/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr, for all classes except some elements of class CostGraph and GraphOperations
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for some elements of class CostGraph and GraphOperations
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.

***************************************************************/

#pragma once

/*!
\file graph.h
\brief Declaration of various classes representing the cost graph and the skeleton graph.
\author Augustin GOUY
*/

#include <stack>
#include <limits>       // for numeric_limits
#include <set>          // for set
#include <utility>      // for pair
#include <algorithm> // for std::sort
#include <fstream>	 // for saving sample file
#include <math.h>
#include <string>		// MSVC2017 requires this
#include <map>
#include <random>
#include <stdlib.h>
#include <queue>
#include <unordered_map>
#include <KarstNSim/nanoflann.hpp>
#include "KarstNSim/surface_sampling.h"
#include "KarstNSim/basics.h"
#include "KarstNSim/vec.h"
#include "KarstNSim/geology.h"
#include "KarstNSim/write_files.h"
#include "KarstNSim/randomgenerator.h"
#include "KarstNSim/simplex_noise.h"

namespace KarstNSim {

	/*!
	\class KarsticConnection
	\brief Represents a connection between nodes in the karst network.

	\details
	A `KarsticConnection` defines the destination node index (`destindex`) for the connection
	and includes an operator to compare the connection's destination index with an integer.

	\public
	- `destindex`: The index of the destination node.
	- `final_branch_id`: The branch ID associated with this connection.
	*/
	class KarsticConnection
	{
	public:
		int destindex; /*!< Index of the destination node. */
		int final_branch_id; /*!< Identifier for the final branch connected by this connection. */
		inline bool operator==(int i) const
		{
			return (destindex == i);
		}

	};

	/*!
	\class KarsticNode
	\brief Represents a node in the karst network.

	\details
	A `KarsticNode` contains information about its position, cost, zone (vadose or phreatic),
	connections to other nodes, and branch identifiers. It also includes utilities for adding
	branch IDs and connections, as well as comparison operators for equality checks.

	\public
	- `index`: The unique identifier of the node.
	- `p`: Position of the node in 3D space.
	- `cost`: A vector of costs associated with the node.
	- `vadose`: A vector indicating the zone of the node (vadose/phreatic).
	- `eq_radius`: Equivalent radius associated with the node.
	- `connections`: A list of `KarsticConnection` objects representing neighboring nodes.
	- `distance`: For debugging; holds a temporary distance value.
	- `branch_id`: Identifier for the branch the node belongs to.
	- `branch_id_ascend`: List of branch IDs from parent nodes, used for backtracking network branches.
	*/
	class KarsticNode
	{
	public:
		int index; /*!< Unique index of the node. */
		Vector3 p; /*!< Position of the node in 3D space. */
		std::vector<float> cost; /*!< Vector of cost values associated with the node, one per water table. */
		std::vector<int> vadose; /*!< Vector of Vadose/phreatic status (one per water table): 1 (vadose), 0 (phreatic), -1 (default). */
		float eq_radius = 0.; /*!< Equivalent radius of the node. */
		std::vector<KarsticConnection> connections; /*!< Connections to neighboring nodes. */
		float distance = 0.; /*!< Debugging aid; holds a temporary distance value. */
		int branch_id; /*!< Identifier for the current branch. */
		std::vector<int> branch_id_ascend; /*!< Branch IDs of parent nodes for backtracking. */

		/*!
		\brief Adds a new branch ID to the node.

		\details
		Combines the current branch ID list (`branch_id_ascend`) with a new set of branch IDs,
		removes duplicates, and ensures the resulting list is sorted.

		\param[in] new_branch_id A vector of new branch IDs to add.
		*/
		inline void add_branch_id(const std::vector<int>& new_branch_id) {
			// Combine both vectors
			branch_id_ascend.insert(branch_id_ascend.end(), new_branch_id.begin(), new_branch_id.end());
			// Sort the combined vector
			std::sort(branch_id_ascend.begin(), branch_id_ascend.end());
			// Remove duplicates
			branch_id_ascend.erase(std::unique(branch_id_ascend.begin(), branch_id_ascend.end()), branch_id_ascend.end());
		}


		inline bool operator==(int i) const
		{
			return (index == i);
		}

		inline bool operator==(Vector3 v) const
		{
			return (p == v);
		}

		/*!
		\brief Adds a connection to another node if it doesn't already exist.

		\details
		This method prevents duplicate connections by first checking if a connection
		to the given destination index already exists.

		\param[in] destindex The index of the destination node.
		\return True if the connection was added, false if it already exists.
		*/
		inline bool add_connection(int destindex)
		{
			for (const KarsticConnection& connection : connections) {
				if (connection.destindex == destindex) {
					return false;
				}
			}
			connections.push_back(KarsticConnection{ destindex });
			return true;
		}
	};

	/*!
	\class CostGraph
	\brief Class representing a nearest neighbor directed cost graph, using an adjacency graph `adj` that consists of a vector of vector of GraphEdge objects.
		   This class also defines methods for updating and pruning the graph, with Dijkstra's algorithm being a key component.
	*/
	class CostGraph
	{
	protected:
		/*!
		\class GraphEdge
		\brief Class representing an edge in the graph, characterized by a target node index, edge weight for both vertices, and fracture flags for both vertices.
		\warning Edges are DIRECTED !! (The cost(s) of a->b may differ from the cost(s) of b->a due to vadose gravity influence (anisotropic))

		*/
		class GraphEdge
		{
		public:
			int target;			//!< The target node index.
			std::vector<float> weight;		//!< The costs of the edge, one per water table.
			std::vector<bool> frac_flag; //!< Flags for indicating whether the edge is part of a fracture's main orientation (UNUSED/DEPRECATED).
		public:
			GraphEdge() : target(-1), weight(), frac_flag() {} // Constructeur par défaut

			/*!
			\brief Constructs a GraphEdge with the specified target node, edge weight, and fracture flags (UNUSED/DEPRECATED).
			\param t The target node index.
			\param w A vector containing the edge's weight for both vertices.
			\param ff A vector indicating the fracture flags for both vertices.
			*/
			inline explicit GraphEdge(int t, std::vector<float> w, std::vector<bool> ff) : target(t), weight(w), frac_flag(ff) {}

			/*!
			\brief Constructs a GraphEdge with the specified target node and edge weight (without fracture flags). (USED VERSION)
			\param t The target node index.
			\param w A vector containing the edge's weight for both vertices.
			*/
			inline explicit GraphEdge(int t, std::vector<float> w) : target(t), weight(w), frac_flag(std::vector<bool>(0)) {}
		};

	protected:

	public:
		/*!
		\brief Constructs a CostGraph with a specified number of nodes.
		\param n The number of nodes in the graph.
		*/
		explicit CostGraph(int n);

		/*!
		\brief The adjacency graph represented as a vector of vector of GraphEdge objects.
		*/
		Array2D<GraphEdge> adj;

	protected:

		/*!
		\brief Normalizes the edge weights in the graph.
		*/
		void NormalizeWeights();

		/*!
		\brief Sets an edge between two nodes with a given weight.
		\param s The source node index.
		\param neigh The neighbor node index.
		\param t The target node index.
		\param w A vector containing the edge weights for both vertices.
		*/
		void SetEdge(const int& s, const int& neigh, const int& t, const std::vector<float>& w);

		/*!
		\brief Sets an edge between two nodes with a given weight and fracture flags.
		\param s The source node index.
		\param neigh The neighbor node index.
		\param t The target node index.
		\param w A vector containing the edge weights for both vertices.
		\param ff A vector indicating the fracture flags for both vertices.
		*/
		void SetEdge(const int& s, const int& neigh, const int& t, const std::vector<float>& w, const std::vector<bool>& ff);

		/*!
		\brief Updates the weight of an edge between two nodes.
		\param s The source node index.
		\param t The target node index.
		\param w A vector containing the updated edge weights for both vertices.
		*/
		void UpdateEdgeWeight(const int& s, const int& t, const std::vector<float>& w);

		/*!
		\brief Retrieves the index of a neighboring node in the adjacency graph.
		\param s The source node index.
		\param t The target node index.
		\return The index of the neighbor node.
		*/
		int GetIdxNeighbor(const int& s, const int& t);

		/*!
		\brief Computes the shortest paths from a source node to all other nodes using Dijkstra's algorithm. In practice, it's used for phreatic shortest path computation between the vadose ending point (source) and the spring (target).
		\param outlet_count The index of the water table associated with the outlet currently considered.
		\param source The index of the inlet node (here, the vadose ending point) in the sampling cloud currently tested.
		\param distance A vector to store the computed shortest distances.
		\param previous A vector to store the previous nodes along the shortest path (can be used to retrieve a path from the source to any other point of the graph).
		\param target Optional parameter to specify a specific target node (the spring) (instead of computing all other nodes, the algorithm will stop as soon as a shortest path is found to this target node) (optional, default value is -1).
		*/
		void DijkstraComputePaths(int outlet_count, int source, std::vector<float>& distance, std::vector<int>& previous, int target = -1) const;

		/*!
		\brief Computes the shortest paths from a source node to the closest point below a given triangulated surface, using Dijkstra's algorithm. In practice, it's used for vadose shortest path computation between the inlet and a vadose ending point (below water table).
		\param outlet_count The index of the water table associated with the outlet currently considered.
		\param source The index of the inlet node in the sampling cloud currently tested.
		\param reach A reference to a variable indicating the first point below the water table level that was reached from the source.
		\param distance A vector to store the computed shortest distances.
		\param previous A vector to store the previous nodes along the shortest path (can be used to retrieve a path from the source to any other point of the graph).
		\param samples_surf_flags Flags indicating whether a node is on the (water table) surface.
		\param target The target node index (the spring).
		\param already_reached Reference to a boolean flag indicating whether the target has already been reached. This corresponds to cases where the phreatic path doesnt have to be computed, since the vadose ending point IS the spring directly.
		*/
		void DijkstraComputePathsSurface(int outlet_count, int source, int& reach, std::vector<float>& distance, std::vector<int>& previous, const  Array2D<char>& samples_surf_flags, int target, bool &already_reached) const;

		/*!
		\brief Retrieves the shortest path and distances from a source node to a target node. Must be used after DijkstraComputePathsSurface or DijkstraComputePaths.
		\param start The index of the start node.
		\param dist A vector containing the computed shortest distances.
		\param prev A vector containing the previous nodes along the shortest path.
		\param cost The total cost of the shortest path.
		\return A pair consisting of the shortest path (as a vector of node indices) and the corresponding distances.
		*/
		std::pair<std::vector<int>, std::vector<float>> DijkstraGetShortestPathTo(int, const std::vector<int>&, const std::vector<float>&, float&) const;
	};

	/*!
	\class GraphOperations
	\brief A class that inherits from `CostGraph` and contains methods for generating a sampling cloud, creating the nearest neighbor graph,
		   assigning costs to graph edges, and pruning the graph to generate the karstic network.
	*/
	class GraphOperations : public CostGraph
	{
	protected:
		/*!
		\struct InternalKeyPoint
		\brief Structure defining internal key points with their index, position, and type.
		*/
		struct InternalKeyPoint
		{
			int index;     //!< The index of the internal key point.
			Vector3 p;     //!< The position of the internal key point.
			KeyPointType type; //!< The type of the key point.
			int wt_idx;    //!< The index of the water table (relevant only for Spring keypoints).
		};

	protected:

		std::vector<float> noise_vector; //!< A vector storing the noise values at each sample point's coordinates.
		Array2D<char> samples_surf_flags; //!< Flags indicating whether each sample is in the phreatic or vadose zone for each water table.
		Array2D<char> samples_on_wt_flags; //!< Flags indicating whether each sample is located exactly on each water table.
		Array2D<float> samples_wt_dist; //!< Distances from each sample to each water table.
		std::vector<float> samples_layer_kp; //!< Intrinsic karstification potential (IKP) value for each sample, scanned from the grid cell containing each point.
		std::vector<float> samples_surf_dist; //!< Distance to the closest inception surface for each sample.
		GeologicalParameters params; //!< Geological parameters for the environment.

	public:

		std::vector<Vector3> samples;	//!< The graph nodes (sampling points). Many other attributes and parameters use indices of points in this list to identify them.

		/*!
		\brief Default constructor for GraphOperations.
		*/
		explicit GraphOperations();

		/*!
		\brief Computes the index of a node given its position.
		\param p The position of the point.
		\return The index of the node corresponding to the given position.
		*/
		int NodeIndex(const Vector3& p) const;

		/*!
		\brief Computes the normalized edge distance cost between two points. It is computed after the other costs.
		\param p The source point.
		\param pn The destination point.
		\param dmin The minimum distance between two points in the whole graph. Used to normalize cost.
		\param dmax The maximum distance between two points in the whole graph. Used to normalize cost.
		\return The computed normalized distance cost.
		*/
		float ComputeEdgeDistanceCost(const Vector3& p, const Vector3& pn, const float& dmin, const float& dmax);

		/*!
		\brief Computes the cost and fracture flags (fracture flags are not used) for an edge between two nodes.
		\param index The index of the node.
		\param p The source point.
		\param pn The destination point.
		\param dmin Reference to the minimum distance between two points in the whole graph. Used to normalize distance cost later.
		\param dmax Reference to the maximum distance between two points in the whole graph. Used to normalize distance cost later.
		\param dist_max The maximum distance of influence of inception surfaces Dmax.
		\return A pair consisting of a vector of computed edge costs and a vector of fracture flags (fracture flags are not used later).
		*/
		std::pair<std::vector<float>, std::vector<bool>> ComputeEdgeCost(const int& index, const Vector3& p, const Vector3& pn, float& dmin, float& dmax, const float&  dist_max) const;


		/*!
		\brief Checks if a point is located below a surface.
		\param pt The point to check.
		\param surface The inception surface to check against.
		\param centers2D The 2D coordinates of the centers of each triangle of the inception surface, stored in a PointCloud.
		\return True if the point is below the surface, otherwise false.
		*/
		static bool CheckBelowSurf(const Vector3&, const Surface*, const PointCloud& centers2D);


		/*!
		\brief Computes the Euclidean distance from a point to a surface.
		\param pt The point to measure the distance from.
		\param surface The surface to measure the distance to.
		\param max_inception_surface_distance Dmax, the maximum inception surface distance of influence.
		\param centers2D The 2D coordinates of the centers of each triangle of the inception surface, stored in a PointCloud.
		\return The computed distance to the surface.
		*/
		float distsurf(const Vector3& pt, Surface* surface, const float& max_inception_surface_distance, const PointCloud& centers2D);

		/*!
		\brief Computes the Euclidean distance of each sample to its closest inception surface (fills samples_surf_dist)
		\param surfaces A vector of all the inception surfaces to compute the distance to.
		\param max_inception_surface_distance Dmax, the maximum inception surface distance of influence.
		*/
		void compute_euclidian_distance_from_surfaces(std::vector<Surface> *surfaces, const float max_inception_surface_distance);

		/*!
		\brief Computes the Euclidean distance of each sample to each water table (fills samples_wt_dist)
		\param surfaces A vector of all the water table surfaces to compute the distance to.
		\param distance Dmax, the maximum inception surface distance of influence.
		*/
		void compute_euclidian_distance_from_wt(std::vector<Surface> *surfaces, const float max_inception_surface_distance);

		/*!
		\brief Saves the nearest neighbor graph (Warning: HEAVY).
		\param geologicalparams The geological parameters for the simulation.
		\param create_nghb_graph_property A flag indicating whether to also save graph properties (and notably, costs of all edges). Warning: VERY HEAVY.
		*/
		void save_nghb_graph(const GeologicalParameters& geologicalparams, const bool& create_nghb_graph_property);

		/*!
		\brief Processes nodes of triangulated surfaces (e.g., inceptions surfaces and water tables), adds them to the samples list and marks the corresponding flag accordingly (e.g., on_surf_flags_idx and on_wt_flags_idx)
		\param node The node to process.
		\param flags_idx A reference to a vector for storing the flags.
		*/
		void process_node_in_samples(const Vector3& node, std::vector<int>& flags_idx);

		/*!
		\brief Initializes the cost graph, sampling the points and processing the key points, creating the nearest neighbor cost graph and properties, and saving it if required.
		\param create_nghb_graph A flag indicating whether to save the nearest neighbor graph.
		\param create_nghb_graph_property A flag indicating whether to save nearest neighbor graph properties.
		\param keypts A vector of key points (inlets, outlets and waypoints) for the graph.
		\param geologicalparams Geological parameters for the environment.
		\param nodes_on_inception_surfaces Nodes on the inception surfaces.
		\param nodes_on_wt_surfaces Nodes on the water table surfaces.
		\param inception_horizons A vector of all the inception surfaces.
		\param water_tables A vector of all the water table surfaces.
		\param use_sampling_points Flag indicating whether to use pre-sampled points instead of generating them.
		\param Box A bounding box for the simulation.
		\param max_inception_surface_distance Maximum distance of influence Dmax to the inception surface.
		\param sampling_points Sampling points (used if use_sampling_points is True).
		\param create_surf_sampling Flag indicating whether to save the sampling cloud.
		\param use_density_property Flag indicating whether to use a density property defined on the background grid.
		\param k_pts Parameter k used in Dwork's algorithm.
		\param fraction_karst_perm Pred, cohesion factor.
		\param propdensity Property density.
		\param propikp Property intrinsic karstification potential.
		\param topo_surface The topographic surface.
		*/
		void InitializeCostGraph(const bool create_nghb_graph, const bool create_nghb_graph_property, const std::vector<KeyPoint>& keypts,
			const GeologicalParameters& geologicalparams, std::vector<Vector3>& nodes_on_inception_surfaces, std::vector<std::vector<Vector3>>& nodes_on_wt_surfaces,
			std::vector<Surface>* inception_horizons,
			std::vector<Surface>* water_tables, const bool use_sampling_points, Box* Box, const float max_inception_surface_distance, std::vector<Vector3>* sampling_points,
			const bool create_surf_sampling, const bool use_density_property, int k_pts, const float fraction_karst_perm, std::vector<float> propdensity,
			std::vector<float> propikp, Surface* topo_surface);


		/*!
		\brief Computes the karstic skeleton for the given key points, and a previously generated cost graph, generating paths for the karst network. To be used after InitializeCostGraph.
		\param pts The key points for the skeleton.
		\param fraction_karst_perm The fraction of karstic permeability.
		\param pathsFinal The final paths of the skeleton.
		\param costsFinal The costs associated with the final paths.
		\param vadoseFinal Vadose zone flags for the final paths (true if in the vadose zone, false in the phreatic zone).
		\param springidxFinal The index of the spring associated with each generated path of the skeleton.
		*/
		void ComputeKarsticSkeleton(const std::vector<KeyPoint>& pts, const float fraction_karst_perm, std::vector<std::vector<int>>& pathsFinal, std::vector<std::vector<float>>& costsFinal, std::vector<std::vector<char>>& vadoseFinal,
			std::vector<int>& springidxFinal);

		/*!
		\brief Finds the shortest path between two nodes in the graph.
		\param u The starting node index.
		\param v The target node index.
		\param idx_specific_wt The index of a specific water table.
		\return A pair consisting of the shortest path as a vector of node indices and the corresponding edge costs.
		*/
		std::pair< std::vector<int>, std::vector<float>> shortest_path(int u, int v, int idx_specific_wt);

		/*!
		\brief amplify the graph from a new set of key points. The new points are linked to the existing base key points.
		\param basekeypts the existing network
		\param newkeypts the new set of key point to connect to the existing network
		\return the set of added paths
		*/
		std::vector<std::vector<int>> amplifykarstickseleton(const std::vector<KarsticNode>& basekeypts, const std::vector<KeyPoint>& newkeypts);

		/*!
		\brief Retrieves the sample at a specified index.
		\param i The index of the sample.
		\return The sample at the specified index.
		*/
		inline Vector3 get_sample(int i) const { return samples[i]; }

		/*!
		\brief Saves the samples to a specified path.
		\param path The path where to save the samples.
		*/
		void save_samples(const std::string& path) const;

		/*!
		\brief Retrieves the number of sampling points in the graph.
		\return The number of sampling points.
		*/
		inline int get_nb_sampling_pts() { return(int)samples.size(); };

		/*!
		\brief Retrieves the indices of the nodes of a previously simulated karst skeleton (used for polyphasic karsts).
		\return A reference to the indices of the nodes of the previously simulated karst skeleton.
		*/
		inline  Array2D<int>& get_IdxOldGraph() { return params.IdxOldGraph; };

		/*!
		\brief Retrieves the flags indicating whether each point in the sampling cloud is in the vadose or phreatic zone, for each water table.
		\return A reference to the the water table flags.
		*/
		inline  Array2D<char>& get_samples_surf_flags() { return samples_surf_flags; };

		/*!
		\brief Adds a ball cost to the graph based on the specified radius. Used during cycle amplification to perturbate costs and force the new path to be different from the previous one, and hence create a cycle (see 2024 article for more details).
		\param center The center of the ball.
		\param ball_radius The radius of the ball.
		\return A vector of pairs, where each pair consists of a node index and its associated cost.
		*/
		std::vector<std::pair<int, float>> add_ball_cost(Vector3& center, float ball_radius);

		/*!
		\brief Adds noise to the sampling points in the graph.
		*/
		void add_noise();

	protected:

		/*!
		\brief Samples the space with the algorithm of Dwork et al. (2021) based on a defined bounding box.
		\param Box The bounding box for the sampling.
		\param use_density_property Flag indicating whether to use a density property or not (if not, density is set as a constant in the whole domain).
		\param k The number of points to sample at each iteration.
		\param propdensity Property density, which is used if use_density_property is True.
		\param topo_surface The topographic surface.
		*/
		void SampleSpaceDwork(Box* Box, const bool use_density_property, const int k, std::vector<float> propdensity, Surface* topo_surface);

		/*!
		\brief Defines flags for the water table surfaces.
		\param water_tables A vector of all the water table surfaces.
		*/
		void DefineWSurfaceFlags(std::vector<Surface>* water_tables);

		/*!
		\brief Builds the nearest neighbor direct cost graph based on the provided parameters.
		\param Box A pointer to the background grid.
		\param fraction_old_karst_perm Ppoly, the factor of cohesion for polyphasic karsts.
		\param max_dist_surf Dmax, the maximum distance of influence of inception surfaces.
		\param keypts A vector of all KeyPoint objects (inlets, outlets, waypoints).
		*/
		void BuildNearestNeighbourGraph(Box*, const float fraction_old_karst_perm, float max_dist_surf, const std::vector<KeyPoint>& keypts);

		/*!
		\brief Saves the noise vector for the given set of points as a new box property.
		\param Points A vector of Vector3 points for which the noise values are to be computed (in principle, the whole samples cloud).
		\return A pair of vectors. The first vector contains the noise values, and the second contains the name of the new box property.
		*/
		std::pair<std::vector<float>, std::vector<std::string>> Save_noise_vector(std::vector<Vector3> Points);
	};

	/*!
	\class KarsticSkeleton
	\brief Class representing a karstic skeleton, defined by its KarsticNode objects, and providing methods for graph amplification and saving.
	\details This class is designed to model a karstic network of interconnected nodes (KarsticNode), compute network properties like branch lengths and valence, and facilitate operations
	such as amplifying the graph with new nodes, detecting intersection points, and saving the network to files.
	*/
	class KarsticSkeleton
	{
	public:
		std::vector<KarsticNode> nodes; //!< The karstic nodes of the skeleton, representing the network.
		Array2D<float> distance_mat; //!< A matrix storing the Euclidean curvilinear distances between all pairs of nodes in the skeleton.
		std::vector<int> valence; //!< A vector representing the number of neighbors (valence) for each node in the skeleton.
		std::vector<int> branch_sizes; //!< A vector storing the sizes (number of nodes) of branches within the karstic network.
		std::vector<int> intersection_points; //!< Flags indicating whether a node is an intersection point (-1 for intersection, 1 for non-intersection).
		int nb_of_intersections = 0; //!< The total number of intersection nodes in the skeleton.

	public:

		/*!
		\brief Default constructor for KarsticSkeleton.
		\details Initializes an empty KarsticSkeleton object with no nodes or associated data.
		*/
		inline explicit KarsticSkeleton() { }

		/*!
		\brief Constructor that initializes the KarsticSkeleton using the GraphOperations sampling cloud and node information.
		\param graph A pointer to GraphOperations containing the full simulation context (and notably, the sampling cloud).
		\param paths A vector of paths representing the connections between nodes.
		\param costs A matrix of costs associated with each path.
		\param vadose_flags A matrix of flags indicating vadose zone status (True if in vadose, False if in phreatic) for each path.
		\param spring_indices A vector of indices for the springs associated with each path.
		*/
		KarsticSkeleton(const GraphOperations* graph, const std::vector<std::vector<int>>& paths, const std::vector<std::vector<float>> costs, const  std::vector<std::vector<char>> vadose, const std::vector<int>& springidx);

		/*!
		\brief Quick constructor when the graph has already been created, useful when the simulation only includes amplification or equivalent section simulation.
		\param pts_graph The points of the edges of the graph.
		\param costs_graph The costs of the edges of the graph.
		\param vadose_flags_graph The vadose flags of the edges of the graph.
		\param spring_indices A vector of indices identifying which springs are connected to each path.
		*/
		KarsticSkeleton(const Array2D<Vector3>& pts_graph, const std::vector<std::vector<float>>& costs_graph, const  std::vector<std::vector<char>>& vadoseflags_graph, const std::vector<int>& springidx); // quick constructor for amplification-only or section generation-only simulations

		/*!
		\brief Computes the number of cycles in the karstic skeleton.
		\return The number of cycles detected in the network.
		*/
		int compute_nb_cycles();

		/*!
		\brief Detects the intersection points in the network.
		\param paths A vector of paths that are evaluated for intersection points.
		*/
		void detect_intersection_points(const std::vector<std::vector<int>>& paths);

		/*!
		\brief Computes the mean branch length in the network.
		\return The average length of all branches in the skeleton.
		*/
		float compute_mean_branch_length();

		/*!
		\brief Computes the mean deviation of branch lengths.
		\return The average deviation in branch lengths across the network.
		*/
		float compute_mean_deviation();

		/**
		* @brief compute the number of points between two water tables (but not on it) over the total number of nodes
		*
		* @param graph
		* @param water_tables the selected water tables
		* @return a ratio between 0 and 1
		*/
		float compute_wt_ratio(GraphOperations* graph, std::vector<Surface>* water_tables);

		/*!
		\brief Counts the number of vadose zone nodes for a specific spring water table.
		\param spring The index of the spring node.
		\return The number of vadose nodes for a specific spring water table.
		*/
		int count_vadose_nodes(int spring);

		/*!
		\brief Computes the average number of vadose nodes in the network.
		\return The average number of vadose nodes.
		*/
		int count_average_vadose_nodes();

		/*!
		\brief Checks if a node is always in the vadose zone for all water tables.
		\param u The node index to check.
		\return True if the node is always in the vadose zone for all water tables, false otherwise.
		*/
		bool always_in_vadose(int u);

		/*!
		\brief Checks if there is an edge between two nodes.
		\param u The first node index.
		\param v The second node index.
		\return True if the edge exists, false otherwise.
		*/
		bool edge_exists(int u, int v);

		/*!
		\brief Amplify an existing network with cycles in the vadose zone (UNUSED).
		\param graph A pointer to GraphOperations used for graph manipulation.
		\param params The geological parameters defining the parameters of cycle generation.
		*/
		void amplify_vadose(GraphOperations* graph, const GeologicalParameters& params);

		/*!
		\brief Amplify an existing network with cycles in the phreatic zone (UNUSED).
		\param graph A pointer to GraphOperations used for graph manipulation.
		\param params The geological parameters defining the parameters of cycle generation.
		*/
		void amplify_phreatic(GraphOperations* graph, const GeologicalParameters& params);

		/*!
		\brief Amplifies the skeleton with deadend points to create new branches.
		\param graph A pointer to GraphOperations used for graph manipulation.
		\param max_distance_of_deadend_pts The maximum distance to consider for deadend points.
		\param nb_deadend_points The number of deadend points to create.
		\param params The geological parameters affecting the amplification.
		*/
		void amplify_deadend(GraphOperations* graph, float max_distance_of_deadend_pts, int nb_deadend_points, const GeologicalParameters& params);

		/*!
		\brief Amplify an existing network with cycles in both vadose and preatic zones (USED VERSION).
		\param graph A pointer to GraphOperations used for graph manipulation.
		\param params The geological parameters defining the parameters of cycle generation.
		\return A pair with the mean and std of the lengths of the cycles created.
		*/
		std::pair<float, float> amplify_noise(GraphOperations* graph, const GeologicalParameters& params);

		/*!
		\brief Saves the karstic skeleton to a file.
		\param file The file path where the skeleton should be saved.
		\param save_directory The directory where the file should be stored.
		*/
		void save(const std::string& file, const std::string& save_directory) const;


		/*!
		\brief
		\param u The parent node index.
		\param v The child node index.
		\return True if u is the parent of v, false otherwise.
		*/
		bool is_parent(int u, int v) const;


		/*!
		\brief Creates a Line object using the karst skeleton. Line objects can the be saved using save_line (in write_files.h).
		\param geologicalparams The geological parameters of the simulation.
		\param network_name The name of the generated network.
		*/
		void create_line(const GeologicalParameters& geologicalparams, std::string network_name) const;

		/*!
		\brief Generates a set of deadend points in the karstic skeleton bounding box. Used in amplify_deadend.
		\param graph A pointer to GraphOperations used for generating deadend points.
		\param max_distance_of_deadend_pts The maximum distance for creating deadend points.
		\param nb_deadend_points The number of deadend points to generate.
		\return A pair of a vector of the simulated deadend points and of a pair of the indices of both these new points and their closest neighbor in the already existing sampling cloud.
		*/
		std::pair<std::vector<Vector3>, std::pair<std::vector<int>, std::vector<int>>> generate_deadend_points(GraphOperations* graph, float max_distance_of_deadend_pts, int nb_deadend_points);

		/*!
		\brief Cleans up the karst skeleton object (e.g., by removing duplicates).
		*/
		void prepare_graph();


		/*!
		\brief Computes the distance matrix between nodes using Johnson's algorithm
		*/
		void compute_distance_matrix();

		/*!
		\brief Computes the size of each branch in the skeleton.
		*/
		void compute_branch_sizes();

		/*!
		\brief Computes the valence (number of neighbors) for each node in the skeleton.
		*/
		void compute_valence();

		/**
		 * @brief Updates the branch IDs for all nodes in the skeleton using a depth-first search (DFS) traversal.
		 *
		 * This method performs the following steps:
		 * - Resets all branch IDs to a default value (-10) for initialization.
		 * - Uses DFS to traverse the skeleton starting from unvisited nodes.
		 * - Assigns a unique branch ID to each branch detected during the traversal.
		 *
		 * A branch is defined as a sequence of connected nodes between intersections or endpoints.
		 * The `branch_id` is incremented for each new branch detected during the traversal.
		 *
		 * This method assumes that the skeleton is represented as a connected graph, where each node has
		 * a list of connections to its neighbors.
		 */
		void update_branch_ID();

		/**
		 * @brief Updates the branch IDs for nodes in the skeleton based on the given paths.
		 *
		 * This method performs the following operations:
		 * - Resets the branch labels for all nodes.
		 * - Refines two-way connections between nodes by removing redundant or "floating" connections.
		 * - Identifies and processes endpoint nodes to prevent incorrect labeling in intersections.
		 * - Traverses the skeleton using a modified depth-first search to assign unique branch IDs, ensuring
		 *   that intersections are correctly handled.
		 *
		 * The branch IDs are assigned based on hierarchical numbering. For instance, at an intersection,
		 * different branches originating from the intersection will have IDs derived from the parent branch ID. See 2025 thesis for more details.
		 *
		 * @param paths A vector of vectors containing paths that define the skeleton's structure.
		 *              Each path is represented as a sequence of node indices.
		 */
		void update_branch_ID(const std::vector<std::vector<int>>& paths);

	private:

		/**
		 * @brief Determines whether one branch ID is a prefix of another. used in is_parent.
		 *
		 * This method checks if the branch ID `id2` is a prefix of the branch ID `id1`.
		 * A prefix relationship is defined such that when `id1` is divided by 10 repeatedly,
		 * it eventually equals `id2`. The method ensures that the larger of the two IDs
		 * is always treated as the starting point of the comparison.
		 *
		 * @param id1 The first branch ID to compare.
		 * @param id2 The second branch ID to compare.
		 * @return `true` if `id2` is a prefix of `id1`, `false` otherwise.
		 *
		 * @note This method swaps the IDs if `id2` is larger than `id1` to simplify the comparison.
		 *       The IDs are treated as numeric sequences, where removing the least significant
		 *       digits from `id1` may reveal a match with `id2`.
		 */
		bool idx_is_prefix(int id1, int id2) const;


		/*!
		\brief Utility function for finding the index of a graph node in the nodes attribute
		\param index index in the sampling cloud
		*/
		int get_internal_index(int index) const;


		/*!
		\brief Utility function for finding the index of a graph node in the nodes attribute
		\param node position of the graph node
		*/
		int get_internal_index(Vector3 node) const;


		/*!
		\brief Appends paths to the skeleton graph.
		\param graph A pointer to the GraphOperations object.
		\param paths A vector of paths to append.
		\param costs A vector of costs associated with the paths.
		\param vadose_flags A matrix of vadose flags for each path.
		\param spring_indices A vector of spring indices for each path.
		*/
		void append_paths(const GraphOperations* graph, const std::vector<std::vector<int>>& paths, const std::vector<std::vector<float>>, const  std::vector<std::vector<char>>, const std::vector<int>& springidx);

		/**
		* @brief Performs a depth-first search (DFS) to traverse and label nodes belonging to the same branch.
		*
		* This method starts from a given node and explores all connected nodes that are part of the same branch,
		* assigning them the specified branch ID (`branch_id_counter`). Nodes are considered part of a branch
		* if they are not intersections (i.e., nodes with more than two connections). Intersection nodes are
		* skipped and assigned a branch ID of `-1`.
		*
		* @param node_index The index of the starting node for the DFS traversal.
		* @param branch_id_counter The branch ID to assign to all nodes within the same branch.
		* @return The number of branch points (nodes assigned the branch ID) encountered during traversal.
		*/
		int dfs(int node_index, int branch_id_counter);
	};
}