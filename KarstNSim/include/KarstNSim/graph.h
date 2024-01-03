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
#include "KarstNSim/KarsticNetwork.h"
#include "KarstNSim/write_files.h"
#include "KarstNSim/randomgenerator.h"

namespace KarstNSim {
	/*!
	\class KarsticConnection
	\brief Class that defines the index of a node destIndex
	*/
	class KarsticConnection
	{
	public:
		int destIndex;

		inline bool operator==(int i) const
		{
			return (destIndex == i);
		}
	};

	/*!
	\class KarsticNode
	\brief Class that defines nodes of the karst network, with their index, position, cost, neighbor indices and whether they are in the vadose or phreatic zone
	*/
	class KarsticNode
	{
	public:
		int index;
		Vector3 p;
		double cost;
		bool vadose; // 1 if in vadose, 0 if in phreatic
		std::vector<KarsticConnection> connections;


		inline bool operator==(int i) const
		{
			return (index == i);
		}
	};

	class GraphOperations;

	/*!
	\class KarsticSkeleton
	\brief Class containing definition of a karstic skeleton object, represented by its KarsticNode. It also contains several methods to amplify the karst network and save it
	*/
	class KarsticSkeleton
	{
	public:
		std::vector<KarsticNode> nodes;

	public:
		inline explicit KarsticSkeleton() { }
		KarsticSkeleton(const GraphOperations*, const std::vector<std::vector<int>>&, const std::vector<std::vector<double>>, const std::vector<std::vector<bool>>);

		int count_vadose_nodes();
		bool edge_exists(int u, int v);
		void Amplify_vadose(GraphOperations* graph, GeologicalParameters params);
		void Amplify_phreatic(GraphOperations* graph, GeologicalParameters params);
		void save(const std::string& file, const std::string& save_directory) const;
		void CreateLine(const GeologicalParameters& geologicalparams, std::string network_name) const;

	private:
		int GetInternalIndex(int index) const;
		void AppendPaths(const GraphOperations* graph, const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>, const std::vector<std::vector<bool>>);
	};

	/*!
	\class CostGraph
	\brief Class that defines the nearest neighbor cost graph, represented by the adjcency graph adj, which is a vector of vector of GraphEdge. Methods pertaining to update and prune this graph are also defined (Dijsktra algorithm).
	*/
	class CostGraph
	{
	protected:
		/*!
		\class GraphEdge
		\brief Class that defines a graph edge as a target node index, a weight for both vertices of the edge and a fracture flag for both vertices of the edge, that tells if the edge is in a fracture main orientation or not
		*/
		class GraphEdge
		{
		public:
			int target;			//!< Target node index.
			std::vector<double> weight;		//!< Weight of the edge.
			std::vector<bool> frac_flag; //!< 
		public:
			inline explicit GraphEdge(int t, std::vector<double> w, std::vector<bool> ff) : target(t), weight(w), frac_flag(ff)
			{
			}
			inline explicit GraphEdge(int t, std::vector<double> w) : target(t), weight(w), frac_flag(std::vector<bool>(0))
			{
			}
		};

	protected:

	public:
		explicit CostGraph(int n);
		std::vector<std::vector<GraphEdge>> adj; //!< Adjacency graph. Attributes : target, weight
		std::vector<std::vector<GraphEdge>> adj_no_cost_red;

	protected:
		void NormalizeWeights();
		void SetEdge(int s, int t, std::vector<double> w);
		void SetEdge(int s, int t, std::vector<double> w, std::vector<bool> ff);
		void UpdateEdgeWeight(int s, int t, std::vector<double> w);
		int GetIdxNeighbor(int s, int t);
		void DijkstraComputePaths(int, int, std::vector<double>&, std::vector<int>&) const;
		void DijkstraComputePathsNoCostRed(int, int, std::vector<double>&, std::vector<int>&) const;
		void DijkstraComputePathsSurface(int, int, int&, std::vector<double>&, std::vector<int>&, std::vector<std::vector<bool>>, int, int, bool&) const;
		std::pair<std::vector<int>, std::vector<double>> DijkstraGetShortestPathTo(int, const std::vector<int>&, const std::vector<double>&, double&) const;
	};

	struct Neighbour {
		int i;
		double d;
		bool operator<(const Neighbour& nei) const {
			return d < nei.d;
		}
		bool operator>(const Neighbour& nei) const {
			return d > nei.d;
		}
	};

	struct PointCloud {

		const std::vector<Vector3>& cloud;
		using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3>;
		std::unique_ptr<my_kd_tree_t> index;

		PointCloud(const std::vector<Vector3>& cloud);

		// Must return the number of data points
		inline size_t kdtree_get_point_count() const { return cloud.size(); }

		// Returns the dim'th component of the idx'th point in the class:
		// Since this is inlined and the "dim" argument is typically an immediate
		// value, the
		//  "if/else's" are actually solved at compile time.
		inline double kdtree_get_pt(const size_t idx, const size_t dim) const
		{
			if (dim == 0)
				return cloud[idx].x;
			else if (dim == 1)
				return cloud[idx].y;
			else
				return cloud[idx].z;
		}

		// Optional bounding-box computation: return false to default to a standard
		// bbox computation loop.
		//   Return true if the BBOX was already computed by the class and returned
		//   in "bb" so it can be avoided to redo it again. Look at bb.size() to
		//   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
		template <class BBOX>
		bool kdtree_get_bbox(BBOX& /* bb */) const
		{
			return false;
		}

		// Function to find nearest neighbors
		void findNearestNeighbors(const Vector3& queryPoint, int querryidx, int N, double distmax, std::vector<Neighbour>& result);

	};

	/*!
	\class GraphOperations
	\brief Class inheriting from CostGraph, containing methods allowing to create the sampling cloud, generate the nearest neighbor graph onto it, apply costs to each edge and prune the graph to generate the karstic network
	*/
	class GraphOperations : public CostGraph
	{
	protected:
		/*!
		\struct InternalKeyPoint
		\brief Structure defining internal key points, by their index, position and type
		*/
		struct InternalKeyPoint
		{
			int index;
			Vector3 p;
			KeyPointType type;
		};
		int NodeIndex(const Vector3& p) const;

	protected:
		std::vector<Vector3> samples;	//!< Graph nodes.
		std::vector<std::vector<bool>> samples_surf_flags;
		std::vector<std::vector<double>> samples_wt_dist;
		std::vector<double> samples_layer_kp;
		std::vector<double> samples_surf_dist;
		GeologicalParameters params;	//!< All geological parameters.

	public:
		explicit GraphOperations();
		double ComputeEdgeDistanceCost(const Vector3& p, const Vector3& pn, const double dmin, const double dmax);
		std::pair<std::vector<double>, std::vector<bool>> ComputeEdgeCostAdapted(const int index, const Vector3& p, const Vector3& pn, double& dmin, double& dmax, const double  dist_max, std::vector<double> max_dist_wt) const;
		bool CheckBelowSurf(Vector3, Surface*);
		double distsurf(Vector3 pt, Surface* surface, const double max_inception_surface_distance);
		double compute_euclidian_distance_from_surfaces(std::vector<Surface>*, const double);
		std::vector<double> compute_euclidian_distance_from_wt(std::vector<Surface>*, const double);
		void InitializeCostGraphAdapted(const bool create_nghb_graph, const bool create_nghb_graph_property, const std::vector<KeyPoint>& keypts, const GeologicalParameters& geologicalparams, const std::vector<Vector3>& Points, std::vector<Surface>* inception_horizons,
		std::vector<Surface>* water_tables, const bool use_sampling_points, Box* Box, const double max_inception_surface_distance, std::vector<Vector3>* sampling_points,
		const bool create_surf_sampling, const bool use_density_property, int k_pts, const double fraction_karst_perm, std::vector<double> propdensity, std::vector<double> propikp, Surface* topo_surface);
		KarsticSkeleton ComputeKarsticSkeleton(const std::vector<KeyPoint>& keypts, const double fraction_old_karst_perm);
		std::pair< std::vector<int>, std::vector<double>> shortest_path(bool under_wt, int u, int v);
		std::vector<std::vector<int>> AmplifyKarsticSkeleton(const std::vector<KarsticNode>& basekeypts, const std::vector<KeyPoint>& newkeypts);
		std::vector<InternalKeyPoint> AddNewSamples(const std::vector<KeyPoint>& samples);
		inline Vector3 GetSample(int i) const { return samples[i]; }
		void saveSamples(const std::string& path) const;
		inline int get_nb_sampling_pts() { return(int)samples.size(); };

	protected:
		void SampleSpaceDwork(std::vector<Vector3>& Points, Box* Box, const bool use_density_property, const int k, std::vector<double> propdensity, Surface* topo_surface);
		void DefineWSurfaceFlags(std::vector<Surface>* water_tables);
		void BuildNearestNeighbourGraphAdapted(Box*, const double fraction_old_karst_perm, const std::vector<double> max_dist_wt, double max_dist_surf);
	};
}
