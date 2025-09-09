/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for classes KeyPointType, KeyPoint and CostTerm
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.

***************************************************************/

#pragma once

/*!
\file geology.h
\brief Declaration of various classes representing simple geological objects and parameters
\author Augustin GOUY
*/

#include "basics.h"
#include "graph.h"
#include <string>

namespace KarstNSim {
	/*!
	\class KeyPointType
	\brief Class that defines different keypoint types
	*/
	enum class KeyPointType
	{
		Sink,		//!< Network inlet
		Spring,		//!< Network outlet
		Waypoint,	//!< Known passage in the bedrock
		Deadend, //!< Deadend point sampled during amplification
		sampling, //!< Sampling point
	};

	/*!
	\class KeyPoint
	\brief Class that defines a keypoint, by its position and KeyPointType
	*/
	class KeyPoint
	{
	public:
		Vector3 p; //!< position of KeyPoint
		KeyPointType type; //!< KeyPointType of KeyPoint
		int wt_idx; //!< Only relevant for Spring KeyPointType. Keeps track of the index of each spring keypoint (to couple them with their water table)
	public:

		/*!
			\brief Default constructor of the KeyPoint class.
			Creates an uninitialized KeyPoint object. The attributes must be set explicitly after using this constructor.
		*/
		inline explicit KeyPoint() { }
		/*!
	\brief Constructor of the KeyPoint class with position and type.

	\param p The position of the KeyPoint in 3D space.
	\param t The type of the KeyPoint, as defined by the KeyPointType enum.
*/
		inline explicit KeyPoint(const Vector3& p, KeyPointType t) : p(p), type(t) {}
		/*!
	\brief Constructor of the KeyPoint class with position, type, and index.

	This constructor is primarily used for spring-type keypoints that require an associated index.

	\param p The position of the KeyPoint in 3D space.
	\param t The type of the KeyPoint, as defined by the KeyPointType enum.
	\param idx The index associated with spring-type keypoints for coupling with a water table.
*/
		inline explicit KeyPoint(const Vector3& p, KeyPointType t, const int& idx) : p(p), type(t), wt_idx(idx) {}
	};

	/*!
		\class CostTerm
		\brief Class that defines cost terms, representing the subcomponents of a cost function.

		The CostTerm class is used to represent individual subcosts of a cost function.
		Each cost term can be enabled or disabled using the `used` attribute, and its relative importance
		in the cost function is specified by the `weight` attribute.

		\details Cost terms are building blocks of the overall cost function.
		By setting the `used` flag, users can selectively include or exclude specific cost terms.
		The `weight` determines the contribution of the cost term to the overall cost calculation when it is enabled.
	*/
	class CostTerm
	{
	public:
		bool used; //!< Indicates whether the cost term is active or not.
		float weight; //!< The weight of the cost term, influencing its contribution to the overall cost function.

	public:
		/*!
	\brief Default constructor of the CostTerm class.

	Initializes the cost term with default values. By default, the cost term is not used (`used = false`)
	and its weight is set to zero (`weight = 0.0`).
*/
		inline explicit CostTerm() : used(false), weight(0.0) { }
		/*!
	\brief Parameterized constructor of the CostTerm class.

	Allows the initialization of a cost term with specific values for `used` and `weight`.

	\param u A boolean value indicating whether the cost term is active.
	\param w A float value specifying the weight of the cost term.
*/
		inline explicit CostTerm(bool u, float w) : used(u), weight(w) { }
	};

	/*!
	\class GeologicalParameters
	\brief Class that defines all geological and simulation parameters, used throughout the code.
	*/
	class GeologicalParameters
	{
	public:
		std::string scenename; //!< Name of the simulation (will appear as a prefix in all output files).
		std::string directoryname; //!< Directory where simulation data is stored.

		bool use_amplification; //!< Flag to enable amplification step.
		float max_distance_amplification; //!< Maximum amplification distance to create a cycle.
		float min_distance_amplification; //!< Minimum amplification distance to create a cycle.
		int nb_cycles; //!< Number of cycles to simulate during amplification.
		bool use_noise; //!< Flag to enable noise generation in the simulation.
		bool use_noise_on_all; //!< Flag to apply noise globally across the entire simulation area.
		int noise_frequency; //!< Frequency of the noise pattern.
		int noise_octaves; //!< Number of octaves used in the noise generation.
		float noise_weight; //!< Weighting factor for noise influence on the simulation.

		float graphPoissonRadius; //!< Poisson sphere sampling radius for graph construction.
		float graphNeighbourRadius; //!< Max radius to consider neighboring points in the graph.
		float maxsize; //!< Average of the dx and dy sizes of the background grid.
		float stretch_factor; //!< Stretching factor of the grid, corresponding to the ratio of maxsize on the dz size of the grid. Used for Nghb computations.
		bool graphuse_max_nghb_radius; //!< Flag to use a maximum neighbor radius constraint graphNeighbourRadius in graph construction.
		int graphNeighbourCount; //!< Number of neighbors to consider for graph connections.
		int nb_springs; //!< Number of springs in the simulation.
		int nb_wt; //!< Number of water tables in the simulation.
		int nb_inception_surf = 0; //!< Number of inception surfaces in the simulation.

		bool multiply_costs; //!< Flag to multiply cost terms during cost function computation instead of adding them.
		bool allow_single_outlet; //!< Allow a connection to a single outlet for each inlet, not more. This will use the ``closest'' spring algorithm (see Thesis for details)
		bool vadose_cohesion; //!< Flag to enable vadose zone cohesion in the simulation (cohesion only in phreatic zone if set to false).

		Array2D<Vector3> PtsOldGraph; //!< Coordinates of points from a previously simulated karst network. Used for polyphasic karstification.
		Array2D<int> IdxOldGraph; //!< Indices of points from a previously simulated karst network. Used for polyphasic karstification.
		std::vector<int> sinks_index; //!< Indices of sink points in the samples vector (main attribute of GraphOperations).
		Array2D<int> connectivity_matrix; //!< Connectivity matrix for inlet/outlet relationships.

		std::vector<float> fractures_orientations; //!< Orientations of fractures in the geological structure.
		std::vector<float> fractures_tolerances; //!< Tolerances for fracture orientations.
		std::vector<float> fractures_max_lengths; //!< Maximum lengths of fractures (CURRENTLY UNUSED).

		float max_dist_loops_vadose; //!< Maximum loop distance in the vadose zone (CURRENTLY UNUSED).
		float loop_density_vadose; //!< Density of loops in the vadose zone (CURRENTLY UNUSED).
		float max_dist_loops_phreatic; //!< Maximum loop distance in the phreatic zone (CURRENTLY UNUSED).
		float loop_density_phreatic; //!< Density of loops in the phreatic zone (CURRENTLY UNUSED).

		bool use_ghost_rocks = false; //!< Flag to enable ghost rocks in the simulation.
		float length; //!< Length of the ghost-rock corridor.
		float width; //!< Width of the ghost-rock corridor.
		Line polyline; //!< Polyline representing the surface alteration lines.
		bool use_max_depth_constraint; //!< Flag to enforce a maximum depth constraint in substratum_surf for ghost-rocks corridors.
		Surface substratum_surf; //!< Surface representing the substratum for ghost-rocks corridors.

		CostTerm distanceCost; //!< Cost term for Euclidean distance.
		CostTerm fractureCost; //!< Cost term for fracture influence.
		CostTerm horizonCost; //!< Cost term for inception horizons (and structural faults) influence.
		CostTerm waterTable1; //!< Cost term associated with the vadose zone (in the 2024 article and the thesis, waterTable1 and waterTable2 are grouped together as a single cost. In practice, users can define them separately if needed).
		CostTerm waterTable2; //!< Cost term associated with the phreatic zone (in the 2024 article and the thesis, waterTable1 and waterTable2 are grouped together as a single cost. In practice, users can define them separately if needed).
		CostTerm karstificationCost; //!< Cost term for Intrinsic Karstification Potential (factors permeability, porosity, solubility of rocks as well as presence of ghost-rock weathering).

		float gamma; //!< factor used to apply rules of gamma-graph in a simulated network. See publication of Paris et al. (2021) for more details. It was not used in the 2024 article and the thesis. Values should in principle be between 0 and 2, 0 excluded and 2 included.

		std::vector<Sphere> spheres; //!< Collection of no_karst_spheres.

		struct Propidx {
			float prop;
			int index; // index in KeyPts vector
		};

		float waypoints_weight; //!< Weighting factor for waypoint contributions in the simulation.
		std::vector<Propidx> waypointsimpactradius; //!< List of impact radii of each waypoint.
		std::vector<Propidx> z_list; //!< List of Z-coordinates of each spring.
		std::vector<Propidx> propspringswtindex; //!< List of index of water table associated with each spring.

	};
}