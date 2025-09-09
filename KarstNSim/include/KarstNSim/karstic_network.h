/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

/*!
\file karstic_network.h
\brief Declaration of the class that stores input parameters and allows to easily set up the karstic network simulation
\author Benoît THEBAULT, Augustin GOUY
*/

#pragma once

#include "KarstNSim/graph.h"
#include "KarstNSim/vec.h"
#include "KarstNSim/write_files.h"
#include "KarstNSim/randomgenerator.h"
#include "KarstNSim/surface_sampling.h"
#include "KarstNSim/ghost_rocks.h"
#include "KarstNSim/geostats.h"
#include <cmath>
#include <math.h>
#include <iostream>
#include <string>
#include <algorithm>    // std::random_shuffle
#include <fstream>
#include <assert.h>
#include <sstream>
#include <vector>
#include <cstring>

namespace KarstNSim {

	/*!
\struct Propidx
\brief Structure representing a property and its index in the key points vector
*/
	struct Propidx {
		float prop;
		int index; // index in KeyPts vector
	};

	/*!
	\class KarsticNetwork
	\brief Class that stores input parameters and allows to easily set up the karstic network simulation
	*/
	class KarsticNetwork {

	public:

		/*!
		\brief Constructor from a karstic network name, a Box (as boundary box), an object to store all simulation
		parameters, a list of keypoints and water tables
		\param karstic_network_name name of the network we want to simulate
		\param Box boundary box of the domain where we will run the simulation
		\param params class that stores all simulation parameters
		\param keypts contains all keypoints to constrain the simulation
		\param water_tables water tables associated to each spring
		*/

		KarsticNetwork(const std::string& karstic_network_name, Box* Box, GeologicalParameters& params,
			const std::vector<KeyPoint>& keypts, std::vector<Surface>* water_tables);

		/*!
		\brief Set sink points into the key points vector
		\param sinks List of sink points
		\param propsinksindex Indices of sinks in the connectivity matrix
		\param propsinksorder Order of iteration for sinks
		\param use_sinks_radius Flag to indicate whether to use sink radii
		\param propsinksradius Radii of sinks
		*/
		void set_sinks(const std::vector<Vector3>* sinks, const std::vector<int>& propsinksindex, const std::vector<int>& propsinksorder, bool use_sinks_radius, const std::vector<float>& propsinksradius);

		/*!
		\brief Puts all spring points in the keypoint vector
		\param springs Pointer to the list of all spring points to put in the vector
		\param propspringsindex Property representing indices of springs in the connectivity matrix
		\param allow_single_outlet_connection Option that will force each inlet to be only connected to its nearest spring when set to true
		\param use_sinks_radius Flag to indicate whether to use sink radii
		\param propsinksradius Radii of sinks
		\param propspringswtsurf Associated water table surfaces for springs
		*/
		void set_springs(const std::vector<Vector3>* springs, const std::vector<int>& propspringsindex, const bool& allow_single_outlet_connection, bool use_sinks_radius, const std::vector<float>& propsinksradius, const std::vector<int>& propspringswtsurf);

		/*!
		\brief Set way points into the key points vector
		\param waypoints List of way points
		\param use_sinks_radius Flag to indicate whether to use sink radii
		\param propsinksradius Radii of sinks
		\param propwaypointsimpactradius Impact radii for way points
		\param waypoints_weight Weight for way points in simulation
		*/
		void set_waypoints(const std::vector<Vector3>* waypoints, bool use_sinks_radius, const std::vector<float>& propsinksradius, const std::vector<float>& propwaypointsimpactradius, float waypoints_weight);

		/*!
		\brief Set the number and maximum distance of dead-end points
		\param nb_deadend_points Number of dead-end points
		\param max_distance_of_deadend_pts Maximum distance for dead-end points
		*/
		void set_deadend_points(int nb_deadend_points, float max_distance_of_deadend_pts);

		/*!
		\brief Adds previous networks to the generated graph
		\param previous_networks Pointer to the list of lines representing the previous networks incorporated to the graph
		*/
		void set_previous_networks(const std::vector<Line>* previous_networks);

		/*!
		\brief Add nodes from surfaces for densified sampling
		\param network_name Name of the network
		\param surfaces_used_to_densify Surfaces used for densification
		\param refine_surface_sampling Level of surface sampling refinement
		\param create_vset_sampling Flag to create a vertex set sampling
		*/
		void set_inception_surfaces_sampling(const std::string& network_name, std::vector<Surface>* surfaces_used_to_densify,
			const int& refine_surface_sampling, const bool& create_vset_sampling);

		/*!
		\brief Set the parameters for inception horizons
		\param inception_horizons List of inception horizon surfaces
		\param inception_horizon_constraint_weight Weight of the inception horizon constraint
		*/
		void set_wt_surfaces_sampling(const std::string& network_name, std::vector<Surface>* surfaces_used_to_densify,
			const int& refine_surface_sampling);

		/*!
		\brief Set the topographic surface
		\param topo_surface Pointer to the topographic surface
		*/
		void set_topo_surface(Surface* topo_surface);

		/*! \brief Configures ghost rocks for the simulation
		\param grid Grid box representing the domain
		\param ikp Vector of intrinsic karstification potential values
		\param alteration_lines Pointer to surface alteration lines
		\param interpolate_lines Boolean indicating whether to interpolate lines
		\param ghostrock_max_vertical_size Maximum vertical size of ghost rocks
		\param use_max_depth_constraint Boolean indicating whether to use depth constraint
		\param ghost_rock_weight Weight of ghost rocks in the simulation
		\param max_depth_horizon Pointer to maximum depth horizon surface
		\param ghostrock_width Max width of ghost rock corridors
		*/
		void set_ghost_rocks(const Box& grid, std::vector<float>& ikp, Line* alteration_lines, const bool& interpolate_lines, const float& ghostrock_max_vertical_size, const bool& use_max_depth_constraint, const float& ghost_rock_weight, Surface* max_depth_horizon, const float& ghostrock_width);

		/*!
		\brief Sets the simulation with the inception hoizon parameters if we want to use this constraint
		\param inception_horizons Pointer list to surfaces representing the inception horizons
		\param inception_horizon_constraint_weight Weight of the inception horizon constraint
		*/
		void set_inception_horizons_parameters(std::vector<Surface>* inception_horizons, const float& inception_horizon_constraint_weight);

		/*!
		\brief Sets the simulation so that it will not consider inception horizon constraint
		*/
		void disable_inception_horizon();

		/*!
		\brief Sets the simulation if we want to use karstification potential constraint
		\param karstification_potential_weight Weight of the karstification potential constraint
		*/
		void set_karstification_potential_parameters(const float& karstification_potential_weight);


		/*!
		\brief Sets the simulation if we want to use the fracture constraint
		\param fracture_plane_normal Normal vector of the fracture plane (normalized vector)
		\param fracture_intensity intensity/norm of the fracture orienation vector
		\param fracture_constraint_weight Weight of the fracture constraint
		*/
		void set_fracture_constraint_parameters(const std::vector<float>* fracture_families_orientations, const std::vector<float>* fracture_families_tolerance, const float& fracture_constraint_weight);

		/*!
		\brief Sets the simulation so that it will not consider fracture constraint
		*/
		void disable_fractures();

		/*!
		\brief Sets the simulation parameters for the sphere we do not want to add sample points
		\param sphere_centers centers of the spheres
		\param sphere_radius radius pf the spheres
		*/
		void set_no_karst_spheres_parameters(const std::vector<Vector3>* sphere_centers, const std::vector<float>& sphere_radius);

		/*!
		\brief Sets the simulation parameters for the sampling, neireast neighbor and gamma skeleton computation phases
		\param nghb_count Simulation parameter corresponding to the number of neighboor for the n neirest neighboor algorithm
		\param use_max_nghb_radius Boolean to use maximum neighbor radius
		\param nghb_radius Radius used for the n neirest neighboor algorithm
		\param poisson_radius Radius used for the Poisson sphere distribution (sampling phase)
		\param gamma Parameter used for the beta-skeleton appraoch to build the network
		\param multiply_costs Boolean to indicate whether to multiply costs (instead of summing them)
		\param vadose_cohesion Boolean for vadose cohesion consideration: if true, cohesion is applied everywhere, if false, cohesion is only applied in the phreatic zone
		*/
		void set_simulation_parameters(const int& nghb_count, const bool& use_max_nghb_radius, const float& nghb_radius, const float& poisson_radius, const float& gamma, const bool& multiply_costs, const bool& vadose_cohesion);

		/*!
		\brief Sets the domain geometry and the boundary box dimensions with the Box class attribute
		*/
		void set_domain_geometry();

		/*! \brief Runs the sampling of the domain and returns the number of sampled points
		\return Number of sampled points
		*/
		int just_sampling();

		/*!
		\brief Reads the connectivity matrix file (Rows = sinks, Columns = springs, 0 = no link, 1 = link, 2 = random) and stores info in a vector of vectors
		\param sinks Pointer to the list of sinks
		\param springs Pointer to the list of springs
		*/
		void read_connectivity_matrix(const std::vector<Vector3>* sinks, const std::vector<Vector3>* springs);

		/*! \brief Sets noise parameters for the simulation
		\param use_noise Boolean indicating whether to use noise for the cycle amplification step only
		\param use_noise_on_all Boolean indicating whether to apply noise for the whole simulation
		\param frequency Frequency of the noise
		\param octaves Number of octaves for noise generation
		\param noise_weight Weight of the noise
		\param globalRng Random number generator for noise application
		*/
		void set_noise_parameters(const bool use_noise, const bool use_noise_on_all, const int frequency, const int octaves, const float noise_weight, std::mt19937 globalRng);

		/*! \brief Creates sections of the karstic skeleton using the 1D-curvilinear branchwise SGS algorithm.
		\param skel Reference to the KarsticSkeleton object
		*/
		void create_sections(KarsticSkeleton& skel);

		/*! \brief Runs a KarstNSim simulation, by skipping ALL steps except the equivalent section simulation using the 1D-curvilinear branchwise SGS algorithm.
		\param skel Reference to the KarsticSkeleton object
		\param alteration_lines Pointer to alteration lines
		\param use_ghost_rocks Boolean indicating whether to use ghost rocks
		\param ghostrock_max_vertical_size Maximum vertical size of ghost rocks
		\param use_max_depth_constraint Boolean indicating whether to use maximum depth constraint
		\param max_depth_horizon Pointer to maximum depth horizon surface
		\param ghostrock_width Width of ghost rocks
		*/
		void run_simulation_properties(KarsticSkeleton& skel, Line* alteration_lines, const bool& use_ghost_rocks, const float& ghostrock_max_vertical_size, const bool& use_max_depth_constraint, Surface* max_depth_horizon, const float& ghostrock_width);

		/*! \brief Runs the simulation and returns the time required for the simulation
		\param sections_simulation_only Boolean indicating whether to simulate sections only (and skip everything else)
		\param create_nghb_graph Boolean to create nearest neighbor graph
		\param create_nghb_graph_property Boolean to create nearest neighbor graph property
		\param use_amplification Boolean to use amplification phase
		\param use_sampling_points_ Boolean to use sampling points
		\param fraction_karst_perm Cohesion factor Pred
		\param fraction_old_karst_perm Polyphasic cohesion factor Ppoly
		\param max_inception_surface_distance Maximum distance of influence of inception surfaces Dmax
		\param sampling_points_ Vector to store sampling points
		\param create_vset_sampling_ Boolean to create sampling pointset
		\param use_density_property_ Boolean to use density property
		\param k_pts Number of key points
		\param propdensity Vector of density property of the box
		\param propikp Vector of IKP property of the box
		\return Simulation time
		*/
		float run_simulation(const bool& sections_simulation_only, const bool& create_nghb_graph, const bool& create_nghb_graph_property, const bool& use_amplification, const bool& use_sampling_points_,
			const float& fraction_karst_perm, const float& fraction_old_karst_perm, const float& max_inception_surface_distance, std::vector<Vector3>* sampling_points_, const bool& create_vset_sampling_,
			const bool& use_density_property_, const int& k_pts, const std::vector<float>& propdensity, const std::vector<float>& propikp);

		/*!
		\brief Sets the save repertory if we want to create link and node files
		\param repertory save repertory
		*/
		void set_save_directory(const std::string& repertory);

		/*!
		\brief Saves the box painted with the properties of the simulation. Mainly useful if ghost rocks were painted onto the IKP property, since it modifies it.
		\param propdensity density property
		\param propikp IKP property, potentially updated during simulation by the ghost-rock corridors.
		*/
		void save_painted_box(std::vector<float>& propdensity, const std::vector<float>& propikp);

		/*!
		\brief Sets all the geostatistical parameters needed for the simulation of equivalent section properties on the nodes of the karstic skeleton
		\param geostat_params struct defining all geostatistical parameters
		*/
		void set_geostat_params(const GeostatParams& geostat_params);

		/*!
		\brief Sets the simulation parameters for the cycle generation part of the amplification step
		\param max_distance_amplification maximum distance between two nodes randomly selected to create a cycle
		\param min_distance_amplification minimum distance between two nodes randomly selected to create a cycle
		\param nb_cycles number of cycles to be generated
		*/
		void set_amplification_params(const float& max_distance_amplification, const float& min_distance_amplification, const int& nb_cycles);

		/*!
		\brief Sets the simulation parameters for the cycle generation part of the amplification step for the vadose zone (CURRENTLY UNUSED)
		\param max_distance_loops_vadose maximum distance between two nodes randomly selected to create a cycle in the vadose zone
		\param loop_density_vadose density of loops in the vadose zone
		*/
		void set_amplification_vadose_params(const float& max_dist_loops_vadose, const float& loop_density_vadose);

		/*!
		\brief Sets the simulation parameters for the cycle generation part of the amplification step for the phreatic zone (CURRENTLY UNUSED)
		\param max_distance_loops_phreatic maximum distance between two nodes randomly selected to create a cycle in the phreatic zone
		\param loop_density_phreatic density of loops in the phreatic zone
		*/
		void set_amplification_phreatic_params(const float& max_dist_loops_phreatic, const float& loop_density_phreatic);

		/*!
		\brief Sets water table weight for vadose and phreatic zones.
		\param water_table_constraint_weight_vadose weight of vadose zone subcost
		\param water_table_constraint_weight_phreatic weight of phreatic zone subcost
		*/
		void set_water_table_weight(const float& water_table_constraint_weight_vadose, const float& water_table_constraint_weight_phreatic);

		/*!
		\brief Sets the simulation so that it will not consider water table cost constraints
		*/
		void disable_water_table();


	private:
		std::string karstic_network_name; /*!< Name of the karstic network to be generated */
		std::vector<Vector3> nodes_on_inception_surfaces; /*!< List of nodes on inception surfaces */
		std::vector<std::vector<Vector3>> nodes_on_wt_surfaces; /*!< List of nodes on water table surfaces */
		std::vector<Vector3> pt_sink; /*!< Sinks */
		std::vector<Vector3> pt_spring; /*!< Springs */
		std::vector<Surface>* inception_horizons = nullptr; /*!< Inception horizon surfaces */
		std::vector<Surface>* water_tables = nullptr; /*!< Water table surfaces */
		Surface* topo_surface_ = nullptr; /*!< Topographic surface */
		bool is_simulation_parametrized = false; /*!< Flag indicating whether the simulation is parametrized or not*/
		bool use_deadend_pts_ = false; /*!< Flag indicating whether deadend amplification should be applied */
		int nb_deadend_points_; /*!< Number of dead-end points */
		float max_distance_of_deadend_pts_; /*!< Maximum distance for dead-end points */
		Box* box; /*!< Bounding box of the simulation */
		GeologicalParameters params; /*!< Geological parameters */
		GeostatParams geostatparams; /*!< Geostatistical parameters */
		std::vector<KeyPoint> keypts; /*!< Keypoints for simulation */
		std::string save_directory; /*!< Directory for saving outputs */

		bool use_sinks_radius_; /*!< Flag to use sinks radius property */
		std::vector<Propidx> propsinksradius_; /*!< Radii of sinks property */
		bool use_springs_radius_; /*!< Flag to use springs radius property */
		std::vector<Propidx> propspringsradius_; /*!< Radii of springs property */
		bool use_waypoints_radius_; /*!< Flag to use waypoints radius property */
		std::vector<Propidx> propwaypointsradius_; /*!< Radii of waypoints property */

	};
}
