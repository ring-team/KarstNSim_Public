/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

/*!
\file KarsticNetwork.h
\brief Declaration of the class that stores input parameters and allows to easily set up the karstic network simulation
\author Benoît THEBAULT, Augustin GOUY
*/

#pragma once

#include "KarstNSim/graph.h"
#include "KarstNSim/vec.h"
#include "KarstNSim/write_files.h"
#include "KarstNSim/randomgenerator.h"
#include "KarstNSim/surface_sampling.h"
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
			const std::vector<KeyPoint>& keypts,  std::vector<Surface>* water_tables);

		/*!
		\brief Puts all sink points in the keypoint vector
		\param sinks Pointer to the list of all sink points to put in the vector
		\param propsinksindex Property representing indices of sinks in the connectivity matrix
		\param propsinksorder Property representing order of iteration of sinks.
		*/
		void set_sinks(const std::vector<Vector3>* sinks, const std::vector<int>& propsinksindex, const std::vector<int>& propsinksorder);

		/*!
		\brief Puts all spring points in the keypoint vector
		\param springs Pointer to the list of all spring points to put in the vector
		\param propspringsindex Property representing indices of springs in the connectivity matrix
		\param allow_single_outlet_connection Option that will force each inlet to be only connected to its nearest spring when set to true
		*/
		void set_springs(const std::vector<Vector3>* springs, const std::vector<int>& propspringsindex, const bool& allow_single_outlet_connection);

		/*!
		\brief Precises index of each water tables surface
		\param  surf_water_table all surfaces of water tables
		*/
		void set_water_table_order(const std::vector<Surface>* surf_water_table);

		/*!
		\brief Puts all way points in the keypoint vector
		\param waypoints Pointer to the list of all way points to put in the vector
		*/
		void set_waypoints(const std::vector<Vector3>* waypoints);


		/*!
		\brief Adds previous networks to the generated graph
		\param previous_networks Pointer to the list of lines representing the previous networks incorporated to the graph
		*/
		void set_previous_networks(const std::vector<Line>* previous_networks);

		/*!
		\brief Adds all nodes of the surfaces in the sampling point list given as input where we want to densify sampling points
		\param surfaces_to_densify Pointer list of the surfaces whose nodes will be included in the sampling point list
		*/
		void set_inception_surfaces_sampling(const std::string& network_name, std::vector<Surface>* surfaces_used_to_densify,
			const int& refine_surface_sampling, const bool& create_vset_sampling);

		void set_wt_surfaces_sampling(const std::string& network_name, std::vector<Surface>* surfaces_used_to_densify,
			const int& refine_surface_sampling);


		void set_topo_surface(Surface* topo_surface);


		/*!
		\brief Sets the simulation with the inception hoizon parameters if we want to use this constraint
		\param inception_horizons Pointer list to surfaces representing the inception horizons
		\param inception_horizon_constraint_weight Weight of the inception horizon constraint
		*/
		void set_inception_horizons_parameters(std::vector<Surface>* inception_horizons, const double& inception_horizon_constraint_weight);

		/*!
		\brief Sets the simulation so that it will not consider inception horizon constraint
		*/
		void disable_inception_horizon();

		/*!
		\brief Sets the simulation if we want to use karstification potential constraint
		*/
		void set_karstification_potential_parameters(const double& karstification_potential_weight);


		/*!
		\brief Sets the simulation if we want to use the fracture constraint
		\param fracture_plane_normal Normal vector of the fracture plane (normalized vector)
		\param fracture_intensity intensity/norm of the fracture orienation vector
		\param fracture_constraint_weight Weight of the fracture constraint
		*/
		void set_fracture_constraint_parameters(const std::vector<double>* fracture_families_orientations, const std::vector<double>* fracture_families_tolerance, const double& fracture_constraint_weight);

		/*!
		\brief Sets the simulation so that it will not consider fracture constraint
		*/
		void disable_fractures();

		/*!
		\brief Sets the simulation parameters for the sphere we do not want to add sample points
		\param sphere_centers centers of the spheres
		\param sphere_radius radius pf the spheres
		*/
		void set_no_karst_spheres_parameters(const std::vector<Vector3>* sphere_centers, const std::vector<double>& sphere_radius);

		/*!
		\brief Sets the simulation parameters for the sampling, neireast neighbor and gamma skeleton computation phases
		\param nghb_count Simulation parameter corresponding to the number of neighboor for the n neirest neighboor algorithm
		\param nghb_radius Radius used for the n neirest neighboor algorithm
		\param poisson_radius Radius used for the Poisson sphere distribution (sampling phase)
		\param gamma Parameter used for the beta-skeleton appraoch to build the network
		*/
		void set_simulation_parameters(const int& nghb_count, const bool& use_max_nghb_radius, const double& nghb_radius, const double& poisson_radius, const double& gamma, const bool& multiply_costs, const bool& vadose_cohesion);

		/*!
		\brief Sets the domain geometry and the boundary box dimensions with the Box class attribute
		*/
		void set_domain_geometry();

		/*!
		\brief Runs the sampling of the domain and returns the number of sampled points
		*/
		int just_sampling();

		/*!
		\brief Reads the connectivity matrix file (Rows = sinks, Columns = springs, 0 = no link, 1 = link, 2 = random) and stores info in a vector of vectors
		*/
		void read_connectivity_matrix(const std::vector<Vector3>* sinks, const std::vector<Vector3>* springs);

		/*!
		\brief Runs the simulation and returns the number of sampling points
		\param save_network Boolean allowing to inform if we want to save the simulated network in node and link files
		\param use_amplification_phase Boolean allowing to inform if we want to use the amplification phase during the simulation
		\param create_line Boolean allowing to inform if we want to directely create the network in SKUA-GOCAD as Line
		\param fracture_distribution Boolean allowing to inform if we want to give fracture distributions to create the fracture
		orientation vector(s) instead of directely use an orienation vector (if we use the fracture constraint)
		*/

		void run_simulation(const bool& create_nghb_graph, const bool& create_nghb_graph_property, const bool& use_amplification_phase_vadose, const bool& use_amplification_phase_phreatic, const bool& use_sampling_points_,
			const double& fraction_karst_perm, const double& fraction_old_karst_perm, const double& max_inception_surface_distance, std::vector<Vector3>* sampling_points_, const bool& create_vset_sampling_,
			const bool& use_density_property_, const int& k_pts, const std::vector<double>& propdensity, const std::vector<double>& propikp);

		/*!
		\brief Sets the save repertory if we want to create link and node files
		\param repertory save repertory
		*/
		void set_save_directory(const std::string& repertory);

		void save_painted_box(const std::vector<double>& propdensity, const std::vector<double>& propikp);

		/*!
		\brief Sets the simulation parameters for the amplification phase (number of new points)
		\param nb_new_deadend_points deadend points (see the method explanations in Paris et al. 2021)
		\param nb_interior_nodes deadend points (see the method explanations in Paris et al. 2021). NB: not yet functional,
		\ same effect as deadend points for now.
		*/
		void set_amplification_vadose_params(const double& max_dist_loops_vadose, const double& loop_density_vadose);

		/*!
		\brief Sets the simulation parameters for the amplification phase (number of new points)
		\param nb_new_deadend_points deadend points (see the method explanations in Paris et al. 2021)
		\param nb_interior_nodes deadend points (see the method explanations in Paris et al. 2021). NB: not yet functional,
		\ same effect as deadend points for now.
		*/
		void set_amplification_phreatic_params(const double& max_dist_loops_phreatic, const double& loop_density_phreatic);

		/*!
		\brief Sets water table level to differentiate vadose and phreatic zones.
		\param water_table_level Z coordinate of the water table
		\param inception_horizon_constraint_weight weight of the cost term to differentiate vadose and phreatic zones
		*/
		void set_water_table_level(const double& water_table_constraint_weight_vadose, const double& water_table_constraint_weight_phreatic);

		/*!
		\brief Sets the simulation so that it will not consider water table constraint
		*/
		void disable_water_table();

		///*!
		//\brief test function
		//\param nghb_count_t 
		//\param gamma_t 
		//*/
		//void test(double nghb_count_t, double gamma_t);


	private:
		std::string karstic_network_name;
		std::vector<Vector3> nodes_on_inception_surfaces;
		std::vector<std::vector<Vector3>> nodes_on_wt_surfaces;
		std::vector<Vector3> pt_sink;
		std::vector<Vector3> pt_spring;
		std::vector<Surface>* inception_horizons = nullptr;
		std::vector<Surface>* water_tables = nullptr;
		Surface* topo_surface_ = nullptr;
		double fracture_azimut;
		double fracture_dip;
		const std::vector<Vector3>* permeability_sphere_centers;
		double permeability_sphere_radius;
		double permeability_intensity;
		double permeability_constraint_weight;
		bool is_simulation_parametrized = false;
		Box* box;
		GeologicalParameters params;
		std::vector<KeyPoint> keypts;
		std::string save_directory;
	};
}
