/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#pragma once

/**
* @mainpage KarstNSim
*
* @section intro_sec Introduction
*
* KarstNSim is a C++ code for graph-based and geologically-driven simulation of
* 3D karst networks, mainly developed during Augustin Gouy's PhD thesis (2022-2025).
*
* It is highly advised to read the main publication first:
*
* Gouy, A., Collon, P., Bailly-Comte, V., Galin, E., Antoine, C., Thebault, B., & Landrein, P. (2024). KarstNSim:
*  A graph-based method for 3D geologically-driven simulation of karst networks. Journal of Hydrology.
*  https://doi.org/10.1016/j.jhydrol.2024.130878
*
* as well as the PhD thesis (2025).
*
* @section features_sec Usage
* 
* You can find a description of the instruction file parameters in the KarstNSim::ParamsSource struct.
* 
* The most effective way to read the documentation is to go to Files and read the information file by file. Indeed, while most files contain a same name class, some
* files are utility files with helper functions, and some others contain multiple classes, whereas each file has a well defined purpose in the project, making reading
* file by file easier than class by class.
*
* The main order of use of each file is the following (apart from main.cpp and some main classes, header files are cited here for convenience):
*	* main.cpp: main of the project. Creates an instance of the ParseInputs class using its constructor, and uses it to call 'run_simulation_full'.
*	* parse_inputs.h: contains a declaration of the ParseInputs class, which reads the instructions.txt file as well as the supplementary ASCII files and transforms these parameters into C++ objects.
*		+ read_files.h: To easily read other ASCII files (pointsets, TINs, graphs, box), helper functions of this file are used.
*	* run_code.h: contains the main function ('run_simulation_full') used to run KarstNSim, and acts as an interface between the main and the KarsticNetwork class.
*	* karstic_network.h: contains a declaration of the class KarsticNetwork, which stores input parameters and allows to easily set up the karstic network simulation. It then creates an instance of
*     GraphOperations to create the cost graph, and then of KarsticSkeleton to create the karst skeleton and if needed, amplify it. Finally, if required, it simulates sections on the karst skeleton nodes, and saves the outputs.
*		+ geology.h: in addition to the attributes of KarsticNetwork, other more generic parameters are stored in classes representing geological parameters, which are defined in this file.
*		+ ghost_rocks.h: if used, ghost-rocks are taken into account in the KarsticNetwork class during parameter setup through helper functions present in this utility file.
*		+ surface_sampling.h: declares the surface_sampling class used to sample points onto triangulated surfaces vertices, and refine sampling if needed, both are done in KarsticNetwork as part of the parameter setup.
*		+ graph.h: declares various classes representing the cost graph (CostGraph), operations on this cost graph (GraphOperations, inheriting from CostGraph), and the skeleton graph (KarsticSkeleton). The two last ones are
*		 instanced and used inside KarsticNetwork once parameter setup is done. These 3 classes are implemented in 3 separate source files with eponymous file names.
*			- GraphOperations: contains methods allowing to create the sampling cloud, generate the nearest neighbor graph onto it, compute costs for each edge and prune the graph to generate the karstic network.
*			- CostGraph: it is the class representing the nearest neighbor cost graph, notably through the attribute 'adj'. It also contains methods used to update and prune this cost graph (and, notably, Dijkstra algorithms).
*			- KarsticSkeleton: intervenes after GraphOperations and CostGraph. It is the class that represents the simulated karstic skeleton. It contains several methods to update it and amplify it, as well as other operations.
*		+ geostats.h: contains methods used in KarstNSim for equivalent section simulation through 1D-Curvilinear branchwise SGS, if required. It is a 1:1 translation of the code of Frantz et al. (2021), but devoid of any
*		 dependencies (as opposed to the initial version). These methods are called in KarsticNetwork, once the karstic skeleton is simulated.
*		+ write_files.h: also called in KarsticNetwork, it contains helper functions to write down outputs in ASCII format and save them once everything else is done.
*
* Several other utility files are used throughout the code:
*	* vec.h: declares and defines utility classes and functions for simple vector objects (in 2D and 3D), as well as a structure handling KD-tree partitioning and search, through the class PointCloud.
*		+ nanoflann.hpp: header-only library for building KD-Trees, called by the PointCloud class methods.
*	* basics.h: declares and defines various classes representing geometrical objects used for the (hydro)geological model and karst network representation: Segment, Line, Plane, Sphere, Triangle, Surface, Box...
*	* randomgenerator.h: declares utility functions allowing generation of seeds, random numbers and random processes.
*	* simplex_noise.h: declares a SimplexNoise class allowing to generate simplex noise. It was implemented in the code to add noise to the whole cost graph, either for the full simulation, or just for cycle generation.
*	In the end, it was not used nor presented in the 2024 article and the 2025 thesis due to its limited value in the algorithm.
*
* @section installation_sec Installation
*
* To use this library, follow the instructions presented in
* the associated github repository at https://github.com/ring-team/KarstNSim_Public.
*
* @section licensing_sec Licensing
*
* This project is licensed under the MIT License.
*/

/**
@file run_code.h
@brief Contains main function used to run KarstNSim.
@author Augustin GOUY
**/

#include <iostream>
#include "KarstNSim/karstic_network.h"
#include <string>
#include "KarstNSim/basics.h"
#include "KarstNSim/geology.h"
#include "KarstNSim/surface_sampling.h"
#include <cmath>
#include "KarstNSim/vec.h"
#include <fstream>
#include <map>
#include <algorithm>
#include <string>
#include <random>
#include "KarstNSim/randomgenerator.h"

namespace KarstNSim {
	/*!
	\struct PramsSource
	\brief Structure defining all simulation parameters, in C++ and custom object format
	*/
	struct ParamsSource
	{
		// Names :

		std::string karstic_network_name; //!< Name of the simulation, will be attached to all outputs
		std::string save_repository; //!< Name of the folder in which results should be saved

		// General parameters :

		Box domain; //!< Background grid for the simulation
		int selected_seed = 1; //!< Seed used for the simulation
		int number_of_iterations = 1; //!< Number of iterations used for the simulation
		bool vary_seed = true; //!< Flag defining if seeds should be changed for each new simulation
		Surface topo_surface; //!< Topographic surface

		// Use already sampled points... :

		bool use_sampling_points = false; //!< Flag defining if already generated sampling cloud should be used instead of generating it
		std::vector<Vector3> sampling_points; //!< If use_sampling_points is true, sampling cloud used for the simulation

		// ...Or sample new points (constant density OR density property scanned from grid/box) :

		float poisson_radius = 0.01f; //!< Poisson radius used uniformly in the domain if no density property is defined (use_density_property = false)
		bool use_density_property = false; //!< Flag defining if a density property is used and defined on the background grid
		int k_pts = 10; //!< Value of the parameter k, used in Dwork's algorithm as the number of points to be generated around a random point at each step

		// Re-use previously generated networks to constrain simulation (polyphasic karstification) :

		bool use_previous_networks = false; //!< Flag defining if previous networks have been simulated and should be incorporated to the simulation (polyphasic simulation)
		std::vector<Line> previous_networks; //!< If use_previous_networks is true, list of previous networks to be incorporated in the simulation
		float fraction_old_karst_perm = 0.5f; //!< Corresponds to Ppoly, the polyphasic cost reduction factor: any edge used by a previous network has its cost reduced by this factor
		bool sections_simulation_only = false; //!< Flag defining whether only equivalent sections should be simulated, as opposed to all steps of KarstNSim. If set to true, previous networks should be defined

		// N nearest-neighbor graph parameters :

		int nghb_count = 25; //!< Number of nearest neighbors
		bool use_max_nghb_radius = false; //!< Flag defining if a maximum radius constraint should be used for the nearest neighbor graph
		float nghb_radius = 0.0f; //!< Maximum radius for nearest neighbor search

		// Ghost-rocks :

		bool use_ghostrocks = false; //!< Flag to include ghost-rocks to the simulation.
		Line alteration_lines; //!< Lines defining surface alteration lines for ghost-rock alterations.
		bool interpolate_lines = false; //!< Flag to enable interpolation between alteration lines (NOT DEVELOPED YET).
		float ghostrock_max_vertical_size; //!< Maximum vertical size for ghost-rock areas.
		bool use_max_depth_constraint = false; //!< Flag to apply a max depth constraint to ghost-rock areas.
		float ghost_rock_weight; //!< Weight for ghost-rock constraints in simulation (will be applied to the IKP cost).
		Surface max_depth_horizon; //!< Surface defining the maximum depth for ghost-rocks.
		float ghostrock_width; //!< Maximum width of the ghost-rock areas.

		// Inlets, outlets and waypoints :

		std::vector<Vector3> sinks; //!< List of sinks.
		std::vector<Vector3> springs; //!< List of springs.
		bool allow_single_outlet_connection = true; //!< Flag to impose inlets to have a single spring connection. This will necessitate the use of the "closest" spring algorithm.
		bool use_waypoints = false; //!< Flag to use waypoints as constraints in simulation.
		std::vector<Vector3> waypoints; //!< List of waypoints.
		bool use_springs_radius; //!< Flag to use radii for spring for equivalent section simulation constraints.
		bool use_sinks_radius; //!< Flag to use radii for sink for equivalent section simulation constraints.
		bool use_waypoints_radius; //!< Flag to use radii for waypoint for equivalent section simulation constraints.
		float waypoints_weight = 0.1f; //!< Weight for waypoint constraints.

		// No-karst spheres :

		bool use_no_karst_spheres = false; //!< Flag to define no-karst spheres in the domain.
		std::vector<Vector3> sphere_centers; //!< Centers of no-karst spheres.

		// Inception surfaces :

		bool add_inception_surfaces = false; //!< Flag to add inception surfaces in simulation.
		int refine_surface_sampling = 0; //!< Refinement level for surface sampling.
		std::vector<Surface> inception_surfaces; //!< List of inception surfaces.
		float inception_surface_constraint_weight = 1.0f; //!< Weight for inception surface constraints.
		float max_inception_surface_distance = 50.0f; //!< Maximum distance for inception surface constraints.

		// Karstification potential (property scanned from grid/box) :

		bool use_karstification_potential = false; //!< Flag to use karstification potential in simulation.
		float karstification_potential_weight = 1.0f; //!< Weight for karstification potential in cost function.

		// Fractures (orientations, as angles to the north, and angular tolerance for each family) :

		bool use_fracture_constraints = false; //!< Flag to use fracture constraints in simulation.
		std::vector<float> fracture_families_orientations; //!< Orientations of fracture families.
		std::vector<float> fracture_families_tolerance; //!< Angular tolerances for each fracture family.
		float fracture_constraint_weight = 1.0f; //!< Weight for fracture constraints in the cost function.

		// Water tables (one for each spring) and vadose and phreatic trends (mandatory) :

		std::vector<Surface> surf_wat_table; //!< Surfaces representing water tables.
		float water_table_constraint_weight_vadose = 1.0f; //!< Weight for vadose zone water table constraints.
		float water_table_constraint_weight_phreatic = 1.0f; //!< Weight for phreatic zone water table constraints.

		// Other general cost graph parameters (gamma-graph parameter, cost reduction (cohesion) factor, and option to multiply costs instead of summing them) :

		float gamma = 2.0f; //!< Gamma parameter for the network pruning based on graph-graph rule.
		float fraction_karst_perm = 0.9f; //!< Cost reduction (cohesion) factor Pred for karsts.
		bool vadose_cohesion = true; //!< Flag to enable cohesion in vadose zone (if false, cohesion is only applied in phreatic zone).
		bool multiply_costs = false; //!< Flag to multiply costs instead of summing them.

		// Deadend points amplification :

		bool use_deadend_points = false; //!< Flag to include dead-end points in simulation.
		int nb_deadend_points; //!< Number of dead-end points.
		float max_distance_of_deadend_pts; //!< Maximum distance for dead-end points.

		// Cycle amplification :

		bool use_amplification = true; //!< Flag to use amplification in simulation.
		float max_distance_amplification; //!< Maximum distance between two random nodes used to create a cycle for amplification.
		float min_distance_amplification; //!< Minimum distance between two random nodes used to create a cycle for amplification.
		int nb_cycles; //!< Number of amplification cycles.

		//Noise Amplification parameters :

		bool use_noise = true; //!< Flag to include noise in amplification only.
		bool use_noise_on_all = false; //!< Flag to apply noise for both karst network simulation AND amplification.
		int noise_frequency = 1; //!< Frequency of noise application.
		int noise_octaves = 1; //!< Number of octaves for noise.
		float noise_weight = 1; //!< Weight for noise contribution.

		// Create sections :

		bool simulate_sections = false; //!< Flag to simulate sections for the karst network.
		GeostatParams geostat_params; //!< Parameters for geostatistical simulation.

		// Save parameters (if you want to save sampling points, cost graph and cost property) :

		bool create_vset_sampling = true; //!< Flag to save sampling points.
		bool create_nghb_graph = false; //!< Flag to save nearest neighbor graph (quite heavy object!)
		bool create_nghb_graph_property = false; //!< Flag to save nearest neighbor graph property (very heavy object!!)
		bool create_grid = false; //!< Flag to save grid data.

		// parameters not included formally with a tag in the prompt file, but defined in another file (properties) :

		std::vector<float> propdensity; //!< Property density values extracted from the domain box.
		std::vector<float> propikp; //!< Property intrinsic karstification potential values from the domain box.

		// extracted from springs pointset :
		std::vector<int> propspringsindex; //!< Indices of spring properties (column index in connectivity matrix).
		std::vector<float> propspringsradius; //!< Radii of springs. -> for conduit size
		std::vector<int> propspringssurfindex; //!< Surface indices for springs (index of water table associated with each spring).

		// extracted from sinks pointset :
		std::vector<int> propsinksindex; //!< Indices of sink properties (row index in connectivity matrix).
		std::vector<int> propsinksorder; //!< Order of sinks (see 2024 paper for signification of sink order)
		std::vector<float> propsinksradius; //!< Radii of sinks. -> for conduit size

		// extracted from waypoints pointset :
		std::vector<float> waypoints_radius; //!< Radii of waypoints. -> for conduit size
		std::vector<float> waypoints_impact_radius; //!< Impact radii of waypoints. -> for influence area

		// extracted from no-karst spheres pointset :
		std::vector<float> sphere_radius; //!< Radii of no-karst spheres.

	};


	/*!
	 * @brief Runs the full KarstNSim simulation with the given parameters.
	 * @param parameters Simulation parameters as defined in ParamsSource.
	 */
	void run_simulation_full(ParamsSource parameters);

	/*!
	 * @brief Runs a partial simulation focusing only on equivalent section property calculations on an already simulated karst network.
	 * @param parameters Simulation parameters as defined in ParamsSource.
	 */
	void run_simulation_properties(ParamsSource parameters);
};