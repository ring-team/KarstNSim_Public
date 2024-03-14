/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#pragma once

/**
@file run_code.h
@brief Contains main function used to run KarstNSim.
@author Augustin GOUY
**/

#include <iostream>
#include "KarstNSim/KarsticNetwork.h"
#include "KarstNSim/graph.h"
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
	struct ParamsSource
	{
		// keep those two ALWAYS at the beginning of the prompt file, because they need to be known to be used for the rest
		std::string karstic_network_name;
		std::string save_repertory;

		// Other parameters. Those can be defined in any order. If an optional parameter is not used, you can either add nothing after
		// its tag or don't add the tag altogether.
		// Some parameters are always mandatory, others sometimes mandatory, others completely optional.

		Box domain;
		int selected_seed = 1;
		int number_of_iterations = 1;
		bool vary_seed = true;
		bool allow_single_outlet_connection = true;
		bool use_sampling_points = false;
		std::vector<Vector3> sampling_points;
		bool use_density_property = false;
		bool create_grid = false;
		double poisson_radius = 0.01;
		bool use_karstification_potential = false;
		double karstification_potential_weight = 1;
		int k_pts = 10;
		int nghb_count = 25;
		bool use_max_nghb_radius = false;
		double nghb_radius = 0;
		double gamma = 2;
		bool vadose_cohesion = true;
		std::vector<Vector3> sinks;
		std::vector<Vector3> springs;
		bool use_waypoints = false;
		std::vector<Vector3> waypoints;
		int refine_surface_sampling = 0;
		Surface topo_surface;
		bool add_inception_surfaces = false;
		std::vector<Surface> inception_surfaces;
		double inception_surface_constraint_weight = 1;
		bool use_fracture_constraints = false;
		std::vector<double> fracture_families_orientations;
		std::vector<double> fracture_families_tolerance;
		double fracture_constraint_weight = 1;
		bool use_previous_networks = false;
		std::vector<Line> previous_networks;
		bool use_no_karst_spheres = false;
		std::vector<Vector3> sphere_centers;
		double fraction_karst_perm = 0.9;
		double max_inception_surface_distance = 50;
		double water_table_constraint_weight_vadose = 1;
		double water_table_constraint_weight_phreatic = 1;
		bool use_amplification_phase_vadose = false;
		double max_dist_loops_vadose;
		double loop_density_vadose;
		bool use_amplification_phase_phreatic = false;
		double max_dist_loops_phreatic;
		double loop_density_phreatic;
		bool create_nghb_graph = false;
		bool create_vset_sampling = true;
		double fraction_old_karst_perm = 0.5;
		std::vector<Surface> surf_wat_table;
		bool create_nghb_graph_property = false;
		bool multiply_costs = false;

		// parameters not included formally with a tag in the prompt file, but defined in another file (properties) :

		// extracted from domain box :values.size()
		std::vector<double> propdensity;
		std::vector<double> propikp;

		// extracted from springs pointset :
		std::vector<int> propspringsindex;

		// extracted from sinks pointset :
		std::vector<int> propsinksindex;
		std::vector<int> propsinksorder;

		// extracted from no-karst spheres pointset :
		std::vector<double> sphere_radius;
	};

	void run_simulation_full(ParamsSource parameters);

};
