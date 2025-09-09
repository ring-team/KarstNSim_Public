/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/parse_inputs.h"
#include <iostream>
#include <fstream>
#include <sstream>

ParseInputs::ParseInputs() {
	// Constructor implementation, if needed
}


bool ParseInputs::parse_boolean(const std::string& input) {
	bool flag;

	if (input == "true" || input == "1") {
		flag = true;
	}
	else if (input == "false" || input == "0") {
		flag = false;
	}
	else {
		std::cerr << "Invalid boolean input: " << input << std::endl;
		flag = false; // Default value for an invalid input
	}

	return flag;
}

std::vector<float> ParseInputs::read_distrib(const std::string& directory, const std::string& filename) {
	std::vector<float> values;
	std::string full_path = directory + "/" + filename;
	std::ifstream file(full_path);
	if (!file) {
		std::cerr << "Error opening file: " << full_path << std::endl;
		std::cerr << "Reason: " << strerror(errno) << std::endl;
		throw std::runtime_error("File not found: " + full_path);
	}

	std::string line;
	while (std::getline(file, line)) {
		std::stringstream ss(line);
		float value;
		while (ss >> value) {
			values.push_back(value);
		}
	}

	return values;
}

KarstNSim::ParamsSource ParseInputs::parse(const std::string& filename) {

	KarstNSim::ParamsSource params;

	// parameters not included formally with a tag in the prompt file, but defined in another file (properties) :

	// extracted from domain box :
	std::vector<float> propdensity;
	std::vector<float> propikp;
	// extracted from springs pointset :
	std::vector<int> propspringsindex;
	std::vector<float> propspringsradius; 
	std::vector<int> propspringssurfindex; 
	// extracted from sinks pointset :
	std::vector<int> propsinksindex;
	std::vector<int> propsinksorder;
	std::vector<float> propsinksradius; 
	// extracted from no-karst spheres pointset :
	std::vector<float> sphere_radius;
	// extracted from waypoints pointset :
	std::vector<float> waypoints_radius;
	std::vector<float> waypoints_impact_radius; 

	std::ifstream inputFile(filename);

	if (!inputFile) {
		std::cerr << "Error opening file: " << filename << std::endl;
		std::cerr << "Error: " << std::strerror(errno) << std::endl;
		throw std::runtime_error("Failed to open input file");
	}

	std::string line;
	while (std::getline(inputFile, line)) {

		// Skip lines starting with "//", "#", "%", or "::" (considered as comments)
		if (line.find("//") == 0 || line.find("#") == 0 || line.find("%") == 0 || line.find("::") == 0) {
			continue;
		}

		// Skip empty lines
		if (line.empty()) {
			continue;
		}

		// File name and repository

		std::istringstream iss(line);
		std::string paramType;
		iss >> paramType;
		if (paramType == "karstic_network_name:") {
			iss >> params.karstic_network_name;
		}
		else if (paramType == "main_repository:") {
			iss >> params.save_repository;
		}

		// General parameters

		else if (paramType == "domain:") {
			std::string box_path;
			iss >> box_path;
			KarstNSim::Box box;
			std::vector<std::vector<float>> prop_box;
			KarstNSim::load_box(box_path, params.save_repository, box, prop_box);
			params.domain = box;

			for (size_t i = 0; i < box.get_nu()*box.get_nv()*box.get_nw(); ++i) {
				if (prop_box[0].size() > 0) {
					params.propdensity.push_back(prop_box[i][0]);
					if (prop_box[0].size() > 1) {
						params.propikp.push_back(prop_box[i][1]);
					}
				}
			}
		}
		else if (paramType == "selected_seed:") {
			iss >> params.selected_seed;
		}
		else if (paramType == "number_of_iterations:") {
			iss >> params.number_of_iterations;
		}
		else if (paramType == "vary_seed:") {
			std::string flag;
			iss >> flag;
			params.vary_seed = parse_boolean(flag);
		}
		else if (paramType == "topo_surface:") {
			std::string surf_path;
			iss >> surf_path;
			std::vector<std::vector<float>> prop;
			KarstNSim::Surface surf;
			KarstNSim::load_surface(surf_path, params.save_repository, surf, prop);
			params.topo_surface = surf;
		}

		// Sampling

		else if (paramType == "use_sampling_points:") {
			std::string flag;
			iss >> flag;
			params.use_sampling_points = parse_boolean(flag);
		}
		else if (paramType == "sampling_points:") {
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<float>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repository,ptset, prop);
			params.sampling_points = ptset;
		}
		else if (paramType == "poisson_radius:") {
			iss >> params.poisson_radius;
		}
		else if (paramType == "use_density_property:") {
			std::string flag;
			iss >> flag;
			params.use_density_property = parse_boolean(flag);
		}
		else if (paramType == "k_pts:") {
			iss >> params.k_pts;
		}

		// Nearest neighbor graph

		else if (paramType == "nghb_count:") {
			iss >> params.nghb_count;
		}
		else if (paramType == "use_max_nghb_radius:") {
			std::string flag;
			iss >> flag;
			params.use_max_nghb_radius = parse_boolean(flag);
		}
		else if (paramType == "nghb_radius:") {
			iss >> params.nghb_radius;
		}

		// Sinks springs and waypoints

		else if (paramType == "use_springs_radius:") {
			std::string flag;
			iss >> flag;
			params.use_springs_radius = parse_boolean(flag);
		}
		else if (paramType == "use_sinks_radius:") {
			std::string flag;
			iss >> flag;
			params.use_sinks_radius = parse_boolean(flag);
		}
		else if (paramType == "use_waypoints_radius:") {
			std::string flag;
			iss >> flag;
			params.use_waypoints_radius = parse_boolean(flag);
		}
		else if (paramType == "sinks:") {
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<float>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repository, ptset, prop);
			params.sinks = ptset;
			for (int i = 0; i < params.sinks.size(); i++) {
				params.propsinksindex.push_back(prop[i][0]);
				params.propsinksorder.push_back(prop[i][1]);
				if (params.use_sinks_radius) {
					params.propsinksradius.push_back(prop[i][2]);
				}
			}
		}
		else if (paramType == "springs:") {
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<float>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repository, ptset, prop);
			params.springs = ptset;
			for (int i = 0; i < params.springs.size(); i++) {
				params.propspringsindex.push_back(prop[i][0]);
				params.propspringssurfindex.push_back(prop[i][1]);
				if (params.use_springs_radius) {
					params.propspringsradius.push_back(prop[i][2]);
				}
			}
		}
		else if (paramType == "use_waypoints:") {
		std::string flag;
		iss >> flag;
		params.use_waypoints = parse_boolean(flag);
		}
		else if (paramType == "waypoints:") {
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<float>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repository, ptset, prop);
			params.waypoints = ptset;
			for (int i = 0; i < params.waypoints.size(); i++) {
				params.waypoints_impact_radius.push_back(prop[i][0]);
				if (params.use_waypoints_radius) {
					params.waypoints_radius.push_back(prop[i][1]);
				}
			}
		}
		else if (paramType == "allow_single_outlet_connection:") {
		std::string flag;
		iss >> flag;
		params.allow_single_outlet_connection = parse_boolean(flag);
		}

		// ghost-rocks

		else if (paramType == "use_ghostrocks:") {
		std::string flag;
		iss >> flag;
		params.use_ghostrocks = parse_boolean(flag);
		}
		else if (paramType == "alteration_lines:") {
		std::string line_path;
		std::vector<std::vector<std::vector<float>>> prop;
		KarstNSim::Line line;
		iss >> line_path;
		KarstNSim::load_line(line_path, params.save_repository, line, prop);
		params.alteration_lines = line;
		}
		else if (paramType == "interpolate_lines:") {
		std::string flag;
		iss >> flag;
		params.interpolate_lines = parse_boolean(flag);
		}
		else if (paramType == "ghostrock_max_vertical_size:") {
		iss >> params.ghostrock_max_vertical_size;
		}
		else if (paramType == "use_max_depth_constraint:") {
		std::string flag;
		iss >> flag;
		params.use_max_depth_constraint = parse_boolean(flag);
		}
		else if (paramType == "ghost_rock_weight:") {
		iss >> params.ghost_rock_weight;
		}
		else if (paramType == "max_depth_horizon:") {
		std::string surf_path;
		std::vector<std::vector<float>> prop;
		KarstNSim::Surface surf;
		iss >> surf_path;
		KarstNSim::load_surface(surf_path, params.save_repository, surf, prop);
		params.max_depth_horizon = surf;
		}
		else if (paramType == "ghostrock_width:") {
		iss >> params.ghostrock_width;
		}

		// Inception surfaces

		else if (paramType == "refine_surface_sampling:") {
			iss >> params.refine_surface_sampling;
		}
		else if (paramType == "add_inception_surfaces:") {
		std::string flag;
		iss >> flag;
		params.add_inception_surfaces = parse_boolean(flag);
		}
		else if (paramType == "inception_surfaces:") {
			std::string surf_path;
			std::vector<std::vector<float>> prop;
			std::vector<KarstNSim::Surface> surf;
			while (iss >> surf_path) {
				KarstNSim::Surface surfi;
				KarstNSim::load_surface(surf_path, params.save_repository, surfi, prop);
				surf.push_back(surfi);
			}

			params.inception_surfaces = surf;
		}
		else if (paramType == "inception_surface_constraint_weight:") {
			iss >> params.inception_surface_constraint_weight;
		}
		else if (paramType == "max_inception_surface_distance:") {
		iss >> params.max_inception_surface_distance;
		}

		// Karstification potential

		else if (paramType == "use_karstification_potential:") {
		std::string flag;
		iss >> flag;
		params.use_karstification_potential = parse_boolean(flag);
		}
		else if (paramType == "karstification_potential_weight:") {
		iss >> params.karstification_potential_weight;
		}

		// Fracture constraints

		else if (paramType == "use_fracture_constraints:") {
		std::string flag;
		iss >> flag;
		params.use_fracture_constraints = parse_boolean(flag);
		}
		else if (paramType == "fracture_families_orientations:") {
		float value;
			while (iss >> value) {
				params.fracture_families_orientations.emplace_back(value);
			}
		}
		else if (paramType == "fracture_families_tolerance:") {
		float value;
			while (iss >> value) {
				params.fracture_families_tolerance.emplace_back(value);
			}
		}
		else if (paramType == "fracture_constraint_weight:") {
		iss >> params.fracture_constraint_weight;
		}

		// Previous networks (polyphasic karstification or simulation of dimensions exclusively)

		else if (paramType == "use_previous_networks:") {
		std::string flag;
		iss >> flag;
		params.use_previous_networks = parse_boolean(flag);
		}
		else if (paramType == "previous_networks:") {
			std::string line_path;
			std::vector<std::vector<std::vector<float>>> prop;
			std::vector<KarstNSim::Line> lines;
			while (iss >> line_path) {
				KarstNSim::Line line;
				KarstNSim::load_line(line_path, params.save_repository, line, prop);
				lines.push_back(line);
			}
			params.previous_networks = lines;
		}
		else if (paramType == "fraction_old_karst_perm:") {
		iss >> params.fraction_old_karst_perm;
		}
		else if (paramType == "sections_simulation_only:") {
		std::string flag;
		iss >> flag;
		params.sections_simulation_only = parse_boolean(flag);
		}

		// No-karst spheres

		else if (paramType == "use_no_karst_spheres:") {
		std::string flag;
		iss >> flag;
		params.use_no_karst_spheres = parse_boolean(flag);
		}
		else if (paramType == "sphere_centers:") {
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<float>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repository, ptset, prop);
			params.sphere_centers = ptset;
		}

		// Water tables

		else if (paramType == "water_table_constraint_weight_vadose:") {
		iss >> params.water_table_constraint_weight_vadose;
		}
		else if (paramType == "water_table_constraint_weight_phreatic:") {
		iss >> params.water_table_constraint_weight_phreatic;
		}
		else if (paramType == "surf_wat_table:") {
		std::string surf_path;
		std::vector<std::vector<float>> prop;
		std::vector<KarstNSim::Surface> surf;
		while (iss >> surf_path) {
			KarstNSim::Surface surfi;
			KarstNSim::load_surface(surf_path, params.save_repository, surfi, prop);
			surf.push_back(surfi);
		}
		params.surf_wat_table = surf;
		}

		// Deadend points amplification :

		else if (paramType == "use_deadend_points:") {
		std::string flag;
		iss >> flag;
		params.use_deadend_points = parse_boolean(flag);
		}
		else if (paramType == "nb_deadend_points:") {
		iss >> params.nb_deadend_points;
		}
		else if (paramType == "max_distance_of_deadend_pts:") {
		iss >> params.max_distance_of_deadend_pts;
		}

		// Cycle amplification :

		else if (paramType == "use_cycle_amplification:") {
		std::string flag;
		iss >> flag;
		params.use_amplification = parse_boolean(flag);
		}
		else if (paramType == "max_distance_amplification:") {
		iss >> params.max_distance_amplification;
		}
		else if (paramType == "min_distance_amplification:") {
		iss >> params.min_distance_amplification;
		}
		else if (paramType == "nb_cycles:") {
		iss >> params.nb_cycles;
		}

		//Noise Amplification parameters :

		else if (paramType == "use_noise:") {
		std::string flag;
		iss >> flag;
		params.use_noise = parse_boolean(flag);
		}
		else if (paramType == "use_noise_on_all:") {
		std::string flag;
		iss >> flag;
		params.use_noise_on_all = parse_boolean(flag);
		}
		else if (paramType == "noise_frequency:") {
		iss >> params.noise_frequency;
		}
		else if (paramType == "noise_octaves:") {
		iss >> params.noise_octaves;
		}
		else if (paramType == "noise_weight:") {
		iss >> params.noise_weight;
		}

		// Create sections :

		else if (paramType == "simulate_sections:") {
		std::string flag;
		iss >> flag;
		params.simulate_sections = parse_boolean(flag);
		if (params.simulate_sections) {
			params.geostat_params.is_used = true;
		}
		}
		else if (paramType == "simulation_distribution:") {
		std::string path;
		iss >> path;
		std::vector<float> sim_distrib_discrete = read_distrib(params.save_repository, path);
		params.geostat_params.simulation_distribution = sim_distrib_discrete;
		}
		else if (paramType == "global_vario_range:") {
		iss >> params.geostat_params.global_vario_range;
		}
		else if (paramType == "global_range_of_neighborhood:") {
		iss >> params.geostat_params.global_range_of_neighborhood;
		}
		else if (paramType == "global_vario_sill:") {
		iss >> params.geostat_params.global_vario_sill;
		}
		else if (paramType == "global_vario_nugget:") {
		iss >> params.geostat_params.global_vario_nugget;
		}
		else if (paramType == "global_vario_model:") {
		iss >> params.geostat_params.global_vario_model;
		}
		else if (paramType == "interbranch_vario_range:") {
		iss >> params.geostat_params.interbranch_vario_range;
		}
		else if (paramType == "interbranch_range_of_neighborhood:") {
		iss >> params.geostat_params.interbranch_range_of_neighborhood;
		}
		else if (paramType == "interbranch_vario_sill:") {
		iss >> params.geostat_params.interbranch_vario_sill;
		}
		else if (paramType == "interbranch_vario_nugget:") {
		iss >> params.geostat_params.interbranch_vario_nugget;
		}
		else if (paramType == "interbranch_vario_model:") {
		iss >> params.geostat_params.interbranch_vario_model;
		}
		else if (paramType == "intrabranch_vario_range:") {
		iss >> params.geostat_params.intrabranch_vario_range;
		}
		else if (paramType == "intrabranch_range_of_neighborhood:") {
		iss >> params.geostat_params.intrabranch_range_of_neighborhood;
		}
		else if (paramType == "intrabranch_vario_sill:") {
		iss >> params.geostat_params.intrabranch_vario_sill;
		}
		else if (paramType == "intrabranch_vario_nugget:") {
		iss >> params.geostat_params.intrabranch_vario_nugget;
		}
		else if (paramType == "intrabranch_vario_model:") {
		iss >> params.geostat_params.intrabranch_vario_model;
		}
		else if (paramType == "number_max_of_neighborhood_points:") {
		iss >> params.geostat_params.number_max_of_neighborhood_points;
		}
		else if (paramType == "nb_points_interbranch:") {
		iss >> params.geostat_params.nb_points_interbranch;
		}
		else if (paramType == "proportion_interbranch:") {
		iss >> params.geostat_params.proportion_interbranch;
		}

		// Save objects

		else if (paramType == "create_vset_sampling:") {
		std::string flag;
		iss >> flag;
		params.create_vset_sampling = parse_boolean(flag);
		}
		else if (paramType == "create_nghb_graph:") {
		std::string flag;
		iss >> flag;
		params.create_nghb_graph = parse_boolean(flag);
		}
		else if (paramType == "create_nghb_graph_property:") {
		std::string flag;
		iss >> flag;
		params.create_nghb_graph_property = parse_boolean(flag);
		}
		else if (paramType == "create_grid:") {
		std::string flag;
		iss >> flag;
		params.create_grid = parse_boolean(flag);
		}

		// Other parameters

		else if (paramType == "fraction_karst_perm:") {
		iss >> params.fraction_karst_perm;
		}
		else if (paramType == "multiply_costs:") {
		std::string flag;
		iss >> flag;
		params.multiply_costs = parse_boolean(flag);
		}
		else if (paramType == "gamma:") {
		iss >> params.gamma;
		}
		else if (paramType == "vadose_cohesion:") {
		std::string flag;
		iss >> flag;
		params.vadose_cohesion = parse_boolean(flag);
		}
	}
	inputFile.close();
	return params;
}