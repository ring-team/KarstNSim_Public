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


bool ParseInputs::parseBoolean(const std::string& input) {
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

KarstNSim::ParamsSource ParseInputs::parse(const std::string& filename) {

	KarstNSim::ParamsSource params;

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

	std::ifstream inputFile(filename);

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

		std::istringstream iss(line);
		std::string paramType;
		iss >> paramType;
		if (paramType == "karstic_network_name:") {
			iss >> params.karstic_network_name;
		}
		else if (paramType == "main_repository:") {
			iss >> params.save_repertory;
		}

		else if (paramType == "use_density_property:") {
			std::string flag;
			iss >>flag;
			params.use_density_property = parseBoolean(flag);
		}

		else if (paramType == "use_karstification_potential:") {
			std::string flag;
			iss >> flag;
			params.use_karstification_potential = parseBoolean(flag);
		}

		else if (paramType == "domain:") {
			// box object requires some steps for parsing
			std::string box_path;
			iss >> box_path;
			KarstNSim::Box box;
			std::vector<std::vector<double>> prop_box;
			KarstNSim::load_box(box_path, params.save_repertory, box, prop_box);
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
			params.vary_seed = parseBoolean(flag);
		}
		else if (paramType == "allow_single_outlet_connection:") {
			std::string flag;
			iss >> flag;
			params.allow_single_outlet_connection = parseBoolean(flag);
		}
		else if (paramType == "vadose_cohesion:") {
			std::string flag;
			iss >> flag;
			params.vadose_cohesion = parseBoolean(flag);
		}
		else if (paramType == "use_sampling_points:") {
			std::string flag;
			iss >> flag;
			params.use_sampling_points = parseBoolean(flag);
		}
		else if (paramType == "sampling_points:") {
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<double>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repertory,ptset, prop);
			params.sampling_points = ptset;
		}

		else if (paramType == "poisson_radius:") {
			iss >> params.poisson_radius;
		}

		else if (paramType == "karstification_potential_weight:") {
			iss >> params.karstification_potential_weight;
		}
		else if (paramType == "k_pts:") {
			iss >> params.k_pts;
		}
		else if (paramType == "nghb_count:") {
			iss >> params.nghb_count;
		}
		else if (paramType == "use_max_nghb_radius:") {
			std::string flag;
			iss >> flag;
			params.use_max_nghb_radius = parseBoolean(flag);
		}
		else if (paramType == "nghb_radius:") {
			iss >> params.nghb_radius;
		}
		else if (paramType == "gamma:") {
			iss >> params.gamma;
		}
		else if (paramType == "sinks:") {
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<double>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repertory, ptset, prop);
			params.sinks = ptset;
			for (int i = 0; i < params.sinks.size(); i++) {
				params.propsinksindex.push_back(prop[i][0]);
				params.propsinksorder.push_back(prop[i][1]);
			}
		}
		else if (paramType == "springs:") {
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<double>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repertory, ptset, prop);
			params.springs = ptset;
			for (int i = 0; i < params.springs.size(); i++) {
				params.propspringsindex.push_back(prop[i][0]);
			}
		}
		else if (paramType == "use_waypoints:") {
		std::string flag;
		iss >> flag;
		params.use_waypoints = parseBoolean(flag);
		}
		else if (paramType == "waypoints:") {
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<double>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repertory, ptset, prop);
			params.waypoints = ptset;
		}
		else if (paramType == "refine_surface_sampling:") {
			iss >> params.refine_surface_sampling;
		}
		else if (paramType == "topo_surface:") {
			std::string surf_path;
			iss >> surf_path;
			std::vector<std::vector<double>> prop;
			KarstNSim::Surface surf;
			KarstNSim::load_surface(surf_path, params.save_repertory, surf, prop);
			params.topo_surface = surf;
		}
		else if (paramType == "add_inception_surfaces:") {
		std::string flag;
		iss >> flag;
		params.add_inception_surfaces = parseBoolean(flag);
		}
		else if (paramType == "inception_surfaces:") {
			std::string surf_path;
			std::vector<std::vector<double>> prop;
			std::vector<KarstNSim::Surface> surf;
			while (iss >> surf_path) {
				KarstNSim::Surface surfi;
				KarstNSim::load_surface(surf_path, params.save_repertory, surfi, prop);
				surf.push_back(surfi);
			}

			params.inception_surfaces = surf;
		}
		else if (paramType == "inception_surface_constraint_weight:") {
			iss >> params.inception_surface_constraint_weight;
		}
		else if (paramType == "use_fracture_constraints:") {
		std::string flag;
		iss >> flag;
		params.use_fracture_constraints = parseBoolean(flag);
		}
		else if (paramType == "fracture_families_orientations:") {
		// Assuming double values are space-separated in the line
		double value;
			while (iss >> value) {
				params.fracture_families_orientations.emplace_back(value);
			}
		}
		else if (paramType == "fracture_families_tolerance:") {
		// Assuming double values are space-separated in the line
		double value;
			while (iss >> value) {
				params.fracture_families_tolerance.emplace_back(value);
			}
		}
		else if (paramType == "fracture_constraint_weight:") {
		iss >> params.fracture_constraint_weight;
		}
		else if (paramType == "use_previous_networks:") {
		std::string flag;
		iss >> flag;
		params.use_previous_networks = parseBoolean(flag);
		}
		else if (paramType == "previous_networks:") {
			std::string line_path;
			std::vector<std::vector<std::vector<double>>> prop;
			std::vector<KarstNSim::Line> lines;
			while (iss >> line_path) {
				KarstNSim::Line line;
				KarstNSim::load_line(line_path, params.save_repertory, line, prop);
				lines.push_back(line);
			}

			params.previous_networks = lines;
		}
		else if (paramType == "use_no_karst_spheres:") {
		std::string flag;
		iss >> flag;
		params.use_no_karst_spheres = parseBoolean(flag);
		}
		else if (paramType == "sphere_centers:") {
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<double>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repertory, ptset, prop);
			params.sphere_centers = ptset;
		}
		else if (paramType == "fraction_karst_perm:") {
		iss >> params.fraction_karst_perm;
		}
		else if (paramType == "max_inception_surface_distance:") {
		iss >> params.max_inception_surface_distance;
		}
		else if (paramType == "water_table_constraint_weight_vadose:") {
		iss >> params.water_table_constraint_weight_vadose;
		}
		else if (paramType == "water_table_constraint_weight_phreatic:") {
		iss >> params.water_table_constraint_weight_phreatic;
		}
		else if (paramType == "use_amplification_phase_vadose:") {
		std::string flag;
		iss >> flag;
		params.use_amplification_phase_vadose = parseBoolean(flag);
		}
		else if (paramType == "max_dist_loops_vadose:") {
		iss >> params.max_dist_loops_vadose;
		}
		else if (paramType == "loop_density_vadose:") {
		iss >> params.loop_density_vadose;
		}
		else if (paramType == "use_amplification_phase_phreatic:") {
		std::string flag;
		iss >> flag;
		params.use_amplification_phase_phreatic = parseBoolean(flag);
		}
		else if (paramType == "max_dist_loops_phreatic:") {
		iss >> params.max_dist_loops_phreatic;
		}
		else if (paramType == "loop_density_phreatic:") {
		iss >> params.loop_density_phreatic;
		}
		else if (paramType == "create_nghb_graph:") {
		std::string flag;
		iss >> flag;
		params.create_nghb_graph = parseBoolean(flag);
		}
		else if (paramType == "create_vset_sampling:") {
		std::string flag;
		iss >> flag;
		params.create_vset_sampling = parseBoolean(flag);
		}
		else if (paramType == "fraction_old_karst_perm:") {
		iss >> params.fraction_old_karst_perm;
		}
		else if (paramType == "surf_wat_table:") {
			std::string surf_path;
			std::vector<std::vector<double>> prop;
			std::vector<KarstNSim::Surface> surf;
			while (iss >> surf_path) {
				KarstNSim::Surface surfi;
				KarstNSim::load_surface(surf_path, params.save_repertory, surfi, prop);
				surf.push_back(surfi);
			}

			params.surf_wat_table = surf;
		}
		else if (paramType == "create_nghb_graph_property:") {
		std::string flag;
		iss >> flag;
		params.create_nghb_graph_property = parseBoolean(flag);
		}
		else if (paramType == "multiply_costs:") {
		std::string flag;
		iss >> flag;
		params.multiply_costs = parseBoolean(flag);
		}
	}
	inputFile.close();
	return params;
}