/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/run_code.h"

namespace KarstNSim {
	void run_simulation_full( ParamsSource parameters) {

		// Uncomment to save inputs for special use

		//// save box voxet
		//Box saved_box = parameters.domain;
		//std::vector<std::string> property_names = { "density","karstif_potential" };
		//std::vector<std::vector<double>> properties;
		//for (size_t i = 0; i < parameters.propdensity.size(); ++i) {
		//	properties.push_back({ parameters.propdensity[i], parameters.propikp[i] });
		//}
		//std::string full_name = parameters.karstic_network_name + "_box.txt";
		//std::string full_dir_name = parameters.save_repertory + "/outputs";
		//save_box(full_name, full_dir_name, saved_box, property_names,properties);

		//// save topography
		//Surface topo = parameters.topo_surface;
		//full_name = parameters.karstic_network_name + "_topo_surf.txt";
		//full_dir_name = parameters.save_repertory + "/outputs";
		//save_surface(full_name, full_dir_name, topo);

		//// save inception surfaces
		//std::vector<Surface> inception_surfaces = parameters.inception_surfaces;
		//full_name = parameters.karstic_network_name + "_inception_surf";
		//full_dir_name = parameters.save_repertory + "/outputs";
		//for (int surfi = 0; surfi < inception_surfaces.size(); surfi++) {
		//	std::string full_name_i = full_name + std::to_string(surfi) + ".txt";
		//	save_surface(full_name_i, full_dir_name, inception_surfaces[surfi]);
		//}

		//// save water tables
		//std::vector<Surface> wtsurf = parameters.surf_wat_table;
		//full_name = parameters.karstic_network_name + "_watertable_surf";
		//full_dir_name = parameters.save_repertory + "/outputs";
		//for (int surfi = 0; surfi < wtsurf.size(); surfi++) {
		//	std::string full_name_i = full_name + std::to_string(surfi) + ".txt";
		//	save_surface(full_name_i, full_dir_name, wtsurf[surfi]);
		//}

		//// save sinks
		//std::vector<Vector3> sinks = parameters.sinks;
		//std::vector<std::string> property_names_sinks = { "index","order" };
		//std::vector<std::vector<double>> properties_sinks;
		//for (size_t i = 0; i < parameters.propsinksindex.size(); ++i) {
		//	properties_sinks.push_back({ double(parameters.propsinksindex[i]), double(parameters.propsinksindex[i])});
		//}
		//full_name = parameters.karstic_network_name + "_sinks.txt";
		//full_dir_name = parameters.save_repertory + "/outputs";
		//save_pointset(full_name, full_dir_name, sinks, property_names_sinks, properties_sinks);

		//// save springs
		//std::vector<Vector3> springs = parameters.springs;
		//std::vector<std::string> property_names_springs(1, "index");
		//std::vector<std::vector<double>> properties_springs;
		//for (size_t i = 0; i < parameters.propspringsindex.size(); ++i) {
		//	properties_springs.push_back({ double(parameters.propspringsindex[i]) });
		//}
		//full_name = parameters.karstic_network_name + "_springs.txt";
		//full_dir_name = parameters.save_repertory + "/outputs";
		//save_pointset(full_name, full_dir_name, springs, property_names_springs, properties_springs);

		//// save waypoints
		//std::vector<Vector3> waypoints = parameters.waypoints;
		//full_name = parameters.karstic_network_name + "_waypoints.txt";
		//full_dir_name = parameters.save_repertory + "/outputs";
		//save_pointset(full_name, full_dir_name, waypoints);

		for (int i = 0; i < parameters.number_of_iterations; i++) {

			// Initialize with selected seed
			if (parameters.number_of_iterations != 1 && parameters.vary_seed) {
				unsigned int uintValue = static_cast<unsigned int>(parameters.selected_seed + i);
				std::vector<std::uint32_t> seed = {uintValue};
				initializeRng(seed);
			}
			else {
				unsigned int uintValue = static_cast<unsigned int>(parameters.selected_seed);
				std::vector<std::uint32_t> seed = { uintValue };
				initializeRng(seed);
			}

			// Creation of the object that will contain all parameters for the karstic simulation
			GeologicalParameters params;

			// Creation of the vector that will contain all points crossed by the karst
			std::vector<KeyPoint> keypts;

			std::string sim_name_iter = parameters.karstic_network_name + "_" +std::to_string(i);

			KarsticNetwork karst(sim_name_iter, &parameters.domain, params, keypts, &parameters.surf_wat_table);
			karst.set_save_directory(parameters.save_repertory);
			karst.set_sinks(&parameters.sinks, parameters.propsinksindex, parameters.propsinksorder);
			karst.set_springs(&parameters.springs, parameters.propspringsindex,parameters.allow_single_outlet_connection);
			karst.set_water_table_order(&parameters.surf_wat_table);

			if (parameters.use_waypoints) {
				karst.set_waypoints(&parameters.waypoints);
			}

			if (parameters.use_previous_networks) {
				karst.set_previous_networks(&parameters.previous_networks);
			}

			if (parameters.use_karstification_potential) {
				karst.set_karstification_potential_parameters(parameters.karstification_potential_weight);
			}

			if (parameters.create_grid) {
				karst.save_painted_box(parameters.propdensity, parameters.propikp);
			}

			karst.set_wt_surfaces_sampling(sim_name_iter, &parameters.surf_wat_table, parameters.refine_surface_sampling);

			if ((parameters.add_inception_surfaces)) {
				if (!parameters.use_sampling_points) {
					karst.set_inception_surfaces_sampling( sim_name_iter, &parameters.inception_surfaces, parameters.refine_surface_sampling, parameters.create_vset_sampling);
				}
				karst.set_inception_horizons_parameters(&parameters.inception_surfaces, parameters.inception_surface_constraint_weight);
			}
			else {
				karst.disable_inception_horizon();
			}

			karst.set_topo_surface(&parameters.topo_surface);

			if (parameters.use_fracture_constraints) {
				karst.set_fracture_constraint_parameters(&parameters.fracture_families_orientations, &parameters.fracture_families_tolerance, parameters.fracture_constraint_weight);
			}
			else {
				karst.disable_fractures();
			}

			if (parameters.use_no_karst_spheres) {
				karst.set_no_karst_spheres_parameters(&parameters.sphere_centers, parameters.sphere_radius);
			}

			if (parameters.use_amplification_phase_vadose) {
				karst.set_amplification_vadose_params(parameters.max_dist_loops_vadose, parameters.loop_density_vadose);
			}

			if (parameters.use_amplification_phase_phreatic) {
				karst.set_amplification_phreatic_params(parameters.max_dist_loops_phreatic, parameters.loop_density_phreatic);
			}

			karst.set_water_table_level(parameters.water_table_constraint_weight_vadose, parameters.water_table_constraint_weight_phreatic);
			karst.set_simulation_parameters(parameters.nghb_count, parameters.use_max_nghb_radius, parameters.nghb_radius, parameters.poisson_radius, parameters.gamma, parameters.multiply_costs,parameters.vadose_cohesion);

			karst.read_connectivity_matrix(&parameters.sinks, &parameters.springs);

			karst.run_simulation(parameters.create_nghb_graph, parameters.create_nghb_graph_property, parameters.use_amplification_phase_vadose, parameters.use_amplification_phase_phreatic,
				parameters.use_sampling_points, parameters.fraction_karst_perm, parameters.fraction_old_karst_perm, parameters.max_inception_surface_distance, &parameters.sampling_points, parameters.create_vset_sampling, parameters.use_density_property,
				parameters.k_pts, parameters.propdensity, parameters.propikp);
		}
	}
}
