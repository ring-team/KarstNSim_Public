/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/run_code.h"

namespace KarstNSim {
	void run_simulation_full(ParamsSource parameters) {

		// Uncomment to save inputs for special use
		std::string full_dir_name = parameters.save_repository + "/outputs";
		std::cout << "Will save outputs to: " << full_dir_name << std::endl;

		// save box voxet
		Box saved_box = parameters.domain;
		std::vector<std::string> property_names = { "density","karstif_potential" };
		std::vector<std::vector<float>> properties;
		for (size_t i = 0; i < parameters.propdensity.size(); ++i) {
			properties.push_back({ parameters.propdensity[i], parameters.propikp[i] });
		}
		std::string full_name = parameters.karstic_network_name + "_box.txt";
		
		save_box(full_name, full_dir_name, saved_box, property_names,properties);

		// save topography
		Surface topo = parameters.topo_surface;
		full_name = parameters.karstic_network_name + "_topo_surf.txt";
		full_dir_name = parameters.save_repository + "/outputs";
		save_surface(full_name, full_dir_name, topo);

		// save inception surfaces
		std::vector<Surface> inception_surfaces = parameters.inception_surfaces;
		full_name = parameters.karstic_network_name + "_inception_surf";
		full_dir_name = parameters.save_repository + "/outputs";
		for (int surfi = 0; surfi < inception_surfaces.size(); surfi++) {
			std::string full_name_i = full_name + std::to_string(surfi) + ".txt";
			save_surface(full_name_i, full_dir_name, inception_surfaces[surfi]);
		}

		// save water tables
		std::vector<Surface> wtsurf = parameters.surf_wat_table;
		full_name = parameters.karstic_network_name + "_watertable_surf";
		full_dir_name = parameters.save_repository + "/outputs";
		for (int surfi = 0; surfi < wtsurf.size(); surfi++) {
			std::string full_name_i = full_name + std::to_string(surfi) + ".txt";
			save_surface(full_name_i, full_dir_name, wtsurf[surfi]);
		}

		// save sinks
		std::vector<Vector3> sinks = parameters.sinks;
		std::vector<std::string> property_names_sinks = { "index","order", "radius" };
		std::vector<std::vector<float>> properties_sinks;
		for (size_t i = 0; i < parameters.propsinksindex.size(); ++i) {
			float v = 0;
			if (parameters.use_sinks_radius) {
				v = parameters.propsinksradius[i];
			}
			properties_sinks.push_back({ float(parameters.propsinksindex[i]), float(parameters.propsinksorder[i]), v });
		}
		full_name = parameters.karstic_network_name + "_sinks.txt";
		full_dir_name = parameters.save_repository + "/outputs";
		save_pointset(full_name, full_dir_name, sinks, property_names_sinks, properties_sinks);

		// save springs
		std::vector<Vector3> springs = parameters.springs;
		std::vector<std::string> property_names_springs = { "index",  "surfindex", "radius" };
		std::vector<std::vector<float>> properties_springs;
		for (size_t i = 0; i < parameters.propspringsindex.size(); ++i) {
			float v = 0;
			if (parameters.use_springs_radius) {
				v = parameters.propspringsradius[i];
			}
			properties_springs.push_back({ float(parameters.propspringsindex[i]), float(parameters.propspringssurfindex[i]), v });
		}
		full_name = parameters.karstic_network_name + "_springs.txt";
		full_dir_name = parameters.save_repository + "/outputs";
		save_pointset(full_name, full_dir_name, springs, property_names_springs, properties_springs);

		// save waypoints
		std::vector<Vector3> waypoints = parameters.waypoints;
		std::vector<std::string> property_names_waypoints = { "impact_radius",  "radius" };
		std::vector<std::vector<float>> properties_waypoints;
		for (size_t i = 0; i < parameters.waypoints_impact_radius.size(); ++i) {
			float v = 0;
			if (parameters.use_waypoints_radius) {
				v = parameters.waypoints_radius[i];
			}
			properties_waypoints.push_back({ float(parameters.waypoints_impact_radius[i]), v });
		}
		full_name = parameters.karstic_network_name + "_waypoints.txt";
		full_dir_name = parameters.save_repository + "/outputs";
		save_pointset(full_name, full_dir_name, waypoints, property_names_waypoints, properties_waypoints);

		for (int i = 0; i < parameters.number_of_iterations; i++) {

			const clock_t begin_time = clock();

			// Initialize with selected seed

			if (parameters.number_of_iterations != 1 && parameters.vary_seed) {
				unsigned int uintValue = static_cast<unsigned int>(parameters.selected_seed + i);
				std::vector<std::uint32_t> seed = { uintValue };
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

			std::string sim_name_iter = parameters.karstic_network_name + "_" + std::to_string(i + parameters.selected_seed - 1);

			KarsticNetwork karst(sim_name_iter, &parameters.domain, params, keypts, &parameters.surf_wat_table);
			karst.set_save_directory(parameters.save_repository);
			karst.set_sinks(&parameters.sinks, parameters.propsinksindex, parameters.propsinksorder, parameters.use_sinks_radius, parameters.propsinksradius);
			karst.set_springs(&parameters.springs, parameters.propspringsindex, parameters.allow_single_outlet_connection, parameters.use_springs_radius, parameters.propspringsradius, parameters.propspringssurfindex);

			if (parameters.use_waypoints) {
				karst.set_waypoints(&parameters.waypoints, parameters.use_waypoints_radius, parameters.waypoints_radius, parameters.waypoints_impact_radius, parameters.waypoints_weight);
			}

			if (parameters.use_deadend_points) {
				karst.set_deadend_points(parameters.nb_deadend_points, parameters.max_distance_of_deadend_pts);
			}

			if (parameters.use_previous_networks) {
				karst.set_previous_networks(&parameters.previous_networks);
			}
			const clock_t timeA = clock();
			if (parameters.use_karstification_potential) {
				karst.set_karstification_potential_parameters(parameters.karstification_potential_weight);

				// that comes after considering karstification potential, as ghostrocks NEED a KP property to be used in the first place (and modify said KP property)
				if (parameters.use_ghostrocks) {
					karst.set_ghost_rocks(parameters.domain, parameters.propikp, &parameters.alteration_lines, parameters.interpolate_lines, parameters.ghostrock_max_vertical_size, parameters.use_max_depth_constraint, parameters.ghost_rock_weight, &parameters.max_depth_horizon, parameters.ghostrock_width);
				}
			}
			const clock_t timeB = clock();
			if (parameters.create_grid) {
				karst.save_painted_box(parameters.propdensity, parameters.propikp);
			}

			if (!parameters.sections_simulation_only) {
				//karst.set_water_table_order(&parameters.surf_wat_table);

				karst.set_wt_surfaces_sampling(sim_name_iter, &parameters.surf_wat_table, parameters.refine_surface_sampling);

				if ((parameters.add_inception_surfaces)) {
					if (!parameters.use_sampling_points) {
						karst.set_inception_surfaces_sampling(sim_name_iter, &parameters.inception_surfaces, parameters.refine_surface_sampling, parameters.create_vset_sampling);
					}
					karst.set_inception_horizons_parameters(&parameters.inception_surfaces, parameters.inception_surface_constraint_weight);
				}
				else {
					karst.disable_inception_horizon();
				}
				const clock_t timeC = clock();
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

				if (parameters.use_amplification) {
					karst.set_amplification_params(parameters.max_distance_amplification, parameters.min_distance_amplification, parameters.nb_cycles);
				}

				karst.set_water_table_weight(parameters.water_table_constraint_weight_vadose, parameters.water_table_constraint_weight_phreatic);
				karst.set_simulation_parameters(parameters.nghb_count, parameters.use_max_nghb_radius, parameters.nghb_radius, parameters.poisson_radius, parameters.gamma, parameters.multiply_costs, parameters.vadose_cohesion);
				karst.read_connectivity_matrix(&parameters.sinks, &parameters.springs);
			}

			if (parameters.simulate_sections) {
				karst.set_geostat_params(parameters.geostat_params);
			}

			//Initialization of noise seed
			karst.set_noise_parameters(parameters.use_noise, parameters.use_noise_on_all, parameters.noise_frequency, parameters.noise_octaves, parameters.noise_weight, globalRng);

			const clock_t time1 = clock();
				float init_time = float(time1 - begin_time) / CLOCKS_PER_SEC;

			float time_needed = karst.run_simulation(parameters.sections_simulation_only, parameters.create_nghb_graph, parameters.create_nghb_graph_property, parameters.use_amplification,
				parameters.use_sampling_points, parameters.fraction_karst_perm, parameters.fraction_old_karst_perm, parameters.max_inception_surface_distance, &parameters.sampling_points, parameters.create_vset_sampling, parameters.use_density_property,
				parameters.k_pts, parameters.propdensity, parameters.propikp);
			const clock_t time2 = clock();
			float real_time_needed = float(time2 - begin_time) / CLOCKS_PER_SEC;

				// File writing
				std::string full_file_name = parameters.save_repository + "/simulation_times.txt";
			std::ofstream outfile(full_file_name, std::ios::app); // Open file in append mode
			if (outfile.is_open()) {
				outfile << parameters.karstic_network_name << ": " << real_time_needed << " s for seed " << parameters.selected_seed << std::endl; // Write the time and seed
				outfile.close(); // Close the file after writing
			}
			else {
			}
		}
	}

	void run_simulation_properties(ParamsSource parameters) {

		for (int i = 0; i < parameters.number_of_iterations; i++) {

			const clock_t begin_time = clock();

			// Initialize with selected seed

			if (parameters.number_of_iterations != 1) {
				unsigned int uintValue = static_cast<unsigned int>(parameters.selected_seed + i);
				std::vector<std::uint32_t> seed = { uintValue };
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

			std::string sim_name_iter = parameters.karstic_network_name + "_" + std::to_string(i + parameters.selected_seed - 1);


			KarsticNetwork karst(sim_name_iter, &parameters.domain, params, keypts, &parameters.surf_wat_table);

			karst.set_save_directory(parameters.save_repository);

			if (parameters.use_waypoints) {
				parameters.use_waypoints_radius = true;
				karst.set_waypoints(&parameters.waypoints, parameters.use_waypoints_radius, parameters.waypoints_radius, parameters.waypoints_impact_radius, parameters.waypoints_weight);
			}

			if (parameters.use_previous_networks) {
				karst.set_previous_networks(&parameters.previous_networks);
			}

			karst.set_geostat_params(parameters.geostat_params);

			std::vector<std::vector<float>> costs_graph(params.PtsOldGraph.size(), std::vector<float>(2, 0.0)); // dummy vector
			std::vector<std::vector<char>> vadoseflags_graph(params.PtsOldGraph.size(), std::vector<char>(2, false)); // dummy vector
			std::vector<int> springidx(params.PtsOldGraph.size(), 1);
			KarsticSkeleton skel(params.PtsOldGraph, costs_graph, vadoseflags_graph, springidx);

			const clock_t time1 = clock();

				karst.run_simulation_properties(skel, &parameters.alteration_lines, parameters.use_ghostrocks, parameters.ghostrock_max_vertical_size, parameters.use_max_depth_constraint, &parameters.max_depth_horizon, parameters.ghostrock_width);
			const clock_t time2 = clock();
		}
	}
}
