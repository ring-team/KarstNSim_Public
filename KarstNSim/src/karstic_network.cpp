/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/karstic_network.h"

namespace KarstNSim {

	KarsticNetwork::KarsticNetwork(const std::string& karstic_network_name, Box* box, GeologicalParameters& params,
		const std::vector<KeyPoint>& keypts, std::vector<Surface>* water_tables) :
		karstic_network_name(karstic_network_name), box(box), params(params), keypts(keypts), water_tables(water_tables) {};

	void KarsticNetwork::set_simulation_parameters(const int& nghb_count, const bool& use_max_nghb_radius, const float& nghb_radius, const float& poisson_radius, const float& gamma, const bool& multiply_costs, const bool& vadose_cohesion) {

		params.gamma = gamma;
		params.graphNeighbourCount = nghb_count;
		params.graphuse_max_nghb_radius = use_max_nghb_radius;
		params.graphNeighbourRadius = nghb_radius;
		params.graphPoissonRadius = poisson_radius;
		params.distanceCost = CostTerm(true, 1);
		params.multiply_costs = multiply_costs;
		params.vadose_cohesion = vadose_cohesion;
		is_simulation_parametrized = true;
		set_domain_geometry();
	}

	void KarsticNetwork::set_sinks(const std::vector<Vector3>* sinks, const std::vector<int>& propsinksindex, const std::vector<int>& propsinksorder, bool use_sinks_radius, const std::vector<float>& propsinksradius) {

		//std::vector<int> propsinksindex = get_gobj_property_karstnetwork(sinks, sinks_index, int(sinks->size()));
		//std::vector<int> propsinksorder = get_gobj_property_karstnetwork(sinks, sinks_order, int(sinks->size()));

		int max_order;
		max_order = *std::max_element(propsinksorder.begin(), propsinksorder.end());
		int scanning_order = 1;
		use_sinks_radius_ = use_sinks_radius;
		// sinks are appended in a specific order to the list :
		// 1) they are appended in increasing order of importance (priority to bigger sinks)
		// 2) among sinks of same order, the order is shuffled randomly

		while (scanning_order <= max_order) {
			std::vector<int> propsinks_iter;
			std::vector<int> propsinks_iter_index;
			std::vector<int> indexes;
			for (int i = 0; i < sinks->size(); i++) {
				if (propsinksorder[i] == scanning_order) {
					propsinks_iter.push_back(i);
					propsinks_iter_index.push_back(propsinksindex[i]);
				}
			}
			for (int i = 0; i < int(propsinks_iter.size()); ++i) {
				indexes.push_back(i);
			}
			shuffleContainer(indexes.begin(), indexes.end());
			reorder(propsinks_iter, indexes);
			reorder(propsinks_iter_index, indexes);
			for (int i = 0; i < int(propsinks_iter.size()); i++) {
				params.sinks_index.push_back(propsinks_iter_index[i]); // insert shuffled vector of indices of sinks of order i
				this->keypts.emplace_back(sinks->at(propsinks_iter[i]), KeyPointType::Sink);
				this->pt_sink.push_back(sinks->at(propsinks_iter[i]));

				if (use_sinks_radius_) {
					propsinksradius_.push_back({ propsinksradius[propsinks_iter[i]], int(keypts.size() - 1) });
				}
			}
			scanning_order++;
		}
	}

	void KarsticNetwork::set_springs(const std::vector<Vector3>* springs, const std::vector<int>& propspringsindex, const bool& allow_single_outlet_connection, bool use_springs_radius, const std::vector<float>& propspringsradius, const std::vector<int>& propspringswtindex) {

		params.allow_single_outlet = allow_single_outlet_connection;
		params.nb_springs = int(springs->size());
		//std::vector<int> propspringsindex = get_gobj_property_karstnetwork(springs, springs_index, params.nb_springs);

		// get all z values for the springs
		use_springs_radius_ = use_springs_radius;

		for (int i = 0; i < springs->size(); i++)
		{
			int index_for_i = int(std::find(propspringsindex.begin(), propspringsindex.end(), i + 1) - propspringsindex.begin());

			this->keypts.emplace_back(springs->at(index_for_i), KeyPointType::Spring, propspringswtindex[index_for_i]);
			params.z_list.push_back({ springs->at(index_for_i).z,  int(keypts.size() - 1) });
			this->pt_spring.push_back(springs->at(index_for_i));
			params.propspringswtindex.push_back({ static_cast<float>(propspringswtindex[index_for_i]), int(keypts.size() - 1) });
			if (use_springs_radius_) {
				propspringsradius_.push_back({ propspringsradius[index_for_i], int(keypts.size() - 1) });
			}
		}
	}

	void KarsticNetwork::set_waypoints(const std::vector<Vector3>* waypoints, bool use_waypoints_radius, const std::vector<float>& propwaypointsradius, const std::vector<float>& propwaypointsimpactradius, float waypoints_weight) {

		use_waypoints_radius_ = use_waypoints_radius;
		params.waypoints_weight = waypoints_weight;
		for (int i = 0; i < waypoints->size(); i++)
		{
			this->keypts.emplace_back(waypoints->at(i), KeyPointType::Waypoint);
			if (use_waypoints_radius_) {
				propwaypointsradius_.push_back({ propwaypointsradius[i], int(keypts.size() - 1) });
			}
			params.waypointsimpactradius.push_back({ propwaypointsimpactradius[i], int(keypts.size() - 1) });
		}
	}

	//void KarsticNetwork::set_water_table_order(const std::vector<Surface>* surf_water_table) { // wt have to be ordered correctly

	//	for (int i = 0; i < surf_water_table->size(); i++) {
	//		params.wt_surf_index.push_back(i+1); // last element is index
	//	}
	//}

	void KarsticNetwork::set_previous_networks(const std::vector<Line>* previous_networks) {

		int count_total_nb_segs = 0;
		for (int i = 0; i < previous_networks->size(); i++) {
			Line Line_i = previous_networks->at(i);
			count_total_nb_segs += int(Line_i.get_nb_segs());
		}
		params.PtsOldGraph.resize(count_total_nb_segs, 2); // vector of all segments of the old graph, represented each with their start and end point
		params.IdxOldGraph.resize(count_total_nb_segs, 2); // vector of all point of the old graph, each represented by its index in the samples object of GraphOperations
		count_total_nb_segs = 0;
		for (int i = 0; i < previous_networks->size(); i++) {
			Line Line_i = previous_networks->at(i);
			int segsize_i = int(Line_i.get_nb_segs());
			for (int j = 0; j < segsize_i; j++)
			{
				Segment segj = Line_i.get_seg(j);
				params.PtsOldGraph[j + count_total_nb_segs][0] = segj.start();
				params.PtsOldGraph[j + count_total_nb_segs][1] = segj.end();
			}
			count_total_nb_segs += segsize_i;
		}
	}

	void KarsticNetwork::set_deadend_points(int nb_deadend_points, float max_distance_of_deadend_pts) {
		use_deadend_pts_ = true;
		nb_deadend_points_ = nb_deadend_points;
		max_distance_of_deadend_pts_ = max_distance_of_deadend_pts;
	}

	void KarsticNetwork::set_inception_surfaces_sampling(const std::string& network_name, std::vector<Surface>* surfaces_used_to_densify, const int& refine_surface_sampling, const bool& create_vset_sampling) {
		params.nb_inception_surf = int(surfaces_used_to_densify->size());
		KarstNSim::surface_sampling::multiple_surface_sampling(params.directoryname, network_name, box, surfaces_used_to_densify, nodes_on_inception_surfaces, refine_surface_sampling, create_vset_sampling);
	}

	void KarsticNetwork::set_wt_surfaces_sampling(const std::string& network_name, std::vector<Surface>* surfaces_used_to_densify, const int& refine_surface_sampling) {

		nodes_on_wt_surfaces.resize(surfaces_used_to_densify->size()); // Resize nodes_on_wt_surfaces to match the number of water tables
		params.nb_wt = surfaces_used_to_densify->size();

		for (int i = 0; i < surfaces_used_to_densify->size(); i++) {
			std::vector<Surface> wt_surface;
			wt_surface.push_back(surfaces_used_to_densify->at(i));
			std::vector<Vector3> nodes_on_wt_surface;
			KarstNSim::surface_sampling::multiple_surface_sampling(params.directoryname, network_name, box, &wt_surface, nodes_on_wt_surface, refine_surface_sampling, false);
			nodes_on_wt_surfaces[i] = nodes_on_wt_surface;
		}
	}

	void KarsticNetwork::save_painted_box(std::vector<float>& propdensity, const std::vector<float>& propikp) {

		std::vector<std::string> property_names = { "density","karstif_potential" };
		if (propdensity.empty()) { // if propdensity wasnt provided, we can just fill it with default values (-99999).
			propdensity.resize(propikp.size(), -99999);
		}
		std::vector<std::vector<float>> properties;
		for (size_t i = 0; i < propikp.size(); ++i) {
			properties.push_back({ propdensity.at(i) , propikp.at(i) });
		}
		std::string full_name = karstic_network_name + "_box.txt";
		std::string full_dir_name = params.directoryname + "/outputs";
		save_box(full_name, full_dir_name, *box, property_names, properties);
	}

	void KarsticNetwork::set_topo_surface(Surface* topo_surface) {
		this->topo_surface_ = topo_surface;
	}

	void KarsticNetwork::set_inception_horizons_parameters(std::vector<Surface>* inception_horizons_list, const float& inception_horizon_constraint_weight) {
		this->inception_horizons = inception_horizons_list;
		params.horizonCost = CostTerm(true, inception_horizon_constraint_weight);
	}

	void KarsticNetwork::set_ghost_rocks(const Box& grid, std::vector<float>& ikp, Line* alteration_lines, const bool& interpolate_lines, const float& ghostrock_max_vertical_size, const bool& use_max_depth_constraint, const float& ghost_rock_weight, Surface* max_depth_horizon, const float& ghostrock_width) {

		params.use_ghost_rocks = true;
		params.length = ghostrock_max_vertical_size;
		params.width = ghostrock_width;
		params.polyline = *alteration_lines;
		params.use_max_depth_constraint = use_max_depth_constraint;
		params.substratum_surf = *max_depth_horizon;

		// modify "ikp" object with ghostrocks ("paint" it). Note that density property is NOT changed, hence passed as const ref
		paint_KP_with_ghostrocks(grid, ikp, ghostrock_max_vertical_size, ghostrock_width, *alteration_lines, use_max_depth_constraint, *max_depth_horizon, ghost_rock_weight);
	}

	void KarsticNetwork::disable_inception_horizon() {
		params.horizonCost = CostTerm(false, 0.0);
	}

	void KarsticNetwork::set_karstification_potential_parameters(const float& karstification_potential_weight) {
		params.karstificationCost = CostTerm(true, karstification_potential_weight);
	}

	void KarsticNetwork::set_fracture_constraint_parameters(const std::vector<float>* fracture_families_orientations, const std::vector<float>* fracture_families_tolerance, const float& fracture_constraint_weight) {
		params.fractures_orientations = *fracture_families_orientations;
		params.fractures_tolerances = *fracture_families_tolerance;
		params.fractureCost = CostTerm(true, fracture_constraint_weight);
	}

	void KarsticNetwork::disable_fractures() {
		params.fractureCost = CostTerm(false, 0.0);
	}

	void KarsticNetwork::set_no_karst_spheres_parameters(const std::vector<Vector3>* sphere_centers, const std::vector<float>& sphere_radius) {
		for (int i = 0; i < sphere_centers->size(); i++) {
			//params.spheres.push_back(Sphere(Vector3(0., 0., 0.), 5.));
			params.spheres.push_back(Sphere(sphere_centers->at(i), sphere_radius[0]));
		}
	}

	void KarsticNetwork::set_domain_geometry() {

		const Vector3& min_pt = box->get_basis();
		const Vector3& max_pt = box->get_end();

		float delta_y = max_pt.y - min_pt.y;
		float delta_x = max_pt.x - min_pt.x;
		float delta_z = max_pt.z - min_pt.z;

		float x_center = (max_pt.x + min_pt.x) / 2;
		float y_center = (max_pt.y + min_pt.y) / 2;

		Vector3 min(min_pt.x + delta_x / 4, min_pt.y + delta_y / 4, min_pt.z + delta_z / 4);
		Vector3 max(max_pt.x - delta_x / 4, max_pt.y - delta_y / 4, max_pt.z - delta_z / 4);

		params.maxsize = (std::abs(delta_x) + std::abs(delta_y)) / 2;
		params.stretch_factor = params.maxsize / std::abs(delta_z);
	}

	int KarsticNetwork::just_sampling() {

		params.scenename = karstic_network_name;

		if (is_simulation_parametrized) {

			// Initialize the graph, because the sampling is in this class
			GraphOperations graph;
			//graph.Sampling(nodes_on_densify_surface, Box);

			return graph.get_nb_sampling_pts();
		}
		return 0;
	}

	void create_connectivity_matrix(const std::string& filename, int nb_springs, int nb_sinks) {

		// by default, creates a connectivity matrix filled with '2'.
		// it should be used with the shortest distance option enabled.

		// Remove the file if it already exists
		if (std::remove(filename.c_str()) != 0 && errno != ENOENT) {
			std::cerr << "Error deleting existing file: " << strerror(errno) << std::endl;
				return;
		}

		std::ofstream outfile(filename);
		if (!outfile.is_open()) {
				std::cerr << "Error creating file: " << strerror(errno) << std::endl;
			return;
		}

		// Write default values to the file
		for (int i = 0; i < nb_sinks; ++i) {
			for (int j = 0; j < nb_springs; ++j) {
				outfile << "2"; // Default value
				if (j < nb_springs - 1) {
					outfile << "\t"; // Separate columns by tabs
				}
			}
			if (i < nb_sinks - 1) {
				outfile << "\n"; // Move to the next line
			}
		}
		outfile.close();
	}

	void KarsticNetwork::read_connectivity_matrix(const std::vector<Vector3>* sinks, const std::vector<Vector3>* springs) {

		std::string connectivity_matrix_path = params.directoryname + "/connectivity_matrix.txt";
		int nb_sinks = static_cast<int>(sinks->size());
		int nb_springs = static_cast<int>(springs->size());

		// Check if the file exists, if not, create it

		// ***************************
		//create_connectivity_matrix(connectivity_matrix_path, nb_springs, nb_sinks); // UNCOMMENT TO USE DEFAULT MATRIX (1's EVERYWERE)
		// ***************************

		std::ifstream myfile(connectivity_matrix_path);
		if (!myfile.is_open()) {
			std::cerr << "Error opening file " << connectivity_matrix_path << " after creation: " << strerror(errno) << std::endl;
			return;
		}

		params.connectivity_matrix.resize(nb_sinks, nb_springs, 1);

		for (int row = 0; row < nb_sinks; ++row) {
			std::string line;
			if (std::getline(myfile, line)) {

				std::istringstream iss(line);
				for (int col = 0; col < nb_springs; ++col) {
					std::string value;
					if (std::getline(iss, value, '\t')) {
						try {
							params.connectivity_matrix[row][col] = std::stoi(value);
						}
						catch (const std::invalid_argument& e) {
							std::cerr << "Error converting string to integer: " << e.what() << std::endl;
						}
					}
					else {
						if (iss.eof()) {
							std::cerr << "Error reading value at row " << row << " column " << col
								<< ": End of file reached unexpectedly." << std::endl;
						}
						else {
							std::cerr << "Error reading value at row " << row << " column " << col
								<< ": " << strerror(errno) << std::endl;
						}
					}
				}
			}
			else {
				if (myfile.eof()) {
					std::cerr << "Error reading line " << row << ": End of file reached unexpectedly." << std::endl;
				}
				else {
					std::cerr << "Error reading line " << row << ": " << strerror(errno) << std::endl;
				}
			}
		}
		myfile.close();
	}

	void KarsticNetwork::create_sections(KarsticSkeleton& skel) {

		if (geostatparams.is_used) {
			skel.prepare_graph(); // removes duplicates, and changes format so that each node is connected to all of its neighbors (not the case at base necesarily)
			skel.update_branch_ID(); // create branch id property on the skeleton nodes (-1 if intersection, branch index >=1 otherwise)
			skel.compute_distance_matrix(); // compute distance matrix between each pair of nodes in the graph
			skel.compute_branch_sizes(); // compute the number of nodes in each branch
			skel.compute_valence(); // compute the valence = number of neighbors of each node
		}

		// Uncomment to print the distance matrix (in case of need to debug!) and to save

		//std::string matrixString = "Distance Matrix:\n";
		//int i = 0;
		//for (const auto& row : skel.distance_mat) {
		//	std::string rowString;
		//	if (i == 0) {
		//		int j = 0;
		//		for (float distance : row) {
		//			skel.nodes[j].distance = distance;
		//			rowString += std::to_string(distance) + " ";
		//			j++;
		//		}
		//		matrixString += rowString + "\n";
		//	}
		//	i++;
		//}

		if (geostatparams.is_used) {

			if (use_sinks_radius_ || use_springs_radius_ || use_waypoints_radius_) {
				geostatparams.simulated_property.resize(skel.nodes.size());
				std::fill(geostatparams.simulated_property.begin(), geostatparams.simulated_property.end(), -99999);
			}

			// Set conditioning data (inlets and/or outlets and/or waypoints)

			if (use_sinks_radius_) {
				for (int j = 0; j < propsinksradius_.size(); j++) {
					for (int i = 0; i < skel.nodes.size(); i++) {
						if (skel.nodes[i].p == keypts[propsinksradius_[j].index].p) {
							geostatparams.simulated_property[i] = propsinksradius_[j].prop;
							break;
						}
					}
				}
			}

			if (use_springs_radius_) {
				for (int j = 0; j < propspringsradius_.size(); j++) {
					for (int i = 0; i < skel.nodes.size(); i++) {
						if (skel.nodes[i].p == keypts[propspringsradius_[j].index].p) {
							geostatparams.simulated_property[i] = propspringsradius_[j].prop;
							break;
						}
					}
				}
			}

			if (use_waypoints_radius_) {
				for (int j = 0; j < propwaypointsradius_.size(); j++) {
					for (int i = 0; i < skel.nodes.size(); i++) {
						if (skel.nodes[i].p == keypts[propwaypointsradius_[j].index].p) {
							geostatparams.simulated_property[i] = propwaypointsradius_[j].prop;
							break;
						}
					}
				}
			}

			// Simulate property with Frantz's modified SGS algorithm (with 3 variograms) (2021)
			SGS3(
				&skel,
				geostatparams.simulated_property,
				&geostatparams.simulation_distribution,
				geostatparams.global_vario_range,
				geostatparams.global_range_of_neighborhood,
				geostatparams.global_vario_sill,
				geostatparams.global_vario_nugget,
				geostatparams.global_vario_model,
				geostatparams.interbranch_vario_range,
				geostatparams.interbranch_range_of_neighborhood,
				geostatparams.interbranch_vario_sill,
				geostatparams.interbranch_vario_nugget,
				geostatparams.interbranch_vario_model,
				geostatparams.intrabranch_vario_range,
				geostatparams.intrabranch_range_of_neighborhood,
				geostatparams.intrabranch_vario_sill,
				geostatparams.intrabranch_vario_nugget,
				geostatparams.intrabranch_vario_model,
				geostatparams.number_max_of_neighborhood_points,
				geostatparams.nb_points_interbranch,
				geostatparams.proportion_interbranch);

			// assign property to skeleton :

			for (int i = 0; i < geostatparams.simulated_property.size(); i++) {
				skel.nodes.at(i).eq_radius = geostatparams.simulated_property.at(i);
			}
		}

		//Finally we paint any node that is within a ghostrock to have the ghost rock width

		if (params.use_ghost_rocks) {
			paint_karst_sections_with_ghostrocks(skel, params.length, params.width, params.polyline, params.use_max_depth_constraint, params.substratum_surf);
		}
	}

	void KarsticNetwork::run_simulation_properties(KarsticSkeleton& skel, Line* alteration_lines, const bool& use_ghost_rocks, const float& ghostrock_max_vertical_size, const bool& use_max_depth_constraint, Surface* max_depth_horizon, const float& ghostrock_width) {

		params.use_ghost_rocks = use_ghost_rocks;
		params.length = ghostrock_max_vertical_size;
		params.width = ghostrock_width;
		params.polyline = *alteration_lines;
		params.use_max_depth_constraint = use_max_depth_constraint;
		params.substratum_surf = *max_depth_horizon;

		skel.prepare_graph(); // removes duplicates, and changes format so that each node is connected to all of its neighbors (not the case at base necesarily)
		skel.update_branch_ID(); // create branch id property on the skeleton nodes (-1 if intersection, branch index >=1 otherwise)
		skel.compute_distance_matrix(); // compute distance matrix between each pair of nodes in the graph
		skel.compute_branch_sizes(); // compute the number of nodes in each branch
		skel.compute_valence(); // compute the valence = number of neighbors of each node

		use_waypoints_radius_ = true;

		if (geostatparams.is_used) {

			if (use_waypoints_radius_) {
				geostatparams.simulated_property.resize(skel.nodes.size());
				std::fill(geostatparams.simulated_property.begin(), geostatparams.simulated_property.end(), -99999);
			}

			// Set conditioning data (inlets and/or outlets and/or waypoints)



			if (use_waypoints_radius_) {
				for (int j = 0; j < propwaypointsradius_.size(); j++) {
					for (int i = 0; i < skel.nodes.size(); i++) {
						if (skel.nodes[i].p == keypts[propwaypointsradius_[j].index].p) {
							geostatparams.simulated_property[i] = propwaypointsradius_[j].prop;
							break;
						}
					}
				}
			}

			// Simulate property with Frantz's modified SGS algorithm (with 3 variograms) (2021)
			SGS3(
				&skel,
				geostatparams.simulated_property,
				&geostatparams.simulation_distribution,
				geostatparams.global_vario_range,
				geostatparams.global_range_of_neighborhood,
				geostatparams.global_vario_sill,
				geostatparams.global_vario_nugget,
				geostatparams.global_vario_model,
				geostatparams.interbranch_vario_range,
				geostatparams.interbranch_range_of_neighborhood,
				geostatparams.interbranch_vario_sill,
				geostatparams.interbranch_vario_nugget,
				geostatparams.interbranch_vario_model,
				geostatparams.intrabranch_vario_range,
				geostatparams.intrabranch_range_of_neighborhood,
				geostatparams.intrabranch_vario_sill,
				geostatparams.intrabranch_vario_nugget,
				geostatparams.intrabranch_vario_model,
				geostatparams.number_max_of_neighborhood_points,
				geostatparams.nb_points_interbranch,
				geostatparams.proportion_interbranch);

			// assign property to skeleton :

			for (int i = 0; i < geostatparams.simulated_property.size(); i++) {
				skel.nodes.at(i).eq_radius = geostatparams.simulated_property.at(i);
			}
		}

		//Finally we paint any node that is within a ghostrock to have the ghost rock width

		if (params.use_ghost_rocks) {
			paint_karst_sections_with_ghostrocks(skel, params.length, params.width, params.polyline, params.use_max_depth_constraint, params.substratum_surf);
		}
	}

	float KarsticNetwork::run_simulation(const bool& sections_simulation_only, const bool& create_nghb_graph, const bool& create_nghb_graph_property, const bool& use_amplification, const bool& use_sampling_points,
		const float& fraction_karst_perm, const float& fraction_old_karst_perm, const float& max_inception_surface_distance, std::vector<Vector3>* sampling_points, const bool& create_vset_sampling,
		const bool& use_density_property, const int& k_pts, const std::vector<float>& propdensity, const std::vector<float>& propikp) {

			params.scenename = karstic_network_name;
		// Cost due to distance
		const clock_t time1 = clock();
		float time_needed = 0.0f;

		if (is_simulation_parametrized) { // full simulation

			// Compute 3D cost graph

			GraphOperations graph;
			graph.InitializeCostGraph(create_nghb_graph, create_nghb_graph_property, keypts, params, nodes_on_inception_surfaces, nodes_on_wt_surfaces,
				inception_horizons, water_tables, use_sampling_points, box, max_inception_surface_distance, sampling_points, create_vset_sampling,
				use_density_property, k_pts, fraction_old_karst_perm, propdensity, propikp, topo_surface_);
			const clock_t time2 = clock();

				// Compute karstic skeleton

				std::vector<std::vector<int>> karst_paths;
			std::vector<std::vector<float>> karst_paths_costs;
			std::vector<std::vector<char>> karst_paths_vadose_flag;
			std::vector<int> springidxFinal;
			graph.ComputeKarsticSkeleton(keypts, fraction_karst_perm, karst_paths, karst_paths_costs, karst_paths_vadose_flag, springidxFinal);

			if (!karst_paths.empty()) { // if a network was successfully computed

				// Build karstic skeleton structure
				KarsticSkeleton skel(&graph, karst_paths, karst_paths_costs, karst_paths_vadose_flag, springidxFinal);

				const clock_t time3 = clock();

					// Network preparation
					skel.detect_intersection_points(karst_paths);
				skel.update_branch_ID(karst_paths);

				// Procedural amplification with deadend nodes
				const clock_t time4 = clock();
				if (use_deadend_pts_) {
					skel.amplify_deadend(&graph, max_distance_of_deadend_pts_, nb_deadend_points_, params);
				}
				const clock_t time5 = clock();

					// Amplification
					if (use_amplification) {
						if (params.use_noise && !params.use_noise_on_all) {
							graph.add_noise();
						}
						std::pair<float, float> result = skel.amplify_noise(&graph, params);
						const clock_t time5bis = clock();
					}

				const clock_t time6 = clock();
				skel.detect_intersection_points(karst_paths);
				skel.update_branch_ID(karst_paths);

				create_sections(skel);
				const clock_t time7 = clock();
					// save network
					skel.create_line(params, karstic_network_name);

				const clock_t time8 = clock();

					const clock_t time_end = clock();
				time_needed = float(time_end - time1) / CLOCKS_PER_SEC;
			}
			else {
			}

		}
		else if (sections_simulation_only) { // only generate sections
			std::vector<std::vector<float>> costs_graph(params.PtsOldGraph.size(), std::vector<float>(2, 0.0)); // dummy vector
			std::vector<std::vector<char>> vadoseflags_graph(params.PtsOldGraph.size(), std::vector<char>(2, false)); // dummy vector
			std::vector<int> springidx(params.PtsOldGraph.size(), 1);
			KarsticSkeleton skel(params.PtsOldGraph, costs_graph, vadoseflags_graph, springidx);

			const clock_t time5 = clock();

			create_sections(skel);

			const clock_t time6 = clock();

				// save network
				skel.create_line(params, karstic_network_name);

			const clock_t time7 = clock();
		}
		return time_needed;
	}

	void KarsticNetwork::set_save_directory(const std::string& repertory) {
		save_directory = repertory;
		params.directoryname = repertory;
	}

	void KarsticNetwork::set_amplification_params(const float& max_distance_amplification, const float& min_distance_amplification, const int& nb_cycles) {
		params.max_distance_amplification = max_distance_amplification;
		params.min_distance_amplification = min_distance_amplification;
		params.nb_cycles = nb_cycles;

	}

	void KarsticNetwork::set_geostat_params(const GeostatParams& geostat_params) {
		geostatparams = geostat_params;
	}

	void KarsticNetwork::set_amplification_vadose_params(const float& max_dist_loops_vadose, const float& loop_density_vadose) {
		params.max_dist_loops_vadose = max_dist_loops_vadose;
		params.loop_density_vadose = loop_density_vadose;
	}

	void KarsticNetwork::set_amplification_phreatic_params(const float& max_dist_loops_phreatic, const float& loop_density_phreatic) {
		params.max_dist_loops_phreatic = max_dist_loops_phreatic;
		params.loop_density_phreatic = loop_density_phreatic;
	}

	void KarsticNetwork::set_water_table_weight(const float& water_table_constraint_weight_vadose, const float& water_table_constraint_weight_phreatic) {
		params.waterTable1 = CostTerm(true, water_table_constraint_weight_vadose); // vadose
		params.waterTable2 = CostTerm(true, water_table_constraint_weight_phreatic); // phreatic
	}

	void KarsticNetwork::disable_water_table() {
		params.waterTable1 = CostTerm(false, 0.0);
		params.waterTable2 = CostTerm(false, 0.0);
	}

	void KarsticNetwork::set_noise_parameters(const bool use_noise, const bool use_noise_on_all, const int frequency, const int octaves, const float noise_weight, std::mt19937 globalRng) {
		if (use_noise) {
			SimplexNoise Perlin = SimplexNoise(float(frequency), 1.0f, 2.0f, 0.5f);
			Perlin.initialize_permutation_table(globalRng);
			params.noise_frequency = frequency;
			params.noise_octaves = octaves;
			params.noise_weight = noise_weight;
		}
		params.use_noise = use_noise;
		params.use_noise_on_all = use_noise_on_all;
	}
}