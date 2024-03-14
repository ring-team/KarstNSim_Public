/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/KarsticNetwork.h"

// basic functions :

// get property of given name from gobj (data stored in vector)
//std::vector<double> get_gobj_property_karstnetwork(GObj* gobj, std::string property_name, int fsize) {
//
//	TopologicalObject* topobj = gobj->topological_object();
//	PropertyDB* propdb = gobj->property_db();
//	const PropertyAccessor& accessor(propdb->find_property(property_name)->property_accessor());
//	TopoHandleIterator_var topitr = topobj->handles_itr(topobj->default_region_db_dimension());
//	std::vector<DataValue> vec_data(fsize);
//	int cell_idx = 0;
//	for (; topitr->more(); topitr->next(), cell_idx++) {
//		TopoHandle prop_th = topitr->cur();
//		accessor.get_real_value(prop_th, &vec_data[cell_idx]);
//	}
//	return vec_data;
//}

namespace KarstNSim {

	KarsticNetwork::KarsticNetwork(const std::string& karstic_network_name,  Box* box, GeologicalParameters& params,
		const std::vector<KeyPoint>& keypts, std::vector<Surface>* water_tables) :
		karstic_network_name(karstic_network_name), box(box), params(params), keypts(keypts), water_tables(water_tables) {};

	void KarsticNetwork::set_simulation_parameters(const int& nghb_count, const bool& use_max_nghb_radius, const double& nghb_radius, const double& poisson_radius, const double& gamma, const bool& multiply_costs, const bool& vadose_cohesion) {
		/*!
		Sets the simulation parameters and the simulation domain geometry
		*/

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

	void KarsticNetwork::set_sinks(const std::vector<Vector3>* sinks, const std::vector<int>& propsinksindex, const std::vector<int>& propsinksorder) {
		/*!
		Adds sink points in the keypoint list
		*/

		//std::vector<int> propsinksindex = get_gobj_property_karstnetwork(sinks, sinks_index, int(sinks->size()));
		//std::vector<int> propsinksorder = get_gobj_property_karstnetwork(sinks, sinks_order, int(sinks->size()));

		int max_order;
		max_order = *std::max_element(propsinksorder.begin(), propsinksorder.end());
		int scanning_order = 1;
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
			}
			scanning_order++;
		}
	}

	void KarsticNetwork::set_springs(const std::vector<Vector3>* springs, const std::vector<int>& propspringsindex, const bool& allow_single_outlet_connection) {
		/*!
		Adds spring points in the keypoint list
		*/
		params.allow_single_outlet = allow_single_outlet_connection;
		params.nb_springs = int(springs->size());
		//std::vector<int> propspringsindex = get_gobj_property_karstnetwork(springs, springs_index, params.nb_springs);

		// get all z values for the springs

		for (int i = 0; i < springs->size(); i++)
		{
			int index_for_i = int(std::find(propspringsindex.begin(), propspringsindex.end(), i + 1) - propspringsindex.begin());
			params.z_list.push_back(springs->at(index_for_i).z);
			this->keypts.emplace_back(springs->at(index_for_i), KeyPointType::Spring,i);
			this->pt_spring.push_back(springs->at(index_for_i));
		}
	}

	void KarsticNetwork::set_water_table_order(const std::vector<Surface>* surf_water_table) { // wt have to be ordered correctly

		for (int i = 0; i < surf_water_table->size(); i++) {
			params.wt_surf_index.push_back(i+1); // last element is index
		}
	}

	void KarsticNetwork::set_previous_networks(const std::vector<Line>* previous_networks) {
		/*!
		Adds previous networks to the generated Line
		*/
		int count_total_nb_segs = 0;
		for (int i = 0; i < previous_networks->size(); i++) {
			Line Line_i = previous_networks->at(i);
			count_total_nb_segs += int(Line_i.get_nb_segs());
		}
		params.PtsOldGraph.resize(count_total_nb_segs); // vector of all segments of the old graph, represented each with their start and end point
		params.IdxOldGraph.resize(count_total_nb_segs); // vector of all segments of the old graph, represented each with their start and end point indices
		count_total_nb_segs = 0;
		for (int i = 0; i < previous_networks->size(); i++) {
			Line Line_i = previous_networks->at(i);
			int segsize_i = int(Line_i.get_nb_segs());
			for (int j = 0; j < segsize_i; j++)
			{
				Segment segj = Line_i.get_seg(j);
				params.PtsOldGraph[j + count_total_nb_segs].push_back(segj.start());
				params.PtsOldGraph[j + count_total_nb_segs].push_back(segj.end());
			}
			count_total_nb_segs += segsize_i;
		}
	}

	void KarsticNetwork::set_waypoints(const std::vector<Vector3>* waypoints) {
		/*!
		Adds waypoints in the keypoint list
		*/
		for (int i = 0; i < waypoints->size(); i++)
		{
			this->keypts.emplace_back(waypoints->at(i), KeyPointType::Waypoint);
		}
	}

	void KarsticNetwork::set_inception_surfaces_sampling(const std::string& network_name, std::vector<Surface>* surfaces_used_to_densify, const int& refine_surface_sampling, const bool& create_vset_sampling){
		params.nb_inception_surf = int(surfaces_used_to_densify->size());
		KarstNSim::surface_sampling::multiple_surface_sampling(params.directoryname, network_name, box, surfaces_used_to_densify, nodes_on_inception_surfaces, refine_surface_sampling, create_vset_sampling);
	}

	void KarsticNetwork::set_wt_surfaces_sampling(const std::string& network_name, std::vector<Surface>* surfaces_used_to_densify, const int& refine_surface_sampling) {

		nodes_on_wt_surfaces.resize(surfaces_used_to_densify->size()); // Resize nodes_on_wt_surfaces to match the number of water tables / springs

		for (int i = 0; i < surfaces_used_to_densify->size(); i++) {
			std::vector<Surface> wt_surface;
			wt_surface.push_back(surfaces_used_to_densify->at(i));
			std::vector<Vector3> nodes_on_wt_surface;
			KarstNSim::surface_sampling::multiple_surface_sampling(params.directoryname, network_name, box, &wt_surface,nodes_on_wt_surface , refine_surface_sampling, false);
			nodes_on_wt_surfaces[i] = nodes_on_wt_surface;
		}
	}
	
	void KarsticNetwork::save_painted_box(const std::vector<double>& propdensity, const std::vector<double>& propikp) {

		std::vector<std::string> property_names = { "density","karstif_potential" };
		std::vector<std::vector<double>> properties;
		for (size_t i = 0; i < propdensity.size(); ++i) {
			properties.push_back({ propdensity[i], propikp[i] });
		}
		std::string full_name = params.scenename + "_box.txt";
		std::string full_dir_name = params.directoryname + "/outputs";
		save_box(full_name, full_dir_name, *box, property_names, properties);
	}

	void KarsticNetwork::set_topo_surface(Surface* topo_surface) {
		this->topo_surface_ = topo_surface;
	}

	void KarsticNetwork::set_inception_horizons_parameters( std::vector<Surface>* inception_horizons_list, const double& inception_horizon_constraint_weight) {
		this->inception_horizons = inception_horizons_list;
		params.horizonCost = CostTerm(true, inception_horizon_constraint_weight);

	}

	void KarsticNetwork::disable_inception_horizon() {
		params.horizonCost = CostTerm(false, 0.0);
	}

	void KarsticNetwork::set_karstification_potential_parameters(const double& karstification_potential_weight) {
		params.karstificationCost = CostTerm(true, karstification_potential_weight);
	}

	void KarsticNetwork::set_fracture_constraint_parameters(const std::vector<double>* fracture_families_orientations, const std::vector<double>* fracture_families_tolerance, const double& fracture_constraint_weight) {
		params.fractures_orientations = *fracture_families_orientations;
		params.fractures_tolerances = *fracture_families_tolerance;
		params.fractureCost = CostTerm(true, fracture_constraint_weight);
	}

	void KarsticNetwork::disable_fractures() {
		params.fractureCost = CostTerm(false, 0.0);
	}

	void KarsticNetwork::set_no_karst_spheres_parameters(const std::vector<Vector3>* sphere_centers, const std::vector<double>& sphere_radius) {
		if (sphere_centers->size() == sphere_radius.size()) {
			for (int i = 0; i < sphere_centers->size(); i++) {
				params.spheres.push_back(Sphere(Vector3(0., 0., 0.), 5.));
				params.spheres.push_back(Sphere(sphere_centers->at(i), sphere_radius[i]));
			}
		}
	}

	void KarsticNetwork::set_domain_geometry() {

		/*!
		Sets the scalarfield and maximun and minimum elevation offset
		*/

		const Vector3& min_pt = box->get_basis();
		const Vector3& max_pt = box->get_end();

		double delta_y = max_pt.y - min_pt.y;
		double delta_x = max_pt.x - min_pt.x;
		double delta_z = max_pt.z - min_pt.z;

		double x_center = (max_pt.x + min_pt.x) / 2;
		double y_center = (max_pt.y + min_pt.y) / 2;

		Vector3 min(min_pt.x + delta_x / 4, min_pt.y + delta_y / 4, min_pt.z + delta_z / 4);
		Vector3 max(max_pt.x - delta_x / 4, max_pt.y - delta_y / 4, max_pt.z - delta_z / 4);

		params.heightfield = ScalarField2D(256, 256, Box2D(Vector2(x_center, y_center), delta_x / 2), 0.0);
		params.elevationOffsetMin = sqrt(min.z*min.z);
		params.elevationOffsetMax = sqrt(max.z*max.z);
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

	void createConnectivityMatrix(const std::string& filename, int nb_springs, int nb_sinks) {
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

		std::string connectivity_matrix_path = params.directoryname + "/Input_files/connectivity_matrix.txt";
		int nb_sinks = static_cast<int>(sinks->size());
		int nb_springs = static_cast<int>(springs->size());
		std::ifstream myfile(connectivity_matrix_path);

		// Check if the file exists, if not, create it

		if (!myfile.is_open()) {
			createConnectivityMatrix(connectivity_matrix_path, nb_springs, nb_sinks);
			myfile.open(connectivity_matrix_path);
			if (!myfile.is_open()) {
				std::cerr << "Error opening file after creation: " << strerror(errno) << std::endl;
				return;
			}
		}

		params.connectivity_matrix.resize(nb_sinks, std::vector<int>(nb_springs));

		for (int row = 0; row < nb_sinks; ++row) {
			std::string line;
			if (std::getline(myfile, line)) {

				std::istringstream iss(line);
				for (int col = 0; col < nb_springs; ++col) {
					std::string value;
					if (std::getline(iss, value, '\t')) {
						try {
							params.connectivity_matrix[row][col] = std::stoi(value);
						} catch (const std::invalid_argument& e) {
							std::cerr << "Error converting string to integer: " << e.what() << std::endl;
						}
					} else {
						if (iss.eof()) {
							std::cerr << "Error reading value at row " << row << " column " << col
									  << ": End of file reached unexpectedly." << std::endl;
						} else {
							std::cerr << "Error reading value at row " << row << " column " << col
									  << ": " << strerror(errno) << std::endl;
						}
					}
				}
			} else {
				if (myfile.eof()) {
					std::cerr << "Error reading line " << row << ": End of file reached unexpectedly." << std::endl;
				} else {
					std::cerr << "Error reading line " << row << ": " << strerror(errno) << std::endl;
				}
			}
		}
    myfile.close();
    }

	void KarsticNetwork::run_simulation(const bool& create_nghb_graph, const bool& create_nghb_graph_property, const bool& use_amplification_phase_vadose, const bool& use_amplification_phase_phreatic, const bool& use_sampling_points,
		const double& fraction_karst_perm, const double& fraction_old_karst_perm, const double& max_inception_surface_distance, std::vector<Vector3>* sampling_points, const bool& create_vset_sampling,
		const bool& use_density_property, const int& k_pts, const std::vector<double>& propdensity, const std::vector<double>& propikp) {

		const clock_t begin_time = clock(); // start clock time


		params.scenename = karstic_network_name;
		// Cost due to distance
		const clock_t time1 = clock();
			if (is_simulation_parametrized) {

				// Compute 3D cost graph
				GraphOperations graph;
				graph.InitializeCostGraph(create_nghb_graph, create_nghb_graph_property, keypts, params, nodes_on_inception_surfaces, nodes_on_wt_surfaces,
					inception_horizons, water_tables, use_sampling_points, box,max_inception_surface_distance, sampling_points, create_vset_sampling,
					use_density_property, k_pts, fraction_old_karst_perm, propdensity, propikp, topo_surface_);
				const clock_t time2 = clock();
					// Compute karstic skeleton
					KarsticSkeleton skel = graph.ComputeKarsticSkeleton(keypts, fraction_karst_perm);

				const clock_t time3 = clock();
					// Procedural amplification
					if (use_amplification_phase_vadose) {
						skel.Amplify_vadose(&graph, params);
					}
				const clock_t time4 = clock();

					if (use_amplification_phase_phreatic) {
						skel.Amplify_phreatic(&graph, params);
					}
				const clock_t time5 = clock();

				//skel.create_sections(&graph, params);
				const clock_t time6 = clock();
				//("we generated conduits sections (" << double(time6 - time5) / CLOCKS_PER_SEC << " s)")

				// save network
				skel.CreateLine(params, karstic_network_name);

				const clock_t time7 = clock();

			}
	}

	void KarsticNetwork::set_save_directory(const std::string& repertory) {
		save_directory = repertory;
		params.directoryname = repertory;
	}

	void KarsticNetwork::set_amplification_vadose_params(const double& max_dist_loops_vadose, const double& loop_density_vadose) {
		params.max_dist_loops_vadose = max_dist_loops_vadose;
		params.loop_density_vadose = loop_density_vadose;
	}

	void KarsticNetwork::set_amplification_phreatic_params(const double& max_dist_loops_phreatic, const double& loop_density_phreatic) {
		params.max_dist_loops_phreatic = max_dist_loops_phreatic;
		params.loop_density_phreatic = loop_density_phreatic;
	}

	void KarsticNetwork::set_water_table_level(const double& water_table_constraint_weight_vadose, const double& water_table_constraint_weight_phreatic) {
		params.waterTable1 = CostTerm(true, water_table_constraint_weight_vadose); // vadose
		params.waterTable2 = CostTerm(true, water_table_constraint_weight_phreatic); // phreatic
	}

	void KarsticNetwork::disable_water_table() {
		params.waterTable1 = CostTerm(false, 0.0);
		params.waterTable2 = CostTerm(false, 0.0);
	}
}
//void KarsticNetwork::test(double nghb_count_t, double gamma_t) {
//	const HomogeneousCoordinateSystem* hcs = GeobaseLib::instance()->default_coordinate_system();
//
//	std::vector<Vector3>API::create_from_points(CString("sinkso"), pt_sink, *hcs);
//
//	std::vector<Vector3>API::create_from_points(CString("springso"), pt_spring, *hcs);
//
//	nghb_count_t = params.graphNeighbourCount;
//	gamma_t = params.gamma;
//
//}