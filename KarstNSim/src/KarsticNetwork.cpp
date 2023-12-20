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
	// reorder vector with given vector of indices
	void reorder(std::vector<int> &vA, std::vector<int> vOrder)
	{
		assert(vA.size() == vOrder.size());
		// for all elements to put in place
		for (int i = 0; i < vA.size() - 1; ++i)
		{
			// while the element i is not yet in place 
			while (i != vOrder[i])
			{
				// swap it with the element at its final place
				int alt = vOrder[i];
				std::swap(vA[i], vA[alt]);
				std::swap(vOrder[i], vOrder[alt]);
			}
		}
	}

	//###########################################################################################################################

	KarsticNetwork::KarsticNetwork(const std::string& karstic_network_name, Box* box, GeologicalParameters& params,
		std::vector<KeyPoint>& keypts, std::vector<Surface>* water_tables) :
		karstic_network_name(karstic_network_name), box(box), params(params), keypts(keypts), water_tables(water_tables) {};

	void KarsticNetwork::set_simulation_parameters(int nghb_count, const bool use_max_nghb_radius, double nghb_radius, double poisson_radius, double gamma, const bool multiply_costs) {
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
		is_simulation_parametrized = true;
		set_domain_geometry();
	}

	void KarsticNetwork::set_sinks(std::vector<Vector3>* sinks, std::vector<int> propsinksindex, std::vector<int> propsinksorder) {
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
		propsinksindex.clear();
		propsinksorder.clear();
	}

	void KarsticNetwork::set_springs(std::vector<Vector3>* springs, std::vector<int> propspringsindex, const bool& allow_single_outlet_connection) {
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
			this->keypts.emplace_back(springs->at(index_for_i), KeyPointType::Spring);
			this->pt_spring.push_back(springs->at(index_for_i));
		}
		propspringsindex.clear();
	}

	void KarsticNetwork::set_water_table_order(std::vector<Surface>* surf_water_table) { // wt have to be ordered correctly

		for (int i = 0; i < surf_water_table->size(); i++) {
			params.wt_surf_index.push_back(i+1); // last element is index
		}
	}

	void KarsticNetwork::set_previous_networks(std::vector<Line>* previous_networks) {
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

	void KarsticNetwork::set_waypoints(std::vector<Vector3>* waypoints) {
		/*!
		Adds waypoints in the keypoint list
		*/
		for (int i = 0; i < waypoints->size(); i++)
		{
			this->keypts.emplace_back(waypoints->at(i), KeyPointType::Waypoint);
		}
	}

	void KarsticNetwork::set_surfaces_to_densify(std::string network_name, std::vector<Surface>* surfaces_used_to_densify, const int refine_surface_sampling, const bool create_vset_sampling) {
		KarstNSim::surface_sampling::multiple_surface_sampling(params.directoryname, network_name, box, surfaces_used_to_densify, nodes_on_densify_surface, refine_surface_sampling, create_vset_sampling);
	}

	void KarsticNetwork::set_topo_surface(Surface* topo_surface) {
		this->topo_surface_ = topo_surface;
	}

	void KarsticNetwork::set_inception_horizons_parameters(std::vector<Surface>* inception_horizons_list,
		double inception_horizon_constraint_weight) {
		this->inception_horizons = inception_horizons_list;
		params.horizonCost = CostTerm(true, inception_horizon_constraint_weight);

	}

	void KarsticNetwork::disable_inception_horizon() {
		params.horizonCost = CostTerm(false, 0.0);
	}

	void vectorial_product(Vector3 vec1, Vector3 vec2, Vector3& result) {
		result.x = vec1.y*vec2.z - vec1.z*vec2.y;
		result.y = vec1.z*vec2.x - vec1.x*vec2.z;
		result.z = vec1.x*vec2.y - vec1.y*vec2.x;
	}

	void KarsticNetwork::set_karstification_potential_parameters(const double karstification_potential_weight) {
		params.karstificationCost = CostTerm(true, karstification_potential_weight);
	}

	void KarsticNetwork::set_fracture_constraint_parameters(std::vector<double>* fracture_families_orientations, std::vector<double>* fracture_families_tolerance, const double fracture_constraint_weight) {
		params.fractures_orientations = *fracture_families_orientations;
		params.fractures_tolerances = *fracture_families_tolerance;
		params.fractureCost = CostTerm(true, fracture_constraint_weight);
	}

	void KarsticNetwork::disable_fractures() {
		params.fractureCost = CostTerm(false, 0.0);
	}

	void KarsticNetwork::set_no_karst_spheres_parameters(std::vector<Vector3>* sphere_centers, std::vector<double> sphere_radius) {
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

		Vector3 minp(min_pt.x + delta_x / 4, min_pt.y + delta_y / 4, min_pt.z + delta_z / 4);
		Vector3 maxp(max_pt.x - delta_x / 4, max_pt.y - delta_y / 4, max_pt.z - delta_z / 4);

		params.heightfield = ScalarField2D(256, 256, Box2D(Vector2(x_center, y_center), delta_x / 2), 0.0);
		params.elevationOffsetMin = sqrt(minp.z*minp.z);
		params.elevationOffsetMax = sqrt(maxp.z*maxp.z);
		params.maxsize = (std::abs(max_pt.x - min_pt.x) + std::abs(max_pt.y - min_pt.y)) / 2;
	}

	int KarsticNetwork::just_sampling() {

		params.scenename = karstic_network_name;

		if (is_simulation_parametrized) {

			// Initialize the graph, because the sampling is in this class
			GraphOperations graph;
			//graph.SamplingAdapted(nodes_on_densify_surface, Box);

			return graph.get_nb_sampling_pts();
		}
		return 0;
	}

	void KarsticNetwork::read_connectivity_matrix(std::vector<Vector3>* sinks, std::vector<Vector3>* springs) {

		std::string connectivity_matrix_path = params.directoryname + "/Inputs_files/connectivity_matrix.txt";
		std::ifstream myfile(connectivity_matrix_path);
		int nb_sinks = int(sinks->size());
		int nb_springs = int(springs->size());
		params.connectivity_matrix.resize(nb_sinks, std::vector<int>(nb_springs));
		for (int row = 0; row < nb_sinks; ++row)
			for (int col = 0; col < nb_springs; ++col)
				myfile >> params.connectivity_matrix[row][col];
	}

	void KarsticNetwork::run_simulation(const bool create_nghb_graph, const bool create_nghb_graph_property, const bool use_amplification_phase_vadose, const bool use_amplification_phase_phreatic,
		const bool use_sampling_points, const double fraction_karst_perm, const double fraction_old_karst_perm, const double max_inception_surface_distance,
		std::vector<Vector3>* sampling_points, const bool create_vset_sampling, const bool use_density_property, const int k_pts, std::vector<double> propdensity, std::vector<double> propikp) {

		const clock_t begin_time = clock(); // start clock time

			params.scenename = karstic_network_name;
		// Cost due to distance

			if (is_simulation_parametrized) {

				// Compute 3D cost graph
				std::cout << "Step 1: Initializing cost graph..." << std::endl;
				GraphOperations graph;
				graph.InitializeCostGraphAdapted(create_nghb_graph, create_nghb_graph_property, keypts, params, nodes_on_densify_surface, inception_horizons, water_tables, use_sampling_points, box,
					max_inception_surface_distance, sampling_points, create_vset_sampling, use_density_property, k_pts, fraction_old_karst_perm, propdensity, propikp, topo_surface_);

				// Compute karstic skeleton
				std::cout << "Step 2: Computing karstic skeleton..." << std::endl;
				KarsticSkeleton skel = graph.ComputeKarsticSkeleton(keypts, fraction_karst_perm);

				// Procedural amplification
				if (use_amplification_phase_vadose) {
					std::cout << "Step 3a: Amplifying vadose zone..." << std::endl;
					skel.Amplify_vadose(&graph, params);
				}

				if (use_amplification_phase_phreatic) {
					std::cout << "Step 3b: Amplifying phreatic zone..." << std::endl;
					skel.Amplify_phreatic(&graph, params);
				}

				// Save network
				std::cout << "Saving network..." << std::endl;
				skel.CreateLine(params, karstic_network_name);
			}
	}

	void KarsticNetwork::set_save_directory(const std::string& repertory) {
		save_directory = repertory;
		params.directoryname = repertory;
	}

	void KarsticNetwork::set_amplification_vadose_params(const double max_dist_loops_vadose, const double loop_density_vadose) {
		params.max_dist_loops_vadose = max_dist_loops_vadose;
		params.loop_density_vadose = loop_density_vadose;
	}

	void KarsticNetwork::set_amplification_phreatic_params(const double max_dist_loops_phreatic, const double loop_density_phreatic) {
		params.max_dist_loops_phreatic = max_dist_loops_phreatic;
		params.loop_density_phreatic = loop_density_phreatic;
	}

	void KarsticNetwork::set_water_table_level(const double  water_table_constraint_weight_vadose, const double  water_table_constraint_weight_phreatic) {
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