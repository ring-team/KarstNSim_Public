/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/karstic_network.h"

namespace KarstNSim {

	KarsticNetwork::KarsticNetwork(const std::string& karstic_network_name, const Box& box, GeologicalParameters& params,
		const std::vector<KeyPoint>& keypts, const std::vector<Surface>& water_tables) :
		karstic_network_name(karstic_network_name), box(box), params(params), keypts(keypts), water_tables(water_tables) {};

	void KarsticNetwork::set_simulation_parameters(const int& nghb_count, const bool& use_max_nghb_radius, const float& nghb_radius, const float& poisson_radius,
		const float& gamma, const bool& multiply_costs, const bool& vadose_cohesion, const float& vertical_distance_stretching_factor) {

		params.gamma = gamma;
		params.graphNeighbourCount = nghb_count;
		params.graphuse_max_nghb_radius = use_max_nghb_radius;
		params.graphNeighbourRadius = nghb_radius;
		params.graphPoissonRadius = poisson_radius;
		params.distanceCost = CostTerm(true, 1);
		params.multiply_costs = multiply_costs;
		params.vadose_cohesion = vadose_cohesion;
		is_simulation_parametrized = true;
		params.vertical_distance_stretching_factor = std::max(1.0f,vertical_distance_stretching_factor);
		set_domain_geometry();
	}

	void KarsticNetwork::update_water_table_cost_channels()
	{
		const int nb_physical_wt = static_cast<int>(water_tables.size());

		params.nb_wt_surfaces = nb_physical_wt;

		if (params.has_springs_without_wt) {
			params.no_wt_cost_index = nb_physical_wt;
			params.nb_wt = nb_physical_wt + 1;
		}
		else {
			params.no_wt_cost_index = -1;
			params.nb_wt = nb_physical_wt;
		}
	}

	void KarsticNetwork::set_sinks(const std::vector<Vector3>& sinks, const std::vector<int>& propsinksindex,
		const std::vector<int>& propsinksorder, bool use_sinks_radius, const std::vector<float>& propsinksradius) {

		const int n_sinks = static_cast<int>(sinks.size());

		if (propsinksindex.size() != static_cast<size_t>(n_sinks)) {
			throw std::runtime_error(
				"[sinks] Invalid property 'sinks_index': expected " + std::to_string(n_sinks) +
				" values (one per sink), but got " + std::to_string(propsinksindex.size()) + "."
			);
		}

		if (propsinksorder.size() != static_cast<size_t>(n_sinks)) {
			throw std::runtime_error(
				"[sinks] Invalid property 'sinks_order': expected " + std::to_string(n_sinks) +
				" values (one per sink), but got " + std::to_string(propsinksorder.size()) + "."
			);
		}

		if (use_sinks_radius && propsinksradius.size() != static_cast<size_t>(n_sinks)) {
			throw std::runtime_error(
				"[sinks] Invalid property 'sinks_radius': option 'use_sinks_radius' is enabled, so "
				"expected " + std::to_string(n_sinks) + " values (one per sink), but got " +
				std::to_string(propsinksradius.size()) + "."
			);
		}

		// Validate sink indices: must be exactly a permutation of 1..n, with no duplicate and no hole.
		std::vector<int> seen_index(n_sinks + 1, 0);
		for (int i = 0; i < n_sinks; ++i) {
			const int idx = propsinksindex[i];

			if (idx < 1 || idx > n_sinks) {
				throw std::runtime_error(
					"[sinks] Invalid property 'sinks_index': value " + std::to_string(idx) +
					" found at sink #" + std::to_string(i + 1) +
					", but valid sink indices must be integers in [1, " + std::to_string(n_sinks) + "]."
				);
			}

			if (seen_index[idx] != 0) {
				throw std::runtime_error(
					"[sinks] Invalid property 'sinks_index': duplicate value " + std::to_string(idx) +
					" found for sink #" + std::to_string(i + 1) +
					" and sink #" + std::to_string(seen_index[idx]) +
					". Sink indices must define a permutation of 1.." + std::to_string(n_sinks) +
					" without repetition."
				);
			}

			seen_index[idx] = i + 1;
		}

		for (int idx = 1; idx <= n_sinks; ++idx) {
			if (seen_index[idx] == 0) {
				throw std::runtime_error(
					"[sinks] Invalid property 'sinks_index': missing value " + std::to_string(idx) +
					". Sink indices must define a complete permutation of 1.." + std::to_string(n_sinks) +
					" without hole."
				);
			}
		}

		// Validate sink orders: positive group identifiers starting at 1.
		int max_order = 0;
		for (int i = 0; i < n_sinks; ++i) {
			const int order = propsinksorder[i];

			if (order < 1) {
				throw std::runtime_error(
					"[sinks] Invalid property 'sinks_order': value " + std::to_string(order) +
					" found at sink #" + std::to_string(i + 1) +
					". Sink orders must be group identifiers starting at 1."
				);
			}

			if (order > max_order) {
				max_order = order;
			}
		}

		use_sinks_radius_ = use_sinks_radius;

		// Sinks are appended in a specific order:
		// 1) groups are processed in increasing order of 'sinks_order';
		// 2) within a same group, sinks are shuffled randomly.
		int scanning_order = 1;
		while (scanning_order <= max_order) {
			std::vector<int> propsinks_iter;
			std::vector<int> propsinks_iter_index;
			std::vector<int> indexes;

			for (int i = 0; i < n_sinks; i++) {
				if (propsinksorder[i] == scanning_order) {
					propsinks_iter.push_back(i);
					propsinks_iter_index.push_back(propsinksindex[i]);
				}
			}

			for (int i = 0; i < static_cast<int>(propsinks_iter.size()); ++i) {
				indexes.push_back(i);
			}

			shuffleContainer(indexes.begin(), indexes.end());
			reorder(propsinks_iter, indexes);
			reorder(propsinks_iter_index, indexes);

			for (int i = 0; i < static_cast<int>(propsinks_iter.size()); i++) {
				params.sinks_index.push_back(propsinks_iter_index[i]);
				this->keypts.emplace_back(sinks[propsinks_iter[i]], KeyPointType::Sink);
				this->pt_sink.push_back(sinks[propsinks_iter[i]]);

				if (use_sinks_radius_) {
					propsinksradius_.push_back({ propsinksradius[propsinks_iter[i]], int(keypts.size() - 1) });
				}
			}

			scanning_order++;
		}
	}

void KarsticNetwork::set_springs(
	const std::vector<Vector3>& springs,
	const std::vector<int>& propspringsindex,
	const bool& allow_single_outlet_connection,
	bool use_springs_radius,
	const std::vector<float>& propspringsradius,
	const std::vector<int>& propspringswtindex) {

	const int n_springs = static_cast<int>(springs.size());

	params.allow_single_outlet = allow_single_outlet_connection;
	params.nb_springs = n_springs;
	use_springs_radius_ = use_springs_radius;

	if (propspringsindex.size() != static_cast<size_t>(n_springs)) {
		throw std::runtime_error(
			"[springs] Invalid property 'springs_index': expected " + std::to_string(n_springs) +
			" values (one per spring), but got " + std::to_string(propspringsindex.size()) + "."
		);
	}

	if (propspringswtindex.size() != static_cast<size_t>(n_springs)) {
		throw std::runtime_error(
			"[springs] Invalid property 'springs_wt_index': expected " + std::to_string(n_springs) +
			" values (one per spring), but got " + std::to_string(propspringswtindex.size()) + "."
		);
	}

	if (use_springs_radius && propspringsradius.size() != static_cast<size_t>(n_springs)) {
		throw std::runtime_error(
			"[springs] Invalid property 'springs_radius': option 'use_springs_radius' is enabled, so "
			"expected " + std::to_string(n_springs) + " values (one per spring), but got " +
			std::to_string(propspringsradius.size()) + "."
		);
	}

	const int n_wt = static_cast<int>(water_tables.size());

	// Validate spring indices: they must define a complete permutation of 1..n_springs.
	std::vector<int> seen_index(n_springs + 1, 0);
	for (int i = 0; i < n_springs; ++i) {
		const int idx = propspringsindex[i];

		if (idx < 1 || idx > n_springs) {
			throw std::runtime_error(
				"[springs] Invalid property 'springs_index': value " + std::to_string(idx) +
				" found at spring #" + std::to_string(i + 1) +
				", but valid spring indices must be integers in [1, " + std::to_string(n_springs) + "]."
			);
		}

		if (seen_index[idx] != 0) {
			throw std::runtime_error(
				"[springs] Invalid property 'springs_index': duplicate value " + std::to_string(idx) +
				" found for spring #" + std::to_string(i + 1) +
				" and spring #" + std::to_string(seen_index[idx]) +
				". Spring indices must define a permutation of 1.." + std::to_string(n_springs) +
				" without repetition."
			);
		}

		seen_index[idx] = i + 1;
	}

	for (int idx = 1; idx <= n_springs; ++idx) {
		if (seen_index[idx] == 0) {
			throw std::runtime_error(
				"[springs] Invalid property 'springs_index': missing value " + std::to_string(idx) +
				". Spring indices must define a complete permutation of 1.." + std::to_string(n_springs) +
				" without hole."
			);
		}
	}

	// Validate water-table index associated with each spring.
	// A value of 0 is valid and means that the spring has no associated water table.
	params.has_springs_without_wt = false;

	for (int i = 0; i < n_springs; ++i) {
		const int wt_idx = propspringswtindex[i];

		if (wt_idx == 0) {
			params.has_springs_without_wt = true;
			continue;
		}

		if (wt_idx < 0 || wt_idx > n_wt) {
			throw std::runtime_error(
				"[springs] Invalid property 'springs_wt_index': value " + std::to_string(wt_idx) +
				" found at spring #" + std::to_string(i + 1) +
				". Valid values are 0 for no associated water table, or integers in [1, " +
				std::to_string(n_wt) + "] for springs associated with a water table."
			);
		}
	}

	update_water_table_cost_channels();

	// Append springs in connectivity-matrix order, i.e. index 1 first, then 2, ..., then n.
	for (int i = 0; i < n_springs; ++i) {
		const auto it = std::find(propspringsindex.begin(), propspringsindex.end(), i + 1);

		// This should never happen after validation, but keep a defensive check.
		if (it == propspringsindex.end()) {
			throw std::runtime_error(
				"[springs] Internal error while ordering springs: index " + std::to_string(i + 1) +
				" was validated but could not be found."
			);
		}

		const int index_for_i = static_cast<int>(std::distance(propspringsindex.begin(), it));

		this->keypts.emplace_back(
			springs.at(index_for_i),
			KeyPointType::Spring,
			propspringswtindex[index_for_i]
		);

		params.z_list.push_back({
			springs.at(index_for_i).z,
			int(keypts.size() - 1)
		});

		this->pt_spring.push_back(springs.at(index_for_i));

		params.propspringswtindex.push_back({
			static_cast<float>(propspringswtindex[index_for_i]),
			int(keypts.size() - 1)
		});

		if (use_springs_radius_) {
			propspringsradius_.push_back({
				propspringsradius[index_for_i],
				int(keypts.size() - 1)
			});
		}
	}
}

	void KarsticNetwork::set_sinks_sections_only(const std::vector<Vector3>& sinks, bool use_sinks_radius, const std::vector<float>& propsinksradius) {

		const int n_sinks = static_cast<int>(sinks.size());

		use_sinks_radius_ = use_sinks_radius;

		if (use_sinks_radius && propsinksradius.size() != static_cast<size_t>(n_sinks)) {
			throw std::runtime_error(
				"[sinks] Invalid property 'sinks_radius': option 'use_sinks_radius' is enabled, so "
				"expected " + std::to_string(n_sinks) + " values (one per sink), but got " +
				std::to_string(propsinksradius.size()) + "."
			);
		}

		for (int i = 0; i < n_sinks; ++i) {
			this->keypts.emplace_back(sinks[i], KeyPointType::Sink);
			this->pt_sink.push_back(sinks[i]);

			if (use_sinks_radius_) {
				propsinksradius_.push_back({ propsinksradius[i], int(keypts.size() - 1) });
			}
		}
	}

void KarsticNetwork::set_springs_sections_only(
	const std::vector<Vector3>& springs,
	bool use_springs_radius,
	const std::vector<float>& propspringsradius) {

	const int n_springs = static_cast<int>(springs.size());

	params.nb_springs = n_springs;
	params.has_springs_without_wt = (n_springs > 0);
	use_springs_radius_ = use_springs_radius;

	update_water_table_cost_channels();

	if (use_springs_radius && propspringsradius.size() != static_cast<size_t>(n_springs)) {
		throw std::runtime_error(
			"[springs] Invalid property 'springs_radius': option 'use_springs_radius' is enabled, so "
			"expected " + std::to_string(n_springs) + " values (one per spring), but got " +
			std::to_string(propspringsradius.size()) + "."
		);
	}

	for (int i = 0; i < n_springs; ++i) {
		this->keypts.emplace_back(springs.at(i), KeyPointType::Spring, 0);
		params.z_list.push_back({ springs.at(i).z, int(keypts.size() - 1) });
		this->pt_spring.push_back(springs.at(i));
		params.propspringswtindex.push_back({ 0.0f, int(keypts.size() - 1) });

		if (use_springs_radius_) {
			propspringsradius_.push_back({ propspringsradius[i], int(keypts.size() - 1) });
		}
	}
}


	void KarsticNetwork::set_waypoints(const std::vector<Vector3>& waypoints,
		bool use_waypoints_radius,
		const std::vector<float>& propwaypointsradius,
		const std::vector<float>& propwaypointsimpactradius,
		float waypoints_weight)
	{
		const int n_pts = static_cast<int>(waypoints.size());

		// Mandatory: impact radius size must match
		if (static_cast<int>(propwaypointsimpactradius.size()) != n_pts) {
			std::cout << "[waypoints][error] impact_radius size mismatch: points=" << n_pts
				<< " impact_radius=" << propwaypointsimpactradius.size() << std::endl;
			throw std::runtime_error("Missing/invalid 'impact_radius' for waypoints: size mismatch.");
		}

		// Optional: per-waypoint radius only if sizes match
		use_waypoints_radius_ = use_waypoints_radius
			&& (static_cast<int>(propwaypointsradius.size()) == n_pts);
		if (use_waypoints_radius && !use_waypoints_radius_) {
			std::cout << "[waypoints][warn] per-waypoint 'radius' ignored (size mismatch)." << std::endl;
		}

		params.waypoints_weight = waypoints_weight;

		for (int i = 0; i < n_pts; ++i) {
			this->keypts.emplace_back(waypoints[i], KeyPointType::Waypoint);
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

	void KarsticNetwork::set_previous_networks(const std::vector<Line>& previous_networks) {

		int count_total_nb_segs = 0;
		for (int i = 0; i < previous_networks.size(); i++) {
			Line Line_i = previous_networks.at(i);
			count_total_nb_segs += int(Line_i.get_nb_segs());
		}
		params.PtsOldGraph.resize(count_total_nb_segs, 2); // vector of all segments of the old graph, represented each with their start and end point
		params.IdxOldGraph.resize(count_total_nb_segs, 2); // vector of all point of the old graph, each represented by its index in the samples object of GraphOperations
		count_total_nb_segs = 0;
		for (int i = 0; i < previous_networks.size(); i++) {
			Line Line_i = previous_networks.at(i);
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

	void KarsticNetwork::set_input_nghb_graph(const InputGraph& input_nghb_graph) {
		params.use_input_nghb_graph = true;
		params.input_nghb_graph = input_nghb_graph;
	}

	void KarsticNetwork::set_deadend_points(int nb_deadend_points, float max_distance_of_deadend_pts) {
		use_deadend_pts_ = true;
		nb_deadend_points_ = nb_deadend_points;
		max_distance_of_deadend_pts_ = max_distance_of_deadend_pts;
	}

	void KarsticNetwork::set_inception_surfaces_sampling(const std::string& network_name, const std::vector<Surface>& surfaces_used_to_densify, const int& refine_surface_sampling, const bool& create_vset_sampling) {
		params.nb_inception_surf = int(surfaces_used_to_densify.size());
		KarstNSim::surface_sampling::multiple_surface_sampling(params.directoryname, network_name, box, surfaces_used_to_densify, nodes_on_inception_surfaces, refine_surface_sampling, create_vset_sampling);
	}

	void KarsticNetwork::set_wt_surfaces_sampling(const std::string& network_name, const std::vector<Surface>& surfaces_used_to_densify, const int& refine_surface_sampling) {

		const int nb_physical_wt = static_cast<int>(surfaces_used_to_densify.size());
		nodes_on_wt_surfaces.resize(nb_physical_wt);
		update_water_table_cost_channels();

		for (int i = 0; i < surfaces_used_to_densify.size(); i++) {
			std::vector<Surface> wt_surface;
			wt_surface.push_back(surfaces_used_to_densify.at(i));
			std::vector<Vector3> nodes_on_wt_surface;
			KarstNSim::surface_sampling::multiple_surface_sampling(params.directoryname, network_name, box, wt_surface, nodes_on_wt_surface, refine_surface_sampling, false);
			nodes_on_wt_surfaces[i] = nodes_on_wt_surface;
		}
	}

	void KarsticNetwork::save_painted_box(
		const std::vector<float>& propdensity,
		const std::vector<float>& propikp)
	{

		const size_t expected_cell_count =
			static_cast<size_t>(box.get_nu()) *
			static_cast<size_t>(box.get_nv()) *
			static_cast<size_t>(box.get_nw());

		if (!propdensity.empty() && propdensity.size() != expected_cell_count) {
			throw std::runtime_error(
				"[save_painted_box] Density property size does not match the domain box cell count."
			);
		}

		if (!propikp.empty() && propikp.size() != expected_cell_count) {
			throw std::runtime_error(
				"[save_painted_box] IKP property size does not match the domain box cell count."
			);
		}

		std::vector<std::string> property_names = { "density", "karstif_potential" };
		std::vector<std::vector<float>> properties(
			expected_cell_count,
			std::vector<float>(2, -99999.0f)
		);

		for (size_t i = 0; i < expected_cell_count; ++i) {
			if (!propdensity.empty()) {
				properties[i][0] = propdensity[i];
			}

			if (!propikp.empty()) {
				properties[i][1] = propikp[i];
			}
		}

		std::string full_name = karstic_network_name + "_box.txt";
		const std::string full_dir_name = params.directoryname + "/outputs";
		save_box(full_name, full_dir_name, box, property_names, properties);
	}



	void KarsticNetwork::set_topo_surface(const Surface& topo_surface) {
		this->topo_surface_ = topo_surface;
	}

	void KarsticNetwork::set_inception_horizons_parameters(const std::vector<Surface>& inception_horizons_list, const float& inception_horizon_constraint_weight) {
		this->inception_horizons = inception_horizons_list;
		params.horizonCost = CostTerm(true, inception_horizon_constraint_weight);
	}

	void KarsticNetwork::set_ghost_rocks(const Box& grid, std::vector<float>& ikp, const Line& alteration_lines, const bool& interpolate_lines, const float& ghostrock_max_vertical_size, const bool& use_max_depth_constraint, const float& ghost_rock_weight, Surface* max_depth_horizon, const float& ghostrock_width) {
		(void)interpolate_lines;
		params.use_ghost_rocks = true;
		params.length = ghostrock_max_vertical_size;
		params.width = ghostrock_width;
		params.polyline = alteration_lines;
		params.use_max_depth_constraint = use_max_depth_constraint;
		params.substratum_surf = *max_depth_horizon;

		// modify "ikp" object with ghostrocks ("paint" it). Note that density property is NOT changed, hence passed as const ref
		paint_KP_with_ghostrocks(grid, ikp, ghostrock_max_vertical_size, ghostrock_width, alteration_lines, use_max_depth_constraint, *max_depth_horizon, ghost_rock_weight);
	}

	void KarsticNetwork::disable_inception_horizon() {
		params.horizonCost = CostTerm(false, 0.0);
	}

	void KarsticNetwork::set_karstification_potential_parameters(const float& karstification_potential_weight) {
		params.karstificationCost = CostTerm(true, karstification_potential_weight);
	}

	void KarsticNetwork::set_fracture_constraint_parameters(const std::vector<float>& fracture_families_orientations, const std::vector<float>& fracture_families_tolerance, const float& fracture_constraint_weight) {
		params.fractures_orientations = fracture_families_orientations;
		params.fractures_tolerances = fracture_families_tolerance;
		params.fractureCost = CostTerm(true, fracture_constraint_weight);
	}

	void KarsticNetwork::disable_fractures() {
		params.fractureCost = CostTerm(false, 0.0);
	}

	void KarsticNetwork::set_no_karst_spheres_parameters(const std::vector<Vector3>& sphere_centers,
		const std::vector<float>& sphere_radius)
	{
		// --- Defensive checks to prevent crashes and undefined behavior ---
		if (sphere_centers.empty()) {
			std::cout << "WARNING: no 'no-karst' spheres provided; skipping spheres setup." << std::endl;
			return;
		}

		if (sphere_radius.empty()) {
			return;
		}

		const bool one_radius_for_all = (sphere_radius.size() == 1);
		const bool per_sphere_radius = (sphere_radius.size() == sphere_centers.size());

		if (!one_radius_for_all && !per_sphere_radius) {
			std::cout << "WARNING: sphere_radius size (" << sphere_radius.size()
				<< ") does not match 1 or number of centers (" << sphere_centers.size()
				<< "). Using the first radius for all spheres as fallback." << std::endl;
		}

		int added = 0;
		for (int i = 0; i < static_cast<int>(sphere_centers.size()); ++i) {
			float r = 0.0f;

			if (one_radius_for_all) {
				r = sphere_radius[0];
			}
			else if (per_sphere_radius) {
				r = sphere_radius[i];
			}
			else {
				// Fallback: use first radius for all (already warned above).
				r = sphere_radius[0];
			}

			if (r <= 0.0f || !std::isfinite(r)) {
				std::cout << "WARNING: invalid sphere radius at index " << i
					<< " (r=" << r << "); skipping this sphere." << std::endl;
				continue;
			}

			params.spheres.push_back(Sphere(sphere_centers.at(i), r));
			++added;
		}
	}


	void KarsticNetwork::set_domain_geometry() {

		const Vector3& min_pt = box.get_basis();
		const Vector3& max_pt = box.get_end();

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

	void KarsticNetwork::read_connectivity_matrix(const std::string& simulation_input_dir, const std::vector<Vector3>& sinks, const std::vector<Vector3>& springs) {

		const std::string connectivity_matrix_path = simulation_input_dir + "/connectivity_matrix.txt";

		int nb_sinks = static_cast<int>(sinks.size());
		int nb_springs = static_cast<int>(springs.size());

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

		//if (geostatparams.is_used) {
		//	//skel.compute_distance_matrix(); // compute distance matrix between each pair of nodes in the graph
		//}

		// Compute z_phreatic as the lowest z among all springs
		float z_phreatic = std::numeric_limits<float>::max();
		for (const auto& zp : params.z_list) {
			if (zp.prop < z_phreatic) {
				z_phreatic = zp.prop;
			}
		}

		// Outputs linked to the external drift computation. Unused if no drift is applied.
		std::vector<float> drift_output(skel.nodes.size(), -99999.f);
		std::vector<float> weights_output(skel.nodes.size(), -99999.f);

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
			SGS3_with_external_drift(
				&skel,
				this->pt_spring,
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
				geostatparams.proportion_interbranch,
				z_phreatic,
				geostatparams.use_drift_zwt,
				geostatparams.use_drift_curv,
				drift_output,
				weights_output
			);
			if (geostatparams.use_drift_zwt) params.use_drift_zwt = true;
			if (geostatparams.use_drift_curv) params.use_drift_curv = true;

			// assign property to skeleton :

			for (int i = 0; i < geostatparams.simulated_property.size(); i++) {
				skel.nodes.at(i).eq_radius = geostatparams.simulated_property.at(i);
				if (geostatparams.use_drift_zwt || geostatparams.use_drift_curv) {

					skel.nodes.at(i).drift_value = drift_output[i];
					skel.nodes.at(i).drift_weight = weights_output[i];
				}
			}
		}

		//Finally we paint any node that is within a ghostrock to have the ghost rock width

		if (params.use_ghost_rocks) {
			paint_karst_sections_with_ghostrocks(skel, params.length, params.width, params.polyline, params.use_max_depth_constraint, params.substratum_surf);
		}
	}

	std::optional<KarstNetworkResult> KarsticNetwork::run_simulation(
		const bool &sections_simulation_only, const bool &create_nghb_graph,
		const bool &create_nghb_graph_property,
		const bool &create_solved_connectivity_matrix,
		const bool &use_amplification, const bool &use_sampling_points,
		const float &fraction_karst_perm, const float &fraction_old_karst_perm,
		const float &max_inception_surface_distance,
		const std::vector<Vector3> &sampling_points, const bool &create_vset_sampling,
		const bool &use_density_property, const int &k_pts,
		const std::vector<float> &propdensity, const std::vector<float> &propikp
	) {
		params.scenename = karstic_network_name;
		const clock_t time1 = clock();

		if (is_simulation_parametrized) { // full simulation

			// Compute 3D cost graph

			std::cout << "\nSTEP 1 - Generation of cost graph:\n\n";

			GraphOperations graph;
			graph.InitializeCostGraph(
				create_nghb_graph, create_nghb_graph_property, keypts, params,
				nodes_on_inception_surfaces, nodes_on_wt_surfaces, inception_horizons,
				water_tables, use_sampling_points, box, max_inception_surface_distance,
				sampling_points, create_vset_sampling, use_density_property, k_pts,
				fraction_old_karst_perm, propdensity, propikp, topo_surface_);
			const clock_t time2 = clock();

			std::cout << "-> STEP 1 completed: Cost graph generated (" << std::fixed
					<< std::setprecision(3) << float(time2 - time1) / CLOCKS_PER_SEC
					<< " s)\n";

			std::cout << "\nSTEP 2 - Simulation of the karst network skeleton: \n\n";

			// Compute karstic skeleton
			std::vector<std::vector<int>> karst_paths;
			std::vector<std::vector<float>> karst_paths_costs;
			std::vector<std::vector<char>> karst_paths_vadose_flag;
			std::vector<int> springidxFinal;

			graph.ComputeKarsticSkeleton(keypts, fraction_karst_perm, karst_paths,
										karst_paths_costs, karst_paths_vadose_flag,
										springidxFinal,
										create_solved_connectivity_matrix);

			if (karst_paths.empty()) {
				std::cout << "No path found between inlets and outlets with the "
							"current parameters."
							<< std::endl;
				return std::nullopt; // no path found
			}

			// Build karstic skeleton structure
			KarsticSkeleton skel(&graph, karst_paths, karst_paths_costs,
								karst_paths_vadose_flag, springidxFinal);

			const clock_t time3 = clock();

			std::cout << "-> STEP 2 completed: Skeleton computed (" << std::fixed << std::setprecision(3)
					<< float(time3 - time2) / CLOCKS_PER_SEC << " s)" << std::endl;
			std::cout << "\nSTEP 3 - Network amplification (cycles and deadends): \n\n";

			// Network preparation
			skel.detect_intersection_points(karst_paths);
			skel.update_branch_ID(karst_paths);

			// Procedural amplification with deadend nodes
			const clock_t time4 = clock();
			if (use_deadend_pts_) {
				skel.amplify_deadend(&graph, max_distance_of_deadend_pts_, nb_deadend_points_, params);
			}
			const clock_t time5 = clock();
			if (use_deadend_pts_) {
				std::cout << "Network amplified with deadend points (" << std::fixed
							<< std::setprecision(3) << float(time5 - time4) / CLOCKS_PER_SEC
							<< " s)" << std::endl;
			}

			// Amplification
			if (use_amplification) {
				if (params.use_noise && !params.use_noise_on_all) {
					graph.add_noise();
				}
				std::pair<float, float> result = skel.amplify_noise(&graph, params);
			}
			const clock_t time5bis = clock();
			if (use_amplification) {

					std::cout << "Network amplified with cycles (" << std::fixed
					<< std::setprecision(3) << float(time5bis - time5) / CLOCKS_PER_SEC
					<< " s)" << std::endl;
			}
						
			if (!(use_amplification) && !(use_deadend_pts_)) {
				std::cout << "-> STEP 3 skipped (no amplification)\n";
			} else {
				std::cout << "-> STEP 3 completed: Network amplified (" << std::fixed << std::setprecision(3)
					<< float(time5bis - time3) / CLOCKS_PER_SEC << " s)" << std::endl; 
			}

			const clock_t time6 = clock();
			std::cout << "\nSTEP 4 - Simulation of conduit sections: \n\n";

			skel.refresh_vadose_flags_from_graph(&graph, params);

			skel.prepare_graph(); // removes duplicates, and changes format so that each node is connected to all of its neighbors (not the case at base necesarily)
			skel.update_branch_ID(); // create branch id property on the skeleton nodes (-1 if intersection, branch index >=1 otherwise)
			skel.compute_branch_sizes(); // compute the number of nodes in each branch
			skel.compute_valence(); // compute the valence = number of neighbors of each node

			create_sections(skel);
			const clock_t time7 = clock();
			if (geostatparams.is_used) {
				std::cout << "-> STEP 4 completed: Conduit sections generated (" << std::fixed
					<< std::setprecision(3) << float(time7 - time6) / CLOCKS_PER_SEC
					<< " s)" << std::endl;
			} else {
				std::cout << "-> STEP 4 skipped (no section simulation required)\n";
			}
			// save network
			auto res = skel.create_line(params, karstic_network_name);
			const clock_t time8 = clock();

			std::cout << "\nKarst network saved (" << std::fixed << std::setprecision(3)
					<< float(time8 - time7) / CLOCKS_PER_SEC << " s)" << std::endl;

			return res;
		} else if (sections_simulation_only) { // only generate sections
		
			std::cout << "\nSimulation of properties only..." << std::endl;
		
			std::vector<std::vector<float>> costs_graph(params.PtsOldGraph.size(), std::vector<float>(2, 0.0)); // dummy vector
			std::vector<std::vector<char>> vadoseflags_graph(params.PtsOldGraph.size(), std::vector<char>(2, false)); // dummy vector
			std::vector<int> springidx(params.PtsOldGraph.size(), 1);
			KarsticSkeleton skel(params.PtsOldGraph, costs_graph, vadoseflags_graph, springidx);

			// Build a clean graph topology and compute branch labels even in sections-only mode,
			// so that exported branch IDs are meaningful.
			skel.prepare_graph();
			skel.update_branch_ID();
			skel.compute_branch_sizes();
			skel.compute_valence();

			create_sections(skel);
			
			const clock_t time2 = clock();
			std::cout << "Conduit sections simulated (" << std::fixed << std::setprecision(3)
					<< float(time2 - time1) / CLOCKS_PER_SEC << " s)" << std::endl << std::endl;

			return skel.create_line(params, karstic_network_name);
		}
		return std::nullopt;
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

	void KarsticNetwork::set_gradient_constraint_weight(const float& gradient_constraint_weight) {
		params.gradient_constraint_weight = std::max(0.0f, gradient_constraint_weight);
	}

	void KarsticNetwork::set_outlet_selection_cost_factor(const float& outlet_selection_cost_factor) {
		params.outlet_selection_cost_factor = std::max(1.0f, outlet_selection_cost_factor);
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