/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for new methods and modifications to original methods
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for original versions of KarsticSkeleton constructor, amplify and append_paths
If you use this code, please cite : Paris et al., 2021, Computer Graphic Forum.

***************************************************************/

#include "KarstNSim/graph.h"
#include "KarstNSim/read_files.h"

namespace {
	int select_lowest_spring_cost_channel_1based(const KarstNSim::GeologicalParameters& params)
	{
		if (params.z_list.empty()) {
			throw std::runtime_error(
				"[skeleton] Cannot select a reference spring: params.z_list is empty."
			);
		}

		auto min_it = std::min_element(
			params.z_list.begin(),
			params.z_list.end(),
			[](const KarstNSim::GeologicalParameters::Propidx& a,
				const KarstNSim::GeologicalParameters::Propidx& b) {
			return a.prop < b.prop;
		}
		);

		const int lowest_spring_keypoint_index = min_it->index;

		auto it = std::find_if(
			params.propspringswtindex.begin(),
			params.propspringswtindex.end(),
			[lowest_spring_keypoint_index](const KarstNSim::GeologicalParameters::Propidx& item) {
			return item.index == lowest_spring_keypoint_index;
		}
		);

		if (it == params.propspringswtindex.end()) {
			throw std::runtime_error(
				"[skeleton] Cannot find the water-table association of the lowest spring."
			);
		}

		const int wt_idx = static_cast<int>(it->prop);

		if (wt_idx > 0) {
			return wt_idx;
		}

		if (params.no_wt_cost_index >= 0) {
			return params.no_wt_cost_index + 1;
		}

		throw std::runtime_error(
			"[skeleton] The lowest spring has no associated water table, but the vadose-only cost channel was not initialized."
		);
	}
}

namespace KarstNSim {
	KarsticSkeleton::KarsticSkeleton(const GraphOperations* graph, const std::vector<std::vector<int>>& paths, const std::vector<std::vector<float>> costs, const  std::vector<std::vector<char>> vadose, const std::vector<int>& springidx)
	{
		append_paths(graph, paths, costs, vadose, springidx);
	}

	KarsticSkeleton::KarsticSkeleton(const Array2D<Vector3>& pts_graph, const std::vector<std::vector<float>>& costs_graph, const  std::vector<std::vector<char>>& vadoseflags_graph, const std::vector<int>& springidx) {

		int nb_springs = *std::max_element(springidx.begin(), springidx.end()) + 1;
		int unique_pts = 0;
		for (int i = 0; i < pts_graph.size(); i++)	// Paths
		{
			int path_size = int(pts_graph.row(i).size()); // path size

			for (int j = path_size - 1; j >= 0; j--)
			{
				Vector3 node = pts_graph(i, j);
				std::vector<float> nodeCost(nb_springs, 0.);
				nodeCost[springidx[i]] = 0.; // DEFAULT VALUE BECAUSE NOT COMPUTED BEFORE costs_graph[i][j];
				std::vector<int> vadoseTest(nb_springs, -1);
				vadoseTest[springidx[i]] = -1; // DEFAULT VALUE BECAUSE NOT COMPUTED BEFORE vadoseflags_graph[i][j];
				auto it = std::find(nodes.begin(), nodes.end(), node);
				if (it == nodes.end())
				{
					nodes.push_back({ unique_pts, node ,nodeCost, vadoseTest });
					unique_pts++;
				}
				else {
					int index = std::distance(nodes.begin(), it);
					nodes[index].vadose[springidx[i]] = -1; // DEFAULT VALUE BECAUSE NOT COMPUTED BEFORE vadoseflags_graph[i][j];
					nodes[index].cost[springidx[i]] = 0.; // DEFAULT VALUE BECAUSE NOT COMPUTED BEFORE costs_graph[i][j];
				}
			}
		}
		for (int i = 0; i < nodes.size(); i++)	// Edges
		{
			Vector3 node = nodes[i].p;
			std::vector<int> edges;
			for (int j = 0; j < pts_graph.size(); j++)
			{
				for (int k = int(pts_graph.row(j).size()) - 1; k >= 0; k--)
				{
					Vector3 n = pts_graph(j, k);
					if (node != n)
						continue;
					if (k > 0)
					{
						Vector3 nafter = pts_graph(j, k - 1);
						int neiIndex = get_internal_index(nafter);
						if (std::find(edges.begin(), edges.end(), neiIndex) == edges.end())
							edges.push_back(neiIndex);
					}
				}
			}
			for (int j = 0; j < edges.size(); j++)
			{
				auto new_connection = KarsticConnection({ edges[j] });
				if (std::find(nodes[i].connections.begin(), nodes[i].connections.end(), edges[j]) == nodes[i].connections.end())
					nodes[i].connections.push_back(new_connection);
			}
		}
	}

	// number of vadose nodes for a given spring
	int KarsticSkeleton::count_vadose_nodes(int spring)
	{
		int count = 0;
		for (int i = 0; i < nodes.size(); i++) {
			if (nodes[i].vadose[spring] == 1) {
				count += 1;
			}
		}
		return count;
	}

	// number of vadose nodes on average for all springs (is equal to vadose nodes of spring if only one spring obviously)
	int KarsticSkeleton::count_average_vadose_nodes() {
		int totalNodes = int(nodes.size());
		int totalVadose = 0;
		int totalSprings = 0; // Assuming all nodes have the same number of springs

		if (totalNodes == 0) return 0; // Avoid division by zero

		// Calculate total vadose nodes over all springs
		for (const auto& node : nodes) {
			totalVadose += std::count(node.vadose.begin(), node.vadose.end(), 1);
			if (!node.vadose.empty()) totalSprings = int(node.vadose.size());
		}

		// Calculate the average count of vadose nodes
		int averageVadose = totalVadose / totalSprings;
		return averageVadose;
	}

	bool KarsticSkeleton::edge_exists(int u, int v) {
		return std::find(nodes[u].connections.begin(), nodes[u].connections.end(), v) != nodes[u].connections.end();
	}

	bool KarsticSkeleton::idx_is_prefix(int id1, int id2) const {
		if (id2 > id1) {
			std::swap(id1, id2);
		}

		while (id1 > 0) {
			if (id1 == id2) {
				return true;
			}
			id1 /= 10;
		}
		return false;
	}

	bool KarsticSkeleton::is_parent(int u, int v) const {
		if (u < 0 || v < 0 || u >= nodes.size() || v >= nodes.size()) {
			return false;
		}
		std::vector<int> branch_id1 = nodes[u].branch_id_ascend;
		std::vector<int> branch_id2 = nodes[v].branch_id_ascend;

		for (auto& id1 : branch_id1) {
			for (auto& id2 : branch_id2) {
				if (idx_is_prefix(id1, id2)) {
					return true;
				}
			}
		}
		return false;
	}

	std::string vectorToString(const std::vector<int>& vec, const std::string& delimiter = " \n ") {
		std::string result;
		for (size_t i = 0; i < vec.size(); ++i) {
			// Append the integer converted to string
			result += std::to_string(vec[i]);

			// Append the delimiter except after the last element
			if (i < vec.size() - 1) {
				result += delimiter;
			}
		}
		return result;
	}

	bool KarsticSkeleton::always_in_vadose(int u) {
		int count = 0;
		for (int i = 0; i < nodes.size(); i++) {
			count += std::count(nodes[i].vadose.begin(), nodes[i].vadose.end(), 1);
		}
		if (count == nodes[u].vadose.size()) {
			return true;
		}
		else {
			return false;
		}
	}

	std::pair<float, float> KarsticSkeleton::amplify_noise(GraphOperations* graph, const GeologicalParameters& params)
	{

		std::vector<std::vector<int>> newPaths;
		std::vector<std::vector<float>> newPaths_cost;
		std::vector<std::vector<char>> newPaths_vadose;
		std::vector<int> springidx;
		float mean_cycles_length = 0.;
		float mean_cycles_length_inter_wt = 0.;
		int cycles_created = 0;

		const int idx_specific_wt = select_lowest_spring_cost_channel_1based(params);

		// --- PARAMETER GUARDS: amplification distance window ---
		if (params.min_distance_amplification > params.max_distance_amplification) {
			throw std::runtime_error(
				"Invalid amplification distance window: min_distance_amplification > max_distance_amplification. "
				"Please fix your instruction file."
			);
		}
		if (params.max_distance_amplification <= 0.f) {
			throw std::runtime_error(
				"Invalid amplification distance window: max_distance_amplification must be > 0."
			);
		}
		if (params.min_distance_amplification < 0.f) {
			throw std::runtime_error(
				"Invalid amplification distance window: min_distance_amplification must be >= 0."
			);
		}

		const std::size_t node_count = nodes.size();

		if (params.nb_cycles > 0 && node_count < 2u) {
			throw std::runtime_error(
				"[amplification] Cycle amplification requires at least two skeleton nodes."
			);
		}

		// Ordered pairs are used because both graph connectivity and costs are directional.
		const std::size_t ordered_pair_count =
			node_count >= 2u
			? node_count * (node_count - 1u)
			: 0u;

		std::size_t tested_pairs_for_current_cycle = 0u;
		std::size_t pair_cursor =
			params.nb_cycles > 0
			? generateRandomIndex(ordered_pair_count)
			: 0u;

		while (cycles_created < params.nb_cycles)
		{
			if (tested_pairs_for_current_cycle >= ordered_pair_count) {
				throw std::runtime_error(
					"[amplification] Unable to create the requested cycle: every ordered "
					"node pair was tested, but none satisfied the distance, connectivity, "
					"and shortest-path requirements."
				);
			}

			const std::size_t pair_index = pair_cursor;
			pair_cursor = (pair_cursor + 1u) % ordered_pair_count;
			++tested_pairs_for_current_cycle;

			const std::size_t u_index = pair_index / (node_count - 1u);
			const std::size_t compressed_v_index =
				pair_index % (node_count - 1u);

			// Restore the skipped diagonal element so that u != v.
			const std::size_t v_index =
				compressed_v_index >= u_index
				? compressed_v_index + 1u
				: compressed_v_index;

			const int u = static_cast<int>(u_index);
			const int v = static_cast<int>(v_index);

			// Preserve the original shortest-path orientation.
			const int int_u_index = nodes[v].index;
			const int int_v_index = nodes[u].index;

			// Check if edge already exists or u and v are the same vertex
			if (u != v && !edge_exists(u, v)) {
				const float distance = graph_anisotropic_distance(
					nodes[u].p,
					nodes[v].p,
					params.vertical_distance_stretching_factor
				);

				if (distance >= params.min_distance_amplification && distance <= params.max_distance_amplification) { // also, check if the euclidian distance separating them is less than the threshold value defined by the user

					std::pair< std::vector<int>, std::vector<float>> pair;

					if (is_parent(u, v)) {
						const float radius = magnitude(nodes[u].p - nodes[v].p) / 4.0f;
						const Vector3 center =
							nodes[u].p + (nodes[v].p - nodes[u].p) / 2.0f;

						const std::vector<GraphOperations::EdgeWeightBackup> reset_values =
							graph->add_ball_cost(center, radius);

						const auto restore_ball_costs = [&]() {
							for (const GraphOperations::EdgeWeightBackup& backup : reset_values) {
								graph->adj[backup.source_node_index]
									[backup.edge_slot]
									.weight[backup.cost_channel] = backup.previous_weight;
							}
						};

						try {
							pair = graph->shortest_path(
								int_u_index,
								int_v_index,
								idx_specific_wt
							);
						}
						catch (...) { // catches all exceptions
							// The temporary perturbation must also be reverted if the
							// shortest-path computation raises an exception.
							restore_ball_costs();
							throw;
						}

						restore_ball_costs();
					}
					else {
						pair = graph->shortest_path(int_u_index, int_v_index, idx_specific_wt); // if all the previous tests were ok, then we can compute a path between those two nodes
					}
					std::vector<int> bestPath = pair.first;
					std::vector<float> bestPath_cost = pair.second;

					// An unreachable target is represented by an incomplete path. Such a path
					// must not be counted as a generated cycle.
					if (bestPath.size() < 2u) {
						continue;
					}

					if (bestPath_cost.size() != bestPath.size()) {
						throw std::runtime_error(
							"[amplification] Shortest-path output is inconsistent: the node and "
							"cost vectors do not have the same size."
						);
					}

					bool inter_wt = false;
					if (nodes[0].vadose.size() == 2) {
						if (((nodes[u].vadose[1] == 0 && nodes[v].vadose[1] == 0) && (nodes[u].vadose[0] == 1 && nodes[v].vadose[0] == 1)) ||
							((nodes[v].vadose[1] == 0 && nodes[u].vadose[1] == 0) && (nodes[v].vadose[0] == 1 && nodes[u].vadose[0] == 1))) {
							inter_wt = true;
						}
					}

					for (std::size_t i = 0; i + 1u < bestPath.size(); ++i) {
						mean_cycles_length += magnitude(nodes[i].p - nodes[i + 1].p);
						if (inter_wt) {
							mean_cycles_length += magnitude(nodes[i].p - nodes[i + 1].p);
						}

					}

					newPaths.push_back(std::vector<int>(bestPath.begin(), bestPath.end()));
					newPaths_cost.push_back(std::vector<float>(bestPath_cost.begin(), bestPath_cost.end()));
					newPaths_vadose.push_back(std::vector<char>(bestPath.size(), true)); // all paths are vadose paths here
					springidx.push_back(0);
					++cycles_created; // and then incerment the number of loops still to be generated
					tested_pairs_for_current_cycle = 0u;
					pair_cursor = generateRandomIndex(ordered_pair_count);
				}


			}
		}
		// Update the internal karstic skeleton struture
		if (springidx.size() != 0) {
			append_paths(graph, newPaths, newPaths_cost, newPaths_vadose, springidx);
		}
		float result1 = params.nb_cycles != 0 ? mean_cycles_length / params.nb_cycles : 0.;
		float result2 = params.nb_cycles != 0 ? mean_cycles_length_inter_wt / params.nb_cycles : 0.;
		return std::make_pair(result1, result2);
	}

	void KarsticSkeleton::amplify_vadose(GraphOperations* graph, const GeologicalParameters& params)
	{
		int vadose_nodes = count_average_vadose_nodes();
		int max_loops = int(params.loop_density_vadose * vadose_nodes); // number of loops is number of vadose nodes times the density fraction defined by the user

		std::vector<std::vector<int>> newPaths;
		std::vector<std::vector<float>> newPaths_cost;
		std::vector<std::vector<char>> newPaths_vadose;
		std::vector<int> springidx;

		const int idx_specific_wt = select_lowest_spring_cost_channel_1based(params);

		while (max_loops > 0)
		{
			int u = int(generateRandomFloat(0, float(nodes.size()))); // Select two random vertices of the main graph
			int v = int(generateRandomFloat(0, float(nodes.size())));
			int int_u_index = nodes[u].index;
			int int_v_index = nodes[v].index;

			// Check if edge already exists or u and v are the same vertex
			if (u != v && !edge_exists(u, v) && always_in_vadose(u) && always_in_vadose(v)) {
				float distance = magnitude(nodes[u].p - nodes[v].p);

				if (distance <= params.max_dist_loops_vadose) { // also, check if the euclidian distance separating them is less than the threshold value defined by the user

					std::pair< std::vector<int>, std::vector<float>> pair = graph->shortest_path(int_u_index, int_v_index, idx_specific_wt); // if all the previous tests were ok, then we can compute a path between those two nodes
					std::vector<int> bestPath = pair.first;
					std::vector<float> bestPath_cost = pair.second;

					newPaths.push_back(std::vector<int>(bestPath.begin(), bestPath.end()));
					newPaths_cost.push_back(std::vector<float>(bestPath_cost.begin(), bestPath_cost.end()));
					newPaths_vadose.push_back(std::vector<char>(bestPath.size(), true)); // all paths are vadose paths here
					springidx.push_back(0);
					--max_loops; // and then decrement the number of loops still to be generated

				}
			}
		}
		// Update the internal karstic skeleton struture
		append_paths(graph, newPaths, newPaths_cost, newPaths_vadose, springidx);
	}

	void KarsticSkeleton::amplify_phreatic(GraphOperations* graph, const GeologicalParameters& params)
	{

		int vadose_nodes = count_average_vadose_nodes();
		int max_loops = int(params.loop_density_phreatic * (nodes.size() - vadose_nodes)); // number of loops is number of phreatic nodes times the density fraction defined by the user

		std::vector<std::vector<int>> newPaths;
		std::vector<std::vector<float>> newPaths_cost;
		std::vector<std::vector<char>> newPaths_vadose;
		std::vector<int> springidx;

		const int idx_specific_wt = select_lowest_spring_cost_channel_1based(params);

		while (max_loops > 0)
		{
			int u = int(generateRandomFloat(0, float(nodes.size()))); // Select two random vertices of the main graph
			int v = int(generateRandomFloat(0, float(nodes.size())));
			int int_u_index = nodes[u].index; // Store those vertices' indices in the graph
			int int_v_index = nodes[v].index;

			// Check if edge already exists or u and v are the same vertex
			if (u != v && !edge_exists(u, v) && !always_in_vadose(u) && !always_in_vadose(v)) {
				float distance = magnitude(nodes[u].p - nodes[v].p);
				if (distance <= params.max_dist_loops_phreatic) { // also, check if the euclidian distance separating them is less than the threshold value defined by the user

					std::pair< std::vector<int>, std::vector<float>> pair = graph->shortest_path(int_u_index, int_v_index, idx_specific_wt); // if all the previous tests were ok, then we can compute a path between those two nodes
					std::vector<int> bestPath = pair.first;
					std::vector<float> bestPath_cost = pair.second;
					newPaths.push_back(std::vector<int>(bestPath.begin(), bestPath.end()));
					newPaths_cost.push_back(std::vector<float>(bestPath_cost.begin(), bestPath_cost.end()));
					newPaths_vadose.push_back(std::vector<char>(bestPath.size(), false)); // all paths are phreatic paths here
					springidx.push_back(0);
					--max_loops; // and then decrement the number of loops still to be generated
				}
			}
		}

		// Update the internal karstic skeleton struture
		append_paths(graph, newPaths, newPaths_cost, newPaths_vadose, springidx);
	}

	void KarsticSkeleton::amplify_deadend(
		GraphOperations* graph,
		float max_distance_of_deadend_pts,
		int nb_deadend_points,
		const GeologicalParameters& params
	) {

		if (nb_deadend_points == 0) {
			return;
		}

		const auto new_deadend_pts = generate_deadend_points(
			graph,
			max_distance_of_deadend_pts,
			nb_deadend_points,
			params
		);

		const std::vector<int>& deadend_indices =
			new_deadend_pts.second.first;
		const std::vector<int>& origin_indices =
			new_deadend_pts.second.second;

		if (deadend_indices.size() !=
			static_cast<std::size_t>(nb_deadend_points) ||
			origin_indices.size() !=
			static_cast<std::size_t>(nb_deadend_points)) {
			throw std::runtime_error(
				"[deadend amplification] Dead-end generation returned an "
				"inconsistent number of point indices."
			);
		}

		std::vector<std::vector<int>> newPaths;
		std::vector<std::vector<float>> newPaths_cost;
		std::vector<std::vector<char>> newPaths_vadose;
		std::vector<int> springidx;

		newPaths.reserve(nb_deadend_points);
		newPaths_cost.reserve(nb_deadend_points);
		newPaths_vadose.reserve(nb_deadend_points);
		springidx.reserve(nb_deadend_points);

		const int idx_specific_wt =
			select_lowest_spring_cost_channel_1based(params);

		for (int i = 0; i < nb_deadend_points; ++i) {
			const auto pair = graph->shortest_path(
				deadend_indices[i],
				origin_indices[i],
				idx_specific_wt
			);

			const std::vector<int>& bestPath = pair.first;
			const std::vector<float>& bestPath_cost = pair.second;

			if (bestPath.size() < 2u) {
				throw std::runtime_error(
					"[deadend amplification] No graph path exists between a generated "
					"dead-end point and its selected skeleton origin."
				);
			}

			if (bestPath_cost.size() != bestPath.size()) {
				throw std::runtime_error(
					"[deadend amplification] Shortest-path output is inconsistent: "
					"the node and cost vectors do not have the same size."
				);
			}

			newPaths.push_back(bestPath);
			newPaths_cost.push_back(bestPath_cost);
			newPaths_vadose.emplace_back(bestPath.size(), false);
			springidx.push_back(0);
		}

		append_paths(
			graph,
			newPaths,
			newPaths_cost,
			newPaths_vadose,
			springidx
		);
	}

	std::pair<std::vector<Vector3>,std::pair<std::vector<int>, std::vector<int>>>
		KarsticSkeleton::generate_deadend_points(
			GraphOperations* graph,
			float max_distance_of_deadend_pts,
			int nb_deadend_points,
			const GeologicalParameters& params
		) {
		using DeadendResult = std::pair<
			std::vector<Vector3>,
			std::pair<std::vector<int>, std::vector<int>>
		>;

		if (nb_deadend_points <= 0) {
			return DeadendResult{};
		}

		if (graph == nullptr) {
			throw std::runtime_error(
				"[deadend amplification] Cannot generate dead-end points from a null graph."
			);
		}

		if (!std::isfinite(max_distance_of_deadend_pts) ||
			max_distance_of_deadend_pts <= 0.0f) {
			throw std::runtime_error(
				"[deadend amplification] The maximum dead-end distance must be finite and > 0."
			);
		}

		std::vector<Vector3> samples_stretched;
		samples_stretched.reserve(graph->samples.size());

		for (const Vector3& sample : graph->samples) {
			samples_stretched.emplace_back(
				sample.x,
				sample.y,
				sample.z * params.vertical_distance_stretching_factor
			);
		}

		PointCloud pointcloud(samples_stretched);

		struct EligibleOrigin {
			int graph_node_index;
			std::vector<Neighbour> candidates;
		};

		std::vector<EligibleOrigin> eligible_origins;
		eligible_origins.reserve(nodes.size());

		const float squared_max_distance =
			max_distance_of_deadend_pts *
			max_distance_of_deadend_pts;

		// Rejection sampling over skeleton origins is equivalent to sampling
		// uniformly among origins that actually possess at least one valid candidate.
		for (const KarsticNode& node : nodes) {
			const Vector3 stretched_position(
				node.p.x,
				node.p.y,
				node.p.z * params.vertical_distance_stretching_factor
			);

			std::vector<Neighbour> candidates;

			pointcloud.findNearestNeighbors(
				stretched_position,
				-1,
				10000,
				squared_max_distance,
				candidates
			);

			// Remove the origin itself and any numerically coincident sample.
			candidates.erase(
				std::remove_if(
					candidates.begin(),
					candidates.end(),
					[](const Neighbour& candidate) {
				return candidate.d <= 1e-4f;
			}
				),
				candidates.end()
			);

			if (!candidates.empty()) {
				eligible_origins.push_back({
					node.index,
					std::move(candidates)
					});
			}
		}

		if (eligible_origins.empty()) {
			throw std::runtime_error(
				"[deadend amplification] No sampling point distinct from the skeleton "
				"was found within max_distance_of_deadend_pts. Dead-end generation "
				"cannot proceed with the current graph and distance threshold."
			);
		}

		std::vector<Vector3> new_points;
		std::vector<int> new_point_indices;
		std::vector<int> closest_skeleton_indices;

		new_points.reserve(nb_deadend_points);
		new_point_indices.reserve(nb_deadend_points);
		closest_skeleton_indices.reserve(nb_deadend_points);

		for (int i = 0; i < nb_deadend_points; ++i) {
			const std::size_t origin_position =
				generateRandomIndex(eligible_origins.size());

			const EligibleOrigin& origin =
				eligible_origins[origin_position];

			const std::size_t candidate_position =
				generateRandomIndex(origin.candidates.size());

			const Neighbour& selected_candidate =
				origin.candidates[candidate_position];

			new_points.push_back(
				graph->get_sample(selected_candidate.i)
			);
			new_point_indices.push_back(
				selected_candidate.i
			);
			closest_skeleton_indices.push_back(
				origin.graph_node_index
			);
		}

		return {
			std::move(new_points),
			{
				std::move(new_point_indices),
				std::move(closest_skeleton_indices)
			}
		};
	}

	float KarsticSkeleton::compute_wt_ratio(const GraphOperations& graph, const std::vector<Surface>& water_tables) {
		float ratio = 0;

		float dist_btw_wt = magnitude(water_tables.at(0).get_boundbox_max() - water_tables.at(1).get_boundbox_max());


		for (auto& node : nodes) {

			float dist1 = graph.distsurf(node.p, water_tables[0], std::numeric_limits<float>::infinity(), water_tables[0].get_centers_cloud());
			float dist2 = graph.distsurf(node.p, water_tables[1], std::numeric_limits<float>::infinity(), water_tables[1].get_centers_cloud());

			if (dist1*dist1 + dist2 * dist2 < 0.95*dist_btw_wt*dist_btw_wt) {
				ratio += 1;
			}

		}
		return ratio / nodes.size();
	}

	float KarsticSkeleton::compute_mean_deviation() {
		float mean_dev = 0;
		int count = 0;
		for (auto& node : nodes) {
			for (auto& connection : node.connections) {
				for (int i = 0; i < nodes[connection.destindex].connections.size(); i++) {
					Vector3 vec_nc = nodes[connection.destindex].p - node.p;
					Vector3 vec_ccbis = nodes[nodes[connection.destindex].connections[i].destindex].p - nodes[connection.destindex].p;

					float mag_nc = magnitude(vec_nc);
					float mag_ccbis = magnitude(vec_ccbis);

					// Check for zero magnitude vectors to prevent division by zero
					if (mag_nc != 0.0 && mag_ccbis != 0.0) {
						float dot_product = Dot(vec_nc, vec_ccbis);
						float denom = mag_nc * mag_ccbis;
						// Clamping the fraction to the range [-1, 1] to ensure it is a valid input for acos
						float fraq = static_cast<float>(std::max(-1.0, std::min(1.0, static_cast<double>(dot_product / denom))));
						float dev = acos(fraq);
						mean_dev += dev * dev;
					}

					count++;
				}
			}
		}

		return mean_dev / count;
	}

	void KarsticSkeleton::detect_intersection_points(const std::vector<std::vector<int>>& paths) {
		std::vector<int> stack;
		stack.push_back(0);
		intersection_points.clear();
		intersection_points.resize(nodes.size(), 0); //default flag

		// Refine two-way connections while preserving the selected directed edge.
		for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
			for (int j = 0;
				j < static_cast<int>(nodes[i].connections.size());
				++j) {

				// Store the destination before modifying its connection vector.
				const int destination_index =
					nodes[i].connections[j].destindex;

				auto& reverse_connections =
					nodes[destination_index].connections;

				for (auto it = reverse_connections.begin();
					it != reverse_connections.end();) {

					if (it->destindex == i) {
						it = reverse_connections.erase(it);
					}
					else {
						++it;
					}
				}
			}
		}

		while (!stack.empty()) {
			int current_index = stack.back();
			stack.pop_back();

			KarsticNode& current_node = nodes[current_index];
			if (intersection_points[current_index] == 0) {	// is usfull (rarely)
				if (current_node.connections.size() != 1) {
					intersection_points[current_index] = -1;   //flag intersection
				}
				else {
					intersection_points[current_index] = 1;  //flag visited
				}
			}

			// Traverse neighbors
			for (const KarsticConnection& connection : current_node.connections) {
				int neighbor_index = connection.destindex;
				//KarsticNode& neighbor_node = nodes[neighbor_index];

				// If the neighbor hasn't been visited yet, add it to the stack
				if (intersection_points[neighbor_index] == 0 && std::find(stack.begin(), stack.end(), neighbor_index) == stack.end()) {
					stack.push_back(neighbor_index);
				}
				// If the neighbor has been visited, it's an new intersection
				else {
					intersection_points[neighbor_index] = -1;
				}
			}

		}

		// set dolines and source as intersections
		for (auto& path : paths) {
			intersection_points[get_internal_index(path[0])] = -1;
		}
		intersection_points[0] = -1;
	}

	float KarsticSkeleton::compute_mean_branch_length() {

		std::vector<float> vec_mean_length;
		std::vector<std::vector<int>> vec_seen_branch_id;

		int i = 0;
		for (auto& node : nodes) {
			for (auto& connection : node.connections) {
				if (intersection_points[connection.destindex] == -1) { // if destination is an intersection
					if (intersection_points[i] == -1) {	 // && current is an intersection
						vec_seen_branch_id.push_back({ -99999 });
						vec_mean_length.push_back(magnitude(node.p - nodes[connection.destindex].p));
					}
					else {
						auto it = std::find(vec_seen_branch_id.begin(), vec_seen_branch_id.end(), node.branch_id_ascend);
						if (it == vec_seen_branch_id.end()) {
							vec_seen_branch_id.push_back(node.branch_id_ascend);
							vec_mean_length.push_back(magnitude(node.p - nodes[connection.destindex].p));
						}
						else {
							vec_mean_length[it - vec_seen_branch_id.begin()] += magnitude(node.p - nodes[connection.destindex].p);
						}
					}
				}
				else {	//if destination is on branch

					auto it2 = std::find(vec_seen_branch_id.begin(), vec_seen_branch_id.end(), nodes[connection.destindex].branch_id_ascend);
					if (it2 == vec_seen_branch_id.end()) {
						vec_seen_branch_id.push_back(nodes[connection.destindex].branch_id_ascend);
						vec_mean_length.push_back(magnitude(node.p - nodes[connection.destindex].p));
					}
					else {
						vec_mean_length[it2 - vec_seen_branch_id.begin()] += magnitude(node.p - nodes[connection.destindex].p);
					}
				}
			}
			i++;
		}

		//calul mean lenth
		float mean_length = 0.;
		for (auto& length : vec_mean_length) {
			mean_length += length;
		}

		return mean_length / vec_mean_length.size();
	}

	int KarsticSkeleton::compute_nb_cycles() {

		// Refine two-way connections while preserving the selected directed edge.
		for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
			for (int j = 0;
				j < static_cast<int>(nodes[i].connections.size());
				++j) {

				// Store the destination before modifying its connection vector.
				const int destination_index =
					nodes[i].connections[j].destindex;

				auto& reverse_connections =
					nodes[destination_index].connections;

				for (auto it = reverse_connections.begin();
					it != reverse_connections.end();) {

					if (it->destindex == i) {
						it = reverse_connections.erase(it);
					}
					else {
						++it;
					}
				}
			}
		}

		int nb_segments = 0;
		int i = 0;
		for (auto& node : nodes) {

			nb_segments += int(node.connections.size());
			i++;
		}
		int nb_cycle = nb_segments - int(nodes.size()) + 1;
		return nb_cycle;
	}

	int KarsticSkeleton::get_internal_index(int nodeIndex) const
	{
		for (int i = 0; i < nodes.size(); i++)
		{
			if (nodes[i].index == nodeIndex)
				return i;
		}
		return -1;
	}

	int KarsticSkeleton::get_internal_index(Vector3 node) const
	{
		for (int i = 0; i < nodes.size(); i++)
		{
			if (nodes[i].p == node)
				return i;
		}
		return -1;
	}

	void KarsticSkeleton::append_paths(const GraphOperations* graph, const std::vector<std::vector<int>>& paths, const std::vector<std::vector<float>> costs, const  std::vector<std::vector<char>> vadose, const std::vector<int>& springidx)
	{

		int nb_springs = *std::max_element(springidx.begin(), springidx.end()) + 1;

		for (int i = 0; i < paths.size(); i++)	// Nodes
		{
			int path_size = int(paths[i].size()); // path size
			for (int j = path_size - 1; j >= 0; j--)
			{
				int node_idx = paths[i][j];
				std::vector<float> nodeCost(nb_springs, 0.);
				nodeCost[springidx[i]] = costs[i][j];
				std::vector<int> vadoseTest(nb_springs, -1);
				vadoseTest[springidx[i]] = vadose[i][j];
				auto it = std::find(nodes.begin(), nodes.end(), node_idx);
				if (it == nodes.end()) // if the node hasn't been added yet
				{
					Vector3 p = graph->get_sample(node_idx);
					nodes.push_back({ node_idx,p ,nodeCost, vadoseTest });
				}
				else { // else: don't add it again, but overwrite the vadose flag and the cost
					int index = std::distance(nodes.begin(), it);
					nodes[index].vadose[springidx[i]] = vadose[i][j];
					nodes[index].cost[springidx[i]] = costs[i][j];
				}
			}

		}
		for (int i = 0; i < nodes.size(); i++)	// Edges
		{
			int nodeIndex = nodes[i].index;
			std::vector<int> edges;
			for (int j = 0; j < paths.size(); j++)
			{
				for (int k = int(paths[j].size()) - 1; k >= 0; k--)
				{
					int n = paths[j][k];
					if (nodeIndex != n)
						continue;
					if (k > 0)
					{
						int nAfter = paths[j][k - 1];
						int neiIndex = get_internal_index(nAfter);
						if (std::find(edges.begin(), edges.end(), neiIndex) == edges.end())
							edges.push_back(neiIndex);
					}
				}
			}
			for (int j = 0; j < edges.size(); j++)
			{
				auto new_connection = KarsticConnection({ edges[j] });
				if (std::find(nodes[i].connections.begin(), nodes[i].connections.end(), edges[j]) == nodes[i].connections.end())
					nodes[i].connections.push_back(new_connection);
			}
		}
	}

	void KarsticSkeleton::update_branch_ID() {
		// Reset branch IDs
		for (KarsticNode& node : nodes) {
			node.branch_id = -10; // Resetting to default value
		}

		int branch_id_counter = 0; // Counter for assigning branch IDs

		// Perform DFS traversal to identify branches
		for (int i = 0; i < nodes.size(); ++i) {
			if (nodes[i].branch_id == -10) {

				// Start a new branch exploration

				//if (nodes[i].connections.size() <= 2) {
				//	nodes[i].branch_id = branch_id_counter;
				//	
				//	// iterate on neighbors

				//	for (int j = 0; j < nodes[i].connections.size(); j++) {
				//		int neigh_idx = nodes[i].connections[j].destindex;
				//		
				//	}


				int number_of_branch_pts = dfs(i, branch_id_counter);
				// Increment branch ID counter after finishing exploration of a branch
				if (number_of_branch_pts > 0) {
					branch_id_counter++;
				}
			}
		}
	}

	int KarsticSkeleton::dfs(int node_index, int branch_id_counter) {

		std::stack<int> stack;
		stack.push(node_index);
		int number_of_branch_pts = 0;

		while (!stack.empty()) {
			int current_index = stack.top();
			stack.pop();

			KarsticNode& current_node = nodes[current_index];

			// Skip nodes already visited
			if (current_node.branch_id != -10) {
				continue;
			}

			// If the current node is an intersection, set its branch ID to -1
			if (current_node.connections.size() > 2) {
				current_node.branch_id = -1;
				continue; // Skip exploration for intersection nodes
			}

			// Assign branch ID to current node
			current_node.branch_id = branch_id_counter;
			number_of_branch_pts++;

			// Traverse neighbors
			for (const KarsticConnection& connection : current_node.connections) {
				int neighbor_index = connection.destindex;
				KarsticNode& neighbor_node = nodes[neighbor_index];

				// If the neighbor hasn't been visited yet, add it to the stack
				if (neighbor_node.branch_id == -10) {
					stack.push(neighbor_index);
				}
			}
		}
		return number_of_branch_pts;
	}

	void KarsticSkeleton::update_branch_ID(const  std::vector<std::vector<int>>& paths) {

		//Reset labelling
		for (int i = 0; i < nodes.size(); i++) {
			nodes[i].branch_id_ascend.clear();
		}

		// Refine two-way connections while preserving the selected directed edge.
		for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
			for (int j = 0;
				j < static_cast<int>(nodes[i].connections.size());
				++j) {

				// Store the destination before modifying its connection vector.
				const int destination_index =
					nodes[i].connections[j].destindex;

				auto& reverse_connections =
					nodes[destination_index].connections;

				for (auto it = reverse_connections.begin();
					it != reverse_connections.end();) {

					if (it->destindex == i) {
						it = reverse_connections.erase(it);
					}
					else {
						++it;
					}
				}
			}
		}

		//create end point vector to avoid branch id to continue labeling the same branche in doline intersection
		std::vector<int> end_points;
		for (int i = 0; i < paths.size(); i++) {
			end_points.push_back(paths[i][0]);
		}

		if (nodes[0].branch_id_ascend.empty()) {
			nodes[0].branch_id_ascend.push_back(1);
		}

		std::vector<int> nodes_index = { 0 };
		int acc = 0;
		std::vector<int> visited_intersection_points;
		while (!nodes_index.empty() && acc < 10000) {
			acc++;
			KarsticNode current_node;
			current_node = nodes[nodes_index.back()];
			nodes_index.pop_back();

			if (current_node.connections.size() == 1 && find(end_points.begin(), end_points.end(), current_node.index) == end_points.end() && intersection_points[get_internal_index(current_node.index)] == 1) {  //current node is on a branch
				int next_index = current_node.connections[0].destindex;
				nodes[next_index].add_branch_id(current_node.branch_id_ascend);
				nodes_index.push_back(next_index);
			}

			else {		//current node is on an intersection
				if (find(visited_intersection_points.begin(), visited_intersection_points.end(), current_node.index) == visited_intersection_points.end()) {
					visited_intersection_points.push_back(current_node.index);

					for (int j = 0; j < current_node.connections.size(); j++) {

						int next_index = current_node.connections[j].destindex;
						nodes_index.push_back(next_index);
						std::vector<int> new_branch_id;
						for (int k = 0; k < current_node.branch_id_ascend.size(); k++) {

							new_branch_id.push_back(current_node.branch_id_ascend[k] * 10 + j);
						}
						nodes[next_index].add_branch_id(new_branch_id);
					}
				}
				else
				{

				}
			}

			if (nodes_index.size() > nodes.size() || acc == 10000) {  //in case infinit cycle (not supposed to happen anymore)
					break;
			}

		}

	}

	void KarsticSkeleton::save(const std::string& file, const std::string& save_directory) const
	{
		std::ofstream out;
		out.open(save_directory + "/outputs/" + file + "_nodes.dat");
		if (out.is_open() == false)
		{
			std::cout << "Cannot save skeleton to file (nodes): " << file << std::endl;
			return;
		}
		std::map<int, int> nodeIndexMapping;
		for (int i = 0; i < nodes.size(); i++)
		{
			Vector3 nodePos = nodes[i].p;
			nodeIndexMapping.insert(std::make_pair(nodes[i].index, i + 1));
			out << std::fixed << std::setprecision(3)
				<< nodePos[0] << " " << nodePos[1] << " " << nodePos[2] << "\n";
		}
		out.flush();
		out.close();

		out.open(save_directory + "/outputs/" + file + "_links.dat");
		if (out.is_open() == false)
		{
			std::cout << "Cannot save skeleton to file (links): " << file << std::endl;
			return;
		}
		for (int i = 0; i < nodes.size(); i++)
		{
			int indexA = i + 1;
			for (int j = 0; j < nodes[i].connections.size(); j++)
			{
				int k = nodes[i].connections[j].destindex;
				out << indexA << " " << (k + 1) << "\n";
			}
		}
		out.flush();
		out.close();
	}

	void KarsticSkeleton::prepare_graph() {
		std::vector<std::pair<int, int>> unique_segments;
		std::vector<bool> removed(nodes.size(), false);

		// Iterate through nodes to find unique segments
		for (size_t i = 0; i < nodes.size(); ++i) {
			const KarsticNode& node = nodes[i];
			for (const KarsticConnection& connection : node.connections) {
				int destindex = connection.destindex;
				if (i < destindex) { // Avoid duplicate segments and self-loops
					const KarsticNode& destNode = nodes[destindex];
					// Check if the segment is unique
					bool found = false;
					for (const auto& seg : unique_segments) {
						if ((seg.first == i && seg.second == destindex) || (seg.first == destindex && seg.second == i)) {
							found = true;
							break;
						}
					}
					if (!found) {
						if (node.p == destNode.p) { // Check if the extremities are at the same location
							removed[i] = true; // Mark one of the nodes to be removed
							removed[destindex] = true;
						}
						else {
							unique_segments.emplace_back(static_cast<int>(i), destindex); // Mark the segment as unique
						}
					}
				}
			}
		}

		// Remove duplicate segments by marking nodes to be removed
		for (size_t i = 0; i < nodes.size(); ++i) {
			if (removed[i]) {
				nodes[i].connections.clear(); // Remove connections from the node
			}
		}

		// Reconnect nodes based on unique segments
		for (const auto& seg : unique_segments) {
			int node1_index = seg.first;
			int node2_index = seg.second;
			nodes[node1_index].add_connection(node2_index);
			nodes[node2_index].add_connection(node1_index);
		}
	}

	KarstNetworkResult KarsticSkeleton::create_line(
		const GeologicalParameters& geologicalparams,
		const std::string& network_name) const
	{
		(void)network_name;

		constexpr float NO_DATA_VALUE = -99999.0f;
		constexpr float EPSILON_COST = 1e-10f;

		KarstNetworkResult result;

		const bool add_drift_props =
			geologicalparams.use_drift_zwt || geologicalparams.use_drift_curv;

		// Total number of cost channels, including the optional internal
		// vadose-only no-water-table channel.
		const int nb_cost_channels = std::max(0, geologicalparams.nb_wt);

		// The no-water-table channel is an internal vadose-only cost channel.
		// It does not carry any discriminating physical water-table flag and
		// must not be exported as vadose_flag_wt_*.
		const int nb_physical_wt =
			geologicalparams.no_wt_cost_index >= 0
			? std::min(geologicalparams.no_wt_cost_index, nb_cost_channels)
			: nb_cost_channels;

		std::vector<int> exported_vadose_channels;

		for (int wt_channel = 0; wt_channel < nb_physical_wt; ++wt_channel) {
			const int wt_index_1based = wt_channel + 1;

			const bool is_used_by_at_least_one_spring =
				std::find_if(
					geologicalparams.propspringswtindex.begin(),
					geologicalparams.propspringswtindex.end(),
					[wt_index_1based](const GeologicalParameters::Propidx& spring_wt_index)
			{
				return static_cast<int>(spring_wt_index.prop) == wt_index_1based;
			}
				) != geologicalparams.propspringswtindex.end();

			if (is_used_by_at_least_one_spring) {
				exported_vadose_channels.push_back(wt_channel);
			}
		}

		result.vadose_property_names.clear();
		for (const int wt_channel : exported_vadose_channels) {
			result.vadose_property_names.push_back(
				"vadose_flag_wt_" + std::to_string(wt_channel + 1)
			);
		}
		result.has_drift_properties = add_drift_props;

		for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
			for (int j = 0; j < static_cast<int>(nodes[i].connections.size()); ++j) {

				const int dest = nodes[i].connections[j].destindex;

				if (dest < 0 || dest >= static_cast<int>(nodes.size())) {
					throw std::runtime_error(
						"[skeleton] Cannot export line: invalid connection index."
					);
				}

				auto it_spring = std::find_if(
					nodes[i].cost.begin(),
					nodes[i].cost.end(),
					[EPSILON_COST](float c)
				{
					return std::isfinite(c) && std::fabs(c) > EPSILON_COST;
				}
				);

				int index_spring = -1;
				if (it_spring != nodes[i].cost.end()) {
					index_spring = static_cast<int>(
						std::distance(nodes[i].cost.begin(), it_spring)
						);
				}

				float start_cost = 0.0f;
				float end_cost = 0.0f;

				if (index_spring >= 0) {
					if (index_spring < static_cast<int>(nodes[i].cost.size())) {
						start_cost = nodes[i].cost[index_spring];
					}
					if (index_spring < static_cast<int>(nodes[dest].cost.size())) {
						end_cost = nodes[dest].cost[index_spring];
					}
				}

				if (!std::isfinite(start_cost) || start_cost <= 0.0f) {
					start_cost = 0.0f;
				}
				if (!std::isfinite(end_cost) || end_cost <= 0.0f) {
					end_cost = 0.0f;
				}

				ResultPoint start;
				ResultPoint end;

				start.p = nodes[i].p;
				end.p = nodes[dest].p;

				start.cost = start_cost;
				end.cost = end_cost;

				start.equivalent_radius = static_cast<float>(nodes[i].eq_radius);
				end.equivalent_radius = static_cast<float>(nodes[dest].eq_radius);

				// Export the branch IDs computed by update_branch_ID().
				// Intersections keep branch_id == -1; branch numbering starts at 0.
				start.branch_id = nodes[i].branch_id;
				end.branch_id = nodes[dest].branch_id;

				start.vadose_flags.clear();
				end.vadose_flags.clear();

				for (const int wt_channel : exported_vadose_channels) {
					float start_vadose = NO_DATA_VALUE;
					float end_vadose = NO_DATA_VALUE;

					if (wt_channel < static_cast<int>(nodes[i].vadose.size()) &&
						nodes[i].vadose[wt_channel] != -1) {
						start_vadose = static_cast<float>(nodes[i].vadose[wt_channel]);
					}

					if (wt_channel < static_cast<int>(nodes[dest].vadose.size()) &&
						nodes[dest].vadose[wt_channel] != -1) {
						end_vadose = static_cast<float>(nodes[dest].vadose[wt_channel]);
					}

					start.vadose_flags.push_back(start_vadose);
					end.vadose_flags.push_back(end_vadose);
				}

				// Optional legacy scalar fallback, if ResultPoint::vadose_flag is still
				// used somewhere else in the public code or bindings.
				start.vadose_flag =
					start.vadose_flags.empty() ? NO_DATA_VALUE : start.vadose_flags.front();
				end.vadose_flag =
					end.vadose_flags.empty() ? NO_DATA_VALUE : end.vadose_flags.front();

				if (add_drift_props) {
					start.external_drift = nodes[i].drift_value;
					start.kriging_weight = nodes[i].drift_weight;

					end.external_drift = nodes[dest].drift_value;
					end.kriging_weight = nodes[dest].drift_weight;
				}

				result.add_segment(start, end);
			}
		}

		return result;
	}

	// Method to populate distance matrix using Johnson's algorithm (which finds shortest distance between a pair of nodes in a graph)
	void KarsticSkeleton::compute_distance_matrix() {
		// Initialize distance matrix with infinity
		distance_mat.resize(nodes.size(), nodes.size(), std::numeric_limits<float>::infinity());

		// Populate distance matrix with direct distances where there is an edge
		for (size_t i = 0; i < nodes.size(); ++i) {
			for (size_t j = 0; j < nodes.size(); ++j) {
				if (edge_exists(int(i), int(j))) {
					distance_mat(i, j) = magnitude(nodes[i].p - nodes[j].p);
				}
			}
		}

		int n = int(nodes.size());

		// Step 1: Add a new node to the graph and connect it to all other nodes with zero-weight edges
		std::vector<std::vector<float>> new_graph(n + 1, std::vector<float>(n + 1, std::numeric_limits<float>::infinity()));
		for (int i = 0; i < n; ++i) {
			new_graph[n][i] = 0; // Connect the new node to all other nodes with zero-weight edges
		}

		// Step 2: Run Bellman-Ford algorithm from the new node to get minimum distances from the new node to all other nodes
		std::vector<float> h(n + 1, std::numeric_limits<float>::infinity());
		h[n] = 0;

		for (int k = 0; k < n; ++k) {
			for (int i = 0; i <= n; ++i) {
				for (int j = 0; j <= n; ++j) {
					if (new_graph[i][j] != std::numeric_limits<float>::infinity()) {
						if (h[i] + new_graph[i][j] < h[j]) {
							h[j] = h[i] + new_graph[i][j];
						}
					}
				}
			}
		}

		// Step 3: Reweight the edges only if there's no loop
		bool has_loop = false;
		for (size_t i = 0; i < nodes.size(); ++i) {
			if (edge_exists(int(i), int(i))) {
				has_loop = true;
				break;
			}
		}

		if (!has_loop) {
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					if (distance_mat(i, j) != std::numeric_limits<float>::infinity()) {
						distance_mat(i, j) += h[i] - h[j];
					}
				}
			}
		}

		// Step 4: Run Dijkstra's algorithm for each node to compute shortest paths
		for (int src = 0; src < n; ++src) {
			std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>> pq;
			pq.push({ 0, src });
			std::vector<float> dist(n, std::numeric_limits<float>::infinity());
			dist[src] = 0;

			while (!pq.empty()) {
				int u = pq.top().second;
				pq.pop();

				for (const KarsticConnection& conn : nodes[u].connections) {
					int v = conn.destindex;
					float weight = distance_mat(u, v);
					if (dist[v] > dist[u] + weight) {
						dist[v] = dist[u] + weight;
						pq.push({ dist[v], v });
					}
				}
			}

			// Update the distance matrix with the shortest distances from src to all other nodes
			for (int i = 0; i < n; ++i) {
				distance_mat(src, i) = dist[i];
			}
		}

		// Step 5: Reverse the reweighting only if there's no loop
		if (!has_loop) {
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					if (distance_mat(i, j) != std::numeric_limits<float>::infinity()) {
						distance_mat(i, j) -= h[i] - h[j];
					}
				}
			}
		}
	}

	// Method to compute the number of neighbors for each node (valence)
	void KarsticSkeleton::compute_valence() {
		valence.resize(nodes.size(), 0);

		for (size_t node = 0; node < nodes.size(); node++) {
			for (size_t i = 0; i < nodes.at(node).connections.size(); ++i) {
				++valence.at(node); // Increment valence for the current node
			}
		}
	}

	void KarsticSkeleton::refresh_vadose_flags_from_graph(
		GraphOperations* graph,
		const GeologicalParameters& params)
	{
		if (graph == nullptr) {
			throw std::runtime_error(
				"[skeleton] Cannot refresh vadose flags: graph pointer is null."
			);
		}

		Array2D<char>& graph_wt_flags = graph->get_samples_surf_flags();

		for (KarsticNode& node : nodes) {
			if (node.index < 0 || node.index >= graph->get_nb_sampling_pts()) {
				throw std::runtime_error(
					"[skeleton] Cannot refresh vadose flags: a skeleton node has an invalid graph index."
				);
			}

			node.vadose.assign(params.nb_wt, -1);

			for (int channel = 0; channel < params.nb_wt; ++channel) {
				// In samples_surf_flags, true means below the water table, i.e. phreatic.
				node.vadose[channel] = graph_wt_flags(node.index, channel) ? 0 : 1;
			}
		}
	}

	// Method to compute and populate branch_sizes attribute
	void KarsticSkeleton::compute_branch_sizes() {
		std::unordered_map<int, int> branchSizeMap;

		branch_sizes.clear();
		nb_of_intersections = 0;

		for (const KarsticNode& node : nodes) {
			if (node.branch_id == -1) {
				nb_of_intersections++;
			}
			else {
				++branchSizeMap[node.branch_id];
			}
		}

		std::vector<std::pair<int, int>> sortedBranchSizes(branchSizeMap.begin(), branchSizeMap.end());
		std::sort(sortedBranchSizes.begin(), sortedBranchSizes.end());

		for (const auto& pair : sortedBranchSizes) {
			branch_sizes.push_back(pair.second);
		}
	}
}