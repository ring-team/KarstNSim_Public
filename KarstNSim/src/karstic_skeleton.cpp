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

		// This scans z_list, finds the index in KeyPts associated with the lowest altitude spring, and finds its associated water table index
		auto min_it = std::min_element(params.z_list.begin(), params.z_list.end(),
			[](const GeologicalParameters::Propidx& a, const GeologicalParameters::Propidx& b) {
			return a.prop < b.prop; // Compare based on prop
		});
		int order = min_it->index; // This is the index of the smallest prop
		auto it = std::find_if(params.propspringswtindex.begin(), params.propspringswtindex.end(),
			[order](const GeologicalParameters::Propidx& item) {
			return item.index == order; // Condition to match the index
		});
		int idx_specific_wt = it->prop;

		while (cycles_created < params.nb_cycles)
		{
			int v = int(generateRandomFloat(0, float(nodes.size()))); // Select two random vertices of the main graph
			int u = int(generateRandomFloat(0, float(nodes.size()))); //Be sure to avoid previous points in the tree order

			int int_u_index = nodes[v].index;
			int int_v_index = nodes[u].index;

			// Check if edge already exists or u and v are the same vertex
			if (u != v && !edge_exists(u, v)) {
				float distance = magnitude(nodes[u].p - nodes[v].p);

				if (distance >= params.min_distance_amplification && distance <= params.max_distance_amplification) { // also, check if the euclidian distance separating them is less than the threshold value defined by the user

					std::pair< std::vector<int>, std::vector<float>> pair;

					if (is_parent(u, v)) {
						float radius = magnitude(nodes[u].p - nodes[v].p) / 4;
						Vector3 center = nodes[u].p + (nodes[v].p - nodes[u].p) / 2;
						std::vector<std::pair<int, float>> reset_values = graph->add_ball_cost(center, radius);
						pair = graph->shortest_path(int_u_index, int_v_index, idx_specific_wt); // if all the previous tests were ok, then we can compute a path between those two nodes
						for (auto& sample_pair : reset_values) {
							for (int i = 0; i < graph->adj[sample_pair.first].size(); i++) {
								if (graph->adj[sample_pair.first][i].target >= 0) { // avoid empty neighbors (when n<N) 
									for (int cost_i = 0; cost_i < params.nb_springs; cost_i++) {
										graph->adj[sample_pair.first][i].weight[cost_i] = sample_pair.second;
									}
								}
							}
						}
					}
					else {
						pair = graph->shortest_path(int_u_index, int_v_index, idx_specific_wt); // if all the previous tests were ok, then we can compute a path between those two nodes
					}
					std::vector<int> bestPath = pair.first;
					std::vector<float> bestPath_cost = pair.second;

					bool inter_wt = false;
					if (nodes[0].vadose.size() == 2) {
						if (((nodes[u].vadose[1] == 0 && nodes[v].vadose[1] == 0) && (nodes[u].vadose[0] == 1 && nodes[v].vadose[0] == 1)) ||
							((nodes[v].vadose[1] == 0 && nodes[u].vadose[1] == 0) && (nodes[v].vadose[0] == 1 && nodes[u].vadose[0] == 1))) {
							inter_wt = true;
						}
					}

					for (int i = 0; i < bestPath.size() - 1; i++) {
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

		// This scans z_list, finds the index in KeyPts associated with the lowest altitude spring, and finds its associated water table index
		auto min_it = std::min_element(params.z_list.begin(), params.z_list.end(),
			[](const GeologicalParameters::Propidx& a, const GeologicalParameters::Propidx& b) {
			return a.prop < b.prop; // Compare based on prop
		});
		int order = min_it->index; // This is the index of the smallest prop
		auto it = std::find_if(params.propspringswtindex.begin(), params.propspringswtindex.end(),
			[order](const GeologicalParameters::Propidx& item) {
			return item.index == order; // Condition to match the index
		});
		int idx_specific_wt = it->prop;

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

		// This scans z_list, finds the index in KeyPts associated with the lowest altitude spring, and finds its associated water table index
		auto min_it = std::min_element(params.z_list.begin(), params.z_list.end(),
			[](const GeologicalParameters::Propidx& a, const GeologicalParameters::Propidx& b) {
			return a.prop < b.prop; // Compare based on prop
		});
		int order = min_it->index; // This is the index of the smallest prop
		auto it = std::find_if(params.propspringswtindex.begin(), params.propspringswtindex.end(),
			[order](const GeologicalParameters::Propidx& item) {
			return item.index == order; // Condition to match the index
		});
		int idx_specific_wt = it->prop;

		while (max_loops > 0)
		{
			int u = int(generateRandomFloat(0, float(nodes.size()))); // Select two random vertices of the main graph
			int v = int(generateRandomFloat(0, float(nodes.size())));
			int int_u_index = nodes[u].index; // Store those vertices' indices in the volumetric graph
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

	void KarsticSkeleton::amplify_deadend(GraphOperations* graph, float max_distance_of_deadend_pts, int nb_deadend_points, const GeologicalParameters& params) {

		std::pair<std::vector<Vector3>, std::pair<std::vector<int>, std::vector<int>>> new_deadend_pts = generate_deadend_points(graph, max_distance_of_deadend_pts, nb_deadend_points);
		std::vector<std::vector<int>> newPaths;
		std::vector<std::vector<float>> newPaths_cost;
		std::vector<std::vector<char>> newPaths_vadose;
		std::vector<int> springidx;

		// This scans z_list, finds the index in KeyPts associated with the lowest altitude spring, and finds its associated water table index
		auto min_it = std::min_element(params.z_list.begin(), params.z_list.end(),
			[](const GeologicalParameters::Propidx& a, const GeologicalParameters::Propidx& b) {
			return a.prop < b.prop; // Compare based on prop
		});
		int order = min_it->index; // This is the index of the smallest prop
		auto it = std::find_if(params.propspringswtindex.begin(), params.propspringswtindex.end(),
			[order](const GeologicalParameters::Propidx& item) {
			return item.index == order; // Condition to match the index
		});
		int idx_specific_wt = it->prop;

		for (int i = 0; i < nb_deadend_points; i++)
		{
			// Check if edge already exists or u and v are the same vertex
			std::pair<std::vector<int>, std::vector<float>> pair = graph->shortest_path(new_deadend_pts.second.first[i], new_deadend_pts.second.second[i], idx_specific_wt);
			std::vector<int> bestPath = pair.first;
			std::vector<float> bestPath_cost = pair.second;
			newPaths.push_back(std::vector<int>(bestPath.begin(), bestPath.end()));
			newPaths_cost.push_back(std::vector<float>(bestPath_cost.begin(), bestPath_cost.end()));
			newPaths_vadose.push_back(std::vector<char>(bestPath.size(), false)); // all paths are phreatic paths here
			springidx.push_back(0);
		}

		// Update the internal karstic skeleton struture
		append_paths(graph, newPaths, newPaths_cost, newPaths_vadose, springidx);
	}

	std::pair<std::vector<Vector3>, std::pair<std::vector<int>, std::vector<int>>> KarsticSkeleton::generate_deadend_points(GraphOperations* graph, float max_distance_of_deadend_pts, int nb_deadend_points)
	{
		int z_stretch_factor = 50; // useful to not create too many deadends far in z from the main branches (can be unrealistic)
		std::vector<Vector3> samples_stretched;
		for (int i = 0; i < graph->samples.size(); i++) {
			Vector3 p_stretched{
				graph->samples[i].x,
				graph->samples[i].y,
				graph->samples[i].z*z_stretch_factor
			};
			samples_stretched.push_back(p_stretched);
		}

		PointCloud pointcloud(samples_stretched);

		// Compute random points inside the box
		std::pair<std::vector<Vector3>, std::pair<std::vector<int>, std::vector<int>>> new_pts_and_idx;
		std::vector<Vector3> new_pts;
		std::vector<int> new_pts_idx;
		std::vector<int> closest_pts_idx;
		float distmax = max_distance_of_deadend_pts * max_distance_of_deadend_pts;
		int nb_added_pts = 0;
		while (nb_added_pts < nb_deadend_points) {

			// select a random node of the graph and get its nearest neighbors in the graph to choose a random point among them
			int random_idx = int(generateRandomInt(0, int(nodes.size()) - 1));
			Vector3 vec(nodes[random_idx].p.x, nodes[random_idx].p.y, nodes[random_idx].p.z*z_stretch_factor);
			std::vector<Neighbour> candidates;
			pointcloud.findNearestNeighbors(vec, -1, 10000, distmax, candidates); // find closest (at most) 10000 neighbors within a distance threshold from the random node
			int random_idx_candidates = int(generateRandomInt(0, int(candidates.size()) - 1));
			Vector3 p = graph->get_sample(candidates[random_idx_candidates].i); // get random neighbor
			if (candidates[random_idx_candidates].d > 1e-4) {
				new_pts.push_back(p);
				new_pts_idx.push_back(candidates[random_idx_candidates].i);
				vec.z /= z_stretch_factor; // undo the stretching in z
				closest_pts_idx.push_back(graph->NodeIndex(vec));
				nb_added_pts++;
			}
		}

		new_pts_and_idx.first = new_pts;
		new_pts_and_idx.second.first = new_pts_idx;
		new_pts_and_idx.second.second = closest_pts_idx;
		return new_pts_and_idx;
	}

	float KarsticSkeleton::compute_wt_ratio(GraphOperations* graph, std::vector<Surface>* water_tables) {
		float ratio = 0;

		float dist_btw_wt = magnitude(water_tables->at(0).get_boundbox_max() - water_tables->at(1).get_boundbox_max());


		for (auto& node : nodes) {

			float dist1 = graph->distsurf(node.p, &water_tables->at(0), std::numeric_limits<float>::infinity(), water_tables->at(0).get_centers_cloud());
			float dist2 = graph->distsurf(node.p, &water_tables->at(1), std::numeric_limits<float>::infinity(), water_tables->at(1).get_centers_cloud());

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

		//Refine two way connection
		for (int i = 0; i < nodes.size(); i++) {
			for (int j = 0; j < nodes[i].connections.size(); j++) {

				int k = 0;
				for (auto& el : nodes[nodes[i].connections[j].destindex].connections) {
					if (el.destindex == i) {
						nodes[nodes[i].connections[j].destindex].connections.erase(nodes[nodes[i].connections[j].destindex].connections.begin() + k);
					}
					k++;
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

		//Refine two way connection
		for (int i = 0; i < nodes.size(); i++) {
			for (int j = 0; j < nodes[i].connections.size(); j++) {

				int k = 0;
				for (auto& el : nodes[nodes[i].connections[j].destindex].connections) {
					if (el.destindex == i) {
						nodes[nodes[i].connections[j].destindex].connections.erase(nodes[nodes[i].connections[j].destindex].connections.begin() + k);
					}
					k++;
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

		//Refine two way connection
		for (int i = 0; i < nodes.size(); i++) {
			for (int j = 0; j < nodes[i].connections.size(); j++) {

				int k = 0;
				for (auto& el : nodes[nodes[i].connections[j].destindex].connections) {
					if (el.destindex == i) {
						nodes[nodes[i].connections[j].destindex].connections.erase(nodes[nodes[i].connections[j].destindex].connections.begin() + k);
					}
					k++;
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
			out << nodePos[0] << " " << nodePos[1] << " " << nodePos[2] << "\n";
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

	void KarsticSkeleton::create_line(const GeologicalParameters& geologicalparams, std::string network_name) const
	{

		std::vector<Segment> nghb_graph_line;
		std::vector<float> list_cost;
		std::vector<float> list_vadose;
		std::vector<float> list_eq_radius;
		std::vector<int> list_branch_id;
		//std::vector<float> list_branch_dist; //debugging

		for (int i = 0; i < nodes.size(); i++) {
			for (int j = 0; j < nodes[i].connections.size(); j++) {

				float epsilon = 1e-10;
				// Find the index of the first element in cost that is not effectively zero (it corresponds to the spring index to which the path leads to)
				int index_spring = std::distance(nodes[i].cost.begin(), std::find_if(nodes[i].cost.begin(), nodes[i].cost.end(),
					[epsilon](float c) { return std::fabs(c) > epsilon; }));

				list_cost.push_back(float(nodes[i].cost[index_spring]));
				list_cost.push_back(float(nodes[nodes[i].connections[j].destindex].cost[index_spring]));

				list_eq_radius.push_back(float(nodes[i].eq_radius));
				list_eq_radius.push_back(float(nodes[nodes[i].connections[j].destindex].eq_radius));

				list_branch_id.push_back(nodes[i].branch_id);
				list_branch_id.push_back(nodes[nodes[i].connections[j].destindex].branch_id);

				list_vadose.push_back(float(nodes[i].vadose[index_spring]));
				list_vadose.push_back(float(nodes[nodes[i].connections[j].destindex].vadose[index_spring]));

				//list_branch_dist.push_back(float(nodes[i].distance)); // debugging
				//list_branch_dist.push_back(float(nodes[nodes[i].connections[j].destindex].distance));

				Vector3 a0 = nodes[i].p;
				Vector3 a1 = nodes[nodes[i].connections[j].destindex].p;

				// Use those if you want the points to be defined in the Lambert 93 coordinates
				//Vector3 a0 = Vector3(nodes[i].p.x+870000.,nodes[i].p.y+6820000.,nodes[i].p.z);
				//Vector3 a1 = Vector3(nodes[nodes[i].connections[j].destindex].p.x + 870000., nodes[nodes[i].connections[j].destindex].p.y + 6820000., nodes[nodes[i].connections[j].destindex].p.z);

				Segment seg(a0, a1);
				nghb_graph_line.push_back(seg);
			}
		}
		std::vector<int> list_new_branch_id;
		std::vector<int> used_id = { 0 };
		std::vector<std::vector<int>> seen_br_id = { nodes[0].branch_id_ascend };

		// Changing branch label for gocad visualization
		for (int i = 0; i < nodes.size(); i++) {
			for (int j = 0; j < nodes[i].connections.size(); j++) {
				auto it = find(seen_br_id.begin(), seen_br_id.end(), nodes[nodes[i].connections[j].destindex].branch_id_ascend);

				if (it != seen_br_id.end()) { //id already found

					ptrdiff_t pos = it - seen_br_id.begin();
					list_new_branch_id.push_back(used_id[pos]);
					list_new_branch_id.push_back(used_id[pos]);
				}
				else {	//id not found
					seen_br_id.push_back(nodes[nodes[i].connections[j].destindex].branch_id_ascend);
					list_new_branch_id.push_back(used_id[find(seen_br_id.begin(), seen_br_id.end(), nodes[i].branch_id_ascend) - seen_br_id.begin()]);
					used_id.push_back(used_id.back() + 1);
					list_new_branch_id.push_back(used_id.back());
				}
			}
		}


		Line full_line = nghb_graph_line;

		int number_segs = full_line.get_nb_segs();
		std::vector<std::string> property_names = { "cost","equivalent_radius", "branch_id","vadose_flag" };
		std::vector<std::vector<std::vector<float>>> properties;

		properties.resize(number_segs);
		for (int segi = 0; segi < int(properties.size()); segi++) {
			properties[segi].resize(4); // there is 4 properties in total
		}
		for (int n = 0; n < number_segs; n++) {
			properties[n][0].push_back(list_cost[2 * n]);
			properties[n][0].push_back(list_cost[2 * n + 1]);

			properties[n][1].push_back(list_eq_radius[2 * n]);
			properties[n][1].push_back(list_eq_radius[2 * n + 1]);

			properties[n][2].push_back(list_new_branch_id[2 * n]);
			properties[n][2].push_back(list_new_branch_id[2 * n + 1]);

			properties[n][3].push_back(list_vadose[2 * n]);
			properties[n][3].push_back(list_vadose[2 * n + 1]);

			//properties[n][3].push_back(list_branch_dist[2 * n]); // for debugging
			//properties[n][3].push_back(list_branch_dist[2 * n + 1]);

		}
		std::string full_name = network_name + "_karst.txt";
		std::string full_dir_name = geologicalparams.directoryname + "/outputs";

		save_line(full_name, full_dir_name, full_line, property_names, properties); // version with costs

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

	// Method to compute and populate branch_sizes attribute

	void KarsticSkeleton::compute_branch_sizes() {
		std::unordered_map<int, int> branchSizeMap;

		// Count sizes of branches
		for (const KarsticNode& node : nodes) {
			if (node.branch_id == -1) { // -1 branch id tag, which means intersection, is handled separately
				nb_of_intersections++;
			}
			else {
				++branchSizeMap[node.branch_id];
			}
		}

		// Sort branchSizeMap by keys (branch IDs)
		std::vector<std::pair<int, int>> sortedBranchSizes(branchSizeMap.begin(), branchSizeMap.end());
		std::sort(sortedBranchSizes.begin(), sortedBranchSizes.end());

		// Populate branch_sizes vector
		branch_sizes.clear();
		for (const auto& pair : sortedBranchSizes) {
			branch_sizes.push_back(pair.second);
		}
	}
}