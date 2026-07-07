/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for new methods and modifications to all original methods
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for original versions of methods of Dijkstra_Compute_Path, DijkstraGetShortestPathTo
If you use this code, please cite : Paris et al., 2021, Computer Graphic Forum.

***************************************************************/

#include "KarstNSim/graph.h"

namespace {
	struct DijkstraQueueNode {
		float distance;
		int node;
	};

	struct DijkstraQueueGreater {
		bool operator()(const DijkstraQueueNode& a, const DijkstraQueueNode& b) const {
			return a.distance > b.distance;
		}
	};
}

namespace KarstNSim {
	/*
	\brief
	*/
	CostGraph::CostGraph(int n)
	{
		adj.resize(n);
	}

	void CostGraph::NormalizeWeights() {
		float min = std::numeric_limits<float>::infinity();
		float max = 0.;
		for (int i = 0; i < adj.size(); i++) {
			for (int j = 0; j < adj[i].size(); j++) {
				if (adj[i][j].target >= 0) { // avoid empty neighbors (when n<N) 
					std::vector<float> weight_itr = adj[i][j].weight;
					for (int k = 0; k < weight_itr.size(); k++) {
						if (max < weight_itr[k]) max = weight_itr[k];
						if (min > weight_itr[k]) min = weight_itr[k];
					}
				}
			}
		}
		for (int i = 0; i < adj.size(); i++) {
			for (int j = 0; j < adj[i].size(); j++) {
				if (adj[i][j].target >= 0) { // avoid empty neighbors (when n<N) 
					for (int k = 0; k < adj[i][j].weight.size(); k++) {
						adj[i][j].weight[k] = (adj[i][j].weight[k] - min) / (max - min);
					}
				}
			}
		}
		InvalidateShortestPathPreprocessing();
	}

	void CostGraph::SetEdge(const int& s, const int& neigh, const int& t, const std::vector<float>& w, const std::vector<bool>& ff)
	{
		InvalidateShortestPathPreprocessing();
		adj[s][neigh] = GraphEdge(t, w, ff);
	}

	void CostGraph::SetEdge(const int& s, const int& neigh, const int& t, const std::vector<float>& w)
	{
		InvalidateShortestPathPreprocessing();
		adj[s][neigh] = GraphEdge(t, w);
	}

	void CostGraph::UpdateEdgeWeight(const int& s, const int& t, const std::vector<float>& w)
	{
		InvalidateShortestPathPreprocessing();
		adj[s][t].weight = w;
	}

	int CostGraph::GetIdxNeighbor(const int& s, const int& t)
	{
		for (int i = 0; i < adj[s].size(); i++) {
			if (adj[s][i].target >= 0) { // avoid empty neighbors (when n<N) 
				if (adj[s][i].target == t)
				{
					return i;
				}
			}
		}
		return 0;
	}

	void CostGraph::ClearShortestPathPreprocessing() const
	{
		reverse_adjacency_ready_ = false;
		std::vector<std::size_t>().swap(reverse_offsets_);
		std::vector<int>().swap(reverse_sources_);
		std::vector<std::uint16_t>().swap(reverse_edge_slots_);
	}

	void CostGraph::InvalidateShortestPathPreprocessing() const
	{
		if (!reverse_adjacency_ready_) {
			return;
		}

		ClearShortestPathPreprocessing();
	}

	float CostGraph::GetDirectedEdgeWeight(int source, int target, int outlet_count) const
	{
		constexpr float inf = std::numeric_limits<float>::infinity();

		if (source < 0 || target < 0 || source >= int(adj.size()) || target >= int(adj.size())) {
			return inf;
		}

		const int n2 = int(adj.row(0).size());

		for (int i = 0; i < n2; ++i) {
			const GraphEdge& edge = adj(source, i);

			if (edge.target != target) {
				continue;
			}

			if (outlet_count < 0 || outlet_count >= int(edge.weight.size())) {
				return inf;
			}

			return edge.weight[outlet_count];
		}

		return inf;
	}

	void CostGraph::BuildReverseAdjacency() const
	{
		const int n = int(adj.size());

		if (n == 0) {
			return;
		}

		if (reverse_adjacency_ready_ && reverse_offsets_.size() == std::size_t(n + 1)) {
			return;
		}

		const int n2 = int(adj.row(0).size());

		if (n2 > int(std::numeric_limits<std::uint16_t>::max())) {
			throw std::runtime_error(
				"[dijkstra] Cannot build reverse adjacency: the local edge slot exceeds uint16_t capacity."
			);
		}

		reverse_adjacency_ready_ = false;
		reverse_offsets_.assign(std::size_t(n + 1), 0);

		for (int u = 0; u < n; ++u) {
			for (int k = 0; k < n2; ++k) {
				const int v = adj(u, k).target;

				if (v >= 0) {
					++reverse_offsets_[std::size_t(v + 1)];
				}
			}
		}

		for (int i = 1; i <= n; ++i) {
			reverse_offsets_[std::size_t(i)] += reverse_offsets_[std::size_t(i - 1)];
		}

		const std::size_t edge_count = reverse_offsets_[std::size_t(n)];
		reverse_sources_.assign(edge_count, -1);
		reverse_edge_slots_.assign(edge_count, 0);

		std::vector<std::size_t> cursor = reverse_offsets_;

		for (int u = 0; u < n; ++u) {
			for (int k = 0; k < n2; ++k) {
				const int v = adj(u, k).target;

				if (v < 0) {
					continue;
				}

				const std::size_t pos = cursor[std::size_t(v)]++;
				reverse_sources_[pos] = u;
				reverse_edge_slots_[pos] = static_cast<std::uint16_t>(k);
			}
		}

		reverse_adjacency_ready_ = true;
	}

	void CostGraph::DijkstraComputePathsBidirectional(int outlet_count, int source, int target, std::vector<float>& distance, std::vector<int>& previous) const
	{
		constexpr float inf = std::numeric_limits<float>::infinity();

		const int n = int(adj.size());

		if (n == 0) {
			return;
		}

		distance.assign(std::size_t(n), inf);
		previous.assign(std::size_t(n), -1);

		if (source < 0 || target < 0 || source >= n || target >= n) {
			return;
		}

		if (source == target) {
			distance[std::size_t(source)] = 0.0f;
			return;
		}

		BuildReverseAdjacency();

		const int n2 = int(adj.row(0).size());

		std::vector<float> dist_forward(std::size_t(n), inf);
		std::vector<float> dist_backward(std::size_t(n), inf);
		std::vector<int> previous_forward(std::size_t(n), -1);
		std::vector<int> next_backward(std::size_t(n), -1);
		std::vector<char> settled_forward(std::size_t(n), 0);
		std::vector<char> settled_backward(std::size_t(n), 0);

		std::priority_queue<DijkstraQueueNode, std::vector<DijkstraQueueNode>, DijkstraQueueGreater> queue_forward;
		std::priority_queue<DijkstraQueueNode, std::vector<DijkstraQueueNode>, DijkstraQueueGreater> queue_backward;

		dist_forward[std::size_t(source)] = 0.0f;
		dist_backward[std::size_t(target)] = 0.0f;

		queue_forward.push({ 0.0f, source });
		queue_backward.push({ 0.0f, target });

		float best_distance = inf;
		int meeting_node = -1;

		while (!queue_forward.empty() && !queue_backward.empty()) {
			const float top_forward = queue_forward.top().distance;
			const float top_backward = queue_backward.top().distance;

			if (std::isfinite(best_distance) && top_forward + top_backward >= best_distance) {
				break;
			}

			if (top_forward <= top_backward) {
				const DijkstraQueueNode current = queue_forward.top();
				queue_forward.pop();

				const float dist = current.distance;
				const int u = current.node;

				if (dist > dist_forward[std::size_t(u)] || settled_forward[std::size_t(u)]) {
					continue;
				}

				settled_forward[std::size_t(u)] = 1;

				if (std::isfinite(dist_backward[std::size_t(u)]) &&
					dist + dist_backward[std::size_t(u)] < best_distance) {
					best_distance = dist + dist_backward[std::size_t(u)];
					meeting_node = u;
				}

				for (int i = 0; i < n2; ++i) {
					const GraphEdge& edge = adj(u, i);
					const int v = edge.target;

					if (v < 0) {
						continue;
					}

					if (outlet_count < 0 || outlet_count >= int(edge.weight.size())) {
						continue;
					}

					const float weight = edge.weight[outlet_count];

					if (!std::isfinite(weight) || weight < 0.0f) {
						continue;
					}

					const float alt = dist + weight;

					if (alt < dist_forward[std::size_t(v)]) {
						dist_forward[std::size_t(v)] = alt;
						previous_forward[std::size_t(v)] = u;
						queue_forward.push({ alt, v });

						if (std::isfinite(dist_backward[std::size_t(v)]) &&
							alt + dist_backward[std::size_t(v)] < best_distance) {
							best_distance = alt + dist_backward[std::size_t(v)];
							meeting_node = v;
						}
					}
				}
			}
			else {
				const DijkstraQueueNode current = queue_backward.top();
				queue_backward.pop();

				const float dist = current.distance;
				const int u = current.node;

				if (dist > dist_backward[std::size_t(u)] || settled_backward[std::size_t(u)]) {
					continue;
				}

				settled_backward[std::size_t(u)] = 1;

				if (std::isfinite(dist_forward[std::size_t(u)]) &&
					dist_forward[std::size_t(u)] + dist < best_distance) {
					best_distance = dist_forward[std::size_t(u)] + dist;
					meeting_node = u;
				}

				for (std::size_t pos = reverse_offsets_[std::size_t(u)];
					pos < reverse_offsets_[std::size_t(u + 1)];
					++pos) {

					const int pred = reverse_sources_[pos];
					const int edge_slot = int(reverse_edge_slots_[pos]);
					const GraphEdge& reverse_edge = adj(pred, edge_slot);

					if (outlet_count < 0 || outlet_count >= int(reverse_edge.weight.size())) {
						continue;
					}

					const float weight = reverse_edge.weight[outlet_count];

					if (!std::isfinite(weight) || weight < 0.0f) {
						continue;
					}

					const float alt = dist + weight;

					if (alt < dist_backward[std::size_t(pred)]) {
						dist_backward[std::size_t(pred)] = alt;
						next_backward[std::size_t(pred)] = u;
						queue_backward.push({ alt, pred });

						if (std::isfinite(dist_forward[std::size_t(pred)]) &&
							dist_forward[std::size_t(pred)] + alt < best_distance) {
							best_distance = dist_forward[std::size_t(pred)] + alt;
							meeting_node = pred;
						}
					}
				}
			}
		}

		if (meeting_node < 0 || !std::isfinite(best_distance)) {
			return;
		}

		std::vector<int> left_path;

		for (int u = meeting_node; u != -1; u = previous_forward[std::size_t(u)]) {
			left_path.push_back(u);
		}

		std::reverse(left_path.begin(), left_path.end());

		std::vector<int> full_path = left_path;

		for (int u = next_backward[std::size_t(meeting_node)]; u != -1; u = next_backward[std::size_t(u)]) {
			full_path.push_back(u);

			if (u == target) {
				break;
			}

			if (full_path.size() > std::size_t(n)) {
				throw std::runtime_error("[dijkstra] Invalid bidirectional path reconstruction.");
			}
		}

		if (full_path.empty() || full_path.front() != source || full_path.back() != target) {
			return;
		}

		float cumulative_distance = 0.0f;
		distance[std::size_t(source)] = 0.0f;
		previous[std::size_t(source)] = -1;

		for (std::size_t i = 1; i < full_path.size(); ++i) {
			const int u = full_path[i - 1];
			const int v = full_path[i];
			const float weight = GetDirectedEdgeWeight(u, v, outlet_count);

			if (!std::isfinite(weight)) {
				throw std::runtime_error("[dijkstra] Could not recover a directed edge while reconstructing a bidirectional path.");
			}

			cumulative_distance += weight;
			distance[std::size_t(v)] = cumulative_distance;
			previous[std::size_t(v)] = u;
		}
	}

	void CostGraph::DijkstraComputePaths(int outlet_count, int source, std::vector<float>& distance, std::vector<int>& previous, int target) const
	{
		constexpr float max_weight = std::numeric_limits<float>::infinity();
		size_t n = adj.size();
		size_t n2 = adj.row(0).size();
		distance.resize(n, max_weight);
		distance[source] = 0;
		previous.resize(n, -1);

		std::set<std::pair<float, int>> vertex_queue;
		vertex_queue.insert(std::make_pair(distance[source], source));

		while (!vertex_queue.empty())
		{
			float dist = vertex_queue.begin()->first;
			int u = vertex_queue.begin()->second;
			vertex_queue.erase(vertex_queue.begin());

			if (u == target) {
				//distance[target] = dist; // Ensure target's distance is updated
				break; // Found shortest path to target, no need to explore further
			}

			for (int i = 0; i < n2; ++i)
			{
				if (adj(u, i).target >= 0) { // avoid empty neighbors (when n<N) 
					const GraphEdge& neighbor(adj(u, i));
					int v = neighbor.target;
					float weight = neighbor.weight[outlet_count];
					float distance_through_u = dist + weight;
					if (distance_through_u < distance[v])
					{
						vertex_queue.erase(std::make_pair(distance[v], v));
						distance[v] = distance_through_u;
						previous[v] = u;
						vertex_queue.insert(std::make_pair(distance[v], v));
					}
				}
			}
		}
	}

	
	void CostGraph::DijkstraComputePathsSurfaceBidirectional(
		int outlet_count,
		int source,
		int& reach,
		std::vector<float>& distance,
		std::vector<int>& previous,
		const Array2D<char>& samples_surf_flags,
		int target,
		bool& already_reached) const
	{
		constexpr float inf = std::numeric_limits<float>::infinity();

		const int n = int(adj.size());

		reach = -1;
		already_reached = false;

		if (n == 0) {
			return;
		}

		distance.assign(std::size_t(n), inf);
		previous.assign(std::size_t(n), -1);

		if (source < 0 || source >= n || outlet_count < 0) {
			return;
		}

		if (samples_surf_flags(source, outlet_count)) {
			reach = source;
			distance[std::size_t(source)] = 0.0f;
			already_reached = (source == target);
			return;
		}

		BuildReverseAdjacency();

		const int n2 = int(adj.row(0).size());

		std::vector<float> dist_forward(std::size_t(n), inf);
		std::vector<float> dist_backward(std::size_t(n), inf);
		std::vector<int> previous_forward(std::size_t(n), -1);
		std::vector<int> next_backward(std::size_t(n), -1);
		std::vector<int> surface_root(std::size_t(n), -1);
		std::vector<char> settled_forward(std::size_t(n), 0);
		std::vector<char> settled_backward(std::size_t(n), 0);

		std::priority_queue<DijkstraQueueNode, std::vector<DijkstraQueueNode>, DijkstraQueueGreater> queue_forward;
		std::priority_queue<DijkstraQueueNode, std::vector<DijkstraQueueNode>, DijkstraQueueGreater> queue_backward;

		dist_forward[std::size_t(source)] = 0.0f;
		queue_forward.push({ 0.0f, source });

		int nb_surface_nodes = 0;

		for (int node = 0; node < n; ++node) {
			if (!samples_surf_flags(node, outlet_count)) {
				continue;
			}

			dist_backward[std::size_t(node)] = 0.0f;
			surface_root[std::size_t(node)] = node;
			queue_backward.push({ 0.0f, node });
			++nb_surface_nodes;
		}

		if (nb_surface_nodes == 0) {
			throw std::runtime_error(
				"[dijkstra] Cannot compute a path to surface: no node is flagged on the requested surface."
			);
		}

		float best_distance = inf;
		int meeting_node = -1;
		int best_surface_node = -1;

		while (!queue_forward.empty() && !queue_backward.empty()) {
			const float top_forward = queue_forward.top().distance;
			const float top_backward = queue_backward.top().distance;

			if (std::isfinite(best_distance) && top_forward + top_backward >= best_distance) {
				break;
			}

			// Prefer forward expansion when the backward queue is large because it is
			// initialized with all surface nodes. This keeps the method exact while
			// avoiding an unnecessary full expansion of the target surface.
			const bool expand_forward =
				(queue_backward.size() > queue_forward.size()) ||
				(top_forward <= top_backward);

			if (expand_forward) {
				const DijkstraQueueNode current = queue_forward.top();
				queue_forward.pop();

				const float dist = current.distance;
				const int u = current.node;

				if (dist > dist_forward[std::size_t(u)] || settled_forward[std::size_t(u)]) {
					continue;
				}

				settled_forward[std::size_t(u)] = 1;

				if (std::isfinite(dist_backward[std::size_t(u)]) &&
					dist + dist_backward[std::size_t(u)] < best_distance) {
					best_distance = dist + dist_backward[std::size_t(u)];
					meeting_node = u;
					best_surface_node = surface_root[std::size_t(u)];
				}

				for (int i = 0; i < n2; ++i) {
					const GraphEdge& edge = adj(u, i);
					const int v = edge.target;

					if (v < 0) {
						continue;
					}

					if (outlet_count >= int(edge.weight.size())) {
						continue;
					}

					const float weight = edge.weight[outlet_count];

					if (!std::isfinite(weight) || weight < 0.0f) {
						continue;
					}

					const float alt = dist + weight;

					if (alt < dist_forward[std::size_t(v)]) {
						dist_forward[std::size_t(v)] = alt;
						previous_forward[std::size_t(v)] = u;
						queue_forward.push({ alt, v });

						if (std::isfinite(dist_backward[std::size_t(v)]) &&
							alt + dist_backward[std::size_t(v)] < best_distance) {
							best_distance = alt + dist_backward[std::size_t(v)];
							meeting_node = v;
							best_surface_node = surface_root[std::size_t(v)];
						}
					}
				}
			}
			else {
				const DijkstraQueueNode current = queue_backward.top();
				queue_backward.pop();

				const float dist = current.distance;
				const int u = current.node;

				if (dist > dist_backward[std::size_t(u)] || settled_backward[std::size_t(u)]) {
					continue;
				}

				settled_backward[std::size_t(u)] = 1;

				if (std::isfinite(dist_forward[std::size_t(u)]) &&
					dist_forward[std::size_t(u)] + dist < best_distance) {
					best_distance = dist_forward[std::size_t(u)] + dist;
					meeting_node = u;
					best_surface_node = surface_root[std::size_t(u)];
				}

				for (std::size_t pos = reverse_offsets_[std::size_t(u)];
					pos < reverse_offsets_[std::size_t(u + 1)];
					++pos) {

					const int pred = reverse_sources_[pos];
					const int edge_slot = int(reverse_edge_slots_[pos]);
					const GraphEdge& reverse_edge = adj(pred, edge_slot);

					if (outlet_count >= int(reverse_edge.weight.size())) {
						continue;
					}

					const float weight = reverse_edge.weight[outlet_count];

					if (!std::isfinite(weight) || weight < 0.0f) {
						continue;
					}

					const float alt = dist + weight;

					if (alt < dist_backward[std::size_t(pred)]) {
						dist_backward[std::size_t(pred)] = alt;
						next_backward[std::size_t(pred)] = u;
						surface_root[std::size_t(pred)] = surface_root[std::size_t(u)];
						queue_backward.push({ alt, pred });

						if (std::isfinite(dist_forward[std::size_t(pred)]) &&
							dist_forward[std::size_t(pred)] + alt < best_distance) {
							best_distance = dist_forward[std::size_t(pred)] + alt;
							meeting_node = pred;
							best_surface_node = surface_root[std::size_t(pred)];
						}
					}
				}
			}
		}

		if (meeting_node < 0 || best_surface_node < 0 || !std::isfinite(best_distance)) {
			return;
		}

		std::vector<int> left_path;

		for (int u = meeting_node; u != -1; u = previous_forward[std::size_t(u)]) {
			left_path.push_back(u);
		}

		std::reverse(left_path.begin(), left_path.end());

		std::vector<int> full_path = left_path;

		for (int u = next_backward[std::size_t(meeting_node)]; u != -1; u = next_backward[std::size_t(u)]) {
			full_path.push_back(u);

			if (u == best_surface_node) {
				break;
			}

			if (full_path.size() > std::size_t(n)) {
				throw std::runtime_error("[dijkstra] Invalid bidirectional surface path reconstruction.");
			}
		}

		if (full_path.empty() || full_path.front() != source || full_path.back() != best_surface_node) {
			return;
		}

		reach = best_surface_node;
		already_reached = (reach == target);

		float cumulative_distance = 0.0f;
		distance[std::size_t(source)] = 0.0f;
		previous[std::size_t(source)] = -1;

		for (std::size_t i = 1; i < full_path.size(); ++i) {
			const int u = full_path[i - 1];
			const int v = full_path[i];
			const float weight = GetDirectedEdgeWeight(u, v, outlet_count);

			if (!std::isfinite(weight)) {
				throw std::runtime_error(
					"[dijkstra] Could not recover a directed edge while reconstructing a bidirectional surface path."
				);
			}

			cumulative_distance += weight;
			distance[std::size_t(v)] = cumulative_distance;
			previous[std::size_t(v)] = u;
		}
	}

	// This version computes a path between a point and a surface with Dijkstra algorithm.
	void CostGraph::DijkstraComputePathsSurface(int outlet_count, int source, int& reach, std::vector<float>& distance, std::vector<int>& previous, const  Array2D<char>& samples_surf_flags, int target, bool &already_reached) const
	{
		constexpr float max_weight = std::numeric_limits<float>::infinity();
		size_t n = adj.size();
		size_t n2 = adj.row(0).size();
		distance.resize(n, max_weight);
		distance[source] = 0;
		previous.resize(n, -1);

		std::set<std::pair<float, int>> vertex_queue;
		vertex_queue.insert(std::make_pair(distance[source], source));
		while (!vertex_queue.empty())
		{

			float dist = vertex_queue.begin()->first;
			int u = vertex_queue.begin()->second;
			vertex_queue.erase(vertex_queue.begin());

			if (samples_surf_flags(u, outlet_count)) {
				reach = u;
				if (u == target) {
					already_reached = true;
				}
				//distance[target] = dist; // Ensure target's distance is updated
				break;
			}

			// Visit each edge exiting u
			for (int i = 0; i < n2; ++i)
			{
				if (adj(u, i).target >= 0) { // avoid empty neighbors (when n<N)

					const GraphEdge& neighbor(adj(u, i));
					int v = neighbor.target;
					float weight = neighbor.weight[outlet_count];
					float distance_through_u = dist + weight;
					if (distance_through_u < distance[v])
					{
						vertex_queue.erase(std::make_pair(distance[v], v));
						distance[v] = distance_through_u;
						previous[v] = u;
						vertex_queue.insert(std::make_pair(distance[v], v));
					}
				}
			}
		}
	}

	std::pair<std::vector<int>, std::vector<float>> CostGraph::DijkstraGetShortestPathTo(int target, const std::vector<int>& previous, const std::vector<float>& distances, float& dist) const
	{
		std::vector<int> path;
		std::vector<float> path_cost;
		int n = target;
		dist = distances[n];
		float dist1 = dist;
		float dist2 = dist;
		float eps = 1e-10f;
		int tocopy = 0;
		for (; n != -1; n = previous[n])
		{
			bool prev = true;
			dist1 = dist2;
			dist2 = distances[n];
			path.insert(path.begin(), n);
			if (abs(dist1 - dist2) < eps) {
				prev = false;
				tocopy++;
			}
			if (prev) {
				path_cost.insert(path_cost.begin(), dist1 - dist2);
				if (tocopy > 0) {
					for (int copy = 0; copy < tocopy; copy++) {
						path_cost.insert(path_cost.begin(), dist1 - dist2);
					}
					tocopy = 0;
				}
			}
		}
		if (!path.empty()) {
			const std::size_t desired = path.size();
			if (path_cost.size() < desired) {
				path_cost.insert(path_cost.begin(), desired - path_cost.size(), 0.f);
			}
		}
		return std::make_pair(path, path_cost);
	}

}