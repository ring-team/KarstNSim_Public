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
	}

	void CostGraph::SetEdge(const int& s, const int& neigh, const int& t, const std::vector<float>& w, const std::vector<bool>& ff)
	{
		adj[s][neigh] = GraphEdge(t, w, ff);
	}

	void CostGraph::SetEdge(const int& s, const int& neigh, const int& t, const std::vector<float>& w)
	{
		adj[s][neigh] = GraphEdge(t, w);
	}

	void CostGraph::UpdateEdgeWeight(const int& s, const int& t, const std::vector<float>& w)
	{
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
		return std::make_pair(path, path_cost);
	}

}