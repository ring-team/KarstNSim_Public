/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for new methods and modifications to all original methods
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for original versions of methods of Dijkstra_Compute_Path, DijkstraGetShortestPathTo
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.
Find the corresponding code at https://github.com/aparis69/Karst-Synthesis

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
		double min = std::numeric_limits<double>::infinity();
		double max = 0.;
		for (int i = 0; i < adj.size(); i++) {
			for (int j = 0; j < adj[i].size(); j++) {
				std::vector<double> weight_itr = adj[i][j].weight;
				for (int k = 0; k < weight_itr.size(); k++) {
					if (max < weight_itr[k]) max = weight_itr[k];
					if (min > weight_itr[k]) min = weight_itr[k];
				}
			}
		}
		for (int i = 0; i < adj.size(); i++) {
			for (int j = 0; j < adj[i].size(); j++) {
				for (int k = 0; k < adj[i][j].weight.size(); k++) {
					adj[i][j].weight[k] = (adj[i][j].weight[k] - min) / (max - min);
				}
			}
		}
	}

	/*
	\brief
	*/
	void CostGraph::SetEdge(int s, int t, std::vector<double> w, std::vector<bool> ff)
	{
		adj[s].push_back(GraphEdge(t, w, ff));
	}

	void CostGraph::SetEdge(int s, int t, std::vector<double> w)
	{
		adj[s].push_back(GraphEdge(t, w));
	}

	void CostGraph::UpdateEdgeWeight(int s, int t, std::vector<double> w)
	{
		adj[s][t].weight = w;
	}

	int CostGraph::GetIdxNeighbor(int s, int t)
	{
		for (int i = 0; i < adj[s].size(); i++) {
			if (adj[s][i].target == t)
			{
				return i;
			}
		}
		return 0;
	}

	/*
	\brief
	*/
	void CostGraph::DijkstraComputePaths(int outlet_count, int source, std::vector<double>& distance, std::vector<int>& previous) const
	{
		constexpr double max_weight = std::numeric_limits<double>::infinity();
		size_t n = adj.size();
		distance.clear();
		distance.resize(n, max_weight);
		distance[source] = 0;
		//int nb_frac_families = int(adj[0][0].frac_flag.size());
		previous.clear();
		previous.resize(n, -1);
		std::set<std::pair<double, int>> vertex_queue;
		vertex_queue.insert(std::make_pair(distance[source], source));
		while (!vertex_queue.empty())
		{
			double dist = vertex_queue.begin()->first;
			int u = vertex_queue.begin()->second;
			vertex_queue.erase(vertex_queue.begin());

			// Visit each edge exiting u
			const std::vector<GraphEdge>& neighbors = adj[u];
			for (std::vector<GraphEdge>::const_iterator neighbor_iter = neighbors.begin(); neighbor_iter != neighbors.end(); neighbor_iter++)
			{
				//double cumulated_length = 0.;
				int v = neighbor_iter->target;
				double weight = neighbor_iter->weight[outlet_count];
				double distance_through_u = dist + weight;
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

	/*
	\brief This version doesn't include cost reduction of the edges when they are created, in order to mimic almost simultaneous processes when one conduit doesn't influence other ones nearby
	*/
	void CostGraph::DijkstraComputePathsNoCostRed(int outlet_count, int source, std::vector<double>& distance, std::vector<int>& previous) const
	{
		constexpr double max_weight = std::numeric_limits<double>::infinity();
		size_t n = adj_no_cost_red.size();
		distance.clear();
		distance.resize(n, max_weight);
		distance[source] = 0;
		//int nb_frac_families = int(adj[0][0].frac_flag.size());
		previous.clear();
		previous.resize(n, -1);
		std::set<std::pair<double, int>> vertex_queue;
		vertex_queue.insert(std::make_pair(distance[source], source));
		while (!vertex_queue.empty())
		{
			double dist = vertex_queue.begin()->first;
			int u = vertex_queue.begin()->second;
			vertex_queue.erase(vertex_queue.begin());

			// Visit each edge exiting u
			const std::vector<GraphEdge>& neighbors = adj_no_cost_red[u];
			for (std::vector<GraphEdge>::const_iterator neighbor_iter = neighbors.begin(); neighbor_iter != neighbors.end(); neighbor_iter++)
			{
				//double cumulated_length = 0.;
				int v = neighbor_iter->target;
				double weight = neighbor_iter->weight[outlet_count];
				double distance_through_u = dist + weight;
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


	// This version computes a path between a point and a surface with Dijkstra algorithm.

	void CostGraph::DijkstraComputePathsSurface(int outlet_count, int source, int& reach, std::vector<double>& distance, std::vector<int>& previous, std::vector<std::vector<bool>> samples_surf_flags, int spring_idx, int target, bool &already_reached) const
	{

		constexpr double max_weight = std::numeric_limits<double>::infinity();
		size_t n = adj.size();
		distance.clear();
		distance.resize(n, max_weight);
		distance[source] = 0;
		previous.clear();
		previous.resize(n, -1);
		std::set<std::pair<double, int>> vertex_queue;
		vertex_queue.insert(std::make_pair(distance[source], source));
		bool watertable_was_reached = false;
		while ((!vertex_queue.empty()) && (!watertable_was_reached))
		{
			double dist = vertex_queue.begin()->first;
			int u = vertex_queue.begin()->second;
			vertex_queue.erase(vertex_queue.begin());

			// Visit each edge exiting u
			const std::vector<GraphEdge>& neighbors = adj[u];
			for (std::vector<GraphEdge>::const_iterator neighbor_iter = neighbors.begin(); neighbor_iter != neighbors.end(); neighbor_iter++)
			{
				int v = neighbor_iter->target;
				double weight = neighbor_iter->weight[outlet_count];
				double distance_through_u = dist + weight;
				if (distance_through_u < distance[v])
				{

					vertex_queue.erase(std::make_pair(distance[v], v));

					distance[v] = distance_through_u;
					previous[v] = u;
					vertex_queue.insert(std::make_pair(distance[v], v));

					if (samples_surf_flags[v][spring_idx]) // If we reach a node below the water table, he have found the shortest path in the vadose zone.
					{
						reach = v;
						watertable_was_reached = true;
					}
					if (v == target) {
						reach = v;
						watertable_was_reached = true;
						already_reached = true;
					}
				}
			}
		}
	}

	/*
	\brief
	*/
	std::pair<std::vector<int>, std::vector<double>> CostGraph::DijkstraGetShortestPathTo(int target, const std::vector<int>& previous, const std::vector<double>& distances, double& dist) const
	{
		std::vector<int> path;
		std::vector<double> path_cost;
		int n = target;
		dist = distances[n];
		double dist1 = dist;
		double dist2 = dist;
		double eps = 1e-4;
		int tocopy = 0;
		for (; n != -1; n = previous[n])
		{
			bool prev = true;
			dist1 = dist2;
			dist2 = distances[n];
			path.insert(path.begin(), n);
			if (abs(dist1 - dist2) < eps) {
				prev = false;
				tocopy ++;
			}
			if (prev) {
				path_cost.insert(path_cost.begin(), dist1 - dist2);
				if (tocopy>0) {
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
