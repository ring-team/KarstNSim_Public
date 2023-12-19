/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for new methods and modifications to original methods
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for original versions of KarsticSkeleton constructor, Amplify and AppendPaths 
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.
Find the corresponding code at https://github.com/aparis69/Karst-Synthesis

***************************************************************/

#include "KarstNSim/graph.h"

/*!
\brief Constructor from the volumetric nearest neighbour graph and a set of path.
\param graph nearest neighbour graph
\param paths all paths in the skeleton
*/
namespace KarstNSim {
	KarsticSkeleton::KarsticSkeleton(const GraphOperations* graph, const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>> costs, const std::vector<std::vector<bool>> vadose)
	{
		AppendPaths(graph, paths, costs, vadose);
	}

	int KarsticSkeleton::count_vadose_nodes()
	{
		int count = 0;
		for (int i = 0; i < nodes.size(); i++) {
			if (nodes[i].vadose) count++; // if node is a vadose node
		}
		return count;
	}

	bool KarsticSkeleton::edge_exists(int u, int v) {
		return std::find(nodes[u].connections.begin(), nodes[u].connections.end(), v) != nodes[u].connections.end();
	}

	/*!
	\brief Amplify an existing network from a new set of key points.
	Calling this method with an empty parameter will procedurally place ten key points in the scene, and amplify the skeleton in the same way.
	\param keypts new key points.
	*/
	void KarsticSkeleton::Amplify_vadose(GraphOperations* graph, GeologicalParameters params)
	{
		int vadose_nodes = count_vadose_nodes();
		int max_loops = int(params.loop_density_vadose * vadose_nodes); // number of loops is number of vadose nodes times the density fraction defined by the user

		std::vector<std::vector<int>> newPaths;
		std::vector<std::vector<double>> newPaths_cost;
		std::vector<std::vector<bool>> newPaths_vadose;

		while (max_loops > 0)
		{
			int u = int(generateRandomDouble(0, double(nodes.size()))); // Select two random vertices of the main graph
			int v = int(generateRandomDouble(0, double(nodes.size())));
			int int_u_index = nodes[u].index;
			int int_v_index = nodes[v].index;

			// Check if edge already exists or u and v are the same vertex
			if (u != v && !edge_exists(u, v) && nodes[u].vadose && nodes[v].vadose) {
				double distance = magnitude(nodes[u].p - nodes[v].p);

				if (distance <= params.max_dist_loops_vadose) { // also, check if the euclidian distance separating them is less than the threshold value defined by the user

					std::pair< std::vector<int>, std::vector<double>> pair = graph->shortest_path(false, int_u_index, int_v_index); // if all the previous tests were ok, then we can compute a path between those two nodes
					std::vector<int> bestPath = pair.first;
					std::vector<double> bestPath_cost = pair.second;

					newPaths.push_back(std::vector<int>(bestPath.begin(), bestPath.end()));
					newPaths_cost.push_back(std::vector<double>(bestPath_cost.begin(), bestPath_cost.end()));
					newPaths_vadose.push_back(std::vector<bool>(bestPath.size(), true)); // all paths are vadose paths here
					--max_loops; // and then decrement the number of loops still to be generated

				}
			}
		}

		// Update the internal karstic skeleton struture
		AppendPaths(graph, newPaths, newPaths_cost, newPaths_vadose);

	}

	/*!
	\brief Amplify an existing network from a new set of key points.
	Calling this method with an empty parameter will procedurally place ten key points in the scene, and amplify the skeleton in the same way.
	\param keypts new key points.
	*/
	void KarsticSkeleton::Amplify_phreatic(GraphOperations* graph, GeologicalParameters params)
	{

		int vadose_nodes = count_vadose_nodes();
		int max_loops = int(params.loop_density_phreatic * (nodes.size() - vadose_nodes)); // number of loops is number of phreatic nodes times the density fraction defined by the user

		std::vector<std::vector<int>> newPaths;
		std::vector<std::vector<double>> newPaths_cost;
		std::vector<std::vector<bool>> newPaths_vadose;

		while (max_loops > 0)
		{
			int u = int(generateRandomDouble(0, double(nodes.size()))); // Select two random vertices of the main graph
			int v = int(generateRandomDouble(0, double(nodes.size())));
			int int_u_index = nodes[u].index; // Store those vertices' indices in the volumetric graph
			int int_v_index = nodes[v].index;

			// Check if edge already exists or u and v are the same vertex
			if (u != v && !edge_exists(u, v) && !nodes[u].vadose && !nodes[v].vadose) {
				double distance = magnitude(nodes[u].p - nodes[v].p);
				if (distance <= params.max_dist_loops_phreatic) { // also, check if the euclidian distance separating them is less than the threshold value defined by the user

					std::pair< std::vector<int>, std::vector<double>> pair = graph->shortest_path(true, int_u_index, int_v_index); // if all the previous tests were ok, then we can compute a path between those two nodes
					std::vector<int> bestPath = pair.first;
					std::vector<double> bestPath_cost = pair.second;
					newPaths.push_back(std::vector<int>(bestPath.begin(), bestPath.end()));
					newPaths_cost.push_back(std::vector<double>(bestPath_cost.begin(), bestPath_cost.end()));
					newPaths_vadose.push_back(std::vector<bool>(bestPath.size(), false)); // all paths are phreatic paths here
					--max_loops; // and then decrement the number of loops still to be generated
				}
			}
		}

		// Update the internal karstic skeleton struture
		AppendPaths(graph, newPaths, newPaths_cost, newPaths_vadose);

		//// if keypts is empty, generate some procedurally
		//if (newkeypts.empty())
		//	newkeypts = GenerateKeyPoints();

		//// First, we must add the new key points to the graph and compute their neighbours
		//graph->AddNewSamples(newkeypts);

		//// Then, we link the new samples to the existing network
		//auto newPaths = graph->AmplifyKarsticSkeleton(nodes, newkeypts);

		//// Update the internal karstic skeleton struture
		//AppendPaths(graph, newPaths);
	}

	/*
	\brief Generates a set of key points inside karstic skeleton domain.
	Key points can be interior waypoints, or deadends
	*/
	//std::vector<KeyPoint> KarsticSkeleton::GenerateKeyPoints() const
	//{
	//	// Compute bounding box of the skeleton
	//	Box box = GetBox();
	//
	//	// Compute random points inside the box
	//	std::vector<KeyPoint> keypts;
	//	for (int i = 0; i < 10; i++)
	//	{
	//		Vector3 p = box.RandomInside();
	//		KeyPointType type = Random::Uniform() < 0.5 ? KeyPointType::Deadend : KeyPointType::Waypoint;
	//		keypts.push_back(KeyPoint(p, type));
	//	}
	//	return keypts;
	//}

	/*
	\brief Generates a set of key points inside karstic skeleton domain.
	Key points can be interior waypoints, or deadends
	*/
	//std::vector<KeyPoint> KarsticSkeleton::GenerateKeyPoints_adapted(int nb_deadend_pts, int nb_interior_nodes) const
	//{
	//	// Compute bounding box of the skeleton
	//	Box box = GetBox();
	//
	//	// Compute random points inside the box
	//	std::vector<KeyPoint> keypts;
	//	/*
	//	for (int i = 0; i < 10; i++)
	//	{
	//		Vector3 p = box.RandomInside();
	//		KeyPointType type = Random::Uniform() < 0.5 ? KeyPointType::Deadend : KeyPointType::Waypoint;
	//		keypts.push_back(KeyPoint(p, type));
	//	}
	//	*/
	//	for (int i = 0; i < nb_deadend_pts + nb_interior_nodes; i++) {
	//		Vector3 p = box.RandomInside();
	//		if (i < nb_deadend_pts) {
	//			keypts.push_back(KeyPoint(p, KeyPointType::Deadend));
	//		}
	//		else {
	//			keypts.push_back(KeyPoint(p, KeyPointType::Waypoint));
	//		}
	//		
	//	}
	//	return keypts;
	//}


	///*!
	//\brief Computes the bounding box of the skeleton.
	//*/
	//Box KarsticSkeleton::GetBox() const
	//{
	//	std::vector<Vector3> p;
	//	p.reserve(nodes.size());
	//	for (int i = 0; i < nodes.size(); i++)
	//		p.push_back(nodes[i].p);
	//	return box(p);
	//}

	/*!
	\brief Utility function for finding the index of a graph node in the karstic skeleton array.
	\param nodeIndex index in the graph (index of the sample)
	*/
	int KarsticSkeleton::GetInternalIndex(int nodeIndex) const
	{
		for (int i = 0; i < nodes.size(); i++)
		{
			if (nodes[i].index == nodeIndex)
				return i;
		}
		return -1;
	}

	/*!
	\brief Appends a set of paths to the existing karstic skeleton.
	*/
	void KarsticSkeleton::AppendPaths(const GraphOperations* graph, const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>> costs, const std::vector<std::vector<bool>> vadose)
	{
		for (int i = 0; i < paths.size(); i++)	// Nodes
		{
			for (int j = 0; j < paths[i].size(); j++)
			{
				int nodeIndex = paths[i][j];
				double nodeCost = costs[i][j];
				bool vadoseTest = vadose[i][j];
				auto it = std::find(nodes.begin(), nodes.end(), nodeIndex);
				if (it == nodes.end())
				{
					Vector3 p = graph->GetSample(nodeIndex);
					nodes.push_back({ nodeIndex, p ,nodeCost, vadoseTest });
				}
			}
		}
		for (int i = 0; i < nodes.size(); i++)	// Edges
		{
			int nodeIndex = nodes[i].index;
			std::vector<int> edges;
			for (int j = 0; j < paths.size(); j++)
			{
				for (int k = 0; k < paths[j].size(); k++)
				{
					int n = paths[j][k];
					if (nodeIndex != n)
						continue;
					if (k < paths[j].size() - 1)
					{
						int nAfter = paths[j][k + 1];
						int neiIndex = GetInternalIndex(nAfter);
						if (std::find(edges.begin(), edges.end(), neiIndex) == edges.end())
							edges.push_back(neiIndex);
					}
				}
			}
			for (int j = 0; j < edges.size(); j++)
			{
				auto newSection = KarsticConnection({ edges[j] });
				if (std::find(nodes[i].connections.begin(), nodes[i].connections.end(), edges[j]) == nodes[i].connections.end())
					nodes[i].connections.push_back(newSection);
			}
		}
	}

	/*!
	\brief Export the karstic skeleton using a simple format.
	\param file filename
	*/
	void KarsticSkeleton::save(const std::string& file, const std::string& save_directory) const
	{
		std::ofstream out;
		out.open(save_directory + '/' + file + "_nodes.dat");
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

		out.open(save_directory + '/' + file + "_links.dat");
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
				int k = nodes[i].connections[j].destIndex;
				out << indexA << " " << (k + 1) << "\n";
			}
		}
		out.flush();
		out.close();
	}

	/*!
	\brief Export the karstic skeleton using a simple format.
	*/
	void KarsticSkeleton::CreateLine(const GeologicalParameters& geologicalparams, std::string network_name) const
	{

		std::vector<Segment> nghb_graph_line;
		std::vector<double> list_cost;
		for (int i = 0; i < nodes.size(); i++) {
			for (int j = 0; j < nodes[i].connections.size(); j++) {

				list_cost.push_back(double(nodes[i].cost));
				list_cost.push_back(double(nodes[nodes[i].connections[j].destIndex].cost));

				Vector3 a0 = nodes[i].p;
				Vector3 a1 = nodes[nodes[i].connections[j].destIndex].p;
				Segment seg(a0, a1);
				nghb_graph_line.push_back(seg);
			}
		}

		Line full_line = nghb_graph_line;

		int number_segs = full_line.get_nb_segs();
		std::vector<std::string> property_names(1, "cost");
		std::vector<std::vector<std::vector<double>>> properties;

		properties.resize(number_segs);
		for (int segi = 0; segi < int(properties.size()); segi++) {
			properties[segi].resize(1); // there is 1 property only
		}

		for (int n = 0; n < number_segs; n++) {
			properties[n][0].push_back(list_cost[2*n]);
			properties[n][0].push_back(list_cost[2*n+1]);
		}
		std::string full_name = network_name + "_karst.txt";
		std::string full_dir_name = geologicalparams.directoryname + "/outputs";

		save_line(full_name, full_dir_name, full_line, property_names, properties); // version with costs

	}
}