/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for new methods + modifications to original methods
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for original versions of methods Compute_Edge_Cost, NodeIndex, BuildNearestNeighborGraph, ComputeKarsticSkeleton
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.

***************************************************************/

#include "KarstNSim/graph.h"

namespace KarstNSim {

	/*
	\brief Default Constructor.
	*/
	GraphOperations::GraphOperations() : CostGraph(0)
	{
	}


	// Samples the domain with Dwork's algorithm adapted in 3D
	void GraphOperations::SampleSpaceDwork( Box* box,
		const bool use_density_property, const int k, std::vector<double> propdensity, Surface* topo_surface) {

		//const clock_t time1 = clock();
		PointCloud centers2D = topo_surface->get_centers_cloud(2);

		double min_z_topo = topo_surface->get_boundbox_min().z;

		// 1 //
		//First, we get the Box characteristics (origin and dimensions in each direction as well as density property).
		//The density property must be defined on a grid with int(sqrt(3) / r_min) cells in each direction, with r_min the
		//minimum density defined in the whole grid (0<rmin<1, e.g., 0.01).

		const Vector3& ptmin = box->get_basis();
		const Vector3& ptmax = box->get_end();
		
		const Vector3 u = box->get_u();
		const Vector3 v = box->get_v();
		const Vector3 w = box->get_w();

		double delta_y = ptmax.y - ptmin.y;
		double delta_x = ptmax.x - ptmin.x;
		double delta_z = ptmax.z - ptmin.z;

		double min_xyz = std::min({ std::abs(delta_x), std::abs(delta_y), std::abs(delta_z) });
		double dim1, dim2, dim3;
		if (min_xyz == std::abs(delta_x)) {
			dim1 = 1;
			dim2 = abs(delta_y / delta_x);
			dim3 = abs(delta_z / delta_x);
		}
		else if (min_xyz == std::abs(delta_y)) {
			dim1 = abs(delta_x / delta_y);
			dim2 = 1;
			dim3 = abs(delta_z / delta_y);
		}
		else {
			dim1 = abs(delta_x / delta_z);
			dim2 = abs(delta_y / delta_z);
			dim3 = 1;
		}
	
		Vector3 min(ptmin.x, ptmin.y, ptmin.z);
		double r_min = 0;
		int nu = box->get_nu(), nv = box->get_nv(), nw = box->get_nw(); //number of samples in u,v,w axis


		double r_cst = 0;
		// We open the property if that's what the user asked for
		if (use_density_property) {

			// find min radius value in array
			r_min = 1;
			int minIndex = 0;
			for (int i = 0; i < nu*nv*nw; i++) {
				if ((propdensity[i] < r_min) && (propdensity[i] > 0)) {
					r_min = propdensity[i];
					minIndex = i;
				}
			}
		}
		else { // else we keep a constant density everywhere
			r_cst = params.graphPoissonRadius;
			r_min = r_cst;
		}
		int dimx = (int)(dim1*std::sqrt(3) / r_min); // dimensions in the Bridson/Dwork algorithms must always be equal to sqrt(3)/r_min (rounded down)
		int dimy = (int)(dim2*std::sqrt(3) / r_min);
		int dimz = (int)(dim3*std::sqrt(3) / r_min);

		//2//
		// Now we intialize the active and point list, as well as the background grid that will be used later
		// It must be noted that the algorithm is based on a cubic grid of size 1 and with int(sqrt(3) / r_min) cells in each direction.
		// Whatever the original dimensions of the Box are, they are transformed during the algorithm and back-transformed at the end.
		std::vector<int> active_list({ 1 });
		bool init_pt_in_grid = false;
		Vector3  pos_init;
		int64_t flattened_index = 0;
		int64_t converted_flat_index = 0;
		if (use_density_property) {
			while (!init_pt_in_grid) {
				pos_init.x = generateRandomDouble(0, dim1);
				pos_init.y = generateRandomDouble(0, dim2);
				pos_init.z = generateRandomDouble(0, dim3); // first random point of index 1
				flattened_index = flatten_indices2(dimx, dimy, dimz, dim1, dim2, dim3, pos_init); // we get the flattened index of the new point
				converted_flat_index = convert_index(dimx, dimy, dimz, nu, nv, nw, flattened_index);
				Vector3 real_pt;
				real_pt.x = min.x + (pos_init.x * u.x / dim1 + pos_init.y * v.x / dim2 + pos_init.z * w.x / dim3);
				real_pt.y = min.y + (pos_init.x * u.y / dim1 + pos_init.y * v.y / dim2 + pos_init.z * w.y / dim3);
				real_pt.z = min.z + (pos_init.x * u.z / dim1 + pos_init.y * v.z / dim2 + pos_init.z * w.z / dim3);
				if (propdensity[converted_flat_index] > 0 && CheckBelowSurf(real_pt,topo_surface,centers2D)) { // we make sure here that the 1st random point is indeed a) in the karstifiable zone (propdensity >0) and b) below topography
					init_pt_in_grid = true;
				}
			}
		}
		else {
			while (!init_pt_in_grid) {
				pos_init.x = generateRandomDouble(0, dim1);
				pos_init.y = generateRandomDouble(0, dim2);
				pos_init.z = generateRandomDouble(0, dim3);
				Vector3 real_pt;
				real_pt.x = min.x + (pos_init.x * u.x / dim1 + pos_init.y * v.x / dim2 + pos_init.z * w.x / dim3);
				real_pt.y = min.y + (pos_init.x * u.y / dim1 + pos_init.y * v.y / dim2 + pos_init.z * w.y / dim3);
				real_pt.z = min.z + (pos_init.x * u.z / dim1 + pos_init.y * v.z / dim2 + pos_init.z * w.z / dim3);
				if (CheckBelowSurf(real_pt, topo_surface,centers2D)) { // we make sure here that the 1st random point is indeed below topography
					init_pt_in_grid = true;
				}
			}
		}

		Vector3 new_p; // add the first point to the samples list
		new_p.x = min.x + (pos_init.x * u.x / dim1 + pos_init.y * v.x / dim2 + pos_init.z * w.x / dim3);
		new_p.y = min.y + (pos_init.x * u.y / dim1 + pos_init.y * v.y / dim2 + pos_init.z * w.y / dim3);
		new_p.z = min.z + (pos_init.x * u.z / dim1 + pos_init.y * v.z / dim2 + pos_init.z * w.z / dim3);
		samples.push_back(new_p);
		std::vector<Vector3> point_list;
		point_list.push_back(pos_init);

		int64_t gridSize = static_cast<int64_t>(dimx) * dimy * dimz;

		// VER 1 : unordered map
		//std::unordered_map<int, std::vector<int>> grid;
		// VER 2 : unordered multimap
		//std::multimap<int, int> grid;
		// VER 3 :  map
		//std::map<int, std::vector<int>> grid;

		std::vector<std::vector<int>> grid(gridSize);

		double r_local;
		if (use_density_property) {
			converted_flat_index = convert_index(dimx, dimy, dimz, nu, nv, nw, flattened_index);  // we convert this index to the property grid's own dimensions (probably different from the one we're using)
			r_local = propdensity[converted_flat_index];
		}
		else {
			r_local = r_cst;
		}
		std::vector<double> r_list({ r_local });
		std::vector<int64_t> scan_area_local = scan_area2(r_local, r_min, pos_init, dimx, dimy, dimz, dim1, dim2, dim3); // we get the scan area corresponding to the local density r_local
		for (int64_t x : scan_area_local) {
			grid[x].push_back(1);
			//grid.insert({ x, 1 });
		}
		//const clock_t time2 = clock();
		//double time_2 = 0.;
		//double time_3 = 0.;
		//double time_4 = 0.;
		//double time_5 = 0.;
		//double time_6 = 0.;
		while (!active_list.empty()) { // while we can still add points
			// find a random index in the active list

			std::size_t random_idx_number = generateRandomIndex(active_list.size());
			int random_idx = active_list[random_idx_number];
			double r_active;
			if (use_density_property) {
				r_active = r_list[random_idx - 1];
			}
			else {
				r_active = r_cst;
			}
			Vector3 random_xi = point_list[random_idx - 1];
			int crit_active = 0;
			int iter_k;
			// 3 //
			// We will now create k points around the active point and trim them if needed with the grid local scan

			for (iter_k = 0; iter_k < k; iter_k++) {
				//const clock_t time3 = clock();
				Vector3 vec_dir(generateNormalRandom(0,1), generateNormalRandom(0, 1), generateNormalRandom(0, 1));
				vec_dir = Normalize(vec_dir);// normalize the vector
				double magnitude_circle = std::cbrt(generateRandomDouble(8 * r_active*r_active*r_active, 27 * r_active*r_active*r_active)); // the magnitude corresponds to the cube root of a value distributed uniformly between the two radii
				Vector3 k_pt(vec_dir[0] * magnitude_circle + random_xi.x, vec_dir[1] * magnitude_circle + random_xi.y, vec_dir[2] * magnitude_circle + random_xi.z); // construct the new k point
				flattened_index = 0; // we get the flattened index of the new point
				converted_flat_index = 0;
				if (use_density_property) {
					flattened_index = flatten_indices2(dimx, dimy, dimz, dim1, dim2, dim3, k_pt);
					converted_flat_index = convert_index(dimx, dimy, dimz, nu, nv, nw, flattened_index);

					if ((k_pt.x > dim1) || (k_pt.x < 0) || (k_pt.y > dim2) || (k_pt.y < 0) || (k_pt.z > dim3) || (k_pt.z < 0) || (propdensity[converted_flat_index] < 0)) {
						crit_active++;
						continue; // if we're outside, we can skip this k_pt. We also add 1 to crit_active to let know that this point wasn't generated
					}
				}

				else {
					flattened_index = flatten_indices2(dimx, dimy, dimz, dim1, dim2, dim3, k_pt);
					if ((k_pt.x > dim1) || (k_pt.x < 0) || (k_pt.y > dim2) || (k_pt.y < 0) || (k_pt.z > dim3) || (k_pt.z < 0)) {
						crit_active++;
						continue; // if we're outside, we can skip this k_pt. We also add 1 to crit_active to let know that this point wasn't generated
					}
				}
				double r_k;
				if (use_density_property) {
					r_k = propdensity[converted_flat_index]; // we get local density at position k
				}
				else {
					r_k = r_cst;
				}
				bool crit_dist = true; // this criterion lets us keep track of the distance of the new point to previously generated points
				//const clock_t time4 = clock();

				//auto range = grid.equal_range(flattened_index);
				//for (auto it = range.first; it != range.second; ++it) {
				//	int point_index = it->second;
				//	Vector3 point_to_check = point_list[point_index - 1];
				//	double dist = magnitude(point_to_check - k_pt);
				//	if (dist < 2 * r_k) { // If distance is too small
				//		crit_dist = false; // Don't add the point
				//		crit_active++;
				//		break; // And don't bother checking the other points
				//	}
				//}
				if ((!grid[flattened_index].empty())) {
				//if (grid.count(flattened_index)) { // first check if there's a key associated in the grid (if it's not the case, why even bother checking anything)
					std::vector<int> tested_pts = grid[flattened_index];
					std::vector<int>::iterator iter_pts = tested_pts.begin();
					for (iter_pts; iter_pts < tested_pts.end(); iter_pts++) { // iterate on the points that should be checked
						Vector3 point_to_check = point_list[*iter_pts - 1];
						double dist = magnitude(point_to_check-k_pt);
						if (dist < 2 * r_k) { // If distance is too small
							crit_dist = false; // Don't add the point
							crit_active++;
							break; // And don't bother checking the other points
						}
					}
				}

				if (crit_dist) {
					// If we get there, it means that the point is far enough from any other point
					Vector3 new_pt; // add the new point to the sampled list
					new_pt.x = min.x + (k_pt.x * u.x / dim1 + k_pt.y * v.x / dim2 + k_pt.z * w.x / dim3);
					new_pt.y = min.y + (k_pt.x * u.y / dim1 + k_pt.y * v.y / dim2 + k_pt.z * w.y / dim3);
					new_pt.z = min.z + (k_pt.x * u.z / dim1 + k_pt.y * v.z / dim2 + k_pt.z * w.z / dim3);
					int valid = 1;
					if (new_pt.z > min_z_topo) {
						if (!(CheckBelowSurf(new_pt, topo_surface, centers2D))) {
							crit_active++;
							valid = 0;
						}
					}
					if (valid == 1)
					{
						//const clock_t time5 = clock();

						point_list.push_back(k_pt); // add the point to list
						samples.push_back(new_pt);
						if (use_density_property) {
							r_list.push_back(r_k); // add its density to the density list
						}
						//const clock_t time6 = clock();

						active_list.push_back(int(point_list.size())); // add it as an active point
						std::vector<int64_t> scan_area_k = scan_area2(r_k, r_min, k_pt, dimx, dimy, dimz, dim1, dim2, dim3); // we get the scan area corresponding to the local density r_k
						//const clock_t time7 = clock();

						for (int64_t x : scan_area_k) {
							grid[x].push_back(int(point_list.size()));
							//grid.insert({ x,int(point_list.size())});
						}
						//for (int x : scan_area_k) { // update the grid
						//	grid[x].push_back(int(point_list.size())); // since we just added the new point, its index corresponds to the size of the point list
						//}
						//const clock_t time8 = clock();
						//time_4 += double(time6 - time5);
						//time_5 += double(time7 - time6);
						//time_6 += double(time8 - time7);
					}


				}
				//time_2 += double(time4-time3);

			}
			if (crit_active == k) {  // If no k_pt was successfully added, we remove the root point from the active list
				active_list.erase(active_list.begin() + random_idx_number);
			}
		}

		propdensity.clear();
		grid.clear(); 	// Get rid of the map memory
	}

	bool GraphOperations::CheckBelowSurf(const Vector3& pt, const Surface* surface, const PointCloud& centers2D)
	{
		const int number_of_tries = 5;
		bool found_triangle = false;
		Vector3 P1, P2, P3;
		const Vector3 pt2D(pt.x, pt.y, 0.);
		std::vector<Neighbour> candidates;
		centers2D.findNearestNeighbors(pt2D, -1, number_of_tries, 1e50, candidates);
		for (int test = 0; test < number_of_tries; test++) {
			const Triangle& trgl_i = surface->get_triangle(candidates[test].i);
			const int& pt1 = trgl_i.point(0), pt2 = trgl_i.point(1), pt3 = trgl_i.point(2);
			const Vector3& p1 = surface->get_node(pt1), p2 = surface->get_node(pt2), p3 = surface->get_node(pt3);
			if (isInsideTriangle(p1, p2, p3, pt)) {
				found_triangle = true;
				P1 = p1;
				P2 = p2;
				P3 = p3;
				break;
			}
		}
		if (!found_triangle) return false; // if point is too far from all the triangles in map view, it means that it's outside the bounds of the surface, and therefore in an unsaturated zone
		return (is_point_under_triangle(pt, P1, P2, P3));
	}

	double GraphOperations::distsurf(const Vector3& pt, Surface* surface, const double& max_inception_surface_distance, const PointCloud& centers2D)
	{
		const int number_of_tries = 5;
		bool found_triangle = false;
		Vector3 P1, P2, P3;
		const Vector3 pt2D(pt.x, pt.y, 0.);
		std::vector<Neighbour> candidates;
		centers2D.findNearestNeighbors(pt2D, -1, number_of_tries, 1e50, candidates);

		for (int test = 0; test < number_of_tries; test++) {
			const Triangle& trgl_i = surface->get_triangle(candidates[test].i);
			const int& pt1 = trgl_i.point(0), pt2 = trgl_i.point(1), pt3 = trgl_i.point(2);
			const Vector3& p1 = surface->get_node(pt1), p2 = surface->get_node(pt2), p3 = surface->get_node(pt3);

			if (isInsideTriangle(p1, p2, p3, pt)) {
				found_triangle = true;
				P1 = p1;
				P2 = p2;
				P3 = p3;
				break;
			}
		}

		if (!found_triangle) return max_inception_surface_distance; // if point is too far from all the triangles in map view, it means that it's outside the bounds of the surface, and therefore in an unsaturated zone

		// Compute barycentric coordinates for interpolation
		double distance = barycentric(pt, P1, P2, P3);

		if (distance < 1e-2) {
			distance = 0.;
		}

		return distance;
	}

	void GraphOperations::DefineWSurfaceFlags(std::vector<Surface>* water_tables) {

		// First and foremost, a water table doesn't project on the whole domain all the time.
		// Therefore any point outside of the projected perimeter of the water table can be considered in the vadose zone.
		std::vector<PointCloud> all_centers2D;

		std::vector<std::vector<double>> coord_boundbox_wt(params.nb_springs); // elements are sorted as followed : xmin,ymin,xmax,ymax
		for (int i = 0; i < water_tables->size(); i++) {
			Vector3 min_wt = water_tables->at(i).get_boundbox_min();
			coord_boundbox_wt[i].push_back(min_wt.x);
			coord_boundbox_wt[i].push_back(min_wt.y);
			Vector3 max_wt = water_tables->at(i).get_boundbox_max();
			coord_boundbox_wt[i].push_back(max_wt.x);
			coord_boundbox_wt[i].push_back(max_wt.y);
			PointCloud centers2D = water_tables->at(i).get_centers_cloud(2);
			all_centers2D.push_back(centers2D);
		}
		samples_surf_flags.resize(int(samples.size()),std::vector<bool>(params.nb_springs));
		for (int i = 0; i<int(samples.size()); i++) {

			for (int j = 0; j < water_tables->size(); j++) {
				bool check;
				if ((samples[i].x < coord_boundbox_wt[j][0]) || (samples[i].x > coord_boundbox_wt[j][2]) || (samples[i].y < coord_boundbox_wt[j][1]) || (samples[i].y > coord_boundbox_wt[j][3])) { // out of wt bonds
					check = false;
				}
				else {
					Surface* surf_j = &water_tables->at(j);

					check = CheckBelowSurf(samples[i], surf_j,all_centers2D[j]);
				}
				samples_surf_flags[i][params.wt_surf_index[j]-1] = check;
			}
		}
	}

	/**
	 * .
	 * 
	 * \param surfaces
	 * \param max_inception_surface_distance
	 * \return 
	 */
	void GraphOperations::compute_euclidian_distance_from_wt(std::vector<Surface> *surfaces, const double max_inception_surface_distance) {

		std::vector<PointCloud> all_centers2D;
		for (int i = 0; i < surfaces->size(); i++) {
			PointCloud centers2D = surfaces->at(i).get_centers_cloud(2);
			all_centers2D.push_back(centers2D);
		}

		samples_wt_dist.resize(int(samples.size()), std::vector<double>(params.nb_springs));
		std::vector<double> max_dist_wt(params.nb_springs,0.);
		std::vector<double> min_dist_wt(params.nb_springs, std::numeric_limits<double>::infinity());
		for (int j = 0; j < surfaces->size(); j++) {
			Surface* surf_j = &surfaces->at(j);
			for (int i = 0; i<int(samples.size()); i++) {
				double distance = distsurf(samples[i], surf_j, max_inception_surface_distance, all_centers2D[j]);
				samples_wt_dist[i][params.wt_surf_index[j] - 1] = distance;
			}
		}
	}

	/**
	 * .
	 * 
	 * \param surfaces
	 * \param max_inception_surface_distance
	 * \return 
	 */
	void GraphOperations::compute_euclidian_distance_from_surfaces(std::vector<Surface> *surfaces, const double max_inception_surface_distance) {

		std::vector<PointCloud> all_centers2D;
		for (int i = 0; i < surfaces->size(); i++) {
			PointCloud centers2D = surfaces->at(i).get_centers_cloud(2);
			all_centers2D.push_back(centers2D);
		}

		samples_surf_dist.resize(int(samples.size()));
		for (int i = 0; i<int(samples.size()); i++) {
			double max_dist_surf_iter = std::numeric_limits<double>::infinity();
			for (int j = 0; j < surfaces->size(); j++) {
				Surface* surf_j = &surfaces->at(j);
				double distance = distsurf(samples[i], surf_j, max_inception_surface_distance,all_centers2D[j]);
				max_dist_surf_iter = std::min(max_dist_surf_iter, distance);
			}
			samples_surf_dist[i] = max_dist_surf_iter;
		}
	}

	void GraphOperations::process_node_in_samples(const Vector3& node, std::vector<int>& flags) {
		int test_already_defined = -1;
		for (int k = 0; k < samples.size(); k++) {
			if (samples[k] == node) {
				test_already_defined = k + 1; // if the corresponding point is already in the samples, we keep its current index
				break;
			}
		}
		if (test_already_defined == -1) { // if the point hasn't been added yet
			samples.push_back(node); // we add it
			// Check if flags is not empty
			//if (!flags.empty()) {
			flags.push_back( int(samples.size()) - 1 );
			//}
		}
		else { // if it has already been added, add the corresponding index
			// Check if flags is not empty
			//if (!flags.empty()) {
			flags.push_back( test_already_defined - 1 );
			//}
		}
		samples.push_back(node);
	}

	void GraphOperations::save_nghb_graph(const GeologicalParameters& geologicalparams, const bool& create_nghb_graph_property) {

		std::string full_dir_name = geologicalparams.directoryname + "/outputs";
		std::string full_name = params.scenename + "_nghb_graph.txt";
		std::vector<Segment> nghb_graph_line;
		std::vector<std::vector<double>> list_cost_up;
		std::vector<std::vector<double>> list_cost_down;
		list_cost_up.resize(params.nb_springs);
		list_cost_down.resize(params.nb_springs);
		for (int cost_i = 0; cost_i < params.nb_springs; cost_i++) {
			for (int i = 0; i < adj.size(); i++) {
				for (int j = 0; j < adj[i].size(); j++) {
					if (create_nghb_graph_property) { // important note : a line property can only be defined on its nodes on Skua Gocad, therefore we chose to consider that the cost at a given node is equal at the same time to all edge costs linked to it
						if (samples[i].z > samples[adj[i][j].target].z) {
							list_cost_down[cost_i].push_back(adj[i][j].weight[cost_i]);
							list_cost_down[cost_i].push_back(adj[i][j].weight[cost_i]);
							bool found_it = false;
							for (int neighbor = 0; neighbor < adj[adj[i][j].target].size(); neighbor++) {
								if (adj[adj[i][j].target][neighbor].target == i) {
									found_it = true;
									list_cost_up[cost_i].push_back(adj[adj[i][j].target][neighbor].weight[cost_i]);
									list_cost_up[cost_i].push_back(adj[adj[i][j].target][neighbor].weight[cost_i]);
									break;
								}
							}
							if (!found_it) { // sometimes, due to the rules of the nghb graph, an edge only exists in one direction. In that case we define the other direction as equal to the first one
								list_cost_up[cost_i].push_back(adj[i][j].weight[cost_i]);
								list_cost_up[cost_i].push_back(adj[i][j].weight[cost_i]);
							}
						}
						else {
							list_cost_up[cost_i].push_back(adj[i][j].weight[cost_i]);
							list_cost_up[cost_i].push_back(adj[i][j].weight[cost_i]);
							bool found_it = false;
							for (int neighbor = 0; neighbor < adj[adj[i][j].target].size(); neighbor++) {
								if (adj[adj[i][j].target][neighbor].target == i) {
									found_it = true;
									list_cost_down[cost_i].push_back(adj[adj[i][j].target][neighbor].weight[cost_i]);
									list_cost_down[cost_i].push_back(adj[adj[i][j].target][neighbor].weight[cost_i]);
									break;
								}
							}
							if (!found_it) { // same
								list_cost_down[cost_i].push_back(adj[i][j].weight[cost_i]);
								list_cost_down[cost_i].push_back(adj[i][j].weight[cost_i]);
							}
						}
					}

					if (cost_i == 0) { // This should only be done once in the property loop

						Vector3 a0 = samples[i];
						Vector3 a1 = samples[adj[i][j].target];
						Segment seg(a0, a1);
						nghb_graph_line.push_back(seg);
					}
				}
			}
		}

		Line full_line(nghb_graph_line);

		if (create_nghb_graph_property) {
			int number_segs = full_line.get_nb_segs();
			std::vector<std::string> property_names;
			std::vector<std::vector<std::vector<double>>> properties;
			properties.resize(number_segs);
			for (int segi = 0; segi < int(properties.size()); segi++) {
				properties[segi].resize(params.nb_springs * 2); // there are 2*number of springs properties (one for each cost)
			}
			for (int cost_i = 0; cost_i < params.nb_springs; cost_i++) {

				std::string name_line_up = "cost_down" + std::to_string(cost_i + 1);
				std::string name_line_down = "cost_up" + std::to_string(cost_i + 1);
				property_names.push_back(name_line_up);
				property_names.push_back(name_line_down);

				for (int n = 0; n < number_segs; n++) {
					properties[n][2 * cost_i].push_back(list_cost_up[cost_i][2 * n + 1]);
					properties[n][2 * cost_i].push_back(list_cost_up[cost_i][2 * n + 1]);
				}

				for (int n = 0; n < number_segs; n++) {
					properties[n][2 * cost_i + 1].push_back(list_cost_down[cost_i][2 * n + 1]);
					properties[n][2 * cost_i + 1].push_back(list_cost_down[cost_i][2 * n + 1]);
				}

			}
			save_line(full_name, full_dir_name, full_line, property_names, properties); // version with costs
		}
		else {
			save_line(full_name, full_dir_name, full_line); // version without costs
		}
	}

	void GraphOperations::InitializeCostGraph(const bool create_nghb_graph, const bool create_nghb_graph_property, const std::vector<KeyPoint>& keypts,
		const GeologicalParameters& geologicalparams, const std::vector<Vector3>& nodes_on_inception_surfaces, const std::vector<std::vector<Vector3>>& nodes_on_wt_surfaces,
		std::vector<Surface>* inception_horizons, std::vector<Surface>* water_tables, const bool use_sampling_points,
		Box* box, const double max_inception_surface_distance,
		std::vector<Vector3>* sampling_points, const bool create_vset_sampling, const bool use_density_property, const int k_pts,
		const double fraction_old_karst_perm, std::vector<double> propdensity, std::vector<double> propikp, Surface* topo_surface)
	{
		PointCloud centers2D = topo_surface->get_centers_cloud(2);

		const clock_t time1 = clock();

		params = geologicalparams;
		int nu = box->get_nu(), nv = box->get_nv();

		// Adaptive sampling of space
		std::vector<std::vector<int>> on_wt_flags_idx(params.nb_springs);
		std::vector<std::vector<int>> on_surf_flags_idx(1);
		Vector3 new_pt;

		// Key points are also added to sample set
		for (int i = 0; i < keypts.size(); i++) {
			samples.push_back(keypts[i].p);
			if (keypts[i].type == KeyPointType::Spring) {
				on_wt_flags_idx.at(keypts[i].wt_idx).push_back(int(samples.size()) - 1);
			}
		}

		if (use_sampling_points) { // first we recover previously sampled points if any
			for (int i = 0; i < sampling_points->size(); i++) {
				new_pt.x = sampling_points->at(i).x;
				new_pt.y = sampling_points->at(i).y;
				new_pt.z = sampling_points->at(i).z;
				samples.push_back(new_pt);
			}
		}
		else { // else we sample them with modified Dwork poisson sphere sampling
			SampleSpaceDwork(box, use_density_property, k_pts, propdensity,topo_surface);
		}

		for (int i = 0; i < nodes_on_inception_surfaces.size(); i++) {
			double value_itr = 1.;
			if (!use_sampling_points && use_density_property) {
				int u1, v1, w1;
				box->xyz2uvw_with_limits_conditions(nodes_on_inception_surfaces[i], u1, v1, w1);
				value_itr = propdensity[u1 + nu * v1 + nv * nu*w1];
			}
			bool below_top_test = true;
			if (!topo_surface->is_empty()) {
				below_top_test = CheckBelowSurf(nodes_on_inception_surfaces[i], topo_surface, centers2D);
			}
			if (below_top_test && value_itr > 0) {  // Accept only nodes that are below topo surf AND are associated to a karstifiable zone (ie. propdensity >0)
				process_node_in_samples(nodes_on_inception_surfaces[i], on_surf_flags_idx[0]);
			}
		}

		//Add water table samples
		for (int i = 0; i < nodes_on_wt_surfaces.size(); i++) { // iterate on nb of wt surfaces
			for (int j = 0; j < nodes_on_wt_surfaces[i].size(); j++) { // iterate on points of wt i
				// Accept only nodes that are below topo surf AND are associated to a karstifiable zone (ie. propdensity >0)
				double value_itr = 1.;
				if (!use_sampling_points && use_density_property) {
					int u1, v1, w1;
					box->xyz2uvw_with_limits_conditions(nodes_on_wt_surfaces[i][j], u1, v1, w1);
					value_itr = propdensity[u1 + nu * v1 + nv * nu*w1];
				}
				bool below_top_test = true;
				if (!topo_surface->is_empty()) {
					below_top_test = CheckBelowSurf(nodes_on_wt_surfaces[i][j], topo_surface, centers2D);
				}
				if (below_top_test && value_itr > 0) {  // Accept only nodes that are below topo surf AND are associated to a karstifiable zone (ie. propdensity >0)
					process_node_in_samples(nodes_on_wt_surfaces[i][j],on_wt_flags_idx[i]);
				}
			}
		}

		// Add previously simulated karst networks samples
		for (int i = 0; i < params.PtsOldGraph.size(); i++) {
			for (int j = 0; j < params.PtsOldGraph[i].size(); j++) { // =2 (start and end points of each segment of the old graph)
				process_node_in_samples(params.PtsOldGraph[i][j], params.IdxOldGraph[i]);
			}
		}

		if (create_vset_sampling) {
			std::string full_name = geologicalparams.scenename + "_pts.txt";
			std::string full_dir_name = geologicalparams.directoryname + "/outputs";
			save_pointset(full_name, full_dir_name, samples);
		}

		const clock_t time2 = clock();

		if (params.karstificationCost.used) {
			const Vector3 u = box->get_u();
			const Vector3 v = box->get_v();
			const Vector3 w = box->get_w();
			int nu = box->get_nu(), nv = box->get_nv(); 
			for (int i = 0; i < samples.size(); i++) {
				int u1, v1, w1;
				box->xyz2uvw_with_limits_conditions(samples[i], u1, v1, w1);
				double value_itr = propikp[u1 + nu * v1 + nv * nu*w1];
				if (value_itr < 0) { // a value smaller than  means a No data value, which means a cell above topography
					samples_layer_kp.push_back(1); // maximal cost (this happens because of the discrepancy between the background grid's resolution and the points' position close to the topography)
					continue;
				}
				samples_layer_kp.push_back(1- propikp[u1 + nu * v1 + nv * nu*w1]);
			}
			propikp.clear();
		}
		const clock_t time21 = clock();

		// For each sampling point, this defines whether the point is below or above each water table surface, and also if they are exactly onto the wt surface

		DefineWSurfaceFlags(water_tables); // above or under wt test
		samples_on_wt_flags = get_flags_from_indices(samples, on_wt_flags_idx); // exactly on wt test

		// For each sampling point, this computes the distance of each point to each water table
		const clock_t time22 = clock();


		compute_euclidian_distance_from_wt(water_tables, max_inception_surface_distance);
		//std::vector<double> max_dist_wt = pair.first;
		//std::vector<double> min_dist_wt = pair.second;

		if (inception_horizons != nullptr) {
			compute_euclidian_distance_from_surfaces(inception_horizons, max_inception_surface_distance);
		}
		const clock_t time23 = clock();


		const clock_t time3 = clock();

		// Nearest neighbour graph initialization
		BuildNearestNeighbourGraph(box, fraction_old_karst_perm, max_inception_surface_distance);

		const clock_t time4 = clock();

		//create_nghb_graph
		if (create_nghb_graph) {
			save_nghb_graph( geologicalparams, create_nghb_graph_property);
		}
	}

	/*!
	\brief Builds the nearest neighbour graph from the set of samples.
	*/
	void GraphOperations::BuildNearestNeighbourGraph(Box* box, const double fraction_old_karst_perm, double max_dist_surf)
	{

		double time_neigbhors = 0.;
		double time_costs = 0.;

		const Vector3& ptmin = box->get_basis();
		const Vector3& ptmax = box->get_end();
		const Vector3 u = box->get_u();
		const Vector3 v = box->get_v();
		const Vector3 w = box->get_w();
		double delta_y_abs = std::abs(ptmax.y - ptmin.y);
		double delta_x_abs = std::abs(ptmax.x - ptmin.x);
		double delta_z_abs = std::abs(ptmax.z - ptmin.z);
		double dim_max = std::max({ delta_y_abs, delta_x_abs, delta_z_abs });
		double one_over_delta_y = 1 / (delta_y_abs);
		double one_over_delta_x = 1 / (delta_x_abs);
		double one_over_delta_z = 1 / (delta_z_abs);

		adj.resize(samples.size());
		double R = 0.; // std::numeric_limits<double>::infinity();
		double distmax = 1e50;
		if (params.graphuse_max_nghb_radius) {
			R = params.graphNeighbourRadius * params.graphNeighbourRadius; // Precise the value of R if it has been defined by the user
			distmax = (R / (dim_max*dim_max));
		}
		const int N = params.graphNeighbourCount;
		double dmin = std::numeric_limits<double>::infinity();
		double dmax = 0.;

		// We want to use a stretched version of the point cloud for perfect scaling :

		std::vector<Vector3> samples_stretched;
		for (int i = 0; i < samples.size(); i++) {
			// Stretch the points based on scaling factors
			Vector3 p_stretched{
				samples[i].x * one_over_delta_x,
				samples[i].y * one_over_delta_y,
				samples[i].z *one_over_delta_z
			};
			samples_stretched.push_back(p_stretched);
		}

		PointCloud pointCloud(samples_stretched);
		
		for (int i = 0; i < samples.size(); i++) {
			const clock_t time0 = clock();
			const Vector3 p = samples[i];
			const Vector3 p_stretched = samples_stretched[i];
			std::vector<Neighbour> candidates;

			// Use the k-d tree to find nearest neighbors
			
			pointCloud.findNearestNeighbors(p_stretched,i, N, distmax, candidates);

			const clock_t time1 = clock();
			int n = KarstNSim::Min((int)candidates.size(), N);

			for (int j = 0; j < n; j++) {
				const Vector3 pn = samples[candidates[j].i];
				const std::pair<std::vector<double>, std::vector<bool>> pair = ComputeEdgeCost(i, p, pn, dmin, dmax, max_dist_surf);
				SetEdge(i, candidates[j].i, pair.first, pair.second);
			}
			const clock_t time2 = clock();
			time_neigbhors += double(time1 - time0);
			time_costs += double(time2 - time1);
		}
		const clock_t time3 = clock();
		// add distance cost (has to be added at the end before it depends on the max distance)
		//double min_val = std::numeric_limits<double>::min();
		double min_val = 1e-5;
		double eps_zero = 1e-5;
		for (int i = 0; i < adj.size(); i++) {
			Vector3 p = samples[i];
			for (int j = 0; j<adj[i].size(); j++) {
				Vector3 pn = samples[adj[i][j].target];
				for (int cost_i = 0; cost_i < params.nb_springs; cost_i++) {

					// Very important step : multiplying the distance to the other costs opens the possibility that the full cost might be zero, which would make Dijkstra's algorithm not work.
					// Therefore, in that case, the cost should be thresholded to a minimum value. That way, the cost, albeit very small, will still be proportional to the distance.

					if (adj[i][j].weight[cost_i] < eps_zero) { // if value is zero (or almost zero)
						adj[i][j].weight[cost_i] = min_val; // change it to the "almost zero" default value.
					}
					adj[i][j].weight[cost_i] *= ComputeEdgeDistanceCost(p, pn, dmin, dmax);

				}
			}
		}

		// reduce cost of edges already karstified in previous iterations
		for (int i = 0; i < params.IdxOldGraph.size(); i++) {
			for (int j = 1; j < params.IdxOldGraph[i].size(); j++) {
				for (int k = 0; k < adj[params.IdxOldGraph[i][0]].size(); k++) {
					if (adj[params.IdxOldGraph[i][0]][k].target == params.IdxOldGraph[i][j]) { // if an edge was created between two nodes of the previous karst
						for (int cost_i = 0; cost_i < params.nb_springs; cost_i++) {
							adj[params.IdxOldGraph[i][0]][k].weight[cost_i] *= fraction_old_karst_perm;
						}
						break; // the edge was found, no need to scan further
					}
				}
			}
		}
		const clock_t time4 = clock();
		time_costs += double(time4 - time3);
		// At that point, the adjacency graph is entirely built, and ready to be used.
		// We keep a copy of adj to use later on for amplification, which works better without cost reduction

		adj_no_cost_red = adj; // copy


	}

	/*
	\brief Computes the normalized distance cost of an edge. This step is done separately because normalizing the distance implies to know
	all the distances in advance.
	\param p start position
	\param pn end position
	*/

	double GraphOperations::ComputeEdgeDistanceCost(const Vector3& p, const Vector3& pn,const double& dmin, const double& dmax)
	{
		double cost = 0;
		cost = ((magnitude(p - pn) - dmin) / (dmax - dmin)) * params.distanceCost.weight;
		cost = KarstNSim::Clamp(cost, 0.0, cost);
		return cost;
	}

	/*
	\brief Computes the cost of an edge : the cost is a combination of a fracture term, a inception horizon term, a water table term.
	\param p start position
	\param pn end position
	\param Surfaces Pointer list of all inception horizons
	*/

	std::pair<std::vector<double>, std::vector<bool>> GraphOperations::ComputeEdgeCost(const int& index, const Vector3& p, const Vector3& pn,
		double& dmin, double& dmax, const double& dist_max) const
	{
		//auto time1 = std::chrono::high_resolution_clock::now();

		std::vector<double> cost(params.nb_springs, params.multiply_costs ? 1.0 : 1.0);
		std::vector<bool> frac_flags(params.fractures_orientations.size(), false);

		const Vector3 diff = pn - p;
		const double dist = magnitude(diff);
		dmax = std::max(dmax, dist);
		dmin = std::min(dmin, dist);


		// Horizon cost
		if (params.horizonCost.used)
		{
			const double distance_to_surf = KarstNSim::Clamp(samples_surf_dist[index], 0.0, dist_max);
			double w = (1.0 - KarstNSim::CubicSmooth(distance_to_surf, dist_max));
			for (int cost_i = 0; cost_i < params.nb_springs; cost_i++) {
				if (params.multiply_costs) cost[cost_i] *= w * params.horizonCost.weight;
				else cost[cost_i] += w * params.horizonCost.weight;
			}
		}


		// Permeability
		if (params.karstificationCost.used)
		{
			for (int cost_i = 0; cost_i < params.nb_springs; cost_i++) {
				if (params.multiply_costs) cost[cost_i] *= samples_layer_kp[index] * params.karstificationCost.weight;
				else cost[cost_i] += samples_layer_kp[index] * params.karstificationCost.weight;
			}
		}


		// Fracture orientation
		if (params.fractureCost.used)
		{
			const double azimut = atan(diff[0]/diff[1]) * (180.0 / 3.141592653589793238463);
			int i;
			double costFrac = 1.;
			for (i = 0; i < params.fractures_orientations.size(); i++) {
				double orientation;
				if (params.fractures_orientations[i] >= 90) {
					orientation = params.fractures_orientations[i] - 180;
				}
				else {
					orientation = params.fractures_orientations[i];
				}
				double tolerance = params.fractures_tolerances[i];
				if (std::abs(azimut - orientation) < tolerance) {
					costFrac = 0.;
					frac_flags[i] = true;
					break;
				}
			}
			for (int cost_i = 0; cost_i < params.nb_springs; cost_i++) {
				if (params.multiply_costs) cost[cost_i] *= costFrac * params.fractureCost.weight;
				else cost[cost_i] += costFrac * params.fractureCost.weight;
			}
		}

		// Vadose and phreatic zone

		for (int cost_i = 0; cost_i < params.nb_springs; ++cost_i) {
			if (!samples_surf_flags[index][cost_i]) {
				const double cos_angle = Dot(Vector3(0, 0, 1), Normalize(diff)); // from -1 to 1
				const double e = (cos_angle + 1) / 2; // from 0 to 1
				if (params.multiply_costs) cost[cost_i] *= params.waterTable1.weight * e;
				else cost[cost_i] += params.waterTable1.weight * e;
			}
			else {
				double w_prim = KarstNSim::Clamp(samples_wt_dist[index][cost_i], 0.0, dist_max);
				double w = (1.0 - KarstNSim::CubicSmooth(w_prim, dist_max));
				if (params.multiply_costs) cost[cost_i] *= w * params.waterTable2.weight;
				else cost[cost_i] += w * params.waterTable2.weight;
			}
		}

		// Clamp the cost value to keep it positive for Dijkstra (this might not be necessary)
		for (int cost_i = 0; cost_i < params.nb_springs; cost_i++) {
			cost[cost_i] = KarstNSim::Clamp(cost[cost_i], 0.0, cost[cost_i]); // make sure to keep a positive cost for Dijsktra
		}

		return std::make_pair(cost, frac_flags);
	}

	/*
	\brief Computes the internal index of a sample of the graph, or -1 if the sample does not exist.
	\param p sample position
	*/
	int GraphOperations::NodeIndex(const Vector3& p) const
	{
		for (int i = 0; i < samples.size(); i++)
		{
			if (samples[i] == p)
				return i;
		}
		return -1;
	}

	KarsticSkeleton GraphOperations::ComputeKarsticSkeleton(const std::vector<KeyPoint>& pts, const double fraction_karst_perm)
	{
		// Internal representation of key points, storing their index (with two sublists : one for sinks and one for springs)
		std::vector<InternalKeyPoint> keypts;
		std::vector<InternalKeyPoint> keyptssinks;
		std::vector<InternalKeyPoint> keyptssprings;
		double eps = 1e5;
		for (int i = 0; i < int(pts.size()); i++) {
			keypts.push_back({ NodeIndex(pts[i].p), pts[i].p, pts[i].type });
			if (pts[i].type == KeyPointType::Sink)
			{
				keyptssinks.push_back({ NodeIndex(pts[i].p), pts[i].p, pts[i].type });
			}
			else if (pts[i].type == KeyPointType::Spring) {
				keyptssprings.push_back({ NodeIndex(pts[i].p), pts[i].p, pts[i].type });
			}
		}
		// Allocates all temporary arrays
		std::vector<std::vector<double>> all_distances;
		std::vector<std::vector<std::vector<double>>> all_distances_detailled; // version where each cost for each segment is precised
		std::vector<std::vector<std::vector<int>>> all_paths;
		std::vector<std::vector<std::vector<bool>>> all_vadose;
		all_paths.resize(keyptssinks.size());
		all_vadose.resize(keyptssinks.size());
		all_distances.resize(keyptssinks.size());
		all_distances_detailled.resize(keyptssinks.size());
		for (int i = 0; i < keyptssinks.size(); i++)
		{
			all_distances[i].resize(keyptssprings.size(), std::numeric_limits<double>::infinity());
			all_distances_detailled[i].resize(keyptssprings.size());
			all_paths[i].resize(keyptssprings.size());
			all_vadose[i].resize(keyptssprings.size());
		}
		// Iterate on all sinks points

		for (int i = 0; i < keyptssinks.size(); i++) {
			int source = keyptssinks[i].index;

			// Pick the inlet/outlet connections when val = 2 in connectivity matrix
			std::vector<int> list_of_two;
			for (int j = 0; j < keyptssprings.size(); j++) {
				int connectivity_flag = params.connectivity_matrix[params.sinks_index[i]-1][j];
				if (connectivity_flag == 2)
					list_of_two.push_back(j);
			}
			if (int(list_of_two.size()) != 0) {
				std::size_t random_idx_number = generateRandomIndex(list_of_two.size());
				int random_idx = list_of_two[random_idx_number];
				for (int j = 0; j < keyptssprings.size(); j++) {

					int connectivity_flag = params.connectivity_matrix[params.sinks_index[i]-1][j];
					
						if (connectivity_flag == 2 && j != random_idx)
							params.connectivity_matrix[params.sinks_index[i]-1][j] = 0;
						else if (connectivity_flag == 2 && j == random_idx)
							params.connectivity_matrix[params.sinks_index[i]-1][j] = 1;
				}
			}

			// now iterate on springs, and look for path only if connectivity flag is equal to 1

			for (int j = 0; j < keyptssprings.size(); j++) {

				int connectivity_flag = params.connectivity_matrix[params.sinks_index[i] - 1][j];
				if (connectivity_flag == 0) continue; // skip springs for which flag equals 0 (no connection)

				int target = keyptssprings[j].index;

				int reach = 0; // we try to find this index, corresponding to the first pt below the water table level that was reached

				// 1) GENERATING VADOSE ZONE CONDUIT

				std::vector<double> distances_vadose;
				std::vector<int> previous_vadose;
				bool already_reached = false;
				double pathSize_vadose = 0.;
				DijkstraComputePathsSurface(j, source, reach, distances_vadose, previous_vadose, samples_on_wt_flags,target, already_reached); // we compute the quickest path to the phreatic zone by finding the first node reached on the water table with Dijkstra
				//reach = FindNearestWaterTablePoint(j,source, samples_on_wt_flags, previous_vadose, distances_vadose);
				//if (target == reach) {
				//	already_reached = true;
				//}
				std::pair<std::vector<int>, std::vector<double>> pair_vadose = DijkstraGetShortestPathTo(reach, previous_vadose, distances_vadose, pathSize_vadose);
				std::vector<int> path_vadose = pair_vadose.first;
				std::vector<double> path_cost_vadose = pair_vadose.second;


				// 2) GENERATING PHREATIC ZONE CONDUIT

				std::vector<double> distances;
				std::vector<int> previous;
				double pathSize = 0.;
				std::vector<int> path;
				std::vector<double> path_cost;
				std::pair< std::vector<int>, std::vector<double>> pair;
				if (!(already_reached)) {
					DijkstraComputePaths(j, reach, distances, previous, target);
					pair = DijkstraGetShortestPathTo(target, previous, distances, pathSize);
					path = pair.first;
					path_cost = pair.second;
				}

				// Concatenation of vadose and phreatic conduit graphs
				std::vector<int> path_full;
				std::vector<double> path_cost_full;
				std::vector<bool> vadose_full;

				if (!path_vadose.empty()) {
					path_full.insert(path_full.end(), path_vadose.begin(), path_vadose.end());
					path_cost_full.insert(path_cost_full.end(), path_cost_vadose.begin(), path_cost_vadose.end());
					vadose_full.resize(path_vadose.size(), true); // Resize and initialize to true
				}

				double pathSize_full = pathSize + pathSize_vadose;

				if (!already_reached && !path.empty()) {
					vadose_full.insert(vadose_full.end(), path.size() - 1, false);
					path_full.insert(path_full.end(), path.begin() + 1, path.end());
					path_cost_full.insert(path_cost_full.end(), path_cost.begin() + 1, path_cost.end());
				}

				all_paths.at(i).at(j) = path_full;
				all_vadose.at(i).at(j) = vadose_full;
				all_distances.at(i).at(j) = pathSize_full;
				all_distances_detailled.at(i).at(j) = path_cost_full;
			}

			// Try to find closest spring from the doline, all other paths are removed
			double shortest_distance = std::numeric_limits<double>::infinity();
			for (int j = 0; j < all_paths[i].size(); j++) // iterate on springs
			{
				if (all_paths[i][j].size() <= 1) continue;
				shortest_distance = std::min(shortest_distance, all_distances[i][j]);
			}
			for (int j = 0; j < all_paths[i].size(); j++) // iterate on springs
			{
				if (all_paths[i][j].size() <= 1) continue;
				if (all_distances[i][j] > (shortest_distance+ eps) && params.allow_single_outlet)
				{
					all_paths[i][j].clear();
				}
				if ((all_distances[i][j] - shortest_distance)<eps)
				{
					// UPDATE THE PATH COST FOR NEXT ITERATIONS : segments already karstified, but only in phreatic zone

					for (int node_path = 0; node_path < all_paths[i][j].size() - 1; node_path++) // iterate on springs and waypoints
					{
						for (int cost_i = 0; cost_i < params.nb_springs; cost_i++)
						{
							if (params.vadose_cohesion || samples_surf_flags[all_paths[i][j][node_path]][cost_i]) // if we enabled cohesion everywhere (or, if not, if we're in the phreatic zone specifically)
							{ 								
								adj[all_paths[i][j][node_path]][GetIdxNeighbor(all_paths[i][j][node_path], all_paths[i][j][node_path + 1])].weight[cost_i] *= fraction_karst_perm;
							}
						}
					}
				}
			}
		}
		// Path pruning based on an anisotropic empty region criterion
		std::vector<std::vector<int>> pathsFinal;
		std::vector<std::vector<double>> costsFinal;
		std::vector<std::vector<bool>> vadoseFinal;
		for (int i = 0; i < all_paths.size(); i++)
		{
			for (int j = 0; j < all_paths[i].size(); j++)
			{
				if (all_paths[i][j].size() <= 1)
					continue;
				// Second, check if there is a cheaper path from nodes[i] going through a node[k] to arrive at nodes[j]
				double d_ij = KarstNSim::Pow(all_distances[i][j], params.gamma);
				bool keep = true;
				for (int k = 0; k < all_paths[i].size(); k++)
				{
					if (i == k || j == k)			 continue;
					if (all_paths[i][k].size() <= 1) continue;
					if (all_paths[k][j].size() <= 1) continue;

					double d_ik = KarstNSim::Pow(all_distances[i][k], params.gamma);
					double d_kj = KarstNSim::Pow(all_distances[k][j], params.gamma);
					if (d_ik + d_kj < d_ij)
					{
						keep = false;
						break;
					}
				}			
				if (keep) {
					pathsFinal.push_back(all_paths[i][j]);
					costsFinal.push_back(all_distances_detailled[i][j]);
					vadoseFinal.push_back(all_vadose[i][j]);
				}
			}
		}
		// Build karstic skeleton structure
		return KarsticSkeleton(this, pathsFinal, costsFinal, vadoseFinal);
	}

	std::pair< std::vector<int>, std::vector<double>> GraphOperations::shortest_path(bool under_wt, int u, int v)
	{
		int order;
		if (under_wt) {
			order = int(std::max_element(params.z_list.begin(), params.z_list.end()) - params.z_list.begin());
		}
		else {
			order = int(std::min_element(params.z_list.begin(), params.z_list.end()) - params.z_list.begin());
		}

		std::vector<double> distances;
		std::vector<int> previous;
		std::pair< std::vector<int>, std::vector<double>> pair;
		double pathSize = 0.;
		std::vector<int> path;
		std::vector<double> path_cost;
		DijkstraComputePaths(order, u, distances, previous); // Computes shortest paths from u to other vertices
		pair = DijkstraGetShortestPathTo(v, previous, distances, pathSize); // gets shortest path computed from v to u (back-propagated)
		path = pair.first;
		path_cost = pair.second;

		return std::make_pair(path,path_cost);
	}

	/*!
	\brief Amplify the graph from a new set of key points. The new points are linked to the existing base key points.
	\param baseSkeletonNodes the existing network
	\param newkeypts the new set of key point to connect to the existing network
	\return the set of added paths
	*/
	std::vector<std::vector<int>> GraphOperations::AmplifyKarsticSkeleton(const std::vector<KarsticNode>& baseSkeletonNodes,
		const std::vector<KeyPoint>& newkeypts)
	{
		// Our goal is to connect the new key points to the original network
		std::vector<std::vector<int>> pathsFinal;
		for (int i = 0; i < newkeypts.size(); i++)
		{
			std::vector<double> distances;
			std::vector<int> previous;

			double minDist = std::numeric_limits<double>::infinity();
			std::vector<int> bestPath;
			for (int j = 0; j < baseSkeletonNodes.size(); j++)
			{
				if (i == j)
					continue;
				int endIndex = baseSkeletonNodes[j].index;
				double totalDistance = 0.0;
				std::vector<int> path;
				std::vector<double> path_cost;
				std::pair< std::vector<int>, std::vector<double>> pair;
				pair = DijkstraGetShortestPathTo(endIndex, previous, distances, totalDistance);
				path = pair.first;
				if (path.size() <= 1)
					continue;
				if (totalDistance < minDist)
				{
					minDist = totalDistance;
					bestPath = path;
				}
			}
			pathsFinal.push_back(std::vector<int>(bestPath.begin(), bestPath.end()));
		}
		std::cout << "Additional path count: " << pathsFinal.size() << std::endl;
		return pathsFinal;
	}

	/*!
	\brief Add new samples to the nearest neighbour graph structure. 
	Existing samples are not removed, but their neighbours are modified.
	\param samples new samples passed as key points
	*/
	std::vector<GraphOperations::InternalKeyPoint> GraphOperations::AddNewSamples(const std::vector<KeyPoint>& newSamples)
	{
		struct Neighbour
		{
			int i;
			double d;
			bool operator<(const Neighbour& nei)
			{
				return d < nei.d;
			}
		};
		const double R = params.graphNeighbourRadius * params.graphNeighbourRadius;
		const int N = params.graphNeighbourCount;
		std::vector<InternalKeyPoint> ret;
		for (int i = 0; i < newSamples.size(); i++)
		{
			Vector3 p = newSamples[i].p;
			samples.push_back(p);
			adj.push_back({});
			int index = (int)samples.size() - 1;

			std::vector<Neighbour> candidates;
			for (int j = 0; j < samples.size(); j++)
			{
				if (index == j) continue;
				Vector3 pn = samples[j];
				double d = squaredmagnitude(p - pn);
				if (d < R)
					candidates.push_back({ j, d });
			}

			std::sort(candidates.begin(), candidates.end());
			int n = KarstNSim::Min((int)candidates.size(), N);
			for (int j = 0; j < n; j++)
			{
				Vector3 pn = samples[candidates[j].i];
				Vector3 d = Normalize(pn - p);
			}

			ret.push_back({ index, p, newSamples[i].type });
		}
		return ret;
	}


	/*!
	\brief Export the sample set as a point set file.
	\param path file path
	*/
	void GraphOperations::saveSamples(const std::string& path) const
	{
		std::ofstream out;
		out.open(path);

		// ply header
		out << "ply\n";
		out << "format ascii 1.0\n";
		out << "element vertex " << samples.size() << "\n";
		out << "property double x\n";
		out << "property double y\n";
		out << "property double z\n";
		out << "end_header\n";

		// data
		for (auto p : samples)
			out << p.x << " " << p.y << " " << p.z << "\n";

		out.close();
	}

}