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

namespace {
	struct QuantizedPointKey {
		int ix, iy, iz;

		bool operator==(const QuantizedPointKey& other) const {
			return ix == other.ix && iy == other.iy && iz == other.iz;
		}
	};

	struct QuantizedPointKeyHash {
		std::size_t operator()(const QuantizedPointKey& k) const {
			std::size_t h1 = std::hash<int>{}(k.ix);
			std::size_t h2 = std::hash<int>{}(k.iy);
			std::size_t h3 = std::hash<int>{}(k.iz);
			return h1 ^ (h2 << 1) ^ (h3 << 2);
		}
	};

	inline QuantizedPointKey make_point_key(const Vector3& p, float tol = 1e-1f) {
		return QuantizedPointKey{
			int(std::round(p.x / tol)),
			int(std::round(p.y / tol)),
			int(std::round(p.z / tol))
		};
	}

	
	// --- Shortest-path method switch ---------------------------------------------------------------
	// NormalDijkstra is kept for legacy support. BidirectionalDijkstra is the
	// preferred method for source-target shortest-path computations.

	//constexpr KarstNSim::ShortestPathMethod K_SHORTEST_PATH_METHOD =
	//	KarstNSim::ShortestPathMethod::NormalDijkstra;

	constexpr KarstNSim::ShortestPathMethod K_SHORTEST_PATH_METHOD =
		KarstNSim::ShortestPathMethod::BidirectionalDijkstra;

	constexpr bool K_CLEAR_SHORTEST_PATH_PREPROCESSING_AFTER_SKELETON = true;
}

namespace KarstNSim {

	GraphOperations::GraphOperations() : CostGraph(0)
	{
	}

	// Samples the domain with Dwork's algorithm adapted in 3D
	void GraphOperations::SampleSpaceDwork(const Box& box,
		const bool use_density_property, const int k, const std::vector<float>& propdensity, const Surface& topo_surface) {

		//const clock_t time1 = clock();
		PointCloud centers2D = topo_surface.get_centers_cloud(2);

		float min_z_topo = topo_surface.get_boundbox_min().z;

		// 1 //
		//First, we get the Box characteristics (origin and dimensions in each direction as well as density property).
		//The density property must be defined on a grid with int(sqrt(3) / r_min) cells in each direction, with r_min the
		//minimum density defined in the whole grid (0<rmin<1, e.g., 0.01).

		const Vector3& ptmin = box.get_basis();
		const Vector3& ptmax = box.get_end();

		const Vector3 u = box.get_u();
		const Vector3 v = box.get_v();
		const Vector3 w = box.get_w();

		float delta_y = ptmax.y - ptmin.y;
		float delta_x = ptmax.x - ptmin.x;
		float delta_z = ptmax.z - ptmin.z;

		float min_xyz = std::min({ std::abs(delta_x), std::abs(delta_y), std::abs(delta_z) });
		float dim1, dim2, dim3;
		if (min_xyz == std::abs(delta_x)) {
			dim1 = 1;
			dim2 = std::abs(delta_y / delta_x);
			dim3 = std::abs(delta_z / delta_x);
		}
		else if (min_xyz == std::abs(delta_y)) {
			dim1 = std::abs(delta_x / delta_y);
			dim2 = 1;
			dim3 = std::abs(delta_z / delta_y);
		}
		else {
			dim1 = std::abs(delta_x / delta_z);
			dim2 = std::abs(delta_y / delta_z);
			dim3 = 1;
		}

		Vector3 min(ptmin.x, ptmin.y, ptmin.z);
		float r_min = 0;
		int nu = box.get_nu(), nv = box.get_nv(), nw = box.get_nw(); //number of samples in u,v,w axis

		float r_cst = 0;
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
			std::cout << "* Using FVDPSS algorithm (Fast Variable Density Poisson-Sphere Sampling, adapted from Dwork et al. (2021)) with r_min = " << r_min << std::endl;
		}
		else { // else we keep a constant density everywhere
			r_cst = params.graphPoissonRadius;
			r_min = r_cst;
			std::cout << "* Using Uniform Density Poisson-Sphere Sampling, with r = " << r_min << std::endl;
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
		Vector3 pos_init;
		int64_t flattened_index = 0;
		int64_t converted_flat_index = 0;
		const int max_attempts = 10000;
		int attempt_count = 0;
		if (use_density_property) {
			while (!init_pt_in_grid && attempt_count < max_attempts) {
				pos_init.x = generateRandomFloat(0, dim1);
				pos_init.y = generateRandomFloat(0, dim2);
				pos_init.z = generateRandomFloat(0, dim3); // first random point of index 1
				flattened_index = flatten_indices2(dimx, dimy, dimz, dim1, dim2, dim3, pos_init); // we get the flattened index of the new point
				converted_flat_index = convert_index(dimx, dimy, dimz, nu, nv, nw, flattened_index);
				Vector3 real_pt;
				real_pt.x = min.x + (pos_init.x * u.x / dim1 + pos_init.y * v.x / dim2 + pos_init.z * w.x / dim3);
				real_pt.y = min.y + (pos_init.x * u.y / dim1 + pos_init.y * v.y / dim2 + pos_init.z * w.y / dim3);
				real_pt.z = min.z + (pos_init.x * u.z / dim1 + pos_init.y * v.z / dim2 + pos_init.z * w.z / dim3);
				// we make sure here that the 1st random point is indeed a) in the karstifiable zone (propdensity >0) and b) below topography
				if (propdensity[converted_flat_index] > 0 && CheckBelowSurf(real_pt, topo_surface, centers2D)) {
					init_pt_in_grid = true;
				}
				attempt_count++;
			}
			if (!init_pt_in_grid) {
				throw std::runtime_error("Failed to find a valid starting point in the grid after maximum attempts.");
			}
		}
		else {
			while (!init_pt_in_grid && attempt_count < max_attempts) {
				pos_init.x = generateRandomFloat(0, dim1);
				pos_init.y = generateRandomFloat(0, dim2);
				pos_init.z = generateRandomFloat(0, dim3);
				Vector3 real_pt;
				real_pt.x = min.x + (pos_init.x * u.x / dim1 + pos_init.y * v.x / dim2 + pos_init.z * w.x / dim3);
				real_pt.y = min.y + (pos_init.x * u.y / dim1 + pos_init.y * v.y / dim2 + pos_init.z * w.y / dim3);
				real_pt.z = min.z + (pos_init.x * u.z / dim1 + pos_init.y * v.z / dim2 + pos_init.z * w.z / dim3);
				if (CheckBelowSurf(real_pt, topo_surface, centers2D)) { // we make sure here that the 1st random point is indeed below topography
					init_pt_in_grid = true;
				}
				attempt_count++;
			}
			if (!init_pt_in_grid) {
				throw std::runtime_error("Failed to find a valid starting point below topography after maximum attempts.");
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
		//int gridSize = static_cast<int>(dimx) * dimy * dimz;

		// test to check that grid size is not too large (max due to vector max size):

		const int64_t MAX_SAFE_SIZE = 500000000;
		if (gridSize > MAX_SAFE_SIZE) {
			throw std::runtime_error("Grid size exceeds safe limit for vector allocation. Increase value for r_min.");
		}

		// VER 1 : unordered map
		//std::unordered_map<int, std::vector<int>> grid;
		// VER 2 : unordered multimap
		//std::multimap<int, int> grid;
		// VER 3 :  map
		//std::map<int, std::vector<int>> grid;

		std::vector<std::vector<int>> grid(gridSize);

		float r_local;
		if (use_density_property) {
			converted_flat_index = convert_index(dimx, dimy, dimz, nu, nv, nw, flattened_index);  // we convert this index to the property grid's own dimensions (probably different from the one we're using)
			r_local = propdensity[converted_flat_index];
		}
		else {
			r_local = r_cst;
		}
		std::vector<float> r_list({ r_local });
		std::vector<int64_t> scan_area_local = scan_area2(r_local, r_min, pos_init, dimx, dimy, dimz, dim1, dim2, dim3); // we get the scan area corresponding to the local density r_local
		for (int64_t x : scan_area_local) {
			grid[x].push_back(1);
			//grid.insert({ x, 1 });
		}
		//float time_2 = 0.;
		//float time_3 = 0.;
		//float time_4 = 0.;
		//float time_5 = 0.;
		//float time_6 = 0.;
		while (!active_list.empty()) { // while we can still add points
			// find a random index in the active list

			std::size_t random_idx_number = generateRandomIndex(active_list.size());
			int random_idx = active_list[random_idx_number];
			float r_active;
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
				Vector3 vec_dir(generateNormalRandom(0, 1), generateNormalRandom(0, 1), generateNormalRandom(0, 1));
				vec_dir = Normalize(vec_dir);// normalize the vector
				float magnitude_circle = std::cbrt(generateRandomFloat(r_active*r_active*r_active, 8 * r_active*r_active*r_active)); // the magnitude corresponds to the cube root of a value distributed uniformly between the two radii
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
				float r_k;
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
				//	float dist = magnitude(point_to_check - k_pt);
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
						float dist = magnitude(point_to_check - k_pt);
						if (dist < r_k) { // If distance is too small
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

					// last test before "validating" the point : checking that it is below the topographic surface
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
						//time_4 += float(time6 - time5);
						//time_5 += float(time7 - time6);
						//time_6 += float(time8 - time7);
					}
				}
				//time_2 += float(time4-time3);

			}
			if (crit_active == k) {  // If no k_pt was successfully added, we remove the root point from the active list
				active_list.erase(active_list.begin() + random_idx_number);
			}
		}

		grid.clear(); 	// Get rid of the map memory
	}

	bool GraphOperations::CheckBelowSurf(const Vector3& pt, const Surface& surface, const PointCloud& centers2D)
	{
		const int number_of_tries = 5;
		bool found_triangle = false;
		Vector3 P1, P2, P3;
		const Vector3 pt2D(pt.x, pt.y, 0.);
		std::vector<Neighbour> candidates;
		centers2D.findNearestNeighbors(pt2D, -1, number_of_tries, 1e25f, candidates);
		for (int test = 0; test < number_of_tries; test++) {
			const Triangle& trgl_i = surface.get_triangle(candidates[test].i);
			const int& pt1 = trgl_i.point(0), pt2 = trgl_i.point(1), pt3 = trgl_i.point(2);
			const Vector3& p1 = surface.get_node(pt1), p2 = surface.get_node(pt2), p3 = surface.get_node(pt3);
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

	float GraphOperations::distsurf(const Vector3& pt, const Surface& surface, const float& max_inception_surface_distance, const PointCloud& centers2D) const
	{
		const int number_of_tries = 5;
		bool found_triangle = false;
		Vector3 P1, P2, P3;
		const Vector3 pt2D(pt.x, pt.y, 0.);
		std::vector<Neighbour> candidates;
		centers2D.findNearestNeighbors(pt2D, -1, number_of_tries, 1e25f, candidates);

		for (int test = 0; test < number_of_tries; test++) {
			const Triangle& trgl_i = surface.get_triangle(candidates[test].i);
			const int& pt1 = trgl_i.point(0), pt2 = trgl_i.point(1), pt3 = trgl_i.point(2);
			const Vector3& p1 = surface.get_node(pt1), p2 = surface.get_node(pt2), p3 = surface.get_node(pt3);

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
		float distance = barycentric(pt, P1, P2, P3);

		if (distance < 1e-2f) {
			distance = 0.0f;
		}

		return distance;
	}

	void GraphOperations::DefineWSurfaceFlags(const std::vector<Surface>& water_tables) {

		// First and foremost, a water table doesn't project on the whole domain all the time.
		// Therefore any point outside of the projected perimeter of the water table can be considered in the vadose zone.
		std::vector<PointCloud> all_centers2D;

		std::vector<std::vector<float>> coord_boundbox_wt(params.nb_wt, std::vector<float>(4)); // elements are sorted as followed : xmin,ymin,xmax,ymax
		for (int i = 0; i < water_tables.size(); i++) {
			Vector3 min_wt = water_tables[i].get_boundbox_min();
			coord_boundbox_wt[i][0] = min_wt.x;
			coord_boundbox_wt[i][1] = min_wt.y;
			Vector3 max_wt = water_tables[i].get_boundbox_max();
			coord_boundbox_wt[i][2] = max_wt.x;
			coord_boundbox_wt[i][3] = max_wt.y;
			PointCloud centers2D = water_tables[i].get_centers_cloud(2);
			all_centers2D.push_back(centers2D);
		}
		samples_surf_flags.resize(int(samples.size()), params.nb_wt);

				// Ensure all flags are initialized as vadose by default.
		// Real water-table are filled below. The optional no-water-table
		// field intentionally remains false everywhere.
		for (int i = 0; i < int(samples.size()); ++i) {
			for (int j = 0; j < params.nb_wt; ++j) {
				samples_surf_flags(i, j) = false;
			}
		}
		for (int i = 0; i<int(samples.size()); i++) {

			for (int j = 0; j < water_tables.size(); j++) {
				bool check;
				if ((samples[i].x < coord_boundbox_wt[j][0]) || (samples[i].x > coord_boundbox_wt[j][2]) || (samples[i].y < coord_boundbox_wt[j][1]) || (samples[i].y > coord_boundbox_wt[j][3])) { // out of wt bonds
					check = false;
				}
				else {
					const Surface& surf_j = water_tables[j];

					check = CheckBelowSurf(samples[i], surf_j, all_centers2D[j]);
				}
				samples_surf_flags(i, j) = check;
			}
		}
	}

	void GraphOperations::compute_euclidian_distance_from_wt(const std::vector<Surface>& surfaces, const float max_inception_surface_distance) {

		std::vector<PointCloud> all_centers2D;
		for (int i = 0; i < surfaces.size(); i++) {
			PointCloud centers2D = surfaces[i].get_centers_cloud(2);
			all_centers2D.push_back(centers2D);
		}

		samples_wt_dist.resize(int(samples.size()), params.nb_wt);

		for (int i = 0; i < int(samples.size()); ++i) {
			for (int j = 0; j < params.nb_wt; ++j) {
				samples_wt_dist(i, j) = max_inception_surface_distance;
			}
		}
		std::vector<float> max_dist_wt(params.nb_wt, 0.);
		std::vector<float> min_dist_wt(params.nb_wt, std::numeric_limits<float>::infinity());
		for (int j = 0; j < surfaces.size(); j++) {
			const Surface& surf_j = surfaces[j];
			for (int i = 0; i<int(samples.size()); i++) {
				float distance = distsurf(samples[i], surf_j, max_inception_surface_distance, all_centers2D[j]);
				samples_wt_dist[i][j] = distance;
			}
		}
	}

	void GraphOperations::compute_euclidian_distance_from_surfaces(const std::vector<Surface>& surfaces, const float max_inception_surface_distance) {

		std::vector<PointCloud> all_centers2D;
		for (int i = 0; i < surfaces.size(); i++) {
			PointCloud centers2D = surfaces[i].get_centers_cloud(2);
			all_centers2D.push_back(centers2D);
		}

		samples_surf_dist.resize(int(samples.size()));
		for (int i = 0; i<int(samples.size()); i++) {
			float max_dist_surf_iter = std::numeric_limits<float>::infinity();
			for (int j = 0; j < surfaces.size(); j++) {
				const Surface& surf_j = surfaces[j];
				float distance = distsurf(samples[i], surf_j, max_inception_surface_distance, all_centers2D[j]);
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
			flags.push_back(int(samples.size()) - 1);
			//}
		}
		else { // if it has already been added, add the corresponding index
			// Check if flags is not empty
			//if (!flags.empty()) {
			flags.push_back(test_already_defined - 1);
			//}
		}
		//samples.push_back(node);
	}

	void GraphOperations::save_nghb_graph(const GeologicalParameters& geologicalparams, const bool& create_nghb_graph_property) {

		std::string full_dir_name = geologicalparams.directoryname + "/outputs";
		std::string full_name = params.scenename + "_nghb_graph.txt";
		std::vector<Segment> nghb_graph_line;
		std::vector<std::vector<float>> list_cost_up;
		std::vector<std::vector<float>> list_cost_down;
		list_cost_up.resize(params.nb_wt);
		list_cost_down.resize(params.nb_wt);
		for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
			for (int i = 0; i < adj.size(); i++) {
				for (int j = 0; j < adj[i].size(); j++) {
					if (create_nghb_graph_property 	&& adj[i][j].target >= 0) { // avoid empty neighbors (when n<N)
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
			std::vector<std::vector<std::vector<float>>> properties;
			properties.resize(number_segs);
			for (int segi = 0; segi < int(properties.size()); segi++) {
				properties[segi].resize(params.nb_wt * 2, std::vector<float>(2)); // there are 2*number of springs properties (one for each cost)
			}
			for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {

				std::string name_line_up = "cost_down" + std::to_string(cost_i + 1);
				std::string name_line_down = "cost_up" + std::to_string(cost_i + 1);
				property_names.push_back(name_line_up);
				property_names.push_back(name_line_down);

				for (int n = 0; n < number_segs; n++) {
					properties[n][2 * cost_i][0] = list_cost_up[cost_i][2 * n + 1];
					properties[n][2 * cost_i][1] = list_cost_up[cost_i][2 * n + 1];
				}

				for (int n = 0; n < number_segs; n++) {
					properties[n][2 * cost_i + 1][0] = list_cost_down[cost_i][2 * n + 1];
					properties[n][2 * cost_i + 1][1] = list_cost_down[cost_i][2 * n + 1];
				}

			}
			save_line(full_name, full_dir_name, full_line, property_names, properties); // version with costs
		}
		else {
			save_line(full_name, full_dir_name, full_line); // version without costs
		}
	}


	/**
	 * @brief Add all nodes from an input neighbor graph to the sampling cloud and
	 *        store, for each segment, the indices of its two endpoints in samples.
	 *
	 * @param[in]  input_nghb_graph Input neighbor graph provided by the user.
	 * @param[out] node_to_sample_idx  indices in samples.
	 */
	void GraphOperations::register_input_nghb_graph_nodes(
		const InputGraph& input_nghb_graph,
		std::vector<int>& input_graph_node_idx_in_samples)
	{
		input_graph_node_idx_in_samples.assign(input_nghb_graph.get_nb_nodes(), -1);

		const float eps_x = 1e-1f;
		const float eps_y = 1e-1f;
		const float eps_z = 1e-3f;

		std::unordered_map<QuantizedPointKey, std::vector<int>, QuantizedPointKeyHash> sample_cells;
		sample_cells.reserve(samples.size() * 2);

		for (int i = 0; i < int(samples.size()); ++i) {
			sample_cells[make_point_key(samples[i], eps_x)].push_back(i);
		}

		for (int i = 0; i < input_nghb_graph.get_nb_nodes(); ++i) {
			const Vector3& p = input_nghb_graph.get_node(i);
			const QuantizedPointKey base_key = make_point_key(p, eps_x);

			int matched_idx = -1;

			for (int dx = -1; dx <= 1 && matched_idx < 0; ++dx) {
				for (int dy = -1; dy <= 1 && matched_idx < 0; ++dy) {
					for (int dz = -1; dz <= 1 && matched_idx < 0; ++dz) {

						QuantizedPointKey neighbor_key{
							base_key.ix + dx,
							base_key.iy + dy,
							base_key.iz + dz
						};

						auto it = sample_cells.find(neighbor_key);
						if (it == sample_cells.end()) {
							continue;
						}

						const std::vector<int>& candidate_indices = it->second;
						for (int candidate_idx : candidate_indices) {
							const Vector3& q = samples[candidate_idx];

							if (std::abs(p.x - q.x) < eps_x &&
								std::abs(p.y - q.y) < eps_y &&
								std::abs(p.z - q.z) < eps_z) {
								matched_idx = candidate_idx;
								break;
							}
						}
					}
				}
			}

			input_graph_node_idx_in_samples[i] = matched_idx;

		}
	}

	void GraphOperations::InitializeCostGraph(const bool create_nghb_graph, const bool create_nghb_graph_property, const std::vector<KeyPoint>& keypts,
		const GeologicalParameters& geologicalparams, std::vector<Vector3>& nodes_on_inception_surfaces, std::vector<std::vector<Vector3>>& nodes_on_wt_surfaces,
		const std::vector<Surface>& inception_horizons, const std::vector<Surface>& water_tables, const bool use_sampling_points,
		const Box& box, const float max_inception_surface_distance,
		const std::vector<Vector3>& sampling_points, const bool create_vset_sampling, const bool use_density_property, const int k_pts,
		const float fraction_old_karst_perm, const std::vector<float>& propdensity, const std::vector<float>& propikp, const Surface& topo_surface)
	{
		PointCloud centers2D = topo_surface.get_centers_cloud(2);
		const clock_t time1 = clock();
		params = geologicalparams;
		int nu = box.get_nu(), nv = box.get_nv();
		std::vector<std::vector<int>> on_wt_flags_idx(params.nb_wt);
		std::vector<std::vector<int>> on_surf_flags_idx(1);
		Vector3 new_pt;

		// 1 Sample points in the domain

		std::cout << "1.1 Generating sampling cloud:\n";

		// 1.1 Key points are added to sample set
		for (int i = 0; i < keypts.size(); i++) {
			samples.push_back(keypts[i].p);
		}

		// 1.2a we recover previously sampled points if any
		if (use_sampling_points) {

			const float eps_x = 1e-1f;
			const float eps_y = 1e-1f;
			const float eps_z = 1e-3f;

			std::unordered_map<QuantizedPointKey, std::vector<int>, QuantizedPointKeyHash> sample_cells;
			sample_cells.reserve(samples.size() * 2);

			for (int i = 0; i < int(samples.size()); ++i) {
				sample_cells[make_point_key(samples[i], eps_x)].push_back(i);
			}

			for (int i = 0; i < sampling_points.size(); i++) {
				new_pt.x = sampling_points.at(i).x;
				new_pt.y = sampling_points.at(i).y;
				new_pt.z = sampling_points.at(i).z;

				const QuantizedPointKey base_key = make_point_key(new_pt, eps_x);
				bool already_exists = false;

				// Check the current cell and its 26 neighbors.
				for (int dx = -1; dx <= 1 && !already_exists; ++dx) {
					for (int dy = -1; dy <= 1 && !already_exists; ++dy) {
						for (int dz = -1; dz <= 1 && !already_exists; ++dz) {

							QuantizedPointKey neighbor_key{
								base_key.ix + dx,
								base_key.iy + dy,
								base_key.iz + dz
							};

							auto it = sample_cells.find(neighbor_key);
							if (it == sample_cells.end()) {
								continue;
							}

							const std::vector<int>& candidate_indices = it->second;
							for (int candidate_idx : candidate_indices) {
								const Vector3& p = samples[candidate_idx];

								if (std::abs(p.x - new_pt.x) < eps_x &&
									std::abs(p.y - new_pt.y) < eps_y &&
									std::abs(p.z - new_pt.z) < eps_z) {
									already_exists = true;
									break;
								}
							}
						}
					}
				}

				if (!already_exists) {
					samples.push_back(new_pt);
					sample_cells[base_key].push_back(int(samples.size()) - 1);
				}
			}
			std::cout << "* Using " << samples.size() << " user-defined sampling points" << std::endl;
		}

		// 1.2b else we sample them with modified Dwork poisson sphere sampling
		else {
			SampleSpaceDwork(box, use_density_property, k_pts, propdensity, topo_surface);
			std::cout << "* Using " << samples.size() << " automatically sampled points" << std::endl;
		}

		// 1.3 remove points sampled inside the "no-karst" spheres
		// /!\ NEVER make a no-karst sphere and a keypoint intersect (makes code crash)
		size_t old_sample_size = samples.size();
		PointCloud centers_sampling_pts(samples);
		std::vector<Neighbour> all_candidates;
		for (int sphere = 0; sphere < params.spheres.size(); sphere++) {
			std::vector<Neighbour> candidates;
			float max_dist = params.spheres[sphere].radius * params.spheres[sphere].radius;
			centers_sampling_pts.findNearestNeighbors(params.spheres[sphere].center, -1, 10000, max_dist, candidates);
			all_candidates.insert(all_candidates.end(), candidates.begin(), candidates.end());
		}
		remove_neighbors(samples, all_candidates);
		// remove inception surface points inside "no-karst" spheres
		PointCloud centers_inception_pts(nodes_on_inception_surfaces);
		std::vector<Neighbour> all_candidates_inception;
		for (int sphere = 0; sphere < params.spheres.size(); sphere++) {
			std::vector<Neighbour> candidates;
			float max_dist = params.spheres[sphere].radius * params.spheres[sphere].radius;
			centers_inception_pts.findNearestNeighbors(params.spheres[sphere].center, -1, 1000, max_dist, candidates);
			all_candidates_inception.insert(all_candidates_inception.end(), candidates.begin(), candidates.end());
		}
		remove_neighbors(nodes_on_inception_surfaces, all_candidates_inception);
		// remove water table surface points inside "no-karst" spheres
		for (int wt = 0; wt < nodes_on_wt_surfaces.size(); wt++) {
			PointCloud centers_wt_pts(nodes_on_wt_surfaces[wt]);
			std::vector<Neighbour> all_candidates_wt;
			for (int sphere = 0; sphere < params.spheres.size(); sphere++) {
				std::vector<Neighbour> candidates;
				float max_dist = params.spheres[sphere].radius * params.spheres[sphere].radius;
				centers_wt_pts.findNearestNeighbors(params.spheres[sphere].center, -1, 1000, max_dist, candidates);
				all_candidates_wt.insert(all_candidates_wt.end(), candidates.begin(), candidates.end());
			}
			remove_neighbors(nodes_on_wt_surfaces[wt], all_candidates_wt);
		}
		std::cout << "* Removed " << old_sample_size - samples.size() << " points from the sampling set because they were inside no-karst spheres" << std::endl;

		// 1.4 Add inception surface points to samples
		old_sample_size = samples.size();
		for (int i = 0; i < nodes_on_inception_surfaces.size(); i++) {
			float value_itr = 1.;
			if (!use_sampling_points && use_density_property) {
				int u1, v1, w1;
				box.xyz2uvw_with_limits_conditions(nodes_on_inception_surfaces[i], u1, v1, w1);
				value_itr = propdensity[u1 + nu * v1 + nv * nu*w1];
			}
			bool below_top_test = true;
			if (!topo_surface.is_empty()) {
				below_top_test = CheckBelowSurf(nodes_on_inception_surfaces[i], topo_surface, centers2D);
			}
			if (below_top_test && value_itr > 0) {  // Accept only nodes that are below topo surf AND are associated to a karstifiable zone (ie. propdensity >0)
				process_node_in_samples(nodes_on_inception_surfaces[i], on_surf_flags_idx[0]);
			}
		}
		std::cout << "* Added " << samples.size() - old_sample_size << " points from the inception surfaces to the sampling set" << std::endl;

		// 1.5 Add spring points to physical water-table flags.
		// Springs with wt_idx == 0 have no associated water table and must not be
		// inserted into any physical water-table flag list.
		old_sample_size = samples.size();
		for (int i = 0; i < keypts.size(); i++) { // first the keypoints
			if (keypts[i].type == KeyPointType::Spring) {
				auto it = std::find(samples.begin(), samples.end(), keypts[i].p);
				int index = (it != samples.end()) ? std::distance(samples.begin(), it) : -1;

				if (index < 0) {
					continue;
				}

				if (keypts[i].wt_idx == 0) {
					continue;
				}

				const int wt_channel = keypts[i].wt_idx - 1;

				if (wt_channel < 0 || wt_channel >= static_cast<int>(on_wt_flags_idx.size())) {
					continue;
				}

				on_wt_flags_idx.at(wt_channel).push_back(index);
			}
		}
		for (int i = 0; i < nodes_on_wt_surfaces.size(); i++) { // iterate on nb of wt surfaces
			for (int j = 0; j < nodes_on_wt_surfaces[i].size(); j++) { // iterate on points of wt i
				// Accept only nodes that are below topo surf AND are associated to a karstifiable zone (ie. propdensity >0)
				float value_itr = 1.;
				if (!use_sampling_points && use_density_property) {
					int u1, v1, w1;
					box.xyz2uvw_with_limits_conditions(nodes_on_wt_surfaces[i][j], u1, v1, w1);
					value_itr = propdensity[u1 + nu * v1 + nv * nu*w1];
				}
				bool below_top_test = true;
				if (!topo_surface.is_empty()) {
					below_top_test = CheckBelowSurf(nodes_on_wt_surfaces[i][j], topo_surface, centers2D);
				}
				if (below_top_test && value_itr > 0) {  // Accept only nodes that are below topo surf AND are associated to a karstifiable zone (ie. propdensity >0)
					process_node_in_samples(nodes_on_wt_surfaces[i][j], on_wt_flags_idx[i]);
				}
			}
		}
		std::cout << "* Added " << samples.size() - old_sample_size << " points from the water table surfaces to the sampling set" << std::endl;

		// 1.6 Add previously simulated karst networks samples
		old_sample_size = samples.size();
		for (int i = 0; i < params.PtsOldGraph.size(); i++) {
			std::vector<int> idx_i;
			for (int j = 0; j < 2; j++) { // =2 (start and end points of each segment of the old graph)
				process_node_in_samples(params.PtsOldGraph[i][j], idx_i);
			}
			params.IdxOldGraph.set_row(i, idx_i);
		}
		std::cout << "* Added " << samples.size() - old_sample_size << " points from the previously simulated karst networks to the sampling set" << std::endl;
		//IdxOldGraph = params.IdxOldGraph;
		std::vector<std::string> noise_prop_name;
		if (params.use_noise) {
			std::pair<std::vector<float>, std::vector<std::string>>(noise_vector, noise_prop_name) = Save_noise_vector(samples);
		}
		if (create_vset_sampling) {
			std::string full_name = geologicalparams.scenename + "_pts.txt";
			std::string full_dir_name = geologicalparams.directoryname + "/outputs";

			if (params.use_noise) {
				std::vector<std::vector<float>> noise_prop;
				for (const float& value : noise_vector) {
					noise_prop.push_back(std::vector<float>{value});
				}
				save_pointset(full_name, full_dir_name, samples, noise_prop_name, noise_prop);
			}
			else {
				save_pointset(full_name, full_dir_name, samples);
			}
		}

		// 1.7 Add user-defined neighbor graph samples
		std::vector<int> input_graph_node_idx_in_samples;
		if (params.use_input_nghb_graph) {
			register_input_nghb_graph_nodes(params.input_nghb_graph, input_graph_node_idx_in_samples);
		}

		const clock_t time2 = clock();
		std::cout << "Points were successfully sampled / recovered ("
			<< std::fixed << std::setprecision(3)
			<< float(time2 - time1) / CLOCKS_PER_SEC << " s)" << std::endl;
		std::cout << "1.2 Storing properties on the sampling cloud:"<<std::endl;
		// 2 Find the karstification potential for each point of the sampling cloud
		if (params.karstificationCost.used) {
			const Vector3 u = box.get_u();
			const Vector3 v = box.get_v();
			const Vector3 w = box.get_w();
			int nu = box.get_nu(), nv = box.get_nv();
			for (int i = 0; i < samples.size(); i++) {
				int u1, v1, w1;
				box.xyz2uvw_with_limits_conditions(samples[i], u1, v1, w1);
				float value_itr = propikp[u1 + nu * v1 + nv * nu*w1];
				if (value_itr < 0) { // a value smaller than  means a No data value, which means a cell above topography
					samples_layer_kp.push_back(1); // maximal cost (this happens because of the discrepancy between the background grid's resolution and the points' position close to the topography)
					continue;
				}
				samples_layer_kp.push_back(1 - propikp[u1 + nu * v1 + nv * nu*w1]);
			}
		}

		const clock_t time21 = clock();

		std::cout << "* Karstification potential was recovered if needed ("
			<< std::fixed << std::setprecision(3)
			<< float(time21 - time2) / CLOCKS_PER_SEC << " s)" << std::endl;

		// 3 For each sampling point, this defines whether the point is below or above each water table surface, and also if they are exactly onto the wt surface
		DefineWSurfaceFlags(water_tables); // above or under wt test
		get_flags_from_indices(samples, on_wt_flags_idx, samples_on_wt_flags); // exactly on wt test

		const clock_t time22 = clock();

		std::cout << "* Points were flagged in vadose and phreatic zones ("
			<< std::fixed << std::setprecision(3)
			<< float(time22 - time21) / CLOCKS_PER_SEC << " s)" << std::endl;

		// 4 Compute the distance from each point to the inception surfaces and the water table surfaces
		compute_euclidian_distance_from_wt(water_tables, max_inception_surface_distance);
		if (!inception_horizons.empty()) {
			compute_euclidian_distance_from_surfaces(inception_horizons, max_inception_surface_distance);
		}
		const clock_t time23 = clock();

		std::cout << "* Distance from points to horizons and water tables was computed ("
			<< std::fixed << std::setprecision(3)
			<< float(time23 - time22) / CLOCKS_PER_SEC << " s)" << std::endl;

		const clock_t time3 = clock();

		std::cout << "Sampling points geometry fully analyzed to create graph ("
			<< float(time3 - time2) / CLOCKS_PER_SEC << " s)" << std::endl;

		std::cout <<"1.3 Creation of nearest neighbor cost graph:"<<std::endl;

		// 5 Nearest neighbour graph creation (initialization + cost computation on each edge)

		if (params.use_input_nghb_graph) {
			BuildNearestNeighbourGraphFromInput(
				box,
				params.input_nghb_graph,
				input_graph_node_idx_in_samples,
				fraction_old_karst_perm,
				max_inception_surface_distance,
				keypts);
		}
		else {
			BuildNearestNeighbourGraph(
				box,
				fraction_old_karst_perm,
				max_inception_surface_distance,
				keypts);
		}

		const clock_t time4 = clock();

		std::cout << "Nearest neighbor & cost graph generated ("
			<< std::fixed << std::setprecision(3)
			<< float(time4 - time3) / CLOCKS_PER_SEC << " s)" << std::endl;

		// 6 Save the nearest neighbour graph
		if (create_nghb_graph) {
			save_nghb_graph(geologicalparams, create_nghb_graph_property);
		}
	}

	std::pair<std::vector<float>, std::vector<std::string>> GraphOperations::Save_noise_vector(std::vector<Vector3> Points) {
		std::vector<std::string> noise_prop_name;
		noise_prop_name.push_back("Noise");
		SimplexNoise Perlin = SimplexNoise(float(params.noise_frequency), 1.0f, 2.0f, 0.5f);
		noise_vector.resize(Points.size());
		float normalization_term = 100.0f;

		for (int i = 0; i < Points.size(); i++) {
			noise_vector[i] = Perlin.fractal(params.noise_octaves, Points[i].x / normalization_term, Points[i].y / normalization_term, Points[i].z / normalization_term);
		}
		std::pair<std::vector<float>, std::vector<std::string>> pair;
		pair.first = noise_vector;
		pair.second = noise_prop_name;
		return pair;
	}

	void GraphOperations::add_noise() {

		//use normalization to avoid too high values in the noise_weight parameter insteaed of value between 0 and 2
		//NormalizeWeights();

		if (params.use_noise) {
			for (int i = 0; i < adj.size(); i++) {
				for (int j = 0; j < adj[i].size(); j++) {
					if (adj[i][j].target >= 0) { // avoid empty neighbors (when n<N) 
					//Vector3 pn = samples[adj[i][j].target];
						float weight1 = 2 + noise_vector[i];
						float weight2 = 2 + noise_vector[adj[i][j].target];
						for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
							adj[i][j].weight[cost_i] += params.noise_weight*(weight1 + weight2) / 2; //supposed to divide by adj[i][j].size() but noise_weight needs to be adjusted to higher values
						}
					}
				}
			}
			//NormalizeWeights();
		}
	}

	void GraphOperations::BuildNearestNeighbourGraph(const Box& box, const float fraction_old_karst_perm, float max_dist_surf, const std::vector<KeyPoint>& keypts)
	{

		float time_neigbhors = 0.;
		float time_costs = 0.;

		const Vector3& ptmin = box.get_basis();
		const Vector3& ptmax = box.get_end();
		const Vector3 u = box.get_u();
		const Vector3 v = box.get_v();
		const Vector3 w = box.get_w();
		float delta_y_abs = std::abs(ptmax.y - ptmin.y);
		float delta_x_abs = std::abs(ptmax.x - ptmin.x);
		float delta_z_abs = std::abs(ptmax.z - ptmin.z);
		float dim_max = std::max({ delta_y_abs, delta_x_abs, delta_z_abs });
		float one_over_delta_y = 1 / (delta_y_abs);
		float one_over_delta_x = 1 / (delta_x_abs);
		float one_over_delta_z = 1 / (delta_z_abs);
		const int N = params.graphNeighbourCount;
		float dmin = std::numeric_limits<float>::infinity();
		float dmax = 0.;

		adj.resize(samples.size(), N); // note that all elements take default constructor GraphEdge as initial value, which has for target -1
		float R = 0.; // std::numeric_limits<float>::infinity();
		float distmax = 1e25f;
		if (params.graphuse_max_nghb_radius) {
			R = params.graphNeighbourRadius * params.graphNeighbourRadius; // Precise the value of R if it has been defined by the user
			distmax = (R / (dim_max*dim_max));
		}

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

			pointCloud.findNearestNeighbors(p_stretched, i, N, distmax, candidates);

			const clock_t time1 = clock();
			int n = KarstNSim::Min((int)candidates.size(), N);

			for (int j = 0; j < n; j++) { // note here that only the n (<=N) first elements of adj[i] will be computed. The others keep default value -1 (meaning not a neighbor)
				const Vector3 pn = samples[candidates[j].i];
				const std::pair<std::vector<float>, std::vector<bool>> pair = ComputeEdgeCost(i, p, pn, dmin, dmax, max_dist_surf);
				SetEdge(i, j, candidates[j].i, pair.first, pair.second);
			}
			const clock_t time2 = clock();
			time_neigbhors += float(time1 - time0);
			time_costs += float(time2 - time1);
		}
		const clock_t time3 = clock();
		for (int i = 0; i < adj.size(); i++) {
			Vector3 p = samples[i];

			for (int j = 0; j < adj[i].size(); j++) {
				if (adj[i][j].target >= 0) {
					Vector3 pn = samples[adj[i][j].target];

					const float distance_cost = ComputeEdgeDistanceCost(p, pn, dmin, dmax);

					if (distance_cost <= 0.0f) {
						throw std::runtime_error(
							"[cost] Zero edge length encountered during graph-cost computation. "
							"This usually indicates duplicated sampling points or a degenerate edge."
						);
					}

					for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
						adj[i][j].weight[cost_i] *= distance_cost;
					}
				}
			}
		}

		if (params.use_noise_on_all) {
			add_noise();
		}

		// reduce cost of edges already karstified in previous iterations
		for (int i = 0; i < params.IdxOldGraph.size(); i++) {
			for (int j = 1; j < 2; j++) { // go only once on each edge
				for (int k = 0; k < adj[params.IdxOldGraph[i][0]].size(); k++) {
					if (adj[params.IdxOldGraph[i][0]][k].target == params.IdxOldGraph[i][j]) { // if an edge was created between two nodes of the previous karst
						for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
							adj[params.IdxOldGraph[i][0]][k].weight[cost_i] *= fraction_old_karst_perm;
						}
						break; // the edge was found, no need to scan further
					}
				}
			}
		}

		if (!params.waypointsimpactradius.empty()) {
			PointCloud pointCloud_not_stretched(samples);
			// reduce cost of edges inside spheres around waypoints
			for (int i = 0; i < params.waypointsimpactradius.size(); i++) {
				float radius = params.waypointsimpactradius[i].prop;
				Vector3 waypt_pos = keypts[params.waypointsimpactradius[i].index].p; // get position of waypoint i
				std::vector<Neighbour> candidates;
				float max_dist = radius * radius;
				pointCloud_not_stretched.findNearestNeighbors(waypt_pos, -1, 10000, max_dist, candidates);
				for (int j = 0; j < candidates.size(); j++) {
					for (int k = 0; k < adj[candidates[j].i].size(); k++) { // iterate on edges starting from point candidates[j].i
						if (adj[candidates[j].i][k].target >= 0) { // avoid empty neighbors (when n<N) 
							for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
								adj[candidates[j].i][k].weight[cost_i] *= params.waypoints_weight; // reduce the corresponding costs
							}
						}
					}
				}
			}
		}

		const clock_t time4 = clock();
		time_costs += float(time4 - time3);
		// At that point, the adjacency graph is entirely built, and ready to be used.
		// We keep a copy of adj to use later on for amplification, which works better without cost reduction

		//adj_no_cost_red = adj; // copy

	}

	void GraphOperations::BuildNearestNeighbourGraphFromInput(
		const Box& box,
		const InputGraph& input_nghb_graph,
		const std::vector<int>& input_graph_node_idx_in_samples,
		const float fraction_old_karst_perm,
		float max_dist_surf,
		const std::vector<KeyPoint>& keypts)
	{
		(void)box;

		float time_neigbhors = 0.;
		float time_costs = 0.;

		float dmin = std::numeric_limits<float>::infinity();
		float dmax = 0.;

		// ---------------------------------------------------------------------
		// 1) Rebuild a clean undirected edge set from the imported input graph
		// ---------------------------------------------------------------------

		auto make_edge_key = [](int a, int b) -> std::uint64_t {
			if (a > b) {
				std::swap(a, b);
			}
			return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32)
				| static_cast<std::uint32_t>(b);
		};

		std::unordered_set<std::uint64_t> unique_edges;
		std::vector<std::pair<int, int>> edges_idx;
		std::vector<int> degree(samples.size(), 0);

		const int nb_input_edges = input_nghb_graph.get_nb_edges();

		for (int i = 0; i < nb_input_edges; ++i) {

			const Vector2i& edge = input_nghb_graph.get_edge(i);
			const int node_a = edge.x;
			const int node_b = edge.y;

			if (node_a < 0 || node_b < 0) {
				continue;
			}
			if (node_a >= int(input_graph_node_idx_in_samples.size()) ||
				node_b >= int(input_graph_node_idx_in_samples.size())) {
				continue;
			}

			const int a = input_graph_node_idx_in_samples[node_a];
			const int b = input_graph_node_idx_in_samples[node_b];

			if (a < 0 || b < 0) {
				continue;
			}
			if (a >= int(samples.size()) || b >= int(samples.size())) {
				continue;
			}
			if (a == b) {
				continue;
			}
			if (samples[a] == samples[b]) {
				continue;
			}

			const std::uint64_t edge_key = make_edge_key(a, b);
			if (!unique_edges.insert(edge_key).second) {
				continue;
			}

			edges_idx.push_back({ a, b });
			degree[a]++;
			degree[b]++;
		}

		int max_degree = 0;
		for (int i = 0; i < int(degree.size()); ++i) {
			max_degree = std::max(max_degree, degree[i]);
		}
		max_degree = std::max(1, max_degree);

		adj.resize(samples.size(), max_degree);

		// For each node, stores the next free adjacency slot
		std::vector<int> next_slot(samples.size(), 0);

		// ---------------------------------------------------------------------
		// 2) Build adjacency and compute geological costs exactly like in the
		//    standard nearest-neighbor graph constructor
		// ---------------------------------------------------------------------

		const clock_t time0 = clock();

		for (int e = 0; e < int(edges_idx.size()); ++e) {
			const int ia = edges_idx[e].first;
			const int ib = edges_idx[e].second;

			const Vector3 p = samples[ia];
			const Vector3 pn = samples[ib];

			// Direction ia -> ib
			{
				const std::pair<std::vector<float>, std::vector<bool>> pair_ab =
					ComputeEdgeCost(ia, p, pn, dmin, dmax, max_dist_surf);
				SetEdge(ia, next_slot[ia], ib, pair_ab.first, pair_ab.second);
				next_slot[ia]++;
			}

			// Direction ib -> ia
			{
				const std::pair<std::vector<float>, std::vector<bool>> pair_ba =
					ComputeEdgeCost(ib, pn, p, dmin, dmax, max_dist_surf);
				SetEdge(ib, next_slot[ib], ia, pair_ba.first, pair_ba.second);
				next_slot[ib]++;
			}
		}

		const clock_t time1 = clock();
		time_neigbhors += float(time1 - time0);
		time_costs += float(time1 - time0);

		// ---------------------------------------------------------------------
		// 3) Add distance cost exactly as in BuildNearestNeighbourGraph()
		// ---------------------------------------------------------------------

		const clock_t time3 = clock();

		for (int i = 0; i < adj.size(); i++) {
			Vector3 p = samples[i];

			for (int j = 0; j < adj[i].size(); j++) {
				if (adj[i][j].target >= 0) {
					Vector3 pn = samples[adj[i][j].target];

					const float distance_cost = ComputeEdgeDistanceCost(p, pn, dmin, dmax);

					if (distance_cost <= 0.0f) {
						throw std::runtime_error(
							"[cost] Zero edge length encountered during graph-cost computation. "
							"This usually indicates duplicated sampling points or a degenerate edge."
						);
					}

					for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
						adj[i][j].weight[cost_i] *= distance_cost;
					}
				}
			}
		}

		// ---------------------------------------------------------------------
		// 4) Optional global noise, exactly as in BuildNearestNeighbourGraph()
		// ---------------------------------------------------------------------

		if (params.use_noise_on_all) {
			add_noise();
		}

		// ---------------------------------------------------------------------
		// 5) Reduce cost of edges already karstified in previous iterations
		//    exactly as in BuildNearestNeighbourGraph()
		// ---------------------------------------------------------------------

		for (int i = 0; i < params.IdxOldGraph.size(); i++) {
			for (int j = 1; j < 2; j++) { // go only once on each edge
				for (int k = 0; k < adj[params.IdxOldGraph[i][0]].size(); k++) {
					if (adj[params.IdxOldGraph[i][0]][k].target == params.IdxOldGraph[i][j]) {
						for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
							adj[params.IdxOldGraph[i][0]][k].weight[cost_i] *= fraction_old_karst_perm;
						}
						break;
					}
				}
			}
		}

		// ---------------------------------------------------------------------
		// 6) Reduce cost of edges inside waypoint influence spheres exactly as
		//    in BuildNearestNeighbourGraph()
		// ---------------------------------------------------------------------

		if (!params.waypointsimpactradius.empty()) {
			PointCloud pointCloud_not_stretched(samples);

			for (int i = 0; i < params.waypointsimpactradius.size(); i++) {
				float radius = params.waypointsimpactradius[i].prop;
				Vector3 waypt_pos = keypts[params.waypointsimpactradius[i].index].p;
				std::vector<Neighbour> candidates;
				float max_dist = radius * radius;
				pointCloud_not_stretched.findNearestNeighbors(waypt_pos, -1, 10000, max_dist, candidates);

				for (int j = 0; j < candidates.size(); j++) {
					for (int k = 0; k < adj[candidates[j].i].size(); k++) {
						if (adj[candidates[j].i][k].target >= 0) {
							for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
								adj[candidates[j].i][k].weight[cost_i] *= params.waypoints_weight;
							}
						}
					}
				}
			}
		}

		const clock_t time4 = clock();
		time_costs += float(time4 - time3);


		// At that point, the adjacency graph is entirely built, and ready to be used.
	}

	float GraphOperations::ComputeEdgeDistanceCost(
		const Vector3& p,
		const Vector3& pn,
		const float& dmin,
		const float& dmax)
	{
		(void)dmin;
		(void)dmax;

		const float distance = graph_anisotropic_distance(
			p,
			pn,
			params.vertical_distance_stretching_factor
		);

		if (!std::isfinite(distance)) {
			throw std::runtime_error(
				"[cost] Non-finite edge distance encountered during graph-cost computation."
			);
		}

		return distance * params.distanceCost.weight;
	}

	std::pair<std::vector<float>, std::vector<bool>> GraphOperations::ComputeEdgeCost(const int& index, const Vector3& p, const Vector3& pn,
		float& dmin, float& dmax, const float& dist_max) const
	{
		//auto time1 = std::chrono::high_resolution_clock::now();
		std::vector<float> cost(params.nb_wt, params.multiply_costs ? 1.0 : 1.0);
		std::vector<bool> frac_flags(params.fractures_orientations.size(), false);

		const Vector3 diff = pn - p;
		const float dist = graph_anisotropic_distance(
			p,
			pn,
			params.vertical_distance_stretching_factor
		);
		dmax = std::max(dmax, dist);
		dmin = std::min(dmin, dist);

		// Horizon cost
		if (params.horizonCost.used)
		{
			const float distance_to_surf = KarstNSim::Clamp(samples_surf_dist[index], 0.0, dist_max);
			float w = (1.0 - KarstNSim::CubicSmooth(distance_to_surf, dist_max));
			for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
				if (params.multiply_costs) cost[cost_i] *= w * params.horizonCost.weight;
				else cost[cost_i] += w * params.horizonCost.weight;
			}
		}

		// Permeability
		if (params.karstificationCost.used)
		{
			for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
				if (params.multiply_costs) cost[cost_i] *= samples_layer_kp[index] * params.karstificationCost.weight;
				else cost[cost_i] += samples_layer_kp[index] * params.karstificationCost.weight;
			}
		}

		// Fracture orientation
		if (params.fractureCost.used)
		{

			const float azimut = atan(diff[0] / diff[1]) * (180.0 / 3.141592653589793238463);
			int i;
			float costFrac = 1.;
			for (i = 0; i < params.fractures_orientations.size(); i++) {
				float orientation;
				if (params.fractures_orientations[i] >= 90) {
					orientation = params.fractures_orientations[i] - 180;
				}
				else {
					orientation = params.fractures_orientations[i];
				}
				float tolerance = params.fractures_tolerances[i];
				if (std::abs(azimut - orientation) < tolerance) {
					costFrac = 0.;
					frac_flags[i] = true;
					break;
				}
			}
			for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
				if (params.multiply_costs) cost[cost_i] *= costFrac * params.fractureCost.weight;
				else cost[cost_i] += costFrac * params.fractureCost.weight;
			}
		}

		// Vadose and phreatic zone

		for (int cost_i = 0; cost_i < params.nb_wt; ++cost_i) {
			if (!samples_surf_flags(index, cost_i)) {
				const float cos_angle = Dot(Vector3(0, 0, 1), Normalize(diff)); // from -1 to 1
				const float e = (cos_angle + 1) / 2; // from 0 to 1
				if (params.multiply_costs) cost[cost_i] *= params.waterTable1.weight * e;
				else cost[cost_i] += params.waterTable1.weight * e;
			}
			else {
				float w_prim = KarstNSim::Clamp(samples_wt_dist(index, cost_i), 0.0, dist_max);
				float w = (1.0 - KarstNSim::CubicSmooth(w_prim, dist_max));
				if (params.multiply_costs) cost[cost_i] *= w * params.waterTable2.weight;
				else cost[cost_i] += w * params.waterTable2.weight;
			}
		}

		// Clamp the cost value to keep it positive for Dijkstra (this might not be necessary)
		for (int cost_i = 0; cost_i < params.nb_wt; cost_i++) {
			cost[cost_i] = KarstNSim::Clamp(cost[cost_i], 0.0, cost[cost_i]); // make sure to keep a positive cost for Dijsktra
		}

		return std::make_pair(cost, frac_flags);
	}

	int GraphOperations::NodeIndex(const Vector3& p) const
	{
		for (int i = 0; i < samples.size(); i++)
		{
			if (samples[i] == p)
				return i;
		}
		return -1;
	}

	void GraphOperations::ComputeKarsticSkeleton(const std::vector<KeyPoint>& pts, const float fraction_karst_perm, std::vector<std::vector<int>>& pathsFinal, std::vector<std::vector<float>>& costsFinal, std::vector<std::vector<char>>& vadoseFinal,
		std::vector<int>& springidxFinal, bool save_new_connectivity_matrix)
	{
		float step1 = 0.0f, step2 = 0.0f, step3 = 0.0f, step4 = 0.0f, step5 = 0.0f, step6 = 0.0f, step7 = 0.0f, step8 = 0.0f, step_sink = 0.0f;
		float step_shortest_path = 0.0f;
		const clock_t time1 = clock();

		Array2D<int> new_connectivity_matrix = params.connectivity_matrix; // initially, we copy the old matrix. We will then update it to solve the uncertain '2' connections

		// Internal representation of key points, storing their index (with two sublists : one for sinks and one for springs)
		std::vector<InternalKeyPoint> keypts;
		std::vector<InternalKeyPoint> keyptssinks;
		std::vector<InternalKeyPoint> keyptssprings;
		float eps = 1e-10f;
		for (int i = 0; i < int(pts.size()); i++) {
			keypts.push_back({ NodeIndex(pts[i].p), pts[i].p, pts[i].type });
			if (pts[i].type == KeyPointType::Sink)
			{
				keyptssinks.push_back({ NodeIndex(pts[i].p), pts[i].p, pts[i].type });
			}
			else if (pts[i].type == KeyPointType::Spring) {
				keyptssprings.push_back({ NodeIndex(pts[i].p), pts[i].p, pts[i].type,pts[i].wt_idx });
			}
		}
		// Allocates all temporary arrays
		std::vector<std::vector<float>> all_distances;
		std::vector<std::vector<float>> all_selection_distances;
		std::vector<std::vector<std::vector<float>>> all_distances_detailled; // version where each cost for each segment is precised
		std::vector<std::vector<std::vector<int>>> all_paths;
		std::vector< std::vector<std::vector<char>>> all_vadose;
		all_paths.resize(keyptssinks.size());
		all_vadose.resize(keyptssinks.size());
		all_distances.resize(keyptssinks.size(), std::vector<float>(keyptssprings.size(), std::numeric_limits<float>::infinity()));
		all_selection_distances.resize(keyptssinks.size(), std::vector<float>(keyptssprings.size(), std::numeric_limits<float>::infinity()));
		all_distances_detailled.resize(keyptssinks.size());
		for (int i = 0; i < keyptssinks.size(); i++)
		{
			all_distances_detailled[i].resize(keyptssprings.size());
			all_paths[i].resize(keyptssprings.size());
			all_vadose[i].resize(keyptssprings.size());
		}
		// Iterate on all sinks points
		clock_t time1point5 = clock();

		step1 = float(time1point5 - time1) / CLOCKS_PER_SEC;

		const int total_sinks = static_cast<int>(keyptssinks.size());
		int next_progress_percent = 10;

		for (int i = 0; i < total_sinks; i++) {
			const clock_t time_sink = clock();

			int source = keyptssinks[i].index;

			const float z_inlet = keyptssinks[i].p.z;
			bool retained_path_for_sink = false;

			// Pick the inlet/outlet connections when val = 2 in connectivity matrix if the random connectivity option was chosen

			if (!params.allow_single_outlet) {
				std::vector<int> list_of_two;
				for (int j = 0; j < keyptssprings.size(); j++) {
					const int connectivity_flag = params.connectivity_matrix[params.sinks_index[i] - 1][j];
					if (connectivity_flag != 2) {
						continue;
					}

					// Reject outlets located above the inlet before random outlet selection.
					// This prevents single-outlet random mode from selecting a hydraulically
					// inadmissible outlet while discarding admissible lower outlets.
					if (z_inlet < keyptssprings[j].p.z) {
						new_connectivity_matrix[params.sinks_index[i] - 1][j] = 0;
						continue;
					}
						list_of_two.push_back(j);
				}
				if (int(list_of_two.size()) != 0) {
					std::size_t random_idx_number = generateRandomIndex(list_of_two.size());
					int random_idx = list_of_two[random_idx_number];
					for (int j = 0; j < keyptssprings.size(); j++) {

						int connectivity_flag = params.connectivity_matrix[params.sinks_index[i] - 1][j];

						if (connectivity_flag == 2 && j != random_idx) {
							params.connectivity_matrix[params.sinks_index[i] - 1][j] = 0;
							new_connectivity_matrix[params.sinks_index[i] - 1][j] = 0;
						}
						else if (connectivity_flag == 2 && j == random_idx) {
							params.connectivity_matrix[params.sinks_index[i] - 1][j] = 1;
							new_connectivity_matrix[params.sinks_index[i] - 1][j] = 1;
						}
					}
				}
			}

			// now iterate on springs, and look for path only if connectivity flag is equal to 1 (or 2 if shortest distance option was chosen)

			for (int j = 0; j < keyptssprings.size(); j++) {

				int connectivity_flag = params.connectivity_matrix[params.sinks_index[i] - 1][j];
				if (connectivity_flag == 0) continue; // skip springs for which flag equals 0 (no connection)
				// Reject hydraulically inadmissible inlet-outlet pairs before any shortest-path computation.
				// The comparison is strict: equal elevations are still allowed as a zero-gradient limit case.
				if (z_inlet < keyptssprings[j].p.z) {
					new_connectivity_matrix[params.sinks_index[i] - 1][j] = 0;
					continue;
				}

				int target = keyptssprings[j].index;

				const bool has_associated_wt = (keyptssprings[j].wt_idx > 0);
				const int cost_channel = has_associated_wt
					? keyptssprings[j].wt_idx - 1
					: params.no_wt_cost_index;

				if (cost_channel < 0 || cost_channel >= params.nb_wt) {
					throw std::runtime_error(
						"[skeleton] Invalid internal cost channel for spring " + std::to_string(j + 1) +
						". This means that a spring without water table was detected but the vadose-only cost channel was not initialized."
					);
				}

				
				auto skip_unreachable_connection = [&](const std::string& reason)
				{
					std::cout << "[skeleton][warn] No path could be computed for inlet "
						<< params.sinks_index[i]
						<< " toward outlet " << (j + 1)
						<< ". The connection will be ignored. "
						<< "Reason: " << reason
						<< ". inlet_order_index=" << (i + 1)
						<< ", inlet_connectivity_index=" << params.sinks_index[i]
						<< ", outlet_connectivity_index=" << (j + 1)
						<< ", source_node=" << source
						<< ", target_node=" << target
						<< ", source_xyz=(" << keyptssinks[i].p.x << ", " << keyptssinks[i].p.y << ", " << keyptssinks[i].p.z << ")"
						<< ", target_xyz=(" << keyptssprings[j].p.x << ", " << keyptssprings[j].p.y << ", " << keyptssprings[j].p.z << ")"
						<< ", wt_idx=" << keyptssprings[j].wt_idx
						<< ", cost_channel=" << cost_channel
						<< ", connectivity_flag=" << connectivity_flag
						<< std::endl;

					new_connectivity_matrix[params.sinks_index[i] - 1][j] = 0;
				};

				const clock_t time2 = clock();

				std::vector<int> path_full;
				std::vector<float> path_cost_full;
				std::vector<char> vadose_full;
				float pathSize_full = 0.0f;

				if (!has_associated_wt) {
					// Springs without an associated water table are handled as overflow/relict
					// outlets. They are connected directly from the inlet to the spring through
					// the vadose-only cost channel, and the whole path is flagged as vadose.
					std::vector<float> distances_direct;
					std::vector<int> previous_direct;
					float pathSize_direct = 0.0f;

										const clock_t time_shortest_path_begin = clock();

					switch (K_SHORTEST_PATH_METHOD) {
					case ShortestPathMethod::NormalDijkstra:
						DijkstraComputePaths(
							cost_channel,
							source,
							distances_direct,
							previous_direct,
							target
						);
						break;

					case ShortestPathMethod::BidirectionalDijkstra:
						DijkstraComputePathsBidirectional(
							cost_channel,
							source,
							target,
							distances_direct,
							previous_direct
						);
						break;
					}

					const clock_t time_shortest_path_end = clock();
					step_shortest_path += float(time_shortest_path_end - time_shortest_path_begin) / CLOCKS_PER_SEC;

					std::pair<std::vector<int>, std::vector<float>> pair_direct =
					DijkstraGetShortestPathTo(target, previous_direct, distances_direct, pathSize_direct);

					path_full = pair_direct.first;
					path_cost_full = pair_direct.second;
					vadose_full.resize(path_full.size(), true);
					pathSize_full = pathSize_direct;

					if (path_full.size() <= 1 || !std::isfinite(pathSize_full)) {
						skip_unreachable_connection("no direct Dijkstra path from inlet to outlet");
						continue;
					}
				}
				else {
					int reach = -1; // first reachable node located on the associated water table

				// 1) GENERATING VADOSE ZONE CONDUIT

					std::vector<float> distances_vadose;
					std::vector<int> previous_vadose;
					bool already_reached = false;
					float pathSize_vadose = 0.0f;

					const clock_t time_shortest_path_begin_vadose = clock();

					switch (K_SHORTEST_PATH_METHOD) {
					case ShortestPathMethod::NormalDijkstra:
					DijkstraComputePathsSurface(
						cost_channel,
						source,
						reach,
						distances_vadose,
						previous_vadose,
						samples_on_wt_flags,
						target,
						already_reached
					);
					break;

					case ShortestPathMethod::BidirectionalDijkstra:
						DijkstraComputePathsSurfaceBidirectional(
							cost_channel,
							source,
							reach,
							distances_vadose,
							previous_vadose,
							samples_on_wt_flags,
							target,
							already_reached
						);
						break;
					}

					const clock_t time_shortest_path_end_vadose = clock();
					step_shortest_path += float(time_shortest_path_end_vadose - time_shortest_path_begin_vadose) / CLOCKS_PER_SEC;

					if (reach < 0) {
						skip_unreachable_connection("no reachable water-table node from inlet");
						continue;
					}

					std::pair<std::vector<int>, std::vector<float>> pair_vadose =
						DijkstraGetShortestPathTo(reach, previous_vadose, distances_vadose, pathSize_vadose);

					std::vector<int> path_vadose = pair_vadose.first;
					std::vector<float> path_cost_vadose = pair_vadose.second;

					if (path_vadose.empty() || !std::isfinite(pathSize_vadose)) {
						skip_unreachable_connection("invalid vadose path from inlet to reached water-table node");
						continue;
					}

					if (source != reach && path_vadose.size() <= 1) {
						skip_unreachable_connection("degenerate vadose path from inlet to reached water-table node");
						continue;
					}

					// 2) GENERATING PHREATIC ZONE CONDUIT
					std::vector<float> distances;
					std::vector<int> previous;
					float pathSize = 0.0f;
					std::vector<int> path;
					std::vector<float> path_cost;

					if (!already_reached) {
						const clock_t time_shortest_path_begin_phreatic = clock();

						switch (K_SHORTEST_PATH_METHOD) {
						case ShortestPathMethod::NormalDijkstra:
							DijkstraComputePaths(
								cost_channel,
								reach,
								distances,
								previous,
								target
							);
							break;

						case ShortestPathMethod::BidirectionalDijkstra:
							DijkstraComputePathsBidirectional(
								cost_channel,
								reach,
								target,
								distances,
								previous
							);
							break;
						}

						const clock_t time_shortest_path_end_phreatic = clock();
						step_shortest_path += float(time_shortest_path_end_phreatic - time_shortest_path_begin_phreatic) / CLOCKS_PER_SEC;

						std::pair<std::vector<int>, std::vector<float>> pair =
							DijkstraGetShortestPathTo(target, previous, distances, pathSize);

						path = pair.first;
						path_cost = pair.second;

						if (path.size() <= 1 || !std::isfinite(pathSize)) {
							skip_unreachable_connection("no phreatic Dijkstra path from reached water-table node to outlet");
							continue;
						}
					}

					if (!path_vadose.empty()) {
						path_full.insert(path_full.end(), path_vadose.begin(), path_vadose.end());
						path_cost_full.insert(path_cost_full.end(), path_cost_vadose.begin(), path_cost_vadose.end());
						vadose_full.resize(path_vadose.size(), true);
					}

					pathSize_full = pathSize + pathSize_vadose;

					if (!already_reached && !path.empty()) {
						vadose_full.insert(vadose_full.end(), path.size() - 1, false);
						path_full.insert(path_full.end(), path.begin() + 1, path.end());
						path_cost_full.insert(path_cost_full.end(), path_cost.begin() + 1, path_cost.end());
					}
				}

				if (path_full.size() <= 1 || !std::isfinite(pathSize_full)) {
					skip_unreachable_connection("final inlet/outlet path is empty, degenerate, or non-finite");
					continue;
				}

				all_paths[i][j] = path_full;
				all_vadose[i][j] = vadose_full;
				all_distances[i][j] = pathSize_full;
				all_selection_distances[i][j] = pathSize_full;
				all_distances_detailled[i][j] = path_cost_full;
				clock_t time7 = clock();
				step6 += float(time7 - time2) / CLOCKS_PER_SEC;
			}
			clock_t time8 = clock();
			// Connectivity choice: select the outlet with the lowest corrected cumulative path cost.
			// The correction is applied only to ambiguous links (connectivity flag == 2) and only
			// for the selection step. Stored path costs remain unchanged.
			std::vector<int> ambiguous_candidates;
			std::vector<float> candidate_dz;

			ambiguous_candidates.reserve(all_paths[i].size());
			candidate_dz.reserve(all_paths[i].size());

			for (int j = 0; j < all_paths[i].size(); ++j) {
				const int connectivity_flag = params.connectivity_matrix[params.sinks_index[i] - 1][j];

				if (connectivity_flag != 2) continue;
				if (all_paths[i][j].size() <= 1) continue;
				if (!std::isfinite(all_distances[i][j])) continue;

				ambiguous_candidates.push_back(j);

				// Positive values mean that the outlet is below the inlet. Larger values
				// represent a larger elevation drop and should therefore be less penalized.
				const float dz = keyptssinks[i].p.z - keyptssprings[j].p.z;
				candidate_dz.push_back(dz);
			}

			if (!ambiguous_candidates.empty()) {
				const auto minmax_dz = std::minmax_element(candidate_dz.begin(), candidate_dz.end());
				const float min_dz = *minmax_dz.first;
				const float max_dz = *minmax_dz.second;
				const float dz_range = max_dz - min_dz;

				for (int k = 0; k < ambiguous_candidates.size(); ++k) {
					const int j = ambiguous_candidates[k];

					float gradient_penalty = 0.0f;

					if (params.gradient_constraint_weight > 0.0f && dz_range > eps) {
						// Penalty convention:
						//   dz == max_dz -> 0, most favorable outlet elevation drop
						//   dz == min_dz -> 1, least favorable outlet elevation drop
						//   intermediate dz values are linearly interpolated.
						gradient_penalty = (max_dz - candidate_dz[k]) / dz_range;
						gradient_penalty = KarstNSim::Clamp(gradient_penalty, 0.0f, 1.0f);
					}

					const float gradient_factor = 1.0f + params.gradient_constraint_weight * gradient_penalty;
					all_selection_distances[i][j] = all_distances[i][j] * gradient_factor;
				}
			}

			float shortest_distance = std::numeric_limits<float>::infinity();
			for (int j = 0; j < all_paths[i].size(); ++j) {
				const int connectivity_flag = params.connectivity_matrix[params.sinks_index[i] - 1][j];

				if (connectivity_flag != 2) continue;
				if (all_paths[i][j].size() <= 1) continue;

				shortest_distance = std::min(shortest_distance, all_selection_distances[i][j]);
			}

			for (int j = 0; j < all_paths[i].size(); ++j) {
				const int connectivity_flag = params.connectivity_matrix[params.sinks_index[i] - 1][j];

				if (all_paths[i][j].size() <= 1) continue;

				const float outlet_selection_threshold =
					std::isfinite(shortest_distance)
					? shortest_distance * (params.outlet_selection_cost_factor)
					: std::numeric_limits<float>::infinity();

				if (connectivity_flag == 2 &&
					all_selection_distances[i][j] > (outlet_selection_threshold + eps))
				{
					std::vector<int>().swap(all_paths[i][j]);
					std::vector<float>().swap(all_distances_detailled[i][j]);
					std::vector<char>().swap(all_vadose[i][j]);

					new_connectivity_matrix[params.sinks_index[i] - 1][j] = 0;
				}
				else
				{
					new_connectivity_matrix[params.sinks_index[i] - 1][j] = 1;
					retained_path_for_sink = true;

					// Update edge costs for subsequent inlet iterations.
					// Segments already karstified have their cost reduced, but only in the phreatic
					// zone unless vadose cohesion is enabled. No shortest-path preprocessing
					// invalidation is required here because bidirectional Dijkstra only caches the
					// reverse topology and reads edge weights directly from adj.
					for (int node_path = 0; node_path < all_paths[i][j].size() - 1; ++node_path)
					{
						const int u = all_paths[i][j][node_path];
						const int v = all_paths[i][j][node_path + 1];
						const int edge_slot = GetIdxNeighbor(u, v);

						if (adj[u][edge_slot].target != v) {
							throw std::runtime_error(
								"[skeleton] Internal error: could not recover an edge while applying cohesion."
							);
						}

						for (int cost_i = 0; cost_i < params.nb_wt; ++cost_i)
						{
							if (params.vadose_cohesion || samples_surf_flags(u, cost_i))
							{
								adj[u][edge_slot].weight[cost_i] *= fraction_karst_perm;
							}
						}
					}
				}
			}

			if (!retained_path_for_sink) {
				std::cout << "[skeleton][warn] No valid outlet could be retained for inlet "
					<< params.sinks_index[i]
					<< ". All candidate outlets were disconnected, located above the inlet, or unreachable."
					<< std::endl;
			}

			clock_t time9 = clock();
			step7 += float(time9 - time8) / CLOCKS_PER_SEC;
			step_sink = float(time9 - time_sink) / CLOCKS_PER_SEC;

			const int completed_sinks = i + 1;
			const int completed_percent = static_cast<int>(
				100.0f * static_cast<float>(completed_sinks) / static_cast<float>(std::max(1, total_sinks))
				);

			if (completed_percent >= next_progress_percent) {
				int logged_progress_percent = next_progress_percent;

				while (next_progress_percent <= 100 && completed_percent >= next_progress_percent) {
					logged_progress_percent = next_progress_percent;
					next_progress_percent += 10;
				}

				std::cout << "Skeleton progress: "
					<< completed_sinks << " / " << total_sinks
					<< " inlets processed ("
					<< logged_progress_percent << "%), last inlet time = "
					<< step_sink << " s"
					<< std::endl;
			}
		}
		// Path pruning based on an anisotropic empty region criterion

		clock_t time10 = clock();
		for (int i = 0; i < all_paths.size(); i++)
		{
			for (int j = 0; j < all_paths[i].size(); j++)
			{
				if (all_paths[i][j].size() <= 1)
					continue;
				// Second, check if there is a cheaper path from nodes[i] going through a node[k] to arrive at nodes[j]
				// float d_ij = KarstNSim::Pow(all_distances[i][j], params.gamma);
				bool keep = true;
				//for (int k = 0; k < all_paths[i].size(); k++)
				//{
				//	if (i == k || j == k)			 continue;
				//	if (all_paths[i][k].size() <= 1) continue;
				//	if (all_paths[k][j].size() <= 1) continue;

				//	float d_ik = KarstNSim::Pow(all_distances[i][k], params.gamma);
				//	float d_kj = KarstNSim::Pow(all_distances[k][j], params.gamma);
				//	if (d_ik + d_kj < d_ij)
				//	{
				//		//keep = false; // uncomment to enable gamma-graph constraint (use at your own risk)
				//		break;
				//	}
				//}
				if (keep) {
					pathsFinal.push_back(std::move(all_paths[i][j]));
					costsFinal.push_back(std::move(all_distances_detailled[i][j]));
					vadoseFinal.push_back(std::move(all_vadose[i][j]));
					springidxFinal.push_back(j);
				}
			}
		}

		if (save_new_connectivity_matrix) {
			std::string full_dir_name = params.directoryname + "/outputs";
			std::string full_name = params.scenename + "_connectivity_matrix.txt";
			save_connectivity_matrix(full_name, full_dir_name, new_connectivity_matrix);
		}

		clock_t time11 = clock();
		step8 += float(time11 - time10) / CLOCKS_PER_SEC;

		if (K_CLEAR_SHORTEST_PATH_PREPROCESSING_AFTER_SKELETON) {
			ClearShortestPathPreprocessing();
		}
	}

	std::pair< std::vector<int>, std::vector<float>> GraphOperations::shortest_path(int u, int v, int idx_specific_wt)
	{
		std::vector<float> distances;
		std::vector<int> previous;
		std::pair< std::vector<int>, std::vector<float>> pair;
		float pathSize = 0.;
		std::vector<int> path;
		std::vector<float> path_cost;
		DijkstraComputePaths(idx_specific_wt - 1, u, distances, previous); // Computes shortest paths from u to other vertices
		pair = DijkstraGetShortestPathTo(v, previous, distances, pathSize); // gets shortest path computed from v to u (back-propagated)
		path = pair.first;
		path_cost = pair.second;

		return std::make_pair(path, path_cost);
	}

	std::vector<std::vector<int>> GraphOperations::amplifykarstickseleton(const std::vector<KarsticNode>& baseSkeletonNodes,
		const std::vector<KeyPoint>& newkeypts)
	{
		// Our goal is to connect the new key points to the original network
		std::vector<std::vector<int>> pathsFinal;
		for (int i = 0; i < newkeypts.size(); i++)
		{
			std::vector<float> distances;
			std::vector<int> previous;

			float minDist = std::numeric_limits<float>::infinity();
			std::vector<int> bestPath;
			for (int j = 0; j < baseSkeletonNodes.size(); j++)
			{
				if (i == j)
					continue;
				int endIndex = baseSkeletonNodes[j].index;
				float totalDistance = 0.0;
				std::vector<int> path;
				std::vector<float> path_cost;
				std::pair< std::vector<int>, std::vector<float>> pair;
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

	std::vector<GraphOperations::EdgeWeightBackup>
		GraphOperations::add_ball_cost(
			const Vector3& center,
			float ball_radius
		) {
		std::vector<EdgeWeightBackup> reset_values;

		if (!std::isfinite(ball_radius) || ball_radius <= 0.0f) {
			return reset_values;
		}

		for (int source_index = 0;
			source_index < static_cast<int>(samples.size());
			++source_index) {

			const Vector3& point = samples[source_index];
			const float distance = magnitude(center - point);

			if (distance >= ball_radius) {
				continue;
			}

			for (int edge_slot = 0;
				edge_slot < static_cast<int>(adj[source_index].size());
				++edge_slot) {

				auto& edge = adj[source_index][edge_slot];

				// Ignore unused adjacency slots.
				if (edge.target < 0) {
					continue;
				}

				// Use the actual number of channels (=wt) stored by this edge
				for (int cost_channel = 0;
					cost_channel < static_cast<int>(edge.weight.size());
					++cost_channel) {

					reset_values.push_back({
						source_index,
						edge_slot,
						cost_channel,
						edge.weight[cost_channel]
						});

					const float relative_distance = distance / ball_radius;
					edge.weight[cost_channel] +=
						(1.0f - relative_distance) *
						2.0f *
						edge.weight[cost_channel];
				}
			}
		}

		return reset_values;
	}

	void GraphOperations::save_samples(const std::string& path) const
	{
		std::ofstream out;
		out.open(path);

		// ply header
		out << "ply\n";
		out << "format ascii 1.0\n";
		out << "element vertex " << samples.size() << "\n";
		out << "property float x\n";
		out << "property float y\n";
		out << "property float z\n";
		out << "end_header\n";

		// data
		for (auto p : samples)
			out << p.x << " " << p.y << " " << p.z << "\n";

		out.close();
	}

}