/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/surface_sampling.h"

// Allow use for unordered_set of the hash operator (which lets it know when there are duplicates)

namespace {
	struct myhash {
		size_t operator()(const Vector3 &p) const {
			return std::hash<float>{}(p.x) ^ std::hash<float>{}(p.x) ^ std::hash<float>{}(p.x);
		}
	};
}

using MySet = std::unordered_set<Vector3, myhash>;

namespace KarstNSim {
	void surface_sampling::multiple_surface_sampling(std::string directoryname, std::string network_name, Box* box, std::vector<Surface>* surface, std::vector<Vector3>& points, const int refine_surface_sampling, const bool create_vset_sampling)
	{
		// we iterate on the surfaces
		MySet res; // unordered sets don't allow for duplicates, which is perfect for our usage here
		for (int k = 0; k < surface->size(); k++) {
			// we iterate on the triangles of the surface
			Surface surfk = surface->at(k);
			for (int j = 0; j < surfk.get_nb_trgls(); j++) {
				bool allow_refinement = false;
				std::vector<Vector3> trgl_pts;
				Triangle trgl = surfk.get_triangle(j);
				// we iterate on the nodes of the triangles
				for (int i = 0; i < 3; i++) {
					// we get Vector3 and put it in the points list if it is not already in this list
					Vector3 p = surfk.get_node(trgl.point(i));
					trgl_pts.push_back(p);
					// continue only if the point is in the grid
					if (box->contains(p) && !allow_refinement) {
						allow_refinement = true;
					}
				}

				if (allow_refinement && refine_surface_sampling != -1) {
					for (int i = 0; i < 3; i++) {
						if (box->contains(trgl_pts[i])) {
							res.insert(trgl_pts[i]);
						}
					}
				}

				// refinement :
				// At each refinement iteration, we consider all non-overlapping possible triplets of points inside a triangle, which describe a sub-triangle.
				// The level of refinement is considered at :
				// -> the side points : a point is created at the middle of each segment of all the previously created sub-triangles
				// This simple algorithm allows to very quickly densify points inside the triangles uniformly, using a triforce pattern

				if (allow_refinement && refine_surface_sampling > 0) { // if triangle is at least partly inside the Box, try to refine it

					std::vector<std::vector<Vector3>> prev_tri;
					prev_tri.push_back(trgl_pts);

					for (int count = 0; count < refine_surface_sampling; count++) {

						std::vector<std::vector<Vector3>> new_tri;
						int prev_tri_size = int(prev_tri.size());
						for (int pt = 0; pt < prev_tri_size; pt++) {
							Vector3 a = prev_tri[pt][0];
							Vector3 b = prev_tri[pt][1];
							Vector3 c = prev_tri[pt][2];
							Vector3 u = (a + b) / 2;
							Vector3 v = (c + b) / 2;
							Vector3 w = (a + c) / 2;
							if (box->contains(u)) {
								res.insert(u);
							}
							if (box->contains(v)) {
								res.insert(v);
							}
							if (box->contains(w)) {
								res.insert(w);
							}
							new_tri.push_back(std::vector<Vector3>({ a,u,w }));
							new_tri.push_back(std::vector<Vector3>({ b,u,v }));
							new_tri.push_back(std::vector<Vector3>({ c,w,v }));
							new_tri.push_back(std::vector<Vector3>({ u,v,w }));
						}

						prev_tri = new_tri;
					}
				}
			}
		}

		MySet::iterator itr_set;
		for (itr_set = res.begin(); itr_set != res.end(); ++itr_set) {
			points.push_back(*itr_set);
		}

		if (create_vset_sampling) {
			std::string name = network_name + "_s_pts.txt";
			std::string full_dir_name = directoryname + "/outputs";
			save_pointset(name, full_dir_name, points);
		}
	}
}