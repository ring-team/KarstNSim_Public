/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/ghost_rocks.h"

namespace KarstNSim {

	float ellipsis_width(float z, float z0, float dw, float dz, int power) {
		// Calculate w(z) using the modified ellipsis formula where ellipsis is of semi-axes dw and dz, and centered at z=z0
		// Note that power = 1 results in a normal ellipsis. Increasing power makes w evolution closer to that of a rectangular section (ie. w = cst = dw/2 everywhere except for z<=z0-dz/2 and z>=z0+dz/2 where w=0)
		float w = dw / 2 * std::sqrt(1 - std::pow((z - z0) / (dz / 2), 2 * power));
		return w;
	}

	// NB : this is done in 2D because the z distance doesnt change anything
	void closest_point_on_closest_segment_to_point_in_polyline(const Vector3& pt, const Line& polyline, Vector3& closest_point, float& min_distance_sq) {
		// Initialize variables to keep track of the closest distance and segment
		min_distance_sq = std::numeric_limits<float>::max();
		float distance_sq;
		Vector3 closest_pt_seg;
		Vector3 p1;
		Vector3 p2;
		// Iterate through each segment in the polyline
		for (int i = 0; i < polyline.get_nb_segs(); ++i) {
			// Get the start and end points of the segment
			p1 = polyline.get_seg(i).start();
			p2 = polyline.get_seg(i).end();

			// Calculate the squared distance from the point to the segment
			squaredistance_to_segment2D(pt, p1, p2, distance_sq, closest_pt_seg);
			// Update the closest segment if this segment is closer
			if (distance_sq < min_distance_sq) {
				min_distance_sq = distance_sq;
				closest_point = closest_pt_seg;
			}
		}
	}

	bool is_pt_in_ghostrock(const Vector3& pt, float length, float width, const Line& polyline, const bool& use_max_depth_constraint, const Surface& substratum_surf, int power, const PointCloud& centers2D, float& width_z) {

		Vector3 closest_point;
		float min_distance_sq;

		// find the closest point on the closest segment in the alteration polyline to the cell
		closest_point_on_closest_segment_to_point_in_polyline(pt, polyline, closest_point, min_distance_sq);
		// check that the z of the cell pt is within boundaries of the ghost-rock (and, if used, above the substratum horizon)
		if (use_max_depth_constraint) {
			if (pt.z < closest_point.z - length || pt.z > closest_point.z || GraphOperations::CheckBelowSurf(pt, &substratum_surf, centers2D)) {
				return false;
			}
		}
		else {
			if (pt.z < closest_point.z - length || pt.z > closest_point.z) {
				return false;
			}
		}

		// check that the 2D distance between grid cell and ghost rock point is smaller than width
		width_z = ellipsis_width(pt.z, closest_point.z - length / 2, width, length, power);
		if (min_distance_sq < width_z*width_z) {
			return true;
		}
		return false;
	}

	void paint_karst_sections_with_ghostrocks(KarsticSkeleton& skel, float length, float width, const Line& polyline, const bool& use_max_depth_constraint, const Surface& substratum_surf) {

		PointCloud centers2D = substratum_surf.get_centers_cloud(2);
		float width_z;
		int power = 2;
		bool is_inside;
		int compt = 0;

		int nb_nodes = int(skel.nodes.size());

		for (int i = 0; i < nb_nodes; ++i) {
			// Calculate the coordinates of the voxel center
			Vector3 pt = skel.nodes[i].p;
			is_inside = is_pt_in_ghostrock(pt, length, width, polyline, use_max_depth_constraint, substratum_surf, power, centers2D, width_z);
			if (!is_inside) {
				continue;
			}
			compt++;
			// if we get there, it means that the node of the karst skeleton IS inside a ghost rock, and we will therefore increase its associated section by the width of the ghostrock corridor
			skel.nodes[i].eq_radius = width_z;
		}
	}

	void paint_KP_with_ghostrocks(const Box& grid, std::vector<float>& ikp, float length, float width, const Line& polyline, const bool& use_max_depth_constraint, const Surface& substratum_surf, int ghost_rock_weight) {

		// Iterate through each voxel in the grid
		PointCloud centers2D = substratum_surf.get_centers_cloud(2);
		float width_z;
		int power = 2;
		int idx;
		bool is_inside;
		int compt = 0;

		for (int u = 0; u < grid.get_nu(); ++u) {
			for (int v = 0; v < grid.get_nv(); ++v) {
				for (int w = 0; w < grid.get_nw(); ++w) {
					// Calculate the coordinates of the voxel center
					Vector3 pt = grid.uvw2xyz(u, v, w);

					// Bounding Box Check: Quick reject cells outside potential ghost rock influence
					Vector2 bbox_min = polyline.get_bbox_min();
					Vector2 bbox_max = polyline.get_bbox_max();
					if (pt.x < bbox_min.x || pt.x > bbox_max.x || pt.y < bbox_min.y || pt.y > bbox_max.y) {
						continue;
					}

					is_inside = is_pt_in_ghostrock(pt, length, width, polyline, use_max_depth_constraint, substratum_surf, power, centers2D, width_z);
					if (!is_inside) {
						continue;
					}
					grid.ravel(u, v, w, idx);
					if (ikp[idx] > -10000) {
						ikp[idx] += ghost_rock_weight; // IKP increased when cell intersected by ghost rock
						compt++;
					}
				}
			}
		}
		// Standardize the values in ikp between 0 and 1
		standardize_to_range(ikp, 0.0, 1.0);
	}
}