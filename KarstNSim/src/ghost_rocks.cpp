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
			if (pt.z < closest_point.z - length || pt.z > closest_point.z || GraphOperations::CheckBelowSurf(pt, substratum_surf, centers2D)) {
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

	void paint_KP_with_ghostrocks(
		const Box& grid,
		std::vector<float>& ikp,
		float length,
		float width,
		const Line& polyline,
		const bool& use_max_depth_constraint,
		const Surface& substratum_surf,
		float ghost_rock_weight)
	{
		const int nu = grid.get_nu();
		const int nv = grid.get_nv();
		const int nw = grid.get_nw();

		const int nb_segments = polyline.get_nb_segs();

		// Keep the existing final normalization behavior even if the alteration
		// polyline unexpectedly contains no segment.
		if (nb_segments <= 0) {
			standardize_to_range(ikp, 0.0f, 1.0f);
			return;
		}

		// Build the substratum spatial index once. It is only queried for cells
		// that have already passed all cheaper ghost-rock geometric tests.
		const PointCloud centers2D = substratum_surf.get_centers_cloud(2);

		const Vector2 bbox_min = polyline.get_bbox_min();
		const Vector2 bbox_max = polyline.get_bbox_max();

		// The plan-view rejection must include the maximum lateral radius of the
		// ghost-rock corridor.
		const float bbox_margin = 0.5f * width;

		const float bbox_x_min = bbox_min.x - bbox_margin;
		const float bbox_x_max = bbox_max.x + bbox_margin;
		const float bbox_y_min = bbox_min.y - bbox_margin;
		const float bbox_y_max = bbox_max.y + bbox_margin;

		// Compute global elevation bounds of the alteration lines once. Any voxel
		// outside [minimum alteration elevation - length, maximum alteration
		// elevation] cannot belong to any ghost-rock corridor and can therefore be
		// rejected before searching for its closest alteration-line segment.
		float alteration_min_z = std::numeric_limits<float>::max();
		float alteration_max_z = std::numeric_limits<float>::lowest();

		for (int segment_index = 0; segment_index < nb_segments; ++segment_index) {
			Segment segment = polyline.get_seg(segment_index);
			const Vector3 p1 = segment.start();
			const Vector3 p2 = segment.end();

			alteration_min_z = std::min(
				alteration_min_z,
				std::min(p1.z, p2.z)
			);

			alteration_max_z = std::max(
				alteration_max_z,
				std::max(p1.z, p2.z)
			);
		}

		const float ghostrock_global_min_z = alteration_min_z - length;
		const float ghostrock_global_max_z = alteration_max_z;

		// When a substratum constraint is used, its map-view bounding box can be
		// used to avoid unnecessary KD-tree queries. Outside this bounding box,
		// CheckBelowSurf() would not find a containing triangle and would return
		// false, so skipping the query preserves the existing behavior.
		Vector3 substratum_bbox_min;
		Vector3 substratum_bbox_max;

		if (use_max_depth_constraint) {
			substratum_bbox_min = substratum_surf.get_boundbox_min();
			substratum_bbox_max = substratum_surf.get_boundbox_max();
		}

		// The IKP array is flattened as:
		// idx = u + nu * (v + nv * w).
		//
		// Keeping u as the innermost loop therefore traverses IKP contiguously in
		// memory and improves cache locality.
		for (int w = 0; w < nw; ++w) {
			for (int v = 0; v < nv; ++v) {
				for (int u = 0; u < nu; ++u) {

					const Vector3 pt = grid.uvw2xyz(u, v, w);

					// Cheap plan-view rejection before any segment-distance
					// calculation.
					if (pt.x < bbox_x_min ||
						pt.x > bbox_x_max ||
						pt.y < bbox_y_min ||
						pt.y > bbox_y_max) {
						continue;
					}

					// Cheap global vertical rejection before searching all alteration
					// segments.
					if (pt.z < ghostrock_global_min_z ||
						pt.z > ghostrock_global_max_z) {
						continue;
					}

					// Find the closest point on the closest alteration-line segment.
					// This logic is kept locally here so that the subsequent tests can
					// be reordered from cheapest to most expensive without adding a
					// new helper function.
					Vector3 closest_point;
					float min_distance_sq = std::numeric_limits<float>::max();

					for (int segment_index = 0;
						segment_index < nb_segments;
						++segment_index) {

						Segment segment =
							polyline.get_seg(segment_index);

						const Vector3 p1 = segment.start();
						const Vector3 p2 = segment.end();

						float distance_sq = 0.0f;
						Vector3 closest_point_on_segment;

						squaredistance_to_segment2D(
							pt,
							p1,
							p2,
							distance_sq,
							closest_point_on_segment
						);

						if (distance_sq < min_distance_sq) {
							min_distance_sq = distance_sq;
							closest_point = closest_point_on_segment;
						}
					}

					// The voxel must lie vertically between the alteration line and
					// the maximum ghost-rock depth.
					if (pt.z < closest_point.z - length ||
						pt.z > closest_point.z) {
						continue;
					}

					// Compute the local lateral radius of the ghost-rock corridor.
					// Keeping ellipsis_width() here preserves the existing corridor
					// geometry exactly.
					const float width_z = ellipsis_width(
						pt.z,
						closest_point.z - length / 2.0f,
						width,
						length,
						2
					);

					// Reject laterally distant voxels before querying the substratum
					// surface. This avoids the comparatively expensive KD-tree search
					// performed by CheckBelowSurf() for voxels that are not actually
					// inside the ghost-rock corridor.
					if (min_distance_sq >= width_z * width_z) {
						continue;
					}

					if (use_max_depth_constraint) {

						// Outside the substratum map-view bounding box,
						// CheckBelowSurf() would return false because no containing
						// triangle can exist there.
						const bool inside_substratum_bbox =
							pt.x >= substratum_bbox_min.x &&
							pt.x <= substratum_bbox_max.x &&
							pt.y >= substratum_bbox_min.y &&
							pt.y <= substratum_bbox_max.y;

						if (inside_substratum_bbox &&
							GraphOperations::CheckBelowSurf(
								pt,
								substratum_surf,
								centers2D
							)) {
							continue;
						}
					}

					const std::size_t idx =
						static_cast<std::size_t>(u) +
						static_cast<std::size_t>(nu) *
						(
							static_cast<std::size_t>(v) +
							static_cast<std::size_t>(nv) *
							static_cast<std::size_t>(w)
							);

					if (ikp[idx] > -10000.0f) {
						ikp[idx] += ghost_rock_weight;
					}
				}
			}
		}

		standardize_to_range(ikp, 0.0f, 1.0f);
	}
}