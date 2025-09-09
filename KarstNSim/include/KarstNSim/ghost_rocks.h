/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#pragma once

/**
@file ghost_rocks.h
@brief Contains methods used to paint ghost-rocks onto the background grid's intrinsic karstification potential property
@author Augustin GOUY
**/

#include <limits>       // for numeric_limits
#include <set>          // for set
#include <utility>      // for pair
#include <algorithm> // for std::sort
#include <fstream>	 // for saving sample file
#include <math.h>
#include <string>		// MSVC2017 requires this
#include <map>
#include <random>
#include <queue>
#include <unordered_map>
#include <stdlib.h>
#include <KarstNSim/nanoflann.hpp>
#include "KarstNSim/basics.h"
#include "KarstNSim/vec.h"
#include "KarstNSim/graph.h"


namespace KarstNSim {

	/**
 * @brief Computes the width of an ellipsis at a given z-coordinate.
 * @param z The z-coordinate.
 * @param z0 The center z-coordinate of the ellipsis.
 * @param dw The width of the ellipsis.
 * @param dz The height of the ellipsis.
 * @param power The power controlling the shape of the ellipsis.
 * @return The width of the ellipsis at z.
 */
	float ellipsis_width(float z, float z0, float dw, float dz, int power);

	/**
 * @brief Finds the closest point on the closest segment in a polyline to a given point.
 * @param pt The point for which the closest segment is to be found.
 * @param polyline The polyline in which to search for the closest segment.
 * @param closest_point The closest point on the segment (output).
 * @param min_distance_sq The squared distance to the closest point (output).
 */
	void closest_point_on_closest_segment_to_point_in_polyline(const Vector3& pt, const Line& polyline, Vector3& closest_point, float& min_distance_sq);

	/**
 * @brief Determines if a point lies within a ghost rock defined by a polyline.
 * @param pt The point to check.
 * @param length The length of the ghost rock.
 * @param width The width of the ghost rock.
 * @param polyline The polyline defining the ghost rock.
 * @param use_max_depth_constraint Whether to enforce a maximum depth constraint.
 * @param substratum_surf The substratum surface for depth constraint checks.
 * @param power The power controlling the ellipsis shape.
 * @param centers2D The 2D centers of the substratum surface.
 * @param width_z The width of the ghost rock at the z-coordinate of the point (output).
 * @return True if the point is inside the ghost rock; false otherwise.
 */
	bool is_pt_in_ghostrock(const Vector3& pt, float length, float width, const Line& polyline, const bool& use_max_depth_constraint, const Surface& substratum_surf, int power, const PointCloud& centers2D, float& width_z);

	/**
 * @brief Paints karstic sections with ghost rocks based on their proximity to a polyline.
 * @param skel The karstic skeleton to paint.
 * @param length The length of the ghost rock.
 * @param width The width of the ghost rock.
 * @param polyline The polyline defining the ghost rock.
 * @param use_max_depth_constraint Whether to enforce a maximum depth constraint.
 * @param substratum_surf The substratum surface for depth constraint checks.
 */
	void paint_karst_sections_with_ghostrocks(KarsticSkeleton& skel, float length, float width, const Line& polyline, const bool& use_max_depth_constraint, const Surface& substratum_surf);


	/**
	 * @brief Paints intrinsic karstification potential (IKP) on a grid using ghost rocks.
	 * @param grid The grid to paint.
	 * @param ikp The IKP property values for the grid (output).
	 * @param length The length of the ghost rock.
	 * @param width The width of the ghost rock.
	 * @param polyline The polyline defining the ghost rock.
	 * @param use_max_depth_constraint Whether to enforce a maximum depth constraint.
	 * @param substratum_surf The substratum surface for depth constraint checks.
	 * @param ghost_rock_weight The weight to apply to IKP values within ghost rocks.
	 */
	void paint_KP_with_ghostrocks(const Box& grid, std::vector<float>& ikp, float length, float width, const Line& polyline, const bool& use_max_depth_constraint, const Surface& substratum_surf, int ghost_rock_weight);

}
