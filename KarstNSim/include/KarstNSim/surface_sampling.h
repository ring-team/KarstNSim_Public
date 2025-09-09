/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#pragma once

/**
@file surface_sampling.h
@brief surface_sampling class used to sample points onto triangulated surfaces vertices, and refine sampling if needed.
@author Augustin GOUY
**/


#include "KarstNSim/write_files.h"
#include "KarstNSim/graph.h"
#include <unordered_set>
#include <set>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <fstream>

namespace KarstNSim {

	/*!
	\class surface_sampling
	\brief Class that contains method allowing triangulated surfaces sampling and refinement
	*/
	class surface_sampling {
	public:
		using pointList = std::vector<Vector3>;
		/*!
		\brief Samples points on several surfaces, with or without surface refinement and adds those points to the point cloud.
		\param directoryname path of the project
		\param network_name name of the karst network
		\param Box box containing the geological model and surfaces
		\param Surface List of surfaces to sample
		\param Points Sampling cloud which will be updated by the function
		\param refine_surface_sampling option to use if surface vertices should be refined.
		A value of 0 means no refinement (only 1 node sampled per vertex).
		A value of 1 means that a new point is sampled at the middle of each edge of each triangle of the surface, effectively forming new triangles.
		A value of n means that this fractal pattern operation is operated n times.
		\param create_surface_sampling option which will create an ASCII file with only the points sampled on the surfaces if set to true.
		*/
		static void multiple_surface_sampling(std::string directoryname, std::string network_name, KarstNSim::Box* Box, std::vector<Surface>* Surface, std::vector<Vector3>& Points, const int refine_surface_sampling, const bool create_surface_sampling);

	};
}