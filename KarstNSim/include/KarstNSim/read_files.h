/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#pragma once

/**
@file read_files.h
@brief Functions used to create standard C++ library objects from ASCII files.
@author Augustin GOUY
**/

#include <string>
#include <cmath>
#include <iostream>
#include "KarstNSim/basics.h"
#include "KarstNSim/vec.h"
#include <sstream>
#include <fstream>
#include <map>
#include <algorithm>
#include <random>
#include <vector>

namespace KarstNSim {

	/**
	 * @brief Loads a single point and its associated properties from an ASCII file.
	 *
	 * @param file_name The path to the input ASCII file.
	 * @param save_directory The directory where input ASCII file should be searched.
	 * @param u A Vector3 object where the point's coordinates will be stored.
	 * @param properties A vector where the associated properties of the point will be stored.
	 */
	void load_point(const std::string& file_name, const std::string& save_directory, Vector3& u, std::vector<float>& properties);

	/**
	 * @brief Loads a surface and its associated properties from an ASCII file.
	 *
	 * @param file_name The path to the input ASCII file.
	 * @param save_directory The directory where input ASCII file should be searched.
	 * @param s A Surface object where the surface data will be stored.
	 * @param properties A 2D vector to store the associated properties of the surface.
	 */
	void load_surface(const std::string& file_name, const std::string& save_directory, Surface& s, std::vector<std::vector<float>>& properties);

	/**
	 * @brief Loads a point set and its associated properties from an ASCII file.
	 *
	 * @param file_name The path to the input ASCII file.
	 * @param save_directory The directory where input ASCII file should be searched.
	 * @param pset A vector of Vector3 objects where the point set will be stored.
	 * @param properties A 2D vector to store the associated properties of the point set.
	 */
	void load_pointset(const std::string& file_name, const std::string& save_directory, std::vector<Vector3>& pset, std::vector<std::vector<float>>& properties);

	/**
	 * @brief Loads a line and its associated properties from an ASCII file.
	 *
	 * @param file_name The path to the input ASCII file.
	 * @param save_directory The directory where input ASCII file should be searched.
	 * @param pline A Line object where the line data will be stored.
	 * @param properties A 3D vector to store the associated properties of the line.
	 */
	void load_line(const std::string& file_name, const std::string& save_directory, Line& pline, std::vector<std::vector<std::vector<float>>>& properties);

	/**
	 * @brief Loads a box and its associated properties from an ASCII file.
	 *
	 * @param file_name The path to the input ASCII file.
	 * @param save_directory The directory where input ASCII file should be searched.
	 * @param box A Box object where the box data will be stored.
	 * @param properties A 2D vector to store the associated properties of the box.
	 */
	void load_box(const std::string& file_name, const std::string& save_directory, Box& box, std::vector<std::vector<float>>& properties);
}
