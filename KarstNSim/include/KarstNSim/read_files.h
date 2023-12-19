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
#include <string>
#include <random>
namespace KarstNSim {
	void load_point(const std::string& file_name, const std::string& save_directory, Vector3& u, std::vector<double>& 
	);
	void load_surface(const std::string& file_name, const std::string& save_directory, Surface& s, std::vector<std::vector<double>>& properties);
	void load_pointset(const std::string& file_name, const std::string& save_directory, std::vector<Vector3>& pset, std::vector<std::vector<double>>& properties);
	void load_line(const std::string& file_name, const std::string& save_directory, Line& pline, std::vector<std::vector<std::vector<double>>>& properties);
	void load_box(const std::string& file_name, const std::string& save_directory, Box& box, std::vector<std::vector<double>>& properties);
}
