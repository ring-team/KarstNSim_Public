/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/


#pragma once

/**
@file write_files.h
@brief Functions used to create ASCII files storing output objects.
@author Augustin GOUY
**/


#include <string>
#include <cmath>
#include <iostream>
#include "KarstNSim/basics.h"
#include "KarstNSim/vec.h"
#include <fstream>
#include <map>
#include <algorithm>
#include <string>
#include <random>

namespace KarstNSim {
	/*!
	 \brief Saves point in ASCII file
	
	 This function saves a given point, represented by its vector in 3D, in an ASCII file
	 with its corresponding properties.
	 
	 \param file_name Name of the created file.
	 \param save_directory Location of the created file.
	 \param u Vector3 representing the coordinates of the point to be saved.
	 \param property_names Vector of string with names of the properties of points.
	 \param properties Vector of double with values of the point for each property.
	!*/
	void save_point(const std::string& file_name, const std::string& save_directory, Vector3 u, std::vector<std::string> property_names = {}, std::vector<double> properties = {});
	/**
	* @brief Saves triangulated surface in ASCII file
	*
	* This function saves a given surface, represented by its individual vertices and
	* its triangles. The output file will first display the vertices coordinates and property values,
	* under the prefix 'VRTX', and then display the indices of the vertices of each triangle, under
	* the prefix 'TRGL'.
	*
	* @param file_name Name of the created file.
	* @param save_directory Location of the created file.
	* @param s Surface to be saved.
	* @param property_names Vector of string with names of the properties of the surface.
	* @param properties Vector of vector of double with values of the nodes of the surface for each property. Format : properties[i][j] is the property j of node i.
	**/
	void save_surface(const std::string& file_name, const std::string& save_directory, Surface s, std::vector<std::string> property_names = {}, std::vector<std::vector<double>> properties = {});
	/**
	* @brief Saves pointset in ASCII file
	*
	* This function saves a given pointset, represented by its nodes each with its
	* coordinated and its properties if any.
	*
	* @param file_name Name of the created file.
	* @param save_directory Location of the created file.
	* @param pset Pointset to be saved.
	* @param property_names Vector of string with names of the properties of the pointset.
	* @param properties Vector of vector of double with values of the nodes of the pointset for each property. Format : properties[i][j] is the property j of node i.
	*/
	void save_pointset(const std::string& file_name, const std::string& save_directory, std::vector<Vector3> pset, std::vector<std::string> property_names = {}, std::vector<std::vector<double>> properties = {});
	/**
	* @brief Saves line in ASCII file
	*
	* This function saves a given line, represented by its edges, each with two points (starting and
	* ending point of the edge) with coordinates and properties if any.
	*
	* @param file_name Name of the created file.
	* @param save_directory Location of the created file.
	* @param pline Line to be saved.
	* @param property_names Vector of string with names of the properties of the line.
	* @param properties Vector of vector of double with values of the nodes of the surface for each property. Format : properties[i][j] is the property j of node i.
	*/
	void save_line(const std::string& file_name, const std::string& save_directory, Line pline, std::vector<std::string> property_names = {}, std::vector<std::vector<std::vector<double>>> properties = {});
	/**
	* @brief Saves triangulated surface in ASCII file
	*
	* This function saves a given surface, represented by its individual vertices and
	* its triangles. The output file will first display the vertices coordinates and property values,
	* under the prefix 'VRTX', and then display the indices of the vertices of each triangle, under
	* the prefix 'TRGL'.
	*
	* @param file_name Name of the created file.
	* @param save_directory Location of the created file.
	* @param s Surface to be saved
	* @param property_names Vector of string with names of the properties of the surface
	* @param properties Vector of vector of double with values of the nodes of the surface for each property. Format : properties[i][j] is the property j of node i.
	*/
	void save_box(const std::string& file_name, const std::string& save_directory, Box box, std::vector<std::string> property_names = {}, std::vector<std::vector<double>> properties = {});
}
