/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#pragma once

/**
* @file parse_inputs.h
* @brief Header file for the ParseInputs class, which provides methods to parse input files for the KarstNSim simulation.
*
* This file declares the ParseInputs class, including methods for parsing boolean values
* and reading simulation parameters from an input file.
*/

#include <string>
#include <vector>
#include "KarstNSim/vec.h"  // header file for Vector3
#include "KarstNSim/basics.h" // header file for Box, Surface, Line, etc.
#include "KarstNSim/run_code.h"
#include "KarstNSim/read_files.h"

/**
* @class ParseInputs
* @brief A utility class for parsing input files and managing simulation parameters.
*
* The ParseInputs class contains methods for reading and interpreting simulation parameters
* from a structured input file. It also includes utility functions for parsing boolean values.
*/
class ParseInputs {
public:

	/**
	 * @brief Constructor for the ParseInputs class.
	 *
	 * Initializes a new instance of the ParseInputs class. Currently, it does not perform
	 * any specific initialization tasks.
	 */
	ParseInputs();  // Constructor declaration

	/**
	 * @brief Parses a boolean value from a string.
	 *
	 * This method interprets a string as a boolean value. The accepted values are:
	 * - `"true"` or `"1"` for true.
	 * - `"false"` or `"0"` for false.
	 *
	 * If the input string is invalid, the method returns `false` and prints an error message.
	 *
	 * @param input A string representation of a boolean value.
	 * @return A boolean value corresponding to the input string.
	 */
	bool parse_boolean(const std::string& input);

	/**
	 * @brief Reads a file containing a discrete distribution (separator = space) and stores it in a vector.
	 * @param directory The directory containing the file.
	 * @param filename The name of the file to read.
	 * @return The vector of floats with the values of the distribution.
	 */
	std::vector<float> read_distrib(const std::string& directory, const std::string& filename);

	/**
	 * @brief Parses simulation parameters from an input file.
	 *
	 * This method reads and parses a structured input file containing parameters for the
	 * KarstNSim simulation. The file may include parameters such as the domain, seed values,
	 * physical properties, and constraints. It supports various data types, including boolean
	 * flags, numerical values, and paths to external files.
	 *
	 * @param filename The path to the input file containing simulation parameters.
	 * @return A ParamsSource object populated with the parsed parameters.
	 */
	KarstNSim::ParamsSource parse(const std::string& filename);

};