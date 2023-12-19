/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#pragma once

#include <string>
#include <vector>
#include "KarstNSim/vec.h"  // header file for Vector3
#include "KarstNSim/basics.h" // header file for Box, Surface, Line, etc.
#include "KarstNSim/run_code.h"
#include "KarstNSim/read_files.h"

class ParseInputs {
public:
	ParseInputs();  // Constructor declaration, if needed
	bool parseBoolean(const std::string& input);
	KarstNSim::ParamsSource parse(const std::string& filename);
	// Other member function declarations, if needed
};