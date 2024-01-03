/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/parse_inputs.h"
#include "KarstNSim/run_code.h"
#include <chrono>
#include <conio.h>

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <instruction_file>" << std::endl;
		return 1;
	}

	std::string instructionFile = argv[1];

	std::cout << "Parsing instruction file: " << instructionFile << std::endl;

	ParseInputs inputParser;
	KarstNSim::ParamsSource params = inputParser.parse(instructionFile);

	std::cout << "Parsing completed. Running simulation..." << std::endl;

	// Measure the start time
	auto startTime = std::chrono::high_resolution_clock::now();

	// Now 'params' contains the parsed parameters
	KarstNSim::run_simulation_full(params);

	// Measure the end time
	auto endTime = std::chrono::high_resolution_clock::now();

	// Calculate the elapsed time in seconds
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();

	std::cout << "Simulation completed successfully. Elapsed time: " << duration << " seconds." << std::endl;

    std::cout << "\nPress any key to exit..." << std::endl;
    _getch();  // Wait for any key without Enter
    return 0;
}