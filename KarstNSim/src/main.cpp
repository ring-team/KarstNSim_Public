/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

/**
* @file main.cpp
* @brief Main entry point for KarstNSim. Parses input instructions and runs the simulation.
*
* This file handles the parsing of the instruction file and executes the main simulation
* using the KarstNSim framework. The elapsed time for the simulation is also measured.
*/

#include "KarstNSim/parse_inputs.h"
#include "KarstNSim/run_code.h"
#include <chrono>

/**
* @brief Entry point for the KarstNSim application.
*
* This function takes an optional command-line argument specifying the instruction file.
* If no argument is provided, it defaults to a predefined file location. It parses the input file,
* runs the simulation, and outputs the elapsed time.
*
* @param argc The number of command-line arguments.
* @param argv The array of command-line arguments.
* @return int Returns 0 if the program runs successfully; otherwise, returns 1.
*
* @note If no command-line argument is provided, the program defaults to "../../../Input_files/instructions.txt".
* @note The simulation parameters are parsed using the `ParseInputs` class, and the simulation is run using the `run_simulation_full` function.
*/
int main(int argc, char* argv[]) {
	//if (argc != 2) {
	//	std::cerr << "Usage: " << argv[0] << " <instruction_file>" << std::endl;
	//	return 1;
	//}
	std::string instructionFile;
	
    // Check if a command-line argument is provided
    if (argc > 1) {
        // Use the provided file name
        instructionFile = argv[1];
    } else {
        // Use a default file name
        instructionFile = "../../../Input_files/instructions.txt";
    }

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


    std::cout << "\nPress Enter to exit...";
    std::cin.get();  // Wait for Enter key
	
    return 0;
}