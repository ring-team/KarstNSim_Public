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


// --- Standard headers ---
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>

// --- Filesystem detection (robust across MSVC & standards) ---
// Silence MSVC's hard error on <experimental/filesystem>
#if defined(_MSC_VER) && !defined(_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#endif

// Determine active language level (MSVC vs others)
#if defined(_MSC_VER)
#ifndef KNS_LANG
#define KNS_LANG _MSVC_LANG
#endif
#else
#ifndef KNS_LANG
#define KNS_LANG __cplusplus
#endif
#endif

#if __has_include(<filesystem>) && (KNS_LANG >= 201703L)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#error "No <filesystem> support found. Compile with C++17 or ensure <experimental/filesystem> exists."
#endif

// --- Your project headers ---
#include "KarstNSim/parse_inputs.h"
#include "KarstNSim/run_code.h"
#include <chrono>


/**
* @brief Discover example simulations under the given root path.
*
* Only subfolders that contain a file named exactly "instructions.txt" are considered.
* The function returns a list of pairs: (display_name, absolute_path_to_instruction_file).
*/
static std::vector<std::pair<std::string, fs::path>>
discover_examples(const fs::path& root)
{
std::vector<std::pair<std::string, fs::path>> out;

if (!fs::exists(root) || !fs::is_directory(root)) {
    return out;
}

for (const auto& entry : fs::directory_iterator(root)) {
    // Use free functions to stay compatible with std::experimental::filesystem
    if (!fs::is_directory(entry.path())) continue;

    const fs::path subdir = entry.path();
    const fs::path preferred = subdir / "instructions.txt";

    if (fs::exists(preferred) && fs::is_regular_file(preferred)) {
        const std::string name = subdir.filename().string();
        out.emplace_back(name, fs::absolute(preferred));
    }
}

// Stable order for the menu
std::sort(out.begin(), out.end(),
    [](const auto& a, const auto& b) { return a.first < b.first; });

return out;
}

/**
 * @brief Show a numeric menu and return the chosen index in [0, n-1].
 *
 * If the user presses Enter without typing a number, default_choice is used
 * (clamped to the valid range). If the input is invalid, the function will
 * re-prompt once and then fall back to default_choice.
 */
static std::size_t prompt_choice(std::size_t n, std::size_t default_choice = 0)
{
    auto clamp = [&](std::size_t v) { return (v < n ? v : std::min(default_choice, n ? n - 1 : 0)); };

    std::cout << "\nEnter a number and press Enter "
        "(blank = " << (default_choice + 1) << "): ";
    std::string line;
    if (!std::getline(std::cin, line)) {
        return clamp(default_choice);
    }

    if (line.empty()) return clamp(default_choice);

    try {
        int v = std::stoi(line);
        if (v >= 1 && static_cast<std::size_t>(v) <= n) {
            return static_cast<std::size_t>(v - 1);
        }
    }
    catch (...) {
        // fall through
    }

    std::cout << "Invalid input. Please enter a number between 1 and " << n << ".\n";
    std::cout << "Choice: ";
    if (!std::getline(std::cin, line) || line.empty()) {
        return clamp(default_choice);
    }
    try {
        int v = std::stoi(line);
        if (v >= 1 && static_cast<std::size_t>(v) <= n) {
            return static_cast<std::size_t>(v - 1);
        }
    }
    catch (...) {
        // ignore
    }
    return clamp(default_choice);
}


/**
* @brief Entry point for the KarstNSim application.
*
* This function takes an optional command-line argument specifying the instruction file.
* If no argument is provided, it scans Input_files/ for subfolders that contain
* "instructions.txt", shows a numeric menu, and lets the user choose.
*/
int main(int argc, char* argv[]) {
    std::string instructionFile;

    // Root folder containing example subdirectories
    const fs::path inputRoot = fs::path("../../../Input_files");

    if (argc > 1) {
        // CLI mode: explicit path passed by the user
        instructionFile = argv[1];
    }
    else {
        // Interactive mode: discover examples with an actual "instructions.txt"
        auto examples = discover_examples(inputRoot);

        if (!examples.empty()) {
            std::cout << "Detected example simulations under: "
                << fs::absolute(inputRoot).string() << "\n\n";
            for (std::size_t i = 0; i < examples.size(); ++i) {
                std::cout << "  [" << (i + 1) << "] " << examples[i].first
                    << "  ->  " << examples[i].second.filename().string() << "\n";
            }

            const std::size_t idx = prompt_choice(examples.size(), /*default_choice=*/0);
            instructionFile = examples[idx].second.string();
        }
        else {
            // Fallback to legacy behavior only if no subfolders were found
            const fs::path legacy = inputRoot / "instructions.txt";
            instructionFile = legacy.string();
            std::cout << "No example subfolders found. Falling back to: "
                << instructionFile << "\n";
        }
    }

    // --- Validate selected instruction file before parsing ---
    {
        const fs::path p = fs::path(instructionFile);
        if (!fs::exists(p)) {
            std::cerr << "Selected instruction file does not exist: " << p.string() << std::endl;
            std::cout << "\nPress Enter to exit...";
            std::cin.get();
            return 1;
        }
        if (fs::is_directory(p)) {
            std::cerr << "Selected instruction path is a directory, not a file: " << p.string() << std::endl;
            std::cout << "\nPress Enter to exit...";
            std::cin.get();
            return 1;
        }
    }

    std::cout << "Parsing instruction file: " << instructionFile << std::endl;

    ParseInputs inputParser;
    KarstNSim::ParamsSource params = inputParser.parse(instructionFile);

    std::cout << "Parsing completed (success). Running simulation..." << std::endl;
	static const char* KARSTNSIM_LOGO = R"ASCII(
	 _  __               _    _   _  ____   _           
	| |/ /__ _ _ __ ___ | |_ | \ | |/ ___| (_) _ __ ___ 
	| ' // _` | '__/ __|| __||  \| |\___ \ | || '_ ` _ \
	| . \ (_| | |  \__ \| |_ | |\  | ___) || || | | | | |
	|_|\_\__,_|_|  |___/ \__||_| \_||____/ |_||_| |_| |_|
			KarstNSim v1.3
	)ASCII";
	
	std::cout << "\n" << KARSTNSIM_LOGO << "\n\n";

    // Measure the start time
    auto startTime = std::chrono::high_resolution_clock::now();

    // Run simulation
    KarstNSim::run_simulation_full(params);

    // Measure the end time
    auto endTime = std::chrono::high_resolution_clock::now();

    // Calculate the elapsed time in seconds
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();

    std::cout << "Simulation completed successfully. Elapsed time: " << duration << " seconds." << std::endl;

    std::cout << "\n\nPress Enter to exit...";
    std::cin.get();  // Wait for Enter key

    return 0;
}
