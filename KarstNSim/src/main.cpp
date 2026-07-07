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
#include <exception>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <memory>
#include <system_error>

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
 * @brief Stream buffer duplicating all console output into two destination buffers.
 *
 * This helper is used to mirror std::cout and std::cerr into a persistent log file
 * while keeping the normal console output visible to the user.
 */
class TeeStreambuf : public std::streambuf {
public:
    /**
     * @brief Construct a stream buffer writing to two existing buffers.
     * @param first First destination buffer, typically the original console buffer.
     * @param second Second destination buffer, typically a file buffer.
     */
    TeeStreambuf(std::streambuf* first, std::streambuf* second)
        : first_(first), second_(second) {
    }

protected:
    int overflow(int ch) override {
        if (traits_type::eq_int_type(ch, traits_type::eof())) {
            return traits_type::not_eof(ch);
        }

        const char c = traits_type::to_char_type(ch);
        const int r1 = first_->sputc(c);
        const int r2 = second_->sputc(c);

        if (traits_type::eq_int_type(r1, traits_type::eof()) ||
            traits_type::eq_int_type(r2, traits_type::eof())) {
            return traits_type::eof();
        }

        return ch;
    }

    std::streamsize xsputn(const char* s, std::streamsize n) override {
        const std::streamsize n1 = first_->sputn(s, n);
        const std::streamsize n2 = second_->sputn(s, n);
        return std::min(n1, n2);
    }

    int sync() override {
        const int r1 = first_->pubsync();
        const int r2 = second_->pubsync();
        return (r1 == 0 && r2 == 0) ? 0 : -1;
    }

private:
    std::streambuf* first_;
    std::streambuf* second_;
};

/**
 * @brief Scoped redirection duplicating std::cout and std::cerr into a log file.
 *
 * The original stream buffers are restored automatically when this object is destroyed.
 * Output streams are put in unit-buffered mode to preserve the most recent messages
 * in case of an abrupt crash.
 */
class ScopedConsoleLog {
public:
    /**
     * @brief Open the log file and start mirroring std::cout and std::cerr.
     * @param log_path Path of the persistent log file.
     */
    explicit ScopedConsoleLog(const fs::path& log_path)
        : old_cout_(std::cout.rdbuf()), old_cerr_(std::cerr.rdbuf()) {

        log_file_.open(log_path.string(), std::ios::out);

        if (!log_file_.is_open()) {
            std::cerr << "Warning: cannot create console log file: "
                << log_path.string() << std::endl;
            return;
        }

        cout_tee_ = std::make_unique<TeeStreambuf>(old_cout_, log_file_.rdbuf());
        cerr_tee_ = std::make_unique<TeeStreambuf>(old_cerr_, log_file_.rdbuf());

        std::cout.rdbuf(cout_tee_.get());
        std::cerr.rdbuf(cerr_tee_.get());

        // Flush after each insertion to preserve the latest diagnostic messages.
        std::cout << std::unitbuf;
        std::cerr << std::unitbuf;
        log_file_ << std::unitbuf;
    }

    /**
     * @brief Restore the original console buffers.
     */
    ~ScopedConsoleLog() {
        std::cout.rdbuf(old_cout_);
        std::cerr.rdbuf(old_cerr_);
    }

    /**
     * @brief Return whether the log file was opened successfully.
     */
    bool is_open() const {
        return log_file_.is_open();
    }

private:
    std::ofstream log_file_;
    std::streambuf* old_cout_;
    std::streambuf* old_cerr_;
    std::unique_ptr<TeeStreambuf> cout_tee_;
    std::unique_ptr<TeeStreambuf> cerr_tee_;
};

/**
 * @brief Return a compact timestamp usable in a filename.
 */
static std::string current_timestamp_for_filename()
{
    const auto now = std::chrono::system_clock::now();
    const std::time_t t = std::chrono::system_clock::to_time_t(now);

    std::tm tm_snapshot;
#if defined(_WIN32)
    localtime_s(&tm_snapshot, &t);
#else
    localtime_r(&t, &tm_snapshot);
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm_snapshot, "%Y%m%d_%H%M%S");
    return oss.str();
}

/**
 * @brief Reads the main_repository entry from an instruction file without fully parsing it.
 *
 * This lightweight scan is only used to place the console log in the output directory
 * before the full parser starts. If the tag is absent or unreadable, an empty path is returned.
 *
 * @param instruction_file Path to the instruction file.
 * @return The path found after the "main_repository:" tag, or an empty path.
 */
static fs::path read_main_repository_for_log(const fs::path& instruction_file)
{
    std::ifstream in(instruction_file.string());
    if (!in.is_open()) {
        return fs::path();
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.find("//") == 0 || line.find("#") == 0 ||
            line.find("%") == 0 || line.find("::") == 0 ||
            line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        std::string tag;
        std::string value;
        iss >> tag >> value;

        if (tag == "main_repository:" && !value.empty()) {
            return fs::path(value);
        }
    }

    return fs::path();
}

/**
 * @brief Builds the path of the persistent console log file.
 *
 * The preferred location is main_repository/outputs, because it is already the
 * standard output location used by the simulation. If this cannot be determined,
 * the log is written next to the instruction file.
 *
 * @param instruction_file Path to the selected instruction file.
 * @return Path to the console log file.
 */
static fs::path make_console_log_path(const fs::path& instruction_file)
{
    fs::path log_dir;

    const fs::path main_repository = read_main_repository_for_log(instruction_file);
    if (!main_repository.empty()) {
        log_dir = main_repository / "outputs";
    }
    else {
        log_dir = instruction_file.parent_path();
    }

    std::error_code ec;
    fs::create_directories(log_dir, ec);

    if (ec) {
        log_dir = instruction_file.parent_path();
    }

    const std::string filename =
        "karstnsim_console_" + current_timestamp_for_filename() + ".log";

    return log_dir / filename;
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
            if (argc <= 1) {
                std::cout << "\n\nPress Enter to exit...";
                std::cin.get();
            }
            return 1;
        }
        if (fs::is_directory(p)) {
            std::cerr << "Selected instruction path is a directory, not a file: " << p.string() << std::endl;
            if (argc <= 1) {
                std::cout << "\n\nPress Enter to exit...";
                std::cin.get();
            }
            return 1;
        }
    }

    const fs::path log_path = make_console_log_path(fs::path(instructionFile));
    ScopedConsoleLog console_log(log_path);

    if (console_log.is_open()) {
        std::cout << "Console log file: "
            << fs::absolute(log_path).string()
            << std::endl;
    }

    try {
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
        KarstNSim v2.0
)ASCII";

        std::cout << "\n" << KARSTNSIM_LOGO << "\n\n";

        // Measure the start time.
        auto startTime = std::chrono::high_resolution_clock::now();

        // Run simulation.
        KarstNSim::run_simulation_full(params);

        // Measure the end time.
        auto endTime = std::chrono::high_resolution_clock::now();

        // Calculate the elapsed time in seconds.
        auto duration =
            std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();

        std::cout << "Simulation completed successfully. Elapsed time: "
            << duration << " seconds." << std::endl;

        if (argc <= 1) {
            std::cout << "\n\nPress Enter to exit...";
            std::cin.get();
        }

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\nKarstNSim stopped because an exception was thrown:\n"
            << e.what() << std::endl;

        std::cerr << "Console log file: "
            << fs::absolute(log_path).string()
            << std::endl;

        if (argc <= 1) {
            std::cout << "\n\nPress Enter to exit...";
            std::cin.get();
        }

        return 1;
    }
    catch (...) {
        std::cerr << "\nKarstNSim stopped because an unknown exception was thrown."
            << std::endl;

        std::cerr << "Console log file: "
            << fs::absolute(log_path).string()
            << std::endl;

        if (argc <= 1) {
            std::cout << "\n\nPress Enter to exit...";
            std::cin.get();
        }

        return 1;
    }
}
