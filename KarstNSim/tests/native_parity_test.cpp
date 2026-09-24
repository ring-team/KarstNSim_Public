// Parity between the in-memory API and the legacy executable on the regression corpus.
//
// Usage: native_parity_test <regression-results-dir>
//
// For every case work root written by tests/regression/run_regression.py (whose exports
// were just checked against the golden hashes), the instruction file is parsed with the
// legacy parser, the connectivity matrix file is loaded into
// ParamsSource::connectivity_matrix, filesystem-only exports are switched off, and
// run_simulation_memory() is called. Every iteration's KarstNetworkResult::to_string()
// must equal the legacy <name>_karst.txt byte for byte, and the in-memory solved
// connectivity matrix must equal the legacy <name>_connectivity_matrix.txt export.

#include "KarstNSim/library.h"
#include "KarstNSim/parse_inputs.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

    std::string read_file(const fs::path& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) {
            throw std::runtime_error("cannot read " + path.string());
        }
        std::ostringstream content;
        content << in.rdbuf();
        return content.str();
    }

    std::vector<std::vector<int>> read_matrix(const fs::path& path) {
        std::vector<std::vector<int>> rows;
        std::istringstream in(read_file(path));
        std::string line;
        while (std::getline(in, line)) {
            if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
            std::istringstream values(line);
            std::vector<int> row;
            int value = 0;
            while (values >> value) row.push_back(value);
            rows.push_back(row);
        }
        return rows;
    }

    std::string matrix_text(const std::vector<std::vector<int>>& matrix) {
        // Same layout as write_files.cpp save_connectivity_matrix.
        std::ostringstream out;
        for (const std::vector<int>& row : matrix) {
            for (std::size_t j = 0; j < row.size(); ++j) {
                out << row[j] << (j + 1 < row.size() ? "\t" : "");
            }
            out << "\n";
        }
        return out.str();
    }

    int check_case(const fs::path& work) {
        const fs::path original = fs::current_path();
        fs::current_path(work);
        int compared = 0;
        try {
            std::ostringstream parse_log;
            std::streambuf* const cout_buffer = std::cout.rdbuf(parse_log.rdbuf()); // silence the legacy parser only
            KarstNSim::ParamsSource params;
            try {
                ParseInputs parser;
                params = parser.parse("instructions.txt");
            }
            catch (...) {
                std::cout.rdbuf(cout_buffer);
                throw;
            }
            std::cout.rdbuf(cout_buffer);

            // Sections-only runs do no routing, so neither runner produces a solved matrix.
            const bool solved_matrix = params.create_solved_connectivity_matrix && !params.sections_simulation_only;
            params.create_vset_sampling = false;
            params.create_nghb_graph = false;
            params.create_nghb_graph_property = false;
            params.create_grid = false;
            if (params.use_user_connectivity_matrix && !params.sections_simulation_only) {
                params.connectivity_matrix = read_matrix(fs::path(params.simulation_input_dir) / "connectivity_matrix.txt");
            }

            KarstNSim::RunOptions options;
            options.work_limit = UINT64_MAX;
            options.max_points = UINT64_MAX;
            options.max_edges = UINT64_MAX;
            const std::vector<KarstNSim::KarstNetworkResult> results =
                KarstNSim::run_simulation_memory(params, options);

            if (results.size() != static_cast<std::size_t>(params.number_of_iterations)) {
                throw std::runtime_error("unexpected number of results");
            }
            for (int i = 0; i < params.number_of_iterations; ++i) {
                const std::string name = KarstNSim::detail::iteration_name(params, i);
                const std::string legacy = read_file(fs::path("outputs") / (name + "_karst.txt"));
                if (results[i].to_string() != legacy) {
                    throw std::runtime_error(name + "_karst.txt differs from the in-memory result");
                }
                ++compared;
                for (const KarstNSim::ResultSegment& segment : results[i].segments) {
                    if (segment.start.node_id == UINT32_MAX || segment.end.node_id == UINT32_MAX) {
                        throw std::runtime_error(name + ": missing node_id");
                    }
                }
                const fs::path legacy_matrix = fs::path("outputs") / (name + "_connectivity_matrix.txt");
                if (solved_matrix) {
                    if (matrix_text(results[i].solved_connectivity_matrix) != read_file(legacy_matrix)) {
                        throw std::runtime_error(name + "_connectivity_matrix.txt differs from the in-memory matrix");
                    }
                    ++compared;
                }
                else if (!results[i].solved_connectivity_matrix.empty() || fs::exists(legacy_matrix)) {
                    throw std::runtime_error(name + ": unexpected solved connectivity matrix");
                }
            }
        }
        catch (...) {
            fs::current_path(original);
            throw;
        }
        fs::current_path(original);
        return compared;
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: native_parity_test <regression-results-dir>" << std::endl;
        return 2;
    }
    const fs::path root = fs::absolute(argv[1]);
    std::vector<fs::path> cases;
    for (const auto& entry : fs::directory_iterator(root)) {
        if (entry.is_directory() && fs::exists(entry.path() / "instructions.txt")) {
            cases.push_back(entry.path());
        }
    }
    std::sort(cases.begin(), cases.end());
    if (cases.empty()) {
        std::cerr << "no regression case under " << root << std::endl;
        return 1;
    }

    int failures = 0;
    int compared = 0;
    for (const fs::path& work : cases) {
        try {
            const int count = check_case(work);
            compared += count;
            std::cout << "[PASS] " << work.filename().string() << " (" << count << " export(s))" << std::endl;
        }
        catch (const std::exception& error) {
            ++failures;
            std::cout << "[FAIL] " << work.filename().string() << ": " << error.what() << std::endl;
        }
    }
    std::cout << (failures == 0 ? "PASS: " : "FAIL: ") << cases.size() << " case(s), " << compared
        << " legacy export(s) reproduced in memory" << std::endl;
    return failures == 0 ? 0 : 1;
}
