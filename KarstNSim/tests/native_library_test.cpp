// Focused tests of the in-memory native backend (library.h).
//
// Covers: result/node identity, sequential, concurrent and nested seed isolation, noisy
// generation, invalid input, cancellation and work/point/edge limits with checkpoints
// inside the sampling, routing and section loops, and absence of filesystem and
// stdout/stderr side effects.

#include "KarstNSim/library.h"
#include "KarstNSim/job_context.h"
#include "native_region_factory.h"

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <unistd.h>

namespace fs = std::filesystem;
using KarstNSimTest::RegionRequest;
using KarstNSimTest::make_region_params;

namespace {

    int failures = 0;

    void require(bool condition, const std::string& message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    void run_test(const char* name, const std::function<void()>& test) {
        try {
            test();
            std::cout << "[PASS] " << name << std::endl;
        }
        catch (const std::exception& error) {
            ++failures;
            std::cout << "[FAIL] " << name << ": " << error.what() << std::endl;
        }
        if (KarstNSim::detail::current_job() != nullptr) {
            ++failures;
            std::cout << "[FAIL] " << name << ": a job context leaked on the test thread" << std::endl;
        }
    }

    // Exact serialization of everything a result carries, including fields that the
    // legacy text export omits (node_id, solved matrix) and full float bit patterns.
    std::string fingerprint(const std::vector<KarstNSim::KarstNetworkResult>& results) {
        std::ostringstream out;
        out << std::hexfloat;
        for (const KarstNSim::KarstNetworkResult& result : results) {
            out << "R " << result.segments.size() << ' ' << result.has_drift_properties << '\n';
            for (const std::string& name : result.vadose_property_names) out << name << ' ';
            out << '\n';
            for (const KarstNSim::ResultSegment& segment : result.segments) {
                for (const KarstNSim::ResultPoint* point : { &segment.start, &segment.end }) {
                    out << point->node_id << ' ' << point->p.x << ' ' << point->p.y << ' ' << point->p.z
                        << ' ' << point->cost << ' ' << point->equivalent_radius << ' ' << point->branch_id
                        << ' ' << point->vadose_flag << ' ' << point->external_drift << ' '
                        << point->kriging_weight;
                    for (float flag : point->vadose_flags) out << ' ' << flag;
                    out << '\n';
                }
            }
            for (const std::vector<int>& row : result.solved_connectivity_matrix) {
                for (int value : row) out << value << ' ';
                out << '\n';
            }
        }
        return out.str();
    }

    std::vector<KarstNSim::KarstNetworkResult> run(const RegionRequest& request,
        const KarstNSim::RunOptions& options = {}) {
        return KarstNSim::run_simulation_memory(make_region_params(request), options);
    }

    template <typename ErrorType>
    ErrorType expect_error(const std::function<void()>& call, KarstNSim::ErrorCode code) {
        try {
            call();
        }
        catch (const ErrorType& error) {
            require(error.code() == code, std::string("unexpected error code ") +
                KarstNSim::to_string(error.code()) + ": " + error.what());
            return error;
        }
        catch (const std::exception& error) {
            throw std::runtime_error(std::string("unexpected exception type: ") + error.what());
        }
        throw std::runtime_error(std::string("expected error ") + KarstNSim::to_string(code));
    }

    // Redirects file descriptors 1 and 2 into a temporary file for the lifetime of the object.
    class DescriptorCapture {
    public:
        DescriptorCapture() {
            std::cout.flush();
            std::cerr.flush();
            std::fflush(nullptr);
            file_ = std::tmpfile();
            require(file_ != nullptr, "tmpfile failed");
            saved_out_ = dup(STDOUT_FILENO);
            saved_err_ = dup(STDERR_FILENO);
            dup2(fileno(file_), STDOUT_FILENO);
            dup2(fileno(file_), STDERR_FILENO);
        }
        std::string finish() {
            std::cout.flush();
            std::cerr.flush();
            std::fflush(nullptr);
            dup2(saved_out_, STDOUT_FILENO);
            dup2(saved_err_, STDERR_FILENO);
            close(saved_out_);
            close(saved_err_);
            std::string content;
            std::rewind(file_);
            char buffer[512];
            std::size_t count = 0;
            while ((count = std::fread(buffer, 1, sizeof buffer, file_)) > 0) {
                content.append(buffer, count);
            }
            std::fclose(file_);
            file_ = nullptr;
            return content;
        }
        ~DescriptorCapture() {
            if (file_ != nullptr) {
                finish();
            }
        }
    private:
        std::FILE* file_ = nullptr;
        int saved_out_ = -1;
        int saved_err_ = -1;
    };

    // Records the log state each time the job polls cancellation.
    struct PollProbe {
        const std::ostringstream* log = nullptr;
        std::map<std::string, int> polls_by_stage;
        std::string cancel_stage; //!< Cancel at the second poll seen in this stage.
        int cancel_after = 2;
        bool cancelled = false;

        static std::string stage_of(const std::string& text) {
            struct Marker { const char* begin; const char* end; const char* stage; };
            static const Marker markers[] = {
                { "STEP 4 - Simulation of conduit sections", "STEP 4 completed", "sections" },
                { "STEP 2 - Simulation of the karst network skeleton", "STEP 2 completed", "routing" },
                { "1.1 Generating sampling cloud", "automatically sampled points", "sampling" },
            };
            for (const Marker& marker : markers) {
                if (text.find(marker.begin) != std::string::npos) {
                    return text.find(marker.end) == std::string::npos ? marker.stage : "other";
                }
            }
            return "other";
        }

        static bool poll(void* user) {
            PollProbe& probe = *static_cast<PollProbe*>(user);
            const std::string stage = stage_of(probe.log->str());
            const int count = ++probe.polls_by_stage[stage];
            if (!probe.cancel_stage.empty() && stage == probe.cancel_stage && count >= probe.cancel_after) {
                probe.cancelled = true;
                return true;
            }
            return false;
        }
    };

    RegionRequest sections_request() {
        RegionRequest request;
        request.seed = 3;
        request.sink_count = 10;
        request.cycles = 20;
        request.poisson_radius = 0.06f;
        request.sections = true; // large enough for several in-loop polls in every stage
        return request;
    }
}

int main() {
    // Reference results computed once, sequentially.
    std::vector<std::string> reference(9);
    RegionRequest base;
    const auto reference_for = [&](int seed, bool noise) {
        RegionRequest request = base;
        request.seed = seed;
        request.noise = noise;
        return fingerprint(run(request));
    };

    run_test("region factory produces a network with native node identity", [&] {
        const auto results = run(base);
        require(results.size() == 1, "one result per iteration");
        const KarstNSim::KarstNetworkResult& result = results[0];
        require(!result.segments.empty(), "network has segments");

        std::map<std::uint32_t, ::Vector3> node_position;
        std::uint32_t max_id = 0;
        for (const KarstNSim::ResultSegment& segment : result.segments) {
            for (const KarstNSim::ResultPoint* point : { &segment.start, &segment.end }) {
                require(point->node_id != UINT32_MAX, "every exported point carries a node_id");
                max_id = std::max(max_id, point->node_id);
                auto inserted = node_position.emplace(point->node_id, point->p);
                const ::Vector3& known = inserted.first->second;
                require(known.x == point->p.x && known.y == point->p.y && known.z == point->p.z,
                    "one node_id always maps to the same exact coordinates");
            }
            require(segment.start.node_id != segment.end.node_id, "segments join two distinct nodes");
        }
        require(max_id < node_position.size() * 4 + 16, "node ids are compact skeleton indices");

        // Legacy text serialization is unchanged: no node_id column, 5-decimal coordinates.
        const std::string text = result.to_string();
        require(text.rfind("Index\tX\tY\tZ\tcost\tequivalent_radius\tbranch_id\n", 0) == 0,
            "legacy header unchanged");
        std::istringstream lines(text);
        std::string header, first;
        std::getline(lines, header);
        std::getline(lines, first);
        std::istringstream fields(first);
        int index = -1;
        double x = 0;
        fields >> index >> x;
        require(index == 0 && std::fabs(x - result.segments[0].start.p.x) <= 5e-6,
            "text export is a rounding of the unrounded in-memory coordinate");

        require(result.solved_connectivity_matrix.size() == static_cast<std::size_t>(base.sink_count),
            "solved connectivity matrix returned in memory");
        for (const auto& row : result.solved_connectivity_matrix) {
            require(row.size() == static_cast<std::size_t>(base.spring_count), "matrix column count");
        }
    });

    run_test("sequential runs are isolated by seed", [&] {
        reference[1] = reference_for(1, false);
        reference[2] = reference_for(2, false);
        require(reference[1] != reference[2], "different seeds give different networks");
        require(reference_for(1, false) == reference[1], "seed 1 unchanged after running seed 2");
        reference[5] = reference_for(5, true);
        reference[6] = reference_for(6, true);
        require(reference_for(2, false) == reference[2], "seed 2 unchanged after noisy runs");
    });

    run_test("noisy generation is seeded per job", [&] {
        require(reference[5] != reference_for(5, false), "noise changes the network");
        require(reference_for(5, true) == reference[5], "noisy run is reproducible");
        require(reference[5] != reference[6], "noise permutation follows the job seed");
    });

    run_test("concurrent jobs match sequential references", [&] {
        struct Job { int seed; bool noise; std::string result{}; std::string error{}; };
        std::vector<Job> jobs = { {1, false, {}, {}}, {2, false, {}, {}}, {5, true, {}, {}}, {6, true, {}, {}}, {1, false, {}, {}}, {5, true, {}, {}},
            {2, false, {}, {}}, {6, true, {}, {}} };
        std::vector<std::thread> threads;
        for (Job& job : jobs) {
            threads.emplace_back([&job, &base] {
                try {
                    RegionRequest request = base;
                    request.seed = job.seed;
                    request.noise = job.noise;
                    job.result = fingerprint(run(request));
                }
                catch (const std::exception& error) {
                    job.error = error.what();
                }
            });
        }
        for (std::thread& thread : threads) thread.join();
        for (const Job& job : jobs) {
            require(job.error.empty(), "concurrent job failed: " + job.error);
            require(job.result == reference[job.seed], "concurrent seed " + std::to_string(job.seed) +
                " differs from its sequential reference");
        }
    });

    run_test("nested jobs keep the outer job state, including after an inner exception", [&] {
        struct Nested {
            RegionRequest inner;
            std::string inner_result;
            bool inner_error_seen = false;
            int calls = 0;
            static bool poll(void* user) {
                Nested& nested = *static_cast<Nested*>(user);
                if (++nested.calls == 2) {
                    nested.inner_result = fingerprint(run(nested.inner));
                    try {
                        KarstNSim::ParamsSource invalid = make_region_params(nested.inner);
                        invalid.nghb_count = 0;
                        KarstNSim::run_simulation_memory(invalid);
                    }
                    catch (const KarstNSim::InvalidInputError&) {
                        nested.inner_error_seen = true;
                    }
                }
                return false;
            }
        } nested;
        nested.inner = base;
        nested.inner.seed = 2;
        KarstNSim::RunOptions options;
        options.is_cancelled = &Nested::poll;
        options.user = &nested;
        RegionRequest outer = base;
        outer.seed = 1;
        const std::string outer_result = fingerprint(run(outer, options));
        require(nested.calls > 2, "outer job polled cancellation repeatedly");
        require(nested.inner_error_seen, "inner invalid job raised InvalidInputError");
        require(nested.inner_result == reference[2], "inner job matches its own reference");
        require(outer_result == reference[1], "outer job unaffected by the nested jobs");
    });

    run_test("invalid input is rejected with InvalidInputError", [&] {
        using KarstNSim::ErrorCode;
        const auto invalid = [&](const std::function<void(KarstNSim::ParamsSource&)>& change) {
            KarstNSim::ParamsSource params = make_region_params(base);
            change(params);
            expect_error<KarstNSim::InvalidInputError>(
                [&] { KarstNSim::run_simulation_memory(params); }, ErrorCode::InvalidInput);
        };
        invalid([](KarstNSim::ParamsSource& p) { p.create_vset_sampling = true; });
        invalid([](KarstNSim::ParamsSource& p) { p.create_nghb_graph = true; });
        invalid([](KarstNSim::ParamsSource& p) { p.create_nghb_graph_property = true; });
        invalid([](KarstNSim::ParamsSource& p) { p.create_grid = true; });
        invalid([](KarstNSim::ParamsSource& p) { p.use_user_connectivity_matrix = true; });
        invalid([](KarstNSim::ParamsSource& p) {
            p.use_user_connectivity_matrix = true;
            p.connectivity_matrix = { {1, 1}, {1, 1} }; // 2 rows for 4 sinks
        });
        invalid([](KarstNSim::ParamsSource& p) {
            p.use_user_connectivity_matrix = true;
            p.connectivity_matrix.assign(p.sinks.size(), std::vector<int>(p.springs.size(), 3));
        });
        invalid([](KarstNSim::ParamsSource& p) { p.connectivity_matrix.assign(p.sinks.size(), std::vector<int>(p.springs.size(), 2)); });
        invalid([](KarstNSim::ParamsSource& p) { p.sinks.clear(); p.propsinksindex.clear(); p.propsinksorder.clear(); });
        invalid([](KarstNSim::ParamsSource& p) { p.poisson_radius = 0.0f; });
        invalid([](KarstNSim::ParamsSource& p) { p.nghb_count = 0; });
        invalid([](KarstNSim::ParamsSource& p) { p.number_of_iterations = 0; });
        invalid([](KarstNSim::ParamsSource& p) { p.propsinksindex[0] = 7; });
        invalid([](KarstNSim::ParamsSource& p) { p.use_amplification = true; p.min_distance_amplification = 9; p.max_distance_amplification = 1; });
    });

    run_test("explicit in-memory connectivity matrix is used", [&] {
        KarstNSim::ParamsSource params = make_region_params(base);
        params.use_user_connectivity_matrix = true;
        params.connectivity_matrix.assign(params.sinks.size(), std::vector<int>(params.springs.size(), 0));
        for (auto& row : params.connectivity_matrix) row[1] = 1; // every sink drains to spring 2
        const auto results = KarstNSim::run_simulation_memory(params);
        require(results.size() == 1 && !results[0].segments.empty(), "network generated");
        for (const auto& row : results[0].solved_connectivity_matrix) {
            require(row[0] == 0 && row[1] == 1, "solved matrix keeps the user connectivity");
        }
        require(fingerprint(results) != reference[1], "user connectivity changes the network");
    });

    run_test("no route raises NoRouteError", [&] {
        KarstNSim::ParamsSource params = make_region_params(base);
        for (::Vector3& spring : params.springs) spring.z = base.origin_z + 0.95f * base.size_z;
        for (::Vector3& sink : params.sinks) sink.z = base.origin_z + 0.3f * base.size_z;
        const auto error = expect_error<KarstNSim::NoRouteError>(
            [&] { KarstNSim::run_simulation_memory(params); }, KarstNSim::ErrorCode::NoRoute);
        require(error.iteration() == 0, "iteration index reported");
    });

    run_test("work, point and edge limits are enforced deterministically", [&] {
        using KarstNSim::ErrorCode;
        KarstNSim::RunOptions options;
        options.work_limit = 20000;
        std::uint64_t work_used = 0;
        options.work_used = &work_used;
        const auto first = expect_error<KarstNSim::LimitExceededError>([&] { run(base, options); },
            ErrorCode::WorkLimitExceeded);
        const auto second = expect_error<KarstNSim::LimitExceededError>([&] { run(base, options); },
            ErrorCode::WorkLimitExceeded);
        require(first.limit() == 20000 && first.requested() == 20001, "work limit reported at the first excess unit");
        require(std::string(first.what()) == second.what(), "work limit failure is reproducible");
        require(work_used == first.requested(), "actual work is reported on exceptional exit");

        options = {};
        options.max_points = 200;
        const auto points = expect_error<KarstNSim::LimitExceededError>([&] { run(base, options); },
            ErrorCode::PointLimitExceeded);
        require(points.requested() == 201, "point admission checked before insertion");

        options = {};
        options.max_edges = 1000;
        const auto edges = expect_error<KarstNSim::LimitExceededError>([&] { run(base, options); },
            ErrorCode::EdgeLimitExceeded);
        require(edges.requested() > 1000, "edge admission checked before adjacency allocation");

        options = {};
        options.work_limit = 0;
        expect_error<KarstNSim::LimitExceededError>([&] { run(base, options); }, ErrorCode::WorkLimitExceeded);

        options = {};
        options.work_used = &work_used;
        require(fingerprint(run(base, options)) == reference[1], "state is clean after limit failures");
        require(work_used > 20000 && work_used <= options.work_limit, "actual work is reported on successful exit");
    });

    run_test("cancellation is polled inside sampling, routing and section loops", [&] {
        std::ostringstream log;
        PollProbe probe;
        probe.log = &log;
        KarstNSim::RunOptions options;
        options.log = &log;
        options.is_cancelled = &PollProbe::poll;
        options.user = &probe;
        const auto results = run(sections_request(), options);
        require(!results.empty(), "section request succeeds");
        bool has_radius = false;
        for (const auto& segment : results[0].segments) {
            has_radius = has_radius || segment.start.equivalent_radius != -99999.0f;
        }
        require(has_radius, "equivalent radii simulated in memory");
        // Each stage has one boundary checkpoint (two for sections); more polls prove in-loop checkpoints.
        require(probe.polls_by_stage["sampling"] >= 2, "in-loop poll during sampling");
        require(probe.polls_by_stage["routing"] >= 2, "in-loop poll during routing");
        require(probe.polls_by_stage["sections"] >= 3, "in-loop poll during section simulation");
    });

    run_test("cancellation stops the job inside each stage", [&] {
        for (const char* stage : { "sampling", "routing", "sections" }) {
            std::ostringstream log;
            PollProbe probe;
            probe.log = &log;
            probe.cancel_stage = stage;
            probe.cancel_after = std::string(stage) == "sections" ? 3 : 2;
            KarstNSim::RunOptions options;
            options.log = &log;
            options.is_cancelled = &PollProbe::poll;
            options.user = &probe;
            expect_error<KarstNSim::CancelledError>([&] { run(sections_request(), options); },
                KarstNSim::ErrorCode::Cancelled);
            require(probe.cancelled, "probe requested cancellation");
            require(PollProbe::stage_of(log.str()) == stage, std::string("job stopped inside ") + stage);
        }
        KarstNSim::RunOptions options;
        options.is_cancelled = [](void*) { return true; };
        expect_error<KarstNSim::CancelledError>([&] { run(base, options); }, KarstNSim::ErrorCode::Cancelled);
    });

    run_test("no filesystem access and no stdout/stderr output", [&] {
        const fs::path original = fs::current_path();
        const fs::path sandbox = fs::temp_directory_path() /
            ("karstnsim-native-" + std::to_string(::getpid()));
        fs::remove_all(sandbox);
        fs::create_directories(sandbox / "inputs");
        {
            // A connectivity_matrix.txt that would be rejected if it were read.
            std::ofstream bad(sandbox / "inputs" / "connectivity_matrix.txt");
            bad << "not a matrix\n";
        }
        const auto before = std::distance(fs::recursive_directory_iterator(sandbox), fs::recursive_directory_iterator());
        fs::current_path(sandbox);

        KarstNSim::ParamsSource params = make_region_params(base);
        params.save_repository = sandbox.string();
        params.simulation_input_dir = (sandbox / "inputs").string();
        params.use_user_connectivity_matrix = true;
        params.connectivity_matrix.assign(params.sinks.size(), std::vector<int>(params.springs.size(), 2));

        std::streambuf* const cout_buffer = std::cout.rdbuf();
        std::streambuf* const cerr_buffer = std::cerr.rdbuf();
        const auto cout_flags = std::cout.flags();
        std::string captured;
        std::string outcome;
        {
            DescriptorCapture capture;
            try {
                const auto quiet = KarstNSim::run_simulation_memory(params);
                RegionRequest sections = sections_request();
                KarstNSim::run_simulation_memory(make_region_params(sections));
                KarstNSim::ParamsSource invalid = params;
                invalid.nghb_count = 0;
                try { KarstNSim::run_simulation_memory(invalid); }
                catch (const KarstNSim::InvalidInputError&) {}
                outcome = quiet.empty() ? "empty" : "ok";
            }
            catch (const std::exception& error) {
                outcome = error.what();
            }
            captured = capture.finish();
        }
        fs::current_path(original);
        const auto after = std::distance(fs::recursive_directory_iterator(sandbox), fs::recursive_directory_iterator());
        fs::remove_all(sandbox);

        require(outcome == "ok", "memory jobs succeeded: " + outcome);
        require(captured.empty(), "nothing written to file descriptors 1/2: " + captured.substr(0, 200));
        require(std::cout.rdbuf() == cout_buffer && std::cerr.rdbuf() == cerr_buffer, "std stream buffers untouched");
        require(std::cout.flags() == cout_flags, "std::cout formatting untouched");
        require(before == after, "no file or directory created in the working or output directories");
    });

    run_test("log stream receives engine messages and keeps its formatting", [&] {
        std::ostringstream log;
        log << std::scientific << std::setprecision(2);
        const auto flags = log.flags();
        KarstNSim::RunOptions options;
        options.log = &log;
        const auto results = run(base, options);
        require(log.str().find("STEP 1 - Generation of cost graph") != std::string::npos, "engine log captured");
        require(log.flags() == flags && log.precision() == 2, "caller stream formatting restored");
        require(fingerprint(results) == reference[1], "logging does not change results");
    });

    run_test("filesystem guard throws inside an in-memory job", [&] {
        KarstNSim::detail::JobContext job;
        job.allow_filesystem = false;
        bool threw = false;
        {
            KarstNSim::detail::ScopedJob scope(job);
            try {
                KarstNSim::save_pointset("guard.txt", "/nonexistent-karstnsim-dir", {});
            }
            catch (const std::logic_error&) {
                threw = true;
            }
        }
        require(threw, "save_pointset refused to write");
        require(!fs::exists("/nonexistent-karstnsim-dir"), "no directory created");
    });

    std::cout << (failures == 0 ? "All native library tests passed." : "Native library tests FAILED.")
        << std::endl;
    return failures == 0 ? 0 : 1;
}
