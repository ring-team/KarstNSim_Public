#pragma once
// Bridge between the regional generator and the in-memory KarstNSim engine.
// This is the only FinalBuildCaves translation unit that includes KarstNSim headers.

#include "fbs/caves.hpp"

#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

namespace fbs::caves::detail {

// Bounded routing graph in region-local meters. nodes[0] is the inlet (native
// sink) and nodes[1] the outlet (native spring); the remaining nodes are support
// points. Every edge has already passed the regional slope and clearance tests,
// so any route KarstNSim selects on it is admissible.
struct RouteGraph {
    std::vector<Vec3> nodes;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
    Bounds domain;
};

struct RouteSettings {
    std::int32_t seed = 1;
    std::uint32_t neighbors = 20;
    bool noise = true;
    float noise_weight = 2.0f;
    int noise_frequency = 4;
    float vertical_stretch = 2.0f;
};

struct RouteLimits {
    std::uint64_t work_limit = 0, max_points = 0, max_edges = 0;
    const std::function<bool()>* cancelled = nullptr;
    // Receives the native engine work counter, also when route_karst throws.
    std::uint64_t* work_used = nullptr;
};

struct RouteResult {
    bool found = false;
    std::vector<std::uint32_t> path;             // node indices from 0 to 1
    std::vector<std::uint32_t> native_node_ids;  // KarstNSim skeleton node id of each path node
};

// Runs KarstNSim::run_simulation_memory on the graph. Returns found=false when the
// engine reports no inlet-outlet route. Throws Error(cancelled) or Error(budget)
// when the engine stops on the caller's limits, Error(internal) on engine failures
// or when the returned skeleton is not a path over the supplied graph.
RouteResult route_karst(const RouteGraph& graph, const RouteSettings& settings, const RouteLimits& limits);

} // namespace fbs::caves::detail
