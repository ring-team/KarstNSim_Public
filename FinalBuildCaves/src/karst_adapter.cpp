// KarstNSim adapter: one geological route per call.
//
// The regional generator supplies an explicit input neighbour graph (sampling
// points plus indexed edges) whose every edge already satisfies the slope limit
// and all reserved-volume clearances. KarstNSim computes its geological edge costs
// (anisotropic distance, vadose gradient term, optional per-job simplex noise) on
// that graph and returns the karst skeleton. Topology is rebuilt from the native
// skeleton node ids. Skeleton node coordinates are exact copies of the supplied
// float sample coordinates, so each node id is resolved to the supplied node by
// exact bit equality (the supplied float coordinates are checked to be unique);
// no rounding, tolerance or merging is involved.

#include "karst_adapter.hpp"

#include "KarstNSim/library.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <new>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace fbs::caves::detail {
namespace {

struct FloatKey {
    std::uint32_t x, y, z;
    bool operator==(const FloatKey& o) const { return x == o.x && y == o.y && z == o.z; }
};
struct FloatKeyHash {
    std::size_t operator()(const FloatKey& k) const {
        std::uint64_t h = (static_cast<std::uint64_t>(k.x) << 32) ^ k.y;
        h ^= static_cast<std::uint64_t>(k.z) * 0x9e3779b97f4a7c15ULL;
        h ^= h >> 29;
        h *= 0xbf58476d1ce4e5b9ULL;
        h ^= h >> 32;
        return static_cast<std::size_t>(h);
    }
};
FloatKey key_of(const ::Vector3& v) {
    FloatKey k{};
    float x = v.x, y = v.y, z = v.z;
    if (x == 0) x = 0;
    if (y == 0) y = 0;
    if (z == 0) z = 0;
    std::memcpy(&k.x, &x, 4);
    std::memcpy(&k.y, &y, 4);
    std::memcpy(&k.z, &z, 4);
    return k;
}
::Vector3 to_native(Vec3 v) {
    return ::Vector3(static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z));
}
std::uint64_t edge_key(std::uint32_t a, std::uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<std::uint64_t>(a) << 32) | b;
}
bool cancelled_trampoline(void* user) {
    const auto* fn = static_cast<const std::function<bool()>*>(user);
    return fn && *fn && (*fn)();
}

KarstNSim::Surface flat_roof(const Bounds& d) {
    // A flat topographic surface just above the routing domain: every support
    // point is underground, so the vadose (gravity) cost term applies everywhere.
    const float mx = static_cast<float>(0.1 * (d.maximum.x - d.minimum.x) + 1.0);
    const float my = static_cast<float>(0.1 * (d.maximum.y - d.minimum.y) + 1.0);
    const float z = static_cast<float>(d.maximum.z + 1.0);
    const float x0 = static_cast<float>(d.minimum.x) - mx, x1 = static_cast<float>(d.maximum.x) + mx;
    const float y0 = static_cast<float>(d.minimum.y) - my, y1 = static_cast<float>(d.maximum.y) + my;
    std::vector<::Vector3> nodes{{x0, y0, z}, {x1, y0, z}, {x0, y1, z}, {x1, y1, z}};
    std::vector<KarstNSim::Triangle> triangles{{0, 1, 2}, {1, 3, 2}};
    return KarstNSim::Surface(nodes, triangles, "fbs_region_roof");
}

} // namespace

RouteResult route_karst(const RouteGraph& graph, const RouteSettings& settings, const RouteLimits& limits) {
    RouteResult out;
    const std::size_t n = graph.nodes.size();
    if (n < 2 || n > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        graph.edges.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw Error(ErrorCode::internal, "route graph size outside the native index range");

    // Exact float identities must be unique, otherwise the engine would merge nodes.
    std::unordered_map<FloatKey, std::uint32_t, FloatKeyHash> by_key;
    by_key.reserve(n * 2);
    std::vector<::Vector3> native_nodes;
    native_nodes.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        native_nodes.push_back(to_native(graph.nodes[i]));
        if (!by_key.emplace(key_of(native_nodes.back()), static_cast<std::uint32_t>(i)).second)
            throw Error(ErrorCode::internal, "route graph has coincident native node coordinates");
    }
    std::unordered_set<std::uint64_t> edge_set;
    edge_set.reserve(graph.edges.size() * 2);
    std::vector<::Vector2i> native_edges;
    native_edges.reserve(graph.edges.size());
    for (const auto& e : graph.edges) {
        if (e.first >= n || e.second >= n || e.first == e.second)
            throw Error(ErrorCode::internal, "route graph edge references an invalid node");
        edge_set.insert(edge_key(e.first, e.second));
        native_edges.emplace_back(static_cast<int>(e.first), static_cast<int>(e.second));
    }

    KarstNSim::ParamsSource p;
    p.karstic_network_name = "fbs_region_route";
    p.selected_seed = settings.seed;
    p.number_of_iterations = 1;
    p.vary_seed = false;
    const Bounds& d = graph.domain;
    p.domain = KarstNSim::Box(to_native(d.minimum),
                              ::Vector3(static_cast<float>(d.maximum.x - d.minimum.x), 0.0f, 0.0f),
                              ::Vector3(0.0f, static_cast<float>(d.maximum.y - d.minimum.y), 0.0f),
                              ::Vector3(0.0f, 0.0f, static_cast<float>(d.maximum.z - d.minimum.z)), 1, 1, 1);
    p.topo_surface = flat_roof(d);

    p.use_sampling_points = true;
    p.sampling_points.assign(native_nodes.begin() + 2, native_nodes.end());
    p.use_input_nghb_graph = true;
    p.input_nghb_graph = KarstNSim::InputGraph(native_nodes, native_edges);
    p.nghb_count = static_cast<int>(settings.neighbors);

    p.sinks = {native_nodes[0]};
    p.springs = {native_nodes[1]};
    p.propsinksindex = {1};
    p.propsinksorder = {1};
    p.propspringsindex = {1};
    p.propspringssurfindex = {0}; // no water table: vadose outlet channel
    p.allow_single_outlet_connection = true;
    p.use_user_connectivity_matrix = false;

    p.fraction_karst_perm = 0.9f;
    p.gamma = 2.0f;
    p.vertical_distance_stretching_factor = settings.vertical_stretch;
    p.use_amplification = false;
    p.nb_cycles = 0;
    p.use_noise = settings.noise;
    p.use_noise_on_all = settings.noise;
    p.noise_frequency = settings.noise_frequency;
    p.noise_octaves = 2;
    p.noise_weight = settings.noise_weight;
    p.simulate_sections = false;

    p.create_vset_sampling = false;
    p.create_nghb_graph = false;
    p.create_nghb_graph_property = false;
    p.create_grid = false;
    p.create_solved_connectivity_matrix = false;

    KarstNSim::RunOptions options;
    options.work_limit = limits.work_limit;
    options.max_points = limits.max_points;
    options.max_edges = limits.max_edges;
    options.is_cancelled = limits.cancelled && *limits.cancelled ? &cancelled_trampoline : nullptr;
    options.user = const_cast<std::function<bool()>*>(limits.cancelled);
    options.log = nullptr;
    options.work_used = limits.work_used;

    std::vector<KarstNSim::KarstNetworkResult> results;
    try {
        results = KarstNSim::run_simulation_memory(std::move(p), options);
    } catch (const KarstNSim::Error& e) {
        switch (e.code()) {
        case KarstNSim::ErrorCode::Cancelled: throw Error(ErrorCode::cancelled, "generation cancelled in KarstNSim");
        case KarstNSim::ErrorCode::WorkLimitExceeded:
        case KarstNSim::ErrorCode::PointLimitExceeded:
        case KarstNSim::ErrorCode::EdgeLimitExceeded:
            throw Error(ErrorCode::budget, std::string("KarstNSim limit: ") + e.what());
        case KarstNSim::ErrorCode::NoRoute: return out;
        default: throw Error(ErrorCode::internal, std::string("KarstNSim failure: ") + e.what());
        }
    } catch (const std::bad_alloc&) {
        throw Error(ErrorCode::budget, "KarstNSim ran out of memory");
    }
    if (results.empty()) return out;

    // Rebuild the skeleton from native node ids.
    std::unordered_map<std::uint32_t, std::uint32_t> native_to_local;
    std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> adjacency;
    auto resolve = [&](const KarstNSim::ResultPoint& rp) {
        const auto it = by_key.find(key_of(rp.p));
        if (it == by_key.end())
            throw Error(ErrorCode::internal, "KarstNSim returned a node that is not in the supplied graph");
        if (rp.node_id != std::numeric_limits<std::uint32_t>::max()) {
            const auto [slot, inserted] = native_to_local.emplace(rp.node_id, it->second);
            if (!inserted && slot->second != it->second)
                throw Error(ErrorCode::internal, "KarstNSim node id maps to two supplied nodes");
        }
        return it->second;
    };
    std::vector<std::uint32_t> native_of(n, std::numeric_limits<std::uint32_t>::max());
    for (const auto& segment : results.front().segments) {
        const std::uint32_t a = resolve(segment.start), b = resolve(segment.end);
        native_of[a] = segment.start.node_id;
        native_of[b] = segment.end.node_id;
        if (a == b) continue;
        if (!edge_set.count(edge_key(a, b)))
            throw Error(ErrorCode::internal, "KarstNSim route uses an edge absent from the input graph");
        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
    }
    if (!adjacency.count(0) || !adjacency.count(1)) return out;

    // Shortest hop path from inlet to outlet over the returned skeleton.
    std::unordered_map<std::uint32_t, std::uint32_t> parent;
    std::queue<std::uint32_t> queue;
    parent[0] = 0;
    queue.push(0);
    while (!queue.empty() && !parent.count(1)) {
        const std::uint32_t u = queue.front();
        queue.pop();
        auto& next = adjacency[u];
        std::sort(next.begin(), next.end());
        for (const std::uint32_t v : next)
            if (parent.emplace(v, u).second) queue.push(v);
    }
    if (!parent.count(1)) return out;
    std::vector<std::uint32_t> reversed;
    for (std::uint32_t v = 1;; v = parent[v]) {
        reversed.push_back(v);
        if (v == 0) break;
        if (reversed.size() > n) throw Error(ErrorCode::internal, "cyclic route reconstruction");
    }
    out.path.assign(reversed.rbegin(), reversed.rend());
    for (const std::uint32_t v : out.path) out.native_node_ids.push_back(native_of[v]);
    out.found = true;
    return out;
}

} // namespace fbs::caves::detail
