#include "KarstNSim/graph.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
class TestGraph : public KarstNSim::CostGraph {
public:
    TestGraph() : CostGraph(6) { adj.reset(6, 3, 2); }
    using CostGraph::SetEdge;
    using CostGraph::UpdateEdgeWeight;
    using CostGraph::DijkstraComputePaths;
    using CostGraph::DijkstraComputePathsBidirectional;
    using CostGraph::DijkstraComputePathsSurface;
    using CostGraph::DijkstraComputePathsSurfaceBidirectional;
    using CostGraph::GetDirectedEdgeWeight;
};

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void compare_routes(const TestGraph& graph) {
    for (int channel = 0; channel < 2; ++channel) {
        for (int source = 0; source < 6; ++source) {
            for (int target = 0; target < 6; ++target) {
                std::vector<float> reference, packed;
                std::vector<int> prev_reference, prev_packed;
                graph.DijkstraComputePaths(channel, source, reference, prev_reference, target);
                graph.DijkstraComputePathsBidirectional(channel, source, target, packed, prev_packed);
                require(reference[target] == packed[target], "directed shortest-path cost changed");
                if (!std::isfinite(packed[target])) continue;
                int node = target;
                float path_cost = 0;
                int steps = 0;
                while (node != source) {
                    const int predecessor = prev_packed[node];
                    require(predecessor >= 0 && ++steps < 6, "invalid predecessor chain");
                    path_cost += graph.GetDirectedEdgeWeight(predecessor, node, channel);
                    node = predecessor;
                }
                require(path_cost == packed[target], "reconstructed route cost changed");
            }
        }
    }
}
}

int main() {
    try {
        using KarstNSim::CostGraph;
        require(CostGraph::ReverseEdgeIndexFitsUint32(65536, 65536), "uint32 last-index boundary rejected");
        require(!CostGraph::ReverseEdgeIndexFitsUint32(65536, 65537), "wide graph incorrectly uses uint32");
        require(CostGraph::ReverseEdgeIndexFitsUint32(0, std::numeric_limits<std::size_t>::max()), "empty shape rejected");
        require(!CostGraph::ReverseEdgeIndexFitsUint32(std::numeric_limits<std::size_t>::max(), 2), "width check overflowed");
        TestGraph graph;
        // Directed edges, asymmetric costs, distinct per-channel best routes,
        // an isolated node, and unoccupied slots all occur in imported graphs.
        graph.SetEdge(0, 0, 1, {1, 10});
        graph.SetEdge(0, 1, 2, {7, 1});
        graph.SetEdge(1, 0, 2, {1, 8});
        graph.SetEdge(1, 1, 3, {4, 2});
        graph.SetEdge(2, 0, 3, {1, 1});
        graph.SetEdge(3, 0, 4, {2, 2});
        graph.SetEdge(4, 0, 0, {9, 9});
        compare_routes(graph);

        // Exercise the real wide-index decoder without allocating a >2^32-slot graph.
        graph.reverse_edges64_.assign(graph.reverse_edges32_.begin(), graph.reverse_edges32_.end());
        std::vector<std::uint32_t>().swap(graph.reverse_edges32_);
        graph.reverse_edges_wide_ = true;
        compare_routes(graph);

        // Cohesion/ball-cost edits change weights without rebuilding topology.
        graph.adj(0, 1).weight[0] = 0.5f;
        compare_routes(graph);
        graph.UpdateEdgeWeight(1, 0, {20, 0.5f});
        compare_routes(graph);

        // SetEdge must invalidate the cached reverse topology.
        graph.SetEdge(0, 2, 4, {0.25f, 0.25f});
        compare_routes(graph);

        Array2D<char> surface(6, 2, 0);
        surface(2, 0) = surface(2, 1) = 1;
        surface(4, 0) = surface(4, 1) = 1;
        for (int channel = 0; channel < 2; ++channel) {
            for (int source = 0; source < 5; ++source) {
                std::vector<float> reference, packed;
                std::vector<int> prev_reference, prev_packed;
                int reach_reference = -1, reach_packed = -1;
                bool at_target_reference = false, at_target_packed = false;
                graph.DijkstraComputePathsSurface(channel, source, reach_reference, reference,
                    prev_reference, surface, 4, at_target_reference);
                graph.DijkstraComputePathsSurfaceBidirectional(channel, source, reach_packed, packed,
                    prev_packed, surface, 4, at_target_packed);
                require(reach_reference >= 0 && reach_packed >= 0, "surface was not reached");
                require(reference[reach_reference] == packed[reach_packed], "surface route cost changed");
                require(at_target_reference == at_target_packed, "surface/outlet classification changed");
            }
        }
        std::cout << "Directed routes, channels, mutations, and surface paths passed.\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
