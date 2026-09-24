#pragma once

/**
@file library.h
@brief Reusable in-memory KarstNSim API (no filesystem access, no process-wide state).

run_simulation_memory() runs the same engine as the karstnsim executable on a ParamsSource
built in memory and returns one KarstNetworkResult per iteration. Each call owns its RNG,
simplex-noise permutation, logger, work counter and limits on its own stack frame, so
independent calls may run concurrently on different threads or be nested on one thread.
For the same parameters and seed, results do not depend on other calls.

Input rules (violations throw InvalidInputError before any sampling):
- Parameters are checked with validate_simulation_parameters(); karstic_network_name,
  save_repository and simulation_input_dir are not used and may be empty.
- Filesystem-only exports are rejected: create_vset_sampling, create_nghb_graph,
  create_nghb_graph_property and create_grid must be false. ParamsSource defaults
  create_vset_sampling to true, so in-memory callers must clear it.
- create_solved_connectivity_matrix is satisfied in memory through
  KarstNetworkResult::solved_connectivity_matrix (empty in sections-only mode, which does
  no routing and for which the legacy runner writes no matrix either).
- When use_user_connectivity_matrix is true (the ParamsSource default) and the run is not
  sections-only, ParamsSource::connectivity_matrix must hold one row per sink and one column
  per spring with values 0, 1 or 2; connectivity_matrix.txt is never read. Set
  use_user_connectivity_matrix to false (and leave the matrix empty) for automatic all-2
  connectivity.

Output: results[i] belongs to iteration i and holds the same segments, attribute values and
order that the legacy runner serializes to <name>_karst.txt, but with unrounded float
coordinates and with ResultPoint::node_id set to the native skeleton node index. Each
skeleton edge is listed once per direction; rebuild topology from node_id, never from
rounded or merged coordinates.

Limits (RunOptions) are checked at deterministic checkpoints, so for identical inputs a
limit is hit at the same point on every run:
- work_limit counts deterministic engine work units: one per Poisson candidate, surface
  triangle and refined sub-triangle, supplied sampling point, sample classified against
  water tables or measured against surfaces, neighbour-graph node and edge slot filled,
  imported edge, settled priority-queue entry in every routing, amplification and section
  distance search, tested amplification pair, dead-end origin, skeleton node visit and
  simulated section node. Linear scans in legacy deduplication loops add one unit per
  256 scanned items. The limit is exceeded when the running total passes work_limit.
- max_points bounds the sampling cloud (keypoints, generated or supplied points, surface
  points, previous-network points) before each insertion, and each refined surface point set.
- max_edges bounds the directed edge slots of the cost graph (points x neighbour slots, or
  points x maximum degree for an imported graph), checked before the adjacency storage is
  allocated, and the deduplicated imported edge list while it is built.
These are admission and progress limits, not an allocator byte cap. Memory also depends on
the number of water-table cost channels, input surfaces and grids, the Poisson bucket grid
(4 bytes per 64 background cells, within the engine's 500M-cell limit), surface refinement (bounded only through work_limit while a level is built),
skeleton size, section simulation and temporary search buffers; see OPTIMIZATION.md for the
measured graph payload formula.

Cancellation: is_cancelled(user) is polled at stage boundaries and every
JobContext::kPollInterval (4096) work units, from the calling thread.
When work_used is non-null, the final engine work counter is written on return
or exceptional exit. The pointer must remain valid until the call ends and must
not be shared between concurrent jobs without caller synchronization.

Logging: log == nullptr is quiet. Otherwise engine messages are written to *log from the
calling thread; share one stream between concurrent jobs only if it is synchronized. The
stream's format flags, precision and fill are restored when the call returns. std::cout and
std::cerr are never written or redirected by this API.

Errors (all derive from KarstNSim::Error, see errors.h):
- InvalidInputError: rejected parameters/options/connectivity, and failures while the
  network is being configured from the inputs.
- CancelledError: is_cancelled returned true.
- LimitExceededError with code WorkLimitExceeded, PointLimitExceeded or EdgeLimitExceeded.
- NoRouteError: an iteration produced no inlet-to-outlet path (no partial result vector is
  returned). Inlets that individually have no route are only logged, as in the legacy runner.
- Error with code SimulationFailed: any other engine exception after configuration; the
  message is the engine message.
std::bad_alloc propagates unchanged.
**/

#include <cstdint>
#include <ostream>
#include <vector>

#include "KarstNSim/errors.h"
#include "KarstNSim/models/results.h"
#include "KarstNSim/run_code.h"

namespace KarstNSim {

	struct RunOptions {
		uint64_t work_limit = 100000000;
		uint64_t max_points = 1000000;
		uint64_t max_edges = 100000000;
		bool (*is_cancelled)(void*) = nullptr;
		void* user = nullptr;
		std::ostream* log = nullptr; // null means quiet, NEVER redirect process cout
		std::uint64_t* work_used = nullptr; // optional output, including exceptional exit
	};

	std::vector<KarstNetworkResult> run_simulation_memory(ParamsSource parameters, const RunOptions& options = {});

}
