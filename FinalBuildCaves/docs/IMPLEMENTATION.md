# Standalone cave backend implementation plan

Goal: complete the generator-side logical/geometric map package approved in
Studio/Research/KarstNSim-2026-09-24/handoff/README.md. Another developer owns
MOOCoW. Its existing stored-data format and existing records must not change.

Architecture: an in-memory KarstNSim C++17 backend feeds an independent
FinalBuildCaves coordinator. The coordinator owns region boundaries, geological
route requests, logical cave geometry and authored constraints. A JSON/C ABI
and CLI expose owned results, plus an adapter to MOOCoW's current sectioned CSV.
No renderer, database, UI, dense global voxel volume or third-party game module
is required by the core. Detailed 3D wall meshes remain outside this delivery.

The public contract is include/fbs/caves.hpp. Meters, XY horizontal, +Z up,
region-local coordinates; double origins and signed 64-bit region keys.
Persistent IDs are UUID strings, with generation metadata in existing string
parameter maps or the separate manifest. Adjacent jobs use the same world
seed/settings and canonical shared-face keys. Each has a real cavern/port on
its own side; stitching creates a normal tunnel without database schema changes.

## Work lanes and acceptance

- [x] Native backend: owned job RNG/noise/logging, in-memory input/result,
  native node identity, cancellation/work/capacity checks, existing CLI parity.
  Worker owns KarstNSim/ only. New entry point is run_simulation_memory with
  RunOptions as recorded in sdk/agents/native-task.md.
- [x] Region coordinator: implement generate/validate/stitch/stable_id in
  src/region*.cpp and its native adapter. Canonical boundaries, deterministic
  chambers/floors/ports/profiles, deliberate loops, copied authored records,
  conservative exclusion routing. Worker owns these files and regional tests.
- [x] Transport/package: coordinator owns public headers, C ABI, JSON transport,
  existing MOOCoW CSV adapter, CLI, CMake install/export and consumer examples.
- [x] Acceptance: fresh installed C/C++ and .NET consumers; unchanged MOOCoW
  importer/validator on a temporary copy; concurrent seed isolation; adjacent
  seams; invalid/cancelled/budget failures; exact authored preservation; full
  scientific export regression; ASan/UBSan; process RSS and a multi-region sweep.
- [x] Independent review and corrections; pin final source, package archives,
  handoff, exact commands and known limits; publish no application changes.

Review focus: preserve coincident but distinct topology; reject invalid numeric
and size inputs before allocation; no partial accepted outputs on failure;
authored ports must be connected without altering authored geometry; unrelated
region generation order must not change IDs, seams or output. Reject infeasible
constraints explicitly. Never silently weaken clearance or slope settings.

Existing room/maze/progression/navigation libraries remain optional adapters.
Their rectilinear/tree limits do not replace natural cave or region topology.
Memory budgets remain measured workload targets, not universal malloc caps.

Execution continues under Max's explicit instruction to finish this side and
his existing autonomous-work preference. No additional permission gate is added.
