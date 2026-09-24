# KarstNSim memory optimization

Implemented 2026-09-24 on upstream `ff55d3969842f456908a6946c8bc8160a5d65965`, branch `perf/compact-cave-graphs`. C++17, GCC 13.3, Release `-O3 -DNDEBUG`, Linux on Ryzen 7 5800X3D. Three local Claude Opus 5.5 agents (`claude-opus-5-5`, verified in runner metadata) handled packed storage, sampler analysis, and independent regression references. The coordinator integrated and measured the final implementation.

## Measured result

Whole-process peak resident memory from GNU `/usr/bin/time -v`. MB and GB here are decimal; MiB and GiB are binary. These are complete command-line runs, including input loading, sampling, graph construction, routing, and the outputs enabled by the input file.

| Workload | Support points | Original peak | Optimized peak | Original wall time | Optimized wall time |
|---|---:|---:|---:|---:|---:|
| Unchanged `Input_files/1_base`, seed 1 | 549,999 | 7,847.57 MB (7.31 GiB) | **723.55 MB (690.03 MiB)** | 41.65 s | 18.72 s |
| Game preset, density radii multiplied by 3 | 24,278 | 359.19 MB | **44.80 MB (42.72 MiB)** | 1.66 s | 0.93 s |

The full example meets a **1,000 MB** budget with identical parameters. Its four exports (point cloud, inception-surface points, cave network, connectivity matrix) are byte-identical to upstream. The cave-network SHA-256 remains `c345582d7df7a1e78684f3ebb562499667f9a579f583b075cba468afe90fa096`.

The game preset meets **100 MB**. It uses the same domain, all the original geological constraints, and 100 neighbor slots per point; positive Poisson radii are multiplied by three. This yields a different, coarser network. All four exports match the original implementation when given that same game preset. This does not establish a 100 MB budget for the original 550,000-point graph, an entire game process, or a world containing millions of miles of tunnels.

Wall times are individual observations, not a statistical benchmark. The memory reductions are about 90.8% on the original example and 87.5% on the game preset. Timing and memory evidence is retained in the Studio research directory described below.

## What changed

- **Packed adjacency:** contiguous 32-bit target IDs and edge-major float costs replace a 72-byte object plus separately allocated vectors per edge. The arbitrary number of water-table cost channels is preserved. Deprecated fracture flags are removed from storage; the fracture penalty remains part of the cost. Costs are calculated directly into their final storage.
- **Compact sampling index:** each accepted Poisson point is stored once in a bucket instead of copied into hundreds or thousands of cell lists. The original covering-cell test and strict distance predicate remain intact. There are no extra random draws. The query can visit points in a different order because its sole result is a boolean rejection.
- **Smaller reverse index:** each incoming edge stores one 32-bit flat index into the existing directed graph, saving two bytes per arc. A 64-bit fallback preserves support for larger addressable graphs. Forward and reverse edge order and live weight access are preserved.
- **Graph-export bug fix:** unused slots are skipped before both geometry and cost export. Upstream attempted to read `samples[-1]` and then read past the cost arrays. ASan reproduced the error in unmodified upstream. The corrected exporter is compared to a separate upstream build containing only that one-line fix.

For N points, K allocated neighbor slots, C cost channels, and E actual directed arcs, the main graph payload is approximately `N*K*(4+4*C) + E*4 + (N+1)*8` bytes while 32-bit reverse indices fit. Sampling, geology, path queues, surfaces, and optional exports add to that. Input density, number of channels, and degree still determine memory use. The budget is measured, not a hard allocator cap.

## Build and verify

From the repository root:

```sh
cmake -S KarstNSim -B KarstNSim/build/release -DCMAKE_BUILD_TYPE=Release
cmake --build KarstNSim/build/release -j4
ctest --test-dir KarstNSim/build/release --output-on-failure

cmake -S KarstNSim -B KarstNSim/build/sanitize \
  -DCMAKE_BUILD_TYPE=Debug -DKARSTNSIM_SANITIZERS=ON
cmake --build KarstNSim/build/sanitize -j4
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  ctest --test-dir KarstNSim/build/sanitize --output-on-failure
```

CTest includes packed storage, compact sampler membership, directed and surface routing, and a Linux/Python regression suite of eight scenarios and 33 exports. The scenarios cover multiple water tables, a spring with no water table, fractures, variable and constant density, noise, waypoints, cohesion, imported graphs, previous networks, loops, dead ends, and section simulation. Original golden hashes are retained alongside the two corrected exporter hashes and the independent reference patch. See [the regression README](KarstNSim/tests/regression/README.md).

Final validation: all four CTest entries passed in both Release and ASan/UBSan builds. The sanitizer run enabled leak detection and halted on undefined behavior; all eight scenarios completed and all 33 exports matched the independent expectations. No sanitizer errors were reported.

The sampler agent also compared rejection decisions against the original index on over 36 million candidate queries across four configurations, with zero differences. Unit tests exercise the 64-bit reverse decoder with converted small graphs; a real graph exceeding 2^32 slots was not allocated.

Compiler warnings are enabled. The standalone packed-storage tests also pass with GCC and Clang under `-Werror`. The wider upstream code still emits existing warnings such as signed/unsigned comparisons and unused declarations; this is not a warning-clean rewrite.

## Reproduce the memory runs

Each output directory must be new. The helper generates the inputs, runs the binary, records full-process peak RSS and output hashes, and returns failure when the selected memory budget is exceeded.

```sh
python3 KarstNSim/tools/benchmark_memory.py \
  --binary KarstNSim/build/release/karstnsim \
  --case full --out /tmp/karst-full

python3 KarstNSim/tools/benchmark_memory.py \
  --binary KarstNSim/build/release/karstnsim \
  --case game --out /tmp/karst-game
```

`--compare /path/to/reference/metrics.json` additionally enforces export parity for the same workload. `--budget-mb 0` is intended only to measure an unoptimized reference without enforcing a memory threshold. Defaults are 1000 for full and 100 for game. The generated game input is `<out>/instructions.txt` with `<out>/game-box.txt`; it can be reused as a concrete preset.

## Scope and evidence

Durable measurements and reports live in the private Studio folder `Research/KarstNSim-2026-09-24/optimized/`: `full/`, `game/`, `game-baseline/`, agent reviews, test logs, and a validation summary. The preserved original full-run measurements are one directory above. The optimized executable SHA-256 for the reported Release measurements is `13cc6e26b7a055614aaa27a4e507c5cf85327476029ed4ce2bfafebc09efce51`.

No scientific features or CLI parameters were removed. The measured 1 GB result is for the original base configuration, which has one cost channel and disables amplification, section simulation, and full neighbor-graph export. Those features remain available and are exercised in smaller regression cases. Arbitrary higher-resolution inputs, many water tables, large graph exports, ghost rocks and external SGS drift are not covered by this memory guarantee. A bounded regional scheduler, persistent world representation, mesh generation and FinalBuildSystems integration remain separate work.

This is a local optimization branch, not an upstream release or a completed game-world streaming system. The upstream MIT license and scientific attribution remain in place.

## Standalone library follow-up

The subsequent [FinalBuildCaves SDK](FinalBuildCaves/README.md) adds in-memory,
independent native jobs and regional logical geometry with canonical shared
boundaries. The scientific CLI and its full parameter model remain available.
Its validation and measurements are recorded separately from the original
optimization measurements above. Application scheduling, storage policy,
rendered wall meshes and navigation integration remain application work.
