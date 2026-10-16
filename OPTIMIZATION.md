# Scientific simulation memory optimization

This change reduces the memory required by KarstNSim's existing scientific
simulator. It preserves the CLI, geological inputs, sampling resolution,
water-table cost channels, directed routing and scientific model. The reference
is upstream commit `ff55d3969842f456908a6946c8bc8160a5d65965`.

## Storage changes

- **Packed adjacency:** contiguous 32-bit target IDs and edge-major float costs
  replace a 72-byte edge object and separately allocated per-edge vectors.
  Costs are calculated directly into their final storage. Deprecated fracture
  flags are removed from storage; fracture penalties remain part of the costs.
- **Compact sampling index:** each accepted Poisson point is stored once in a
  bucket instead of copied into many cell lists. Queries preserve the original
  covering-cell membership and strict distance test, with no extra random draws.
- **Smaller reverse index:** incoming arcs use a 32-bit flat index into the
  existing directed graph, with a 64-bit fallback. Forward/reverse edge order
  and live weight access are preserved.

For N points, K neighbor slots, C cost channels and E directed arcs, the main
graph payload is approximately `N*K*(4+4*C) + E*4 + (N+1)*8` bytes when 32-bit
reverse indices fit. Sampling, geology, path queues, surfaces and exports add
to this. Resolution, graph degree and the number of cost channels still govern
memory use; there is no fixed allocator cap.

The neighbor-graph exporter also skips unused slots before accessing geometry
and costs. Upstream dereferenced `samples[-1]` and misaligned edge properties
for these slots. The two affected regression hashes come from a separate
upstream build with only the documented one-line guard, rather than from the
optimized implementation. See [the regression corpus](KarstNSim/tests/regression/README.md)
for the retained original hashes, corrected reference patch and semantic checks.

## Full-resolution scientific example

The measured workload is unmodified `Input_files/1_base`, seed 1: 549,999
support points and 100 neighbor slots per point. The benchmark changes only
`main_repository` to route outputs into a new directory. No geological input
or sampling radius is changed.

Original baseline and a fresh run of this isolated optimization branch on
2026-09-24, Linux, Ryzen 7 5800X3D, GCC 13.3, CMake Release
(`-O3 -DNDEBUG`, C++17):

| Full example | Original | Optimized |
| --- | ---: | ---: |
| Whole-process peak RSS | 7,847.57 MB | 721.99 MB |
| Wall time | 41.65 s | 24.32 s |

MB is decimal. These are individual observations, not a statistical timing
benchmark. Peak RSS includes input loading, sampling, graph construction,
routing and enabled exports. All four exports match the original exactly:
sampling points, inception-surface points, cave network and connectivity matrix.
The cave-network SHA-256 remains
`c345582d7df7a1e78684f3ebb562499667f9a579f583b075cba468afe90fa096`.

Portable [baseline measurements](KarstNSim/tests/benchmarks/scientific-base-baseline.json),
[optimized measurements](KarstNSim/tests/benchmarks/scientific-base-optimized.json)
and [validation results](KarstNSim/tests/benchmarks/validation.json) are included.
Peak RSS fell by approximately 90.8% for this workload, without reducing its
scientific resolution or changing its inputs.

## Build and verify

```sh
cmake -S KarstNSim -B KarstNSim/build/release -DCMAKE_BUILD_TYPE=Release
cmake --build KarstNSim/build/release -j4
ctest --test-dir KarstNSim/build/release --output-on-failure

cmake -S KarstNSim -B KarstNSim/build/sanitize \
  -DCMAKE_BUILD_TYPE=Debug -DKARSTNSIM_SANITIZERS=ON
cmake --build KarstNSim/build/sanitize -j4
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  ctest --test-dir KarstNSim/build/sanitize --output-on-failure

python3 KarstNSim/tools/benchmark_memory.py \
  --binary KarstNSim/build/release/karstnsim \
  --out /tmp/karst-full \
  --compare KarstNSim/tests/benchmarks/scientific-base-baseline.json
```

The benchmark refuses existing output directories and records full-process RSS
and export hashes. It measures without a memory threshold by default. An
optional `--budget-mb` can enforce a limit selected by the researcher. `--compare`
requires the export names and hashes to match the supplied reference.

CTest covers packed storage, compact sampler membership, directed and surface
routing, and eight reduced scientific regression scenarios with 33 exports.
They exercise multiple water tables, perched springs, fractures, density
variation, noise, waypoints, cohesion, imported graphs, previous networks,
amplification, loops, dead ends and section simulation. The full example is
measured separately above. Expected exports use the original upstream result,
except for the two explicitly corrected neighbor-graph exports.

All four CTest entries passed in fresh Release and ASan/UBSan builds of this
branch. Leak detection was enabled and no sanitizer checks were disabled.

## Limits

The [additional workload measurements](SCALING.md) cover denser support graphs,
higher neighbor counts, multiple water tables and combined scientific features.
They include a direct feature-rich comparison against upstream with identical
exports and approximately 87% lower peak RSS, plus bounded graph-export tests.

The measured full example has one cost channel and disables amplification,
section simulation and full neighbor-graph export. Those modes remain available
and are exercised in smaller regression cases, but the full-example memory
result is not a bound for arbitrary scientific inputs. Ghost rocks and external
SGS drift are not covered by the supplied regression corpus.

The 64-bit reverse-index decoder is tested using small converted graphs; no
graph exceeding 2^32 slots was allocated. Existing compiler warnings remain.
Validation is on Linux with GCC 13.3; Windows/macOS builds and cross-compiler
byte-identical output have not been verified. The upstream MIT license and
scientific attribution remain unchanged.
