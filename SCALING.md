# Additional scientific workload measurements

These experiments assess higher-resolution support graphs and optional
scientific features of the existing simulator. The geological domain remains
the original `Input_files/1_base` volume. Denser cases multiply positive
Poisson radii by 0.8 or 0.6; all other density-grid properties remain unchanged.
They do not enlarge the physical domain.

## Method

Measurements were collected on 2026-09-24 using Linux x86_64, a Ryzen 7 5800X3D,
GCC 13.3 and CMake Release (`-O3 -DNDEBUG`, C++17). The optimized executable
SHA-256 is `13cc6e26b7a055614aaa27a4e507c5cf85327476029ed4ce2bfafebc09efce51`.
Its C++ implementation is the one in commit
`e88f7cf47fb5f7a9aa0a912096df0d331f33e354`. Later documentation and benchmark
harness changes do not change that executable.

Cases ran serially at nice level 5 on a shared machine. GNU `/usr/bin/time -v`
measured whole-process peak RSS and elapsed time; output hashing happened after
the measured process exited. MB is decimal. Timings are individual observations,
not a controlled statistical estimate of speedup.

Each test process had a 5 GiB virtual-address-space ceiling, a 4 GiB individual
output-file ceiling and a 600-second wall-time ceiling. A host-memory guard
would stop a case if available physical memory fell below 2 GiB. These are
limits of the measurement apparatus. A run stopped by them is incomplete;
it does not establish that the same workload fails without those limits.

The two-water-table cases assign three inlets to one outlet and two to the
other, using two edge-cost channels. The combined-feature cases additionally
enable waypoints, both supplied previous networks, 20 requested dead ends,
30 requested cycles, noise during amplification, and conduit-section
simulation. Original inception surfaces, fractures, karstification potential
and no-karst spheres remain enabled.

## Completed workload sweep

| Workload | Support points | Peak RSS | Elapsed time |
| --- | ---: | ---: | ---: |
| Original-resolution base | 549,999 | 722.64 MB | 19.01 s |
| Original resolution, 200 neighbors | 549,999 | 1,382.54 MB | 31.26 s |
| Original resolution, two water tables | 555,068 | 954.77 MB | 22.39 s |
| Original resolution, combined features | 555,333 | 979.03 MB | 63.02 s |
| Density radii x0.8 | 1,068,440 | 1,390.80 MB | 46.43 s |
| Density radii x0.6 | 2,519,423 | 3,187.09 MB | 129.15 s |
| Radii x0.8, 150 neighbors, combined features | 1,073,773 | 2,678.96 MB | 174.65 s |

Increasing point count from 549,999 to 1,068,440 increased peak RSS by about
1.92 times. The 2,519,423-point case used about 4.41 times the base memory for
4.58 times as many points. These observations are consistent with the packed
graph storage scaling with point count, neighbor count and cost channels.

The combined-feature full-resolution run spent 39.34 CPU seconds in
amplification, out of its roughly 63-second elapsed run. Conduit-section
simulation itself took 0.023 CPU seconds for that generated network. Feature
selection can increase runtime substantially without the same proportional
increase in peak memory.

## Direct comparison against upstream

The paired comparison uses all combined features above, with positive density
radii multiplied by 1.5 in both executables. This produces 171,571 support
points with 100 neighbors and two water-table cost channels.

| Implementation | Peak RSS | Elapsed time |
| --- | ---: | ---: |
| Unmodified upstream `ff55d39` | 2,461.90 MB | 35.63 s |
| Optimized implementation | 319.61 MB | 19.02 s |

All four exports are byte-identical between these runs: sampling points,
inception-surface points, final cave network including section properties, and
solved connectivity matrix. Peak RSS is approximately 87.0% lower. The unchanged
full-resolution base control also matches the previously recorded upstream
hashes for all four exports.

The larger optimized-only cases above establish completion and resource use.
They were not rerun against upstream and therefore do not establish additional
byte-for-byte equivalence at those larger sizes. Matching upstream output
checks numerical compatibility; it is not independent geological validation.

## Neighbor-graph export

The full-resolution cost-property graph export did not complete under the
5 GiB virtual-address-space ceiling. It threw `std::bad_alloc`, with
4,019.42 MB observed peak RSS after 31.95 seconds. Only the two point exports
were present, so this is recorded as an incomplete run, not a successful
simulation. The corresponding run without graph export completed normally.

A smaller, 72,785-point graph-export case did complete:

| Same scientific configuration | Peak RSS | Elapsed time |
| --- | ---: | ---: |
| Without graph export | 111.02 MB | 2.47 s |
| With graph and cost-property export | 2,704.63 MB | 293.05 s |

All four common outputs matched exactly. The additional neighbor-graph file
was 864,987,777 bytes. This demonstrates that optional graph export remains
a substantial memory and runtime cost even after the core graph is compacted.

Source inspection identifies existing work separate from packed graph storage:
`GraphOperations::save_nghb_graph` builds complete segment and nested property
buffers, and `Line::Line` uses repeated linear searches to discover unique
nodes. Those paths remain candidates for further memory and runtime work.
No export streaming or unique-node-index optimization is included here.

## Reproduction

Build the optimized executable as documented in `OPTIMIZATION.md`. Build a
reference at the pinned upstream revision, then run the supplied harness:

```sh
git worktree add --detach /tmp/karst-scientific-reference ff55d3969842f456908a6946c8bc8160a5d65965
cmake -S /tmp/karst-scientific-reference/KarstNSim \
  -B /tmp/karst-scientific-reference/KarstNSim/build/release \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/karst-scientific-reference/KarstNSim/build/release -j4

python3 KarstNSim/tools/benchmark_scaling.py \
  --source . \
  --binary KarstNSim/build/release/karstnsim \
  --reference /tmp/karst-scientific-reference/KarstNSim/build/release/karstnsim \
  --out /tmp/karst-scientific-scaling
```

The output directory must be new. It retains exact generated inputs, parameter
changes, process logs, time reports, export hashes and comparison results.
Portable measurements and case settings are included in
[`scaling-results.json`](KarstNSim/tests/benchmarks/scaling-results.json).

Ghost rocks and external SGS drift remain outside this corpus. Windows/macOS
behavior and cross-compiler byte-identical output have not been verified.
