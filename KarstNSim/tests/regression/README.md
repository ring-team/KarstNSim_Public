# KarstNSim reduced regression corpus

This corpus checks that a changed `karstnsim` binary still produces byte-identical exports to the upstream baseline (commit `ff55d3969842f456908a6946c8bc8160a5d65965`, "Patch 2.1") on small, reproducible inputs.

**These are reduced fixtures, not the original examples.** Most cases use a coarsened copy of `Input_files/1_base/example_box.txt` and smaller `nghb_count` values, so they run in seconds. The golden hashes apply only to these reduced configurations. Checking the full-size `1_base` example is a separate step. For reference, its `base_0_karst.txt` SHA-256 is `c345582d7df7a1e78684f3ebb562499667f9a579f583b075cba468afe90fa096`.

## Usage

```sh
# Compare a candidate that includes the save_nghb_graph fix (see "Golden files") (exit status 1 on any difference)
python3 run_regression.py --binary /path/to/karstnsim --out /tmp/karst-rg-candidate \
    --golden golden/expected_optimized.json \
    [--reference-out /tmp/karst-opus-validation/baseline-run1]   # optional: line-level diff summaries

# Re-record the golden hashes (baseline binary only)
python3 run_regression.py --binary /path/to/baseline/karstnsim --out /tmp/karst-rg-base \
    --write-golden golden/baseline_ff55d39.json
```

`--cases c01_base_coarse c05_import_graph` runs a subset. `--assets` overrides the `Input_files` location, which defaults to `../../../Input_files` relative to this directory. Only the Python 3 standard library is required. When `/usr/bin/time` is available, the runner uses it to record peak RSS.

The run writes `<out>/result.json`, which holds the exit code, wall/user/sys time, peak RSS and the SHA-256, byte count and line count of each export. Each case's work root also keeps `stdout.txt` and `time.txt`.

## How it works

- `make_fixtures.py` derives the reduced inputs into `<out>/_gen/`. It is deterministic, uses a fixed-seed PRNG and runs in about 2 s. Its output hashes are recorded in the golden file. The derived files are:
  - `box_coarse.txt`: every positive density radius ×3 (0.03→0.09, 0.007→0.021). NDV cells and IKP are unchanged, so the density stays variable.
  - `box_medium.txt`: the same with a ×2 factor.
  - `box_gradient.txt`: positive density replaced by a smooth u/w gradient between 0.05 and 0.12, so the sampler sees many distinct local radii.
  - `graph_pts.txt` and `input_graph.txt`: a jittered 40 m lattice below the topography and outside the no-karst spheres, plus a symmetric 12-NN graph in `NODES/EDGES` format. The graph's nodes are exactly the sinks, springs, waypoints and lattice points.
- Each case directory `cases/<name>/` contains `instructions.txt` and, when needed, `connectivity_matrix.txt`. The runner copies them into a fresh work root, `<out>/<name>/`, and links three directories there: `assets` → `Input_files`, `fixtures` → `fixtures/` and `gen` → `<out>/_gen`. The binary runs with that work root as cwd, reading `main_repository: .`. KarstNSim reads path values as single whitespace-delimited tokens, so the instruction files never contain the repository path, which has spaces in it.
- `fixtures/` contains two small spring files that keep the original coordinates but change `surfindex`. In `springs_wt12.txt`, spring 1 uses water table 1 and spring 2 uses water table 2. In `springs_wt1_perched.txt`, spring 1 uses water table 1 and spring 2 has `surfindex` 0, meaning no water table.

## What is compared

Every file in `<work>/outputs/` is compared byte for byte by SHA-256, and a missing or extra file counts as a failure. Exit codes are compared too. Two files are excluded:

- `karstnsim_console_<timestamp>.log`, which contains wall-clock timings and absolute paths.
- `simulation_times.txt`, which is written to the work root rather than `outputs/`.

The compared exports contain no timing values or timestamps. They are:

| Export | Content |
| --- | --- |
| `*_karst.txt` | Final network: route nodes, cost, equivalent radius, branch id |
| `*_pts.txt` | Full sampling cloud; with noise enabled it also carries a noise property |
| `*_s_pts.txt` | Inception-surface sampling |
| `*_connectivity_matrix.txt` | Solved inlet/outlet associations |
| `*_nghb_graph.txt` | Exported neighbour graph (c05 also has per-edge cost properties). The upstream export is buggy; see "Golden files" |

The exported graph's segment order follows the internal adjacency order. A change to the graph storage that keeps the same edge set and costs, but visits neighbours in a different order, will therefore fail this byte check even when routes are unchanged. Treat an `_nghb_graph.txt`-only mismatch as needing semantic review with `check_nghb_export.py` or by comparing sorted undirected edges and costs. It should not be read as an automatic route regression. Route, point, matrix and section exports should match exactly: the code is single-threaded, and two baseline runs were byte-identical.

The golden hashes were recorded with a baseline binary built by GCC 13.3.0 using CMake `Release` (`-O3 -DNDEBUG`), binary SHA-256 `da488944…b255`. Build the candidate with the same compiler and floating-point flags. `-ffast-math`, `-march=native`, or a different FP-contraction setting can change the last bits of float output even when the logic is unchanged.

## Golden files

- `golden/baseline_ff55d39.json` records the original upstream exports, unchanged since it was first written. Use it to check a binary that should reproduce ff55d39 exactly, bugs included.
- `golden/expected_optimized.json` is identical except for the two neighbour-graph exports listed under `bug_fix_overrides`: `c05_import_graph/c05_0_nghb_graph.txt` and `c07_const_density_loops_deadends/c07_10_nghb_graph.txt`. **These two hashes are a bug fix. They are not an identical original export.**

### The upstream bug these overrides fix

`GraphOperations::save_nghb_graph` in ff55d39 has two problems:

- **Out-of-bounds read.** It iterates every slot of the rectangular adjacency and appends a segment even for unused slots, where `target == -1`. Such a segment's end point is `samples[-1]`, an out-of-bounds read.
- **Misaligned costs.** The cost properties are appended only for valid targets. Every cost after the first unused slot is therefore attached to the wrong segment, and the property array is read past its end.

The original export therefore contains junk segments whose end point is heap metadata. Examples are `0.000 0.000 0.000` and `538.738 0.000 0.000`. Their bytes depend on the heap layout: an earlier trial baseline run of c05 produced 10597018 bytes, while the golden run produced 10597064. The trial files were not kept.

AddressSanitizer on unmodified ff55d39 reports this on c05 and c07:

```
heap-buffer-overflow ... READ of size 12 ... located 12 bytes before <samples> region
SUMMARY: AddressSanitizer: heap-buffer-overflow graph_operations.cpp:591 in GraphOperations::save_nghb_graph
```

### How the corrected hashes were derived

- **Independent reference build.** The corrected hashes come from a separate build, not from the optimized binary: ff55d39 plus only `golden/save_nghb_graph_unused_slot_guard.patch`. That patch adds one line, `if (adj[i][j].target < 0) continue;` at the top of the `j` loop.
- **Build configuration.** The reference was built with GCC 13.3 in CMake `Release` (`-O3 -DNDEBUG`). The reference binary's SHA-256 is recorded in the JSON.
- **Other exports unchanged.** In c05 and c07 on the reference build, all other exports stayed byte-identical to `baseline_ff55d39.json`.
- **Semantic check.** `check_nghb_export.py ORIGINAL CORRECTED PTS [--edges INPUT_GRAPH]` confirmed three things:
  - the corrected geometry equals the original geometry with its junk segments removed (c05: 93020 → 62366 segments, exactly 2 × 31183 imported edges; c07: 195780 → 189254);
  - every segment end point lies in the sampling cloud;
  - the c05 cost sequence is the original sequence, no longer shifted. In the original, 62353 of the 62366 valid segments carried a wrong cost.
- **Sanitizer check.** An ASan/UBSan build of the patched reference ran c05 and c07 with no report and gave the same outputs as the Release reference.

## Cases

| Case | Main coverage |
| --- | --- |
| c01_base_coarse | `1_base` settings on the ×3 variable-density box; 2 iterations with `vary_seed` (graph rebuilt in one process); fractures; inception surfaces; IKP; no-karst spheres; user matrix; no water table |
| c02_gradient_noise_waypoints | Gradient density; noise on the whole simulation (2 octaves); waypoints; generated all-`2` matrix with closest-outlet selection, `gradient_constraint_weight` 0.5 and `outlet_selection_cost_factor` 1.1; two water tables with one perched spring (`surfindex` 0); cohesion only in phreatic areas |
| c03_multi_wt_random_outlet | Two water tables, one per spring; mixed 0/1/2 user matrix with random outlet selection; `use_max_nghb_radius`; `multiply_costs`; vertical stretching ×3 |
| c04_polyphasic_amplified_sections | Previous networks (`3_amplification` base + polyphasic); waypoints; dead-end amplification (20) and cycle amplification (30 loops); noise only during amplification; SGS section simulation; two water tables |
| c05_import_graph | Imported sampling points and imported neighbour graph; neighbour-graph export with cost properties; waypoints |
| c06_sections_only | Example 4 (`sections_simulation_only`) on the original `amplification_0_karst.txt`. Not reduced; it is cheap. It does not load the example's `base_0_pts.txt`, which sections-only mode does not use |
| c07_const_density_loops_deadends | Constant-density Poisson sampling (`poisson_radius` 0.05); dead-ends and cycles without previous networks; sections with Spherical/Gaussian variogram models; neighbour-graph export |
| c08_base_medium | `1_base` settings with the original `nghb_count` 100 on the ×2 box: about 1 GB peak RSS and 5 s. This is the largest graph in the corpus and gives a memory signal |

Not covered: ghost rocks (the examples have no alteration lines), external SGS drift (the examples have no radius observations), `create_grid` and the full-size configurations.
