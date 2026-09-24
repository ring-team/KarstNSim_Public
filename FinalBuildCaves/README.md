# FinalBuildCaves

Independent C++17 logical cave generation backed by the optimized KarstNSim
engine. It supplies regional cavern outlines, floors, transitions, ports and
tunnel profiles; native C++ and C APIs; JSON transport; a CLI; and an adapter to
MOOCoW's existing sectioned CSV format. No engine, editor, database or renderer
is a dependency. Detailed 3D wall meshes and navigation baking are outside this
logical-map package.

This is the generator side of the parallel MOOCoW effort. Application code,
stored-data format, database schema and existing records remain unchanged.
Generated bookkeeping is carried in region JSON and the completion manifest.
See `docs/MOOCOW-ADAPTER.md` for existing-format constraints and import policy.
See [validation and measured memory](docs/VALIDATION.md) for the tested workloads
and the limits of this delivery.

## Build

From the enclosing KarstNSim source root:

```sh
cmake -S FinalBuildCaves -B FinalBuildCaves/build/release -DCMAKE_BUILD_TYPE=Release
cmake --build FinalBuildCaves/build/release --parallel 4
ctest --test-dir FinalBuildCaves/build/release --output-on-failure
cmake --install FinalBuildCaves/build/release --prefix "$PWD/FinalBuildCaves/stage"
```

Use `-DBUILD_SHARED_LIBS=ON` for a dynamic library and the .NET P/Invoke consumer.
Use `-DFBS_CAVES_SANITIZERS=ON -DCMAKE_BUILD_TYPE=Debug` for ASan/UBSan. The legacy
scientific executable is built under `build/release/karst/karstnsim` with the
original instruction-file interface. The higher-level regional API is additive.
The native core and optional tools do not depend on Python or .NET. Existing
scientific regression tooling uses Python; the application compatibility test
uses an available .NET SDK in an isolated temporary copy.

## Generate an adjacent pair

The destination directory must not exist. Use a fresh map UUID when importing
into an application: MOOCoW's existing importer can replace a map with the same
ID. This command creates files only and never opens a database.

```sh
FinalBuildCaves/build/release/fbs-caves \
  --request FinalBuildCaves/examples/region.json \
  --request FinalBuildCaves/examples/neighbor.json \
  --map-id e8cf283d-6cb9-49e6-8ad4-2325f0ee7519 \
  --out /tmp/cave-pair
```

The output contains `region-0.json`, `region-1.json`, `dataset.csv` and
`manifest.json`. Only accept a CLI output directory with a completion manifest;
it is written last through a temporary file and atomic rename. CLI failure removes its newly created directory. Existing
destinations are rejected. SIGINT/SIGTERM requests cooperative cancellation.

## Native consumer

```cmake
find_package(FinalBuildCaves 1 CONFIG REQUIRED)
target_link_libraries(my_generator PRIVATE fbs::caves)
```

```cpp
#include <fbs/caves.hpp>
fbs::caves::Request request;
request.world_seed = 42;
request.key = {0, 0, 0};
auto region = fbs::caves::generate(request);
auto json = fbs::caves::to_json(region);
```

See `examples/installed-consumer` for separately compiled C and C++ consumers,
and `examples/dotnet` for P/Invoke without a MOOCoW dependency. Scientific users
may link `fbs::karst` and include `KarstNSim/library.h` for the lower-level
in-memory API and full geological parameter model.

## Contract

- Public schema and ABI version are 1; algorithm identity is `fbs-caves-1`.
  Unknown versions are rejected. JSON request keys and enums are strict.
- Coordinates are local meters, XY horizontal and +Z vertical. Region keys use
  signed 64-bit integers; origin is kept separately in double precision. The
  regional validator bounds supported values before float conversion.
- Generation uses the world seed, region key and settings. Independent calls
  own their state. Reproducibility is tested on the pinned native build; binary
  equality across different compilers/standard libraries is not promised.
- Adjacent regions share canonical face keys and seam geometry. Each owns its
  own cavern/port; `stitch` produces ordinary cross-region tunnel records.
  An incompatible seed, configuration or boundary is an error.
- Native skeleton IDs remain provenance. Persistent UUIDs are generated from
  domain-separated identities, not rounded coordinates. Projected line
  intersections do not implicitly join tunnels.
- Authored records and metadata are preserved, and required connection ports
  and reserved space constrain new content. An infeasible request returns an
  error; it must not be treated as a partially successful cave.
- Work, point and edge limits apply per region. The work counter covers geometry,
  validation and every native route attempt together. These are not
  a universal byte allocator cap. Memory measurements cover the named test
  workload. The application controls concurrent jobs, storage and resident
  regions for large worlds.

The existing room, maze, progression and navigation modules remain optional.
Their rectilinear/tree constraints are not imposed on natural cave topology.

`caves_profile 64` exercises sequential generation and stitching while retaining
only two regions. Use an external process RSS tool for memory measurements.
`tools/plot_regions.py region-0.json region-1.json --out preview.png` plots actual
logical geometry and requires matplotlib only for that optional visualization.

## Ownership and C ABI

`fbs_caves_generate_json` accepts UTF-8 and returns an owned opaque result.
Independent jobs can run concurrently. Copy/view its JSON, then call
`fbs_caves_result_destroy`. Failed calls leave the output handle unchanged.
The optional cancellation callback runs synchronously on the calling thread.
Do not throw across the C boundary, mutate a request during a call, destroy a
result while another caller reads it, or mix allocation/free functions.

The C API catches C++ exceptions and reports status plus an optional error
message. JSON input is limited to 16 MiB and 64 nesting levels. An export call
accepts at most 1024 regions; stream larger worlds as separate application jobs.
`fbs_caves_export_moocow_csv` returns an owned text buffer, never modifies a DB.

## Sources and licenses

The new SDK code is MIT licensed. KarstNSim is copyright Université de Lorraine,
ANDRA and BRGM under MIT; its original notice is retained. Its authors request
citation of Gouy et al., 2024, Journal of Hydrology,
DOI `10.1016/j.jhydrol.2024.130878`. The vendored JSON parser is nlohmann/json
3.11.3 under MIT; `vendor/nlohmann/SOURCE.json` records exact source URLs and
SHA-256 hashes. No licensed game art or Studio private account data is included.
