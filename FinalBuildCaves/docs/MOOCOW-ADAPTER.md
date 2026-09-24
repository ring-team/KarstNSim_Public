# MOOCoW sectioned-CSV adapter

```cpp
std::string fbs::caves::to_moocow_csv(const std::vector<Region>& regions,
                                      const std::string& map_id,
                                      const std::string& map_name = "Generated caves");
```

`to_moocow_csv` turns one or more generated regions into a single new map in
MOOCoW's existing sectioned CSV interchange format. It targets the format that
`CsvImport.ImportDataset` reads and `CsvExport.ExportDataset` writes, as of
MOOCoW commit `eb3318a`. It adds no section, column, table, migration or record
kind. MOOCoW is not modified. The adapter has no .NET or JSON-library dependency
and builds as C++17.

The call either returns the complete CSV text or throws `fbs::caves::Error`.
It never returns partial output. It writes no files and never reads or writes
a database. Before returning, it checks every rule that MOOCoW's CSV parser,
domain constructors and `CaveDatasetValidator` apply. Returned text therefore
imports and validates in unmodified MOOCoW, and MOOCoW does not change any value
while importing it. If a value cannot pass through MOOCoW unchanged, the adapter
rejects it and explains why instead of changing it.

## Using the output as a new dataset

- `map_id` must be a canonical, non-nil UUID that the target database does not
  already use. MOOCoW's `import` command replaces an existing map that has the
  same ID. The adapter cannot see the database, so the caller must check.
  `MapRepository.GetById` and `GetByCode` should both return null.
- The map code is derived as `fbs-caves-<map_id>`. `maps.code` is unique across
  the whole database (NOCASE), and a derived code cannot collide with another map.
- Record IDs are primary keys across the whole MOOCoW database, not per map.
  Importing the same generated world, with the same `world_seed` and settings,
  as a second new map collides on these keys. MOOCoW then rolls back the whole
  import. This is verified in the acceptance run. In one database, one world
  maps to one MOOCoW map.
- The CSV equals MOOCoW's own export of the stored map. MOOCoW's revision
  `content_hash` is therefore reproducible from the package bytes:
  `sha256(csv + "\nrevision=1")`. The acceptance run verifies this.

## Format written

The sections and their order are exactly the ones `CsvExport` writes:
`MAP`, `ELEVATION_BANDS`, `CAVERNS`, `CAVERN_OUTLINE_POINTS`, `FLOOR_REGIONS`,
`FLOOR_REGION_POINTS`, `CLIFF_EDGES`, `CLIFF_EDGE_POINTS`, `CAVERN_PORTS`,
`TUNNELS` and `TUNNEL_POINTS`. The adapter uses the same header lines, the same
column order and a single blank line between sections. Lines end in LF. The file
has no BOM, and text is UTF-8.

- **Numbers.** `CsvText.Number` uses `ToString("G17", InvariantCulture)`. The
  adapter writes the equivalent `%.17g` with an upper-case `E`, for example
  `1E+20`, `1.0000000000000001E-05` or `-0`. The values round-trip bit-exactly.
- **Integers.** Generation and variation seeds are written as signed decimal
  int32 values. Point indexes start at 0 and are contiguous.
- **Text cells.** These follow `CsvText.Escape`: a cell containing `,`, `"`,
  CR or LF is quoted, and embedded `"` is doubled. Empty optional IDs are
  written as empty cells.
- **Tags.** Tags are a JSON array of strings.
- **`generation_json`.** This is a JSON object with string values. JSON strings
  use the escaping of System.Text.Json's default encoder: printable ASCII stays
  literal except `" & ' + < > \``, which become upper-case `\uXXXX` escapes.
  `\\` and the short escapes `\b \t \n \f \r` are kept, and non-ASCII
  characters become UTF-16 `\uXXXX` escapes.
- **Order.** Order follows MOOCoW's storage order in
  `CaveDatasetRepository.LoadDatasetWithRevision`. Caverns and tunnels are
  sorted by ID. Floors, cliffs and ports are grouped by cavern, in cavern ID
  order, and sorted by ID within each group. Bands are sorted by display order.
  Generation parameters are sorted by ordinal UTF-16 order. Tags are sorted
  with `OrdinalIgnoreCase`. IDs use the lower-case `D` form.
- **Region order.** Input order does not affect the output. Regions are sorted
  by key in `(z, y, x)` order.

## Field mapping

In this table, `o` is the region's `origin`.

| CSV section.column | Source |
| --- | --- |
| MAP.id / code / name | `map_id` in lower-case form / `fbs-caves-<map_id>` / `map_name` |
| MAP.min/max xyz | Union of every region box `[o, o+size]` and all written geometry, including tunnel profiles ± width/2 and height/2, port openings, floor vertices at their sloped elevation ± variation, cavern bounds, outlines and cliff paths |
| ELEVATION_BANDS | One band per region Z layer (`key.z`), described below |
| CAVERNS.* | `Cavern`: `center + o`, `bounds + o`, `authoring`, `locked` (1/0), `reserved_clearance`, `generation.seed`, `tags`, `generation.parameters` |
| CAVERN_OUTLINE_POINTS | `Cavern::boundary + o.xy`. Written only when the boundary is non-empty |
| FLOOR_REGIONS.* | `Floor`: `base_z + o.z`, and `slope_x`, `slope_y`, `variation_amplitude`, `variation_seed`, `material`, `traversal` and `tags`, all unchanged |
| FLOOR_REGION_POINTS | `Floor::boundary + o.xy` |
| CLIFF_EDGES.* | `Cliff`: `from_floor_id`, `to_floor_id` (empty means a void edge), `code`, `transition`, `height` |
| CLIFF_EDGE_POINTS | `Cliff::boundary + o.xy` |
| CAVERN_PORTS.* | `Port`: `floor_id` (may be empty), `position + o`, `facing` (unchanged), `width`, `height`, `state`, `traversal` |
| TUNNELS.* | `Tunnel`: `map_id`, ports, `maximum_slope_degrees`, `minimum_clearance`, `traversal`, `locked`, `generation`, `tags` |
| TUNNEL_POINTS | `Profile`: `center + o`, and `width` and `height` for each profile, so variable profiles are kept |

Enumerations are written in the lower-case snake-case form MOOCoW parses, for
example `constrained_procedural`, `climbable_cliff` or `sealed`.

These values are not written to the CSV. They stay in the Region JSON, which
the coordinator's transport produces:

- `Region::version`, `algorithm`, `id`, `settings_id`, `world_seed`, `key`,
  `size`, `origin`
- `boundaries`
- `work_used`, `support_points`
- `Tunnel::source_node_ids`

No database field is invented for them. Band `style_json` is not part of the
CSV, so it keeps its database default of `{}`.

### Elevation bands

Each distinct `key.z` becomes one band:

| Field | Value |
| --- | --- |
| Code | `fbs-z<key.z>` |
| Name | `Region layer z=<key.z>` |
| Display order | 0, 1, … from bottom to top |
| ID | `stable_id("fbs.caves.moocow.elevation_band\n<map_id>\n<key.z>")` |

The band range starts at `[origin.z, origin.z + size.z]`. The lowest band then
extends down to the map's minimum Z and the highest band extends up to the map's
maximum Z. Any gap between layers is closed by raising the lower band's maximum.
Together, the bands cover the full map Z range. Neighbouring bands touch but
never overlap, which MOOCoW's `ELEVATION_BAND_OVERLAP` rule allows.

## Coordinates and precision

- Both formats use meters, with X and Y horizontal and +Z up, so no axis is
  swapped. Every point is translated with `world = origin + local`. Directions
  (`facing`), sizes, slopes, gradients and angles are not translated.
- Tunnels returned by `stitch(first, second)` are local to `first`, so the
  adapter offsets them by `first.origin`.
- **Floor slope convention.** `slope_x` and `slope_y` are gradients in m/m,
  anchored at the floor's first boundary vertex:
  `z(x, y) = base_z + slope_x·(x − x0) + slope_y·(y − y0)`. A translation
  therefore changes only `base_z`. MOOCoW itself does not define an anchor.
  Generated floors are currently flat.
- **Exactness.** The written world value is the correctly rounded sum
  `fl(origin + local)`. G17 then reproduces that double bit-exactly in MOOCoW.
  Two cases apply when recovering a local value:
  - When `local` and `origin` are exactly representable at the world value's
    precision, for example dyadic authored values with integer origins, the
    subtraction `world − origin` gives back the exact original local value. The
    constructed test checks this bit for bit.
  - Otherwise, the recovered value is within half an ULP of the world value.
    For example, at |world| < 2^20 m that is ≤ 1.2e-10 m.

  The Region JSON keeps the local values exactly.
- Regions exported together must lie on one lattice: `origin − key·size` must
  be equal for all of them, within 1e-9 plus 8 ULP. They must also have equal
  `size`, `version`, `algorithm`, `world_seed` and `settings_id`.

## Multi-region export

1. Regions are sorted by key. A duplicate key or a repeated region ID is
   rejected. Each region must pass `validate()`.
2. For every face-adjacent pair `(a, b)` with `b.key = a.key + e_axis`, the
   adapter checks whether either region has a `Boundary` on the shared face
   (`a` on the positive face, `b` on the negative face). If so, it calls
   `stitch(a, b)` with the lower key first, so the result is deterministic, and
   adds the ordinary tunnels it returns.
3. Every boundary port on that shared face must be the start or end of a
   returned tunnel. Otherwise the export fails with `incompatible_boundary`.
   Errors from `stitch`, such as a different seed or settings, or a one-sided or
   inconsistent seam, are passed through unchanged.
4. Duplicate entity IDs anywhere in the export are rejected, never merged. This
   covers regions, stitched tunnels, bands and `map_id`. So are duplicate codes
   in a MOOCoW scope, compared ignoring case.

## Rejections

`Error::code` is set as follows:

- **`invalid_input`**: malformed values.
- **`incompatible_version`**: a version or algorithm mismatch.
- **`incompatible_boundary`**: inconsistent regions or seams.
- **`constraint`**: valid SDK data that the existing MOOCoW format or validator
  cannot hold without changing it.

Violations of MOOCoW validator rules are collected and reported together, up to
25 lines.

### Identity and input

- `map_id`, and every record or reference ID, must be a canonical
  8-4-4-4-12 hex UUID and must not be nil. Braced and hyphenless forms are
  rejected.
- The input must not be empty, and there must be at least one cavern in total.

### Text and numbers

- Every number in every record must be finite, and it must stay finite after
  the offset is added.
- Code, name, material and `map_name`:
  - must be valid UTF-8 and must not be empty;
  - must not have leading or trailing .NET whitespace, because
    `DomainGuard.RequiredText` trims it; this includes U+00A0 and U+2000–U+200A;
  - must not contain CR or LF, because `CsvImport` splits lines before it reads
    quotes.

### Tags and generation parameters

- **Tags** follow the same trimming rule. Tags that differ only in case are
  rejected, because MOOCoW would silently merge them. Non-ASCII tags are also
  rejected (see the gaps below).
- **Generation parameter keys** follow the trimming rule. A key named exactly
  `parameters` is rejected:
  - `CsvText.ParseGenerationParameters` and `Cavern/TunnelRepository` treat a
    root property named `parameters` as a legacy wrapper.
  - The plain form fails `CsvImport`. The wrapped form
    `{"parameters":{...}}` imports and saves, but the map can then no longer be
    loaded. The acceptance probe verifies this.
  - Compatible mapping: rename the key, for example to `fbs.parameters`.

### Port facing

- A facing must be a bit-exact fixed point of MOOCoW's `UnitVector3`
  normalization: `x/m` with `m = sqrt(x*x + y*y + z*z)` and no fused
  multiply-add. MOOCoW would otherwise store a rescaled value.
- For horizontal directions, a single identical normalization step always
  reached a fixed point on 2M sampled angles. The regional lane was given this
  recipe in `region-review.md`.

### Records and geometry

The adapter mirrors the domain constructors:

| Rule | Detail |
| --- | --- |
| Bounds | Minimum ≤ maximum, with finite extents |
| Cavern center | Inside the cavern bounds |
| Clearance, amplitude, cliff height | Must not be negative |
| Polygons | At least 3 distinct vertices |
| Cliff paths | At least 2 points; a cliff must not go from a floor to itself |
| Port width and height | Greater than 0 |
| Tunnel ports | Two different ports |
| Tunnel slope | `0 ≤ maximum_slope_degrees < 90` |
| Tunnel clearance | Greater than 0 |
| Tunnel profiles | At least 2; every profile height ≥ clearance and width > 0 |

### `CaveDatasetValidator` rules

These rules are checked in world coordinates, with the same arithmetic MOOCoW
uses:

- **IDs.** Duplicate IDs are not allowed across all record kinds.
- **Codes.** Codes must be unique, ignoring case, within each scope:
  - map scope covers bands, caverns and tunnels together;
  - cavern scope covers floors, cliffs and ports together.
- **References.** Each reference must point to a record that exists and has the
  right owner.
- **Polygons.**
  - Consecutive points must not repeat.
  - Edges must not cross (with a 1e-9 tolerance on orientation).
  - Outlines and floors must lie inside the cavern's XY bounds.
- **Ports.**
  - A port must be inside its cavern's 3D bounds.
  - A port with a floor must be inside that floor's polygon.
  - A required port must be used by a tunnel. A sealed port must not be used.
- **Tunnel shape.**
  - Each end of a tunnel must be within `max(width, height)` of its port.
  - Each segment's slope must be ≤ the tunnel's maximum + 1e-9.
- **Transitions.** A `one_way_drop` must have a destination floor. The
  transitions `slope`, `ramp`, `stairs`, `climbable_cliff`, `impassable_cliff`
  and `one_way_drop` need a height greater than 1e-9.
- **Connectivity.**
  - Each cavern needs a port that a tunnel uses, or an `isolated` tag.
  - The floors of a cavern must be connected through transitions that are not
    `unconnected` or `impassable_cliff`. A floor outside the main component
    needs an `isolated` tag.
- **Exact caverns.** An `exact` cavern needs an outline.

## Known gaps

| Gap | Cause | Proposed compatible mapping |
| --- | --- | --- |
| Tag case folding outside ASCII | MOOCoW deduplicates tags with .NET `OrdinalIgnoreCase`, which follows Unicode simple case mapping. The adapter has no Unicode tables, so it rejects non-ASCII tags rather than risk MOOCoW silently merging them. | Use ASCII tags, or keep non-ASCII labels in names or generation parameters. |
| Code case folding outside ASCII | The adapter's duplicate-code check folds ASCII only. For non-ASCII codes, a collision that differs only in case is caught by MOOCoW's validator instead. That check is loud, not silent. | Use ASCII codes. |
| Generation key `parameters` | This is a MOOCoW defect: its own editor and exporter accept the key, but the repository cannot reload it. | Rename the key. On the MOOCoW side, only unwrap `parameters` when it is a JSON object. That fix belongs to the MOOCoW developer; the adapter does not change it. |
| Precision of the translated coordinates | This is inherent to storing world doubles. | The exact local values stay in the Region JSON. |
| Replacing existing maps and custom input | These are out of scope. MOOCoW's A04-P4 defines replacement; this adapter creates new maps only. | — |

## Tests

These are the owned files:

| File | Purpose |
| --- | --- |
| `tests/adapter_tests.cpp` | Native adapter tests. The coordinator's CMake builds them as `caves_adapter_test`, linked against the real library. |
| `tests/adapter_region_stub.cpp` | Test-only `validate`, `stitch` and `stable_id` for a standalone build. Never link it together with `src/region*.cpp`. |
| `tests/adapter_acceptance.py` | Driver for the unmodified-MOOCoW acceptance run. It is optional: it exits 77 if `dotnet` or the MOOCoW checkout is missing. |
| `tests/adapter_acceptance.cs`, `tests/adapter_acceptance.csproj` | The .NET consumer compiled against a temporary copy of MOOCoW's own Domain and Infrastructure projects. |

```sh
# Native, real library (coordinator CMake):
ctest -R caves_adapter
# Native, standalone with the stub (no region library needed):
g++ -std=c++17 -IFinalBuildCaves/include -DFBS_ADAPTER_STUB_REGION \
    FinalBuildCaves/src/moocow_csv.cpp FinalBuildCaves/tests/adapter_region_stub.cpp \
    FinalBuildCaves/tests/adapter_tests.cpp -o adapter_stub_test && ./adapter_stub_test
# Unmodified MOOCoW acceptance (.NET optional; everything in a temp dir):
python3 FinalBuildCaves/tests/adapter_acceptance.py --standalone
python3 FinalBuildCaves/tests/adapter_acceptance.py --adapter-test build/caves_adapter_test
```

### What the acceptance run does

1. It copies MOOCoW's `Code/src/MOOCoW.Domain`, `Code/src/MOOCoW.Infrastructure`
   and `Directory.Build.props` from a commit with `git archive`. This is
   read-only. The checkout's `global.json` is not used.
2. It restores NuGet packages from the local package cache only, and puts all
   build output in the temporary directory.
3. For each fixture, it runs:
   - `CsvImport`, then `CaveDatasetValidator`, which must report 0 issues;
   - `CsvExport`, whose output must be byte-identical to the fixture;
   - `CaveDatasetRepository.ImportDatasetWithBackup` into a throwaway SQLite
     file that already holds MOOCoW's seed map;
   - a reload, whose export must be byte-identical;
   - a check of the revision `content_hash`.
4. It then checks that the seed map, `sqlite_master` and `schema_migrations` did
   not change, and that re-importing records whose IDs collide rolls back
   completely.
5. Probes show what MOOCoW does to values the adapter rejects:
   - it trims whitespace;
   - it merges tags that differ only in case;
   - it rescales facings;
   - it cannot store and reload a `parameters` key;
   - it rejects quoted line breaks.
