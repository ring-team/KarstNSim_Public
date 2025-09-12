# KarstNSim — Instruction File Parameter Reference

This document explains every parameter accepted in the KarstNSim instruction file (parsed into `ParamsSource` and then propagated to the simulation).  
Lines starting with `//`, `%`, `::`, or `#` are treated as comments and ignored by the parser.  
Unless stated otherwise, distances/lengths are in the same units as your model coordinates (typically meters).

---

## 1) Project name & repository paths

### `karstic_network_name: <string>`
- **Type**: `string`
- **What it does**: Scenario/run name used as a prefix for outputs and logs.
- **Practical effect**: Change it to partition outputs per scenario or batch.

### `main_repository: <path>`
- **Type**: `path`
- **What it does**: Root directory containing your `Input_files/` and where `outputs/` will be created.
- **Practical effect**: Set this to the repository root of the archive.


---

## 2) General parameters

### `domain: <path>`
- **Type**: `path` (text voxet/box description with density property and optionally IKP property)
- **What it does**: Defines the 3D domain (origin, axes, grid size) used for sampling and grid-based properties (density and IKP (Intrinsic Karstification Potential), the second being optional). Density values range from 0 to 1 both excluded, and IKP from 0 to 1 both included. Any negative value is considered a NDV, but -99999 is the default NDV for both.

### `selected_seed: <int>`
- **Type**: `int` (≥ 0)
- **What it does**: Random seed for reproducibility.

### `number_of_iterations: <int>`
- **Type**: `int` (≥ 1)
- **What it does**: Number of simulation iterations.

### `vary_seed: <bool>`
- **Type**: `bool`
- **What it does**: If true, different iterations will use different seeds. Intended to be kept true.

### `topo_surface: <path>`
- **Type**: `path` (triangulated surface)
- **What it does**: Topographic surface used to validate sample points (reject above-ground points).


---

## 3) Sampling points

You can **reuse an existing point set** or **sample new points**.

### Reuse existing points
- `use_sampling_points: <bool>` — **Type**: `bool`  
- `sampling_points: <path>` — **Type**: `path` (point set)  
**Behavior**: If `true`, the provided points are used as the sampling cloud from which the neighbor graph is built.

### Sample new points
Two modes are available; choose **one**.

1) **Constant-density Poisson sampling**
- `poisson_radius: <float>` — **Type**: `float` (> 0 and typically < 0.1)  
  **Behavior**: Sets the (approximate) minimum inter-point spacing. A value of 0.01 means that the minimum distance between two sampled points in the smallest domain dimension will be at least 1% of that dimension size. Smaller values produce denser clouds.  
  **Optionality**: Used when `use_density_property = false`.

2) **Property-driven density (from box/grid)**  
- `use_density_property: <bool>` — **Type**: `bool`  
  **Behavior**: A density property defined on the domain grid controls local spacing.
  **Practical effects**: Denser sampling increases graph connectivity and resolution, but raises memory/CPU. Use property-driven density to concentrate nodes where karstification is expected.

- `k_pts: <int>` — **Type**: `int` (≥ 1)
  **Behavior**: `k_pts` is the number of candidate trials around each active sample (higher → better sampling quality, but higher computation time). A good compromise is between 10 and 30.

---

## 4) N-nearest neighbor graph

- `nghb_count: <int>` — **Type**: `int` (≥ 1)  
- `use_max_nghb_radius: <bool>` — **Type**: `bool`  
- `nghb_radius: <float>` — **Type**: `float` (> 0)  
**Behavior**: Each point connects to up to `nghb_count` neighbors, optionally capped by `nghb_radius`.  
**Practical effects**: Larger `nghb_count` increases simulation precision by allowing more edge directions but produces higher computation times; a radius cap avoids unrealistic long edges. Typical nghb_count should be between 100 and 200.

---

## 5) Inlets, outlets & waypoints

### Sinks & springs

- `sinks: <path>` — **Type**: `path` (point set with mandatory "Index" and "Order" properties) 
- `use_sinks_radius: <bool>` — **Type**: `bool` (read per-sink radius from properties if true)  
- `springs: <path>` — **Type**: `path` (point set with mandatory "Index" and "Surfindex" properties)  
- `use_springs_radius: <bool>` — **Type**: `bool`  
**Behavior**: Define inlet(s) and outlet(s). Inlets MUST be associated with two properties: Index, which allows to pair them in the connectivity matrix with springs, and Order, which controls the order of iteration of sinks (order 1 first, then 2 etc.). Sinks of the same order are iterated in a random order within that group. Springs must be associated with an Index property, to pair them with inlets in the connectivity matrix, and with a Surfindex property which links each spring to its watertable (watertable file names end with a number ; that number is their index). Radii can be optionally read per point when the respective `use_*_radius` flag is true. These radii will be used as observation data to constrain the conduit size simulation.

### Outlet selection for ambiguous links
- `allow_single_outlet_connection: <bool>` — **Type**: `bool`  
**Behavior**: When the connectivity matrix indicates uncertainty (value `2`), this toggles between choosing the **closest** spring (in the sense of the shortest cumulative cost path, not the euclidean shortest path) (`true`) or **randomly** selecting one (`false`).  
**Practical effect**: Set `true` for more deterministic routing to the nearest outlet. This option can also help in delineating probable catchment areas of springs, based on the cost function field.

### Waypoints
- `use_waypoints: <bool>` — **Type**: `bool` 
- `waypoints: <path>` — **Type**: `path` (point set with mandatory "waypoints_impact_radius" property)
- `waypoints_weight: <float>` — **Type**: `float` (≥ 0 and < 1)  
- `use_waypoints_radius: <bool>` — **Type**: `bool` (read per-waypoint radius if true)  
**Behavior**: Waypoints locally reduce costs based on their weight (costs are locally multiplied by the waypoint weight). An Impact Radius property must be defined for each waypoint (in the pointset file), which defines the radius of impact of each waypoint to reduce the cost. This is different from the optionally added Radius property, which can be added to the waypoints and which will locally constrain the conduit size simulation.
**Practical effect**: Waypoints can be seen as soft-data constraints. They are used to honor surveyed conduits/passages or enforce transit through checkpoints.

---

## 6) Ghost-rocks

- `use_ghostrocks: <bool>` — **Type**: `bool`  
- `alteration_lines: <path>` — **Type**: `path` (polyline set)  
- `interpolate_lines: <bool>` — **Type**: `bool` (Not implemented yet, keep false)  
- `ghostrock_max_vertical_size: <float>` — **Type**: `float` (≥ 0)  
- `ghostrock_width: <float>` — **Type**: `float` (≥ 0)  
- `use_max_depth_constraint: <bool>` — **Type**: `bool`  
- `max_depth_horizon: <path>` — **Type**: `path` (surface; required if previous is true)  
- `ghost_rock_weight: <float>` — **Type**: `float` (≥ 0)  
**Behavior**: “Paints” karstification potential property along alteration polylines, limited by vertical size, lateral width and optional max-depth horizon. Requires `use_karstification_potential = true` to be used. The ghost-rock weight will be the weight of the ghost-rock in the overall karstification potential cost 
**Practical effects**: Increases "karstifiability" in corridors aligned with the lines.

---

## 7) Inception surfaces

- `add_inception_surfaces: <bool>` — **Type**: `bool`  
- `refine_surface_sampling: <int>` — **Type**: `int` (≥ 0, or -1 to avoid adding TIN nodes to the sampling cloud)  
- `inception_surfaces: <list<path>>` — **Type**: list of paths (one or more surfaces)  
- `inception_surface_constraint_weight: <float>` — **Type**: `float` (≥ 0)  
- `max_inception_surface_distance: <float>` — **Type**: `float` (> 0)  
**Behavior**: Enrich sampling on each surface by using the internal discretization of TINs, and add a sub-cost term pulling paths towards them, limited by a maximum 3D connection distance. The user can control the level of surface refinement with refine_surface_sampling. A value of 0 means no refinement, -1 means not adding the TIN nodes, 1+ means adding mid-points on each triangle edge of the TIN, in a fractal pattern.
**Practical effects**: Higher weight locks paths to structural planes (bedding, fault interfaces). Use the max inception surface impact distance to fine-tune the "magnet" effect of inception surfaces.

---

## 8) Intrinsic arstification potential (IKP)

- `use_karstification_potential: <bool>` — **Type**: `bool`  
- `karstification_potential_weight: <float>` — **Type**: `float` (≥ 0)  
**Behavior**: Enables a cost term favoring high-IKP rock. IKP values (ranging from 0 to 1 included) are read from the domain grid property, that MUST be provided if the karstification potential option is enabled.  
**Practical effects**: Larger weight increases the tendency to traverse high-IKP volumes.

---

## 9) Fracture constraints

- `use_fracture_constraints: <bool>` — **Type**: `bool`  
- `fracture_families_orientations: <list<float>>` — **Type**: list of angles separated by single spaces (° clockwise from North, between 0 and 180)  
- `fracture_families_tolerance: <list<float>>` — **Type**: list of angles separated by single spaces (°; same length as orientations)  
- `fracture_constraint_weight: <float>` — **Type**: `float` (≥ 0)  
**Behavior**: Penalizes edges deviating from target fracture orientation sets; tolerance acts as an angular half-width around each set. Tolerance can typically range between 5 and 20°. A smaller tolerance can lead to the local impossibility to find an edge exactly in the right direction.
**Practical effects**: Higher weight and tighter tolerance produce more angular, fracture-aligned networks.

---

## 10) Previous networks / polyphasic karstification

- `use_previous_networks: <bool>` — **Type**: `bool`  
- `previous_networks: <path>` — **Type**: `path` (existing network polylines)  
- `fraction_old_karst_perm: <float>` — **Type**: `float` in (0, 1]  
**Behavior**: Reimports existing network segments into the sampling and reduces their costs by multiplying with `fraction_old_karst_perm` (i.e., cheaper edges along old karst).  
**Practical effects**: Encourages reuse of paleokarst pathways; smaller fractions give stronger “prefer old path” bias. Can be used to simulate polyphasic networks.

### `sections_simulation_only: <bool>`
- **Type**: `bool`  
- **Behavior**: If true, skip graph generation and only simulate section properties (equivalent radii) on a provided network. use_previous_networks MUST be set to true to be used.

---

## 11) No-karst spheres

- `use_no_karst_spheres: <bool>` — **Type**: `bool`  
- `sphere_centers: <path>` — **Type**: `path` (point set with radius property)
**Behavior**: Defines spherical exclusion zones. The pointset defines the center of the spheres and the radius property defines the radius of each sphere. Edges intersecting spheres are removed from the simulation.  
**Practical effects**: Deflects paths away from forbidden volumes (e.g., non-karst units or protected areas). In practice you can also achieve the same effect by defining a NDV (-99999) density in the density property.

---

## 12) Water tables & vadose/phreatic trends (mandatory)

- `surf_wat_table: <list<path>>` — **Type**: list of surfaces (text files)
- `water_table_constraint_weight_vadose: <float>` — **Type**: `float` (≥ 0)  
- `water_table_constraint_weight_phreatic: <float>` — **Type**: `float` (≥ 0)  
**Behavior**: Classifies nodes as vadose/phreatic relative to their associated water table and adds dedicated cost terms with independent weights for vadose and phreatic zones.
**Important**: The format of the watertable surface file name must very precisely be "(name)_k" with k the index of the water table. The spring pointset has a mandatory property called surfindex which links each spring to its watertable k. There is currently no option to avoid defining any water table in the simulation ; this option will be developed in the future (perched springs).
**Practical effects**: Larger weights strengthen vertically downward conduits in the vadose zone (and penalize any horizontal or vertically upward conduit) and base level (water table) control in the phreatic zone.


---

## 13) Deadend points amplification

- `use_deadend_points: <bool>` — **Type**: `bool`  
- `nb_deadend_points: <int>` — **Type**: `int` (≥ 0)  
- `max_distance_of_deadend_pts: <float>` — **Type**: `float` (≥ 0)  
**Behavior**: Adds dead-end points during the amplification phase. They are points randomly sampled around the main network, and then connected to it using a shortest path computation.
**Practical effects**: One can see dead-end branches as smaller scale branches of the network that are not connected to a surface sink. Using this option will densify the network.


---

## 14) Cycle amplification

- `use_cycle_amplification: <bool>` — **Type**: `bool`  
- `max_distance_amplification: <float>` — **Type**: `float` (≥ 0)  
- `min_distance_amplification: <float>` — **Type**: `float` (≥ 0; ≤ max)  
- `nb_cycles: <int>` — **Type**: `int` (≥ 0)  
**Behavior**: Adds loops (cycles) to the skeleton within the specified distance range.s
**Practical effects**: Larger distances and more cycles produce more "mazey" networks, which is ideal to represent either epiphreatic mazes or epikarst fracture mazes (or even spongelike paleokarst).

---

## 15) Noise amplification

- `use_noise: <bool>` — **Type**: `bool`  
- `use_noise_on_all: <bool>` — **Type**: `bool` (apply noise for whole simulation (true) or only for amplification (false))
- `noise_frequency: <int>` — **Type**: `int` (≥ 1)  
- `noise_octaves: <int>` — **Type**: `int` (≥ 1)  
- `noise_weight: <float>` — **Type**: `float` (≥ 0)  
**Behavior**: Applies simplex-noise–based heterogeneity to cost field. Frequency/octaves set spatial texture; weight sets intensity.  
**Practical effects**: In all cases, use moderate values to avoid overpowering the other costs.


---

## 16) Sections (equivalent radii) simulation

- `simulate_sections: <bool>` — **Type**: `bool`  
- `simulation_distribution: <path>` — **Type**: `path` (file with initial distribution ==> a set of values sampled from the initial distribution (50 min. values is better))  
- **Global variogram**:  
  - `global_vario_range`, `global_range_of_neighborhood`, `global_vario_sill`, `global_vario_nugget`, `global_vario_model`
- **Inter-branch variogram**:  
  - `interbranch_vario_range`, `interbranch_range_of_neighborhood`, `interbranch_vario_sill`, `interbranch_vario_nugget`, `interbranch_vario_model`
- **Intra-branch variogram**:  
  - `intrabranch_vario_range`, `intrabranch_range_of_neighborhood`, `intrabranch_vario_sill`, `intrabranch_vario_nugget`, `intrabranch_vario_model`
- **Neighborhood & mixing**:  
  - `number_max_of_neighborhood_points`, `nb_points_interbranch`, `proportion_interbranch`
  
**Behavior**: Performs curvilinear SGS along the skeleton to simulate `eq_radius` at nodes, combining a global field with structured inter- and intra-branch covariance. It is advised to use log(radius) for more linearity in computations. You can also define on inlets, outlets, but also waypoints (only if they are traversed directly by a conduit) a radius property that will act as observation data in the SGS.
**Practical effects**:  
- Larger `range` ⇒ smoother variations along branches; higher `nugget` ⇒ noisier radii; `sill` sets total variance.  
- `nb_points_interbranch` and `proportion_interbranch` tune how much cross-branch correlation is allowed.  
- `*_range_of_neighborhood` and `number_max_of_neighborhood_points` control kriging neighborhood max size.
**Allowed Kriging Models**: `"Exponential"`, `"Spherical"`, `"Gaussian"`
**Note**: A functionality to incorporate external drifts in the curvilinear branchwise SGS algorithm is currently being developed.

---

## 17) Cost graph combination & cohesion

- `gamma: <float>` — **Type**: `float` in (0, 2]  
  **Behavior**: gamma-graph pruning parameter (anisotropic empty-region rule). Deprecated (unused) parameter. Might be reincorporated in the future for special cases.

- `fraction_karst_perm: <float>` — **Type**: `float` in (0, 1]  
  **Behavior**: After selecting a path/iteration, edge costs along it are multiplied by this fraction (cohesion). Lower = stronger cohesion.  
  **Practical effects**: Smaller values produce more hierarchical, cohesive networks over successive iterations. Since the cohesion factor is applied at EACH iterated inlet, the chosen value should be somewhat proportional to the number of iterated inlets. A value ranging from 0.9 for 1-10 inlets to 0.99 for 500+ inlets yields good results.

- `vadose_cohesion: <bool>` — **Type**: `bool`  
  **Behavior**: If true, cohesion applies everywhere; if false, restrict cohesion to phreatic areas only (this option was added because it is debatable whether cohesion (ie. backwards carving influence) only appears in saturated media or not).

- `multiply_costs: <bool>` — **Type**: `bool`  
  **Behavior**: Combine sub-costs by multiplication instead of summation. Not recommended, but might be useful in some cases (cost penalties become much more predominant, even if only in one sub-cost).

---

## 18) Connectivity matrix (separate file)

- **File**: typically `Input_files/connectivity_matrix.txt` (or similar).  
- **Meaning (row = sink i, column = spring j)**:  
  - `0` → **forbid** connection `i → j`  
  - `1` → **force** connection `i → j`  
  - `2` → **let the algorithm decide** among candidates (subject to `allow_single_outlet_connection`)  
**Notes**:  
- This parameter is mandatory because it controls which inlets are connected to which outlets. The inlet and outlet index in the matrix (line and column number resp.) are defined by the Index property in the two pointset files.


---

## 19) Saving extra artifacts

- `create_vset_sampling: <bool>` — save the sampling points as a Pointset.  
- `create_nghb_graph: <bool>` — save the neighbor graph as a polyline object. Warning: Heavy object!
- `create_nghb_graph_property: <bool>` — add per-segment properties (cost components) to the saved neighbor graph. Warning: Very Heavy object!
- `create_grid: <bool>` — export a voxet with useful scalar properties (e.g., density, IKP). This is useful because ghost-rocks printed onto IKP modify the IKP property during the simulation.


---