# KarstNSim — Instruction File Parameter Reference

This document explains every parameter accepted in the KarstNSim instruction file (parsed into `ParamsSource` and then propagated to the simulation).  
Lines starting with `//`, `%`, `::`, or `#` are treated as comments and ignored by the parser. 
Each active declaration must use the exact syntax `parameter_name: value`. Unknown tags, duplicated declarations, missing values, invalid booleans, malformed numerical values, and missing conditionally mandatory parameters raise an explicit exception before external simulation objects are prepared. 
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

### `simulation_input_dir: <path>`
- **Type**: `path` (optional)
- **What it does**: Overrides the directory used to locate input files that are resolved relative to the instruction file, including `connectivity_matrix.txt`.
- **Default behavior**: If omitted, KarstNSim uses the directory containing the active instruction file.
- **Practical effect**: Keep this parameter omitted for self-contained simulation folders; define it only when the instruction file and its auxiliary ASCII inputs are intentionally stored in different directories.


---

## 2) General parameters

### `domain: <path>`
- **Type**: `path` (text voxet/box description with optional density property and IKP property)
- **What it does**: Defines the 3D domain (origin, axes, grid size) used for sampling and optional grid-based properties (density and IKP (Intrinsic Karstification Potential)). Density values range from 0 to 1 both excluded, and IKP from 0 to 1 both included. Any negative value is considered a NDV, but -99999 is the default NDV for both.
- **Indexing reminder**: KarstNSim expects grid properties (`density`, `IKP`, etc.) to be flattened with the **u axis as the fastest varying index**, followed by **v**, then **w**. In other words, the linear index must be computed as `idx = u + nu * v + nu * nv * w` where `0 ≤ u < nu`, `0 ≤ v < nv`, `0 ≤ w < nw`. When exporting or building voxets from external tools, make sure the serialization order matches this convention; mismatches will map sampled points to the wrong cells.

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

You have the option to import a graph onto which all the simulation steps will be computed (simulation support), or create a new one based on an N-nearest neighbor logic.
Be aware that if you import an input graph, this graph should include exactly all the sampling points, the inlets, the outlets, the waypoints, and the surface points as its nodes, and no other point.

1) **Import a graph**

- `use_input_nghb_graph: <bool>` — **Type**: `bool`  
- `input_nghb_graph: <path>` — **Type**: `path` (existing graph in a graph format)  
**Behavior**: Imports the input graph to the simulation. This is especially useful if (i) you've already simulated a network and saved the previously generated nghb graph, and want to reuse exactly the same for subsequent simulations, but this does not really save computation time, or (ii) you need the KarstNSim network to be conformal to a given grid (ie. a hydrogeological model grid) in order to perform groundwater simulation afterwards.  

2) ** Create a NGHB graph during simulation **

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
**Behavior**: Define inlet(s) and outlet(s).
Inlets MUST be associated with two or three properties: 
- First property: Index (1-based), which allows to pair them in the connectivity matrix with springs
- Second property: Order (1-based too), which controls the order of iteration of sinks (order 1 first, then 2 etc.). Sinks of the same order are iterated in a random order within that group.
- Third property (optional): Radii (or any other section property) can be optionally read per sink node when the respective `use_*_radius` flag is true. These radii will be used as observation data to constrain the conduit size simulation if activated. They do not have another effect.
Springs must be associated with two or three properties:
- First property: Index property (1-based), to pair them with inlets in the connectivity matrix
- Second property: Surfindex property which links each spring to its watertable (watertable file names end with a number ; that number is their index, it is 1-based). It is possible to define Surfindex=0 for springs that are perched / not associated with any perennial water body.
- Third property (optional): Radii (or any other section property) can be optionally read per spring node when the respective `use_*_radius` flag is true. These radii will be used as observation data to constrain the conduit size simulation if activated. They do not have another effect.

Note: in a sections-only simulation, only the radii properties are (optionally) required, but the 3-property format is still accepted.

### Connectivity source and outlet selection for ambiguous links

- `use_user_connectivity_matrix: <bool>` — **Type**: `bool`  
  **Mandatory status**: Required for a full network simulation; ignored in `sections_simulation_only` mode.  
  **Behavior**:
  - `true`: read and strictly validate `connectivity_matrix.txt` from `simulation_input_dir`;
  - `false`: ignore any matrix file supplied by the user and create a sink-by-spring matrix filled exclusively with `2` values directly in memory.
  
  In generated mode, the initial all-`2` matrix is not written to disk. The independent `create_solved_connectivity_matrix` option may still export the matrix after ambiguous associations have been resolved.

- `allow_single_outlet_connection: <bool>` — **Type**: `bool`  
**Behavior**: For connectivity entries equal to `2`, this toggles between choosing the **closest** spring in terms of corrected cumulative path cost (`true`) and **randomly** selecting an admissible spring (`false`).  
**Practical effect**: When `use_user_connectivity_matrix = false`, every entry is `2`; consequently this parameter controls the routing strategy for the entire generated connectivity matrix. Set it to `true` for cost-based, more deterministic routing. This option can also help delineate probable spring catchments from the cost field.

### Waypoints
- `use_waypoints: <bool>` — **Type**: `bool` 
- `waypoints: <path>` — **Type**: `path` (point set with mandatory "waypoints_impact_radius" property)
- `waypoints_weight: <float>` — **Type**: `float` (≥ 0 and < 1)  
- `use_waypoints_radius: <bool>` — **Type**: `bool` (read per-waypoint radius if true)  
**Behavior**: Waypoints locally reduce costs based on their weight (costs are locally multiplied by the waypoint weight). An Impact Radius property must be defined for each waypoint (in the pointset file),
 which defines the radius of impact of each waypoint to reduce the cost. This is different from the optionally added Radius property,
 which can be added to the waypoints if `use_waypoints_radius` is true, and which will locally constrain the conduit size simulation in the same way as sinks and springs.
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

## 8) Intrinsic karstification potential (IKP)

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
**Important**: The format of the watertable surface file name must very precisely be "(name)_k" with k the index of the water table. The spring pointset has a mandatory property called surfindex which links each spring to its watertable k.
The user can define perched/relict springs without associated water bodies by defining surfindex = 0 and not associating any water table.
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
**Important**: Variogram `sill` and `nugget` parameters are expected to be defined in **Gaussian normal-score space**, since SGS is performed after normal-score transformation of the simulated property. Variogram `range` parameters remain curvilinear distances expressed in the model coordinate units.
- **Neighborhood & mixing**:  
  - `number_max_of_neighborhood_points`, `nb_points_interbranch`, `proportion_interbranch`
  
**Behavior**: Performs curvilinear SGS along the skeleton to simulate `eq_radius` at nodes, combining a global field with structured inter- and intra-branch covariance. It is advised to use log(radius) for more linearity in computations. You can also define on inlets, outlets, but also waypoints (only if they are traversed directly by a conduit) a radius property that will act as observation data in the SGS.
**Practical effects**:  
- Larger `range` ⇒ smoother variations along branches; higher `nugget` ⇒ noisier radii; `sill` sets total variance.  
- `nb_points_interbranch` and `proportion_interbranch` tune how much cross-branch correlation is allowed.  
- `*_range_of_neighborhood` and `number_max_of_neighborhood_points` control kriging neighborhood max size.
**Allowed Kriging Models**: `"Exponential"`, `"Spherical"`, `"Gaussian"`

** External drift **: These options allow to add an external drift encompassing long-range variations in conduit size. The drift is separated between a vadose-phreatic trend and an upstream-downstream trend (conduits becoming larger in the phreatic zone and in the downstream part of the network, respectively). Both can be applied at the same time. The drift values are directly inferred from observational data.
** Warning **: It is not possible to add a drift if no data (inlets, outlets, waypoints radii) are provided. At least two data are needed, but preferably min. 10 to 20.
- `use_drift_zwt: <bool>` — **Type**: `bool` (to add a vadose/phreatic trend)
- `use_drift_curv: <bool>` — **Type**: `bool` (to add an upstream/downstream trend)
  
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

- `vertical_distance_stretching_factor: <float>` — **Type**: `float` (> 0)
  **Behavior**: Vertical anisotropy factor used for selected graph-related distances. Any computed Euclidean distance will have its vertical component scaled by the provided factor. This includes cost computation and distances for amplification (dead-end points and cycle distances). A large factor (eg. 10) can be used if conduits are expected to propagate as horizontally as possible.

- `gradient_constraint_weight: ` — **Type**: `float` (≥ 0, default `0`)  
  **Behavior**: When `allow_single_outlet_connection = true`, ambiguous inlet/outlet links
  (`2` values in the connectivity matrix) are selected using a corrected cumulative path cost:

  `corrected_cost = cumulative_cost * (1 + gradient_constraint_weight * normalized_gradient_penalty)`

  where `normalized_gradient_penalty` is computed, for each inlet, from the vertical drop
  between this inlet and each candidate outlet. The outlet with the largest vertical drop
  receives a penalty of `0`; the outlet with the smallest vertical drop receives a penalty
  of `1`; intermediate outlets are linearly interpolated.

  **Practical effect**: Higher values favor outlets located at lower elevation, i.e. paths
  with stronger inlet-to-outlet hydraulic gradient, even if their raw cumulative graph cost
  is larger.
  
- `outlet_selection_cost_factor: ` — **Type**: `float` (≥ 1, default `1`)  
  **Behavior**: When `allow_single_outlet_connection = true`, ambiguous inlet/outlet links
  (`2` values in the connectivity matrix) are filtered using the corrected cumulative
  path cost.

  Let `best_corrected_cost` be the minimum corrected cumulative cost among the candidate
  outlets for a given inlet. KarstNSim keeps every ambiguous candidate path satisfying:

  `corrected_cost <= best_corrected_cost * outlet_selection_cost_factor`

  The corrected cost includes the optional `gradient_constraint_weight` correction.
  Therefore, this factor is applied after the elevation-drop correction, not on the raw
  cumulative graph cost. This will typically be used to link an inlet to outlets part of the same
  system, which would typically be close together and yield a very similar cumulative cost from the inlet.

  **Examples**:
  - `1`: keep only the minimum corrected-cost path.
  - `1.1`: keep paths up to 10% more expensive than the best corrected path.
  - `2`: keep paths up to 100% more expensive than the best corrected path.


---

## 18) Connectivity matrix modes

The active mode is selected by `use_user_connectivity_matrix`.

### User-matrix mode (`use_user_connectivity_matrix: true`)

- **Required file**: `simulation_input_dir/connectivity_matrix.txt`.
- **Shape**: rows = sinks, columns = springs.
- **Allowed integer values**:
  - `0` → **forbid** connection `i → j`;
  - `1` → **force** connection `i → j`;
  - `2` → **leave the association ambiguous** and resolve it according to `allow_single_outlet_connection`.

The row and column identifiers are defined by the one-based `Index` properties of the sink and spring point sets.

### Generated all-ambiguous mode (`use_user_connectivity_matrix: false`)

KarstNSim does not read `connectivity_matrix.txt`, even if such a file is present. It creates an in-memory matrix with dimensions `number_of_sinks × number_of_springs` and initializes every entry to `2`. This initial generated matrix is not saved as an input artifact. Ambiguous links are then resolved by the closest-cost strategy when `allow_single_outlet_connection: true`, or by random admissible selection when it is `false`.

No connectivity matrix is required in `sections_simulation_only` mode because the network topology is read from `previous_networks`.

---

## 19) Saving extra artifacts

- `create_vset_sampling: <bool>` — save the sampling points as a Pointset.  
- `create_nghb_graph: <bool>` — save the neighbor graph as a polyline object. Warning: Heavy object!
- `create_nghb_graph_property: <bool>` — add per-segment properties (cost components) to the saved neighbor graph. Warning: Very Heavy object!
- `create_solved_connectivity_matrix: <bool>` — save the connectivity matrix after all ambiguous (`2`) associations have been resolved. This output option is independent of `use_user_connectivity_matrix`; it can therefore be enabled for both user-matrix and generated all-ambiguous modes.
- `create_grid: <bool>` — export a voxet with useful scalar properties (e.g., density, IKP). This is useful because ghost-rocks printed onto IKP modify the IKP property during the simulation.


---