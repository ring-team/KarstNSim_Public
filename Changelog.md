# KarstNSimPublic


## Version 2.1

08/07/2026

### Modifications

This update adds a new option:

- Added a new parameter `use_user_connectivity_matrix` to explicitly select the connectivity source: if true, the user must provide the connectivity matrix. If false, the matrix is populated entirely with ambiguous (`2`) connections.

Several code improvements:

- Added a stricter instruction file parser that detects errors in the file (parameter range, non-provided mandatory parameter, missing values, unknown tags...).
- Improved the ASCII input readers: they now detect and fail explicitly when malformed ASCII files are provided.
- Several fixes. Notably:
	- Fixed nearest-neighbor graph operations to safely handle empty point clouds.
	- Fixed cycle amplification robustness.
	- Fixed dead-end amplification to avoid invalid empty inputs and potentially infinite candidate-generation loops.
	- Added consistency and validity checks across graph amplification methods.
- Changed several functionalities of the CB-SGS with external drift algorithm:
	- Variogram parameters are now explicitly expected to be in Gaussian space, not in the simulated-property space. An internal switch (not accessible to the user) allows to change this behavior if needed. In that case, an empirical normal-score variogram converter is used (should be used with caution).
	- The matrix singularity detection is now more robust by being scale-aware (instead of checking for an exactly zero pivot in the matrix).
	- The explicit inversion of the kriging system has been replaced by a direct Cholesky decomposition, reducing computational cost and increasing numerical robustness.
	- Reworked the MAD-based outlier filtering (during the regression) to work on more robust leave-one-out (LOO) residuals.
	- Changed redundancy weighting from a single predictor axis to a joint predictor space density using both drift variables (when both are active).
	- Changed the behavior of outlier removal when only one inlet or outlet datum is provided: in that case the single inlet/outlet is preserved and not removed independently of its value.

## Version 2.0

07/07/2026.

### Modifications

This update adds several new options:

- A graph can be given as an input, instead of generating the NGHB graph. This is especially useful if the KarstNSim karst skeleton should be conformal to another grid, for instance a groundwater model grid. The grid should simply be transformed into an edges & nodes graph format and imported as an input to KarstNSim.
- An external drift can now be superimposed onto the CB-SGS simulation of conduit sizes. It it automatically inferred based on observations.
- New shortest path computation options are available :
	- vertical_distance_stretching_factor: Vertical anisotropy factor used for selected graph-related distances. Any computed Euclidean distance will have its vertical component scaled by the provided factor. This includes cost computation and distances for amplification (dead-end points and cycle distances). A large factor (eg. 10) can be used if conduits are expected to propagate as horizontally as possible.
	- gradient_constraint_weight: Factor multiplied to cumulative cost distances to reflect the influence of the gradient (difference in altitude) between inlets and outlets.
	- outlet_selection_cost_factor: Tolerance factor to keep paths to outlets almost as close as the closest outlet to a given inlet. Useful for outlets close together and part of the same system.

And several code improvements:

- Possibility to not define a water table for certain (or all) springs. These springs are therefore relict or overflow springs that the user cannot associate with a particular perennial water body.
- Improved shortest-path computation by using an exact bidirectional heap-based Dijkstra algorithm for inlet-outlet path searches (~2x speed).
- Several bug fixes.
- We provide access to a launcher for parallel simulations (parallel_launcher.py). Each job can be comprised of one KarstNSim simulation or a sequence of simulation (each with specific parameters).


## Version 1.3

09/12/2025.

### Modifications

Minor update in terms of code functionalities:
- Added an input configuration reference manual (config_reference.md) which precises the type, theoretical and typical range of values, meaning, behavior and practical effect of all input parameters. Users can now refer to it for correct KarstNSim use, which avoids the requirement to generate the doxygen documentation for most users.
- Updated the "base" example provided with some new functionalities, and added three other examples:
	- A "Polygenic Karst" example which simulated the carving of new conduits following a base level drop from the base case.
	- An "Amplification" example which uses the amplification step on the network generated in the "polygenic karst" example.
	- A "Karst section generation" example, which, based on the amplified network, simulates an equivalent radius on the nodes of the network using a curvilinear SGS algorithm (Frantz et al. 2021).
- Added more log lines during simulation so that the user is kept informed of the completion and time duration of each step. In case of a crash, the user can also, as a result, provide more precise information to the developer for a fix.
- Some minor bug fixes.

## Version 1.2

02/14/2025.

### Modifications

Major update with new features from the PhD thesis (Gouy, 2025), including:
- Amplification step to increase the network density and generate maze patterns (cycles). There is also an option to add noise to the cost graph, either just for the amplification step, or for the whole simulation.
- New algorithm to pick inlet/outlet pairing when uncertain (no tracer tests and/or no strong geological argument): either by picking the pairings randomly, or only keeping the path from each inlet to the closest spring (by computing the cumulative shortest path cost). This allows to draw probable catchment areas of each spring.
- Possibility to incorporate ghost-rocks in the simulation through surface alteration lines: the code will first interpolate an alterite volume beneath the line (with elliptic cross-section) and then reduce the cost of edges traversing this volume (through the intrinsic karstification potential subcost).
- New key points available : waypoints and karst-free points. Both affect the path chosen during shortest path computations. Waypoints can be seen as soft data constraints, by reducing cost of all edges less than a given radius away from a waypoint (can be used when some conduits positions are known) and no-karst points also have a spherical effect zone around them which simply deletes any sampling point in it, hence blocking any path close to them.
- Conduit dimensions simulation step using the 1D-Curvilinear Branchwise modified SGS algorithm of Frantz et al. (2021). Allows to generate any property on the skeleton nodes, such as Re (equivalent radius) and WH (widht-height ratio).

## Version 1.1

03/14/2024.

### Modifications

- Optimization of sampling, surface detection, and Dijkstra implementation, mostly with the use of a K-D tree structure.
- Debugging of the public version of the code:
  - Point-surface detection is now much more precise (no resolution artifacts of the karst network at the borders of triangulated surfaces, as described in the article (last paragraph of section 4.3)).
  - The two-step shortest path computation between inlets and outlets was changed slightly for more realistic results in complex geometries: before, the first Dijkstra computation (vadose shortest path) was made between the inlet and any point below the water table. Now, the computation is between the inlet and the closest point precisely onto the water table (not below, not above).
  - Corrected bug that would only compute costs based on the position of the water table n°1, never using other water tables when necesary.
- Added 2 new options:
  - Apply cost reduction only in the phreatic zone (instead of everywhere).
  - Option to keep only the path to the closest outlet for each inlet if more than one "1" is present in a single connectivity matrix line.
- Changed the cost function from C = L(αICI + αFCF + αWTCWT + αPCP) to C = L(1 + αICI + αFCF + αWTCWT + αPCP). The added "+1" simplifies counterbalancing of costs, particularly the fracture cost, to control the level of fracturation of the network more easily, while maintaining complete proportionality between distance L and full cost C (and thus avoiding resolution artifacts).
- A property "equivalent_radius" appears in the karst network output file. It is not currently implemented, hence returning a default value for now.

## Version 1.0

02/14/2024.
Original version, as used in the 2024 [article](https://doi.org/10.1016/j.jhydrol.2024.130878).