# KarstNSimPublic

## Version 1.0

02/14/2024.
Original version, as used in the 2024 [article](https://doi.org/10.1016/j.jhydrol.2024.130878).

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

## Version 1.2

02/14/2025.

### Modifications

Major update with new features from the PhD thesis (Gouy, 2025), including:
- Amplification step to increase the network density and generate maze patterns (cycles). There is also an option to add noise to the cost graph, either just for the amplification step, or for the whole simulation.
- New algorithm to pick inlet/outlet pairing when uncertain (no tracer tests and/or no strong geological argument): either by picking the pairings randomly, or only keeping the path from each inlet to the closest spring (by computing the cumulative shortest path cost). This allows to draw probable catchment areas of each spring.
- Possibility to incorporate ghost-rocks in the simulation through surface alteration lines: the code will first interpolate an alterite volume beneath the line (with elliptic cross-section) and then reduce the cost of edges traversing this volume (through the intrinsic karstification potential subcost).
- New key points available : waypoints and karst-free points. Both affect the path chosen during shortest path computations. Waypoints can be seen as soft data constraints, by reducing cost of all edges less than a given radius away from a waypoint (can be used when some conduits positions are known) and no-karst points also have a spherical effect zone around them which simply deletes any sampling point in it, hence blocking any path close to them.
- Conduit dimensions simulation step using the 1D-Curvilinear Branchwise modified SGS algorithm of Frantz et al. (2021). Allows to generate any property on the skeleton nodes, such as Re (equivalent radius) and WH (widht-height ratio).
