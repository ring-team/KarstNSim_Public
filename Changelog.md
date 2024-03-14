# KarstNSimPublic

## Version 1.0

02/14/2024.
Original version, as used in the [article](https://doi.org/10.1016/j.jhydrol.2024.130878).

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