/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for new methods + modifications to original methods
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Rewritten entirely (except the variogram_value fonction, written by G. Rongier in 2015) from the pseudo-code algorithm (Algorithm 2) in the paper of Frantz et al., 2021 "Analysis and stochastic simulation of geometrical properties of conduits in karstic networks".
This work was performed in the frame of the RING project at Université de Lorraine.

***************************************************************/

#pragma once

/**
@file geostats.h
@brief Contains methods used in KarstNSim for 1D-Curvilinear branchwise SGS. It is a 1:1 translation of the code of Frantz et al. (2021), but devoid of any dependencies (as opposed to the initial version).
@author Augustin GOUY
**/

#include <queue>
#include <limits>
#include <algorithm>
#include <string>
#include <vector>
#include "KarstNSim/graph.h"
#include "KarstNSim/randomgenerator.h"
#include "KarstNSim/vec.h"
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
/**
 * \brief Hydraulic role of a conduit-size conditioning datum.
 *
 * The role is used by the leave-one-out external-drift diagnostics to avoid
 * testing the sole inlet or outlet observation of a hydraulic boundary class.
 */
enum class ConditioningDataRole {
	None,
	Inlet,
	Outlet,
	Waypoint
};

/*!
	\struct GeostatParams
	\brief Struct to group all input parameters for geostatistical simulation of conduit shapes.

	This structure encapsulates parameters required for geostatistical simulations, including
	variogram parameters, neighborhood properties, and proportions for inter-branch.

	\details
	- Parameters for global, inter-branch, and intra-branch variogram models are included.
	- Includes options to control neighborhood ranges and proportions for simulation.

*/
struct GeostatParams {
	bool is_used = false; //!< Flag indicating if geostatistical simulation is enabled.
	std::vector<float> simulated_property; //!< Simulated property values for the conduits.
	std::vector<float> simulation_distribution; //!< Input distribution used for the simulation.
	float global_vario_range; //!< Range parameter for the global variogram model.
	float global_range_of_neighborhood; //!< Range for neighborhood selection in global variogram.
	float global_vario_sill; //!< Sill of the junction-node variogram in the input-property space.
	float global_vario_nugget; //!< Nugget of the junction-node variogram in the input-property space.
	std::string global_vario_model; //!< Type of global variogram model (e.g., Spherical, Exponential, Gaussian).
	float interbranch_vario_range; //!< Range for inter-branch variogram model.
	float interbranch_range_of_neighborhood; //!< Neighborhood range for inter-branch model.
	float interbranch_vario_sill; //!< Sill of the inter-branch variogram in the input-property space.
	float interbranch_vario_nugget; //!< Nugget of the inter-branch variogram in the input-property space.
	std::string interbranch_vario_model; //!< Type of inter-branch variogram model.
	float intrabranch_vario_range; //!< Range for intra-branch variogram model.
	float intrabranch_range_of_neighborhood; //!< Neighborhood range for intra-branch model.
	float intrabranch_vario_sill; //!< Sill of the intra-branch variogram in the input-property space.
	float intrabranch_vario_nugget; //!< Nugget of the intra-branch variogram in the input-property space.
	std::string intrabranch_vario_model; //!< Type of intra-branch variogram model.
	int number_max_of_neighborhood_points; //!< Maximum number of points in a neighborhood used for covariance matrix.
	int nb_points_interbranch; //!< Number of points considered per branch for inter-branch variogram.
	float proportion_interbranch; //!< Proportion of points considered per branch for inter-branch variogram.
	bool use_drift_zwt = false; //!< Flag indicating if an external drift based on the distance to the water table should be enabled.
	bool use_drift_curv = false; //!< Flag indicating if an external drift based on an upstream/downstream increasing drain size trend should be enabled.
	std::vector<float> drift_output; //!< Simulated drift term for each node of the karst skeleton. Only filled if a drift is applied.
	std::vector<float> weights_output; //!< Simulated drift weight for each node of the karst skeleton. Only filled if a drift is applied.
};

/**
\brief Finds the neighborhood of a given node in a KarsticSkeleton graph based on curvilinear distances.

\param current_node_index Index of the current node.
\param curve Pointer to the KarsticSkeleton graph.
\param number_max_of_neighborhood_points Maximum number of points in the neighborhood.
\param range_of_neighborhood Range for the neighborhood selection.
\param type_neighborhood Type of neighborhood ("base", always adds node (used for Interbranch and Global steps), or "branch", only adds nodes within the same branch (used for Intrabranch step))
\param node_values Vector of property values associated with nodes in the graph.
\return A vector of indices representing the neighborhood nodes.
*/
std::vector<int> find_neighborhood(
	int current_node_index,
	const KarstNSim::KarsticSkeleton* const curve,
	int number_max_of_neighborhood_points,
	float range_of_neighborhood,
	std::string type_neighborhood,
	const std::vector<float>& node_values
);

/**
\brief Returns a variogram value at a given distance using variogram parameters and model type.

\param distance Distance at which the variogram is evaluated.
\param sill Sill parameter of the variogram.
\param nugget Nugget parameter of the variogram.
\param range Range parameter of the variogram.
\param vario_model Type of variogram model (e.g., Spherical, Exponential, Gaussian).
\return The computed variogram value.
*/
float variogram_value(
	const float& distance,
	const float& sill,
	const float& nugget,
	const float& range,
	const std::string& vario_model);

/**
\brief Inverts a given matrix.

\param input Matrix to be inverted.
\param output Inverted matrix is stored here.
\return True if the inversion is successful, false otherwise.
*/
bool invert_matrix(const std::vector<std::vector<float>>& input, std::vector<std::vector<float>>& output);

/**
\brief Computes shortest-path distances from a source node to a list of target nodes using a truncated Dijkstra search.

This helper function is used internally by kriging_in_point_on_the_fly() to avoid
relying on the global dense distance matrix. It computes the minimum path length
between the source and each target node in the graph, but only explores nodes within
a maximum distance threshold (range_cap). For unreachable nodes or distances beyond
the threshold, the returned distance is set to infinity.

\param src Index of the source node.
\param curve Pointer to the KarsticSkeleton containing the graph structure.
\param targets Indices of target nodes for which distances are required.
\param range_cap Maximum search radius; nodes beyond this distance are ignored.
\return Vector of distances (same order as targets). Values are infinity for unreachable nodes.
*/
std::vector<float> dijkstra_to_targets_truncated(
	int src,
	const KarstNSim::KarsticSkeleton* curve,
	const std::vector<int>& targets,
	float range_cap
);

/**
\brief Estimates the value and variance of a property at a given node using kriging,
	   computing all required distances on the fly from the graph structure.

This version removes the dependency on the dense distance matrix (mat_distance) and
instead computes local distances dynamically using truncated Dijkstra searches.
It provides identical kriging results as the standard function for the same
neighborhood and variogram parameters, but with much lower memory footprint
and better scalability for large graphs.

\param current_node_index Index of the current node.
\param neighborhood Indices of neighborhood nodes.
\param kriging_distribution Pointer to the random distribution used for SGS kriging.
\param curve Pointer to the KarsticSkeleton providing graph connectivity and node coordinates.
\param vario_range Range parameter of the variogram.
\param vario_sill Sill parameter of the variogram.
\param vario_nugget Nugget parameter of the variogram.
\param vario_model Type of variogram model ("Exponential", "Spherical", etc.).
\param node_values Values of the nodes in the graph (NDV = -99999).
\param[out] var_estimation Estimated variance (used for SGS simulation).
\param[out] val_estimation Estimated value (used for SGS simulation).
\param range_cap Maximum distance (graph path length) considered for neighborhood influence.
*/
void kriging_in_point_on_the_fly(
	int current_node_index,
	const std::vector<int>& neighborhood,
	const std::vector<float>* kriging_distribution,
	const KarstNSim::KarsticSkeleton* curve,
	const float& vario_range,
	const float& vario_sill,
	const float& vario_nugget,
	const std::string& vario_model,
	const std::vector<float>& node_values,
	float& var_estimation,
	float& val_estimation,
	float range_cap
);

/**
\brief Estimates the value and variance of a property at a given node using kriging.

\param current_node_index Index of the current node.
\param neighborhood Indices of neighborhood nodes.
\param kriging_distribution Pointer to the distribution used in kriging.
\param mat_distance Matrix of distances between nodes.
\param vario_range Range parameter of the variogram.
\param vario_sill Sill parameter of the variogram.
\param vario_nugget Nugget parameter of the variogram.
\param vario_model Type of variogram model.
\param node_values Values of the nodes in the graph.
\param[out] var_estimation Estimated variance (used for SGS simulation).
\param[out] val_estimation Estimated value (used for SGS simulation).
*/
void kriging_in_point(
	const int& current_node_index,
	const std::vector<int>& neighborhood,
	const std::vector<float>* kriging_distribution,
	const Array2D<float>& mat_distance,
	const float& vario_range,
	const float& vario_sill,
	const float& vario_nugget,
	const std::string& vario_model,
	const std::vector<float>& node_values, // node values, not updated here
	float& var_estimation, // variance estimated, used for SGS simulation
	float& val_estimation // value estimated, used for SGS simulation
);

/**
\brief Saves a vector of data to a file (useful for debugging).

\param data Vector of data to save.
\param filename Name of the output file.
*/
void save_data(const std::vector<float>& data, const std::string& filename);

/**
\brief Performs the SGS3/"branchwise" algorithm (3-variogram curvilinear graph-based SGS, Frantz et al., 2021).

\param curve Pointer to the KarsticSkeleton graph.
\param simulated_property Vector where simulated properties will be stored.
\param simulation_distribution Pointer to the distribution used for simulation,
expressed in the same property space as the conditioning data.
\param global_vario_range Range parameter of the global variogram.
\param global_range_of_neighborhood Range for neighborhood in the global variogram.
\param global_vario_sill Sill parameter of the global variogram.
\param global_vario_nugget Nugget parameter of the global variogram.
\param global_vario_model Type of global variogram model.
\param interbranch_vario_range Range parameter of the inter-branch variogram.
\param interbranch_range_of_neighborhood Range for neighborhood in the inter-branch variogram.
\param interbranch_vario_sill Sill parameter of the inter-branch variogram.
\param interbranch_vario_nugget Nugget parameter of the inter-branch variogram.
\param interbranch_vario_model Type of inter-branch variogram model.
\param intrabranch_vario_range Range parameter of the intra-branch variogram.
\param intrabranch_range_of_neighborhood Range for neighborhood in the intra-branch variogram.
\param intrabranch_vario_sill Sill parameter of the intra-branch variogram.
\param intrabranch_vario_nugget Nugget parameter of the intra-branch variogram.
\param intrabranch_vario_model Type of intra-branch variogram model.
\param number_max_of_neighborhood_points Maximum number of points in the neighborhood.
\param nb_points_interbranch Number of points considered per branch for inter-branch variogram.
\param proportion_interbranch Proportion of points considered per branch for inter-branch variogram.
*/
void SGS3(
	const KarstNSim::KarsticSkeleton* curve,
	std::vector<float>& simulated_property,
	const std::vector<float>* simulation_distribution,
	const float& global_vario_range,
	const float& global_range_of_neighborhood,
	const float& global_vario_sill,
	const float& global_vario_nugget,
	const std::string& global_vario_model,
	const float& interbranch_vario_range,
	const float& interbranch_range_of_neighborhood,
	const float& interbranch_vario_sill,
	const float& interbranch_vario_nugget,
	const std::string& interbranch_vario_model,
	const float& intrabranch_vario_range,
	const float& intrabranch_range_of_neighborhood,
	const float& intrabranch_vario_sill,
	const float& intrabranch_vario_nugget,
	const std::string& intrabranch_vario_model,
	const int& number_max_of_neighborhood_points,
	const int& nb_points_interbranch,
	const float& proportion_interbranch);


/**
\brief Computes the actual number of nodes in a branch that should be simulated in a given step.

\param total_nb_nodes_branch Total number of nodes in the branch.
\param nb_points_interbranch Number of points considered per branch for inter-branch variogram.
\param proportion_interbranch Proportion of points to consider for inter-branch variogram.
\return Number of nodes to be simulated in the branch.
*/
static int compute_prop_branch(const int& total_nb_nodes_branch, const int& nb_points_interbranch, const float& proportion_interbranch);

/**
\brief Computes the total upstream curvilinear length for each node in the karstic network.

This function recursively accumulates distances along all upstream branches of a given node
by traversing the branch ancestry defined in branch_id_ascend. It is used to capture the cumulative
position of the node in the network, from an upstream-downstream perspective.

\param curve Pointer to the KarsticSkeleton graph.
\return A vector of curvilinear lengths for each node in the graph.
*/
std::vector<float> compute_upstream_curvilinear_length(const KarstNSim::KarsticSkeleton* curve);

/**
\brief Computes an external drift field using a robust weighted linear regression.

The drift is modeled as a linear combination of user-specified explanatory variables:
- vertical distance above the water table (z - z_phreatic, clamped to zero),
- upstream curvilinear length.
Outliers are screened from leave-one-out residuals using the MAD (Median Absolute
Deviation) criterion. A unique inlet or outlet datum is protected from rejection.

\param curve Pointer to the KarsticSkeleton graph.
\param springs_xyz Outlet positions used to orient the curvilinear predictor and balance regression classes.
\param z_phreatic Reference water table elevation.
\param use_drift_zwt Whether to include the water table trend.
\param use_drift_curv Whether to include the upstream length trend.
\param eq_radius_values Observed equivalent radius values at each node.
\param conditioning_roles Hydraulic role of each conditioning datum.
\param[out] weights_out Regression weight for retained observations and zero otherwise.
\return External drift values m(u) for each node.
*/
std::vector<float> compute_external_drift(
	const KarstNSim::KarsticSkeleton* curve,
	const std::vector<Vector3>& springs_xyz,
	const float& z_phreatic,
	const bool& use_drift_zwt,
	const bool& use_drift_curv,
	const std::vector<float>& eq_radius_values,
	const std::vector<ConditioningDataRole>& conditioning_roles,
	std::vector<float>& weights_out);

/**
\brief Performs SGS simulation using an external drift field to account for large-scale trends.

This version enhances the classical branchwise SGS by decoupling large-scale deterministic trends
from local spatial variability. The drift m(u) is inferred by robust regression and added back
to the simulated residuals ε(u). Both the drift and the final weights are returned for diagnostics.

\param curve Pointer to the KarsticSkeleton graph.
\param springs_xyz Outlet positions used by the external-drift regression.
\param simulated_property Vector to store the final simulated values (drift + residuals).
\param conditioning_roles Hydraulic role of each conditioning datum.
\param simulation_distribution Pointer to the input distribution in the simulated-property space.
\param global_vario_range Range for the global variogram model.
\param global_range_of_neighborhood Neighborhood search radius for global step.
\param global_vario_sill Sill for the global variogram model.
\param global_vario_nugget Nugget for the global variogram model.
\param global_vario_model Type of global variogram (e.g., Spherical, Gaussian).
\param interbranch_vario_range Range for inter-branch variogram model.
\param interbranch_range_of_neighborhood Neighborhood range for inter-branch.
\param interbranch_vario_sill Sill for inter-branch variogram.
\param interbranch_vario_nugget Nugget for inter-branch variogram.
\param interbranch_vario_model Type of inter-branch variogram.
\param intrabranch_vario_range Range for intra-branch variogram model.
\param intrabranch_range_of_neighborhood Neighborhood range for intra-branch.
\param intrabranch_vario_sill Sill for intra-branch variogram.
\param intrabranch_vario_nugget Nugget for intra-branch variogram.
\param intrabranch_vario_model Type of intra-branch variogram.
\param number_max_of_neighborhood_points Maximum number of conditioning nodes for kriging.
\param nb_points_interbranch Maximum number of inter-branch points per branch.
\param proportion_interbranch Fraction of nodes to simulate in inter-branch step.
\param z_phreatic Constant elevation of the water table used in the drift.
\param use_drift_zwt Whether to activate the vadose/phreatic trend in the drift.
\param use_drift_curv Whether to activate the upstream length trend in the drift.
\param[out] drift_output Vector storing the external drift values m(u) for each node.
\param[out] weights_output Vector storing the weights assigned to each observation during regression.
*/
void SGS3_with_external_drift(
	const KarstNSim::KarsticSkeleton* curve,
	const std::vector<Vector3>& springs_xyz,
	std::vector<float>& simulated_property,
	const std::vector<ConditioningDataRole>& conditioning_roles,
	const std::vector<float>* simulation_distribution,
	const float& global_vario_range,
	const float& global_range_of_neighborhood,
	const float& global_vario_sill,
	const float& global_vario_nugget,
	const std::string& global_vario_model,
	const float& interbranch_vario_range,
	const float& interbranch_range_of_neighborhood,
	const float& interbranch_vario_sill,
	const float& interbranch_vario_nugget,
	const std::string& interbranch_vario_model,
	const float& intrabranch_vario_range,
	const float& intrabranch_range_of_neighborhood,
	const float& intrabranch_vario_sill,
	const float& intrabranch_vario_nugget,
	const std::string& intrabranch_vario_model,
	const int& number_max_of_neighborhood_points,
	const int& nb_points_interbranch,
	const float& proportion_interbranch,
	const float& z_phreatic,
	const bool& use_drift_zwt,
	const bool& use_drift_curv,
	std::vector<float>& drift_output,
	std::vector<float>& weights_output);
