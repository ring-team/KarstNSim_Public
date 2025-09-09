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

#include "KarstNSim/graph.h"
#include "KarstNSim/randomgenerator.h"

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
	float global_vario_sill; //!< Sill parameter for the global variogram model.
	float global_vario_nugget; //!< Nugget parameter for the global variogram model.
	std::string global_vario_model; //!< Type of global variogram model (e.g., Spherical, Exponential, Gaussian).
	float interbranch_vario_range; //!< Range for inter-branch variogram model.
	float interbranch_range_of_neighborhood; //!< Neighborhood range for inter-branch model.
	float interbranch_vario_sill; //!< Sill parameter for inter-branch variogram.
	float interbranch_vario_nugget; //!< Nugget parameter for inter-branch variogram.
	std::string interbranch_vario_model; //!< Type of inter-branch variogram model.
	float intrabranch_vario_range; //!< Range for intra-branch variogram model.
	float intrabranch_range_of_neighborhood; //!< Neighborhood range for intra-branch model.
	float intrabranch_vario_sill; //!< Sill parameter for intra-branch variogram.
	float intrabranch_vario_nugget; //!< Nugget parameter for intra-branch variogram.
	std::string intrabranch_vario_model; //!< Type of intra-branch variogram model.
	int number_max_of_neighborhood_points; //!< Maximum number of points in a neighborhood used for covariance matrix.
	int nb_points_interbranch; //!< Number of points considered per branch for inter-branch variogram.
	float proportion_interbranch; //!< Proportion of points considered per branch for inter-branch variogram.
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
\param simulation_distribution Pointer to the distribution used for simulation.
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
	std::vector<float>&  simulated_property,
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


