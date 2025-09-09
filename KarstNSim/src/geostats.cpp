/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for new methods + modifications to original methods
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Rewritten entirely (except the variogram_value fonction, written by G. Rongier in 2015) from the pseudo-code algorithm (Algorithm 2) in the paper of Frantz et al., 2021 "Analysis and stochastic simulation of geometrical properties of conduits in karstic networks".
This work was performed in the frame of the RING project at Université de Lorraine.

***************************************************************/

#include "KarstNSim/geostats.h"

std::vector<int> find_neighborhood(
	int current_node_index,
	const KarstNSim::KarsticSkeleton* const curve,
	int number_max_of_neighborhood_points,
	float range_of_neighborhood,
	std::string type_neighborhood,
	const std::vector<float>& node_values
) {

	int node_id = curve->nodes.at(current_node_index).branch_id;
	const Array2D<float>& distance_mat = curve->distance_mat;
	int num_nodes = int(curve->nodes.size());
	std::vector<float> distances(num_nodes, std::numeric_limits<float>::max());

	std::vector<int> neighborhood;

	// Calculate distances from the current node to all other nodes
	for (int i = 0; i < num_nodes; ++i) {
		distances[i] = distance_mat(current_node_index, i);
	}

	// Find the closest neighbors within the specified range
	for (int i = 0; i < num_nodes; ++i) {
		if (i != current_node_index && distances[i] <= range_of_neighborhood) { // check that the other node is not too far
			if ((node_values[i] - (-99999)) > 1e-12) { // if node is not ndv (-99999)
				if (type_neighborhood == "base") { // in the case of the base neighborhood search, always add node (used for Interbranch and Global steps)
					neighborhood.push_back(i);
				}
				else if (type_neighborhood == "branch") { // this type of search only uses nodes within the same branch (used for Intrabranch step)
					int other_node_id = curve->nodes.at(i).branch_id;
					if (other_node_id == node_id) { // same branch
						neighborhood.push_back(i);
					}
				}
			}

		}
	}

	// Sort the neighborhood based on distances
	std::sort(neighborhood.begin(), neighborhood.end(), [&](int a, int b) {
		return distances[a] < distances[b];
	});

	// Keep only the closest number_max_of_neighborhood_points neighbors
	if (neighborhood.size() > number_max_of_neighborhood_points) {
		neighborhood.resize(number_max_of_neighborhood_points);
	}
	return neighborhood;
}

// variogram_value (Rongier-2015)
float variogram_value(
	const float& distance,
	const float& sill,
	const float& nugget,
	const float& range,
	const std::string& vario_model) {

	float vario_value = -1;
	// Check the input parameters.

	if (nugget >= 0. && sill >= nugget && range >= 0. && vario_model != "") {
		// Get the variogram value following the model of variogram.
		if (vario_model == "Gaussian") {
			vario_value =
				nugget +
				(sill - nugget)*
				(1 - exp(-(3 * distance / range)*(distance / range)));
		}
		else if (vario_model == "Spherical") {
			if (distance <= range) {
				vario_value =
					nugget +
					(sill - nugget)*
					((3 / 2)*(distance / range) - (1 / 2)*pow((distance / range), 3.));
			}
			else if (distance > range) {
				vario_value = sill;
			}
		}
		else if (vario_model == "Exponential") {
			vario_value =
				nugget +
				(sill - nugget)*
				(1 - exp(-3 * distance / range));
		}
		else if (vario_model == "Nugget") {
			vario_value = nugget;
		}
	}
	return vario_value;
}

// Function to calculate the inverse of a square matrix using Gaussian elimination
bool invert_matrix(const std::vector<std::vector<float>>& input, std::vector<std::vector<float>>& output) {
	// Check if the input matrix is square
	int size = int(input.size());
	if (size != input[0].size()) {
		std::cerr << "Input matrix is not square." << std::endl;
		return false;
	}

	// Initialize the output matrix as the identity matrix
	output = std::vector<std::vector<float>>(size, std::vector<float>(size, 0.0));
	for (int i = 0; i < size; ++i) {
		output[i][i] = 1.0;
	}

	// Copy the input matrix to avoid modifying the original
	std::vector<std::vector<float>> A = input;

	// Gaussian elimination with partial pivoting
	for (int i = 0; i < size; ++i) {
		// Find the pivot row
		int max_row = i;
		for (int k = i + 1; k < size; ++k) {
			if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
				max_row = k;
			}
		}

		// Swap the current row with the pivot row
		if (max_row != i) {
			A[i].swap(A[max_row]);
			output[i].swap(output[max_row]);
		}

		// Check if the matrix is singular
		if (A[i][i] == 0.0) {
			std::cerr << "Matrix is singular." << std::endl;
			return false;
		}

		// Scale the current row to make the diagonal element 1
		float scale = 1.0 / A[i][i];
		for (int j = 0; j < size; ++j) {
			A[i][j] *= scale;
			output[i][j] *= scale;
		}

		// Eliminate the current column
		for (int k = 0; k < size; ++k) {
			if (k != i) {
				float factor = A[k][i];
				for (int j = 0; j < size; ++j) {
					A[k][j] -= factor * A[i][j];
					output[k][j] -= factor * output[i][j];
				}
			}
		}
	}
	return true;
}

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
) {

	float average = std::accumulate(kriging_distribution->begin(), kriging_distribution->end(), 0.0) / kriging_distribution->size();

	val_estimation = -99999;

	// If the neighborhood is empty, perform a random draw in the initial distribution.
	if (neighborhood.empty()) {
		val_estimation = select_random_element(*kriging_distribution);
		return;
	}
	else {
		int nb_neigh = int(neighborhood.size());

		// Create matrix of distances that can be used for calculation
		std::vector<std::vector<float>> local_mat_distance(nb_neigh + 1, std::vector<float>(nb_neigh + 1));
		std::vector<std::vector<float>> local_mat_vario_values(nb_neigh + 1, std::vector<float>(nb_neigh + 1));
		std::vector<std::vector<float>> local_mat_cov(nb_neigh + 1, std::vector<float>(nb_neigh + 1));

		local_mat_distance[0][0] = 0.;
		for (int i = 0; i < nb_neigh; ++i) {
			local_mat_distance[i + 1][0] = mat_distance(current_node_index, neighborhood[i]);
			local_mat_distance[0][i + 1] = mat_distance(current_node_index, neighborhood[i]);
		}

		// Fill the rest of the matrix
		for (int i = 0; i < nb_neigh; ++i) {
			for (int j = 0; j < nb_neigh; ++j) {
				local_mat_distance[i + 1][j + 1] = mat_distance(neighborhood[i], neighborhood[j]);
			}
		}

		float a = vario_range;
		float C = vario_sill;
		float C0 = vario_nugget;

		// Fill matrix of variogram values
		for (int j = 0; j < nb_neigh + 1; ++j) {
			float h = local_mat_distance[0][j];
			if (h == 0.) { // It should not be the case
				local_mat_vario_values[0][j] = 0.;
			}
			else {
				local_mat_vario_values[0][j] = variogram_value(h, C, C0, a, vario_model);
			}
			local_mat_vario_values[j][0] = local_mat_vario_values[0][j];
		}

		for (int i = 1; i < nb_neigh + 1; ++i) {
			for (int j = i; j < nb_neigh + 1; ++j) {
				float h = local_mat_distance[i][j]; // Matrix of values in variogram
				if (h == 0.) {
					local_mat_vario_values[i][j] = 0.;
				}
				else {
					local_mat_vario_values[i][j] = variogram_value(h, C, C0, a, vario_model);
				}
				local_mat_vario_values[j][i] = local_mat_vario_values[i][j];
			}
		}

		// Matrix of covariance
		for (int i = 0; i < nb_neigh + 1; ++i) {
			for (int j = 0; j < nb_neigh + 1; ++j) {
				local_mat_cov[i][j] = C - local_mat_vario_values[i][j];
			}
		}

		// Perform simple kriging
		std::vector<std::vector<float>> k0(nb_neigh, std::vector<float>(1));

		for (int i = 0; i < nb_neigh; ++i) {
			k0[i][0] = local_mat_cov[0][i + 1];
		}

		std::vector<std::vector<float>> K(nb_neigh, std::vector<float>(nb_neigh));

		for (int i = 0; i < nb_neigh; ++i) {
			for (int j = 0; j < nb_neigh; ++j) {
				K[i][j] = local_mat_cov[i + 1][j + 1];
			}
		}

		std::vector<std::vector<float>> lambda(nb_neigh, std::vector<float>(1));
		std::vector<std::vector<float>> inv_K(nb_neigh, std::vector<float>(nb_neigh));

		bool invert_success = invert_matrix(K, inv_K);

		if (!invert_success) {
			// Handle inversion failure
			return;
		}

		for (int i = 0; i < nb_neigh; ++i) {
			lambda[i][0] = 0;
			for (int j = 0; j < nb_neigh; ++j) {
				lambda[i][0] += inv_K[i][j] * k0[j][0];
			}
		}

		val_estimation = average;
		float sum_lambda = 0;

		for (int i = 0; i < nb_neigh; ++i) {
			val_estimation += lambda[i][0] * node_values[neighborhood[i]];
			sum_lambda += lambda[i][0];
		}

		val_estimation -= sum_lambda * average;

		float var_temp = 0;

		for (int i = 0; i < nb_neigh; ++i) {
			var_temp += k0[i][0] * lambda[i][0];
		}

		var_estimation = local_mat_cov[0][0] - var_temp;

		return;
	}
}



void save_data(const std::vector<float>& data, const std::string& filename) {
	std::ofstream file(filename);

	for (const auto& value : data) {
		file << value << std::endl;
	}

	file.close();

}

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
	const float& proportion_interbranch) {

	// 1) Perform Normal Score Transform of initial distrib AND initial data vector
	std::vector<float> simulated_prop_gauss;
	std::vector<float> simulated_property_copy = simulated_property;
	if (!simulated_property.empty()) {
		simulated_prop_gauss = nst_data_with_nodata(simulated_property, *simulation_distribution, 0., 1.); // gaussianize the data values if any (do not change the no data values though)
	}
	else {
		simulated_prop_gauss.resize(curve->nodes.size());
		std::fill(simulated_prop_gauss.begin(), simulated_prop_gauss.end(), -99999);
	}

	std::vector<float> sim_distrib_gauss(simulation_distribution->size());
	if (!simulation_distribution->empty()) {
		sim_distrib_gauss = nst(*simulation_distribution, 0., 1.);
	}

	//save_data(*simulation_distribution, "distrib.txt");
	//save_data(sim_distrib_gauss, "gauss_distrib.txt");

	// 2) Interbranch simulation with interbranch variogram

	// iterate on branches:
	for (int branch_id = 0; branch_id < curve->branch_sizes.size(); branch_id++) {

		// find number of nodes to simulate in given branch
		int nb_nodes_to_simulate = compute_prop_branch(curve->branch_sizes.at(branch_id), nb_points_interbranch, proportion_interbranch);

		std::vector<int> all_branch_nodes;
		// create list of indices of nodes of that branch in nodes, and shuffle that list to get the random order of iteration
		for (int i = 0; i < curve->nodes.size(); i++) {
			if (curve->nodes.at(i).branch_id == branch_id && (simulated_prop_gauss[i] - (-99999)) < 1e-12) { // check that the node is in the same branch and isnt assigned to a value yet (no data value)
				all_branch_nodes.push_back(i);
			}
		}
		std::vector<int> nodes_to_simulate = select_random_elements(all_branch_nodes, nb_nodes_to_simulate);

		for (int i = 0; i < nodes_to_simulate.size(); i++) {

			// determine neighborhood (conditional nodes nearby)

			std::vector<int> neighbors_interbranch = find_neighborhood(nodes_to_simulate.at(i), curve, number_max_of_neighborhood_points, interbranch_range_of_neighborhood, "base", simulated_prop_gauss);

			// if empty neighborhood, simulate directly by sampling from distrib, otherwise, simulate value (kriging + sample from gaussian estimate)
			float val_estimation;
			float var_estimation;

			kriging_in_point(nodes_to_simulate.at(i), neighbors_interbranch, &sim_distrib_gauss, curve->distance_mat, interbranch_vario_range, interbranch_vario_sill, interbranch_vario_nugget, interbranch_vario_model, simulated_prop_gauss, var_estimation, val_estimation);

			if (!neighbors_interbranch.empty()) { // if we had a neighborhood, we simulate value from estimated (kriged) value
				var_estimation = (var_estimation < 0.) ? 0. : var_estimation; // due to numerical uncertainties, var estimation can sometimes be slightly negative.
				float simulated_val = generateNormalRandom(val_estimation, std::sqrt(var_estimation));
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = simulated_val;
			}
			else { // if the neighborhood was empty, we simply sampled a value from the distribution
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = val_estimation;
			}
		}
	}

	// 3) Intrabranch simulation with intrabranch variogram

	for (int branch_id = 0; branch_id < curve->branch_sizes.size(); branch_id++) {

		int nb_nodes_to_simulate = compute_prop_branch(curve->branch_sizes.at(branch_id), curve->branch_sizes.at(branch_id), 1.); // compute ALL remaining nodes

		std::vector<int> all_branch_nodes;
		// create list of indices of nodes of that branch in nodes, and shuffle that list to get the random order of iteration
		for (int i = 0; i < curve->nodes.size(); i++) {
			if (curve->nodes.at(i).branch_id == branch_id && (simulated_prop_gauss[i] - (-99999)) < 1e-12) { // check that the node is in the same branch and isnt assigned to a value yet (no data value)
				all_branch_nodes.push_back(i);
			}
		}
		std::vector<int> nodes_to_simulate = select_random_elements(all_branch_nodes, nb_nodes_to_simulate);

		for (int i = 0; i < nodes_to_simulate.size(); i++) {
			std::vector<int> neighbors_intrabranch = find_neighborhood(nodes_to_simulate.at(i), curve, number_max_of_neighborhood_points, intrabranch_range_of_neighborhood, "branch", simulated_prop_gauss);

			float val_estimation;
			float var_estimation;

			kriging_in_point(nodes_to_simulate.at(i), neighbors_intrabranch, &sim_distrib_gauss, curve->distance_mat, intrabranch_vario_range, intrabranch_vario_sill, intrabranch_vario_nugget, intrabranch_vario_model, simulated_prop_gauss, var_estimation, val_estimation);

			if (!neighbors_intrabranch.empty()) {
				var_estimation = (var_estimation < 0.) ? 0. : var_estimation;// due to numerical uncertainties, var estimation can sometimes be slightly negative.
				float simulated_val = generateNormalRandom(val_estimation, std::sqrt(var_estimation));
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = simulated_val;
			}
			else {
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = val_estimation;
			}
		}
	}

	// 4) Intersections simulation with global variogram

	std::vector<int> all_intersection_nodes;
	for (int i = 0; i < curve->nodes.size(); i++) {
		if (curve->nodes.at(i).branch_id == -1 && (simulated_prop_gauss[i] - (-99999)) < 1e-12) // if intersection and value not already assigned
			all_intersection_nodes.push_back(i);
	}
	std::vector<int> nodes_to_simulate_global = select_random_elements(all_intersection_nodes, int(all_intersection_nodes.size())); // we compute ALL intersections, no exceptions

	for (int i = 0; i < nodes_to_simulate_global.size(); i++) {

		std::vector<int> neighbors_global = find_neighborhood(nodes_to_simulate_global.at(i), curve, number_max_of_neighborhood_points, global_range_of_neighborhood, "base", simulated_prop_gauss);

		float val_estimation_global;
		float var_estimation_global;

		kriging_in_point(nodes_to_simulate_global.at(i), neighbors_global, &sim_distrib_gauss, curve->distance_mat, global_vario_range, global_vario_sill, global_vario_nugget, global_vario_model, simulated_prop_gauss, var_estimation_global, val_estimation_global);

		if (!neighbors_global.empty()) {
			var_estimation_global = (var_estimation_global < 0.) ? 0. : var_estimation_global;// due to numerical uncertainties, var estimation can sometimes be slightly negative.
			float simulated_val_global = generateNormalRandom(val_estimation_global, std::sqrt(var_estimation_global));
			simulated_prop_gauss.at(nodes_to_simulate_global.at(i)) = simulated_val_global;
		}
		else {
			simulated_prop_gauss.at(nodes_to_simulate_global.at(i)) = val_estimation_global;
		}
	}

	// 5) Back-transform of generated distribution.

	// if we had data at the beginning, use it as basis for back transformation
	//if (!simulated_property.empty()) {
	//	simulated_property = back_transform(simulated_prop_gauss, simulated_property_copy);
	//}
	//else { // else, use the initial distribution instead
	simulated_property = back_transform(simulated_prop_gauss, *simulation_distribution);
	//}
}

int compute_prop_branch(const int& total_nb_nodes_branch, const int& nb_points_interbranch, const float& proportion_interbranch) {

	float true_proportion_interbranch;
	int needed_value_points_for_this_branch;

	if (proportion_interbranch > 1.) {
		true_proportion_interbranch = 1.;
	}
	else {
		true_proportion_interbranch = proportion_interbranch;
	}

	int proportion_value = round(total_nb_nodes_branch * true_proportion_interbranch);

	if (proportion_value <= nb_points_interbranch) {
		needed_value_points_for_this_branch = proportion_value;
	}
	else {
		needed_value_points_for_this_branch = nb_points_interbranch;
	}

	return needed_value_points_for_this_branch;
}