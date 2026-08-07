/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/parse_inputs.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>


namespace {
	struct InstructionOccurrence {
		std::size_t line_number;
		std::string value;
	};

	using InstructionIndex = std::unordered_map<std::string, InstructionOccurrence>;

	std::string trim_copy(const std::string& input) {
		const std::size_t first = input.find_first_not_of(" \t\r\n");
		if (first == std::string::npos) {
			return std::string();
		}
		const std::size_t last = input.find_last_not_of(" \t\r\n");
		return input.substr(first, last - first + 1);
	}

	bool is_instruction_comment(const std::string& trimmed_line) {
		return trimmed_line.compare(0, 2, "//") == 0 ||
			trimmed_line.compare(0, 2, "::") == 0 ||
			trimmed_line.compare(0, 1, "#") == 0 ||
			trimmed_line.compare(0, 1, "%") == 0;
	}

	const std::unordered_set<std::string>& known_parameter_tags() {
		static const std::unordered_set<std::string> tags = {
			"add_inception_surfaces:",
			"allow_single_outlet_connection:",
			"alteration_lines:",
			"create_grid:",
			"create_nghb_graph:",
			"create_nghb_graph_property:",
			"create_solved_connectivity_matrix:",
			"create_vset_sampling:",
			"domain:",
			"fraction_karst_perm:",
			"fraction_old_karst_perm:",
			"fracture_constraint_weight:",
			"fracture_families_orientations:",
			"fracture_families_tolerance:",
			"gamma:",
			"ghost_rock_weight:",
			"ghostrock_max_vertical_size:",
			"ghostrock_width:",
			"global_range_of_neighborhood:",
			"global_vario_model:",
			"global_vario_nugget:",
			"global_vario_range:",
			"global_vario_sill:",
			"gradient_constraint_weight:",
			"inception_surface_constraint_weight:",
			"inception_surfaces:",
			"input_nghb_graph:",
			"interbranch_range_of_neighborhood:",
			"interbranch_vario_model:",
			"interbranch_vario_nugget:",
			"interbranch_vario_range:",
			"interbranch_vario_sill:",
			"interpolate_lines:",
			"intrabranch_range_of_neighborhood:",
			"intrabranch_vario_model:",
			"intrabranch_vario_nugget:",
			"intrabranch_vario_range:",
			"intrabranch_vario_sill:",
			"k_pts:",
			"karstic_network_name:",
			"karstification_potential_weight:",
			"main_repository:",
			"max_depth_horizon:",
			"max_distance_amplification:",
			"max_distance_of_deadend_pts:",
			"max_inception_surface_distance:",
			"min_distance_amplification:",
			"multiply_costs:",
			"nb_cycles:",
			"nb_deadend_points:",
			"nb_points_interbranch:",
			"nghb_count:",
			"nghb_radius:",
			"noise_frequency:",
			"noise_octaves:",
			"noise_weight:",
			"number_max_of_neighborhood_points:",
			"number_of_iterations:",
			"outlet_selection_cost_factor:",
			"poisson_radius:",
			"previous_networks:",
			"proportion_interbranch:",
			"refine_surface_sampling:",
			"sampling_points:",
			"sections_simulation_only:",
			"selected_seed:",
			"simulate_sections:",
			"simulation_distribution:",
			"simulation_input_dir:",
			"sinks:",
			"sphere_centers:",
			"springs:",
			"surf_wat_table:",
			"topo_surface:",
			"use_cycle_amplification:",
			"use_deadend_points:",
			"use_density_property:",
			"use_drift_curv:",
			"use_drift_zwt:",
			"use_fracture_constraints:",
			"use_ghostrocks:",
			"use_input_nghb_graph:",
			"use_karstification_potential:",
			"use_max_depth_constraint:",
			"use_max_nghb_radius:",
			"use_no_karst_spheres:",
			"use_noise:",
			"use_noise_on_all:",
			"use_previous_networks:",
			"use_sampling_points:",
			"use_sinks_radius:",
			"use_springs_radius:",
			"use_user_connectivity_matrix:",
			"use_waypoints:",
			"use_waypoints_radius:",
			"vadose_cohesion:",
			"vary_seed:",
			"vertical_distance_stretching_factor:",
			"water_table_constraint_weight_phreatic:",
			"water_table_constraint_weight_vadose:",
			"waypoints:",
			"waypoints_weight:"
		};
		return tags;
	}

	bool parse_boolean_token_strict(
		const std::string& raw_value,
		const std::string& parameter_name,
		std::size_t line_number) {

		std::istringstream stream(raw_value);
		std::string token;
		std::string extra;
		stream >> token;
		if (token.empty()) {
			throw std::runtime_error(
				"[instructions] Missing boolean value for '" + parameter_name +
				"' at line " + std::to_string(line_number) + "."
			);
		}
		if (stream >> extra) {
			throw std::runtime_error(
				"[instructions] Unexpected extra token after boolean parameter '" +
				parameter_name + "' at line " + std::to_string(line_number) + "."
			);
		}
		if (token == "true" || token == "1") {
			return true;
		}
		if (token == "false" || token == "0") {
			return false;
		}
		throw std::runtime_error(
			"[instructions] Invalid boolean value '" + token + "' for parameter '" +
			parameter_name + "' at line " + std::to_string(line_number) +
			". Accepted values are true, false, 1, and 0."
		);
	}

	InstructionIndex scan_instruction_file(std::ifstream& input_file) {
		InstructionIndex index;
		std::string line;
		std::size_t line_number = 0;

		while (std::getline(input_file, line)) {
			++line_number;
			const std::string trimmed_line = trim_copy(line);
			if (trimmed_line.empty() || is_instruction_comment(trimmed_line)) {
				continue;
			}

			std::istringstream line_stream(trimmed_line);
			std::string parameter_name;
			line_stream >> parameter_name;

			if (known_parameter_tags().count(parameter_name) == 0) {
				const std::string colon_candidate = parameter_name + ":";
				if (!parameter_name.empty() && parameter_name.back() != ':' &&
					known_parameter_tags().count(colon_candidate) != 0) {
					throw std::runtime_error(
						"[instructions] Invalid declaration at line " +
						std::to_string(line_number) + ": parameter '" + parameter_name +
						"' must be followed by ':' (expected '" + colon_candidate + " value')."
					);
				}
				throw std::runtime_error(
					"[instructions] Unknown parameter '" + parameter_name +
					"' at line " + std::to_string(line_number) + "."
				);
			}

			std::string raw_value;
			std::getline(line_stream, raw_value);
			raw_value = trim_copy(raw_value);
			if (raw_value.empty()) {
				throw std::runtime_error(
					"[instructions] Parameter '" + parameter_name +
					"' has no value at line " + std::to_string(line_number) + "."
				);
			}

			const auto previous = index.find(parameter_name);
			if (previous != index.end()) {
				throw std::runtime_error(
					"[instructions] Parameter '" + parameter_name + "' is defined more than once "
					"(lines " + std::to_string(previous->second.line_number) + " and " +
					std::to_string(line_number) + ")."
				);
			}

			index.emplace(parameter_name, InstructionOccurrence{ line_number, raw_value });
		}

		input_file.clear();
		input_file.seekg(0, std::ios::beg);
		return index;
	}

	bool has_parameter(const InstructionIndex& index, const std::string& name) {
		return index.find(name) != index.end();
	}

	bool boolean_parameter(
		const InstructionIndex& index,
		const std::string& name,
		bool default_value) {
		const auto it = index.find(name);
		if (it == index.end()) {
			return default_value;
		}
		return parse_boolean_token_strict(it->second.value, name, it->second.line_number);
	}

	std::string first_parameter_token(
		const InstructionIndex& index,
		const std::string& name,
		const std::string& default_value = std::string()) {
		const auto it = index.find(name);
		if (it == index.end()) {
			return default_value;
		}
		std::istringstream stream(it->second.value);
		std::string token;
		stream >> token;
		return token;
	}


	void validate_single_token(
		const InstructionOccurrence& occurrence,
		const std::string& parameter_name) {
		std::istringstream stream(occurrence.value);
		std::string value;
		std::string extra;
		if (!(stream >> value)) {
			throw std::runtime_error(
				"[instructions] Missing value for parameter '" + parameter_name +
				"' at line " + std::to_string(occurrence.line_number) + "."
			);
		}
		if (stream >> extra) {
			throw std::runtime_error(
				"[instructions] Parameter '" + parameter_name + "' at line " +
				std::to_string(occurrence.line_number) +
				" accepts exactly one value, but additional token '" + extra + "' was found."
			);
		}
	}

	template <typename NumericType>
	void validate_single_numeric_value(
		const InstructionOccurrence& occurrence,
		const std::string& parameter_name,
		const std::string& expected_type) {
		std::istringstream stream(occurrence.value);
		NumericType value{};
		if (!(stream >> value)) {
			throw std::runtime_error(
				"[instructions] Invalid " + expected_type + " value for parameter '" +
				parameter_name + "' at line " + std::to_string(occurrence.line_number) + "."
			);
		}
		stream >> std::ws;
		if (!stream.eof()) {
			throw std::runtime_error(
				"[instructions] Invalid trailing character or extra value for parameter '" +
				parameter_name + "' at line " + std::to_string(occurrence.line_number) +
				". Expected exactly one " + expected_type + " value."
			);
		}
	}

	void validate_float_list(
		const InstructionOccurrence& occurrence,
		const std::string& parameter_name) {
		std::istringstream stream(occurrence.value);
		bool found_value = false;
		while (true) {
			stream >> std::ws;
			if (stream.eof()) {
				break;
			}
			float value = 0.0f;
			if (!(stream >> value)) {
				throw std::runtime_error(
					"[instructions] Invalid floating-point list for parameter '" +
					parameter_name + "' at line " + std::to_string(occurrence.line_number) +
					". Values must be separated by whitespace and must not contain trailing punctuation."
				);
			}
			found_value = true;
		}
		if (!found_value) {
			throw std::runtime_error(
				"[instructions] Parameter '" + parameter_name + "' at line " +
				std::to_string(occurrence.line_number) + " requires at least one value."
			);
		}
	}

	void validate_token_list(
		const InstructionOccurrence& occurrence,
		const std::string& parameter_name) {
		std::istringstream stream(occurrence.value);
		std::string token;
		if (!(stream >> token)) {
			throw std::runtime_error(
				"[instructions] Parameter '" + parameter_name + "' at line " +
				std::to_string(occurrence.line_number) + " requires at least one value."
			);
		}
	}

	void validate_instruction_value_syntax(const InstructionIndex& index) {
		static const std::unordered_set<std::string> boolean_parameters = {
			"vary_seed:", "use_sampling_points:", "use_density_property:",
			"use_max_nghb_radius:", "use_springs_radius:", "use_sinks_radius:",
			"use_waypoints_radius:", "allow_single_outlet_connection:",
			"use_user_connectivity_matrix:", "use_waypoints:", "use_ghostrocks:",
			"interpolate_lines:", "use_max_depth_constraint:", "add_inception_surfaces:",
			"use_karstification_potential:", "use_fracture_constraints:",
			"use_previous_networks:", "sections_simulation_only:",
			"use_input_nghb_graph:", "use_no_karst_spheres:", "use_deadend_points:",
			"use_cycle_amplification:", "use_noise:", "use_noise_on_all:",
			"simulate_sections:", "use_drift_zwt:", "use_drift_curv:",
			"create_vset_sampling:", "create_nghb_graph:",
			"create_nghb_graph_property:", "create_solved_connectivity_matrix:",
			"create_grid:", "multiply_costs:", "vadose_cohesion:"
		};

		static const std::unordered_set<std::string> integer_parameters = {
			"selected_seed:", "number_of_iterations:", "k_pts:", "nghb_count:",
			"refine_surface_sampling:", "nb_deadend_points:", "nb_cycles:",
			"noise_frequency:", "noise_octaves:",
			"number_max_of_neighborhood_points:", "nb_points_interbranch:"
		};

		static const std::unordered_set<std::string> floating_parameters = {
			"poisson_radius:", "nghb_radius:", "waypoints_weight:",
			"ghostrock_max_vertical_size:", "ghost_rock_weight:", "ghostrock_width:",
			"inception_surface_constraint_weight:", "max_inception_surface_distance:",
			"karstification_potential_weight:", "fracture_constraint_weight:",
			"fraction_old_karst_perm:", "water_table_constraint_weight_vadose:",
			"water_table_constraint_weight_phreatic:", "max_distance_of_deadend_pts:",
			"max_distance_amplification:", "min_distance_amplification:",
			"noise_weight:", "global_vario_range:", "global_range_of_neighborhood:",
			"global_vario_sill:", "global_vario_nugget:",
			"interbranch_vario_range:", "interbranch_range_of_neighborhood:",
			"interbranch_vario_sill:", "interbranch_vario_nugget:",
			"intrabranch_vario_range:", "intrabranch_range_of_neighborhood:",
			"intrabranch_vario_sill:", "intrabranch_vario_nugget:",
			"proportion_interbranch:", "fraction_karst_perm:", "gamma:",
			"vertical_distance_stretching_factor:", "gradient_constraint_weight:",
			"outlet_selection_cost_factor:"
		};

		static const std::unordered_set<std::string> floating_list_parameters = {
			"fracture_families_orientations:", "fracture_families_tolerance:"
		};

		static const std::unordered_set<std::string> token_list_parameters = {
			"inception_surfaces:", "previous_networks:", "surf_wat_table:"
		};

		for (const auto& item : index) {
			const std::string& parameter_name = item.first;
			const InstructionOccurrence& occurrence = item.second;

			if (boolean_parameters.count(parameter_name) != 0) {
				parse_boolean_token_strict(
					occurrence.value, parameter_name, occurrence.line_number);
			}
			else if (integer_parameters.count(parameter_name) != 0) {
				validate_single_numeric_value<int>(
					occurrence, parameter_name, "integer");
			}
			else if (floating_parameters.count(parameter_name) != 0) {
				validate_single_numeric_value<float>(
					occurrence, parameter_name, "floating-point");
			}
			else if (floating_list_parameters.count(parameter_name) != 0) {
				validate_float_list(occurrence, parameter_name);
			}
			else if (token_list_parameters.count(parameter_name) != 0) {
				validate_token_list(occurrence, parameter_name);
			}
			else {
				validate_single_token(occurrence, parameter_name);
			}
		}
	}


	void validate_required_parameter_presence(const InstructionIndex& index) {
		std::vector<std::string> missing;
		std::unordered_set<std::string> already_reported;
		auto require = [&](const std::string& name, const std::string& reason) {
			if (!has_parameter(index, name) && already_reported.insert(name).second) {
				missing.push_back(name + " — " + reason);
			}
		};

		// Parameters that define the global execution context are required in every workflow.
		require("karstic_network_name:", "required for output naming");
		require("main_repository:", "required to resolve input and output paths");
		require("selected_seed:", "required for reproducibility");
		require("number_of_iterations:", "required to define the number of realizations");
		require("vary_seed:", "required to define multi-realization seed behavior");
		require("sections_simulation_only:", "required to select the simulation workflow");
		require("simulate_sections:", "required to explicitly enable or disable conduit-size simulation");
		require("use_previous_networks:", "required to explicitly enable or disable previous-network reuse");

		const bool sections_only = boolean_parameter(index, "sections_simulation_only:", false);
		const bool simulate_sections = boolean_parameter(index, "simulate_sections:", false);
		const bool use_previous_networks = boolean_parameter(index, "use_previous_networks:", false);
		const bool use_waypoints = boolean_parameter(index, "use_waypoints:", false);
		const bool use_ghostrocks = boolean_parameter(index, "use_ghostrocks:", false);
		const bool use_karstification_potential = boolean_parameter(
			index, "use_karstification_potential:", false);
		const bool use_sinks_radius = boolean_parameter(
			index, "use_sinks_radius:", false);
		const bool use_springs_radius = boolean_parameter(
			index, "use_springs_radius:", false);
		const bool use_drift_zwt = boolean_parameter(
			index, "use_drift_zwt:", false);
		const bool use_drift_curv = boolean_parameter(
			index, "use_drift_curv:", false);
		const bool create_grid = boolean_parameter(index, "create_grid:", false);

		if (sections_only) {
			require("previous_networks:",
				"required when sections_simulation_only is true");

			// The domain is otherwise unused in sections-only mode. It becomes necessary
			// only for operations that explicitly rasterize or export the background grid.
			if (create_grid) {
				require("domain:", "required when create_grid is true");
			}
		}
		else {
			require("domain:", "required for full network simulation");
			require("sinks:", "at least one inlet is required for full network simulation");
			require("use_sinks_radius:", "required to define the expected sink property layout");
			require("springs:", "at least one outlet is required for full network simulation");
			require("use_springs_radius:", "required to define the expected spring property layout");
			require("allow_single_outlet_connection:",
				"required to select the ambiguous-link resolution strategy");
			require("use_user_connectivity_matrix:",
				"required to select file-based or generated connectivity");

			// Sampling is an exclusive workflow. Topography is required only by the
			// internal Poisson-sphere sampler. When a complete sampling cloud is supplied,
			// it remains an optional filter for added surface nodes.
			require("use_sampling_points:", "required to select the sampling workflow");
			if (has_parameter(index, "use_sampling_points:")) {
				const bool use_sampling_points = boolean_parameter(
					index, "use_sampling_points:", false);
				if (use_sampling_points) {
					require("sampling_points:", "required when use_sampling_points is true");
				}
				else {
					require("topo_surface:",
						"required when KarstNSim generates the sampling cloud");
					require("use_density_property:",
						"required when a sampling cloud must be generated");
					require("k_pts:",
						"required when a sampling cloud must be generated");
					if (has_parameter(index, "use_density_property:") &&
						!boolean_parameter(index, "use_density_property:", false)) {
						require("poisson_radius:",
							"required for constant-density Poisson sampling");
					}
				}
			}

			require("use_input_nghb_graph:",
				"required to select the graph-construction workflow");
			if (has_parameter(index, "use_input_nghb_graph:")) {
				const bool use_input_graph = boolean_parameter(
					index, "use_input_nghb_graph:", false);
				if (use_input_graph) {
					require("input_nghb_graph:",
						"required when use_input_nghb_graph is true");
				}
				else {
					require("nghb_count:",
						"required when the neighbor graph is generated internally");
					require("use_max_nghb_radius:",
						"required when the neighbor graph is generated internally");
					if (has_parameter(index, "use_max_nghb_radius:") &&
						boolean_parameter(index, "use_max_nghb_radius:", false)) {
						require("nghb_radius:",
							"required when use_max_nghb_radius is true");
					}
				}
			}

			for (const char* selector : {
				"use_waypoints:", "use_ghostrocks:", "add_inception_surfaces:",
				"use_karstification_potential:", "use_fracture_constraints:",
				"use_no_karst_spheres:", "use_deadend_points:",
				"use_cycle_amplification:", "use_noise:" }) {
				require(selector,
					"required to explicitly enable or disable this simulation component");
			}

			require("water_table_constraint_weight_vadose:",
				"required for full network simulation");
			require("water_table_constraint_weight_phreatic:",
				"required for full network simulation");
			require("fraction_karst_perm:", "required for path-cohesion updates");
			require("multiply_costs:", "required to define cost-term aggregation");
			require("gamma:", "required for the graph parameter set");
			require("vadose_cohesion:", "required to define cohesion behavior");
			require("vertical_distance_stretching_factor:",
				"required to define graph-distance anisotropy");

			if (has_parameter(index, "allow_single_outlet_connection:") &&
				boolean_parameter(index, "allow_single_outlet_connection:", false)) {
				require("gradient_constraint_weight:",
					"required for corrected closest-outlet selection");
				require("outlet_selection_cost_factor:",
					"required for corrected closest-outlet selection");
			}
		}

		if (use_previous_networks) {
			require("previous_networks:",
				"required when use_previous_networks is true");
			if (!sections_only) {
				require("fraction_old_karst_perm:",
					"required for polyphasic cost reduction");
			}
		}

		// Waypoints can provide section-conditioning radii in both workflows.
		// Their cost-reduction weight is used only during full graph simulation.
		if (use_waypoints) {
			require("waypoints:", "required when use_waypoints is true");
			require("use_waypoints_radius:",
				"required when use_waypoints is true");
			if (!sections_only) {
				require("waypoints_weight:",
					"required for waypoint cost reduction in full network simulation");
			}
		}

		// Ghost rocks can also overwrite conduit radii in sections-only mode.
		if (use_ghostrocks) {
			require("domain:", "required to rasterize ghost-rock corridors");
			require("use_karstification_potential:",
				"ghost rocks require the IKP workflow");
			require("alteration_lines:", "required when use_ghostrocks is true");
			require("ghostrock_max_vertical_size:",
				"required when use_ghostrocks is true");
			require("use_max_depth_constraint:",
				"required when use_ghostrocks is true");
			require("ghost_rock_weight:",
				"required when use_ghostrocks is true");
			require("ghostrock_width:",
				"required when use_ghostrocks is true");
			if (boolean_parameter(index, "use_max_depth_constraint:", false)) {
				require("max_depth_horizon:",
					"required when use_max_depth_constraint is true");
			}
		}

		// The IKP weight contributes to edge costs only during full network simulation.
		if (!sections_only && use_karstification_potential) {
			require("karstification_potential_weight:",
				"required when use_karstification_potential is true");
		}

		if (!sections_only) {
			if (boolean_parameter(index, "add_inception_surfaces:", false)) {
				require("refine_surface_sampling:",
					"required when add_inception_surfaces is true");
				require("inception_surfaces:",
					"required when add_inception_surfaces is true");
				require("inception_surface_constraint_weight:",
					"required when add_inception_surfaces is true");
				require("max_inception_surface_distance:",
					"required when add_inception_surfaces is true");
			}

			if (boolean_parameter(index, "use_fracture_constraints:", false)) {
				require("fracture_families_orientations:",
					"required when use_fracture_constraints is true");
				require("fracture_families_tolerance:",
					"required when use_fracture_constraints is true");
				require("fracture_constraint_weight:",
					"required when use_fracture_constraints is true");
			}

			if (boolean_parameter(index, "use_no_karst_spheres:", false)) {
				require("sphere_centers:",
					"required when use_no_karst_spheres is true");
			}

			if (boolean_parameter(index, "use_deadend_points:", false)) {
				require("nb_deadend_points:",
					"required when use_deadend_points is true");
				require("max_distance_of_deadend_pts:",
					"required when use_deadend_points is true");
			}

			if (boolean_parameter(index, "use_cycle_amplification:", false)) {
				require("max_distance_amplification:",
					"required when use_cycle_amplification is true");
				require("min_distance_amplification:",
					"required when use_cycle_amplification is true");
				require("nb_cycles:",
					"required when use_cycle_amplification is true");
			}

			if (boolean_parameter(index, "use_noise:", false)) {
				require("use_noise_on_all:", "required when use_noise is true");
				require("noise_frequency:", "required when use_noise is true");
				require("noise_octaves:", "required when use_noise is true");
				require("noise_weight:", "required when use_noise is true");
			}
		}

		if (simulate_sections) {
			for (const char* parameter : {
				"simulation_distribution:",
				"global_vario_range:", "global_range_of_neighborhood:",
				"global_vario_sill:", "global_vario_nugget:", "global_vario_model:",
				"interbranch_vario_range:", "interbranch_range_of_neighborhood:",
				"interbranch_vario_sill:", "interbranch_vario_nugget:",
				"interbranch_vario_model:",
				"intrabranch_vario_range:", "intrabranch_range_of_neighborhood:",
				"intrabranch_vario_sill:", "intrabranch_vario_nugget:",
				"intrabranch_vario_model:",
				"number_max_of_neighborhood_points:", "nb_points_interbranch:",
				"proportion_interbranch:", "use_drift_zwt:", "use_drift_curv:" }) {
				require(parameter, "required when simulate_sections is true");
			}

			if (has_parameter(index, "sinks:")) {
				require("use_sinks_radius:",
					"required to interpret optional sink conditioning radii");
			}
			if (has_parameter(index, "springs:")) {
				require("use_springs_radius:",
					"required to interpret optional spring conditioning radii");
			}
			if (use_sinks_radius) {
				require("sinks:",
					"required when use_sinks_radius is true");
			}
			if (use_springs_radius) {
				require("springs:",
					"required when use_springs_radius is true");
			}
			if (use_drift_zwt || use_drift_curv) {
				require("springs:",
					"required to define the outlet reference for external drift");
			}
		}

		if (!missing.empty()) {
			std::ostringstream message;
			message << "[instructions] Missing required parameter(s):";
			for (const std::string& item : missing) {
				message << "\n - " << item;
			}
			throw std::runtime_error(message.str());
		}
	}

	int checked_integer_property(
		float value,
		const std::string& object_name,
		const std::string& property_name,
		std::size_t object_index) {
		if (!std::isfinite(value) || std::fabs(value - std::round(value)) > 1e-5f) {
			throw std::runtime_error(
				"[" + object_name + "] Property '" + property_name + "' for item #" +
				std::to_string(object_index + 1) + " must be a finite integer."
			);
		}
		if (value < static_cast<float>(std::numeric_limits<int>::min()) ||
			value > static_cast<float>(std::numeric_limits<int>::max())) {
			throw std::runtime_error(
				"[" + object_name + "] Property '" + property_name + "' for item #" +
				std::to_string(object_index + 1) + " is outside the supported integer range."
			);
		}
		return static_cast<int>(std::round(value));
	}

	std::string strip_optional_quotes(const std::string& value) {
		if (value.size() >= 2 &&
			((value.front() == '"' && value.back() == '"') ||
			 (value.front() == '\'' && value.back() == '\''))) {
			return value.substr(1, value.size() - 2);
		}
		return value;
	}

	void validate_variogram(
		const std::string& label,
		float range,
		float neighborhood_range,
		float sill,
		float nugget,
		const std::string& model) {

		if (!std::isfinite(range) || range <= 0.0f) {
			throw std::runtime_error("[parameters][sections] '" + label + "_vario_range' must be > 0.");
		}
		if (!std::isfinite(neighborhood_range) || neighborhood_range <= 0.0f) {
			throw std::runtime_error("[parameters][sections] '" + label + "_range_of_neighborhood' must be > 0.");
		}
		if (!std::isfinite(sill) || sill < 0.0f) {
			throw std::runtime_error("[parameters][sections] '" + label + "_vario_sill' must be >= 0.");
		}
		if (!std::isfinite(nugget) || nugget < 0.0f || nugget > sill) {
			throw std::runtime_error("[parameters][sections] '" + label + "_vario_nugget' must satisfy 0 <= nugget <= sill.");
		}
		static const std::unordered_set<std::string> allowed_models = {
			"Exponential", "Spherical", "Gaussian", "Nugget"
		};
		if (allowed_models.count(model) == 0) {
			throw std::runtime_error(
				"[parameters][sections] Invalid variogram model '" + model +
				"' for " + label + ". Allowed models are Exponential, Spherical, Gaussian, and Nugget."
			);
		}
	}

	void validate_parameter_values(const KarstNSim::ParamsSource& params) {
		if (params.karstic_network_name.empty()) {
			throw std::runtime_error("[parameters] 'karstic_network_name' cannot be empty.");
		}
		if (params.save_repository.empty()) {
			throw std::runtime_error("[parameters] 'main_repository' cannot be empty.");
		}
		if (params.simulation_input_dir.empty()) {
			throw std::runtime_error("[parameters] 'simulation_input_dir' cannot be empty.");
		}
		if (params.selected_seed < 0) {
			throw std::runtime_error("[parameters] 'selected_seed' must be >= 0.");
		}
		if (params.number_of_iterations < 1) {
			throw std::runtime_error("[parameters] 'number_of_iterations' must be >= 1.");
		}

		const int nu = params.domain.get_nu();
		const int nv = params.domain.get_nv();
		const int nw = params.domain.get_nw();
		const bool has_valid_domain = nu > 0 && nv > 0 && nw > 0;
		const std::size_t cell_count = has_valid_domain
			? static_cast<std::size_t>(nu) * static_cast<std::size_t>(nv) *
				static_cast<std::size_t>(nw)
			: 0u;

		if (params.sections_simulation_only) {
			if (params.create_grid && !has_valid_domain) {
				throw std::runtime_error(
					"[parameters] 'create_grid: true' requires a valid domain grid, "
					"including in sections-only mode."
				);
			}
			if (!params.use_previous_networks || params.previous_networks.empty()) {
				throw std::runtime_error(
					"[parameters] 'sections_simulation_only: true' requires "
					"'use_previous_networks: true' and at least one previous network."
				);
			}
			if (!params.simulate_sections) {
				throw std::runtime_error(
					"[parameters] 'sections_simulation_only: true' requires 'simulate_sections: true'."
				);
			}
		}
		else {
			if (!has_valid_domain) {
				throw std::runtime_error(
					"[parameters] 'domain' must define strictly positive nu, nv, and nw dimensions."
				);
			}
			if (!params.use_sampling_points && params.topo_surface.is_empty()) {
				throw std::runtime_error(
					"[parameters] 'topo_surface' is required when KarstNSim generates "
					"the sampling cloud (use_sampling_points is false)."
				);
			}
			if (params.sinks.empty()) {
				throw std::runtime_error(
					"[parameters] At least one sink is required for full network simulation."
				);
			}
			if (params.springs.empty()) {
				throw std::runtime_error(
					"[parameters] At least one spring is required for full network simulation."
				);
			}

			if (!params.use_sampling_points) {
				if (params.k_pts < 1) {
					throw std::runtime_error(
						"[parameters] 'k_pts' must be >= 1 when sampling points are generated."
					);
				}
				if (params.use_density_property) {
					if (params.propdensity.size() != cell_count) {
						throw std::runtime_error(
							"[parameters] Density-property size mismatch: expected " +
							std::to_string(cell_count) + " values, but found " +
							std::to_string(params.propdensity.size()) + "."
						);
					}
					bool has_positive_density = false;
					for (float value : params.propdensity) {
						if (value < 0.0f) {
							continue;
						}
						if (!std::isfinite(value) || value <= 0.0f || value >= 1.0f) {
							throw std::runtime_error(
								"[parameters] Every valid density value must satisfy 0 < density < 1; "
								"negative values are reserved for no-data cells."
							);
						}
						has_positive_density = true;
					}
					if (!has_positive_density) {
						throw std::runtime_error(
							"[parameters] The density property contains no valid positive cell."
						);
					}
				}
				else if (!std::isfinite(params.poisson_radius) || params.poisson_radius <= 0.0f) {
					throw std::runtime_error(
						"[parameters] 'poisson_radius' must be > 0 for constant-density sampling."
					);
				}
			}
			else if (params.sampling_points.empty()) {
				throw std::runtime_error(
					"[parameters] 'sampling_points' is empty while use_sampling_points is true."
				);
			}

			if (params.use_input_nghb_graph) {
				if (params.input_nghb_graph.get_nb_nodes() <= 0 ||
					params.input_nghb_graph.get_nb_edges() <= 0) {
					throw std::runtime_error(
						"[parameters] 'input_nghb_graph' must contain at least one node and one edge."
					);
				}
			}
			else {
				if (params.nghb_count < 1) {
					throw std::runtime_error("[parameters] 'nghb_count' must be >= 1.");
				}
				if (params.use_max_nghb_radius &&
					(!std::isfinite(params.nghb_radius) || params.nghb_radius <= 0.0f)) {
					throw std::runtime_error(
						"[parameters] 'nghb_radius' must be > 0 when use_max_nghb_radius is true."
					);
				}
			}

			if (!std::isfinite(params.water_table_constraint_weight_vadose) ||
				params.water_table_constraint_weight_vadose < 0.0f ||
				!std::isfinite(params.water_table_constraint_weight_phreatic) ||
				params.water_table_constraint_weight_phreatic < 0.0f) {
				throw std::runtime_error(
					"[parameters] Water-table constraint weights must be finite and >= 0."
				);
			}

			int highest_water_table_index = 0;
			for (int index : params.propspringssurfindex) {
				if (index < 0) {
					throw std::runtime_error(
						"[parameters] Spring water-table indices must be >= 0."
					);
				}
				highest_water_table_index = std::max(highest_water_table_index, index);
			}
			if (highest_water_table_index > static_cast<int>(params.surf_wat_table.size())) {
				throw std::runtime_error(
					"[parameters] A spring references water-table index " +
					std::to_string(highest_water_table_index) + ", but only " +
					std::to_string(params.surf_wat_table.size()) +
					" water-table surface(s) were provided."
				);
			}
			for (const KarstNSim::Surface& surface : params.surf_wat_table) {
				if (surface.is_empty()) {
					throw std::runtime_error(
						"[parameters] At least one provided water-table surface is empty."
					);
				}
			}

			if (!std::isfinite(params.fraction_karst_perm) ||
				params.fraction_karst_perm <= 0.0f || params.fraction_karst_perm > 1.0f) {
				throw std::runtime_error(
					"[parameters] 'fraction_karst_perm' must satisfy 0 < value <= 1."
				);
			}
			if (!std::isfinite(params.gamma) || params.gamma <= 0.0f || params.gamma > 2.0f) {
				throw std::runtime_error("[parameters] 'gamma' must satisfy 0 < gamma <= 2.");
			}
			if (!std::isfinite(params.vertical_distance_stretching_factor) ||
				params.vertical_distance_stretching_factor < 1.0f) {
				throw std::runtime_error(
					"[parameters] 'vertical_distance_stretching_factor' must be >= 1."
				);
			}
			if (params.allow_single_outlet_connection) {
				if (!std::isfinite(params.gradient_constraint_weight) ||
					params.gradient_constraint_weight < 0.0f) {
					throw std::runtime_error(
						"[parameters] 'gradient_constraint_weight' must be >= 0."
					);
				}
				if (!std::isfinite(params.outlet_selection_cost_factor) ||
					params.outlet_selection_cost_factor < 1.0f) {
					throw std::runtime_error(
						"[parameters] 'outlet_selection_cost_factor' must be >= 1."
					);
				}
			}
		}

		if (params.use_previous_networks) {
			if (params.previous_networks.empty()) {
				throw std::runtime_error(
					"[parameters] 'previous_networks' is empty while use_previous_networks is true."
				);
			}
			if (!params.sections_simulation_only &&
				(!std::isfinite(params.fraction_old_karst_perm) ||
				 params.fraction_old_karst_perm <= 0.0f ||
				 params.fraction_old_karst_perm > 1.0f)) {
				throw std::runtime_error(
					"[parameters] 'fraction_old_karst_perm' must satisfy 0 < value <= 1."
				);
			}
		}

		// Section-conditioning properties may be supplied either in the original
		// physical domain or in a transformed domain, such as log10(radius).
		// Consequently, negative and zero values are valid; only non-finite values
		// are rejected here.
		auto validate_finite_section_properties = [](
			const std::vector<float>& values,
			const std::string& label) {

			for (std::size_t i = 0; i < values.size(); ++i) {
				if (!std::isfinite(values[i])) {
					throw std::runtime_error(
						"[parameters] Invalid " + label + " section-conditioning value " +
						std::to_string(values[i]) + " at item #" +
						std::to_string(i + 1) +
						". Section-conditioning properties must be finite, but may be "
						"negative or zero when expressed in a transformed domain."
					);
				}
			}
		};

		// Geometric radii and distances are expressed in model-coordinate units.
		// They must therefore be finite and strictly positive.
		auto validate_positive_geometric_radii = [](
			const std::vector<float>& radii,
			const std::string& label) {

			for (std::size_t i = 0; i < radii.size(); ++i) {
				if (!std::isfinite(radii[i]) || radii[i] <= 0.0f) {
					throw std::runtime_error(
						"[parameters] Invalid " + label + " radius " +
						std::to_string(radii[i]) + " at item #" +
						std::to_string(i + 1) +
						". Geometric radii must be finite and strictly greater than zero."
					);
				}
			}
		};

		if (params.use_sinks_radius) {
			if (params.propsinksradius.size() != params.sinks.size()) {
				throw std::runtime_error(
					"[parameters] use_sinks_radius is true, but the number of sink radii "
					"does not match the number of sinks."
				);
			}
			validate_finite_section_properties(
				params.propsinksradius,
				"sink"
			);
		}
		if (params.use_springs_radius) {
			if (params.propspringsradius.size() != params.springs.size()) {
				throw std::runtime_error(
					"[parameters] use_springs_radius is true, but the number of spring radii "
					"does not match the number of springs."
				);
			}
			validate_finite_section_properties(
				params.propspringsradius,
				"spring"
			);
		}

		if (params.use_waypoints) {
			if (params.waypoints.empty()) {
				throw std::runtime_error(
					"[parameters] 'waypoints' is empty while use_waypoints is true."
				);
			}
			if (params.waypoints_impact_radius.size() != params.waypoints.size()) {
				throw std::runtime_error(
					"[parameters] Waypoints require one impact radius per point."
				);
			}
			if (!std::isfinite(params.waypoints_weight) ||
				params.waypoints_weight < 0.0f || params.waypoints_weight >= 1.0f) {
				throw std::runtime_error(
					"[parameters] 'waypoints_weight' must satisfy 0 <= value < 1."
				);
			}
			validate_positive_geometric_radii(
				params.waypoints_impact_radius,
				"waypoint impact"
			);
			if (params.use_waypoints_radius) {
				if (params.waypoints_radius.size() != params.waypoints.size()) {
					throw std::runtime_error(
						"[parameters] use_waypoints_radius is true, but the number of radii "
						"does not match the number of waypoints."
					);
				}
				validate_finite_section_properties(
					params.waypoints_radius,
					"waypoint"
				);
			}
		}

		if (params.use_ghostrocks) {
			if (!params.use_karstification_potential) {
				throw std::runtime_error(
					"[parameters] 'use_ghostrocks: true' requires "
					"'use_karstification_potential: true'."
				);
			}
			if (!has_valid_domain) {
				throw std::runtime_error(
					"[parameters] Ghost-rock processing requires a valid domain grid."
				);
			}
			if (params.propikp.size() != cell_count) {
				throw std::runtime_error(
					"[parameters] Ghost-rock processing requires one IKP value per domain cell: expected " +
					std::to_string(cell_count) + ", found " +
					std::to_string(params.propikp.size()) + "."
				);
			}
			if (!std::isfinite(params.ghostrock_max_vertical_size) ||
				params.ghostrock_max_vertical_size <= 0.0f ||
				!std::isfinite(params.ghostrock_width) || params.ghostrock_width <= 0.0f ||
				!std::isfinite(params.ghost_rock_weight) || params.ghost_rock_weight < 0.0f) {
				throw std::runtime_error(
					"[parameters] Ghost-rock dimensions must be > 0 and ghost_rock_weight must be >= 0."
				);
			}
			if (params.use_max_depth_constraint && params.max_depth_horizon.is_empty()) {
				throw std::runtime_error(
					"[parameters] 'max_depth_horizon' is empty while use_max_depth_constraint is true."
				);
			}
		}

		if (params.use_karstification_potential &&
			(!params.sections_simulation_only || params.use_ghostrocks)) {
			if (params.propikp.size() != cell_count) {
				throw std::runtime_error(
					"[parameters] IKP-property size mismatch: expected " +
					std::to_string(cell_count) + " values, but found " +
					std::to_string(params.propikp.size()) + "."
				);
			}
			for (float value : params.propikp) {
				if (value < 0.0f) {
					continue;
				}
				if (!std::isfinite(value) || value > 1.0f) {
					throw std::runtime_error(
						"[parameters] Every valid IKP value must satisfy 0 <= IKP <= 1."
					);
				}
			}
			if (!std::isfinite(params.karstification_potential_weight) ||
				params.karstification_potential_weight < 0.0f) {
				throw std::runtime_error(
					"[parameters] 'karstification_potential_weight' must be >= 0."
				);
			}
		}

		if (!params.sections_simulation_only) {
			if (params.add_inception_surfaces) {
				if (params.refine_surface_sampling < -1) {
					throw std::runtime_error(
						"[parameters] 'refine_surface_sampling' must be >= -1."
					);
				}
				if (params.inception_surfaces.empty()) {
					throw std::runtime_error(
						"[parameters] At least one inception surface is required when "
						"add_inception_surfaces is true."
					);
				}
				if (!std::isfinite(params.inception_surface_constraint_weight) ||
					params.inception_surface_constraint_weight < 0.0f ||
					!std::isfinite(params.max_inception_surface_distance) ||
					params.max_inception_surface_distance <= 0.0f) {
					throw std::runtime_error(
						"[parameters] Inception-surface weight must be >= 0 and maximum distance must be > 0."
					);
				}
			}

			if (params.use_fracture_constraints) {
				if (params.fracture_families_orientations.empty() ||
					params.fracture_families_orientations.size() !=
					params.fracture_families_tolerance.size()) {
					throw std::runtime_error(
						"[parameters] Fracture orientation and tolerance lists must be non-empty "
						"and have equal sizes."
					);
				}
				for (float orientation : params.fracture_families_orientations) {
					if (!std::isfinite(orientation) || orientation < 0.0f || orientation > 180.0f) {
						throw std::runtime_error(
							"[parameters] Fracture orientations must lie in [0, 180] degrees."
						);
					}
				}
				for (float tolerance : params.fracture_families_tolerance) {
					if (!std::isfinite(tolerance) || tolerance < 0.0f || tolerance > 90.0f) {
						throw std::runtime_error(
							"[parameters] Fracture tolerances must lie in [0, 90] degrees."
						);
					}
				}
				if (!std::isfinite(params.fracture_constraint_weight) ||
					params.fracture_constraint_weight < 0.0f) {
					throw std::runtime_error(
						"[parameters] 'fracture_constraint_weight' must be >= 0."
					);
				}
			}

			if (params.use_no_karst_spheres) {
				if (params.sphere_centers.empty() ||
					params.sphere_radius.size() != params.sphere_centers.size()) {
					throw std::runtime_error(
						"[parameters] No-karst spheres require one radius per sphere center."
					);
				}
				validate_positive_geometric_radii(
					params.sphere_radius,
					"no-karst sphere"
				);
			}

			if (params.use_deadend_points) {
				if (params.nb_deadend_points < 0 ||
					!std::isfinite(params.max_distance_of_deadend_pts) ||
					params.max_distance_of_deadend_pts <= 0.0f) {
					throw std::runtime_error(
						"[parameters] Dead-end count must be >= 0 and maximum distance must be > 0."
					);
				}
			}

			if (params.use_amplification) {
				if (params.nb_cycles < 0 ||
					!std::isfinite(params.min_distance_amplification) ||
					!std::isfinite(params.max_distance_amplification) ||
					params.min_distance_amplification < 0.0f ||
					params.max_distance_amplification <= 0.0f ||
					params.min_distance_amplification > params.max_distance_amplification) {
					throw std::runtime_error(
						"[parameters] Cycle amplification requires nb_cycles >= 0 and "
						"0 <= min_distance_amplification <= max_distance_amplification with max > 0."
					);
				}
			}

			if (params.use_noise_on_all && !params.use_noise) {
				throw std::runtime_error(
					"[parameters] 'use_noise_on_all: true' requires 'use_noise: true'."
				);
			}
			if (params.use_noise) {
				if (params.noise_frequency < 1 || params.noise_octaves < 1 ||
					!std::isfinite(params.noise_weight) || params.noise_weight < 0.0f) {
					throw std::runtime_error(
						"[parameters] Noise frequency and octaves must be >= 1, and noise_weight must be >= 0."
					);
				}
			}
		}

		if (params.simulate_sections) {
			if (params.geostat_params.simulation_distribution.size() < 2u) {
				throw std::runtime_error(
					"[parameters][sections] 'simulation_distribution' must contain at least "
					"two numerical values for normal-score transformation and back-transformation."
				);
			}

			for (float value : params.geostat_params.simulation_distribution) {
				if (!std::isfinite(value)) {
					throw std::runtime_error(
						"[parameters][sections] 'simulation_distribution' contains a non-finite value."
					);
				}
			}

			validate_variogram(
				"global", params.geostat_params.global_vario_range,
				params.geostat_params.global_range_of_neighborhood,
				params.geostat_params.global_vario_sill,
				params.geostat_params.global_vario_nugget,
				params.geostat_params.global_vario_model);
			validate_variogram(
				"interbranch", params.geostat_params.interbranch_vario_range,
				params.geostat_params.interbranch_range_of_neighborhood,
				params.geostat_params.interbranch_vario_sill,
				params.geostat_params.interbranch_vario_nugget,
				params.geostat_params.interbranch_vario_model);
			validate_variogram(
				"intrabranch", params.geostat_params.intrabranch_vario_range,
				params.geostat_params.intrabranch_range_of_neighborhood,
				params.geostat_params.intrabranch_vario_sill,
				params.geostat_params.intrabranch_vario_nugget,
				params.geostat_params.intrabranch_vario_model);

			if (params.geostat_params.number_max_of_neighborhood_points < 1) {
				throw std::runtime_error(
					"[parameters][sections] 'number_max_of_neighborhood_points' must be >= 1."
				);
			}
			if (params.geostat_params.nb_points_interbranch < 0) {
				throw std::runtime_error(
					"[parameters][sections] 'nb_points_interbranch' must be >= 0."
				);
			}
			if (!std::isfinite(params.geostat_params.proportion_interbranch) ||
				params.geostat_params.proportion_interbranch < 0.0f ||
				params.geostat_params.proportion_interbranch > 1.0f) {
				throw std::runtime_error(
					"[parameters][sections] 'proportion_interbranch' must lie in [0, 1]."
				);
			}

			if ((params.geostat_params.use_drift_zwt ||
				 params.geostat_params.use_drift_curv) && params.springs.empty()) {
				throw std::runtime_error(
					"[parameters][sections] External drift requires at least one spring "
					"to define the outlet reference."
				);
			}
			if (params.geostat_params.use_drift_zwt || params.geostat_params.use_drift_curv) {
				const std::size_t observation_count =
					(params.use_sinks_radius ? params.propsinksradius.size() : 0u) +
					(params.use_springs_radius ? params.propspringsradius.size() : 0u) +
					(params.use_waypoints_radius ? params.waypoints_radius.size() : 0u);
				if (observation_count < 2u) {
					throw std::runtime_error(
						"[parameters][sections] External drift requires at least two radius observations "
						"from sinks, springs, and/or waypoints."
					);
				}
			}
		}
	}
}

ParseInputs::ParseInputs() {
	// Constructor implementation, if needed
}


bool ParseInputs::parse_boolean(const std::string& input) {
	return parse_boolean_token_strict(input, "boolean parameter", 0);
}

std::vector<float> ParseInputs::read_distrib(const std::string& directory, const std::string& filename) {
	std::vector<float> values;

	// --- Defensive guard: allow empty/none filenames gracefully ---
	// Treat empty or sentinel values as "no distribution provided".
	if (filename.empty() ||
		filename == "-" ||
		filename == "none" || filename == "None" || filename == "NONE" ||
		filename == "null" || filename == "Null" || filename == "NULL") {
		// Not an error: just return an empty vector. The caller can check emptiness if needed.
		return values;
	}

	const std::string full_path = directory + "/" + filename;

	std::ifstream file(full_path);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << full_path << std::endl;
		std::cerr << "Reason: " << std::strerror(errno) << std::endl;
		throw std::runtime_error("File not found: " + full_path);
	}

	std::string line;
	std::size_t line_number = 0;
	while (std::getline(file, line)) {
		++line_number;
		std::stringstream ss(line);
		while (true) {
			ss >> std::ws;
			if (ss.eof()) {
				break;
			}
			float value = 0.0f;
			if (!(ss >> value)) {
				throw std::runtime_error(
					"[simulation_distribution] Invalid floating-point value in " +
					full_path + " at line " + std::to_string(line_number) + "."
				);
			}
			values.push_back(value);
		}
	}

	if (values.empty()) {
		throw std::runtime_error(
			"[simulation_distribution] No numerical value was read from " + full_path + "."
		);
	}

	// Normal-score transforms and quantile lookup require an ordered empirical
	// distribution, independently of the order used in the input file.
	std::sort(values.begin(), values.end());

	return values;
}


KarstNSim::ParamsSource ParseInputs::parse(const std::string& filename) {

	KarstNSim::ParamsSource params;

	// parameters not included formally with a tag in the prompt file, but defined in another file (properties) :

	// extracted from domain box :
	std::vector<float> propdensity;
	std::vector<float> propikp;
	// extracted from springs pointset :
	std::vector<int> propspringsindex;
	std::vector<float> propspringsradius; 
	std::vector<int> propspringssurfindex; 
	// extracted from sinks pointset :
	std::vector<int> propsinksindex;
	std::vector<int> propsinksorder;
	std::vector<float> propsinksradius; 
	// extracted from no-karst spheres pointset :
	std::vector<float> sphere_radius;
	// extracted from waypoints pointset :
	std::vector<float> waypoints_radius;
	std::vector<float> waypoints_impact_radius; 

	std::ifstream inputFile(filename);

	if (!inputFile) {
		std::cerr << "Error opening file: " << filename << std::endl;
		std::cerr << "Error: " << std::strerror(errno) << std::endl;
		throw std::runtime_error("Failed to open input file");
	}
	// Default: input files are located in the same directory as the main instruction file.
	params.simulation_input_dir = std::filesystem::absolute(filename).parent_path().string();

	// First pass: validate instruction syntax and determine which parameters are mandatory
	// before any external ASCII object is loaded. This prevents missing switches or file
	// parameters from degenerating into unrelated downstream crashes.
	const InstructionIndex instruction_index = scan_instruction_file(inputFile);
	validate_instruction_value_syntax(instruction_index);
	validate_required_parameter_presence(instruction_index);

	// Resolve path and property-layout switches before the second pass so that parameter
	// declaration order does not affect file loading.
	params.save_repository = first_parameter_token(instruction_index, "main_repository:");
	if (has_parameter(instruction_index, "simulation_input_dir:")) {
		params.simulation_input_dir = first_parameter_token(instruction_index, "simulation_input_dir:");
	}
	params.sections_simulation_only = boolean_parameter(
		instruction_index, "sections_simulation_only:", false);
	params.simulate_sections = boolean_parameter(
		instruction_index, "simulate_sections:", false);
	params.geostat_params.is_used = params.simulate_sections;
	params.geostat_params.use_drift_zwt = boolean_parameter(
		instruction_index, "use_drift_zwt:", false);
	params.geostat_params.use_drift_curv = boolean_parameter(
		instruction_index, "use_drift_curv:", false);

	params.use_previous_networks = boolean_parameter(
		instruction_index, "use_previous_networks:", false);
	params.use_sampling_points = boolean_parameter(
		instruction_index, "use_sampling_points:", false);
	params.use_density_property = boolean_parameter(
		instruction_index, "use_density_property:", false);
	params.use_input_nghb_graph = boolean_parameter(
		instruction_index, "use_input_nghb_graph:", false);

	params.use_sinks_radius = boolean_parameter(
		instruction_index, "use_sinks_radius:", false);
	params.use_springs_radius = boolean_parameter(
		instruction_index, "use_springs_radius:", false);
	params.use_waypoints = boolean_parameter(
		instruction_index, "use_waypoints:", false);
	params.use_waypoints_radius = boolean_parameter(
		instruction_index, "use_waypoints_radius:", false);

	params.use_ghostrocks = boolean_parameter(
		instruction_index, "use_ghostrocks:", false);
	params.use_max_depth_constraint = boolean_parameter(
		instruction_index, "use_max_depth_constraint:", false);
	params.add_inception_surfaces = boolean_parameter(
		instruction_index, "add_inception_surfaces:", false);
	params.use_karstification_potential = boolean_parameter(
		instruction_index, "use_karstification_potential:", false);
	params.use_fracture_constraints = boolean_parameter(
		instruction_index, "use_fracture_constraints:", false);
	params.use_no_karst_spheres = boolean_parameter(
		instruction_index, "use_no_karst_spheres:", false);
	params.use_deadend_points = boolean_parameter(
		instruction_index, "use_deadend_points:", false);
	params.use_amplification = boolean_parameter(
		instruction_index, "use_cycle_amplification:", false);
	params.use_noise = boolean_parameter(
		instruction_index, "use_noise:", false);
	params.create_grid = boolean_parameter(
		instruction_index, "create_grid:", false);

	std::string line;
	while (std::getline(inputFile, line)) {

		line = trim_copy(line);
		if (line.empty() || is_instruction_comment(line)) {
			continue;
		}

		// File name and repository

		std::istringstream iss(line);
		std::string paramType;
		iss >> paramType;
		if (paramType == "karstic_network_name:") {
			iss >> params.karstic_network_name;
		}
		else if (paramType == "main_repository:") {
			iss >> params.save_repository;
		}
		else if (paramType == "simulation_input_dir:") {
			// allow changing the input directory from the default (same as main input file)
			iss >> params.simulation_input_dir;
		}

		// General parameters

		else if (paramType == "domain:") {
			const bool domain_is_used = !params.sections_simulation_only ||
				params.use_ghostrocks || params.create_grid;
			if (!domain_is_used) {
				continue;
			}
			std::string box_path;
			iss >> box_path;
			KarstNSim::Box box;
			std::vector<std::vector<float>> prop_box;
			KarstNSim::load_box(box_path, params.save_repository, box, prop_box);
			params.domain = box;

			const std::size_t cell_count = static_cast<std::size_t>(box.get_nu()) *
				static_cast<std::size_t>(box.get_nv()) * static_cast<std::size_t>(box.get_nw());
			if (prop_box.size() != cell_count) {
				throw std::runtime_error(
					"[domain] Property-row count mismatch: expected " +
					std::to_string(cell_count) + " rows, but found " +
					std::to_string(prop_box.size()) + "."
				);
			}
			for (std::size_t i = 0; i < cell_count; ++i) {
				if (!prop_box[i].empty()) {
					params.propdensity.push_back(prop_box[i][0]);
				}
				if (prop_box[i].size() > 1) {
					params.propikp.push_back(prop_box[i][1]);
				}
			}
		}
		else if (paramType == "selected_seed:") {
			iss >> params.selected_seed;
		}
		else if (paramType == "number_of_iterations:") {
			iss >> params.number_of_iterations;
		}
		else if (paramType == "vary_seed:") {
			std::string flag;
			iss >> flag;
			params.vary_seed = parse_boolean(flag);
		}
		else if (paramType == "topo_surface:") {
			if (params.sections_simulation_only) {
				continue;
			}
			std::string surf_path;
			iss >> surf_path;
			std::vector<std::vector<float>> prop;
			KarstNSim::Surface surf;
			KarstNSim::load_surface(surf_path, params.save_repository, surf, prop);
			params.topo_surface = surf;
		}

		// Sampling

		else if (paramType == "use_sampling_points:") {
			std::string flag;
			iss >> flag;
			params.use_sampling_points = parse_boolean(flag);
		}
		else if (paramType == "sampling_points:") {
			if (params.sections_simulation_only || !params.use_sampling_points) {
				continue;
			}
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<float>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repository,ptset, prop);
			params.sampling_points = ptset;
		}
		else if (paramType == "poisson_radius:") {
			iss >> params.poisson_radius;
		}
		else if (paramType == "use_density_property:") {
			std::string flag;
			iss >> flag;
			params.use_density_property = parse_boolean(flag);
		}
		else if (paramType == "k_pts:") {
			iss >> params.k_pts;
		}

		// Nearest neighbor graph

		else if (paramType == "nghb_count:") {
			iss >> params.nghb_count;
		}
		else if (paramType == "use_max_nghb_radius:") {
			std::string flag;
			iss >> flag;
			params.use_max_nghb_radius = parse_boolean(flag);
		}
		else if (paramType == "nghb_radius:") {
			iss >> params.nghb_radius;
		}

		// Sinks springs and waypoints

		else if (paramType == "use_springs_radius:") {
			std::string flag;
			iss >> flag;
			params.use_springs_radius = parse_boolean(flag);
		}
		else if (paramType == "use_sinks_radius:") {
			std::string flag;
			iss >> flag;
			params.use_sinks_radius = parse_boolean(flag);
		}
		else if (paramType == "use_waypoints_radius:") {
			std::string flag;
			iss >> flag;
			params.use_waypoints_radius = parse_boolean(flag);
		}
		else if (paramType == "sinks:") {
			if (params.sections_simulation_only && !params.use_sinks_radius) {
				continue;
			}
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<float>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repository, ptset, prop);
			params.sinks = ptset;

			if (prop.size() != ptset.size()) {
				throw std::runtime_error(
					"[sinks] Property-row count does not match the number of sink points."
				);
			}

			for (std::size_t i = 0; i < ptset.size(); ++i) {
				if (params.sections_simulation_only) {
					// In sections-only mode, connectivity metadata is optional. A compact
					// one-property point set may store the radius directly in the first
					// property column. The full index/order/radius layout remains accepted.
					if (params.use_sinks_radius) {
						if (prop[i].size() >= 3u) {
							params.propsinksradius.push_back(prop[i][2]);
						}
						else if (!prop[i].empty()) {
							params.propsinksradius.push_back(prop[i][0]);
						}
						else {
							throw std::runtime_error(
								"[sinks] Sink #" + std::to_string(i + 1) +
								" is missing its radius property in sections-only mode."
							);
						}
					}
					continue;
				}

				const std::size_t required_columns = params.use_sinks_radius ? 3u : 2u;
				if (prop[i].size() < required_columns) {
					throw std::runtime_error(
						"[sinks] Sink #" + std::to_string(i + 1) + " requires at least " +
						std::to_string(required_columns) + " properties (index, order" +
						(params.use_sinks_radius ? ", radius" : "") + ")."
					);
				}
				params.propsinksindex.push_back(
					checked_integer_property(prop[i][0], "sinks", "index", i));
				params.propsinksorder.push_back(
					checked_integer_property(prop[i][1], "sinks", "order", i));
				if (params.use_sinks_radius) {
					params.propsinksradius.push_back(prop[i][2]);
				}
			}
		}
		else if (paramType == "springs:") {
			const bool springs_are_used = !params.sections_simulation_only ||
				params.use_springs_radius || params.geostat_params.use_drift_zwt ||
				params.geostat_params.use_drift_curv;
			if (!springs_are_used) {
				continue;
			}
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<float>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repository, ptset, prop);
			params.springs = ptset;

			if (prop.size() != ptset.size()) {
				throw std::runtime_error(
					"[springs] Property-row count does not match the number of spring points."
				);
			}

			for (std::size_t i = 0; i < ptset.size(); ++i) {
				if (params.sections_simulation_only) {
					// Sections-only mode accepts either a compact radius-only property
					// layout or the full index/surfindex/radius layout.
					if (params.use_springs_radius) {
						if (prop[i].size() >= 3u) {
							params.propspringsradius.push_back(prop[i][2]);
						}
						else if (!prop[i].empty()) {
							params.propspringsradius.push_back(prop[i][0]);
						}
						else {
							throw std::runtime_error(
								"[springs] Spring #" + std::to_string(i + 1) +
								" is missing its radius property in sections-only mode."
							);
						}
					}
					continue;
				}

				const std::size_t required_columns = params.use_springs_radius ? 3u : 2u;
				if (prop[i].size() < required_columns) {
					throw std::runtime_error(
						"[springs] Spring #" + std::to_string(i + 1) + " requires at least " +
						std::to_string(required_columns) + " properties (index, surfindex" +
						(params.use_springs_radius ? ", radius" : "") + ")."
					);
				}
				params.propspringsindex.push_back(
					checked_integer_property(prop[i][0], "springs", "index", i));
				params.propspringssurfindex.push_back(
					checked_integer_property(prop[i][1], "springs", "surfindex", i));
				if (params.use_springs_radius) {
					params.propspringsradius.push_back(prop[i][2]);
				}
			}
		}
		else if (paramType == "use_waypoints:") {
		std::string flag;
		iss >> flag;
		params.use_waypoints = parse_boolean(flag);
		}
		else if (paramType == "waypoints:") {
			if (!params.use_waypoints) {
				continue;
			}
			std::string pts_path;
			iss >> pts_path;
			std::vector<std::vector<float>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repository, ptset, prop);
			params.waypoints = ptset;
			const std::size_t required_columns = params.use_waypoints_radius ? 2u : 1u;
			if (prop.size() != ptset.size()) {
				throw std::runtime_error("[waypoints] Property-row count does not match the number of waypoint points.");
			}
			for (std::size_t i = 0; i < ptset.size(); ++i) {
				if (prop[i].size() < required_columns) {
					throw std::runtime_error(
						"[waypoints] Waypoint #" + std::to_string(i + 1) +
						" requires an impact_radius property" +
						(params.use_waypoints_radius ? " and a radius property." : ".")
					);
				}
				params.waypoints_impact_radius.push_back(prop[i][0]);
				if (params.use_waypoints_radius) {
					params.waypoints_radius.push_back(prop[i][1]);
				}
			}
		}
		else if (paramType == "waypoints_weight:") {
			iss >> params.waypoints_weight;
		}
		else if (paramType == "allow_single_outlet_connection:") {
		std::string flag;
		iss >> flag;
		params.allow_single_outlet_connection = parse_boolean(flag);
		}
		else if (paramType == "use_user_connectivity_matrix:") {
			std::string flag;
			iss >> flag;
			params.use_user_connectivity_matrix = parse_boolean(flag);
		}

		// ghost-rocks

		else if (paramType == "use_ghostrocks:") {
		std::string flag;
		iss >> flag;
		params.use_ghostrocks = parse_boolean(flag);
		}
		else if (paramType == "alteration_lines:") {
			if (!params.use_ghostrocks) {
				continue;
			}
		std::string line_path;
		std::vector<std::vector<std::vector<float>>> prop;
		KarstNSim::Line line;
		iss >> line_path;
		KarstNSim::load_line(line_path, params.save_repository, line, prop);
		params.alteration_lines = line;
		}
		else if (paramType == "interpolate_lines:") {
		std::string flag;
		iss >> flag;
		params.interpolate_lines = parse_boolean(flag);
		}
		else if (paramType == "ghostrock_max_vertical_size:") {
		iss >> params.ghostrock_max_vertical_size;
		}
		else if (paramType == "use_max_depth_constraint:") {
		std::string flag;
		iss >> flag;
		params.use_max_depth_constraint = parse_boolean(flag);
		}
		else if (paramType == "ghost_rock_weight:") {
		iss >> params.ghost_rock_weight;
		}
		else if (paramType == "max_depth_horizon:") {
			if (!params.use_ghostrocks || !params.use_max_depth_constraint) {
				continue;
			}
		std::string surf_path;
		std::vector<std::vector<float>> prop;
		KarstNSim::Surface surf;
		iss >> surf_path;
		KarstNSim::load_surface(surf_path, params.save_repository, surf, prop);
		params.max_depth_horizon = surf;
		}
		else if (paramType == "ghostrock_width:") {
		iss >> params.ghostrock_width;
		}

		// Inception surfaces

		else if (paramType == "refine_surface_sampling:") {
			iss >> params.refine_surface_sampling;
		}
		else if (paramType == "add_inception_surfaces:") {
		std::string flag;
		iss >> flag;
		params.add_inception_surfaces = parse_boolean(flag);
		}
		else if (paramType == "inception_surfaces:") {
			if (params.sections_simulation_only || !params.add_inception_surfaces) {
				continue;
			}
			std::string surf_path;
			std::vector<std::vector<float>> prop;
			std::vector<KarstNSim::Surface> surf;
			while (iss >> surf_path) {
				KarstNSim::Surface surfi;
				KarstNSim::load_surface(surf_path, params.save_repository, surfi, prop);
				surf.push_back(surfi);
			}

			params.inception_surfaces = surf;
		}
		else if (paramType == "inception_surface_constraint_weight:") {
			iss >> params.inception_surface_constraint_weight;
		}
		else if (paramType == "max_inception_surface_distance:") {
		iss >> params.max_inception_surface_distance;
		}

		// Karstification potential

		else if (paramType == "use_karstification_potential:") {
		std::string flag;
		iss >> flag;
		params.use_karstification_potential = parse_boolean(flag);
		}
		else if (paramType == "karstification_potential_weight:") {
		iss >> params.karstification_potential_weight;
		}

		// Fracture constraints

		else if (paramType == "use_fracture_constraints:") {
		std::string flag;
		iss >> flag;
		params.use_fracture_constraints = parse_boolean(flag);
		}
		else if (paramType == "fracture_families_orientations:") {
		float value;
			while (iss >> value) {
				params.fracture_families_orientations.emplace_back(value);
			}
		}
		else if (paramType == "fracture_families_tolerance:") {
		float value;
			while (iss >> value) {
				params.fracture_families_tolerance.emplace_back(value);
			}
		}
		else if (paramType == "fracture_constraint_weight:") {
		iss >> params.fracture_constraint_weight;
		}

		// Previous networks (polyphasic karstification or simulation of dimensions exclusively)

		else if (paramType == "use_previous_networks:") {
		std::string flag;
		iss >> flag;
		params.use_previous_networks = parse_boolean(flag);
		}
		else if (paramType == "previous_networks:") {
			if (!params.use_previous_networks) {
				continue;
			}
			std::string line_path;
			std::vector<std::vector<std::vector<float>>> prop;
			std::vector<KarstNSim::Line> lines;
			while (iss >> line_path) {
				KarstNSim::Line line;
				KarstNSim::load_line(line_path, params.save_repository, line, prop);
				lines.push_back(line);
			}
			params.previous_networks = lines;
		}
		else if (paramType == "fraction_old_karst_perm:") {
		iss >> params.fraction_old_karst_perm;
		}
		else if (paramType == "sections_simulation_only:") {
		std::string flag;
		iss >> flag;
		params.sections_simulation_only = parse_boolean(flag);
		}

		// Previous nghb graph (if needed)

		else if (paramType == "use_input_nghb_graph:") {
			std::string flag;
			iss >> flag;
			params.use_input_nghb_graph = parse_boolean(flag);
		}
		else if (paramType == "input_nghb_graph:") {
			if (params.sections_simulation_only || !params.use_input_nghb_graph) {
				continue;
			}
			std::string graph_path;
			KarstNSim::InputGraph graph;
			iss >> graph_path;
			graph = KarstNSim::translate_input_graph(graph_path, params.save_repository);
			params.input_nghb_graph = graph;
		}

		// No-karst spheres

		else if (paramType == "use_no_karst_spheres:") {
		std::string flag;
		iss >> flag;
		params.use_no_karst_spheres = parse_boolean(flag);
		}
		else if (paramType == "sphere_centers:") {
			if (params.sections_simulation_only || !params.use_no_karst_spheres) {
				continue;
			}
			std::string pts_path;
			iss >> pts_path;

			std::vector<std::vector<float>> prop;
			std::vector<Vector3> ptset;
			KarstNSim::load_pointset(pts_path, params.save_repository, ptset, prop);
			params.sphere_centers = ptset;

			if (prop.size() != ptset.size()) {
				throw std::runtime_error("[sphere_centers] Property-row count does not match the number of sphere centers.");
			}
			sphere_radius.reserve(ptset.size());
			for (std::size_t i = 0; i < ptset.size(); ++i) {
				if (prop[i].empty()) {
					throw std::runtime_error(
						"[sphere_centers] Sphere #" + std::to_string(i + 1) +
						" is missing its mandatory radius property."
					);
				}
				sphere_radius.push_back(prop[i][0]);
			}
			params.sphere_radius = sphere_radius;
		}

		// Water tables

		else if (paramType == "water_table_constraint_weight_vadose:") {
		iss >> params.water_table_constraint_weight_vadose;
		}
		else if (paramType == "water_table_constraint_weight_phreatic:") {
		iss >> params.water_table_constraint_weight_phreatic;
		}
		else if (paramType == "surf_wat_table:") {
			if (params.sections_simulation_only) {
				continue;
			}
		std::string surf_path;
		std::vector<std::vector<float>> prop;
		std::vector<KarstNSim::Surface> surf;
		while (iss >> surf_path) {
			KarstNSim::Surface surfi;
			KarstNSim::load_surface(surf_path, params.save_repository, surfi, prop);
			surf.push_back(surfi);
		}
		params.surf_wat_table = surf;
		}

		// Deadend points amplification :

		else if (paramType == "use_deadend_points:") {
		std::string flag;
		iss >> flag;
		params.use_deadend_points = parse_boolean(flag);
		}
		else if (paramType == "nb_deadend_points:") {
		iss >> params.nb_deadend_points;
		}
		else if (paramType == "max_distance_of_deadend_pts:") {
		iss >> params.max_distance_of_deadend_pts;
		}

		// Cycle amplification :

		else if (paramType == "use_cycle_amplification:") {
		std::string flag;
		iss >> flag;
		params.use_amplification = parse_boolean(flag);
		}
		else if (paramType == "max_distance_amplification:") {
		iss >> params.max_distance_amplification;
		}
		else if (paramType == "min_distance_amplification:") {
		iss >> params.min_distance_amplification;
		}
		else if (paramType == "nb_cycles:") {
		iss >> params.nb_cycles;
		}

		//Noise Amplification parameters :

		else if (paramType == "use_noise:") {
		std::string flag;
		iss >> flag;
		params.use_noise = parse_boolean(flag);
		}
		else if (paramType == "use_noise_on_all:") {
		std::string flag;
		iss >> flag;
		params.use_noise_on_all = parse_boolean(flag);
		}
		else if (paramType == "noise_frequency:") {
		iss >> params.noise_frequency;
		}
		else if (paramType == "noise_octaves:") {
		iss >> params.noise_octaves;
		}
		else if (paramType == "noise_weight:") {
		iss >> params.noise_weight;
		}

		// Create sections :

		else if (paramType == "simulate_sections:") {
		std::string flag;
		iss >> flag;
		params.simulate_sections = parse_boolean(flag);
		if (params.simulate_sections) {
			params.geostat_params.is_used = true;
		}
		}
		else if (paramType == "simulation_distribution:") {
			if (!params.simulate_sections) {
				continue;
			}
		std::string path;
		iss >> path;
		std::vector<float> sim_distrib_discrete = read_distrib(params.save_repository, path);
		params.geostat_params.simulation_distribution = sim_distrib_discrete;
		}
		else if (paramType == "global_vario_range:") {
		iss >> params.geostat_params.global_vario_range;
		}
		else if (paramType == "global_range_of_neighborhood:") {
		iss >> params.geostat_params.global_range_of_neighborhood;
		}
		else if (paramType == "global_vario_sill:") {
		iss >> params.geostat_params.global_vario_sill;
		}
		else if (paramType == "global_vario_nugget:") {
		iss >> params.geostat_params.global_vario_nugget;
		}
		else if (paramType == "global_vario_model:") {
			std::string model;
			iss >> model;
			params.geostat_params.global_vario_model = strip_optional_quotes(model);
		}
		else if (paramType == "interbranch_vario_range:") {
		iss >> params.geostat_params.interbranch_vario_range;
		}
		else if (paramType == "interbranch_range_of_neighborhood:") {
		iss >> params.geostat_params.interbranch_range_of_neighborhood;
		}
		else if (paramType == "interbranch_vario_sill:") {
		iss >> params.geostat_params.interbranch_vario_sill;
		}
		else if (paramType == "interbranch_vario_nugget:") {
		iss >> params.geostat_params.interbranch_vario_nugget;
		}
		else if (paramType == "interbranch_vario_model:") {
			std::string model;
			iss >> model;
			params.geostat_params.interbranch_vario_model = strip_optional_quotes(model);
		}
		else if (paramType == "intrabranch_vario_range:") {
		iss >> params.geostat_params.intrabranch_vario_range;
		}
		else if (paramType == "intrabranch_range_of_neighborhood:") {
		iss >> params.geostat_params.intrabranch_range_of_neighborhood;
		}
		else if (paramType == "intrabranch_vario_sill:") {
		iss >> params.geostat_params.intrabranch_vario_sill;
		}
		else if (paramType == "intrabranch_vario_nugget:") {
		iss >> params.geostat_params.intrabranch_vario_nugget;
		}
		else if (paramType == "intrabranch_vario_model:") {
			std::string model;
			iss >> model;
			params.geostat_params.intrabranch_vario_model = strip_optional_quotes(model);
		}
		else if (paramType == "number_max_of_neighborhood_points:") {
		iss >> params.geostat_params.number_max_of_neighborhood_points;
		}
		else if (paramType == "nb_points_interbranch:") {
		iss >> params.geostat_params.nb_points_interbranch;
		}
		else if (paramType == "proportion_interbranch:") {
		iss >> params.geostat_params.proportion_interbranch;
		}
		else if (paramType == "use_drift_zwt:") {
			std::string flag;
			iss >> flag;
			params.geostat_params.use_drift_zwt = parse_boolean(flag);
		}
		else if (paramType == "use_drift_curv:") {
			std::string flag;
			iss >> flag;
			params.geostat_params.use_drift_curv = parse_boolean(flag);
		}

		// Save objects

		else if (paramType == "create_vset_sampling:") {
		std::string flag;
		iss >> flag;
		params.create_vset_sampling = parse_boolean(flag);
		}
		else if (paramType == "create_nghb_graph:") {
		std::string flag;
		iss >> flag;
		params.create_nghb_graph = parse_boolean(flag);
		}
		else if (paramType == "create_nghb_graph_property:") {
		std::string flag;
		iss >> flag;
		params.create_nghb_graph_property = parse_boolean(flag);
		}
		else if (paramType == "create_solved_connectivity_matrix:") {
		std::string flag;
		iss >> flag;
		params.create_solved_connectivity_matrix = parse_boolean(flag);
		}
		else if (paramType == "create_grid:") {
		std::string flag;
		iss >> flag;
		params.create_grid = parse_boolean(flag);
		}

		// Other parameters

		else if (paramType == "fraction_karst_perm:") {
		iss >> params.fraction_karst_perm;
		}
		else if (paramType == "multiply_costs:") {
		std::string flag;
		iss >> flag;
		params.multiply_costs = parse_boolean(flag);
		}
		else if (paramType == "gamma:") {
		iss >> params.gamma;
		}
		else if (paramType == "vadose_cohesion:") {
		std::string flag;
		iss >> flag;
		params.vadose_cohesion = parse_boolean(flag);
		}
		else if (paramType == "vertical_distance_stretching_factor:") {
		iss >> params.vertical_distance_stretching_factor;
		}
		else if (paramType == "gradient_constraint_weight:") {
		iss >> params.gradient_constraint_weight;
		}
		else if (paramType == "outlet_selection_cost_factor:") {
		iss >> params.outlet_selection_cost_factor;
		}
	}
	inputFile.close();
	validate_parameter_values(params);
	return params;
}