/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/read_files.h"
#include <utility>

namespace KarstNSim {
	void load_point(
		const std::string& file_name,
		const std::string& save_directory,
		Vector3& u,
		std::vector<float>& properties)
	{
		const std::string full_name = save_directory + "/" + file_name;
		std::ifstream input(full_name);
		if (!input.is_open()) {
			throw std::runtime_error("[load_point] Cannot open file: " + full_name);
		}

		std::string line;
		if (!std::getline(input, line)) {
			throw std::runtime_error("[load_point] Empty file: " + full_name);
		}

		std::istringstream header_stream(line);
		std::vector<std::string> header_tokens;
		std::string token;
		while (header_stream >> token) {
			header_tokens.push_back(token);
		}
		if (header_tokens.size() < 4u) {
			throw std::runtime_error(
				"[load_point] Invalid header in " + full_name +
				": expected at least Index, X, Y, and Z."
			);
		}
		const std::size_t property_count = header_tokens.size() - 4u;

		bool point_read = false;
		std::size_t line_number = 1;
		while (std::getline(input, line)) {
			++line_number;
			if (line.find_first_not_of(" \t\r\n") == std::string::npos) {
				continue;
			}
			if (point_read) {
				throw std::runtime_error(
					"[load_point] More than one data row was found in " + full_name +
					"; a point file must contain exactly one point."
				);
			}

			std::istringstream data_stream(line);
			int index = 0;
			float x = 0.0f;
			float y = 0.0f;
			float z = 0.0f;
			if (!(data_stream >> index >> x >> y >> z)) {
				throw std::runtime_error(
					"[load_point] Invalid point row in " + full_name + " at line " +
					std::to_string(line_number) + "."
				);
			}

			std::vector<float> parsed_properties(property_count, 0.0f);
			for (std::size_t i = 0; i < property_count; ++i) {
				if (!(data_stream >> parsed_properties[i])) {
					throw std::runtime_error(
						"[load_point] Missing property value in " + full_name +
						" at line " + std::to_string(line_number) + "."
					);
				}
			}
			std::string extra;
			if (data_stream >> extra) {
				throw std::runtime_error(
					"[load_point] Unexpected extra token '" + extra + "' in " +
					full_name + " at line " + std::to_string(line_number) + "."
				);
			}

			u = Vector3(x, y, z);
			properties = std::move(parsed_properties);
			point_read = true;
		}

		if (!point_read) {
			throw std::runtime_error("[load_point] No point row was found in " + full_name + ".");
		}
	}

	void load_surface(
		const std::string& file_name,
		const std::string& save_directory,
		Surface& s,
		std::vector<std::vector<float>>& properties)
	{
		const std::string full_name = save_directory + "/" + file_name;
		std::ifstream input(full_name);
		if (!input.is_open()) {
			throw std::runtime_error("[load_surface] Cannot open file: " + full_name);
		}

		std::string line;
		if (!std::getline(input, line)) {
			throw std::runtime_error("[load_surface] Empty file: " + full_name);
		}

		std::istringstream header_stream(line);
		std::vector<std::string> header_tokens;
		std::string token;
		while (header_stream >> token) {
			header_tokens.push_back(token);
		}
		if (header_tokens.size() < 5u) {
			throw std::runtime_error(
				"[load_surface] Invalid header in " + full_name +
				": expected at least Type, Index, X, Y, and Z."
			);
		}
		const std::size_t property_count = header_tokens.size() - 5u;

		std::vector<Vector3> nodes;
		std::vector<Triangle> triangles;
		properties.clear();

		std::size_t line_number = 1;
		while (std::getline(input, line)) {
			++line_number;
			if (line.find_first_not_of(" \t\r\n") == std::string::npos) {
				continue;
			}

			std::istringstream data_stream(line);
			std::string tag;
			data_stream >> tag;

			if (tag == "VRTX") {
				int index = 0;
				float x = 0.0f;
				float y = 0.0f;
				float z = 0.0f;
				if (!(data_stream >> index >> x >> y >> z)) {
					throw std::runtime_error(
						"[load_surface] Invalid VRTX row in " + full_name + " at line " +
						std::to_string(line_number) + "."
					);
				}

				std::vector<float> vertex_properties(property_count, 0.0f);
				for (std::size_t i = 0; i < property_count; ++i) {
					if (!(data_stream >> vertex_properties[i])) {
						throw std::runtime_error(
							"[load_surface] Missing vertex property in " + full_name +
							" at line " + std::to_string(line_number) + "."
						);
					}
				}
				std::string extra;
				if (data_stream >> extra) {
					throw std::runtime_error(
						"[load_surface] Unexpected extra token '" + extra + "' in " +
						full_name + " at line " + std::to_string(line_number) + "."
					);
				}

				nodes.emplace_back(x, y, z);
				properties.push_back(std::move(vertex_properties));
			}
			else if (tag == "TRGL") {
				int index = 0;
				int id1 = 0;
				int id2 = 0;
				int id3 = 0;
				if (!(data_stream >> index >> id1 >> id2 >> id3)) {
					throw std::runtime_error(
						"[load_surface] Invalid TRGL row in " + full_name + " at line " +
						std::to_string(line_number) + "."
					);
				}
				std::string extra;
				if (data_stream >> extra) {
					throw std::runtime_error(
						"[load_surface] Unexpected extra token '" + extra + "' in " +
						full_name + " at line " + std::to_string(line_number) + "."
					);
				}
				triangles.emplace_back(id1, id2, id3);
			}
			else {
				throw std::runtime_error(
					"[load_surface] Unknown row tag '" + tag + "' in " + full_name +
					" at line " + std::to_string(line_number) + "."
				);
			}
		}

		if (nodes.empty()) {
			throw std::runtime_error("[load_surface] No VRTX row was found in " + full_name + ".");
		}
		if (triangles.empty()) {
			throw std::runtime_error("[load_surface] No TRGL row was found in " + full_name + ".");
		}

		for (std::size_t i = 0; i < triangles.size(); ++i) {
			for (int vertex = 0; vertex < 3; ++vertex) {
				const int node_index = triangles[i].point(vertex);
				if (node_index < 0 || node_index >= static_cast<int>(nodes.size())) {
					throw std::runtime_error(
						"[load_surface] Triangle #" + std::to_string(i + 1) +
						" in " + full_name + " references invalid node index " +
						std::to_string(node_index) + "."
					);
				}
			}
		}

		Surface loaded_surface(nodes, triangles);

		// Surface construction can remove degenerate triangles. A file may therefore
		// contain TRGL rows while still producing no usable triangle.
		if (loaded_surface.get_nb_trgls() == 0) {
			throw std::runtime_error(
				"[load_surface] All triangles in " + full_name +
				" are degenerate in map view; the file does not define "
				"a usable triangulated surface."
			);
		}

		s = std::move(loaded_surface);
	}

	void load_pointset(
		const std::string& file_name,
		const std::string& save_directory,
		std::vector<Vector3>& pset,
		std::vector<std::vector<float>>& properties)
	{
		const std::string full_name = save_directory + "/" + file_name;
		std::ifstream input(full_name);
		if (!input.is_open()) {
			throw std::runtime_error("[load_pointset] Cannot open file: " + full_name);
		}

		std::string line;
		if (!std::getline(input, line)) {
			throw std::runtime_error("[load_pointset] Empty file: " + full_name);
		}

		std::istringstream header_stream(line);
		std::vector<std::string> header_tokens;
		std::string token;
		while (header_stream >> token) {
			header_tokens.push_back(token);
		}
		if (header_tokens.size() < 4u) {
			throw std::runtime_error(
				"[load_pointset] Invalid header in " + full_name +
				": expected at least Index, X, Y, and Z."
			);
		}
		const std::size_t property_count = header_tokens.size() - 4u;

		pset.clear();
		properties.clear();
		std::size_t line_number = 1;
		while (std::getline(input, line)) {
			++line_number;
			if (line.find_first_not_of(" \t\r\n") == std::string::npos) {
				continue;
			}

			std::istringstream data_stream(line);
			int index = 0;
			float x = 0.0f;
			float y = 0.0f;
			float z = 0.0f;
			if (!(data_stream >> index >> x >> y >> z)) {
				throw std::runtime_error(
					"[load_pointset] Invalid point row in " + full_name + " at line " +
					std::to_string(line_number) + "."
				);
			}

			std::vector<float> point_properties(property_count, 0.0f);
			for (std::size_t i = 0; i < property_count; ++i) {
				if (!(data_stream >> point_properties[i])) {
					throw std::runtime_error(
						"[load_pointset] Missing property value in " + full_name +
						" at line " + std::to_string(line_number) + "."
					);
				}
			}
			std::string extra;
			if (data_stream >> extra) {
				throw std::runtime_error(
					"[load_pointset] Unexpected extra token '" + extra + "' in " +
					full_name + " at line " + std::to_string(line_number) + "."
				);
			}

			pset.emplace_back(x, y, z);
			properties.push_back(std::move(point_properties));
		}

		if (pset.empty()) {
			throw std::runtime_error("[load_pointset] No point row was found in " + full_name + ".");
		}
	}



	void load_line(
		const std::string& file_name,
		const std::string& save_directory,
		Line& pline,
		std::vector<std::vector<std::vector<float>>>& properties,
		std::vector<std::string>& property_names)
	{
		const std::string full_name = save_directory + "/" + file_name;
		std::ifstream input(full_name);
		if (!input.is_open()) {
			throw std::runtime_error("[load_line] Cannot open file: " + full_name);
		}

		std::string line;
		if (!std::getline(input, line)) {
			throw std::runtime_error("[load_line] Empty file: " + full_name);
		}

		std::istringstream header_stream(line);
		std::vector<std::string> header_tokens;
		std::string token;
		while (header_stream >> token) {
			header_tokens.push_back(token);
		}
		if (header_tokens.size() < 4u) {
			throw std::runtime_error(
				"[load_line] Invalid header in " + full_name +
				": expected at least Index, X, Y, and Z."
			);
		}

		property_names.assign(header_tokens.begin() + 4, header_tokens.end());
		const std::size_t property_count = property_names.size();

		struct EndpointRow {
			Vector3 point;
			std::vector<float> properties;
		};
		std::vector<EndpointRow> endpoints;

		std::size_t line_number = 1;
		while (std::getline(input, line)) {
			++line_number;
			if (line.find_first_not_of(" \t\r\n") == std::string::npos) {
				continue;
			}

			std::istringstream data_stream(line);
			int index = 0;
			float x = 0.0f;
			float y = 0.0f;
			float z = 0.0f;
			if (!(data_stream >> index >> x >> y >> z)) {
				throw std::runtime_error(
					"[load_line] Invalid endpoint row in " + full_name + " at line " +
					std::to_string(line_number) + "."
				);
			}

			std::vector<float> endpoint_properties(property_count, 0.0f);
			for (std::size_t i = 0; i < property_count; ++i) {
				if (!(data_stream >> endpoint_properties[i])) {
					throw std::runtime_error(
						"[load_line] Missing property value in " + full_name +
						" at line " + std::to_string(line_number) + "."
					);
				}
			}
			std::string extra;
			if (data_stream >> extra) {
				throw std::runtime_error(
					"[load_line] Unexpected extra token '" + extra + "' in " +
					full_name + " at line " + std::to_string(line_number) + "."
				);
			}

			endpoints.push_back({ Vector3(x, y, z), std::move(endpoint_properties) });
		}

		if (endpoints.empty()) {
			throw std::runtime_error("[load_line] No endpoint row was found in " + full_name + ".");
		}
		if (endpoints.size() % 2u != 0u) {
			throw std::runtime_error(
				"[load_line] Invalid endpoint count in " + full_name +
				": each segment requires exactly two consecutive endpoint rows."
			);
		}

		std::vector<Segment> segments;
		segments.reserve(endpoints.size() / 2u);
		properties.clear();
		properties.resize(endpoints.size() / 2u);

		for (std::size_t segment_index = 0; segment_index < endpoints.size() / 2u;
			++segment_index) {
			const EndpointRow& first = endpoints[2u * segment_index];
			const EndpointRow& second_endpoint = endpoints[2u * segment_index + 1u];
			segments.emplace_back(first.point, second_endpoint.point);

			std::vector<std::vector<float>> segment_properties(
				property_count, std::vector<float>(2u, 0.0f));
			for (std::size_t property_index = 0; property_index < property_count;
				++property_index) {
				segment_properties[property_index][0] = first.properties[property_index];
				segment_properties[property_index][1] =
					second_endpoint.properties[property_index];
			}
			properties[segment_index] = std::move(segment_properties);
		}

		pline = Line(segments);
	}

	void load_line(
		const std::string& file_name,
		const std::string& save_directory,
		Line& pline,
		std::vector<std::vector<std::vector<float>>>& properties)
	{
		std::vector<std::string> property_names;

		load_line(
			file_name,
			save_directory,
			pline,
			properties,
			property_names
		);
	}

	void load_box(
		const std::string& file_name,
		const std::string& save_directory,
		Box& box,
		std::vector<std::vector<float>>& properties)
	{
		std::string line;
		const std::string full_name = save_directory + "/" + file_name;
		std::ifstream in(full_name);

		if (!in.is_open()) {
			throw std::runtime_error("[load_box] Error opening file: " + full_name);
		}

		// Skip the format header.
		if (!std::getline(in, line)) {
			throw std::runtime_error("[load_box] Empty box file: " + full_name);
		}

		std::string name;
		int prop_size = 0;
		Vector3 basis, u, v, w;
		int nu = 0;
		int nv = 0;
		int nw = 0;
		float x = 0.0f;
		float y = 0.0f;
		float z = 0.0f;

		// Property number.
		if (!std::getline(in, line)) {
			throw std::runtime_error("[load_box] Missing property count in file: " + full_name);
		}
		std::stringstream iss1(line);
		if (!(iss1 >> name >> prop_size) || prop_size <= 0) {
			throw std::runtime_error("[load_box] Invalid property count in file: " + full_name);
		}

		// Basis.
		if (!std::getline(in, line)) {
			throw std::runtime_error("[load_box] Missing basis vector in file: " + full_name);
		}
		std::stringstream iss2(line);
		if (!(iss2 >> name >> x >> y >> z)) {
			throw std::runtime_error("[load_box] Invalid basis vector in file: " + full_name);
		}
		basis = Vector3(x, y, z);

		// U vector.
		if (!std::getline(in, line)) {
			throw std::runtime_error("[load_box] Missing u vector in file: " + full_name);
		}
		std::stringstream iss3(line);
		if (!(iss3 >> name >> x >> y >> z)) {
			throw std::runtime_error("[load_box] Invalid u vector in file: " + full_name);
		}
		u = Vector3(x, y, z);

		// V vector.
		if (!std::getline(in, line)) {
			throw std::runtime_error("[load_box] Missing v vector in file: " + full_name);
		}
		std::stringstream iss4(line);
		if (!(iss4 >> name >> x >> y >> z)) {
			throw std::runtime_error("[load_box] Invalid v vector in file: " + full_name);
		}
		v = Vector3(x, y, z);

		// W vector.
		if (!std::getline(in, line)) {
			throw std::runtime_error("[load_box] Missing w vector in file: " + full_name);
		}
		std::stringstream iss5(line);
		if (!(iss5 >> name >> x >> y >> z)) {
			throw std::runtime_error("[load_box] Invalid w vector in file: " + full_name);
		}
		w = Vector3(x, y, z);

		// Grid dimensions.
		if (!std::getline(in, line)) {
			throw std::runtime_error("[load_box] Missing nu value in file: " + full_name);
		}
		std::stringstream iss6(line);
		if (!(iss6 >> name >> nu) || nu <= 0) {
			throw std::runtime_error("[load_box] Invalid nu value in file: " + full_name);
		}

		if (!std::getline(in, line)) {
			throw std::runtime_error("[load_box] Missing nv value in file: " + full_name);
		}
		std::stringstream iss7(line);
		if (!(iss7 >> name >> nv) || nv <= 0) {
			throw std::runtime_error("[load_box] Invalid nv value in file: " + full_name);
		}

		if (!std::getline(in, line)) {
			throw std::runtime_error("[load_box] Missing nw value in file: " + full_name);
		}
		std::stringstream iss8(line);
		if (!(iss8 >> name >> nw) || nw <= 0) {
			throw std::runtime_error("[load_box] Invalid nw value in file: " + full_name);
		}

		box = Box(basis, u, v, w, nu, nv, nw);

		// Skip the data column header.
		if (!std::getline(in, line)) {
			throw std::runtime_error("[load_box] Missing data header in file: " + full_name);
		}

		const int expected_data_lines = nu * nv * nw;
		properties.clear();
		properties.resize(static_cast<size_t>(expected_data_lines));

		int line_count = 0;
		while (std::getline(in, line) && line_count < expected_data_lines) {
			if (line.empty()) {
				continue;
			}

			std::istringstream iss(line);
			int idx = -1;
			std::vector<float> prop(static_cast<size_t>(prop_size), 0.0f);

			if (!(iss >> idx)) {
				throw std::runtime_error("[load_box] Invalid cell index in file: " + full_name);
			}

			for (int k = 0; k < prop_size; ++k) {
				if (!(iss >> prop[static_cast<size_t>(k)])) {
					throw std::runtime_error("[load_box] Invalid property value in file: " + full_name);
				}
			}

			properties[static_cast<size_t>(line_count)] = prop;
			++line_count;
		}

		if (line_count != expected_data_lines) {
			throw std::runtime_error(
				"[load_box] Box property count mismatch in file: " + full_name +
				". Expected " + std::to_string(expected_data_lines) +
				" cells, read " + std::to_string(line_count) + "."
			);
		}

		in.close();
	}

	KarstNSim::InputGraph translate_input_graph(const std::string& file_name, const std::string& save_directory) {

		std::string full_name = save_directory + "/" + file_name;
		std::ifstream in(full_name);

		if (!in.is_open()) {
			throw std::runtime_error("Error opening input graph file: " + full_name);
		}

		std::string line;
		std::string tag;

		int nb_nodes = -1;
		int nb_edges = -1;

		std::vector<Vector3> nodes;
		std::vector<Vector2i> edges;

		// -----------------------------
		// Read NODES section header
		// -----------------------------
		while (std::getline(in, line)) {
			if (line.empty() || line[0] == '#') {
				continue;
			}
			std::istringstream iss(line);
			iss >> tag;
			if (tag == "NODES") {
				iss >> nb_nodes;
				break;
			}
		}

		if (nb_nodes < 0) {
			throw std::runtime_error("Invalid input graph file: missing 'NODES <count>' section.");
		}

		nodes.resize(nb_nodes);

		// -----------------------------
		// Read node list
		// -----------------------------
		for (int i = 0; i < nb_nodes; ++i) {
			if (!std::getline(in, line)) {
				throw std::runtime_error("Invalid input graph file: unexpected end of file while reading nodes.");
			}

			if (line.empty() || line[0] == '#') {
				--i;
				continue;
			}

			std::istringstream iss(line);
			int node_id;
			float x, y, z;
			iss >> node_id >> x >> y >> z;

			if (!iss) {
				throw std::runtime_error("Invalid input graph file: malformed node line: " + line);
			}
			if (node_id < 0 || node_id >= nb_nodes) {
				throw std::runtime_error("Invalid input graph file: node id out of range: " + std::to_string(node_id));
			}

			nodes[node_id] = Vector3(x, y, z);
		}

		// -----------------------------
		// Read EDGES section header
		// -----------------------------
		while (std::getline(in, line)) {
			if (line.empty() || line[0] == '#') {
				continue;
			}
			std::istringstream iss(line);
			iss >> tag;
			if (tag == "EDGES") {
				iss >> nb_edges;
				break;
			}
		}

		if (nb_edges < 0) {
			throw std::runtime_error("Invalid input graph file: missing 'EDGES <count>' section.");
		}

		edges.reserve(nb_edges);

		// -----------------------------
		// Read edge list
		// -----------------------------
		for (int i = 0; i < nb_edges; ++i) {
			if (!std::getline(in, line)) {
				throw std::runtime_error("Invalid input graph file: unexpected end of file while reading edges.");
			}

			if (line.empty() || line[0] == '#') {
				--i;
				continue;
			}

			std::istringstream iss(line);
			int a, b;
			iss >> a >> b;

			if (!iss) {
				throw std::runtime_error("Invalid input graph file: malformed edge line: " + line);
			}
			if (a < 0 || a >= nb_nodes || b < 0 || b >= nb_nodes) {
				throw std::runtime_error("Invalid input graph file: edge endpoint out of range: " + line);
			}

			edges.push_back(Vector2i(a, b));
		}

		in.close();

		return KarstNSim::InputGraph(nodes, edges);
	}
}