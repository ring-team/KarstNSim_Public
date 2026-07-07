/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/read_files.h"

namespace KarstNSim {
	void load_point(const std::string& file_name, const std::string& save_directory, Vector3& u, std::vector<float>& properties)
	{
		//(std::vector<Vector3>, std::vector<Triangle>)
		std::string line;
		std::string full_name = save_directory + "/" + file_name;
		std::ifstream in(full_name);

		//get property size

		std::getline(in, line); // skip first line, which is a header !
		std::getline(in, line);
		std::stringstream iss;
		iss << line;
		int nb_columns = 0;
		std::string value;
		while (iss >> value) nb_columns++;
		int prop_size = nb_columns - 4;
		// don't forget to reset to beginning, to not skip first line after header!
		in.clear();                 // clear fail and eof bits
		in.seekg(0, std::ios::beg); // back to the start!

		// go through lines

		std::getline(in, line); // skip first line, which is a header !
		while (std::getline(in, line)) { // get line
			std::istringstream iss2(line);
			// initialize all params of each line
			int idx;
			float x, y, z;
			std::vector<float> prop(prop_size, 0.);
			// and get their values
			iss2 >> idx >> x >> y >> z;
			// iterate on property size to get the rest
			for (int k = 0; k < prop_size; k++) {
				iss2 >> prop[k];
			}
			u = Vector3(x, y, z);
			properties = prop;
		}
		in.close();
	}

	void load_surface(const std::string& file_name, const std::string& save_directory, Surface& s, std::vector<std::vector<float>>& properties)
	{
		std::string line;
		std::string full_name = save_directory + "/" + file_name;
		std::ifstream in(full_name);

		if (!in.is_open()) {
			std::cerr << "Error opening file: " << full_name << std::endl;
			std::cerr << "Reason: " << strerror(errno) << std::endl;
			throw std::runtime_error("File not found: " + full_name);
		}

		// Skip headers and get the number of columns
		std::getline(in, line); // skip first line (header)
		std::getline(in, line); // read the second line for property size calculation
		std::stringstream iss;
		iss << line;
		int nb_columns = 0;
		std::string value;
		while (iss >> value) nb_columns++;
		int prop_size = nb_columns - 5; // Adjust according to actual properties

		in.clear(); // clear eof flag
		in.seekg(0, std::ios::beg); // back to the start

		// Prepare to read data
		std::vector<Vector3> nodes;
		std::vector<Triangle> trgls;
		std::string tag;

		std::getline(in, line); // Skip the header line

		int line_count = 0;
		while (std::getline(in, line)) {
			std::istringstream iss2(line);
			iss2 >> tag;

			if (tag == "VRTX") { // If we're reading the vertices
				int idx;
				float x, y, z;
				std::vector<float> prop(prop_size, 0.);
				iss2 >> idx >> x >> y >> z;

				for (int k = 0; k < prop_size; k++) {
					iss2 >> prop[k];
				}
				Vector3 v(x, y, z);
				nodes.push_back(v);
				properties.push_back(prop);
			}
			else if (tag == "TRGL") { // When vertices have been read, we now have to read the triangles
				int idx;
				int id1, id2, id3;
				iss2 >> idx >> id1 >> id2 >> id3;
				Triangle tgl(id1, id2, id3);
				trgls.push_back(tgl);
			}
			line_count++;
		}

		s = Surface(nodes, trgls);
		in.close();
	}

	void load_pointset(const std::string& file_name, const std::string& save_directory, std::vector<Vector3>& pset, std::vector<std::vector<float>>& properties)
	{
		std::string line;
		std::string full_name = save_directory + "/" + file_name;
		std::ifstream in(full_name);

		if (!in.is_open()) {
			std::cerr << "Error opening file: " << full_name << std::endl;
			std::cerr << "Reason: " << strerror(errno) << std::endl;
			throw std::runtime_error("File not found: " + full_name);
		}

		// Count total number of lines in the file
		int total_lines = 0;
		while (std::getline(in, line)) {
			total_lines++;
		}
		in.clear(); // clear eof flag
		in.seekg(0, std::ios::beg); // back to the start

		// Skip headers and get the number of columns
		std::getline(in, line); // skip first line (header)
		std::getline(in, line); // read the second line for property size calculation
		std::stringstream iss;
		iss << line;
		int nb_columns = 0;
		std::string value;
		while (iss >> value) nb_columns++;
		int prop_size = nb_columns - 4; // Adjust according to actual properties

		// Calculate the number of lines to be processed
		int data_lines = total_lines - 1; // minus 1 header

		in.clear(); // clear eof flag
		in.seekg(0, std::ios::beg); // back to the start

		// Resize properties based on the number of data lines
		properties.resize(data_lines);

		// Prepare to read data
		std::getline(in, line); // Skip the header line

		int line_count = 0;
		while (std::getline(in, line) && line_count < data_lines) {
			std::istringstream iss2(line);
			// Initialize all params of each line
			int idx;
			float x, y, z;
			std::vector<float> prop(prop_size, 0.);
			// Get their values
			iss2 >> idx >> x >> y >> z;

			// Iterate on property size to get the rest
			for (int k = 0; k < prop_size; k++) {
				iss2 >> prop[k];
			}
			Vector3 u(x, y, z);
			pset.push_back(u);
			properties[line_count] = prop; // Store the properties in the array

			line_count++;
		}

		in.close();
	}



	void load_line(
		const std::string& file_name,
		const std::string& save_directory,
		Line& pline,
		std::vector<std::vector<std::vector<float>>>& properties,
		std::vector<std::string>& property_names)
	{
		std::string line;
		std::string full_name = save_directory + "/" + file_name;
		std::ifstream in(full_name);

		if (!in.is_open()) {
			std::cerr << "Error opening file: " << full_name << std::endl;
			return;
		}

		int total_lines = 0;
		while (std::getline(in, line)) {
			total_lines++;
		}

		in.clear();
		in.seekg(0, std::ios::beg);

		if (total_lines <= 0) {
			throw std::runtime_error(
				"[load_line] Invalid line file: the file is empty."
			);
		}

		std::getline(in, line);

		std::stringstream header_stream(line);
		std::string token;
		std::vector<std::string> header_tokens;

		while (header_stream >> token) {
			header_tokens.push_back(token);
		}

		if (header_tokens.size() < 4) {
			throw std::runtime_error(
				"[load_line] Invalid line file header: expected at least Index, X, Y, and Z."
			);
		}

		property_names.assign(header_tokens.begin() + 4, header_tokens.end());

		const int prop_size = static_cast<int>(property_names.size());
		const int data_lines = total_lines - 1;

		if (data_lines < 0 || data_lines % 2 != 0) {
			throw std::runtime_error(
				"[load_line] Invalid line file: each segment must be represented by two endpoint rows."
			);
		}

		properties.clear();
		properties.resize(data_lines / 2);

		std::vector<Segment> segs;
		segs.reserve(data_lines / 2);

		Vector3 u1, u2;
		std::vector<float> prop(prop_size, 0.0f);
		std::vector<float> prop2(prop_size, 0.0f);

		bool first_endpoint = true;
		int line_count = 0;

		while (std::getline(in, line) && line_count < data_lines) {
			std::istringstream data_stream(line);

			int idx;
			float x, y, z;

			if (!(data_stream >> idx >> x >> y >> z)) {
				throw std::runtime_error(
					"[load_line] Invalid line file: endpoint row does not contain Index, X, Y, and Z."
				);
			}

			if (first_endpoint) {
				for (int k = 0; k < prop_size; k++) {
					if (!(data_stream >> prop[k])) {
						throw std::runtime_error(
							"[load_line] Invalid line file: missing property value on first segment endpoint."
						);
					}
				}

				u1 = Vector3(x, y, z);
			}
			else {
				for (int k = 0; k < prop_size; k++) {
					if (!(data_stream >> prop2[k])) {
						throw std::runtime_error(
							"[load_line] Invalid line file: missing property value on second segment endpoint."
						);
					}
				}

				u2 = Vector3(x, y, z);

				std::vector<std::vector<float>> props(prop_size, std::vector<float>(2));

				for (int k = 0; k < prop_size; k++) {
					props[k] = { prop[k], prop2[k] };
				}

				properties[line_count / 2] = props;
				segs.push_back(Segment(u1, u2));
			}

			first_endpoint = !first_endpoint;
			line_count++;
		}

		pline = Line(segs);
		in.close();
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