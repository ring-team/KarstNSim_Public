/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/read_files.h"
#include <cerrno>
#include <cstring>

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



	void load_line(const std::string& file_name, const std::string& save_directory, Line& pline, std::vector<std::vector<std::vector<float>>>& properties)
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

		in.clear(); // clear eof flag
		in.seekg(0, std::ios::beg); // back to the start

		// Calculate the number of lines to be processed
		int data_lines = total_lines - 1; // minus 1 header

		// Resize properties based on the number of data lines
		properties.resize(data_lines / 2); // Each segment has 2 points

		// Prepare to read data
		std::getline(in, line); // Skip the header line

		bool count = false;
		std::vector<Segment> segs;
		Vector3 u1, u2;
		std::vector<float> prop(prop_size, 0.);
		std::vector<float> prop2(prop_size, 0.);

		int line_count = 0;
		while (std::getline(in, line) && line_count < data_lines) {
			std::istringstream iss2(line);
			int idx;
			float x, y, z;
			iss2 >> idx >> x >> y >> z;

			if (!count) {
				for (int k = 0; k < prop_size; k++) {
					iss2 >> prop[k];
				}
				u1 = Vector3(x, y, z);
			}
			else {
				for (int k = 0; k < prop_size; k++) {
					iss2 >> prop2[k];
				}
				u2 = Vector3(x, y, z);

				// Create properties array for this segment
				std::vector<std::vector<float>> props(prop_size, std::vector<float>(2));
				for (int k = 0; k < prop_size; k++) {
					props[k] = { prop[k], prop2[k] };
				}
				properties[line_count / 2] = props; // Store properties for this segment
				segs.push_back(Segment(u1, u2));
			}

			count = !count; // Toggle count
			line_count++;
		}

		pline = Line(segs);
		in.close();
	}


	void load_box(const std::string& file_name, const std::string& save_directory, Box& box, std::vector<std::vector<float>>& properties)
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
		in.clear(); // Clear eof flag
		in.seekg(0, std::ios::beg); // Back to the start

		// Skip headers and get property size and other characteristics
		std::getline(in, line); // Skip first line (header)

		std::string name;
		int prop_size;
		Vector3 basis, u, v, w;
		int nu, nv, nw;
		float x, y, z;

		// Property number
		std::getline(in, line);
		std::stringstream iss1(line);
		iss1 >> name >> prop_size;

		// Basis
		std::getline(in, line);
		std::stringstream iss2(line);
		iss2 >> name >> x >> y >> z;
		basis = Vector3(x, y, z);

		// u
		std::getline(in, line);
		std::stringstream iss3(line);
		iss3 >> name >> x >> y >> z;
		u = Vector3(x, y, z);

		// v
		std::getline(in, line);
		std::stringstream iss4(line);
		iss4 >> name >> x >> y >> z;
		v = Vector3(x, y, z);

		// w
		std::getline(in, line);
		std::stringstream iss5(line);
		iss5 >> name >> x >> y >> z;
		w = Vector3(x, y, z);

		// nu
		std::getline(in, line);
		std::stringstream iss6(line);
		iss6 >> name >> nu;

		// nv
		std::getline(in, line);
		std::stringstream iss7(line);
		iss7 >> name >> nv;

		// nw
		std::getline(in, line);
		std::stringstream iss8(line);
		iss8 >> name >> nw;

		box = Box(basis, u, v, w, nu, nv, nw);

		// Skip next line (header)
		std::getline(in, line);

		// Calculate the number of data lines
		int data_lines = total_lines - 9; // Total lines minus 9 header lines (1+8)

		// Resize properties based on the number of data lines
		properties.resize(data_lines);

		// Read the property values
		int line_count = 0;
		while (std::getline(in, line) && line_count < data_lines) {
			std::istringstream iss(line);
			int idx;
			std::vector<float> prop(prop_size, 0.);
			iss >> idx;
			for (int k = 0; k < prop_size; k++) {
				iss >> prop[k];
			}
			properties[line_count] = prop; // Store the properties in the array
			line_count++;
		}
		in.close();
	}
}