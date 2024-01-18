/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/read_files.h"

namespace KarstNSim {
	void load_point(const std::string& file_name, const std::string& save_directory, Vector3& u, std::vector<double>& properties)
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
		double value;
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
			double x, y, z;
			std::vector<double> prop(prop_size, 0.);
			// and get their values
			iss2 >> idx >> x >> y >> z;
			// iterate on property size to get the rest
			for (int k = 0; k < prop_size; k++) {
				iss2 >> prop[k];
			}
			u = Vector3(x, y, z);
			properties = prop;
		}
	}

	void load_surface(const std::string& file_name, const std::string& save_directory, Surface& s, std::vector<std::vector<double>>& properties)
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
		int prop_size = nb_columns - 5;
		// don't forget to reset to beginning, to not skip first line after header!
		in.clear();                 // clear fail and eof bits
		in.seekg(0, std::ios::beg); // back to the start!
		// go through lines

		std::vector<Vector3> nodes;
		std::vector<Triangle> trgls;
		std::string tag;

		std::getline(in, line); // skip first line, which is a header !
		while (std::getline(in, line)) { // get line
			std::istringstream iss2(line);
			iss2 >> tag;
			if (tag == "VRTX") { 			// if we're reading the vertices
				// initialize all params of each line
				int idx;
				double x, y, z;
				std::vector<double> prop(prop_size, 0.);
				// and get their values
				iss2 >> idx >> x >> y >> z;
				// iterate on property size to get the rest
				for (int k = 0; k < prop_size; k++) {
					iss2 >> prop[k];
				}
				Vector3 v(x, y, z);
				nodes.push_back(v);
				properties.push_back(prop);
			}
			else if (tag == "TRGL") { 		// when vertices have been read, we now have to read the triangles
				int idx;
				int id1, id2, id3;
				iss2 >> idx >> id1 >> id2 >> id3;
				Triangle tgl(id1, id2, id3);
				trgls.push_back(tgl);
			}
		}
		s = Surface(nodes, trgls);
	}

	void load_pointset(const std::string& file_name, const std::string& save_directory, std::vector<Vector3>& pset, std::vector<std::vector<double>>& properties)
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
		double value;
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
			double x, y, z;
			std::vector<double> prop(prop_size, 0.);
			// and get their values
			iss2 >> idx >> x >> y >> z;
			// iterate on property size to get the rest
			for (int k = 0; k < prop_size; k++) {
				iss2 >> prop[k];
			}
			Vector3 u(x, y, z);
			pset.push_back(u);
			properties.push_back(prop);
		}
	}

	void load_line(const std::string& file_name, const std::string& save_directory, Line& pline, std::vector<std::vector<std::vector<double>>>& properties)
	{
		//(std::vector<Vector3>, std::vector<Triangle>)
		std::string line;
		std::string full_name = save_directory + "/" + file_name;
		std::ifstream in(full_name);
		std::vector<Segment> segs;
		//get property size

		std::getline(in, line); // skip first line, which is a header !
		std::getline(in, line);
		std::stringstream iss;
		iss << line;
		int nb_columns = 0;
		double value;
		while (iss >> value) nb_columns++;
		int prop_size = nb_columns - 4;
		// don't forget to reset to beginning, to not skip first line after header!
		in.clear();                 // clear fail and eof bits
		in.seekg(0, std::ios::beg); // back to the start!

		// go through lines
		bool count = false;
		std::getline(in, line); // skip first line, which is a header !
		Vector3 u1;
		Vector3 u2;
		std::vector<double> prop(prop_size, 0.);
		std::vector<double> prop2(prop_size, 0.);

		while (std::getline(in, line)) { // get line
			std::istringstream iss2(line);
			// initialize all params of each line
			int idx;
			double x, y, z;

			// and get their values
			iss2 >> idx >> x >> y >> z;
			// iterate on property size to get the rest
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
			}

			if (count) {
				std::vector<std::vector<double>> props;
				for (int k = 0; k < prop_size; k++) {
					props.push_back({ prop[k],prop2[k] });
				}
				properties.push_back(props);
				Segment seg(u1, u2);
				segs.push_back(seg);
				count = false;
			}
			else {
				count = true;
			}
		}
		pline = Line(segs);
	}

	void load_box(const std::string& file_name, const std::string& save_directory, Box& box, std::vector<std::vector<double>>& properties)
	{
		//(std::vector<Vector3>, std::vector<Triangle>)
		std::string line;
		std::string full_name = save_directory + "/" + file_name;
		std::ifstream in(full_name);

		//get property size and other characteristics

		std::getline(in, line); // skip first line, which is a header !

		std::string name;
		int prop_size;
		Vector3 basis;
		Vector3 u;
		Vector3 v;
		Vector3 w;
		int nu;
		int nv;
		int nw;
		double x, y, z;

		// prop number
		std::getline(in, line);
		std::stringstream iss1(line);
		iss1 >> name >> prop_size;

		// basis
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

		// skip next line (header)
		std::getline(in, line);

		// go through next lines for property values

		while (std::getline(in, line)) { // get line
			std::istringstream iss(line);
			// initialize all params of each line
			int idx;
			std::vector<double> prop(prop_size, 0.);
			// and get their values
			iss >> idx;
			// iterate on property size to get the rest
			for (int k = 0; k < prop_size; k++) {
				iss >> prop[k];
			}
			properties.push_back(prop);
		}
	}
}