/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/write_files.h"

namespace KarstNSim {

	void save_point(const std::string& file_name, const std::string& save_directory, Vector3 u, std::vector<std::string> property_names, std::vector<double> properties)
	{

		std::ofstream out;
		out.open(save_directory + '/' + file_name);
		if (out.is_open() == false)
		{
			std::cout << "Cannot save skeleton to file (nodes): " << file_name << std::endl;
			return;
		}

		out << "Index	X	Y	Z";
		for (int j = 0; j < property_names.size(); j++) { // iterate on properties (if any)
			out << "	" << property_names[j]; // get name of property j
		}
		out << "\n";

		out << 1 << "	" << u[0] << "	" << u[1] << "	" << u[2];
		if (properties.size()) {
			for (int j = 0; j < int(properties.size()); j++) { // iterate on properties (if any)
				out << "	" << std::fixed << std::setprecision(10) <<  properties[j];
			}
			out << "\n";
		}
		out.flush();
		out.close();
	}

	void save_surface(const std::string& file_name, const std::string& save_directory, Surface s, std::vector<std::string> property_names, std::vector<std::vector<double>> properties)
	{
		int nb_pts = s.get_nb_pts();
		int nb_trgls = s.get_nb_trgls();

		std::ofstream out;
		out.open(save_directory + '/' + file_name);
		if (out.is_open() == false)
		{
			std::cout << "Cannot save skeleton to file (nodes): " << file_name << std::endl;
			return;
		}

		out << "Type	Index	X	Y	Z";
		for (int j = 0; j < property_names.size(); j++) { // iterate on properties (if any)
			out << "	" << property_names[j]; // get name of property j
		}
		out << "\n";

		for (int i = 0; i < nb_pts; i++)
		{
			Vector3 v = s.get_node(i);
			out << "VRTX " << i + 1 << "	" << v[0] << "	" << v[1] << "	" << v[2];
			if (properties.size()) {
				for (int j = 0; j < properties[0].size(); j++) { // iterate on properties (if any)
					out << "	" << std::fixed << std::setprecision(10) << properties[i][j];
				}
			}
			out << "\n";
		}

		for (int i = 0; i < nb_trgls; i++) {
			Triangle tr = s.get_triangle(i);
			out << "TRGL " << i + 1 << "	" << tr.point(0) << "	" << tr.point(1) << "	" << tr.point(2);
			out << "\n";
		}
		out.flush();
		out.close();
	}

	void save_pointset(const std::string& file_name, const std::string& save_directory, std::vector<Vector3> pset, std::vector<std::string> property_names, std::vector<std::vector<double>> properties)
	{

		int nb_pts = int(pset.size());
		std::ofstream out;
		out.open(save_directory + '/' + file_name);
		if (out.is_open() == false)
		{
			std::cout << "Cannot save skeleton to file (nodes): " << file_name << std::endl;
			return;
		}
			out << "Index	X	Y	Z";
		for (int j = 0; j < property_names.size(); j++) { // iterate on properties (if any)
			out << "	" << property_names[j]; // get name of property j
		}
		out << "\n";
			for (int i = 0; i < nb_pts; i++)
			{

			Vector3 v = pset[i];
			out << i + 1 << "	" << v[0] << "	" << v[1] << "	" << v[2];
			if (properties.size()) {
				for (int j = 0; j < properties[0].size(); j++) { // iterate on properties (if any)
					out << "	" << std::fixed << std::setprecision(10) << properties[i][j];
				}
			}
			out << "\n";
		}
		out.flush();
		out.close();
	}

	void save_line(const std::string& file_name, const std::string& save_directory, Line pline, std::vector<std::string> property_names, std::vector<std::vector<std::vector<double>>> properties)
	{

		int nb_segs = pline.get_nb_segs();

		std::ofstream out;
		out.open(save_directory + '/' + file_name);
		if (out.is_open() == false)
		{
			std::cout << "Cannot save skeleton to file (nodes): " << file_name << std::endl;
			return;
		}

		out << "Index	X	Y	Z";
		for (int j = 0; j < property_names.size(); j++) { // iterate on properties (if any)
			out << "	" << property_names[j]; // get name of property j
		}
		out << "\n";

		for (int i = 0; i < nb_segs; i++)
		{

			Segment seg = pline.get_seg(i);
			Vector3 start = seg.start();
			Vector3 end = seg.end();
			out << i + 1 << "	" << start[0] << "	" << start[1] << "	" << start[2];
			if (properties.size()) {
				for (int j = 0; j < properties[0].size(); j++) { // iterate on properties (if any)
					out << "	" << std::fixed << std::setprecision(10) <<  properties[i][j][0];
				}
			}
			out << "\n";
			out << i + 1 << "	" << end[0] << "	" << end[1] << "	" << end[2];
			if (properties.size()) {
				for (int j = 0; j < properties[0].size(); j++) { // iterate on properties (if any)
					out << "	" << std::fixed << std::setprecision(10) <<properties[i][j][1];
				}
			}
			out << "\n";
		}
		out.flush();
		out.close();
	}

	void save_box(const std::string& file_name, const std::string& save_directory, Box box, std::vector<std::string> property_names, std::vector<std::vector<double>> properties)
	{

		std::ofstream out;
		out.open(save_directory + '/' + file_name);
		if (out.is_open() == false)
		{
			std::cout << "Cannot save skeleton to file (nodes): " << file_name << std::endl;
			return;
		}
		int nb_prop = int(property_names.size());
		Vector3 basis = box.get_basis();
		Vector3 u = box.get_u();
		Vector3 v = box.get_v();
		Vector3 w = box.get_w();
		int nu = box.get_nu();
		int nv = box.get_nv();
		int nw = box.get_nw();

		out << "Parameters" << "\n";
		out << "number_properties	" << nb_prop << "\n";
		out << "basis	" << basis.x << "	" << basis.y << "	" << basis.z << "\n";
		out << "u	" << u.x << "	" << u.y << "	" << u.z << "\n";
		out << "v	" << v.x << "	" << v.y << "	" << v.z << "\n";
		out << "w	" << w.x << "	" << w.y << "	" << w.z << "\n";
		out << "nu	" << nu << "\n";
		out << "nv	" << nv << "\n";
		out << "nw	" << nw << "\n";
		out << "Index";
		for (int j = 0; j < property_names.size(); j++) { // iterate on properties (if any)
			out << "	" << property_names[j]; // get name of property j
		}
		out << "\n";

		for (int wk = 0; wk < nw; wk++) { // w is always external loop, v middle loop and u internal loop !!!!
			for (int vk = 0; vk < nv; vk++) {
				for (int uk = 0; uk < nu; uk++) {
					int index = uk + vk * nu + wk * nu * nv;
					out << index;
					for (int p = 0; p < int(properties[0].size()); p++) { // iterate on properties (if any)
						out << "	" << std::fixed << std::setprecision(10) << properties[index][p];
					}
					out << "\n";
				}
			}
		}

		out.flush();
		out.close();
	}
}