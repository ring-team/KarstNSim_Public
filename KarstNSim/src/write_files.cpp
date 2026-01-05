/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/write_files.h"
#include <filesystem>

namespace KarstNSim {

	static bool directory_exists(const std::string& directory) {
		std::error_code ec;
		return std::filesystem::is_directory(directory, ec) && !ec;
	}

	static void create_directory(const std::string& directory) {
		std::error_code ec;
		if (!std::filesystem::create_directories(directory, ec) && ec) {
			std::cerr << "Error creating directory: " << directory << " (" << ec.message() << ")" << std::endl;
		}
	}

	static void ensure_directory_exists(const std::string& directory) {
		if (!directory_exists(directory)) {
			create_directory(directory);
		}
	}

	// Generate a unique filename by appending (0), (1), etc. if needed
	std::filesystem::path make_unique_filename(const std::string& base_filename, const std::string& save_directory) {
		std::filesystem::path dirPath(save_directory);
		std::filesystem::path filePath(base_filename);
		std::filesystem::path fullPath = dirPath / filePath;

		if (!std::filesystem::exists(fullPath)) {
			return fullPath;  // Already unique
		}

		std::string baseName = filePath.stem().string();
		std::string extension = filePath.extension().string();

		int copyIndex = 0;
		std::filesystem::path uniquePath;
		do {
			std::string uniqueFilename = baseName + "(" + std::to_string(copyIndex) + ")" + extension;
			uniquePath = dirPath / uniqueFilename;
			copyIndex++;
		} while (std::filesystem::exists(uniquePath));
		
		return uniquePath;
	}

	// Save the point data to a unique file
	void save_point(const std::string& file_name, const std::string& save_directory, Vector3 u, std::vector<std::string> property_names, std::vector<float> properties) {

		ensure_directory_exists(save_directory);

		// Generate a unique file path
		std::filesystem::path filePath = make_unique_filename(file_name, save_directory);
		std::ofstream out(filePath);
		if (!out.is_open()) {
			std::cout << "Cannot save skeleton to file (nodes): " << filePath.filename() << std::endl;
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
				out << "	" << std::fixed << std::setprecision(10) << properties[j];
			}
			out << "\n";
		}
		out.flush();
		out.close();
	}

	void save_surface(const std::string& file_name, const std::string& save_directory, Surface s, std::vector<std::string> property_names, std::vector<std::vector<float>> properties)
	{

		ensure_directory_exists(save_directory);

		// Generate a unique file path
		std::filesystem::path filePath = make_unique_filename(file_name, save_directory);

		int nb_pts = s.get_nb_pts();
		int nb_trgls = s.get_nb_trgls();

		std::ofstream out(filePath);
		if (!out.is_open())
		{
			std::cout << "Cannot save skeleton to file (nodes): " << filePath.filename() << std::endl;
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

	void save_pointset(const std::string& file_name, const std::string& save_directory, std::vector<Vector3> pset, std::vector<std::string> property_names, std::vector<std::vector<float>> properties)
	{

		ensure_directory_exists(save_directory);

		// Generate a unique file path
		std::filesystem::path filePath = make_unique_filename(file_name, save_directory);

		int nb_pts = int(pset.size());
		std::ofstream out(filePath);
		if (!out.is_open())
		{
			std::cout << "Cannot save skeleton to file (nodes): " << filePath.filename() << std::endl;
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

	void save_line(const std::string& file_name, const std::string& save_directory, Line pline, std::vector<std::string> property_names, std::vector<std::vector<std::vector<float>>> properties)
	{

		ensure_directory_exists(save_directory);

		// Generate a unique file path
		std::filesystem::path filePath = make_unique_filename(file_name, save_directory);

		int nb_segs = pline.get_nb_segs();

		std::ofstream out(filePath);
		if (!out.is_open())
		{
			std::cout << "Cannot save skeleton to file (nodes): " << filePath.filename() << std::endl;
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
					out << "	" << std::fixed << std::setprecision(10) << properties[i][j][0];
				}
			}
			out << "\n";
			out << i + 1 << "	" << end[0] << "	" << end[1] << "	" << end[2];
			if (properties.size()) {
				for (int j = 0; j < properties[0].size(); j++) { // iterate on properties (if any)
					out << "	" << std::fixed << std::setprecision(10) << properties[i][j][1];
				}
			}
			out << "\n";
		}
		out.flush();
		out.close();
	}

	void save_connectivity_matrix(const std::string& file_name, const std::string& save_directory, Array2D<int> matrix) {

		ensure_directory_exists(save_directory);

		// Generate a unique file path
		std::filesystem::path filePath = make_unique_filename(file_name, save_directory);

		std::ofstream out(filePath);
		if (!out.is_open())
		{
			std::cout << "Cannot save connectivity matrix " << filePath.filename() << std::endl;
			return;
		}

		for (int i = 0; i < matrix.size(); i++) { // iterate on properties (if any)
			for (int j = 0; j < matrix[i].size(); j++) {
				if (j < matrix[i].size() - 1) {
					out << matrix[i][j] << "\t"; // get name of property j
				}
				else {
					out << matrix[i][j];
				}
			}
			out << "\n";
		}
		out.flush();
		out.close();

	}

	void save_box(const std::string& file_name, const std::string& save_directory, Box box, std::vector<std::string> property_names, std::vector<std::vector<float>> properties)
	{

		ensure_directory_exists(save_directory);

		// Generate a unique file path
		std::filesystem::path filePath = make_unique_filename(file_name, save_directory);

		std::ofstream out(filePath);
		if (!out.is_open())
		{
			std::cout << "Cannot save skeleton to file (nodes): " << filePath.filename() << std::endl;
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
					if (properties.size()) {
						for (int p = 0; p < int(properties[0].size()); p++) { // iterate on properties (if any)
							out << "	" << std::fixed << std::setprecision(10) << properties[index][p];
						}
					}
					out << "\n";
				}
			}
		}

		out.flush();
		out.close();
	}
}