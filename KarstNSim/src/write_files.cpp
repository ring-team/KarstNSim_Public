/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/write_files.h"
#include <cerrno>
#include <cstring>

namespace KarstNSim {

	bool directory_exists(const std::string& directory) {
		struct stat info;
		if (stat(directory.c_str(), &info) != 0) {
			return false;
		}
		return (info.st_mode & S_IFDIR) != 0;
	}

	void create_directory(const std::string& directory) {
#ifdef _WIN32
		int res = _mkdir(directory.c_str());
#else
		int res = mkdir(directory.c_str(), 0777);
#endif
		if (res != 0) {
			std::cerr << "Error creating directory: " << directory << " (" << strerror(errno) << ")" << std::endl;
		}
	}

	// Function to check if a file exists
	bool file_exists(const std::string& file_path) {
		struct stat buffer;
		return (stat(file_path.c_str(), &buffer) == 0);
	}

	// Function to modify the file name by reference to ensure uniqueness
	void make_unique_filename(std::string& fileName, const std::string& saveDirectory) {
		std::string baseName = fileName;
		std::string extension = "";

		// Find the last dot to separate the file name and extension
		size_t dotPos = fileName.find_last_of('.');
		if (dotPos != std::string::npos) {
			baseName = fileName.substr(0, dotPos);
			extension = fileName.substr(dotPos);
		}

		int copyIndex = 0;

		// Keep modifying the file name until it's unique
		while (file_exists(saveDirectory + '/' + fileName)) {
			std::ostringstream oss;
			oss << baseName << '(' << copyIndex << ')' << extension;
			fileName = oss.str();
			copyIndex++;
		}
	}

	void make_unique_filename(std::string& file_name, const std::string& save_directory, Vector3 u, std::vector<std::string> property_names, std::vector<float> properties)
	{

		// Check if the save directory exists, if not, create it
		if (!directory_exists(save_directory)) {
			create_directory(save_directory);
		}

		// Ensure the file name is unique (modifies file_name by reference if needed)
		make_unique_filename(file_name, save_directory);

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
				out << "	" << std::fixed << std::setprecision(10) << properties[j];
			}
			out << "\n";
		}
		out.flush();
		out.close();
	}

	// Adjust the file name to the "last" existing version
	void adjust_file_name_to_last(std::string& file_name, const std::string& save_directory) {
		std::string base_name = file_name;
		std::string extension = "";

		// Separate the file name and extension
		size_t dot_pos = file_name.find_last_of('.');
		if (dot_pos != std::string::npos) {
			base_name = file_name.substr(0, dot_pos);
			extension = file_name.substr(dot_pos);
		}

		int copy_index = 0;
		std::string candidate_name = file_name;
		std::string candidate_name_new = file_name;
		std::ostringstream oss;
		std::string ossstr = candidate_name;
		oss << file_name;

		// Search for the last available file name in the sequence
		while (file_exists(save_directory + '/' + candidate_name_new)) {
			candidate_name = ossstr;
			std::ostringstream oss;
			oss << base_name << '(' << copy_index << ')' << extension;
			ossstr = oss.str();
			candidate_name_new = ossstr;
			copy_index++;
		}

		//// Adjust file_name to the "last existing" version
		//if (copy_index > 0) {
		//	std::ostringstream oss;
		//	oss << base_name << '(' << (copy_index - 1) << ')' << extension;
		//	
		//}
		file_name = candidate_name;
	}

	// Save the point data to a unique file
	void save_point(std::string& file_name, const std::string& save_directory, Vector3 u, std::vector<std::string> property_names, std::vector<float> properties) {

		// Check if the save directory exists, if not, create it
		if (!directory_exists(save_directory)) {
			create_directory(save_directory);
		}

		// Ensure the file name is unique (modifies file_name by reference)
		make_unique_filename(file_name, save_directory);

		std::ofstream out;
		out.open(save_directory + '/' + file_name);
		if (out.is_open() == false) {
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
				out << "	" << std::fixed << std::setprecision(10) << properties[j];
			}
			out << "\n";
		}
		out.flush();
		out.close();
	}

	void save_surface(std::string& file_name, const std::string& save_directory, Surface s, std::vector<std::string> property_names, std::vector<std::vector<float>> properties)
	{

		// Check if the save directory exists, if not, create it
		if (!directory_exists(save_directory)) {
			create_directory(save_directory);
		}

		// Ensure the file name is unique (modifies file_name by reference if needed)
		make_unique_filename(file_name, save_directory);

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

	void save_pointset(std::string& file_name, const std::string& save_directory, std::vector<Vector3> pset, std::vector<std::string> property_names, std::vector<std::vector<float>> properties)
	{

		// Check if the save directory exists, if not, create it
		if (!directory_exists(save_directory)) {
			create_directory(save_directory);
		}

		// Ensure the file name is unique (modifies file_name by reference if needed)
		make_unique_filename(file_name, save_directory);

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

	void save_line(std::string& file_name, const std::string& save_directory, Line pline, std::vector<std::string> property_names, std::vector<std::vector<std::vector<float>>> properties)
	{

		// Check if the save directory exists, if not, create it
		if (!directory_exists(save_directory)) {
			create_directory(save_directory);
		}

		// Ensure the file name is unique (modifies file_name by reference if needed)
		make_unique_filename(file_name, save_directory);

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

	void save_connectivity_matrix(std::string& file_name, const std::string& save_directory, Array2D<int> matrix) {

		// Check if the save directory exists, if not, create it
		if (!directory_exists(save_directory)) {
			create_directory(save_directory);
		}

		// Ensure the file name is unique (modifies file_name by reference if needed)
		make_unique_filename(file_name, save_directory);

		std::ofstream out;
		out.open(save_directory + '/' + file_name);
		if (out.is_open() == false)
		{
			std::cout << "Cannot save connectivity matrix " << file_name << std::endl;
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

	void save_box(std::string& file_name, const std::string& save_directory, Box box, std::vector<std::string> property_names, std::vector<std::vector<float>> properties)
	{

		// Check if the save directory exists, if not, create it
		if (!directory_exists(save_directory)) {
			create_directory(save_directory);
		}

		// Ensure the file name is unique (modifies file_name by reference if needed)
		make_unique_filename(file_name, save_directory);

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