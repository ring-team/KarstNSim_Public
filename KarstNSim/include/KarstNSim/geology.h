/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for classes KeyPointType, KeyPoint and CostTerm
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.
Find the corresponding code at https://github.com/aparis69/Karst-Synthesis

***************************************************************/

#pragma once

/*!
\file geology.h
\brief Declaration of various classes representing simple geological objects and parameters
\author Augustin GOUY
*/

#include "basics.h"
#include "graph.h"
#include <string>

namespace KarstNSim {
	/*!
	\class KeyPointType
	\brief Class that defines different keypoint types
	*/
	enum class KeyPointType
	{
		Sink,		// Entrance
		Spring,		// Exists
		Waypoint,	// waypoints are located inside the bedrock.
		Deadend,
		sampling,
	};

	/*!
	\class KeyPoint
	\brief Class that defines a keypoint, by its position and KeyPointType
	*/
	class KeyPoint
	{
	public:
		Vector3 p;
		KeyPointType type;
	public:
		inline explicit KeyPoint() { }
		inline explicit KeyPoint(const Vector3& p, KeyPointType t) : p(p), type(t)
		{
		}
	};

	/*!
	\class CostTerm
	\brief Class that defines costs terms, which are the different subcosts of the cost function. It defines if they are used or not and if yes what is their weight.
	*/
	class CostTerm
	{
	public:
		bool used;
		double weight;
	public:
		inline explicit CostTerm() : used(false), weight(1.0) { }
		inline explicit CostTerm(bool u, double w) : used(u), weight(w) { }
	};

	/*!
	\class GeologicalParameters
	\brief Class that defines all geological and simulation parameters, used throughout the code.
	*/
	class GeologicalParameters
	{
	public:
		std::string scenename;
		std::string directoryname;

		double graphPoissonRadius;
		double graphNeighbourRadius;
		double maxsize;
		bool graphuse_max_nghb_radius;
		int graphNeighbourCount;
		int nb_springs;
		bool multiply_costs;
		bool allow_single_outlet;
		std::vector<std::vector<Vector3>> PtsOldGraph;
		std::vector<std::vector<int>> IdxOldGraph;
		std::vector<double> z_list;
		std::vector<int> sinks_index;
		std::vector<int> wt_surf_index;

		std::vector<std::vector<int>> connectivity_matrix;

		std::vector<double> horizons;

		std::vector<double> fractures_orientations;
		std::vector<double> fractures_tolerances;
		std::vector<double> fractures_max_lengths;


		ScalarField2D heightfield;
		double elevationOffsetMin;
		double elevationOffsetMax;

		double max_dist_loops_vadose;
		double loop_density_vadose;
		double max_dist_loops_phreatic;
		double loop_density_phreatic;

		CostTerm distanceCost;
		CostTerm fractureCost;
		CostTerm horizonCost;
		CostTerm waterTable1;
		CostTerm waterTable2;
		CostTerm karstificationCost;
		CostTerm permeabilityCost;
		double gamma;

		double fracture_constraint_weight;
		std::vector<Sphere> spheres;
	};
}