/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for classes Triangle, Surface, Box, Segment, Line
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for classes Random, Ray, Plane, Sphere
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.

***************************************************************/

#pragma once

/*!
\file basics.h
\brief Declaration of various classes representing objects used for (hydro)geological model and karst network representation.
\author Augustin GOUY, Axel PARIS
*/

#include "KarstNSim/vec.h"
#include <time.h>
#include <vector>
#include "KarstNSim/randomgenerator.h"
#include <algorithm>
#include <memory>

namespace KarstNSim {


	/*!
	\class Ray
	\brief Class that defines a ray with an origin and direction.
	\param o Origin of the ray, represented as a Vector3.
	\param d Direction of the ray, represented as a Vector3 (should be normalized).
	*/
	class Ray
	{
	public:
		Vector3 o; //!< Origin of the ray.
		Vector3 d; //!< Direction of the ray.

		Ray(); //!< Default constructor for the Ray class.
		Ray(const Vector3& oo, const Vector3& dd); //!< Constructor that initializes the ray with a given origin and direction.
		Vector3 operator()(float t) const; //!< Computes a point along the ray at a given parameter t.
	};

	inline Ray::Ray()
	{
		// Default constructor implementation.
	}

	inline Ray::Ray(const Vector3& oo, const Vector3& dd)
	{
		o = oo;
		d = dd;
	}

	inline Vector3 Ray::operator()(float t) const
	{
		return o + t * d;
	}

	/*!
	\class Plane
	\brief Class that defines a plane in 3D space using a point and a normal vector.
	\param p A point on the plane, represented as a Vector3.
	\param n The normal vector of the plane, represented as a Vector3 (should be normalized).
	\param c The constant term in the plane equation (calculated internally).
	*/
	class Plane
	{
	protected:
		Vector3 p; //!< A point on the plane.
		Vector3 n; //!< The normal vector of the plane.
		float c;   //!< Constant term in the plane equation.

	public:
		Plane(); //!< Default constructor for the Plane class.
		Plane(const Vector3& p, const Vector3& n); //!< Constructor that initializes the plane with a point and normal.

		Vector3 Normal() const; //!< Returns the normal vector of the plane.
		Vector3 Point() const; //!< Returns a point on the plane.
		int Side(const Vector3& p) const; //!< Determines the relative position of a point with respect to the plane.
		float Signed(const Vector3& p) const; //!< Computes the signed distance from a point to the plane.
		static bool Intersection(const Plane& a, const Plane& b, const Plane& c, Vector3& p); //!< Finds the intersection of three planes (if it exists).
		static std::vector<Vector3> ConvexPoints(const std::vector<Plane>& planes); //!< Computes the intersection of a set of half-spaces.
	};

	/*!
	\brief Default constructor that initializes the plane to the origin with a zero normal vector.
	*/
	inline Plane::Plane()
	{
		p = Vector3(0);
		n = Vector3(0);
		c = 0.0;
	}

	/*!
	\brief Constructor from a point and a normal vector.
	\param pp A point on the plane.
	\param nn The normal vector of the plane (should be normalized).
	*/
	inline Plane::Plane(const Vector3& pp, const Vector3& nn)
	{
		p = pp;
		n = nn;
		c = Dot(p, n);
	}

	/*!
	\brief Returns the normal vector of the plane.
	\return A Vector3 representing the normal vector.
	*/
	inline Vector3 Plane::Normal() const
	{
		return n;
	}

	/*!
	\brief Computes the signed distance from a point to the plane.
	\param pp The point to compute the distance from.
	\return The signed distance as a float.
	*/
	inline float Plane::Signed(const Vector3& pp) const
	{
		return Dot(n, pp) - c;
	}

	/*!
	\brief Returns a point on the plane.
	\return A Vector3 representing a point on the plane.
	*/
	inline Vector3 Plane::Point() const
	{
		return p;
	}

	/*!
	\brief Determines the relative position of a point with respect to the plane.
	The epsilon tolerance used for testing is 10<SUP>-6</SUP>.
	\param pp The point to test.
	\return 1 if the point is above the plane, -1 if below, and 0 if on the plane.
	*/
	inline int Plane::Side(const Vector3& pp) const
	{
		float r = Dot(n, pp) - c;

		// Epsilon test
		if (r > 1e-6)
			return 1;
		else if (r < -1e-6)
			return -1;
		else
			return 0;
	}

	/*!
	\brief Finds the intersection point of three planes.
	\param a The first plane.
	\param b The second plane.
	\param c The third plane.
	\param p The intersection point (output parameter).
	\return True if the planes intersect at a point, false otherwise.
	*/
	inline bool Plane::Intersection(const Plane& a, const Plane& b, const Plane& c, Vector3& p)
	{
		float e = Matrix4(Matrix3(a.Normal(), b.Normal(), c.Normal())).Determinant();
		if (e < 1e-06)
			return false;
		p = (Dot(a.Point(), a.Normal()) * (Cross(b.Normal(), c.Normal()))) +
			(Dot(b.Point(), b.Normal()) * (Cross(c.Normal(), a.Normal()))) +
			(Dot(c.Point(), c.Normal()) * (Cross(a.Normal(), b.Normal())));
		p = p / (-e);
		return true;
	}

	/*!
	\brief Computes the intersection of a set of half-spaces.
	This is an O(n^4) algorithm.
	\param planes A set of planes defining the half-spaces.
	\return A vector of Vector3 points representing the minimal convex polygon embedding all half-spaces.
	*/
	inline std::vector<Vector3> Plane::ConvexPoints(const std::vector<Plane>& planes)
	{
		std::vector<Vector3> pts;
		for (int i = 0; i < int(planes.size()); i++)
		{
			for (int j = i + 1; j < int(planes.size()); j++)
			{
				for (int k = j + 1; k < int(planes.size()); k++)
				{
					Vector3 p;
					bool intersect = Intersection(planes[i], planes[j], planes[k], p);
					if (intersect)
					{
						bool isInside = true;
						for (int l = 0; l < int(planes.size()); l++)
						{
							// Do not check point ijk with one of its generating plane
							if (l == i || l == j || l == k)
								continue;
							int s = planes[l].Side(p);
							if (s > 0)
							{
								isInside = false;
								break;
							}
						}
						if (isInside)
							pts.push_back(p);
					}
				}
			}
		}
		return pts;
	}


	/*!
	\class Sphere
	\brief Class that represents a spherical geometric element with a center and radius.
	\param c The center of the sphere, represented as a Vector3.
	\param r The radius of the sphere, represented as a float.
	*/
	class Sphere
	{
	public:
		Vector3 center; //!< The center of the sphere.
		float radius;   //!< The radius of the sphere.

	public:
		Sphere(const Vector3& c, float r); //!< Constructor that initializes the sphere with a center and radius.

		float Distance(const Vector3& p) const; //!< Computes the distance between a point and the surface of the sphere.
		bool Contains(const Vector3& p) const; //!< Checks if a point lies inside the sphere.
		Vector3 RandomSurface() const; //!< Returns a random point on the surface of the sphere.
		Vector3 Center() const; //!< Returns the center of the sphere.
		float Radius() const; //!< Returns the radius of the sphere.
	};

	/*!
	\brief Constructor that initializes a sphere with a center and radius.
	\param c The center of the sphere.
	\param r The radius of the sphere.
	*/
	inline Sphere::Sphere(const Vector3& c, float r) : center(c), radius(r)
	{

	}

	/*!
	\brief Computes the distance between a point and the surface of the sphere.
	If the point lies inside the sphere, the distance is considered to be 0.
	\param p The world point to compute the distance from.
	\return The distance from the point to the sphere's surface (0 if inside the sphere).
	*/
	inline float Sphere::Distance(const Vector3& p) const
	{
		float a = Dot((p - center), (p - center));
		if (a < radius * radius)
		{
			return 0.0;
		}
		a = sqrt(a) - radius;
		a *= a;
		return a;
	}

	/*!
	\brief Checks if a point lies inside the sphere (i.e., if Distance(p) == 0).
	\param p The world point to check.
	\return True if the point lies inside the sphere, false otherwise.
	*/
	inline bool Sphere::Contains(const Vector3& p) const
	{
		return Distance(p) == 0;
	}

	/*!
	\brief Returns a random point on the surface of the sphere.
	\return A random point on the surface of the sphere as a Vector3.
	*/
	inline Vector3 Sphere::RandomSurface() const
	{
		return Vector3(0);
	}

	/*!
	\brief Returns the center of the sphere.
	\return The center of the sphere as a Vector3.
	*/
	inline Vector3 Sphere::Center() const
	{
		return center;
	}

	/*!
	\brief Returns the radius of the sphere.
	\return The radius of the sphere as a float.
	*/
	inline float Sphere::Radius() const
	{
		return radius;
	}

	/*!
	\class Triangle
	\brief Class that represents a 3D triangle defined by three points.
	\param a The index of the first vertex (point) of the triangle.
	\param b The index of the second vertex (point) of the triangle.
	\param c The index of the third vertex (point) of the triangle.
	*/
	class Triangle
	{
	private:
		int pts[3]; //!< The array of three points (vertices) that define the triangle.

	public:
		Triangle(); //!< Default constructor for the Triangle class.
		Triangle(const int& a, const int& b, const int& c); //!< Constructor that initializes the triangle with three point indices.
		const int& point(const int& i) const; //!< Returns the index of one of the points of the triangle.
		bool containsVertex(const int& vertexIndex) const; //!< Checks if a given vertex index is part of the triangle.
	};

	/*!
	\brief Default constructor that initializes an empty triangle.
	*/
	inline Triangle::Triangle()
	{
	}

	/*!
	\brief Constructor that initializes the triangle with three given points.
	\param a The index of the first point.
	\param b The index of the second point.
	\param c The index of the third point.
	*/
	inline Triangle::Triangle(const int& a, const int& b, const int& c)
	{
		pts[0] = a;
		pts[1] = b;
		pts[2] = c;
	}

	/*!
	\brief Returns one of the points of the triangle by index.
	\param i The index of the point (0, 1, or 2).
	\return The index of the point as an integer.
	*/
	inline const int& Triangle::point(const int& i) const
	{
		return pts[i];
	}

	/*!
	\brief Checks if a given vertex index is part of the triangle.
	\param vertexIndex The index of the vertex to check.
	\return True if the vertex index is part of the triangle, false otherwise.
	*/
	inline bool Triangle::containsVertex(const int& vertexIndex) const {
		return (vertexIndex == pts[0] || vertexIndex == pts[1] || vertexIndex == pts[2]);
	}


	/*!
	\class Surface
	\brief Class that represents a 3D surface defined by a collection of points and triangles.
	\param pts A vector of points that define the vertices of the triangles on the surface.
	\param trgls A vector of triangles that define the surface, each triangle consisting of three point indices.
	\param name A string name for the surface (default is "DefaultName").
	*/
	class Surface
	{
	private:
		std::vector<Vector3> pts_; //!< The list of points (vertices) defining the surface.
		std::vector<Triangle> trgls_; //!< The list of triangles that make up the surface.
		Vector3 b_min_; //!< The minimum bound of the surface (bounding box).
		Vector3 b_max_; //!< The maximum bound of the surface (bounding box).
		std::vector<float> circumradii_; //!< A vector of circumradii, one per triangle.
		std::vector<Vector3> centers_; //!< A vector of triangle centers.
		std::string name_; //!< The name of the surface.

	public:
		Surface(const std::string& name = "DefaultName"); //!< Default constructor with optional surface name.
		Surface(std::vector<Vector3>, std::vector<Triangle>, const std::string& name = "DefaultName"); //!< Constructor that initializes the surface with points, triangles, and an optional name.
		int get_nb_trgls(); //!< Returns the number of triangles in the surface.
		int get_nb_pts(); //!< Returns the number of points in the surface.
		bool is_empty() const; //!< Checks if the surface is empty (no points or triangles).
		const Triangle& get_triangle(const int& i); //!< Returns the triangle at index i.
		const Triangle& get_triangle(const int& i) const; //!< Returns the triangle at index i (const version).
		const Vector3& get_node(const int& i) const; //!< Returns the point (node) at index i.
		Vector3 nearest(Vector3 p1, Vector3 p2, Vector3 ptest); //!< Finds the nearest point between two points p1 and p2 to the test point ptest.
		Vector3 nearest(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 ptest); //!< Finds the nearest point among three points p1, p2, and p3 to the test point ptest.
		Vector3 get_boundbox_min(); //!< Returns the minimum bound of the surface (bounding box).
		Vector3 get_boundbox_max(); //!< Returns the maximum bound of the surface (bounding box).
		float get_circumradii(const int& i); //!< Returns the circumradius of the triangle at index i.
		Vector3 get_trgl_center(const int& trgl_idx); //!< Returns the center of the triangle at index trgl_idx.
		PointCloud get_centers_cloud(int dim) const; //!< Returns a PointCloud of the centers of the triangles in the surface.
		int get_nb_valid_trgls(); //!< Returns the number of valid triangles (those with a non-zero circumradius).
	};

	/*!
	\brief Default constructor that initializes the surface with a given name.
	\param name The name of the surface (default is "DefaultName").
	*/
	inline Surface::Surface(const std::string& name) : name_(name) {}

	/*!
	\brief Constructor that initializes the surface with points and triangles.
	\param pts The vector of points that define the surface.
	\param trgls The vector of triangles that define the surface.
	\param name The name of the surface (default is "DefaultName").
	*/
	inline Surface::Surface(std::vector<Vector3> pts, std::vector<Triangle> trgls, const std::string& name)
		: name_(name) {
		std::vector<Vector3> centers;
		pts_ = pts;
		// phase of trimming to remove bad triangles (3 points aligned, or worse...)
		// another parameter that is practical for TINs is the length of the largest circumradius among all its triangles (in map view)
		// We compute it directly when creating the surface.
		//float eps = 1e-5;
		//circumradii_.resize(trgls.size(), 0.);
		for (int trgl = 0; trgl < trgls.size(); trgl++) {
			Triangle trgl_i = trgls[trgl];
			int pt1 = trgl_i.point(0), pt2 = trgl_i.point(1), pt3 = trgl_i.point(2);
			Vector3 p1 = pts[pt1], p2 = pts[pt2], p3 = pts[pt3];
			// 2D dist between p1 and p2
			float dist1 = magnitude2D(p1 - p2);
			// 2D dist between p2 and p3
			float dist2 = magnitude2D(p2 - p3);
			// 2D dist between p1 and p3
			float dist3 = magnitude2D(p1 - p3);
			Vector3 center;
			center.x = (p1.x + p2.x + p3.x) / 3.0;
			center.y = (p1.y + p2.y + p3.y) / 3.0;
			center.z = (p1.z + p2.z + p3.z) / 3.0;
			// compute circumradius (proxy of the max distance Fermat-Torricelli point in the triangle)
			// Simple meaning = if you pick any random point in the domain, if it is less than a circumradius away in aerial view
			// from any vertex of the TIN, then it is very probably inside the TIN, or at least very close from it.
			// From my research it is the quickest way to make this test while keeping (very) good precision.
			if (isTriangleValid(dist1, dist2, dist3)) {
				float circum = calculateCircumradius(dist1, dist2, dist3);
				//if ((circum - dist1) < eps && (circum - dist2) < eps && (circum - dist3) < eps) {
				circumradii_.push_back(circum);
				trgls_.push_back(trgl_i);
				centers.push_back(center);
				//}
			}
		}

		centers_ = centers;

		float x_min = pts[0].x;
		float x_max = pts[0].x;
		float y_min = pts[0].y;
		float y_max = pts[0].y;
		float z_min = pts[0].z;
		float z_max = pts[0].z;
		for (int i = 0; i < int(pts.size()); i++) {
			if (pts[i].x < x_min) x_min = pts[i].x;
			if (pts[i].x > x_max) x_max = pts[i].x;
			if (pts[i].y < y_min) y_min = pts[i].y;
			if (pts[i].y > y_max) y_max = pts[i].y;
			if (pts[i].z < z_min) z_min = pts[i].z;
			if (pts[i].z > z_max) z_max = pts[i].z;
		}
		b_min_ = Vector3(x_min, y_min, z_min); // this is the coordinates of the bottom corner in all 3 directions
		b_max_ = Vector3(x_max, y_max, z_max); // this is the coordinates of the top corner in all 3 directions

	}


	/*!
	\brief Returns a PointCloud of the centers of the triangles in the surface.
	\param dim The dimensionality of the points (default is 3).
	\return A PointCloud of the triangle centers.
	*/
	inline PointCloud Surface::get_centers_cloud(int dim = 3) const {
		std::vector<Vector3> centers;
		for (int trgl_idx = 0; trgl_idx < trgls_.size(); ++trgl_idx) {
			int pt1 = trgls_[trgl_idx].point(0), pt2 = trgls_[trgl_idx].point(1), pt3 = trgls_[trgl_idx].point(2);
			Vector3 p1 = get_node(pt1), p2 = get_node(pt2), p3 = get_node(pt3);
			Vector3 center(0., 0., 0.);
			center.x = (p1.x + p2.x + p3.x) / 3.0;
			center.y = (p1.y + p2.y + p3.y) / 3.0;
			if (dim == 3) {
				center.z = (p1.z + p2.z + p3.z) / 3.0;
			}
			centers.push_back(center);
		}
		return PointCloud(centers);
	}

	/*!
	\brief Returns the number of triangles in the surface.
	\return The number of triangles.
	*/
	inline int Surface::get_nb_trgls() {
		return int(trgls_.size());
	}

	/*!
	\brief Returns the number of points in the surface.
	\return The number of points.
	*/
	inline int Surface::get_nb_pts() {
		return int(pts_.size());
	}

	/*!
	\brief Checks if the surface is empty (no points or triangles).
	\return True if the surface is empty, false otherwise.
	*/
	inline bool Surface::is_empty() const {
		return pts_.size() == 0 && trgls_.size() == 0;
	}

	/*!
	\brief Returns the triangle at index i.
	\param i The index of the triangle.
	\return A reference to the triangle at index i.
	*/
	inline const Triangle& Surface::get_triangle(const int& i) {
		return trgls_[i];
	}

	/*!
	\brief Returns the triangle at index i (const version).
	\param i The index of the triangle.
	\return A constant reference to the triangle at index i.
	*/
	inline const Triangle& Surface::get_triangle(const int& i) const {
		return trgls_[i];
	}

	/*!
	\brief Returns the point (node) at index i.
	\param i The index of the point.
	\return A constant reference to the point at index i.
	*/
	inline const Vector3& Surface::get_node(const int& i) const {
		return pts_[i];
	}

	/*!
	\brief Finds the nearest point between two points p1 and p2 to the test point ptest.
	\param p1 The first point.
	\param p2 The second point.
	\param ptest The test point.
	\return The nearest point to ptest.
	*/
	inline Vector3 Surface::nearest(Vector3 p1, Vector3 p2, Vector3 ptest) {
		float dist = squaredmagnitude(p1 - ptest);
		float dist2 = squaredmagnitude(p2 - ptest);
		return (dist <= dist2) ? p1 : p2;
	}

	/*!
	\brief Returns the number of valid triangles (those with a non-zero circumradius).
	\return The number of valid triangles.
	*/
	inline int Surface::get_nb_valid_trgls() {
		return int(circumradii_.size());
	}

	/*!
	\brief Finds the nearest point among three points p1, p2, and p3 to the test point ptest.
	\param p1 The first point.
	\param p2 The second point.
	\param p3 The third point.
	\param ptest The test point.
	\return The nearest point to ptest.
	*/
	inline Vector3 Surface::nearest(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 ptest) {

		float dist = squaredmagnitude(p1 - ptest);
		float dist2 = squaredmagnitude(p2 - ptest);
		float dist3 = squaredmagnitude(p3 - ptest);
		if (dist < dist2) {
			if (dist < dist3) {
				return p1;
			}
			return p3;
		}
		else if (dist2 < dist3) {
			return p2;
		}
		return p3;
	}

	/*!
	\brief Returns the center of the triangle at index trgl_idx.
	\param trgl_idx The index of the triangle.
	\return The center of the triangle at the specified index.
	*/
	inline Vector3 Surface::get_trgl_center(const int& trgl_idx) {

		return centers_[trgl_idx];
	}

	/*!
	\brief Returns the minimum bound of the surface (bounding box).
	\return The minimum bound vector.
	*/
	inline Vector3 Surface::get_boundbox_min() {
		return b_min_;
	}

	/*!
	\brief Returns the maximum bound of the surface (bounding box).
	\return The maximum bound vector.
	*/
	inline Vector3 Surface::get_boundbox_max() {
		return b_max_;
	}

	/*!
	\brief Returns the circumradius of the triangle at index i.
	\param i The index of the triangle.
	\return The circumradius of the triangle at index i.
	*/
	inline float Surface::get_circumradii(const int& i) {
		return circumradii_[i];
	}

	/*!
	\class Segment
	\brief Class that represents a line segment defined by two points in 3D space.
	\param p1 The first point (start point) of the segment.
	\param p2 The second point (end point) of the segment.
	*/
	class Segment
	{
	protected:
		Vector3 p1_; //!< The start point of the segment.
		Vector3 p2_; //!< The end point of the segment.
	public:
		Segment(); //!< Default constructor.
		Segment(const Vector3& a, const Vector3& b); //!< Constructor that initializes the segment with two points.
		Vector3 start(); //!< Returns the start point of the segment.
		Vector3 end(); //!< Returns the end point of the segment.
	};

	/*!
	\brief Default constructor that initializes the segment with default values.
	*/
	inline Segment::Segment()
	{
	}

	/*!
	\brief Constructor that initializes the segment with two points.
	\param a The start point of the segment.
	\param b The end point of the segment.
	*/
	inline Segment::Segment(const Vector3& a, const Vector3& b)
	{
		p1_ = a;
		p2_ = b;
	}

	/*!
	\brief Returns the start point of the segment.
	\return The start point of the segment (Vector3).
	*/
	inline Vector3 Segment::start()
	{
		return p1_;
	}

	/*!
	\brief Returns the end point of the segment.
	\return The end point of the segment (Vector3).
	*/
	inline Vector3 Segment::end()
	{
		return p2_;
	}

	/*!
	\class Line
	\brief A class representing a polyline made of connected segments.
	\param segs A vector of segments that form the polyline.
	*/
	class Line
	{
	protected:
		std::vector<Segment> segs_; //!< A vector of segments forming the line.
		int nb_unique_nodes_; //!< The number of unique nodes (points) in the line.
		std::vector<Vector3> unique_nodes_; //!< A vector of unique nodes (points).
	public:
		Line(); //!< Default constructor.
		Line(std::vector<Segment> segs); //!< Constructor that initializes the line with a vector of segments.
		void append(Segment seg); //!< Appends a segment to the line.

		int size(); //!< Returns the number of segments in the line.
		int size() const; //!< Const version of the size function, returning the number of segments.

		Segment get_seg(int i) const; //!< Returns the i-th segment (const version).
		Segment& get_seg(int i); //!< Returns the i-th segment (non-const version).

		std::vector<Vector3> get_unique_nodes() const; //!< Returns the list of unique nodes (points) of the line (const version).
		std::vector<Vector3>& get_unique_nodes(); //!< Returns the list of unique nodes (non-const version).

		int get_nb_unique_nodes() const; //!< Returns the number of unique nodes (const version).
		int get_nb_unique_nodes(); //!< Returns the number of unique nodes (non-const version).

		int get_nb_segs() const; //!< Returns the number of segments (const version).
		int get_nb_segs(); //!< Returns the number of segments (non-const version).

		Vector2 get_bbox_min() const; //!< Returns the bounding box minimum (2D) for the line.
		Vector2 get_bbox_max() const; //!< Returns the bounding box maximum (2D) for the line.
	};

	/*!
	\brief Default constructor that initializes the line with no segments and no nodes.
	*/
	inline Line::Line() {}

	/*!
	\brief Constructor that initializes the line with a vector of segments.
	\param segs The vector of segments to initialize the line.
	*/
	inline Line::Line(std::vector<Segment> segs) {
		segs_ = segs;

		for (int i = 0; i<int(segs_.size()); i++) {

			//int idx = std::find(unique_nodes_.begin(), unique_nodes_.end(), segs_[i].start()) - unique_nodes_.end();

			if (std::find(unique_nodes_.begin(), unique_nodes_.end(), segs_[i].start()) == unique_nodes_.end()) { // node is unique
				unique_nodes_.push_back(segs_[i].start());
			}
			if (std::find(unique_nodes_.begin(), unique_nodes_.end(), segs_[i].end()) == unique_nodes_.end()) { // node is unique
				unique_nodes_.push_back(segs_[i].end());
			}
		}
		nb_unique_nodes_ = int(unique_nodes_.size());
	}

	/*!
	\brief Appends a segment to the line.
	\param seg The segment to be appended.
	*/
	inline void Line::append(Segment seg) {
		segs_.push_back(seg);
	}

	/*!
	\brief Returns the number of segments in the line.
	\return The number of segments.
	*/
	inline int Line::size() {
		return int(segs_.size());
	}

	/*!
	\brief Returns the number of segments in the line (const version).
	\return The number of segments.
	*/
	inline int Line::size() const {
		return int(segs_.size());
	}

	/*!
	\brief Returns the i-th segment in the line (const version).
	\param i The index of the segment to retrieve.
	\return The i-th segment.
	*/
	inline Segment Line::get_seg(int i) const {
		return segs_[i];
	}

	/*!
	\brief Returns the i-th segment in the line (non-const version).
	\param i The index of the segment to retrieve.
	\return The i-th segment.
	*/
	inline Segment& Line::get_seg(int i) {
		return segs_[i];
	}

	/*!
	\brief Returns the list of unique nodes (points) of the line (const version).
	\return A vector of unique nodes.
	*/
	inline std::vector<Vector3> Line::get_unique_nodes() const {
		return unique_nodes_;
	}

	/*!
	\brief Returns the list of unique nodes (points) of the line (non-const version).
	\return A vector of unique nodes.
	*/
	inline std::vector<Vector3>& Line::get_unique_nodes() {
		return unique_nodes_;
	}

	/*!
	\brief Returns the number of unique nodes (const version).
	\return The number of unique nodes.
	*/
	inline int Line::get_nb_unique_nodes() const {
		return nb_unique_nodes_;
	}

	/*!
	\brief Returns the number of unique nodes (non-const version).
	\return The number of unique nodes.
	*/
	inline int Line::get_nb_unique_nodes() {
		return nb_unique_nodes_;
	}

	/*!
	\brief Returns the number of segments in the line (const version).
	\return The number of segments.
	*/
	inline int Line::get_nb_segs() const {
		return int(segs_.size());
	}

	/*!
	\brief Returns the number of segments in the line (non-const version).
	\return The number of segments.
	*/
	inline int Line::get_nb_segs() {
		return int(segs_.size());
	}

	/*!
	\brief Returns the minimum (bounding box) 2D point for the line.
	\return The minimum point of the bounding box (Vector2).
	*/
	inline Vector2 Line::get_bbox_min() const {
		float min_x = std::numeric_limits<float>::max();
		float min_y = std::numeric_limits<float>::max();

		for (const auto& node : unique_nodes_) {
			if (node.x < min_x) min_x = node.x;
			if (node.y < min_y) min_y = node.y;
		}

		return Vector2(min_x, min_y);
	}

	/*!
	\brief Returns the maximum (bounding box) 2D point for the line.
	\return The maximum point of the bounding box (Vector2).
	*/
	inline Vector2 Line::get_bbox_max() const {
		float max_x = -std::numeric_limits<float>::max();
		float max_y = -std::numeric_limits<float>::max();

		for (const auto& node : unique_nodes_) {
			if (node.x > max_x) max_x = node.x;
			if (node.y > max_y) max_y = node.y;
		}

		return Vector2(max_x, max_y);
	}

	/**
	 * \brief Represents a 3D Axis-Aligned Bounding Box (AABB).
	 *
	 * This class provides methods to calculate geometric properties such as the center, diagonal,
	 * and distance to a point, as well as methods for transforming points between coordinate systems
	 * and generating random points inside the box.
	 */
	class Box
	{
	protected:
		Vector3 basis; ///< The lower-left corner of the box in world coordinates.
		Vector3 end; ///< The upper-right corner of the box in world coordinates.
		Vector3 u, v, w; ///< The step vectors of the box in the u, v, w directions.
		int nu, nv, nw; ///< Number of steps in each direction (u, v, w).
		std::vector<float> mat; ///< 3x3 matrix representing the box's coordinate system.
		std::vector<float> mat_t; ///< The inverse of the transformation matrix (useful for transforming coordinates).

	public:
		/**
		  * \brief Default constructor for a box.
		  */
		inline Box() { }
		/**
		 * \brief Constructs a box with given basis, step vectors, and grid dimensions.
		 * \param basis The lower-left corner of the box in world coordinates.
		 * \param u Step vector in the u-direction.
		 * \param v Step vector in the v-direction.
		 * \param w Step vector in the w-direction.
		 * \param nu Number of steps in the u-direction.
		 * \param nv Number of steps in the v-direction.
		 * \param nw Number of steps in the w-direction.
		 */
		explicit Box(const Vector3& basis, const Vector3& u, const Vector3& v, const Vector3& w, const int& nu, const int& nv, const int& nw);

		/**
	 * \brief Checks if the point `p` is inside the box.
	 * \param p The point in world coordinates.
	 * \return True if the point is inside the box, false otherwise.
	 */
		bool contains(const Vector3&) const;

		/**
	 * \brief Computes and returns the center of the box.
	 * \return The center of the box as a `Vector3` object.
	 */
		Vector3 Center() const;

		/**
	 * \brief Computes the distance between a point and the box.
	 * \param p The point in world coordinates.
	 * \return The distance from the point to the box.
	 */
		float Distance(const Vector3& p) const;

		/**
	 * \brief Computes and returns the diagonal (size) of the box.
	 * \return The diagonal of the box as a `Vector3` object.
	 */
		Vector3 Diagonal() const;

		/**
	 * \brief Generates a random point inside the box.
	 * \return A random point inside the box.
	 */
		Vector3 RandomInside() const;

		/**
	 * \brief Returns one of the vertices of the box.
	 * \param i The index of the vertex (0 for the `basis` and 1 for the `end`).
	 * \return The vertex of the box corresponding to the index `i`.
	 */
		Vector3 Vertex(int) const;

		/**
		 * \brief Returns the lower-left corner of the box in world coordinates.
		 * \return The `basis` of the box.
		 */
		Vector3 get_basis() const;

		/**
		 * \brief Returns the upper-right corner of the box in world coordinates.
		 * \return The `end` of the box.
		 */
		Vector3 get_end() const;

		/**
		 * \brief Returns the u-direction step vector.
		 * \return The `u` vector of the box.
		 */
		Vector3 get_u() const;

		/**
		 * \brief Returns the v-direction step vector.
		 * \return The `v` vector of the box.
		 */
		Vector3 get_v() const;

		/**
		 * \brief Returns the w-direction step vector.
		 * \return The `w` vector of the box.
		 */
		Vector3 get_w() const;

		/**
		 * \brief Returns the number of steps in the u-direction.
		 * \return The number of steps in the u-direction (`nu`).
		 */
		int get_nu() const;

		/**
		 * \brief Returns the number of steps in the v-direction.
		 * \return The number of steps in the v-direction (`nv`).
		 */
		int get_nv() const;

		/**
		 * \brief Returns the number of steps in the w-direction.
		 * \return The number of steps in the w-direction (`nw`).
		 */
		int get_nw() const;

		/**
		 * \brief Transforms the point from the uvw coordinate system to xyz world coordinates.
		 * \param up The u-coordinate.
		 * \param vp The v-coordinate.
		 * \param wp The w-coordinate.
		 * \param cellcentered Flag indicating whether the point corresponds to the center of a cell (default is false).
		 * \return The transformed point in world coordinates.
		 */
		Vector3 uvw2xyz(int up, int vp, int wp, bool cellcentered = false) const;

		/**
		 * \brief Transforms a point from world coordinates to the uvw coordinate system.
		 * \param pt The point in world coordinates.
		 * \param u The transformed u-coordinate.
		 * \param v The transformed v-coordinate.
		 * \param w The transformed w-coordinate.
		 */
		void xyz2uvw(const Vector3& pt, int &u, int &v, int &w);

		/**
		 * \brief Transforms a point from world coordinates to the uvw coordinate system, with bounds checking.
		 * \param pt The point in world coordinates.
		 * \param u1 The transformed u-coordinate.
		 * \param v1 The transformed v-coordinate.
		 * \param w1 The transformed w-coordinate.
		 */
		void xyz2uvw_with_limits_conditions(const Vector3& pt, int &u1, int &v1, int &w1);

		/**
		 * \brief Converts the cell coordinates (u, v, w) to a 1D index.
		 * \param u The u-coordinate.
		 * \param v The v-coordinate.
		 * \param w The w-coordinate.
		 * \param idx The computed index.
		 */
		void ravel(int u, int v, int w, int& idx) const;

		/**
		 * \brief Converts a 1D index to the corresponding (u, v, w) cell coordinates.
		 * \param u The u-coordinate.
		 * \param v The v-coordinate.
		 * \param w The w-coordinate.
		 * \param idx The index to convert.
		 */
		void unravel(int& u, int& v, int& w, int idx) const;

		/**
		 * \brief Computes a Poisson disk sampling distribution inside the box.
		 * \param p The existing sample points (can be empty).
		 * \param r The radius of the spheres.
		 * \param n The number of candidate points to generate.
		 */
		void Poisson(std::vector<Vector3>& samples, float r, int n) const;

		/**
		 * \brief Computes a Poisson disk sampling distribution inside the box, avoiding certain spheres.
		 * \param p The existing sample points (can be empty).
		 * \param r The radius of the spheres.
		 * \param n The number of candidate points to generate.
		 * \param spheres The list of spheres where points should not be generated.
		 */
		void Poisson_adapted(std::vector<Vector3>& p, float r, int n, std::vector<Sphere>& spheres) const;

		Vector3& operator[](int i);
		Vector3 operator[](int i) const;

	};

	inline Box::Box(const Vector3& basis, const Vector3& u, const Vector3& v, const Vector3& w, const int& nu, const int& nv, const int& nw) : basis(basis), u(u), v(v), w(w), nu(nu), nv(nv), nw(nw)
	{
		end = basis + u + v + w;
		mat = { u.x,v.x,w.x,u.y,v.y,w.y,u.z,v.z,w.z };
		float a = mat[0], b = mat[1], c = mat[2], d = mat[3], e = mat[4], f = mat[5], g = mat[6], h = mat[7], i = mat[8];
		float det = a * (e*i - f * h) - b * (d*i - f * g) + c * (d*h - e * g);
		mat_t.push_back(e*i - f * h);
		mat_t.push_back(-(b*i - c * h));
		mat_t.push_back(b*f - c * e);
		mat_t.push_back(-(d*i - f * g));
		mat_t.push_back(a*i - c * g);
		mat_t.push_back(-(a*f - c * d));
		mat_t.push_back(d*h - e * g);
		mat_t.push_back(-(a*h - b * g));
		mat_t.push_back(a*e - b * d);
		for (int iter = 0; iter < int(mat_t.size()); iter++) {
			mat_t[iter] /= det;
		}
	}

	inline Vector3 Box::RandomInside() const
	{
		Vector3 s = end - basis;
		float randw = generateRandomFloat(-1.0f * s[0] / 2.0f, s[0] / 2.0f);
		float randh = generateRandomFloat(-1.0f * s[1] / 2.0f, s[1] / 2.0f);
		float randl = generateRandomFloat(-1.0f * s[2] / 2.0f, s[2] / 2.0f);
		return (end + basis) / 2.0f + Vector3(randw, randh, randl);
	}

	/*
	\brief Returns true if p is inside the box, false otherwise.
	\param p world point
	*/
	//inline bool Box::contains(const Vector3& p) const
	//{
	//	// we don't know for sure if basis x, y and z are smaller than end x, y and z so we operate this simple test which works in either case:
	//	float dx = (p.x - basis.x)*(p.x - end.x);
	//	float dy = (p.y - basis.y)*(p.y - end.y);
	//	float dz = (p.z - basis.z)*(p.z - end.z);
	//	return ((dx<=0)&&(dy<=0)&&(dz<=0));
	//}
	inline bool Box::contains(const Vector3& p) const {
		// Transform the point to the box's local coordinate system
		float px = p.x - basis.x;
		float py = p.y - basis.y;
		float pz = p.z - basis.z;

		// Apply the inverse transformation matrix mat_t to the point (px, py, pz)
		float local_x = mat_t[0] * px + mat_t[1] * py + mat_t[2] * pz;
		float local_y = mat_t[3] * px + mat_t[4] * py + mat_t[5] * pz;
		float local_z = mat_t[6] * px + mat_t[7] * py + mat_t[8] * pz;

		// Check if the transformed point is within the bounds of the box
		return (local_x >= 0 && local_x <= 1 &&
			local_y >= 0 && local_y <= 1 &&
			local_z >= 0 && local_z <= 1);
	}

	/*!
	\brief Computes and returns the center of the box.
	*/
	inline Vector3 Box::Center() const
	{
		return (basis + end) / 2.0f;
	}

	/*!
	\brief Returns the diagonal/size of the box.
	*/
	inline Vector3 Box::Diagonal() const
	{
		return (end - basis);
	}

	/*!
	DOESNT WORK AS IT IS NOW
	\brief Compute the distance between a point and the box.
	\param p point
	*/
	inline float Box::Distance(const Vector3& p) const
	{
		float r = 0.0;
		for (int i = 0; i < 3; i++)
		{
			if (p[i] < basis[i])
			{
				float s = p[i] - basis[i];
				r += s * s;
			}
			else if (p[i] > end[i])
			{
				float s = p[i] - end[i];
				r += s * s;
			}
		}
		return r;
	}


	/*
	\brief Get one of the vertex of the box.
	*/
	inline Vector3 Box::Vertex(int i) const
	{
		if (i == 0)
			return basis;
		return end;
	}

	/*
	\brief Get bottom left vertex in world coordinates
	*/
	inline Vector3 Box::get_basis() const
	{
		return basis;
	}

	/*
	\brief Get top right vertex in world coordinates
	*/
	inline Vector3 Box::get_end() const
	{
		return end;
	}

	inline Vector3 Box::get_u() const
	{
		return u;
	}

	inline Vector3 Box::get_v() const
	{
		return v;
	}

	inline Vector3 Box::get_w() const
	{
		return w;
	}

	inline int Box::get_nu() const
	{
		return nu;
	}

	inline int Box::get_nv() const
	{
		return nv;
	}

	inline int Box::get_nw() const
	{
		return nw;
	}

	/*
	\brief Access box vertex by reference
	*/
	inline Vector3& Box::operator[](int i)
	{
		if (i == 0)
			return basis;
		return end;
	}

	/*
	\brief Access box vertex by const value
	*/
	inline Vector3 Box::operator[](int i) const
	{
		if (i == 0)
			return basis;
		return end;
	}

	inline void Box::xyz2uvw(const Vector3& pt, int &up, int &vp, int &wp) {

		// 1: normalize vector

		Vector3 ptn;
		ptn.x = mat_t[0] * (pt.x - basis.x) + mat_t[1] * (pt.y - basis.y) + mat_t[2] * (pt.z - basis.z);
		ptn.y = mat_t[3] * (pt.x - basis.x) + mat_t[4] * (pt.y - basis.y) + mat_t[5] * (pt.z - basis.z);
		ptn.z = mat_t[6] * (pt.x - basis.x) + mat_t[7] * (pt.y - basis.y) + mat_t[8] * (pt.z - basis.z);

		// 2: get u/v/w coordinates in normalized domain

		up = static_cast<int>(std::floor(ptn.x * nu));
		vp = static_cast<int>(std::floor(ptn.y * nv));
		wp = static_cast<int>(std::floor(ptn.z * nw));

		// Ensure that the indices are within bounds
		up = std::max(0, std::min(up, nu - 1));
		vp = std::max(0, std::min(vp, nv - 1));
		wp = std::max(0, std::min(wp, nw - 1));
	}

	inline Vector3 Box::uvw2xyz(int up, int vp, int wp, bool cellcentered) const {

		// if cell centered, the returned point is the center of the corresponding cell.
		// if not, the grid is corner centered, and then the returned value is the cell corner (the one closer to the origin)
		float offset = cellcentered ? 0.5 : 0.0;

		// Calculate the position of the cell corner or center in the normalized domain
		float u_normalized = (static_cast<float>(up) + offset) / nu;
		float v_normalized = (static_cast<float>(vp) + offset) / nv;
		float w_normalized = (static_cast<float>(wp) + offset) / nw;

		// Transform back to the original coordinate system
		Vector3 point;
		point.x = u_normalized * u.x + v_normalized * v.x + w_normalized * w.x + basis.x;
		point.y = u_normalized * u.y + v_normalized * v.y + w_normalized * w.y + basis.y;
		point.z = u_normalized * u.z + v_normalized * v.z + w_normalized * w.z + basis.z;

		return point;
	}

	inline void Box::xyz2uvw_with_limits_conditions(const Vector3& pt, int &u1, int &v1, int &w1) {
		xyz2uvw(pt, u1, v1, w1);
		if (u1 < 0) {
			u1 = 0;
		}
		if (u1 >= nu) {
			u1 = nu - 1;
		}
		if (v1 < 0) {
			v1 = 0;
		}
		if (v1 >= nv) {
			v1 = nv - 1;
		}
		if (w1 < 0) {
			w1 = 0;
		}
		if (w1 >= nw) {
			w1 = nw - 1;
		}
	}

	/*!
	\brief Compute a Poisson sphere distribution inside a box.
	This function uses a simple O(n<SUP>3</SUP>) dart throwing algorithm.
	\param array of existing samples (possibly empty)
	\param r Radius of the sphere.
	\param n Number of candidate points.
	*/
	inline void Box::Poisson(std::vector<Vector3>& p, float r, int n) const
	{
		float c = 4.0 * r * r;
		for (int i = 0; i < n; i++)
		{
			Vector3 t = RandomInside();
			bool hit = false;
			for (int j = 0; j < int(p.size()); j++)
			{
				if (squaredmagnitude(t - p.at(j)) < c)
				{
					hit = true;
					break;
				}
			}
			if (!hit)
				p.push_back(t);
		}
	}


	/*!
	\brief Compute a Poisson sphere distribution inside a box.
	This function uses a simple O(n<SUP>3</SUP>) dart throwing algorithm.
	\param array of existing samples (possibly empty)
	\param r Radius of the sphere.
	\param n Number of candidate points.
	\param spheres list of the sphere in which we don't want to add sampling points
	*/
	inline void Box::Poisson_adapted(std::vector<Vector3>& p, float r, int n, std::vector<Sphere>& spheres) const
	{
		float c = 4.0 * r * r;
		for (int i = 0; i < n; i++)
		{
			Vector3 t = RandomInside();
			bool hit = false;
			for (int j = 0; j < int(spheres.size()); j++) {
				if (spheres[j].Contains(t)) {
					hit = true;
					break;
				}
			}
			for (int j = 0; j < int(p.size()); j++)
			{
				if (squaredmagnitude(t - p.at(j)) < c)
				{
					hit = true;
					break;
				}
			}
			if (!hit)
				p.push_back(t);
		}
	}

	inline void Box::ravel(int u, int v, int w, int& idx) const {
		// Compute the index based on the cell indices and grid dimensions
		idx = u + nu * (v + nv * w);
	}

	inline void Box::unravel(int& u, int& v, int& w, int idx) const {
		// Compute the cell indices from the index and grid dimensions
		w = idx / (nu * nv);
		idx %= (nu * nv);
		v = idx / nu;
		u = idx % nu;
	}

}
