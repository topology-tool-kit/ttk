/// \ingroup baseCode
/// \class ttk::Geometry
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2016.
/// 
/// \brief Minimalist class that handles simple geometric computations 
/// (operations on vectors, barycentric coordinates, etc.).
///
#ifndef                 _GEOMETRY_H
#define                 _GEOMETRY_H

#include                <Debug.h>

namespace ttk{

  class Geometry : public Debug{
    
    public:
      
      // 1) constructors, destructors, operators, etc.
      Geometry();
      
      virtual ~Geometry();
     
      /// Compute the angle between two vectors
      /// \param vA0 xyz coordinates of vA's origin
      /// \param vA1 xyz coordinates of vA's destination
      /// \param vB0 xyz coordinates of vB's origin
      /// \param vB1 xyz coordinates of vB's destination
      static double angle(
        const double *vA0, const double *vA1, 
        const double *vB0, const double *vB1);
      
      /// Check if two 3D vectors vA and vB are colinear.
      /// \param vA0 xyz coordinates of vA's origin
      /// \param vA1 xyz coordinates of vA's destination
      /// \param vB0 xyz coordinates of vB's origin
      /// \param vB1 xyz coordinates of vB's destination
      /// \param coefficients Optional output vector of colinearity 
      /// coefficients.
      /// \returns Returns true if the vectors are colinear, false otherwise.
      static bool areVectorsColinear(const double *vA0, const double *vA1,
        const double *vB0, const double *vB1, 
        vector<double> *coefficients = NULL, 
        const double *tolerance = NULL);
      
      /// Compute the barycentric coordinates of point \p p with regard to the 
      /// edge defined by the 3D points \p p0 and \p p1.
      /// \param p0 xyz coordinates of the first vertex of the edge 
      /// \param p1 xyz coordinates of the second vertex of the edge
      /// \param p xyz coordinates of the queried point
      /// \param baryCentrics Output barycentric coordinates (all in [0, 1] if
      /// \p p belongs to the edge).
      /// \param dimension Optional parameter that specifies the dimension of 
      /// the point set (by default 3).
      /// \return Returns 0 upon success, negative values otherwise.
      static int computeBarycentricCoordinates(
        const double *p0, const double *p1, const double *p,
        vector<double> &baryCentrics, const int &dimension = 3);
      
      /// Compute the barycentric coordinates of point \p xyz with regard to
      /// the edge defined by the points \p xy0 and \p xy1.
      /// \param x0 x coordinate of the first vertex of the edge 
      /// \param y0 y coordinate of the first vertex of the edge
      /// \param x1 x coordinate of the second vertex of the edge
      /// \param y1 y coordinate of the second vertex of the edge
      /// \param x x coordinate of the queried point
      /// \param y y coordinate of the queried point
      /// \param baryCentrics Output barycentric coordinates (all in [0, 1] if
      /// \p p belongs to the edge).
      /// \return Returns 0 upon success, negative values otherwise.
      static int computeBarycentricCoordinates(
        const double &x0, const double &y0,
        const double &x1, const double &y1,
        const double &x, const double &y,
        vector<double> &baryCentrics);
      
      /// Compute the barycentric coordinates of point \p p with regard to the 
      /// triangle defined by the 3D points \p p0, \p p1, and \p p2.
      /// \param p0 xyz coordinates of the first vertex of the triangle 
      /// \param p1 xyz coordinates of the second vertex of the triangle
      /// \param p2 xyz coordinates of the third vertex of the triangle
      /// \param p xyz coordinates of the queried point
      /// \param baryCentrics Output barycentric coordinates (all in [0, 1] if
      /// \p p belongs to the triangle).
      /// \return Returns 0 upon success, negative values otherwise.
      static int computeBarycentricCoordinates(
        const double *p0, const double *p1, const double *p2,
        const double *p,
        vector<double> &baryCentrics);
      
      /// Compute the barycentric coordinates of point \p p with regard to the 
      /// triangle defined by the points \p p0, \p p1, and \p p2.
      /// \param p0 xyz coordinates of the first vertex of the triangle 
      /// \param p1 xyz coordinates of the second vertex of the triangle
      /// \param p2 xyz coordinates of the third vertex of the triangle
      /// \param p xyz coordinates of the queried point
      /// \param baryCentrics Output barycentric coordinates (all in [0, 1] if
      /// \p p belongs to the triangle).
      /// \return Returns 0 upon success, negative values otherwise.
      /// \note Backward compatibility.
      static int computeBarycentricCoordinates(
        const float *p0, const float *p1, const float *p2,
        const float *p,
        vector<double> &baryCentrics);
      
      /// Compute the intersection between two 2D segments AB and CD.
      /// \param xA x coordinate of the first vertex of the first segment (AB)
      /// \param yA y coordinate of the first vertex of the first segment (AB)
      /// \param xB x coordinate of the second vertex of the first segment (AB)
      /// \param yB y coordinate of the second vertex of the first segment (AB)
      /// \param xC x coordinate of the first vertex of the second segment (CD)
      /// \param yC y coordinate of the first vertex of the second segment (CD)
      /// \param xD x coordinate of the second vertex of the second segment (CD)
      /// \param yD y coordinate of the second vertex of the second segment (CD)
      /// \param x x coordinate of the output intersection (if any).
      /// \param y y coordinate of the output intersection (if any).
      /// \return Returns true if the segments intersect (false otherwise).
      static bool computeSegmentIntersection(
        const double &xA, const double &yA, const double &xB, const double &yB,
        const double &xC, const double &yC, const double &xD, const double &yD,
        double &x, double &y);
     
      /// Compute the angles of a triangle
      /// \param p0 xyz coordinates of the first vertex of the triangle 
      /// \param p1 xyz coordinates of the second vertex of the triangle
      /// \param p2 xyz coordinates of the third vertex of the triangle
      /// \param angles Angles (p0p1, p1p2) (p1p2, p2p0) (p2p0, p0p1)
      static int computeTriangleAngles(
        const double *p0, const double *p1, const double *p2,
        vector<double> &angles);
      
      /// Compute the area of a 3D triangle.
      /// \param p0 xyz coordinates of the first vertex of the triangle 
      /// \param p1 xyz coordinates of the second vertex of the triangle
      /// \param p2 xyz coordinates of the third vertex of the triangle
      /// \param area Output area.
      /// \return Returns 0 upon success, negative values otherwise.
      static int computeTriangleArea(
        const double *p0, const double *p1, const double *p2,
        double &area);
      
      /// Compute the cross product of two 3D vectors
      /// \param vA0 xyz coordinates of vA's origin
      /// \param vA1 xyz coordinates of vA's destination
      /// \param vB0 xyz coordinates of vB's origin
      /// \param vB1 xyz coordinates of vB's destination
      /// \param crossProduct Output cross product.
      /// \return Returns 0 upon success, negative values otherwise.
      static int crossProduct(
        const double *vA0, const double *vA1, 
        const double *vB0, const double *vB1,
        vector<double> &crossProduct);

			/// Compute the cross product of two 3D vectors
      /// \param vA xyz coordinates of vA vector
      /// \param vB xyz coordinates of vB vector
      /// \param vC Output cross product.
      /// \return Returns 0 upon success, negative values otherwise.
			static int crossProduct(const double *vA, const double *vB,
					double *vC);
  
      /// Compute the Euclidean distance between two points
      /// \param p0 xyz coordinates of the first input point.
      /// \param p1 xyz coordinates of the second input point.
      /// \param dimension Optional parameter that specifies the dimension of 
      /// the point set (by default 3).
      static double distance(const double *p0, const double *p1,
        const int &dimension = 3);
      
      /// Compute the Euclidean distance between two points
      /// \param p0 xyz coordinates of the first input point.
      /// \param p1 xyz coordinates of the second input point.
      /// \param dimension Optional parameter that specifies the dimension of 
      /// the point set (by default 3).
      static double distance(const float *p0, const float *p1,
        const int &dimension = 3);
      
      /// Compute the dot product of two 3D vectors
      /// \param vA0 xyz coordinates of vA's origin
      /// \param vA1 xyz coordinates of vA's destination
      /// \param vB0 xyz coordinates of vB's origin
      /// \param vB1 xyz coordinates of vB's destination
      /// \return Returns Output dot product
      static double dotProduct(
        const double *vA0, const double *vA1, 
        const double *vB0, const double *vB1);

			/// Compute the dot product of two 3D vectors
      /// \param vA0 xyz coordinates of vA vector
      /// \param vB0 xyz coordinates of vB vector
      /// \return Returns Output dot product
			static double dotProduct(const double *vA,
					const double *vB);
    
      /// Compute the bounding box of a point set.
      /// \param points Vector of points. Each entry is a vector whose size is 
      /// equal to the dimension of the space embedding the points.
      /// \param bBox Output bounding box. The number of entries in this vector
      /// is equal to the dimension of the space embedding the points.
      /// \return Returns 0 upon success, negative values otherwise.
      static int getBoundingBox(const vector<vector<float> > &points,
        vector<pair<double, double>> &bBox);
      
      /// Check if the point \p p is inside the triangle (\p p0, \p p1, \p p2).
      /// \param p0 xyz coordinates of the first vertex of the triangle 
      /// \param p1 xyz coordinates of the second vertex of the triangle
      /// \param p2 xyz coordinates of the third vertex of the triangle
      /// \param p xyz coordinates of the queried point
      /// \return Returns true if \p p is in the triangle, false otherwise.
      static bool isPointInTriangle(const double *p0, const double *p1,
        const double *p2, const double *p);
      
      /// Check if the point \p p is inside the triangle (\p p0, \p p1, \p p2).
      /// \param p0 xyz coordinates of the first vertex of the triangle 
      /// \param p1 xyz coordinates of the second vertex of the triangle
      /// \param p2 xyz coordinates of the third vertex of the triangle
      /// \param p xyz coordinates of the queried point
      /// \return Returns true if \p p is in the triangle, false otherwise.
      static bool isPointInTriangle(const float *p0, const float *p1,
        const float *p2, const float *p);
      
      /// Check if a 2D point \p xy lies on a 2D segment \p AB.
      /// \param x x coordinate of the input point.
      /// \param y y coordinate of the input point.
      /// \param xA x coordinate of the first vertex of the segment [AB]
      /// \param yA y coordinate of the first vertex of the segment [AB]
      /// \param xB x coordinate of the second vertex of the segment [AB]
      /// \param yB y coordinate of the second vertex of the segment [AB]
      /// \return Returns true if the point lies on the segment, false 
      /// otherwise.
      static bool isPointOnSegment(const double &x, const double &y,
        const double &xA, const double &yA, const double &xB, const double &yB);
      
      /// Check if a point \p p lies on a segment \p AB.
      /// \param p xyz coordinates of the input point.
      /// \param pA xyz coordinates of the first vertex of the segment [AB]
      /// \param pB xyz coordinate of the second vertex of the segment [AB]
      /// \param dimension Optional parameter that specifies the dimension of 
      /// the point set (by default 3).
      /// \return Returns true if the point lies on the segment, false 
      /// otherwise.
      static bool isPointOnSegment(const double *p,
        const double *pA, const double *pB, const int &dimension = 3);
      
      /// Check if all the edges of a triangle are colinear.
      /// \param p0 xyz coordinates of the first vertex of the triangle.
      /// \param p1 xyz coordinates of the second vertex of the triangle.
      /// \param p2 xyz coordinates of the third vertex of the triangle.
      /// \return Returns true if all the edges are colinear (false otherwise).
      static bool isTriangleColinear(
        const double *p0, const double *p1, const double *p2, 
        const double *tolerance = NULL);
      
      /// Compute the magnitude of a 3D vector \p v.
      /// \param v xyz coordinates of the input vector.
      /// \return Returns the magnitude upon success, negative values otherwise.
      static double magnitude(const double *v);
      
      /// Compute the magnitude of a 3D vector.
      /// \param o xyz coordinates of the vector's origin
      /// \param d xyz coordinates of the vector's destination
      static double magnitude(const double *o, const double *d);
      
    protected:
      
  };
}

#endif
