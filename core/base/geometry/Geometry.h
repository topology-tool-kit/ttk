/// \ingroup base
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2016.
///
/// \brief Minimalist namespace that handles simple geometric computations
/// (operations on std::vectors, barycentric coordinates, etc.).

#pragma once

#include <Debug.h>

namespace ttk {

  namespace Geometry {

    /// Compute the angle between two std::vectors
    /// \param vA0 xyz coordinates of vA's origin
    /// \param vA1 xyz coordinates of vA's destination
    /// \param vB0 xyz coordinates of vB's origin
    /// \param vB1 xyz coordinates of vB's destination
    template <typename T>
    T angle(const T *vA0, const T *vA1, const T *vB0, const T *vB1);

    /// Check if two 3D std::vectors vA and vB are colinear.
    /// \param vA0 xyz coordinates of vA's origin
    /// \param vA1 xyz coordinates of vA's destination
    /// \param vB0 xyz coordinates of vB's origin
    /// \param vB1 xyz coordinates of vB's destination
    /// \param coefficients Optional output std::vector of colinearity
    /// coefficients.
    /// \param tolerance Optional tolerance value (default: PREC_FLT)
    /// \returns Returns true if the std::vectors are colinear, false
    /// otherwise.
    template <typename T>
    bool areVectorsColinear(const T *vA0,
                            const T *vA1,
                            const T *vB0,
                            const T *vB1,
                            std::vector<T> *coefficients = NULL,
                            const T *tolerance = NULL);

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
    template <typename T>
    int computeBarycentricCoordinates(const T *p0,
                                      const T *p1,
                                      const T *p,
                                      std::vector<T> &baryCentrics,
                                      const int &dimension = 3);

    /// Compute the barycentric coordinates of point \p p with regard to the
    /// triangle defined by the 3D points \p p0, \p p1, and \p p2.
    /// \param p0 xyz coordinates of the first vertex of the triangle
    /// \param p1 xyz coordinates of the second vertex of the triangle
    /// \param p2 xyz coordinates of the third vertex of the triangle
    /// \param p xyz coordinates of the queried point
    /// \param baryCentrics Output barycentric coordinates (all in [0, 1] if
    /// \p p belongs to the triangle).
    /// \return Returns 0 upon success, negative values otherwise.
    template <typename T>
    int computeBarycentricCoordinates(const T *p0,
                                      const T *p1,
                                      const T *p2,
                                      const T *p,
                                      std::vector<T> &baryCentrics);

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
    template <typename T>
    bool computeSegmentIntersection(const T &xA,
                                    const T &yA,
                                    const T &xB,
                                    const T &yB,
                                    const T &xC,
                                    const T &yC,
                                    const T &xD,
                                    const T &yD,
                                    T &x,
                                    T &y);

    /// Compute the angles of a triangle
    /// \param p0 xyz coordinates of the first vertex of the triangle
    /// \param p1 xyz coordinates of the second vertex of the triangle
    /// \param p2 xyz coordinates of the third vertex of the triangle
    /// \param angles Angles (p0p1, p1p2) (p1p2, p2p0) (p2p0, p0p1)
    template <typename T>
    int computeTriangleAngles(const T *p0,
                              const T *p1,
                              const T *p2,
                              std::vector<T> &angles);

    // Get the angle opposite to edge s2 using cosine law
    /// \param s0 length of the first side of the triangle
    /// \param s1 length of the second side of the triangle
    /// \param s2 length of the third side of the triangle
    /// \param angle Angle opposite to edge s2
    /// \return Returns 0 upon success, negative values otherwise.
    template <typename T>
    int computeTriangleAngleFromSides(const T s0,
                                      const T s1,
                                      const T s2,
                                      T &angle);

    /// Compute the area of a 3D triangle.
    /// \param p0 xyz coordinates of the first vertex of the triangle
    /// \param p1 xyz coordinates of the second vertex of the triangle
    /// \param p2 xyz coordinates of the third vertex of the triangle
    /// \param area Output area.
    /// \return Returns 0 upon success, negative values otherwise.
    template <typename T>
    int computeTriangleArea(const T *p0, const T *p1, const T *p2, T &area);

    /// Compute the area of a triangle given length of its sides
    /// \param s0 length of the first side of the triangle
    /// \param s1 length of the second side of the triangle
    /// \param s2 length of the third side of the triangle
    /// \param area Output area.
    /// \return Returns 0 upon success, negative values otherwise.
    template <typename T>
    int
      computeTriangleAreaFromSides(const T s0, const T s1, const T s2, T &area);

    /// Compute the cross product of two 3D std::vectors
    /// \param vA0 xyz coordinates of vA's origin
    /// \param vA1 xyz coordinates of vA's destination
    /// \param vB0 xyz coordinates of vB's origin
    /// \param vB1 xyz coordinates of vB's destination
    /// \param crossProduct Output cross product.
    /// \return Returns 0 upon success, negative values otherwise.
    template <typename T>
    int crossProduct(const T *vA0,
                     const T *vA1,
                     const T *vB0,
                     const T *vB1,
                     std::vector<T> &crossProduct);

    /// Compute the cross product of two 3D std::vectors
    /// \param vA xyz coordinates of vA std::vector
    /// \param vB xyz coordinates of vB std::vector
    /// \param crossProduct Output cross product.
    /// \return Returns 0 upon success, negative values otherwise.
    template <typename T>
    int crossProduct(const T *vA, const T *vB, T *crossProduct);

    /// Compute the Euclidean distance between two points
    /// \param p0 xyz coordinates of the first input point.
    /// \param p1 xyz coordinates of the second input point.
    /// \param dimension Optional parameter that specifies the dimension of
    /// the point set (by default 3).
    template <typename T>
    T distance(const T *p0, const T *p1, const int &dimension = 3);

    /// Compute the dot product of two 3D std::vectors
    /// \param vA0 xyz coordinates of vA's origin
    /// \param vA1 xyz coordinates of vA's destination
    /// \param vB0 xyz coordinates of vB's origin
    /// \param vB1 xyz coordinates of vB's destination
    /// \return Returns Output dot product
    template <typename T>
    T dotProduct(const T *vA0, const T *vA1, const T *vB0, const T *vB1);

    /// Compute the dot product of two 3D std::vectors
    /// \param vA xyz coordinates of vA std::vector
    /// \param vB xyz coordinates of vB std::vector
    /// \return Returns Output dot product
    template <typename T>
    T dotProduct(const T *vA, const T *vB);

    /// Compute the bounding box of a point set.
    /// \param points Vector of points. Each entry is a std::vector whose size
    /// is
    /// equal to the dimension of the space embedding the points.
    /// \param bBox Output bounding box. The number of entries in this
    /// std::vector
    /// is equal to the dimension of the space embedding the points.
    /// \return Returns 0 upon success, negative values otherwise.
    template <typename T>
    int getBoundingBox(const std::vector<std::vector<float>> &points,
                       std::vector<std::pair<T, T>> &bBox);

    /// Check if the point \p p is inside the triangle (\p p0, \p p1, \p p2).
    /// \param p0 xyz coordinates of the first vertex of the triangle
    /// \param p1 xyz coordinates of the second vertex of the triangle
    /// \param p2 xyz coordinates of the third vertex of the triangle
    /// \param p xyz coordinates of the queried point
    /// \return Returns true if \p p is in the triangle, false otherwise.
    template <typename T>
    bool isPointInTriangle(const T *p0, const T *p1, const T *p2, const T *p);

    /// Check if a 2D point \p xy lies on a 2D segment \p AB.
    /// \param x x coordinate of the input point.
    /// \param y y coordinate of the input point.
    /// \param xA x coordinate of the first vertex of the segment [AB]
    /// \param yA y coordinate of the first vertex of the segment [AB]
    /// \param xB x coordinate of the second vertex of the segment [AB]
    /// \param yB y coordinate of the second vertex of the segment [AB]
    /// \return Returns true if the point lies on the segment, false
    /// otherwise.
    template <typename T>
    bool isPointOnSegment(const T &x,
                          const T &y,
                          const T &xA,
                          const T &yA,
                          const T &xB,
                          const T &yB);

    /// Check if a point \p p lies on a segment \p AB.
    /// \param p xyz coordinates of the input point.
    /// \param pA xyz coordinates of the first vertex of the segment [AB]
    /// \param pB xyz coordinate of the second vertex of the segment [AB]
    /// \param dimension Optional parameter that specifies the dimension of
    /// the point set (by default 3).
    /// \return Returns true if the point lies on the segment, false
    /// otherwise.
    template <typename T>
    bool isPointOnSegment(const T *p,
                          const T *pA,
                          const T *pB,
                          const int &dimension = 3);

    /// Check if all the edges of a triangle are colinear.
    /// \param p0 xyz coordinates of the first vertex of the triangle.
    /// \param p1 xyz coordinates of the second vertex of the triangle.
    /// \param p2 xyz coordinates of the third vertex of the triangle.
    /// \param tolerance Optional tolerance value (default: PREC_FLT)
    /// \return Returns true if all the edges are colinear (false otherwise).
    template <typename T>
    bool isTriangleColinear(const T *p0,
                            const T *p1,
                            const T *p2,
                            const T *tolerance = NULL);

    /// Compute the magnitude of a 3D std::vector \p v.
    /// \param v xyz coordinates of the input std::vector.
    /// \return Returns the magnitude upon success, negative values otherwise.
    template <typename T>
    T magnitude(const T *v);

    /// Compute the magnitude of a 3D std::vector.
    /// \param o xyz coordinates of the std::vector's origin
    /// \param d xyz coordinates of the std::vector's destination
    template <typename T>
    T magnitude(const T *o, const T *d);

    /// Compute the integer power of a floating-point value
    /// (std::pow is optimised for floating-point exponents)
    template <typename T>
    inline T powInt(const T val, const int n) {
      if(n < 0) {
        return 1.0 / powInt(val, -n);
      } else if(n == 0) {
        return 1;
      } else if(n == 1) {
        return val;
      } else if(n == 2) {
        return val * val;
      } else if(n == 3) {
        return val * val * val;
      } else {
        T ret = val;
        for(int i = 0; i < n - 1; ++i) {
          ret *= val;
        }
        return ret;
      }
    }

    /// Compute the nth power of ten
    template <typename T = double>
    inline T powIntTen(const int n) {
      return powInt(static_cast<T>(10), n);
    }

    /// Compute the power of an arithmetic value
    /// (redirect to std::pow with a floating-point exponent and to
    /// Geometry::powInt with an integer exponent)
    template <typename T1, typename T2>
    inline T1 pow(const T1 val, const T2 n) {
      static_assert(
        std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value,
        "pow can only be applied on arithmetic values");

      if(std::is_integral<T2>::value) {
        return powInt(val, n);
      } else if(std::is_floating_point<T2>::value) {
        return std::pow(val, n);
      }
      // this return statement should be unreachable thanks to the
      // static_assert
      return T1{};
    }

    /// Computes the difference between two vectors.
    /// \param a xyz coordinates of the first vector
    /// \param b xyz coordinates of the second vector
    /// \param out the difference between the two vectors
    /// \return Returns 0 for success, negative otherwise.
    template <typename T>
    int subtractVectors(const T *a, const T *b, T *out);

  } // namespace Geometry
} // namespace ttk
