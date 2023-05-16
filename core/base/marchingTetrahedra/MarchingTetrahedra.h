/// \ingroup base
/// \class ttk::MarchingTetrahedra
/// \author Robin G. C. Maack <maack@rptu.de>
/// \date May 2023.
///
/// \brief TTK processing package for Marching Tetra/Triangles computations
///
/// Given an input scalar field with labels attached to the point data this
/// class executes the marching tetrahedra/triangles algorithm. It has three
/// options that either separate each label with a single separating geometry
/// inbetween two labels, or a separating geometry enclosing each label
/// (detailed and fast mode).
///
/// \b Related \b publication \n
/// "Parallel Computation of Piecewise Linear Morse-Smale Segmentations" \n
/// Robin G. C. Maack, Jonas Lukasczyk, Julien Tierny, Hans Hagen,
/// Ross Maciejewski, Christoph Garth \n
/// IEEE Transactions on Visualization and Computer Graphics \n
///
/// \sa ttkMarchingTetrahedra.cpp %for a usage example.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleSegmentation_at/">Morse-Smale
///   segmentation example</a> \n

#pragma once

// base code includes
#include <Triangulation.h>

#include <queue>

/**
 * Lookup table to retrieve the indices of the edges (0-5) and triangles (6-9)
 * of the tetrahedron used for each generated triangle. The indices are ordered
 * by the combination of vertex indices that span the simplex
 * (e.g. [0,1] = edge 0, [0,2] = edge 1, [1,2] = edge 3, [0,1,3] = triangle 7).
 * Every 3 entries represent one generated triangle by using the middle of the
 * edges and triangles respectively.
 * The comments next to each line first indicate the number of unique vertex
 * labels on the tetrahedron, the 2nd label, 3rd label, 4th label, and its index
 * (the first label is alway considered to be 0).
 * This allows to directly retrieve the triangulation from this lookup table by
 * supplying a bit code, where the 4th label supplies the two LSB
 */
constexpr int tetLookupWall[28][15] = {
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (1) 0,0,0 - 0
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,1 - 1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,2 - 2
  {2, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, //    (2) 0,0,3 - 3
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,0 - 4
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,1 - 5
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,2 - 6
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,3 - 7
  {1, 3, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, //    (2) 0,2,0 - 8
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,2,1 - 9
  {1, 2, 4, 1, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1}, //       (2) 0,2,2 -10
  {5, 8, 9, 1, 9, 8, 1, 9, 3, 2, 9, 4, 2, 9, 8}, //                (3) 0,2,3 -11
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,0 -12
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,1 -13
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,2 -14
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,3 -15
  {0, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, //    (2) 1,0,0 -16
  {0, 2, 3, 2, 3, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1}, //       (2) 1,0,1 -17
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,0,2 -18
  {4, 7, 9, 0, 7, 9, 0, 3, 9, 2, 7, 9, 2, 5, 9}, //                (3) 1,0,3 -19
  {0, 1, 4, 1, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1}, //       (2) 1,1,0 -20
  {0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, //    (2) 1,1,1 -21
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,1,2 -22
  {2, 7, 8, 1, 7, 8, 0, 1, 7, 4, 7, 8, 4, 5, 8}, //                (3) 1,1,3 -23
  {3, 6, 9, 4, 6, 9, 0, 4, 6, 1, 6, 9, 1, 5, 9}, //                (3) 1,2,0 -24
  {1, 6, 8, 2, 6, 8, 2, 6, 0, 3, 8, 5, 3, 8, 6}, //                (3) 1,2,1 -25
  {0, 6, 7, 2, 6, 1, 2, 6, 7, 3, 7, 4, 3, 7, 6}, //                (3) 1,2,2 -26
  {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10} //  (4) 1,2,3 -27
};

/**
 * Lookup table providing the labels of two vertex indexes that are separated by
 * the triangle created by tetLookupWall.
 * Every two entries describe the two labels of both vertices that are separated
 * by each triangle, allowing to calculate a hash for each triangle.
 */
constexpr int tetLookupWallLabel[27][10] = {
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (1) 0,0,0 - 0
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,1 - 1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,2 - 2
  {0, 3, -1, -1, -1, -1, -1, -1, -1, -1}, //   (2) 0,0,3 - 3
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,0 - 4
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,1 - 5
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,2 - 6
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,3 - 7
  {0, 2, -1, -1, -1, -1, -1, -1, -1, -1}, //   (2) 0,2,0 - 8
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,2,1 - 9
  {0, 2, 0, 2, -1, -1, -1, -1, -1, -1}, //     (2) 0,2,2 -10
  {2, 3, 0, 2, 0, 2, 0, 3, 0, 3}, //           (3) 0,2,3 -11
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,0 -12
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,1 -13
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,2 -14
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,3 -15
  {0, 1, -1, -1, -1, -1, -1, -1, -1, -1}, //   (2) 1,0,0 -16
  {0, 1, 0, 1, -1, -1, -1, -1, -1, -1}, //     (2) 1,0,1 -17
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,0,2 -18
  {1, 3, 0, 1, 0, 1, 0, 3, 0, 3}, //           (3) 1,0,3 -19
  {0, 1, 0, 1, -1, -1, -1, -1, -1, -1}, //     (2) 1,1,0 -20
  {0, 1, -1, -1, -1, -1, -1, -1, -1, -1}, //   (2) 1,1,1 -21
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,1,2 -22
  {0, 3, 0, 1, 0, 1, 1, 3, 1, 3}, //           (3) 1,1,3 -23
  {1, 2, 0, 1, 0, 1, 0, 2, 0, 2}, //           (3) 1,2,0 -24
  {0, 2, 0, 1, 0, 1, 1, 2, 1, 2}, //           (3) 1,2,1 -25
  {0, 1, 0, 2, 0, 2, 1, 2, 1, 2} //            (3) 1,2,2 -26
};

/**
 * Lookup table providing number of triangles created by tetLookupWall.
 */
constexpr size_t tetLookupNumWallTriangles[28] = {
  0, // (1) 0,0,0 - 0
  0, // (-) 0,0,1 - 1
  0, // (-) 0,0,2 - 2
  1, // (2) 0,0,3 - 3
  0, // (-) 0,1,0 - 4
  0, // (-) 0,1,1 - 5
  0, // (-) 0,1,2 - 6
  0, // (-) 0,1,3 - 7
  1, // (2) 0,2,0 - 8
  0, // (-) 0,2,1 - 9
  2, // (2) 0,2,2 -10
  5, // (3) 0,2,3 -11
  0, // (-) 0,3,0 -12
  0, // (-) 0,3,1 -13
  0, // (-) 0,3,2 -14
  0, // (-) 0,3,3 -15
  1, // (2) 1,0,0 -16
  2, // (2) 1,0,1 -17
  0, // (-) 1,0,2 -18
  5, // (3) 1,0,3 -19
  2, // (2) 1,1,0 -20
  1, // (2) 1,1,1 -21
  0, // (-) 1,1,2 -22
  5, // (3) 1,1,3 -23
  5, // (3) 1,2,0 -24
  5, // (3) 1,2,1 -25
  5, // (3) 1,2,2 -26
  12 // (4) 1,2,3 -27
};

/**
 * Lookup table providing the number of triangles created by tetLookupWall in
 * each case, using the "detailed boundaries" Mode that splits every triangle
 * into two slightly shifted triangles.
 */
constexpr size_t tetLookupNumTrianglesDetailedBoundary[28] = {
  0, //  (1) 0,0,0 - 0
  0, //  (-) 0,0,1 - 1
  0, //  (-) 0,0,2 - 2
  2, //  (2) 0,0,3 - 3
  0, //  (-) 0,1,0 - 4
  0, //  (-) 0,1,1 - 5
  0, //  (-) 0,1,2 - 6
  0, //  (-) 0,1,3 - 7
  2, //  (2) 0,2,0 - 8
  0, //  (-) 0,2,1 - 9
  4, //  (2) 0,2,2 -10
  10, // (3) 0,2,3 -11
  0, //  (-) 0,3,0 -12
  0, //  (-) 0,3,1 -13
  0, //  (-) 0,3,2 -14
  0, //  (-) 0,3,3 -15
  2, //  (2) 1,0,0 -16
  4, //  (2) 1,0,1 -17
  0, //  (-) 1,0,2 -18
  10, // (3) 1,0,3 -19
  4, //  (2) 1,1,0 -20
  2, //  (2) 1,1,1 -21
  0, //  (-) 1,1,2 -22
  10, // (3) 1,1,3 -23
  10, // (3) 1,2,0 -24
  10, // (3) 1,2,1 -25
  10, // (3) 1,2,2 -26
  24 //  (4) 1,2,3 -27
};

/**
 * Lookup table providing the number of triangles created by tetLookupWall in
 * each case, using the basin separating boundary mode that only draws a
 * triangle if the tetrahedron has three vertices of the same label and one
 * vertex of another label.
 */
constexpr size_t tetLookupNumTrianglesBoundaries[28] = {
  0, // (1) 0,0,0 - 0
  0, // (-) 0,0,1 - 1
  0, // (-) 0,0,2 - 2
  1, // (2) 0,0,3 - 3
  0, // (-) 0,1,0 - 4
  0, // (-) 0,1,1 - 5
  0, // (-) 0,1,2 - 6
  0, // (-) 0,1,3 - 7
  1, // (2) 0,2,0 - 8
  0, // (-) 0,2,1 - 9
  0, // (2) 0,2,2 -10
  0, // (3) 0,2,3 -11
  0, // (-) 0,3,0 -12
  0, // (-) 0,3,1 -13
  0, // (-) 0,3,2 -14
  0, // (-) 0,3,3 -15
  1, // (2) 1,0,0 -16
  0, // (2) 1,0,1 -17
  0, // (-) 1,0,2 -18
  0, // (3) 1,0,3 -19
  0, // (2) 1,1,0 -20
  1, // (2) 1,1,1 -21
  0, // (-) 1,1,2 -22
  0, // (3) 1,1,3 -23
  0, // (3) 1,2,0 -24
  0, // (3) 1,2,1 -25
  0, // (3) 1,2,2 -26
  0 //  (4) 1,2,3 -27
};

/** Lookup table providing, if a lookupIndex is a multi label case. */
constexpr bool tetLookupIsMultiLabel[28] = {
  false, // (1) 0,0,0 - 0
  false, // (-) 0,0,1 - 1
  false, // (-) 0,0,2 - 2
  true, //  (2) 0,0,3 - 3
  false, // (-) 0,1,0 - 4
  false, // (-) 0,1,1 - 5
  false, // (-) 0,1,2 - 6
  false, // (-) 0,1,3 - 7
  true, //  (2) 0,2,0 - 8
  false, // (-) 0,2,1 - 9
  true, //  (2) 0,2,2 -10
  true, //  (3) 0,2,3 -11
  false, // (-) 0,3,0 -12
  false, // (-) 0,3,1 -13
  false, // (-) 0,3,2 -14
  false, // (-) 0,3,3 -15
  true, //  (2) 1,0,0 -16
  true, //  (2) 1,0,1 -17
  false, // (-) 1,0,2 -18
  true, //  (3) 1,0,3 -19
  true, //  (2) 1,1,0 -20
  true, //  (2) 1,1,1 -21
  false, // (-) 1,1,2 -22
  true, //  (3) 1,1,3 -23
  true, //  (3) 1,2,0 -24
  true, //  (3) 1,2,1 -25
  true, //  (3) 1,2,2 -26
  true //   (4) 1,2,3 -27
};

/**  Lookup table providing if a lookupIndex is a two label case. */
constexpr bool tetLookupIs2Label[28] = {
  false, // (1) 0,0,0 - 0
  false, // (-) 0,0,1 - 1
  false, // (-) 0,0,2 - 2
  true, //  (2) 0,0,3 - 3
  false, // (-) 0,1,0 - 4
  false, // (-) 0,1,1 - 5
  false, // (-) 0,1,2 - 6
  false, // (-) 0,1,3 - 7
  true, //  (2) 0,2,0 - 8
  false, // (-) 0,2,1 - 9
  true, //  (2) 0,2,2 -10
  false, // (3) 0,2,3 -11
  false, // (-) 0,3,0 -12
  false, // (-) 0,3,1 -13
  false, // (-) 0,3,2 -14
  false, // (-) 0,3,3 -15
  true, //  (2) 1,0,0 -16
  true, //  (2) 1,0,1 -17
  false, // (-) 1,0,2 -18
  false, // (3) 1,0,3 -19
  true, //  (2) 1,1,0 -20
  true, //  (2) 1,1,1 -21
  false, // (-) 1,1,2 -22
  false, // (3) 1,1,3 -23
  false, // (3) 1,2,0 -24
  false, // (3) 1,2,1 -25
  false, // (3) 1,2,2 -26
  false //  (4) 1,2,3 -27
};

/**  Lookup table providing if a lookupIndex is a three label case. */
constexpr bool tetLookupIs3Label[28] = {
  false, // (1) 0,0,0 - 0
  false, // (-) 0,0,1 - 1
  false, // (-) 0,0,2 - 2
  false, // (2) 0,0,3 - 3
  false, // (-) 0,1,0 - 4
  false, // (-) 0,1,1 - 5
  false, // (-) 0,1,2 - 6
  false, // (-) 0,1,3 - 7
  false, // (2) 0,2,0 - 8
  false, // (-) 0,2,1 - 9
  false, // (2) 0,2,2 -10
  true, //  (3) 0,2,3 -11
  false, // (-) 0,3,0 -12
  false, // (-) 0,3,1 -13
  false, // (-) 0,3,2 -14
  false, // (-) 0,3,3 -15
  false, // (2) 1,0,0 -16
  false, // (2) 1,0,1 -17
  false, // (-) 1,0,2 -18
  true, //  (3) 1,0,3 -19
  false, // (2) 1,1,0 -20
  false, // (2) 1,1,1 -21
  false, // (-) 1,1,2 -22
  true, //  (3) 1,1,3 -23
  true, //  (3) 1,2,0 -24
  true, //  (3) 1,2,1 -25
  true, //  (3) 1,2,2 -26
  false //  (4) 1,2,3 -27
};

/**
 * Lookup table providing the vertex index(local to the tetrahedron) to specify
 * the label for the triangles using the basin separating boundaries mode.
 */
constexpr int tetLookupFastCase[28] = {
  -1, // (1) 0,0,0 - 0
  -1, // (-) 0,0,1 - 1
  -1, // (-) 0,0,2 - 2
  0, //  (2) 0,0,3 - 3
  -1, // (-) 0,1,0 - 4
  -1, // (-) 0,1,1 - 5
  -1, // (-) 0,1,2 - 6
  -1, // (-) 0,1,3 - 7
  1, //  (2) 0,2,0 - 8
  -1, // (-) 0,2,1 - 9
  -1, // (2) 0,2,2 -10
  -1, // (3) 0,2,3 -11
  -1, // (-) 0,3,0 -12
  -1, // (-) 0,3,1 -13
  -1, // (-) 0,3,2 -14
  -1, // (-) 0,3,3 -15
  2, //  (2) 1,0,0 -16
  -1, // (2) 1,0,1 -17
  -1, // (-) 1,0,2 -18
  -1, // (3) 1,0,3 -19
  -1, // (2) 1,1,0 -20
  3, //  (2) 1,1,1 -21
  -1, // (-) 1,1,2 -22
  -1, // (3) 1,1,3 -23
  -1, // (3) 1,2,0 -24
  -1, // (3) 1,2,1 -25
  -1, // (3) 1,2,2 -26
  -1 //  (4) 1,2,3 -27
};

/** Provides if a tetrahedron has three vertices with the same label and one
 * vertex with another label.
 */
constexpr bool tetLookupFast[28] = {
  false, // (1) 0,0,0 - 0
  false, // (-) 0,0,1 - 1
  false, // (-) 0,0,2 - 2
  true, //  (2) 0,0,3 - 3
  false, // (-) 0,1,0 - 4
  false, // (-) 0,1,1 - 5
  false, // (-) 0,1,2 - 6
  false, // (-) 0,1,3 - 7
  true, //  (2) 0,2,0 - 8
  false, // (-) 0,2,1 - 9
  false, // (2) 0,2,2 -10
  false, // (3) 0,2,3 -11
  false, // (-) 0,3,0 -12
  false, // (-) 0,3,1 -13
  false, // (-) 0,3,2 -14
  false, // (-) 0,3,3 -15
  true, //  (2) 1,0,0 -16
  false, // (2) 1,0,1 -17
  false, // (-) 1,0,2 -18
  false, // (3) 1,0,3 -19
  false, // (2) 1,1,0 -20
  true, //  (2) 1,1,1 -21
  false, // (-) 1,1,2 -22
  false, // (3) 1,1,3 -23
  false, // (3) 1,2,0 -24
  false, // (3) 1,2,1 -25
  false, // (3) 1,2,2 -26
  false //  (4) 1,2,3 -27
};

/**
 * Provide the vertex indices, local to the tetrahedron, for the cases in the
 * lookup table tetLookupFastCase. These three vertices have the same label.
 */
constexpr int tetLookupFastTri[4][3]
  = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};

/**
 * Vertex indices for interpolation in the detailed basin boundaries mode where
 * two labels are present at one tetrahedron. Every pair of two vertex indices
 * is used to interpolate along their connecting edge. Therefore, three edges
 * are given by six vertex indices when three vertices have the same label, as
 * one triangle is enough to separate the remaining vertex. In the case of two
 * vertices with the same label, two triangles are needed, so 4 edge centers are
 * given by the 8 vertex indices.
 */
constexpr int tetLookupSplitBasins2Label[22][8] = {
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (1) 0,0,0 - 0
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,1 - 1
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,2 - 2
  {0, 3, 1, 3, 2, 3, -1, -1}, //       (2) 0,0,3 - 3
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,0 - 4
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,1 - 5
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,2 - 6
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,3 - 7
  {0, 2, 1, 2, 3, 2, -1, -1}, //       (2) 0,2,0 - 8
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,2,1 - 9
  {0, 2, 1, 3, 0, 3, 1, 2}, //         (2) 0,2,2 -10
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (3) 0,2,3 -11
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,0 -12
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,1 -13
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,2 -14
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,3 -15
  {0, 1, 2, 1, 3, 1, -1, -1}, //       (2) 1,0,0 -16
  {0, 1, 2, 3, 0, 3, 2, 1}, //         (2) 1,0,1 -17
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,0,2 -18
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (3) 1,0,3 -19
  {0, 1, 3, 2, 0, 2, 3, 1}, //         (2) 1,1,0 -20
  {1, 0, 2, 0, 3, 0, -1, -1}, //       (2) 1,1,1 -21
};

/**
 * Vertex indices for interpolation in the detailed basin boundaries mode where
 * three labels are present at one tetrahedron. In the case of three labels, two
 * vertices have the same label. For those two vertices the vertex index local
 * to the tetrahedron, the two edges connected to the vertices that are needed
 * for the triangulation, and the triangle index are given at positions 0-3 and
 * 4-7 respectively. Position 8 and 9 give the vertex index of the two remaining
 * vertices and Position 10 gives the edge index conntecting them.
 */
constexpr int tetLookupSplitBasisns3Label[27][11] = {
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (1) 0,0,0 - 0
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,1 - 1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,2 - 2
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,0,3 - 3
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,0 - 4
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,1 - 5
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,2 - 6
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,3 - 7
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,2,0 - 8
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,2,1 - 9
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,2,2 -10
  {0, 1, 2, 2, 1, 3, 4, 3, 2, 3, 5}, //            (3) 0,2,3 -11
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,0 -12
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,1 -13
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,2 -14
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,3 -15
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,0,0 -16
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,0,1 -17
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,0,2 -18
  {0, 2, 0, 1, 2, 5, 3, 3, 3, 1, 4}, //            (3) 1,0,3 -19
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,1,0 -20
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,1,1 -21
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,1,2 -22
  {1, 0, 4, 1, 2, 1, 5, 2, 0, 3, 2}, //            (3) 1,1,3 -23
  {0, 0, 1, 0, 3, 4, 5, 3, 1, 2, 3}, //            (3) 1,2,0 -24
  {1, 0, 3, 0, 3, 2, 5, 2, 0, 2, 1}, //            (3) 1,2,1 -25
  {2, 1, 3, 0, 3, 2, 4, 1, 0, 1, 0}, //            (3) 1,2,2 -26
};

/** Multi label triangles cases. */
constexpr bool triangleLookupIsMultiLabel[7] = {
  false, // (1) 0,0 - 0
  false, // (-) 0,1 - 1
  true, //  (2) 0,2 - 2
  false, // (-) 0,3 - 3
  true, //  (2) 1,0 - 4
  true, //  (2) 1,1 - 5
  true //   (3) 1,2 - 6
};

/** Two label triangle cases. */
constexpr bool triangleLookupIs2Label[7] = {
  false, // (1) 0,0 - 0
  false, // (-) 0,1 - 1
  true, //  (2) 0,2 - 2
  false, // (-) 0,3 - 3
  true, //  (2) 1,0 - 4
  true, //  (2) 1,1 - 5
  false //  (3) 1,2 - 6
};

/** Number of triangles in basin separators mode */
constexpr size_t triangleNumberLookup[7] = {
  0, // (1) 0,0 - 0
  0, // (-) 0,1 - 1
  1, // (2) 0,2 - 2
  0, // (-) 0,3 - 3
  1, // (2) 1,0 - 4
  1, // (2) 1,1 - 5
  3 //  (3) 1,2 - 6
};

/** Number of triangles in basin boundaries mode */
constexpr size_t triangleNumberLookupBoundary[7] = {
  0, // (1) 0,0 - 0
  0, // (-) 0,1 - 1
  1, // (2) 0,2 - 2
  0, // (-) 0,3 - 3
  1, // (2) 1,0 - 4
  1, // (2) 1,1 - 5
  0 //  (3) 1,2 - 6
};

/** Number of triangles in detailed basin boundaries mode */
constexpr size_t triangleNumberLookupBoundaryDetailed[7] = {
  0, // (1) 0,0 - 0
  0, // (-) 0,1 - 1
  2, // (2) 0,2 - 2
  0, // (-) 0,3 - 3
  2, // (2) 1,0 - 4
  2, // (2) 1,1 - 5
  6 //  (3) 1,2 - 6
};

/**
 * Two label triangle cases. Every 2 entries represent the edge between
 * the vertex indices, such that two edges are given for each case.
 */
constexpr int triangleLookupEdgeVerts[7][4] = {
  {-1, -1, -1, -1}, // (1) 0,0 - 0
  {-1, -1, -1, -1}, // (-) 0,1 - 1
  {0, 2, 1, 2}, //     (2) 0,2 - 2
  {-1, -1, -1, -1}, // (-) 0,3 - 3
  {0, 1, 2, 1}, //     (2) 1,0 - 4
  {1, 0, 2, 0}, //     (2) 1,1 - 5
};

/**
 * @brief Get a hash value from two keys
 *
 * @param a First Hash key
 * @param b Second hash key
 *
 * @return Hash value
 */
constexpr unsigned long long int getHash(const unsigned long long int a,
                                         const unsigned long long int b) {
  return (a * b + (a * a) + (b * b) + (a * a * a) * (b * b * b))
         % ULONG_LONG_MAX;
}

using ttk::SimplexId;

namespace ttk {
  class MarchingTetrahedra : public virtual Debug {
  public:
    MarchingTetrahedra();

    /** @brief Type of 2-separatrix output */
    enum class SURFACE_MODE {
      SM_SEPARATORS = 0,
      SM_BOUNDARIES = 1,
      SM_BOUNDARIES_DETAILED = 2
    };

    /**
     * Main function for the Marching tetrahedra computation.
     *
     * @pre MarchingTetrahedra::preconditionTriangulation must be
     * called prior to this.
     */
    template <typename dataType, typename triangulationType>
    inline int execute(const dataType *const scalars,
                       const triangulationType &triangulation);

    /**
     * Computes a binary code for each tetrahedron and counts the number of
     * triangles generated, depending on the triangleCounter used.
     *
     * @tparam dataType Scalartype
     * @tparam triangulationType Triangulationtype
     * @param[out] tetCases Binary codes
     * @param[out] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangleCounter Table with binary code to number of triangles
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <
      typename dataType,
      typename triangulationType,
      class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
    int computeMarchingCases_2D(unsigned char *const tetCases,
                                size_t *const numTriangles,
                                const dataType *const scalars,
                                const size_t *const triangleCounter,
                                const triangulationType &triangulation) const;

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 2D basin separators.
     *
     * @tparam dataType Scalartype
     * @tparam triangulationType Triangulationtype
     * @tparam std::enable_if<std::is_unsigned<dataType>::value>::type
     * @param[in] tetCases Binary codes
     * @param[in] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <
      typename dataType,
      typename triangulationType,
      class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
    int writeSeparators_2D(const unsigned char *const tetCases,
                           const size_t *numTriangles,
                           const dataType *const scalars,
                           const triangulationType &triangulation) const;

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 2D basin boundaries.
     *
     * @tparam dataType Scalartype
     * @tparam triangulationType Triangulationtype
     * @tparam std::enable_if<std::is_unsigned<dataType>::value>::type
     * @param[in] tetCases Binary codes
     * @param[in] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <
      typename dataType,
      typename triangulationType,
      class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
    int writeBoundaries_2D(const unsigned char *const tetCases,
                           const size_t *numTriangles,
                           const dataType *const scalars,
                           const triangulationType &triangulation) const;

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 2D detailed basin boundaries.
     *
     * @tparam dataType Scalartype
     * @tparam triangulationType Triangulationtype
     * @tparam std::enable_if<std::is_unsigned<dataType>::value>::type
     * @param[in] tetCases Binary codes
     * @param[in] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <
      typename dataType,
      typename triangulationType,
      class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
    int
      writeBoundariesDetailed_2D(const unsigned char *const tetCases,
                                 const size_t *numTriangles,
                                 const dataType *const scalars,
                                 const triangulationType &triangulation) const;

    /* 3D Datasets */

    /**
     * Computes a binary code for each tetrahedron and counts the number of
     * triangles generated, depending on the triangleCounter used.
     *
     * @tparam dataType Scalartype
     * @tparam triangulationType Triangulation
     * @param[out] tetCases Binary codes
     * @param[out] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangleCounter Table with binary code to number of triangles
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <
      typename dataType,
      typename triangulationType,
      class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
    int computeMarchingCases_3D(unsigned char *const tetCases,
                                size_t *const numTriangles,
                                const dataType *const scalars,
                                const size_t *const triangleCounter,
                                const triangulationType &triangulation) const;

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 3D basin separators.
     *
     * @tparam dataType Scalartype
     * @tparam triangulationType Triangulationtype
     * @tparam std::enable_if<std::is_unsigned<dataType>::value>::type
     * @param[in] tetCases Binary codes
     * @param[in] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <
      typename dataType,
      typename triangulationType,
      class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
    int writeSeparators_3D(const unsigned char *const tetCases,
                           const size_t *numTriangles,
                           const dataType *const scalars,
                           const triangulationType &triangulation) const;

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 3D basin boundaries.
     *
     * @tparam dataType Scalartype
     * @tparam triangulationType Triangulationtype
     * @tparam std::enable_if<std::is_unsigned<dataType>::value>::type
     * @param[in] tetCases Binary codes
     * @param[in] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <
      typename dataType,
      typename triangulationType,
      class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
    int writeBoundaries_3D(const unsigned char *const tetCases,
                           const size_t *numTriangles,
                           const dataType *const scalars,
                           const triangulationType &triangulation) const;

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 3D detailed basin boundaries.
     *
     * @tparam dataType Scalartype
     * @tparam triangulationType Triangulationtype
     * @tparam std::enable_if<std::is_unsigned<dataType>::value>::type
     * @param[in] tetCases Binary codes
     * @param[in] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <
      typename dataType,
      typename triangulationType,
      class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
    int
      writeBoundariesDetailed_3D(const unsigned char *const tetCases,
                                 const size_t *numTriangles,
                                 const dataType *const scalars,
                                 const triangulationType &triangulation) const;

    /**
     * Set all necessary 2-separatrix output variables
     */
    inline int
      setOutput(SimplexId *const output_numberOfPoints,
                std::vector<float> *const output_points,
                SimplexId *const output_numberOfCells,
                std::vector<SimplexId> *const output_cells_connectivity,
                std::vector<unsigned long long> *const output_cells_mscIds) {
      output_numberOfPoints_ = output_numberOfPoints;
      output_points_ = output_points;
      output_numberOfCells_ = output_numberOfCells;
      output_cells_connectivity_ = output_cells_connectivity;
      output_cells_mscIds_ = output_cells_mscIds;
      return 0;
    }

    /**
     * @brief Get the center of two points
     *
     * @param[in] pos0 Position 0
     * @param[in] pos1 Position 1
     * @param[out] incenter Resulting position
     * @return int 0 on success
     */
    inline void
      getCenter(float pos0[3], float pos1[3], float incenter[3]) const {
      incenter[0] = 0.5 * (pos0[0] + pos1[0]);
      incenter[1] = 0.5 * (pos0[1] + pos1[1]);
      incenter[2] = 0.5 * (pos0[2] + pos1[2]);
    }

    /**
     * @brief Get the center of three points
     *
     * @param[in] pos0 Position 0
     * @param[in] pos1 Position 1
     * @param[in] pos2 Position 2
     * @param[out] incenter Resulting position
     * @return int 0 on success
     */
    inline void getCenter(float pos0[3],
                          float pos1[3],
                          float pos2[3],
                          float incenter[3]) const {
      incenter[0] = 0.3333 * (pos0[0] + pos1[0] + pos2[0]);
      incenter[1] = 0.3333 * (pos0[1] + pos1[1] + pos2[1]);
      incenter[2] = 0.3333 * (pos0[2] + pos1[2] + pos2[2]);
    }

    /**
     * @brief Get the center of four points
     *
     * @param[in] pos0 Position 0
     * @param[in] pos1 Position 1
     * @param[in] pos2 Position 2
     * @param[in] pos3 Position 3
     * @param[out] incenter Resulting position
     * @return int 0 on success
     */
    inline void getCenter(float pos0[3],
                          float pos1[3],
                          float pos2[3],
                          float pos3[3],
                          float incenter[3]) const {
      incenter[0] = 0.25 * (pos0[0] + pos1[0] + pos2[0] + pos3[0]);
      incenter[1] = 0.25 * (pos0[1] + pos1[1] + pos2[1] + pos3[1]);
      incenter[2] = 0.25 * (pos0[2] + pos1[2] + pos2[2] + pos3[2]);
    }

    /**
     * @brief Interpolate between two points (lambda = 0 -> pos1 / 1 -> pos0)
     *
     * @param[in] pos0 Position 0
     * @param[in] pos1 Position 1
     * @param[in] lambda Interpolation parameter
     * @param[out] result Resulting position
     * @return int 0 on success
     */
    inline void interpolatePoints(float pos0[3],
                                  float pos1[3],
                                  float lambda,
                                  float result[3]) const {

      result[0] = lambda * pos0[0] + (1 - lambda) * pos1[0];
      result[1] = lambda * pos0[1] + (1 - lambda) * pos1[1];
      result[2] = lambda * pos0[2] + (1 - lambda) * pos1[2];
    }

  protected:
    // Output options
    SURFACE_MODE SurfaceMode{SURFACE_MODE::SM_SEPARATORS};

  private:
    // Output data
    SimplexId *output_numberOfPoints_{};
    std::vector<float> *output_points_{};
    SimplexId *output_numberOfCells_{};
    std::vector<unsigned long long> *output_cells_mscIds_{};
    std::vector<SimplexId> *output_cells_connectivity_{};
  };
} // namespace ttk

template <typename dataType, typename triangulationType>
int ttk::MarchingTetrahedra::execute(const dataType *const scalars,
                                     const triangulationType &triangulation) {

  Timer t;

  if(scalars == nullptr)
    return this->printErr("Input scalar field pointer is null.");

  const SimplexId nV = triangulation.getNumberOfVertices();
  const SimplexId nC = triangulation.getNumberOfCells();
  const int dim = triangulation.getDimensionality();
  int num_threads_current = omp_get_max_threads();

  unsigned char *tetCases = nullptr;
  size_t *numberOfTetCases = nullptr;
  unsigned long long *cScalars = nullptr;

  tetCases = (unsigned char *)malloc(nC * sizeof(unsigned char));
  numberOfTetCases = (size_t *)malloc(num_threads_current * sizeof(size_t));
  cScalars = (unsigned long long *)malloc(nV * sizeof(unsigned long long));

  for(int vert = 0; vert < nV; vert++)
    cScalars[vert] = (unsigned long long)scalars[vert];

  if(dim == 2) {
    if(SurfaceMode == SURFACE_MODE::SM_SEPARATORS) {
      computeMarchingCases_2D(tetCases, numberOfTetCases, cScalars,
                              triangleNumberLookup, triangulation);
      writeSeparators_2D(tetCases, numberOfTetCases, cScalars, triangulation);
    } else if(SurfaceMode == SURFACE_MODE::SM_BOUNDARIES) {
      computeMarchingCases_2D(tetCases, numberOfTetCases, cScalars,
                              triangleNumberLookupBoundary, triangulation);
      writeBoundaries_2D(tetCases, numberOfTetCases, cScalars, triangulation);
    } else if(SurfaceMode == SURFACE_MODE::SM_BOUNDARIES_DETAILED) {
      computeMarchingCases_2D(tetCases, numberOfTetCases, cScalars,
                              triangleNumberLookupBoundaryDetailed,
                              triangulation);
      writeBoundariesDetailed_2D(
        tetCases, numberOfTetCases, cScalars, triangulation);
    }
  } else if(dim == 3) {
    if(SurfaceMode == SURFACE_MODE::SM_SEPARATORS) {
      computeMarchingCases_3D(tetCases, numberOfTetCases, cScalars,
                              tetLookupNumWallTriangles, triangulation);
      writeSeparators_3D(tetCases, numberOfTetCases, cScalars, triangulation);
    } else if(SurfaceMode == SURFACE_MODE::SM_BOUNDARIES) {
      computeMarchingCases_3D(tetCases, numberOfTetCases, cScalars,
                              tetLookupNumTrianglesBoundaries, triangulation);
      writeBoundaries_3D(tetCases, numberOfTetCases, cScalars, triangulation);
    } else if(SurfaceMode == SURFACE_MODE::SM_BOUNDARIES_DETAILED) {
      computeMarchingCases_3D(tetCases, numberOfTetCases, cScalars,
                              tetLookupNumTrianglesDetailedBoundary,
                              triangulation);
      writeBoundariesDetailed_3D(
        tetCases, numberOfTetCases, cScalars, triangulation);
    }
  } else {
    return this->printErr("Data of dimension " + std::to_string(dim)
                          + "not suported");
  }

  this->printMsg("Data-set ("
                   + std::to_string(triangulation.getNumberOfVertices())
                   + " points) processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <
  typename dataType,
  typename triangulationType,
  class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
int ttk::MarchingTetrahedra::computeMarchingCases_2D(
  unsigned char *const tetCases,
  size_t *const numTriangles,
  const dataType *const scalars,
  const size_t *const triangleCounter,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing separator cases", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  const SimplexId numTetra = triangulation.getNumberOfCells();

#ifdef TTK_ENABLE_OPENMP
  if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

    size_t numTris = 0;

    for(SimplexId tet = 0; tet < numTetra; tet++) {
      SimplexId vertices[3];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);

      const unsigned long long msm[3]
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

      const unsigned char index0 = (msm[0] == msm[1]) ? 0x00 : 0x04; // 0 : 1
      const unsigned char index1 = (msm[0] == msm[2]) ? 0x00 : // 0
                                     (msm[1] == msm[2]) ? 0x01
                                                        : 0x02; // 1 : 2

      tetCases[tet] = index0 | index1;
      numTris += triangleCounter[tetCases[tet]];
    }

    numTriangles[0] = numTris;

#ifdef TTK_ENABLE_OPENMP
  } else {
#pragma omp parallel num_threads(threadNumber_)
    {
      SimplexId threadTriangles = 0;
      const int tid = omp_get_thread_num();

#pragma omp for schedule(static)
      for(SimplexId tet = 0; tet < numTetra; ++tet) {
        SimplexId vertices[3];
        triangulation.getCellVertex(tet, 0, vertices[0]);
        triangulation.getCellVertex(tet, 1, vertices[1]);
        triangulation.getCellVertex(tet, 2, vertices[2]);

        const unsigned long long msm[3]
          = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

        const unsigned char index0 = (msm[0] == msm[1]) ? 0x00 : 0x04; // 0 : 1
        const unsigned char index1 = (msm[0] == msm[2]) ? 0x00 : // 0
                                       (msm[1] == msm[2]) ? 0x01
                                                          : 0x02; // 1 : 2

        tetCases[tet] = index0 | index1;
        threadTriangles += triangleCounter[tetCases[tet]];
      }
      numTriangles[tid] = threadTriangles;
    }
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed separator cases", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}

template <
  typename dataType,
  typename triangulationType,
  class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
int ttk::MarchingTetrahedra::writeSeparators_2D(
  const unsigned char *const tetCases,
  const size_t *numEdges,
  const dataType *const scalars,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Writing Boundaries", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
  if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

    output_points_->resize(6 * numEdges[0]);
    output_cells_connectivity_->resize(2 * numEdges[0]);
    output_cells_mscIds_->resize(numEdges[0]);
    *output_numberOfPoints_ = 2 * numEdges[0];
    *output_numberOfCells_ = numEdges[0];

    auto &points = *output_points_;
    auto &cellsConn = *output_cells_connectivity_;
    auto &cellsMSCIds = *output_cells_mscIds_;

    float *p = points.data();
    SimplexId *c = cellsConn.data();
    unsigned long long *m = cellsMSCIds.data();
    SimplexId cellIndex = 0;
    const SimplexId numTets = triangulation.getNumberOfCells();

    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!triangleLookupIsMultiLabel[tetCases[tet]]) {
        continue;
      }

      const int *edgeVerts = triangleLookupEdgeVerts[tetCases[tet]];

      SimplexId vertices[3];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);

      float vertPos[3][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);

      const unsigned long long msm[3]
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

      if(triangleLookupIs2Label[tetCases[tet]]) {
        float eC[2][3];
        getCenter(vertPos[edgeVerts[0]], vertPos[edgeVerts[1]], eC[0]);
        getCenter(vertPos[edgeVerts[2]], vertPos[edgeVerts[3]], eC[1]);

        const unsigned long long sparseID
          = getHash(msm[edgeVerts[0]], msm[edgeVerts[1]]);

        p[0] = eC[0][0];
        p[1] = eC[0][1];
        p[2] = eC[0][2];
        p[3] = eC[1][0];
        p[4] = eC[1][1];
        p[5] = eC[1][2];
        p += 6;

        c[0] = cellIndex;
        c[1] = cellIndex + 1;
        c += 2;
        cellIndex += 2;

        m[0] = sparseID;
        m += 1;

      } else {
        float eC[4][3];
        getCenter(vertPos[0], vertPos[1], eC[0]);
        getCenter(vertPos[0], vertPos[2], eC[1]);
        getCenter(vertPos[1], vertPos[2], eC[2]);
        getCenter(vertPos[0], vertPos[1], vertPos[2], eC[3]);

        const unsigned long long sparseID[3]
          = {getHash(msm[0], msm[1]), getHash(msm[0], msm[2]),
             getHash(msm[1], msm[2])};

        p[0] = eC[0][0];
        p[1] = eC[0][1];
        p[2] = eC[0][2];
        p[3] = eC[3][0];
        p[4] = eC[3][1];
        p[5] = eC[3][2];
        p[6] = eC[1][0];
        p[7] = eC[1][1];
        p[8] = eC[1][2];
        p[9] = eC[3][0];
        p[10] = eC[3][1];
        p[11] = eC[3][2];
        p[12] = eC[2][0];
        p[13] = eC[2][1];
        p[14] = eC[2][2];
        p[15] = eC[3][0];
        p[16] = eC[3][1];
        p[17] = eC[3][2];
        p += 18;

        c[0] = cellIndex + 0;
        c[1] = cellIndex + 1;
        c[2] = cellIndex + 2;
        c[3] = cellIndex + 3;
        c[4] = cellIndex + 4;
        c[5] = cellIndex + 5;
        c += 6;
        cellIndex += 6;

        m[0] = sparseID[0];
        m[1] = sparseID[1];
        m[2] = sparseID[2];
        m += 3;
      }
    }

#ifdef TTK_ENABLE_OPENMP
  } else {
    std::vector<size_t> triangleStartIndex(threadNumber_ + 1);

    triangleStartIndex[0] = 0;

    // Count triangle number and create iterator start indices
    for(int t = 1; t <= threadNumber_; ++t) {
      triangleStartIndex[t] = numEdges[t - 1] + triangleStartIndex[t - 1];
    }

    size_t numTotalTriangles = triangleStartIndex[threadNumber_];

    output_points_->resize(9 * numTotalTriangles);
    output_cells_connectivity_->resize(3 * numTotalTriangles);
    output_cells_mscIds_->resize(numTotalTriangles);
    *output_numberOfPoints_ = 3 * numTotalTriangles;
    *output_numberOfCells_ = numTotalTriangles;

    const SimplexId numTets = triangulation.getNumberOfCells();

#pragma omp parallel num_threads(threadNumber_)
    {
      const int tid = omp_get_thread_num();
      size_t numThreadIndex = triangleStartIndex[tid];

      auto &points = *output_points_;
      auto &cellsConn = *output_cells_connectivity_;
      auto &cellsMSCIds = *output_cells_mscIds_;

      float *p = points.data();
      SimplexId *c = cellsConn.data();
      unsigned long long *m = cellsMSCIds.data();

      p += (numThreadIndex * 6);
      c += (numThreadIndex * 2);
      m += numThreadIndex;

      numThreadIndex = 2 * numThreadIndex;

#pragma omp for schedule(static)
      for(SimplexId tet = 0; tet < numTets; ++tet) {
        if(!triangleLookupIsMultiLabel[tetCases[tet]]) {
          continue;
        }

        const int *edgeVerts = triangleLookupEdgeVerts[tetCases[tet]];

        SimplexId vertices[3];
        triangulation.getCellVertex(tet, 0, vertices[0]);
        triangulation.getCellVertex(tet, 1, vertices[1]);
        triangulation.getCellVertex(tet, 2, vertices[2]);

        float vertPos[3][3];
        triangulation.getVertexPoint(
          vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);

        const unsigned long long msm[3]
          = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

        if(triangleLookupIs2Label[tetCases[tet]]) {
          float eC[2][3];
          getCenter(vertPos[edgeVerts[0]], vertPos[edgeVerts[1]], eC[0]);
          getCenter(vertPos[edgeVerts[2]], vertPos[edgeVerts[3]], eC[1]);

          const unsigned long long sparseID
            = getHash(msm[edgeVerts[0]], msm[edgeVerts[1]]);

          p[0] = eC[0][0];
          p[1] = eC[0][1];
          p[2] = eC[0][2];
          p[3] = eC[1][0];
          p[4] = eC[1][1];
          p[5] = eC[1][2];
          p += 6;

          c[0] = numThreadIndex;
          c[1] = numThreadIndex + 1;
          c += 2;
          numThreadIndex += 2;

          m[0] = sparseID;
          m += 1;

        } else {
          float eC[4][3];
          getCenter(vertPos[0], vertPos[1], eC[0]);
          getCenter(vertPos[0], vertPos[2], eC[1]);
          getCenter(vertPos[1], vertPos[2], eC[2]);
          getCenter(vertPos[0], vertPos[1], vertPos[2], eC[3]);

          const unsigned long long sparseID[3]
            = {getHash(msm[0], msm[1]), getHash(msm[0], msm[2]),
               getHash(msm[1], msm[2])};

          p[0] = eC[0][0];
          p[1] = eC[0][1];
          p[2] = eC[0][2];
          p[3] = eC[3][0];
          p[4] = eC[3][1];
          p[5] = eC[3][2];
          p[6] = eC[1][0];
          p[7] = eC[1][1];
          p[8] = eC[1][2];
          p[9] = eC[3][0];
          p[10] = eC[3][1];
          p[11] = eC[3][2];
          p[12] = eC[2][0];
          p[13] = eC[2][1];
          p[14] = eC[2][2];
          p[15] = eC[3][0];
          p[16] = eC[3][1];
          p[17] = eC[3][2];
          p += 18;

          c[0] = numThreadIndex + 0;
          c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2;
          c[3] = numThreadIndex + 3;
          c[4] = numThreadIndex + 4;
          c[5] = numThreadIndex + 5;
          c += 6;
          numThreadIndex += 6;

          m[0] = sparseID[0];
          m[1] = sparseID[1];
          m[2] = sparseID[2];
          m += 3;
        }
      }
    }
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg(
    "Wrote Boundaries", 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <
  typename dataType,
  typename triangulationType,
  class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
int ttk::MarchingTetrahedra::writeBoundaries_2D(
  const unsigned char *const tetCases,
  const size_t *numEdges,
  const dataType *const scalars,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Writing Boundaries", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
  if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

    output_points_->resize(6 * numEdges[0]);
    output_cells_connectivity_->resize(2 * numEdges[0]);
    output_cells_mscIds_->resize(numEdges[0]);
    *output_numberOfPoints_ = 2 * numEdges[0];
    *output_numberOfCells_ = numEdges[0];

    auto &points = *output_points_;
    auto &cellsConn = *output_cells_connectivity_;
    auto &cellsMSCIds = *output_cells_mscIds_;

    float *p = points.data();
    SimplexId *c = cellsConn.data();
    unsigned long long *m = cellsMSCIds.data();
    SimplexId cellIndex = 0;
    const SimplexId numTets = triangulation.getNumberOfCells();

    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!triangleLookupIsMultiLabel[tetCases[tet]]) {
        continue;
      }

      SimplexId vertices[3];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);

      float vertPos[3][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);

      const unsigned long long msm[3]
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

      if(triangleLookupIs2Label[tetCases[tet]]) {
        const int *edgeVerts = triangleLookupEdgeVerts[tetCases[tet]];

        p[0] = vertPos[edgeVerts[0]][0];
        p[1] = vertPos[edgeVerts[0]][1];
        p[2] = vertPos[edgeVerts[0]][2];
        p[3] = vertPos[edgeVerts[2]][0];
        p[4] = vertPos[edgeVerts[2]][1];
        p[5] = vertPos[edgeVerts[2]][2];
        p += 6;

        c[0] = cellIndex;
        c[1] = cellIndex + 1;
        c += 2;
        cellIndex += 2;

        m[0] = msm[edgeVerts[0]];
        m += 1;
      }
    }

#ifdef TTK_ENABLE_OPENMP
  } else {
    std::vector<size_t> triangleStartIndex(threadNumber_ + 1);

    triangleStartIndex[0] = 0;

    // Count triangle number and create iterator start indices
    for(int t = 1; t <= threadNumber_; ++t) {
      triangleStartIndex[t] = numEdges[t - 1] + triangleStartIndex[t - 1];
    }

    size_t numTotalTriangles = triangleStartIndex[threadNumber_];

    output_points_->resize(9 * numTotalTriangles);
    output_cells_connectivity_->resize(3 * numTotalTriangles);
    output_cells_mscIds_->resize(numTotalTriangles);
    *output_numberOfPoints_ = 3 * numTotalTriangles;
    *output_numberOfCells_ = numTotalTriangles;

    const SimplexId numTets = triangulation.getNumberOfCells();

#pragma omp parallel num_threads(threadNumber_)
    {
      const int tid = omp_get_thread_num();
      size_t numThreadIndex = triangleStartIndex[tid];

      auto &points = *output_points_;
      auto &cellsConn = *output_cells_connectivity_;
      auto &cellsMSCIds = *output_cells_mscIds_;

      float *p = points.data();
      SimplexId *c = cellsConn.data();
      unsigned long long *m = cellsMSCIds.data();

      p += (numThreadIndex * 6);
      c += (numThreadIndex * 2);
      m += numThreadIndex;

      numThreadIndex = 2 * numThreadIndex;

#pragma omp for schedule(static)
      for(SimplexId tet = 0; tet < numTets; ++tet) {
        if(!triangleLookupIsMultiLabel[tetCases[tet]]) {
          continue;
        }

        SimplexId vertices[3];
        triangulation.getCellVertex(tet, 0, vertices[0]);
        triangulation.getCellVertex(tet, 1, vertices[1]);
        triangulation.getCellVertex(tet, 2, vertices[2]);

        float vertPos[3][3];
        triangulation.getVertexPoint(
          vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);

        const unsigned long long msm[3]
          = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

        if(triangleLookupIs2Label[tetCases[tet]]) {
          const int *edgeVerts = triangleLookupEdgeVerts[tetCases[tet]];

          p[0] = vertPos[edgeVerts[0]][0];
          p[1] = vertPos[edgeVerts[0]][1];
          p[2] = vertPos[edgeVerts[0]][2];
          p[3] = vertPos[edgeVerts[2]][0];
          p[4] = vertPos[edgeVerts[2]][1];
          p[5] = vertPos[edgeVerts[2]][2];
          p += 6;

          c[0] = numThreadIndex;
          c[1] = numThreadIndex + 1;
          c += 2;
          numThreadIndex += 2;

          m[0] = msm[edgeVerts[0]];
          m += 1;
        }
      }
    }
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg(
    "Wrote Boundaries", 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <
  typename dataType,
  typename triangulationType,
  class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
int ttk::MarchingTetrahedra::writeBoundariesDetailed_2D(
  const unsigned char *const tetCases,
  const size_t *numEdges,
  const dataType *const scalars,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Writing Boundaries", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  constexpr float diff = 0.02;
  constexpr float d0 = 0.5 + diff;
  constexpr float d1 = 0.5 - diff;
  constexpr float dc = 2 * diff;

#ifdef TTK_ENABLE_OPENMP
  if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

    output_points_->resize(6 * numEdges[0]);
    output_cells_connectivity_->resize(2 * numEdges[0]);
    output_cells_mscIds_->resize(numEdges[0]);
    *output_numberOfPoints_ = 2 * numEdges[0];
    *output_numberOfCells_ = numEdges[0];

    auto &points = *output_points_;
    auto &cellsConn = *output_cells_connectivity_;
    auto &cellsMSCIds = *output_cells_mscIds_;

    float *p = points.data();
    SimplexId *c = cellsConn.data();
    unsigned long long *m = cellsMSCIds.data();
    SimplexId cellIndex = 0;
    const SimplexId numTets = triangulation.getNumberOfCells();

    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!triangleLookupIsMultiLabel[tetCases[tet]]) {
        continue;
      }

      const int *vIds = triangleLookupEdgeVerts[tetCases[tet]];

      SimplexId vertices[3];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);

      float vPos[3][3];
      triangulation.getVertexPoint(
        vertices[0], vPos[0][0], vPos[0][1], vPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vPos[1][0], vPos[1][1], vPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vPos[2][0], vPos[2][1], vPos[2][2]);

      const unsigned long long msm[3]
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

      if(triangleLookupIs2Label[tetCases[tet]]) {
        float vert00[3], vert01[3], vert10[3], vert11[3];

        interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
        interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
        interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d1, vert10);
        interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d1, vert11);

        p[0] = vert00[0];
        p[1] = vert00[1];
        p[2] = vert00[2];
        p[3] = vert01[0];
        p[4] = vert01[1];
        p[5] = vert01[2];
        p[6] = vert10[0];
        p[7] = vert10[1];
        p[8] = vert10[2];
        p[9] = vert11[0];
        p[10] = vert11[1];
        p[11] = vert11[2];
        p += 12;

        c[0] = cellIndex;
        c[1] = cellIndex + 1;
        c[2] = cellIndex + 2;
        c[3] = cellIndex + 3;
        c += 4;
        cellIndex += 4;

        m[0] = msm[vIds[0]];
        m[1] = msm[vIds[1]];
        m += 2;

      } else {
        float vert00[3], vert01[3], vert10[3], vert11[3], vert20[3], vert21[3];
        interpolatePoints(vPos[0], vPos[1], d0, vert00);
        interpolatePoints(vPos[0], vPos[2], d0, vert01);
        interpolatePoints(vPos[1], vPos[0], d0, vert10);
        interpolatePoints(vPos[1], vPos[2], d0, vert11);
        interpolatePoints(vPos[2], vPos[0], d0, vert20);
        interpolatePoints(vPos[2], vPos[1], d0, vert21);

        float triCenter[3];
        getCenter(vPos[0], vPos[1], vPos[2], triCenter);

        float triS0[3], triS1[3], triS2[3];
        interpolatePoints(vPos[0], triCenter, dc, triS0);
        interpolatePoints(vPos[1], triCenter, dc, triS1);
        interpolatePoints(vPos[2], triCenter, dc, triS2);

        p[0] = vert00[0];
        p[1] = vert00[1];
        p[2] = vert00[2];
        p[3] = triS0[0];
        p[4] = triS0[1];
        p[5] = triS0[2];
        p[6] = triS0[0];
        p[7] = triS0[1];
        p[8] = triS0[2];
        p[9] = vert01[0];
        p[10] = vert01[1];
        p[11] = vert01[2];

        p[12] = vert10[0];
        p[13] = vert10[1];
        p[14] = vert10[2];
        p[15] = triS1[0];
        p[16] = triS1[1];
        p[17] = triS1[2];
        p[18] = triS1[0];
        p[19] = triS1[1];
        p[20] = triS1[2];
        p[21] = vert11[0];
        p[22] = vert11[1];
        p[23] = vert11[2];

        p[24] = vert20[0];
        p[25] = vert20[1];
        p[26] = vert20[2];
        p[27] = triS2[0];
        p[28] = triS2[1];
        p[29] = triS2[2];
        p[30] = triS2[0];
        p[31] = triS2[1];
        p[32] = triS2[2];
        p[33] = vert21[0];
        p[34] = vert21[1];
        p[35] = vert21[2];
        p += 36;

        c[0] = cellIndex + 0;
        c[1] = cellIndex + 1;
        c[2] = cellIndex + 2;
        c[3] = cellIndex + 3;

        c[4] = cellIndex + 4;
        c[5] = cellIndex + 5;
        c[6] = cellIndex + 6;
        c[7] = cellIndex + 7;

        c[8] = cellIndex + 8;
        c[9] = cellIndex + 9;
        c[10] = cellIndex + 10;
        c[11] = cellIndex + 11;
        c += 12;
        cellIndex += 12;

        m[0] = msm[0];
        m[1] = msm[0];
        m[2] = msm[1];
        m[3] = msm[1];
        m[4] = msm[2];
        m[5] = msm[2];
        m += 6;
      }
    }

#ifdef TTK_ENABLE_OPENMP
  } else {
    std::vector<size_t> triangleStartIndex(threadNumber_ + 1);

    triangleStartIndex[0] = 0;

    // Count triangle number and create iterator start indices
    for(int t = 1; t <= threadNumber_; ++t) {
      triangleStartIndex[t] = numEdges[t - 1] + triangleStartIndex[t - 1];
    }

    size_t numTotalTriangles = triangleStartIndex[threadNumber_];

    output_points_->resize(9 * numTotalTriangles);
    output_cells_connectivity_->resize(3 * numTotalTriangles);
    output_cells_mscIds_->resize(numTotalTriangles);
    *output_numberOfPoints_ = 3 * numTotalTriangles;
    *output_numberOfCells_ = numTotalTriangles;

    const SimplexId numTets = triangulation.getNumberOfCells();

#pragma omp parallel num_threads(threadNumber_)
    {
      const int tid = omp_get_thread_num();
      size_t numThreadIndex = triangleStartIndex[tid];

      auto &points = *output_points_;
      auto &cellsConn = *output_cells_connectivity_;
      auto &cellsMSCIds = *output_cells_mscIds_;

      float *p = points.data();
      SimplexId *c = cellsConn.data();
      unsigned long long *m = cellsMSCIds.data();

      p += (numThreadIndex * 6);
      c += (numThreadIndex * 2);
      m += numThreadIndex;

      numThreadIndex = 2 * numThreadIndex;

#pragma omp for schedule(static)
      for(SimplexId tet = 0; tet < numTets; ++tet) {
        if(!triangleLookupIsMultiLabel[tetCases[tet]]) {
          continue;
        }

        const int *vIds = triangleLookupEdgeVerts[tetCases[tet]];

        SimplexId vertices[3];
        triangulation.getCellVertex(tet, 0, vertices[0]);
        triangulation.getCellVertex(tet, 1, vertices[1]);
        triangulation.getCellVertex(tet, 2, vertices[2]);

        float vPos[3][3];
        triangulation.getVertexPoint(
          vertices[0], vPos[0][0], vPos[0][1], vPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vPos[1][0], vPos[1][1], vPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vPos[2][0], vPos[2][1], vPos[2][2]);

        const unsigned long long msm[3]
          = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

        if(triangleLookupIs2Label[tetCases[tet]]) {
          float vert00[3], vert01[3], vert10[3], vert11[3];

          interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
          interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
          interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d1, vert10);
          interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d1, vert11);

          p[0] = vert00[0];
          p[1] = vert00[1];
          p[2] = vert00[2];
          p[3] = vert01[0];
          p[4] = vert01[1];
          p[5] = vert01[2];
          p[6] = vert10[0];
          p[7] = vert10[1];
          p[8] = vert10[2];
          p[9] = vert11[0];
          p[10] = vert11[1];
          p[11] = vert11[2];
          p += 12;

          c[0] = numThreadIndex;
          c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2;
          c[3] = numThreadIndex + 3;
          c += 4;
          numThreadIndex += 4;

          m[0] = msm[vIds[0]];
          m[1] = msm[vIds[1]];
          m += 2;

        } else {
          float vert00[3], vert01[3], vert10[3], vert11[3], vert20[3],
            vert21[3];
          interpolatePoints(vPos[0], vPos[1], d0, vert00);
          interpolatePoints(vPos[0], vPos[2], d0, vert01);
          interpolatePoints(vPos[1], vPos[0], d0, vert10);
          interpolatePoints(vPos[1], vPos[2], d0, vert11);
          interpolatePoints(vPos[2], vPos[0], d0, vert20);
          interpolatePoints(vPos[2], vPos[1], d0, vert21);

          float triCenter[3];
          getCenter(vPos[0], vPos[1], vPos[2], triCenter);

          float triS0[3], triS1[3], triS2[3];
          interpolatePoints(vPos[0], triCenter, dc, triS0);
          interpolatePoints(vPos[1], triCenter, dc, triS1);
          interpolatePoints(vPos[2], triCenter, dc, triS2);

          p[0] = vert00[0];
          p[1] = vert00[1];
          p[2] = vert00[2];
          p[3] = triS0[0];
          p[4] = triS0[1];
          p[5] = triS0[2];
          p[6] = triS0[0];
          p[7] = triS0[1];
          p[8] = triS0[2];
          p[9] = vert01[0];
          p[10] = vert01[1];
          p[11] = vert01[2];

          p[12] = vert10[0];
          p[13] = vert10[1];
          p[14] = vert10[2];
          p[15] = triS1[0];
          p[16] = triS1[1];
          p[17] = triS1[2];
          p[18] = triS1[0];
          p[19] = triS1[1];
          p[20] = triS1[2];
          p[21] = vert11[0];
          p[22] = vert11[1];
          p[23] = vert11[2];

          p[24] = vert20[0];
          p[25] = vert20[1];
          p[26] = vert20[2];
          p[27] = triS2[0];
          p[28] = triS2[1];
          p[29] = triS2[2];
          p[30] = triS2[0];
          p[31] = triS2[1];
          p[32] = triS2[2];
          p[33] = vert21[0];
          p[34] = vert21[1];
          p[35] = vert21[2];
          p += 36;

          c[0] = numThreadIndex + 0;
          c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2;
          c[3] = numThreadIndex + 3;

          c[4] = numThreadIndex + 4;
          c[5] = numThreadIndex + 5;
          c[6] = numThreadIndex + 6;
          c[7] = numThreadIndex + 7;

          c[8] = numThreadIndex + 8;
          c[9] = numThreadIndex + 9;
          c[10] = numThreadIndex + 10;
          c[11] = numThreadIndex + 11;
          c += 12;
          numThreadIndex += 12;

          m[0] = msm[0];
          m[1] = msm[0];
          m[2] = msm[1];
          m[3] = msm[1];
          m[4] = msm[2];
          m[5] = msm[2];
          m += 6;
        }
      }
    }
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg(
    "Wrote Boundaries", 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <
  typename dataType,
  typename triangulationType,
  class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
int ttk::MarchingTetrahedra::computeMarchingCases_3D(
  unsigned char *const tetCases,
  size_t *const numTriangles,
  const dataType *const scalars,
  const size_t *const triangleCounter,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing separator cases", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  const SimplexId numTetra = triangulation.getNumberOfCells();

#ifdef TTK_ENABLE_OPENMP
  if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

    size_t numTris = 0;

    for(SimplexId tet = 0; tet < numTetra; tet++) {
      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const unsigned long long msm[4]
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
           scalars[vertices[3]]};

      const unsigned char index1 = (msm[0] == msm[1]) ? 0x00 : 0x10; // 0 : 1
      const unsigned char index2 = (msm[0] == msm[2]) ? 0x00 : // 0
                                     (msm[1] == msm[2]) ? 0x04
                                                        : 0x08; // 1 : 2
      const unsigned char index3 = (msm[0] == msm[3]) ? 0x00 : // 0
                                     (msm[1] == msm[3]) ? 0x01
                                                        : // 1
                                     (msm[2] == msm[3]) ? 0x02
                                                        : 0x03; // 2 : 3

      tetCases[tet] = index1 | index2 | index3;
      numTris += triangleCounter[tetCases[tet]];
    }

    numTriangles[0] = numTris;

#ifdef TTK_ENABLE_OPENMP
  } else {
#pragma omp parallel num_threads(threadNumber_)
    {
      SimplexId threadTriangles = 0;
      const int tid = omp_get_thread_num();

#pragma omp for schedule(static)
      for(SimplexId tet = 0; tet < numTetra; ++tet) {

        SimplexId vertices[4];
        triangulation.getCellVertex(tet, 0, vertices[0]);
        triangulation.getCellVertex(tet, 1, vertices[1]);
        triangulation.getCellVertex(tet, 2, vertices[2]);
        triangulation.getCellVertex(tet, 3, vertices[3]);

        const unsigned long long msm[4]
          = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
             scalars[vertices[3]]};

        const unsigned char index1 = (msm[0] == msm[1]) ? 0x00 : 0x10; // 0 : 1
        const unsigned char index2 = (msm[0] == msm[2]) ? 0x00 : // 0
                                       (msm[1] == msm[2]) ? 0x04
                                                          : 0x08; // 1 : 2
        const unsigned char index3 = (msm[0] == msm[3]) ? 0x00 : // 0
                                       (msm[1] == msm[3]) ? 0x01
                                                          : // 1
                                       (msm[2] == msm[3]) ? 0x02
                                                          : 0x03; // 2 : 3
        tetCases[tet] = index1 | index2 | index3;

        threadTriangles += triangleCounter[tetCases[tet]];
      }
      numTriangles[tid] = threadTriangles;
    }
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed separator cases", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}

template <
  typename dataType,
  typename triangulationType,
  class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
int ttk::MarchingTetrahedra::writeSeparators_3D(
  const unsigned char *const tetCases,
  const size_t *numTriangles,
  const dataType *const scalars,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Writing separators", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
  if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

    output_points_->resize(9 * numTriangles[0]);
    output_cells_connectivity_->resize(3 * numTriangles[0]);
    output_cells_mscIds_->resize(numTriangles[0]);
    *output_numberOfPoints_ = 3 * numTriangles[0];
    *output_numberOfCells_ = numTriangles[0];

    auto &points = *output_points_;
    auto &cellsConn = *output_cells_connectivity_;
    auto &cellsMSCIds = *output_cells_mscIds_;

    float *p = points.data();
    SimplexId *c = cellsConn.data();
    unsigned long long *m = cellsMSCIds.data();
    SimplexId cellIndex = 0;
    const SimplexId numTets = triangulation.getNumberOfCells();

    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!tetLookupIsMultiLabel[tetCases[tet]]) {
        continue;
      }

      const int *tetEdgeIndices = tetLookupWall[tetCases[tet]];
      const int *tetVertLabel = tetLookupWallLabel[tetCases[tet]];

      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const unsigned long long msm[4]
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
           scalars[vertices[3]]};

      float vertPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      float eC[10][3];
      // 6 edge centers
      getCenter(vertPos[0], vertPos[1], eC[0]);
      getCenter(vertPos[0], vertPos[2], eC[1]);
      getCenter(vertPos[0], vertPos[3], eC[2]);
      getCenter(vertPos[1], vertPos[2], eC[3]);
      getCenter(vertPos[1], vertPos[3], eC[4]);
      getCenter(vertPos[2], vertPos[3], eC[5]);

      // 4 triangle centers
      getCenter(vertPos[0], vertPos[1], vertPos[2], eC[6]);
      getCenter(vertPos[0], vertPos[1], vertPos[3], eC[7]);
      getCenter(vertPos[0], vertPos[2], vertPos[3], eC[8]);
      getCenter(vertPos[1], vertPos[2], vertPos[3], eC[9]);

      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        float tetCenter[3];
        getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

        unsigned long long sparseMSIds[6]
          = {getHash(msm[0], msm[1]), getHash(msm[0], msm[2]),
             getHash(msm[0], msm[3]), getHash(msm[1], msm[2]),
             getHash(msm[1], msm[3]), getHash(msm[2], msm[3])};

        p[0] = eC[7][0];
        p[1] = eC[7][1];
        p[2] = eC[7][2];
        p[3] = eC[0][0];
        p[4] = eC[0][1];
        p[5] = eC[0][2];
        p[6] = tetCenter[0];
        p[7] = tetCenter[1];
        p[8] = tetCenter[2];
        p[9] = eC[0][0];
        p[10] = eC[0][1];
        p[11] = eC[0][2];
        p[12] = eC[6][0];
        p[13] = eC[6][1];
        p[14] = eC[6][2];
        p[15] = tetCenter[0];
        p[16] = tetCenter[1];
        p[17] = tetCenter[2];
        p[18] = eC[8][0];
        p[19] = eC[8][1];
        p[20] = eC[8][2];
        p[21] = eC[1][0];
        p[22] = eC[1][1];
        p[23] = eC[1][2];
        p[24] = tetCenter[0];
        p[25] = tetCenter[1];
        p[26] = tetCenter[2];
        p[27] = eC[1][0];
        p[28] = eC[1][1];
        p[29] = eC[1][2];
        p[30] = eC[6][0];
        p[31] = eC[6][1];
        p[32] = eC[6][2];
        p[33] = tetCenter[0];
        p[34] = tetCenter[1];
        p[35] = tetCenter[2];
        p[36] = eC[8][0];
        p[37] = eC[8][1];
        p[38] = eC[8][2];
        p[39] = eC[2][0];
        p[40] = eC[2][1];
        p[41] = eC[2][2];
        p[42] = tetCenter[0];
        p[43] = tetCenter[1];
        p[44] = tetCenter[2];
        p[45] = eC[2][0];
        p[46] = eC[2][1];
        p[47] = eC[2][2];
        p[48] = eC[7][0];
        p[49] = eC[7][1];
        p[50] = eC[7][2];
        p[51] = tetCenter[0];
        p[52] = tetCenter[1];
        p[53] = tetCenter[2];
        p[54] = eC[6][0];
        p[55] = eC[6][1];
        p[56] = eC[6][2];
        p[57] = eC[3][0];
        p[58] = eC[3][1];
        p[59] = eC[3][2];
        p[60] = tetCenter[0];
        p[61] = tetCenter[1];
        p[62] = tetCenter[2];
        p[63] = eC[3][0];
        p[64] = eC[3][1];
        p[65] = eC[3][2];
        p[66] = eC[9][0];
        p[67] = eC[9][1];
        p[68] = eC[9][2];
        p[69] = tetCenter[0];
        p[70] = tetCenter[1];
        p[71] = tetCenter[2];
        p[72] = eC[7][0];
        p[73] = eC[7][1];
        p[74] = eC[7][2];
        p[75] = eC[4][0];
        p[76] = eC[4][1];
        p[77] = eC[4][2];
        p[78] = tetCenter[0];
        p[79] = tetCenter[1];
        p[80] = tetCenter[2];
        p[81] = eC[4][0];
        p[82] = eC[4][1];
        p[83] = eC[4][2];
        p[84] = eC[9][0];
        p[85] = eC[9][1];
        p[86] = eC[9][2];
        p[87] = tetCenter[0];
        p[88] = tetCenter[1];
        p[89] = tetCenter[2];
        p[90] = eC[9][0];
        p[91] = eC[9][1];
        p[92] = eC[9][2];
        p[93] = eC[5][0];
        p[94] = eC[5][1];
        p[95] = eC[5][2];
        p[96] = tetCenter[0];
        p[97] = tetCenter[1];
        p[98] = tetCenter[2];
        p[99] = eC[5][0];
        p[100] = eC[5][1];
        p[101] = eC[5][2];
        p[102] = eC[8][0];
        p[103] = eC[8][1];
        p[104] = eC[8][2];
        p[105] = tetCenter[0];
        p[106] = tetCenter[1];
        p[107] = tetCenter[2];
        p += 108;

        c[0] = cellIndex + 0;
        c[1] = cellIndex + 1;
        c[2] = cellIndex + 2;
        c[3] = cellIndex + 3;
        c[4] = cellIndex + 4;
        c[5] = cellIndex + 5;
        c[6] = cellIndex + 6;
        c[7] = cellIndex + 7;
        c[8] = cellIndex + 8;
        c[9] = cellIndex + 9;
        c[10] = cellIndex + 10;
        c[11] = cellIndex + 11;
        c[12] = cellIndex + 12;
        c[13] = cellIndex + 13;
        c[14] = cellIndex + 14;
        c[15] = cellIndex + 15;
        c[16] = cellIndex + 16;
        c[17] = cellIndex + 17;
        c[18] = cellIndex + 18;
        c[19] = cellIndex + 19;
        c[20] = cellIndex + 20;
        c[21] = cellIndex + 21;
        c[22] = cellIndex + 22;
        c[23] = cellIndex + 23;
        c[24] = cellIndex + 24;
        c[25] = cellIndex + 25;
        c[26] = cellIndex + 26;
        c[27] = cellIndex + 27;
        c[28] = cellIndex + 28;
        c[29] = cellIndex + 29;
        c[30] = cellIndex + 30;
        c[31] = cellIndex + 31;
        c[32] = cellIndex + 32;
        c[33] = cellIndex + 33;
        c[34] = cellIndex + 34;
        c[35] = cellIndex + 35;
        c += 36;
        cellIndex += 36;

        m[0] = sparseMSIds[0];
        m[1] = sparseMSIds[0];
        m[2] = sparseMSIds[1];
        m[3] = sparseMSIds[1];
        m[4] = sparseMSIds[2];
        m[5] = sparseMSIds[2];
        m[6] = sparseMSIds[3];
        m[7] = sparseMSIds[3];
        m[8] = sparseMSIds[4];
        m[9] = sparseMSIds[4];
        m[10] = sparseMSIds[5];
        m[11] = sparseMSIds[5];
        m += 12;
      } else { // 2 or 3 labels on tetraeder
        const size_t numTris = tetLookupNumWallTriangles[tetCases[tet]];
        std::vector<unsigned long long> sparseIds(numTris);
        for(size_t t = 0; t < numTris; ++t) {
          sparseIds[t]
            = getHash(msm[tetVertLabel[t * 2]], msm[tetVertLabel[(t * 2) + 1]]);

          p[0] = eC[tetEdgeIndices[(t * 3)]][0];
          p[1] = eC[tetEdgeIndices[(t * 3)]][1];
          p[2] = eC[tetEdgeIndices[(t * 3)]][2];
          p[3] = eC[tetEdgeIndices[(t * 3) + 1]][0];
          p[4] = eC[tetEdgeIndices[(t * 3) + 1]][1];
          p[5] = eC[tetEdgeIndices[(t * 3) + 1]][2];
          p[6] = eC[tetEdgeIndices[(t * 3) + 2]][0];
          p[7] = eC[tetEdgeIndices[(t * 3) + 2]][1];
          p[8] = eC[tetEdgeIndices[(t * 3) + 2]][2];
          p += 9;

          c[0] = cellIndex + 0;
          c[1] = cellIndex + 1;
          c[2] = cellIndex + 2;
          c += 3;
          cellIndex += 3;

          m[0] = sparseIds[t];
          m += 1;
        }
      }
    }

#ifdef TTK_ENABLE_OPENMP
  } else {
    std::vector<size_t> triangleStartIndex(threadNumber_ + 1);

    triangleStartIndex[0] = 0;

    // Count triangle number and create iterator start indices
    for(int t = 1; t <= threadNumber_; ++t) {
      triangleStartIndex[t] = numTriangles[t - 1] + triangleStartIndex[t - 1];
    }

    size_t numTotalTriangles = triangleStartIndex[threadNumber_];

    output_points_->resize(9 * numTotalTriangles);
    output_cells_connectivity_->resize(3 * numTotalTriangles);
    output_cells_mscIds_->resize(numTotalTriangles);
    *output_numberOfPoints_ = 3 * numTotalTriangles;
    *output_numberOfCells_ = numTotalTriangles;

    const SimplexId numTets = triangulation.getNumberOfCells();

#pragma omp parallel num_threads(threadNumber_)
    {
      const int tid = omp_get_thread_num();
      size_t numThreadIndex = triangleStartIndex[tid];

      auto &points = *output_points_;
      auto &cellsConn = *output_cells_connectivity_;
      auto &cellsMSCIds = *output_cells_mscIds_;

      float *p = points.data();
      SimplexId *c = cellsConn.data();
      unsigned long long *m = cellsMSCIds.data();

      p += (numThreadIndex * 9);
      c += (numThreadIndex * 3);
      m += numThreadIndex;

      numThreadIndex = 3 * numThreadIndex;

#pragma omp for schedule(static)
      for(SimplexId tet = 0; tet < numTets; ++tet) {
        if(!tetLookupIsMultiLabel[tetCases[tet]]) {
          continue;
        }

        const int *tetEdgeIndices = tetLookupWall[tetCases[tet]];
        const int *tetVertLabel = tetLookupWallLabel[tetCases[tet]];

        SimplexId vertices[4];
        triangulation.getCellVertex(tet, 0, vertices[0]);
        triangulation.getCellVertex(tet, 1, vertices[1]);
        triangulation.getCellVertex(tet, 2, vertices[2]);
        triangulation.getCellVertex(tet, 3, vertices[3]);

        const unsigned long long msm[4]
          = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
             scalars[vertices[3]]};

        float vertPos[4][3];
        triangulation.getVertexPoint(
          vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
        triangulation.getVertexPoint(
          vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

        float eC[10][3];
        // 6 edge centers
        getCenter(vertPos[0], vertPos[1], eC[0]);
        getCenter(vertPos[0], vertPos[2], eC[1]);
        getCenter(vertPos[0], vertPos[3], eC[2]);
        getCenter(vertPos[1], vertPos[2], eC[3]);
        getCenter(vertPos[1], vertPos[3], eC[4]);
        getCenter(vertPos[2], vertPos[3], eC[5]);

        // 4 triangle centers
        getCenter(vertPos[0], vertPos[1], vertPos[2], eC[6]);
        getCenter(vertPos[0], vertPos[1], vertPos[3], eC[7]);
        getCenter(vertPos[0], vertPos[2], vertPos[3], eC[8]);
        getCenter(vertPos[1], vertPos[2], vertPos[3], eC[9]);

        if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
          float tetCenter[3];
          getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

          unsigned long long sparseMSIds[6]
            = {getHash(msm[0], msm[1]), getHash(msm[0], msm[2]),
               getHash(msm[0], msm[3]), getHash(msm[1], msm[2]),
               getHash(msm[1], msm[3]), getHash(msm[2], msm[3])};

          p[0] = eC[7][0];
          p[1] = eC[7][1];
          p[2] = eC[7][2];
          p[3] = eC[0][0];
          p[4] = eC[0][1];
          p[5] = eC[0][2];
          p[6] = tetCenter[0];
          p[7] = tetCenter[1];
          p[8] = tetCenter[2];
          p[9] = eC[0][0];
          p[10] = eC[0][1];
          p[11] = eC[0][2];
          p[12] = eC[6][0];
          p[13] = eC[6][1];
          p[14] = eC[6][2];
          p[15] = tetCenter[0];
          p[16] = tetCenter[1];
          p[17] = tetCenter[2];
          p[18] = eC[8][0];
          p[19] = eC[8][1];
          p[20] = eC[8][2];
          p[21] = eC[1][0];
          p[22] = eC[1][1];
          p[23] = eC[1][2];
          p[24] = tetCenter[0];
          p[25] = tetCenter[1];
          p[26] = tetCenter[2];
          p[27] = eC[1][0];
          p[28] = eC[1][1];
          p[29] = eC[1][2];
          p[30] = eC[6][0];
          p[31] = eC[6][1];
          p[32] = eC[6][2];
          p[33] = tetCenter[0];
          p[34] = tetCenter[1];
          p[35] = tetCenter[2];
          p[36] = eC[8][0];
          p[37] = eC[8][1];
          p[38] = eC[8][2];
          p[39] = eC[2][0];
          p[40] = eC[2][1];
          p[41] = eC[2][2];
          p[42] = tetCenter[0];
          p[43] = tetCenter[1];
          p[44] = tetCenter[2];
          p[45] = eC[2][0];
          p[46] = eC[2][1];
          p[47] = eC[2][2];
          p[48] = eC[7][0];
          p[49] = eC[7][1];
          p[50] = eC[7][2];
          p[51] = tetCenter[0];
          p[52] = tetCenter[1];
          p[53] = tetCenter[2];
          p[54] = eC[6][0];
          p[55] = eC[6][1];
          p[56] = eC[6][2];
          p[57] = eC[3][0];
          p[58] = eC[3][1];
          p[59] = eC[3][2];
          p[60] = tetCenter[0];
          p[61] = tetCenter[1];
          p[62] = tetCenter[2];
          p[63] = eC[3][0];
          p[64] = eC[3][1];
          p[65] = eC[3][2];
          p[66] = eC[9][0];
          p[67] = eC[9][1];
          p[68] = eC[9][2];
          p[69] = tetCenter[0];
          p[70] = tetCenter[1];
          p[71] = tetCenter[2];
          p[72] = eC[7][0];
          p[73] = eC[7][1];
          p[74] = eC[7][2];
          p[75] = eC[4][0];
          p[76] = eC[4][1];
          p[77] = eC[4][2];
          p[78] = tetCenter[0];
          p[79] = tetCenter[1];
          p[80] = tetCenter[2];
          p[81] = eC[4][0];
          p[82] = eC[4][1];
          p[83] = eC[4][2];
          p[84] = eC[9][0];
          p[85] = eC[9][1];
          p[86] = eC[9][2];
          p[87] = tetCenter[0];
          p[88] = tetCenter[1];
          p[89] = tetCenter[2];
          p[90] = eC[9][0];
          p[91] = eC[9][1];
          p[92] = eC[9][2];
          p[93] = eC[5][0];
          p[94] = eC[5][1];
          p[95] = eC[5][2];
          p[96] = tetCenter[0];
          p[97] = tetCenter[1];
          p[98] = tetCenter[2];
          p[99] = eC[5][0];
          p[100] = eC[5][1];
          p[101] = eC[5][2];
          p[102] = eC[8][0];
          p[103] = eC[8][1];
          p[104] = eC[8][2];
          p[105] = tetCenter[0];
          p[106] = tetCenter[1];
          p[107] = tetCenter[2];
          p += 108;

          c[0] = numThreadIndex + 0;
          c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2;
          c[3] = numThreadIndex + 3;
          c[4] = numThreadIndex + 4;
          c[5] = numThreadIndex + 5;
          c[6] = numThreadIndex + 6;
          c[7] = numThreadIndex + 7;
          c[8] = numThreadIndex + 8;
          c[9] = numThreadIndex + 9;
          c[10] = numThreadIndex + 10;
          c[11] = numThreadIndex + 11;
          c[12] = numThreadIndex + 12;
          c[13] = numThreadIndex + 13;
          c[14] = numThreadIndex + 14;
          c[15] = numThreadIndex + 15;
          c[16] = numThreadIndex + 16;
          c[17] = numThreadIndex + 17;
          c[18] = numThreadIndex + 18;
          c[19] = numThreadIndex + 19;
          c[20] = numThreadIndex + 20;
          c[21] = numThreadIndex + 21;
          c[22] = numThreadIndex + 22;
          c[23] = numThreadIndex + 23;
          c[24] = numThreadIndex + 24;
          c[25] = numThreadIndex + 25;
          c[26] = numThreadIndex + 26;
          c[27] = numThreadIndex + 27;
          c[28] = numThreadIndex + 28;
          c[29] = numThreadIndex + 29;
          c[30] = numThreadIndex + 30;
          c[31] = numThreadIndex + 31;
          c[32] = numThreadIndex + 32;
          c[33] = numThreadIndex + 33;
          c[34] = numThreadIndex + 34;
          c[35] = numThreadIndex + 35;
          c += 36;
          numThreadIndex += 36;

          m[0] = sparseMSIds[0];
          m[1] = sparseMSIds[0];
          m[2] = sparseMSIds[1];
          m[3] = sparseMSIds[1];
          m[4] = sparseMSIds[2];
          m[5] = sparseMSIds[2];
          m[6] = sparseMSIds[3];
          m[7] = sparseMSIds[3];
          m[8] = sparseMSIds[4];
          m[9] = sparseMSIds[4];
          m[10] = sparseMSIds[5];
          m[11] = sparseMSIds[5];
          m += 12;

        } else { // 2 or 3 labels on tetraeder
          const size_t numTris = tetLookupNumWallTriangles[tetCases[tet]];
          std::vector<unsigned long long> sparseIds(numTris);
          for(size_t t = 0; t < numTris; ++t) {
            sparseIds[t] = getHash(
              msm[tetVertLabel[t * 2]], msm[tetVertLabel[(t * 2) + 1]]);

            p[0] = eC[tetEdgeIndices[(t * 3)]][0];
            p[1] = eC[tetEdgeIndices[(t * 3)]][1];
            p[2] = eC[tetEdgeIndices[(t * 3)]][2];
            p[3] = eC[tetEdgeIndices[(t * 3) + 1]][0];
            p[4] = eC[tetEdgeIndices[(t * 3) + 1]][1];
            p[5] = eC[tetEdgeIndices[(t * 3) + 1]][2];
            p[6] = eC[tetEdgeIndices[(t * 3) + 2]][0];
            p[7] = eC[tetEdgeIndices[(t * 3) + 2]][1];
            p[8] = eC[tetEdgeIndices[(t * 3) + 2]][2];
            p += 9;

            c[0] = numThreadIndex + 0;
            c[1] = numThreadIndex + 1;
            c[2] = numThreadIndex + 2;
            c += 3;
            numThreadIndex += 3;

            m[0] = sparseIds[t];
            m += 1;
          }
        }
      }
    }
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg(
    "Wrote separators", 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <
  typename dataType,
  typename triangulationType,
  class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
int ttk::MarchingTetrahedra::writeBoundaries_3D(
  const unsigned char *const tetCases,
  const size_t *numTriangles,
  const dataType *const scalars,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Writing Boundaries", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
  if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

    output_points_->resize(9 * numTriangles[0]);
    output_cells_connectivity_->resize(3 * numTriangles[0]);
    output_cells_mscIds_->resize(numTriangles[0]);
    *output_numberOfPoints_ = 3 * numTriangles[0];
    *output_numberOfCells_ = numTriangles[0];

    auto &points = *output_points_;
    auto &cellsConn = *output_cells_connectivity_;
    auto &cellsMSCIds = *output_cells_mscIds_;

    float *p = points.data();
    SimplexId *c = cellsConn.data();
    unsigned long long *m = cellsMSCIds.data();
    SimplexId cellIndex = 0;
    const SimplexId numTets = triangulation.getNumberOfCells();

    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!tetLookupFast[tetCases[tet]]) {
        continue;
      }

      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const unsigned long long msm[4]
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
           scalars[vertices[3]]};

      const int id0 = tetLookupFastTri[tetLookupFastCase[tetCases[tet]]][0];
      const int id1 = tetLookupFastTri[tetLookupFastCase[tetCases[tet]]][1];
      const int id2 = tetLookupFastTri[tetLookupFastCase[tetCases[tet]]][2];

      float vertPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      p[0] = vertPos[id0][0];
      p[1] = vertPos[id0][1];
      p[2] = vertPos[id0][2];
      p[3] = vertPos[id1][0];
      p[4] = vertPos[id1][1];
      p[5] = vertPos[id1][2];
      p[6] = vertPos[id2][0];
      p[7] = vertPos[id2][1];
      p[8] = vertPos[id2][2];
      p += 9;

      c[0] = cellIndex + 0;
      c[1] = cellIndex + 1;
      c[2] = cellIndex + 2;
      c += 3;
      cellIndex += 3;

      m[0] = msm[id0];
      m += 1;
    }

#ifdef TTK_ENABLE_OPENMP
  } else {
    std::vector<size_t> triangleStartIndex(threadNumber_ + 1);

    triangleStartIndex[0] = 0;

    // Count triangle number and create iterator start indices
    for(int t = 1; t <= threadNumber_; ++t) {
      triangleStartIndex[t] = numTriangles[t - 1] + triangleStartIndex[t - 1];
    }

    size_t numTotalTriangles = triangleStartIndex[threadNumber_];

    output_points_->resize(9 * numTotalTriangles);
    output_cells_connectivity_->resize(3 * numTotalTriangles);
    output_cells_mscIds_->resize(numTotalTriangles);
    *output_numberOfPoints_ = 3 * numTotalTriangles;
    *output_numberOfCells_ = numTotalTriangles;

    const SimplexId numTets = triangulation.getNumberOfCells();

#pragma omp parallel num_threads(threadNumber_)
    {
      const int tid = omp_get_thread_num();
      size_t numThreadIndex = triangleStartIndex[tid];

      auto &points = *output_points_;
      auto &cellsConn = *output_cells_connectivity_;
      auto &cellsMSCIds = *output_cells_mscIds_;

      float *p = points.data();
      SimplexId *c = cellsConn.data();
      unsigned long long *m = cellsMSCIds.data();

      p += (numThreadIndex * 9);
      c += (numThreadIndex * 3);
      m += numThreadIndex;

      numThreadIndex = 3 * numThreadIndex;

#pragma omp for schedule(static)
      for(SimplexId tet = 0; tet < numTets; ++tet) {
        if(!tetLookupFast[tetCases[tet]]) {
          continue;
        }

        SimplexId vertices[4];
        triangulation.getCellVertex(tet, 0, vertices[0]);
        triangulation.getCellVertex(tet, 1, vertices[1]);
        triangulation.getCellVertex(tet, 2, vertices[2]);
        triangulation.getCellVertex(tet, 3, vertices[3]);

        const unsigned long long msm[4]
          = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
             scalars[vertices[3]]};

        const int id0 = tetLookupFastTri[tetLookupFastCase[tetCases[tet]]][0];
        const int id1 = tetLookupFastTri[tetLookupFastCase[tetCases[tet]]][1];
        const int id2 = tetLookupFastTri[tetLookupFastCase[tetCases[tet]]][2];

        float vertPos[4][3];
        triangulation.getVertexPoint(
          vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
        triangulation.getVertexPoint(
          vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

        p[0] = vertPos[id0][0];
        p[1] = vertPos[id0][1];
        p[2] = vertPos[id0][2];
        p[3] = vertPos[id1][0];
        p[4] = vertPos[id1][1];
        p[5] = vertPos[id1][2];
        p[6] = vertPos[id2][0];
        p[7] = vertPos[id2][1];
        p[8] = vertPos[id2][2];
        p += 9;

        c[0] = numThreadIndex + 0;
        c[1] = numThreadIndex + 1;
        c[2] = numThreadIndex + 2;
        c += 3;
        numThreadIndex += 3;

        m[0] = msm[id0];
        m += 1;
      }
    }
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg(
    "Wrote Boundaries", 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <
  typename dataType,
  typename triangulationType,
  class = typename std::enable_if<std::is_unsigned<dataType>::value>::type>
int ttk::MarchingTetrahedra::writeBoundariesDetailed_3D(
  const unsigned char *const tetCases,
  const size_t *numTriangles,
  const dataType *const scalars,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Writing detailed boundaries", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  constexpr float diff = 0.02;
  constexpr float d0 = 0.5 + diff;
  constexpr float d1 = 0.5 - diff;
  constexpr float dc = diff * 2;

#ifdef TTK_ENABLE_OPENMP
  if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

    output_points_->resize(9 * numTriangles[0]);
    output_cells_connectivity_->resize(3 * numTriangles[0]);
    output_cells_mscIds_->resize(numTriangles[0]);
    *output_numberOfPoints_ = 3 * numTriangles[0];
    *output_numberOfCells_ = numTriangles[0];

    auto &points = *output_points_;
    auto &cellsConn = *output_cells_connectivity_;
    auto &cellsMSCIds = *output_cells_mscIds_;

    float *p = points.data();
    SimplexId *c = cellsConn.data();
    unsigned long long *m = cellsMSCIds.data();
    SimplexId cellIndex = 0;
    const SimplexId numTets = triangulation.getNumberOfCells();

    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!tetLookupIsMultiLabel[tetCases[tet]]) {
        continue;
      }

      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const unsigned long long msm[4]
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
           scalars[vertices[3]]};

      float vPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vPos[0][0], vPos[0][1], vPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vPos[1][0], vPos[1][1], vPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vPos[2][0], vPos[2][1], vPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vPos[3][0], vPos[3][1], vPos[3][2]);

      if(tetLookupIs2Label[tetCases[tet]]) { // 2 labels (eg. AAAB / AABB)
        const int *vIds = tetLookupSplitBasins2Label[tetCases[tet]];

        float vert00[3], vert01[3], vert02[3], vert10[3], vert11[3], vert12[3];

        interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
        interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
        interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d0, vert02);
        interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d1, vert10);
        interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d1, vert11);
        interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d1, vert12);

        p[0] = vert00[0];
        p[1] = vert00[1];
        p[2] = vert00[2];
        p[3] = vert01[0];
        p[4] = vert01[1];
        p[5] = vert01[2];
        p[6] = vert02[0];
        p[7] = vert02[1];
        p[8] = vert02[2];

        p[9] = vert10[0];
        p[10] = vert10[1];
        p[11] = vert10[2];
        p[12] = vert11[0];
        p[13] = vert11[1];
        p[14] = vert11[2];
        p[15] = vert12[0];
        p[16] = vert12[1];
        p[17] = vert12[2];
        p += 18;

        c[0] = cellIndex + 0;
        c[1] = cellIndex + 1;
        c[2] = cellIndex + 2;
        c[3] = cellIndex + 3;
        c[4] = cellIndex + 4;
        c[5] = cellIndex + 5;
        c += 6;
        cellIndex += 6;

        m[0] = msm[vIds[0]];
        m[1] = msm[vIds[1]];
        m += 2;

        if(vIds[6] != -1) { // 2 vertices per label (e.g. AABB)
          float vert03[3], vert13[3];

          interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d0, vert03);
          interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d1, vert13);

          p[0] = vert00[0];
          p[1] = vert00[1];
          p[2] = vert00[2];
          p[3] = vert01[0];
          p[4] = vert01[1];
          p[5] = vert01[2];
          p[6] = vert03[0];
          p[7] = vert03[1];
          p[8] = vert03[2];
          p[9] = vert10[0];
          p[10] = vert10[1];
          p[11] = vert10[2];
          p[12] = vert11[0];
          p[13] = vert11[1];
          p[14] = vert11[2];
          p[15] = vert13[0];
          p[16] = vert13[1];
          p[17] = vert13[2];
          p += 18;

          c[0] = cellIndex + 0;
          c[1] = cellIndex + 1;
          c[2] = cellIndex + 2;
          c[3] = cellIndex + 3;
          c[4] = cellIndex + 4;
          c[5] = cellIndex + 5;
          c += 6;
          cellIndex += 6;

          m[0] = msm[vIds[0]];
          m[1] = msm[vIds[1]];
          m += 2;
        }
      } else if(tetLookupIs3Label[tetCases[tet]]) {
        const int *vIds = tetLookupSplitBasisns3Label[tetCases[tet]];

        float triCenter[4][3];
        getCenter(vPos[0], vPos[1], vPos[2], triCenter[0]);
        getCenter(vPos[0], vPos[1], vPos[3], triCenter[1]);
        getCenter(vPos[0], vPos[2], vPos[3], triCenter[2]);
        getCenter(vPos[1], vPos[2], vPos[3], triCenter[3]);

        float edgeCenters[10][3];

        getCenter(vPos[0], vPos[1], edgeCenters[0]);
        getCenter(vPos[0], vPos[2], edgeCenters[1]);
        getCenter(vPos[0], vPos[3], edgeCenters[2]);
        getCenter(vPos[1], vPos[2], edgeCenters[3]);
        getCenter(vPos[1], vPos[3], edgeCenters[4]);
        getCenter(vPos[2], vPos[3], edgeCenters[5]);

        float edge00[3], edge01[3], edge02[3], edge03[3], tri00[3], tri01[3],
          edge10[3], edge11[3], edge12[3], tri10[3], tri11[3], edge20[3],
          edge21[3], edge22[3], tri20[3], tri21[3];

        interpolatePoints(vPos[vIds[0]], edgeCenters[vIds[1]], dc, edge00);
        interpolatePoints(vPos[vIds[0]], edgeCenters[vIds[2]], dc, edge01);
        interpolatePoints(vPos[vIds[0]], triCenter[vIds[3]], dc, tri00);
        interpolatePoints(vPos[vIds[4]], edgeCenters[vIds[5]], dc, edge02);
        interpolatePoints(vPos[vIds[4]], edgeCenters[vIds[6]], dc, edge03);
        interpolatePoints(vPos[vIds[4]], triCenter[vIds[7]], dc, tri01);

        interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[1]], dc, edge10);
        interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[5]], dc, edge11);
        interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[10]], dc, edge12);
        interpolatePoints(vPos[vIds[8]], triCenter[vIds[3]], dc, tri10);
        interpolatePoints(vPos[vIds[8]], triCenter[vIds[7]], dc, tri11);

        interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[2]], dc, edge20);
        interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[6]], dc, edge21);
        interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[10]], dc, edge22);
        interpolatePoints(vPos[vIds[9]], triCenter[vIds[3]], dc, tri20);
        interpolatePoints(vPos[vIds[9]], triCenter[vIds[7]], dc, tri21);

        // Label 0
        p[0] = edge00[0];
        p[1] = edge00[1];
        p[2] = edge00[2];
        p[3] = edge02[0];
        p[4] = edge02[1];
        p[5] = edge02[2];
        p[6] = tri00[0];
        p[7] = tri00[1];
        p[8] = tri00[2];
        p[9] = edge02[0];
        p[10] = edge02[1];
        p[11] = edge02[2];
        p[12] = tri00[0];
        p[13] = tri00[1];
        p[14] = tri00[2];
        p[15] = tri01[0];
        p[16] = tri01[1];
        p[17] = tri01[2];
        p[18] = edge01[0];
        p[19] = edge01[1];
        p[20] = edge01[2];
        p[21] = edge03[0];
        p[22] = edge03[1];
        p[23] = edge03[2];
        p[24] = tri00[0];
        p[25] = tri00[1];
        p[26] = tri00[2];
        p[27] = edge03[0];
        p[28] = edge03[1];
        p[29] = edge03[2];
        p[30] = tri00[0];
        p[31] = tri00[1];
        p[32] = tri00[2];
        p[33] = tri01[0];
        p[34] = tri01[1];
        p[35] = tri01[2];

        // Label 1
        p[36] = edge10[0];
        p[37] = edge10[1];
        p[38] = edge10[2];
        p[39] = edge11[0];
        p[40] = edge11[1];
        p[41] = edge11[2];
        p[42] = tri10[0];
        p[43] = tri10[1];
        p[44] = tri10[2];
        p[45] = edge11[0];
        p[46] = edge11[1];
        p[47] = edge11[2];
        p[48] = tri10[0];
        p[49] = tri10[1];
        p[50] = tri10[2];
        p[51] = tri11[0];
        p[52] = tri11[1];
        p[53] = tri11[2];
        p[54] = edge12[0];
        p[55] = edge12[1];
        p[56] = edge12[2];
        p[57] = tri10[0];
        p[58] = tri10[1];
        p[59] = tri10[2];
        p[60] = tri11[0];
        p[61] = tri11[1];
        p[62] = tri11[2];

        // Label 2
        p[63] = edge20[0];
        p[64] = edge20[1];
        p[65] = edge20[2];
        p[66] = edge21[0];
        p[67] = edge21[1];
        p[68] = edge21[2];
        p[69] = tri20[0];
        p[70] = tri20[1];
        p[71] = tri20[2];
        p[72] = edge21[0];
        p[73] = edge21[1];
        p[74] = edge21[2];
        p[75] = tri20[0];
        p[76] = tri20[1];
        p[77] = tri20[2];
        p[78] = tri21[0];
        p[79] = tri21[1];
        p[80] = tri21[2];
        p[81] = edge22[0];
        p[82] = edge22[1];
        p[83] = edge22[2];
        p[84] = tri20[0];
        p[85] = tri20[1];
        p[86] = tri20[2];
        p[87] = tri21[0];
        p[88] = tri21[1];
        p[89] = tri21[2];
        p += 90;

        c[0] = cellIndex + 0;
        c[1] = cellIndex + 1;
        c[2] = cellIndex + 2;
        c[3] = cellIndex + 3;
        c[4] = cellIndex + 4;
        c[5] = cellIndex + 5;
        c[6] = cellIndex + 6;
        c[7] = cellIndex + 7;
        c[8] = cellIndex + 8;
        c[9] = cellIndex + 9;
        c[10] = cellIndex + 10;
        c[11] = cellIndex + 11;
        c[12] = cellIndex + 12;
        c[13] = cellIndex + 13;
        c[14] = cellIndex + 14;
        c[15] = cellIndex + 15;
        c[16] = cellIndex + 16;
        c[17] = cellIndex + 17;
        c[18] = cellIndex + 18;
        c[19] = cellIndex + 19;
        c[20] = cellIndex + 20;
        c[21] = cellIndex + 21;
        c[22] = cellIndex + 22;
        c[23] = cellIndex + 23;
        c[24] = cellIndex + 24;
        c[25] = cellIndex + 25;
        c[26] = cellIndex + 26;
        c[27] = cellIndex + 27;
        c[28] = cellIndex + 28;
        c[29] = cellIndex + 29;
        c += 30;
        cellIndex += 30;

        m[0] = msm[vIds[0]];
        m[1] = msm[vIds[0]];
        m[2] = msm[vIds[0]];
        m[3] = msm[vIds[0]];
        m[4] = msm[vIds[8]];
        m[5] = msm[vIds[8]];
        m[6] = msm[vIds[8]];
        m[7] = msm[vIds[9]];
        m[8] = msm[vIds[9]];
        m[9] = msm[vIds[9]];
        m += 10;
      } else { // 4 labels
        float tetCenter[3];
        getCenter(vPos[0], vPos[1], vPos[2], vPos[3], tetCenter);

        // the 4 triangle centers
        float triCenter[4][3];
        getCenter(vPos[0], vPos[1], vPos[2], triCenter[0]);
        getCenter(vPos[0], vPos[1], vPos[3], triCenter[1]);
        getCenter(vPos[0], vPos[2], vPos[3], triCenter[2]);
        getCenter(vPos[1], vPos[2], vPos[3], triCenter[3]);

        float vert00[3], vert01[3], vert02[3], vert0tet[3], vert0t0[3],
          vert0t1[3], vert0t2[3], vert10[3], vert11[3], vert12[3], vert1tet[3],
          vert1t0[3], vert1t1[3], vert1t2[3], vert20[3], vert21[3], vert22[3],
          vert2tet[3], vert2t0[3], vert2t1[3], vert2t2[3], vert30[3], vert31[3],
          vert32[3], vert3tet[3], vert3t0[3], vert3t1[3], vert3t2[3];

        interpolatePoints(vPos[0], vPos[1], d0, vert00);
        interpolatePoints(vPos[0], vPos[2], d0, vert01);
        interpolatePoints(vPos[0], vPos[3], d0, vert02);
        interpolatePoints(vPos[0], tetCenter, dc, vert0tet);
        interpolatePoints(vPos[0], triCenter[0], dc, vert0t0);
        interpolatePoints(vPos[0], triCenter[1], dc, vert0t1);
        interpolatePoints(vPos[0], triCenter[2], dc, vert0t2);

        interpolatePoints(vPos[1], vPos[0], d0, vert10);
        interpolatePoints(vPos[1], vPos[2], d0, vert11);
        interpolatePoints(vPos[1], vPos[3], d0, vert12);
        interpolatePoints(vPos[1], tetCenter, dc, vert1tet);
        interpolatePoints(vPos[1], triCenter[0], dc, vert1t0);
        interpolatePoints(vPos[1], triCenter[1], dc, vert1t1);
        interpolatePoints(vPos[1], triCenter[3], dc, vert1t2);

        interpolatePoints(vPos[2], vPos[0], d0, vert20);
        interpolatePoints(vPos[2], vPos[1], d0, vert21);
        interpolatePoints(vPos[2], vPos[3], d0, vert22);
        interpolatePoints(vPos[2], tetCenter, dc, vert2tet);
        interpolatePoints(vPos[2], triCenter[0], dc, vert2t0);
        interpolatePoints(vPos[2], triCenter[2], dc, vert2t1);
        interpolatePoints(vPos[2], triCenter[3], dc, vert2t2);

        interpolatePoints(vPos[3], vPos[0], d0, vert30);
        interpolatePoints(vPos[3], vPos[1], d0, vert31);
        interpolatePoints(vPos[3], vPos[2], d0, vert32);
        interpolatePoints(vPos[3], tetCenter, dc, vert3tet);
        interpolatePoints(vPos[3], triCenter[1], dc, vert3t0);
        interpolatePoints(vPos[3], triCenter[2], dc, vert3t1);
        interpolatePoints(vPos[3], triCenter[3], dc, vert3t2);

        // Label Vert 0
        p[0] = vert00[0];
        p[1] = vert00[1];
        p[2] = vert00[2];
        p[3] = vert0t0[0];
        p[4] = vert0t0[1];
        p[5] = vert0t0[2];
        p[6] = vert0tet[0];
        p[7] = vert0tet[1];
        p[8] = vert0tet[2];
        p[9] = vert00[0];
        p[10] = vert00[1];
        p[11] = vert00[2];
        p[12] = vert0t1[0];
        p[13] = vert0t1[1];
        p[14] = vert0t1[2];
        p[15] = vert0tet[0];
        p[16] = vert0tet[1];
        p[17] = vert0tet[2];
        p[18] = vert01[0];
        p[19] = vert01[1];
        p[20] = vert01[2];
        p[21] = vert0t0[0];
        p[22] = vert0t0[1];
        p[23] = vert0t0[2];
        p[24] = vert0tet[0];
        p[25] = vert0tet[1];
        p[26] = vert0tet[2];
        p[27] = vert01[0];
        p[28] = vert01[1];
        p[29] = vert01[2];
        p[30] = vert0t2[0];
        p[31] = vert0t2[1];
        p[32] = vert0t2[2];
        p[33] = vert0tet[0];
        p[34] = vert0tet[1];
        p[35] = vert0tet[2];
        p[36] = vert02[0];
        p[37] = vert02[1];
        p[38] = vert02[2];
        p[39] = vert0t2[0];
        p[40] = vert0t2[1];
        p[41] = vert0t2[2];
        p[42] = vert0tet[0];
        p[43] = vert0tet[1];
        p[44] = vert0tet[2];
        p[45] = vert02[0];
        p[46] = vert02[1];
        p[47] = vert02[2];
        p[48] = vert0t1[0];
        p[49] = vert0t1[1];
        p[50] = vert0t1[2];
        p[51] = vert0tet[0];
        p[52] = vert0tet[1];
        p[53] = vert0tet[2];

        // Label Vert 1
        p[54] = vert10[0];
        p[55] = vert10[1];
        p[56] = vert10[2];
        p[57] = vert1t0[0];
        p[58] = vert1t0[1];
        p[59] = vert1t0[2];
        p[60] = vert1tet[0];
        p[61] = vert1tet[1];
        p[62] = vert1tet[2];
        p[63] = vert10[0];
        p[64] = vert10[1];
        p[65] = vert10[2];
        p[66] = vert1t1[0];
        p[67] = vert1t1[1];
        p[68] = vert1t1[2];
        p[69] = vert1tet[0];
        p[70] = vert1tet[1];
        p[71] = vert1tet[2];
        p[72] = vert11[0];
        p[73] = vert11[1];
        p[74] = vert11[2];
        p[75] = vert1t0[0];
        p[76] = vert1t0[1];
        p[77] = vert1t0[2];
        p[78] = vert1tet[0];
        p[79] = vert1tet[1];
        p[80] = vert1tet[2];
        p[81] = vert11[0];
        p[82] = vert11[1];
        p[83] = vert11[2];
        p[84] = vert1t2[0];
        p[85] = vert1t2[1];
        p[86] = vert1t2[2];
        p[87] = vert1tet[0];
        p[88] = vert1tet[1];
        p[89] = vert1tet[2];
        p[90] = vert12[0];
        p[91] = vert12[1];
        p[92] = vert12[2];
        p[93] = vert1t2[0];
        p[94] = vert1t2[1];
        p[95] = vert1t2[2];
        p[96] = vert1tet[0];
        p[97] = vert1tet[1];
        p[98] = vert1tet[2];
        p[99] = vert12[0];
        p[100] = vert12[1];
        p[101] = vert12[2];
        p[102] = vert1t1[0];
        p[103] = vert1t1[1];
        p[104] = vert1t1[2];
        p[105] = vert1tet[0];
        p[106] = vert1tet[1];
        p[107] = vert1tet[2];

        // Label Vert 2
        p[108] = vert20[0];
        p[109] = vert20[1];
        p[110] = vert20[2];
        p[111] = vert2t0[0];
        p[112] = vert2t0[1];
        p[113] = vert2t0[2];
        p[114] = vert2tet[0];
        p[115] = vert2tet[1];
        p[116] = vert2tet[2];
        p[117] = vert20[0];
        p[118] = vert20[1];
        p[119] = vert20[2];
        p[120] = vert2t1[0];
        p[121] = vert2t1[1];
        p[122] = vert2t1[2];
        p[123] = vert2tet[0];
        p[124] = vert2tet[1];
        p[125] = vert2tet[2];
        p[126] = vert21[0];
        p[127] = vert21[1];
        p[128] = vert21[2];
        p[129] = vert2t0[0];
        p[130] = vert2t0[1];
        p[131] = vert2t0[2];
        p[132] = vert2tet[0];
        p[133] = vert2tet[1];
        p[134] = vert2tet[2];
        p[135] = vert21[0];
        p[136] = vert21[1];
        p[137] = vert21[2];
        p[138] = vert2t2[0];
        p[139] = vert2t2[1];
        p[140] = vert2t2[2];
        p[141] = vert2tet[0];
        p[142] = vert2tet[1];
        p[143] = vert2tet[2];
        p[144] = vert22[0];
        p[145] = vert22[1];
        p[146] = vert22[2];
        p[147] = vert2t2[0];
        p[148] = vert2t2[1];
        p[149] = vert2t2[2];
        p[150] = vert2tet[0];
        p[151] = vert2tet[1];
        p[152] = vert2tet[2];
        p[153] = vert22[0];
        p[154] = vert22[1];
        p[155] = vert22[2];
        p[156] = vert2t1[0];
        p[157] = vert2t1[1];
        p[158] = vert2t1[2];
        p[159] = vert2tet[0];
        p[160] = vert2tet[1];
        p[161] = vert2tet[2];

        // Label Vert 3
        p[162] = vert30[0];
        p[163] = vert30[1];
        p[164] = vert30[2];
        p[165] = vert3t0[0];
        p[166] = vert3t0[1];
        p[167] = vert3t0[2];
        p[168] = vert3tet[0];
        p[169] = vert3tet[1];
        p[170] = vert3tet[2];
        p[171] = vert30[0];
        p[172] = vert30[1];
        p[173] = vert30[2];
        p[174] = vert3t1[0];
        p[175] = vert3t1[1];
        p[176] = vert3t1[2];
        p[177] = vert3tet[0];
        p[178] = vert3tet[1];
        p[179] = vert3tet[2];
        p[180] = vert31[0];
        p[181] = vert31[1];
        p[182] = vert31[2];
        p[183] = vert3t0[0];
        p[184] = vert3t0[1];
        p[185] = vert3t0[2];
        p[186] = vert3tet[0];
        p[187] = vert3tet[1];
        p[188] = vert3tet[2];
        p[189] = vert31[0];
        p[190] = vert31[1];
        p[191] = vert31[2];
        p[192] = vert3t2[0];
        p[193] = vert3t2[1];
        p[194] = vert3t2[2];
        p[195] = vert3tet[0];
        p[196] = vert3tet[1];
        p[197] = vert3tet[2];
        p[198] = vert32[0];
        p[199] = vert32[1];
        p[200] = vert32[2];
        p[201] = vert3t2[0];
        p[202] = vert3t2[1];
        p[203] = vert3t2[2];
        p[204] = vert3tet[0];
        p[205] = vert3tet[1];
        p[206] = vert3tet[2];
        p[207] = vert32[0];
        p[208] = vert32[1];
        p[209] = vert32[2];
        p[210] = vert3t1[0];
        p[211] = vert3t1[1];
        p[212] = vert3t1[2];
        p[213] = vert3tet[0];
        p[214] = vert3tet[1];
        p[215] = vert3tet[2];
        p += 216;

        c[0] = cellIndex + 0;
        c[1] = cellIndex + 1;
        c[2] = cellIndex + 2;
        c[3] = cellIndex + 3;
        c[4] = cellIndex + 4;
        c[5] = cellIndex + 5;
        c[6] = cellIndex + 6;
        c[7] = cellIndex + 7;
        c[8] = cellIndex + 8;
        c[9] = cellIndex + 9;
        c[10] = cellIndex + 10;
        c[11] = cellIndex + 11;
        c[12] = cellIndex + 12;
        c[13] = cellIndex + 13;
        c[14] = cellIndex + 14;
        c[15] = cellIndex + 15;
        c[16] = cellIndex + 16;
        c[17] = cellIndex + 17;
        c[18] = cellIndex + 18;
        c[19] = cellIndex + 19;
        c[20] = cellIndex + 20;
        c[21] = cellIndex + 21;
        c[22] = cellIndex + 22;
        c[23] = cellIndex + 23;
        c[24] = cellIndex + 24;
        c[25] = cellIndex + 25;
        c[26] = cellIndex + 26;
        c[27] = cellIndex + 27;
        c[28] = cellIndex + 28;
        c[29] = cellIndex + 29;
        c[30] = cellIndex + 30;
        c[31] = cellIndex + 31;
        c[32] = cellIndex + 32;
        c[33] = cellIndex + 33;
        c[34] = cellIndex + 34;
        c[35] = cellIndex + 35;
        c[36] = cellIndex + 36;
        c[37] = cellIndex + 37;
        c[38] = cellIndex + 38;
        c[39] = cellIndex + 39;
        c[40] = cellIndex + 40;
        c[41] = cellIndex + 41;
        c[42] = cellIndex + 42;
        c[43] = cellIndex + 43;
        c[44] = cellIndex + 44;
        c[45] = cellIndex + 45;
        c[46] = cellIndex + 46;
        c[47] = cellIndex + 47;
        c[48] = cellIndex + 48;
        c[49] = cellIndex + 49;
        c[50] = cellIndex + 50;
        c[51] = cellIndex + 51;
        c[52] = cellIndex + 52;
        c[53] = cellIndex + 53;
        c[54] = cellIndex + 54;
        c[55] = cellIndex + 55;
        c[56] = cellIndex + 56;
        c[57] = cellIndex + 57;
        c[58] = cellIndex + 58;
        c[59] = cellIndex + 59;
        c[60] = cellIndex + 60;
        c[61] = cellIndex + 61;
        c[62] = cellIndex + 62;
        c[63] = cellIndex + 63;
        c[64] = cellIndex + 64;
        c[65] = cellIndex + 65;
        c[66] = cellIndex + 66;
        c[67] = cellIndex + 67;
        c[68] = cellIndex + 68;
        c[69] = cellIndex + 69;
        c[70] = cellIndex + 70;
        c[71] = cellIndex + 71;
        c += 72;
        cellIndex += 72;

        m[0] = msm[0];
        m[1] = msm[0];
        m[2] = msm[0];
        m[3] = msm[0];
        m[4] = msm[0];
        m[5] = msm[0];
        m[6] = msm[1];
        m[7] = msm[1];
        m[8] = msm[1];
        m[9] = msm[1];
        m[10] = msm[1];
        m[11] = msm[1];
        m[12] = msm[2];
        m[13] = msm[2];
        m[14] = msm[2];
        m[15] = msm[2];
        m[16] = msm[2];
        m[17] = msm[2];
        m[18] = msm[3];
        m[19] = msm[3];
        m[20] = msm[3];
        m[21] = msm[3];
        m[22] = msm[3];
        m[23] = msm[3];
        m += 24;
      }
    }

#ifdef TTK_ENABLE_OPENMP
  } else {
    std::vector<size_t> triangleStartIndex(threadNumber_ + 1);

    triangleStartIndex[0] = 0;

    // Count triangle number and create iterator start indices
    for(int t = 1; t <= threadNumber_; ++t) {
      triangleStartIndex[t] = numTriangles[t - 1] + triangleStartIndex[t - 1];
    }

    size_t numTotalTriangles = triangleStartIndex[threadNumber_];

    output_points_->resize(9 * numTotalTriangles);
    output_cells_connectivity_->resize(3 * numTotalTriangles);
    output_cells_mscIds_->resize(numTotalTriangles);
    *output_numberOfPoints_ = 3 * numTotalTriangles;
    *output_numberOfCells_ = numTotalTriangles;

    const SimplexId numTets = triangulation.getNumberOfCells();

#pragma omp parallel num_threads(threadNumber_)
    {
      const int tid = omp_get_thread_num();
      size_t numThreadIndex = triangleStartIndex[tid];

      auto &points = *output_points_;
      auto &cellsConn = *output_cells_connectivity_;
      auto &cellsMSCIds = *output_cells_mscIds_;

      float *p = points.data();
      SimplexId *c = cellsConn.data();
      unsigned long long *m = cellsMSCIds.data();

      p += (numThreadIndex * 9);
      c += (numThreadIndex * 3);
      m += numThreadIndex;

      numThreadIndex = 3 * numThreadIndex;

#pragma omp for schedule(static)
      for(SimplexId tet = 0; tet < numTets; ++tet) {
        if(!tetLookupIsMultiLabel[tetCases[tet]]) {
          continue;
        }

        SimplexId vertices[4];
        triangulation.getCellVertex(tet, 0, vertices[0]);
        triangulation.getCellVertex(tet, 1, vertices[1]);
        triangulation.getCellVertex(tet, 2, vertices[2]);
        triangulation.getCellVertex(tet, 3, vertices[3]);

        const unsigned long long msm[4]
          = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
             scalars[vertices[3]]};

        float vPos[4][3];
        triangulation.getVertexPoint(
          vertices[0], vPos[0][0], vPos[0][1], vPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vPos[1][0], vPos[1][1], vPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vPos[2][0], vPos[2][1], vPos[2][2]);
        triangulation.getVertexPoint(
          vertices[3], vPos[3][0], vPos[3][1], vPos[3][2]);

        if(tetLookupIs2Label[tetCases[tet]]) { // 2 labels (eg. AAAB / AABB)
          const int *vIds = tetLookupSplitBasins2Label[tetCases[tet]];

          float vert00[3], vert01[3], vert02[3], vert10[3], vert11[3],
            vert12[3];

          interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
          interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
          interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d0, vert02);
          interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d1, vert10);
          interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d1, vert11);
          interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d1, vert12);

          p[0] = vert00[0];
          p[1] = vert00[1];
          p[2] = vert00[2];
          p[3] = vert01[0];
          p[4] = vert01[1];
          p[5] = vert01[2];
          p[6] = vert02[0];
          p[7] = vert02[1];
          p[8] = vert02[2];

          p[9] = vert10[0];
          p[10] = vert10[1];
          p[11] = vert10[2];
          p[12] = vert11[0];
          p[13] = vert11[1];
          p[14] = vert11[2];
          p[15] = vert12[0];
          p[16] = vert12[1];
          p[17] = vert12[2];
          p += 18;

          c[0] = numThreadIndex + 0;
          c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2;
          c[3] = numThreadIndex + 3;
          c[4] = numThreadIndex + 4;
          c[5] = numThreadIndex + 5;
          c += 6;
          numThreadIndex += 6;

          m[0] = msm[vIds[0]];
          m[1] = msm[vIds[1]];
          m += 2;

          if(vIds[6] != -1) { // 2 vertices per label (e.g. AABB)
            float vert03[3], vert13[3];

            interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d0, vert03);
            interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d1, vert13);

            p[0] = vert00[0];
            p[1] = vert00[1];
            p[2] = vert00[2];
            p[3] = vert01[0];
            p[4] = vert01[1];
            p[5] = vert01[2];
            p[6] = vert03[0];
            p[7] = vert03[1];
            p[8] = vert03[2];
            p[9] = vert10[0];
            p[10] = vert10[1];
            p[11] = vert10[2];
            p[12] = vert11[0];
            p[13] = vert11[1];
            p[14] = vert11[2];
            p[15] = vert13[0];
            p[16] = vert13[1];
            p[17] = vert13[2];
            p += 18;

            c[0] = numThreadIndex + 0;
            c[1] = numThreadIndex + 1;
            c[2] = numThreadIndex + 2;
            c[3] = numThreadIndex + 3;
            c[4] = numThreadIndex + 4;
            c[5] = numThreadIndex + 5;
            c += 6;
            numThreadIndex += 6;

            m[0] = msm[vIds[0]];
            m[1] = msm[vIds[1]];
            m += 2;
          }
        } else if(tetLookupIs3Label[tetCases[tet]]) {
          const int *vIds = tetLookupSplitBasisns3Label[tetCases[tet]];

          float triCenter[4][3];
          getCenter(vPos[0], vPos[1], vPos[2], triCenter[0]);
          getCenter(vPos[0], vPos[1], vPos[3], triCenter[1]);
          getCenter(vPos[0], vPos[2], vPos[3], triCenter[2]);
          getCenter(vPos[1], vPos[2], vPos[3], triCenter[3]);

          float edgeCenters[10][3];

          getCenter(vPos[0], vPos[1], edgeCenters[0]);
          getCenter(vPos[0], vPos[2], edgeCenters[1]);
          getCenter(vPos[0], vPos[3], edgeCenters[2]);
          getCenter(vPos[1], vPos[2], edgeCenters[3]);
          getCenter(vPos[1], vPos[3], edgeCenters[4]);
          getCenter(vPos[2], vPos[3], edgeCenters[5]);

          float edge00[3], edge01[3], edge02[3], edge03[3], tri00[3], tri01[3],
            edge10[3], edge11[3], edge12[3], tri10[3], tri11[3], edge20[3],
            edge21[3], edge22[3], tri20[3], tri21[3];

          interpolatePoints(vPos[vIds[0]], edgeCenters[vIds[1]], dc, edge00);
          interpolatePoints(vPos[vIds[0]], edgeCenters[vIds[2]], dc, edge01);
          interpolatePoints(vPos[vIds[0]], triCenter[vIds[3]], dc, tri00);
          interpolatePoints(vPos[vIds[4]], edgeCenters[vIds[5]], dc, edge02);
          interpolatePoints(vPos[vIds[4]], edgeCenters[vIds[6]], dc, edge03);
          interpolatePoints(vPos[vIds[4]], triCenter[vIds[7]], dc, tri01);

          interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[1]], dc, edge10);
          interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[5]], dc, edge11);
          interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[10]], dc, edge12);
          interpolatePoints(vPos[vIds[8]], triCenter[vIds[3]], dc, tri10);
          interpolatePoints(vPos[vIds[8]], triCenter[vIds[7]], dc, tri11);

          interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[2]], dc, edge20);
          interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[6]], dc, edge21);
          interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[10]], dc, edge22);
          interpolatePoints(vPos[vIds[9]], triCenter[vIds[3]], dc, tri20);
          interpolatePoints(vPos[vIds[9]], triCenter[vIds[7]], dc, tri21);

          // Label 0
          p[0] = edge00[0];
          p[1] = edge00[1];
          p[2] = edge00[2];
          p[3] = edge02[0];
          p[4] = edge02[1];
          p[5] = edge02[2];
          p[6] = tri00[0];
          p[7] = tri00[1];
          p[8] = tri00[2];
          p[9] = edge02[0];
          p[10] = edge02[1];
          p[11] = edge02[2];
          p[12] = tri00[0];
          p[13] = tri00[1];
          p[14] = tri00[2];
          p[15] = tri01[0];
          p[16] = tri01[1];
          p[17] = tri01[2];
          p[18] = edge01[0];
          p[19] = edge01[1];
          p[20] = edge01[2];
          p[21] = edge03[0];
          p[22] = edge03[1];
          p[23] = edge03[2];
          p[24] = tri00[0];
          p[25] = tri00[1];
          p[26] = tri00[2];
          p[27] = edge03[0];
          p[28] = edge03[1];
          p[29] = edge03[2];
          p[30] = tri00[0];
          p[31] = tri00[1];
          p[32] = tri00[2];
          p[33] = tri01[0];
          p[34] = tri01[1];
          p[35] = tri01[2];

          // Label 1
          p[36] = edge10[0];
          p[37] = edge10[1];
          p[38] = edge10[2];
          p[39] = edge11[0];
          p[40] = edge11[1];
          p[41] = edge11[2];
          p[42] = tri10[0];
          p[43] = tri10[1];
          p[44] = tri10[2];
          p[45] = edge11[0];
          p[46] = edge11[1];
          p[47] = edge11[2];
          p[48] = tri10[0];
          p[49] = tri10[1];
          p[50] = tri10[2];
          p[51] = tri11[0];
          p[52] = tri11[1];
          p[53] = tri11[2];
          p[54] = edge12[0];
          p[55] = edge12[1];
          p[56] = edge12[2];
          p[57] = tri10[0];
          p[58] = tri10[1];
          p[59] = tri10[2];
          p[60] = tri11[0];
          p[61] = tri11[1];
          p[62] = tri11[2];

          // Label 2
          p[63] = edge20[0];
          p[64] = edge20[1];
          p[65] = edge20[2];
          p[66] = edge21[0];
          p[67] = edge21[1];
          p[68] = edge21[2];
          p[69] = tri20[0];
          p[70] = tri20[1];
          p[71] = tri20[2];
          p[72] = edge21[0];
          p[73] = edge21[1];
          p[74] = edge21[2];
          p[75] = tri20[0];
          p[76] = tri20[1];
          p[77] = tri20[2];
          p[78] = tri21[0];
          p[79] = tri21[1];
          p[80] = tri21[2];
          p[81] = edge22[0];
          p[82] = edge22[1];
          p[83] = edge22[2];
          p[84] = tri20[0];
          p[85] = tri20[1];
          p[86] = tri20[2];
          p[87] = tri21[0];
          p[88] = tri21[1];
          p[89] = tri21[2];
          p += 90;

          c[0] = numThreadIndex + 0;
          c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2;
          c[3] = numThreadIndex + 3;
          c[4] = numThreadIndex + 4;
          c[5] = numThreadIndex + 5;
          c[6] = numThreadIndex + 6;
          c[7] = numThreadIndex + 7;
          c[8] = numThreadIndex + 8;
          c[9] = numThreadIndex + 9;
          c[10] = numThreadIndex + 10;
          c[11] = numThreadIndex + 11;
          c[12] = numThreadIndex + 12;
          c[13] = numThreadIndex + 13;
          c[14] = numThreadIndex + 14;
          c[15] = numThreadIndex + 15;
          c[16] = numThreadIndex + 16;
          c[17] = numThreadIndex + 17;
          c[18] = numThreadIndex + 18;
          c[19] = numThreadIndex + 19;
          c[20] = numThreadIndex + 20;
          c[21] = numThreadIndex + 21;
          c[22] = numThreadIndex + 22;
          c[23] = numThreadIndex + 23;
          c[24] = numThreadIndex + 24;
          c[25] = numThreadIndex + 25;
          c[26] = numThreadIndex + 26;
          c[27] = numThreadIndex + 27;
          c[28] = numThreadIndex + 28;
          c[29] = numThreadIndex + 29;
          c += 30;
          numThreadIndex += 30;

          m[0] = msm[vIds[0]];
          m[1] = msm[vIds[0]];
          m[2] = msm[vIds[0]];
          m[3] = msm[vIds[0]];
          m[4] = msm[vIds[8]];
          m[5] = msm[vIds[8]];
          m[6] = msm[vIds[8]];
          m[7] = msm[vIds[9]];
          m[8] = msm[vIds[9]];
          m[9] = msm[vIds[9]];
          m += 10;

        } else { // 4 labels
          float tetCenter[3];
          getCenter(vPos[0], vPos[1], vPos[2], vPos[3], tetCenter);

          // the 4 triangle centers
          float triCenter[4][3];
          getCenter(vPos[0], vPos[1], vPos[2], triCenter[0]);
          getCenter(vPos[0], vPos[1], vPos[3], triCenter[1]);
          getCenter(vPos[0], vPos[2], vPos[3], triCenter[2]);
          getCenter(vPos[1], vPos[2], vPos[3], triCenter[3]);

          float vert00[3], vert01[3], vert02[3], vert0tet[3], vert0t0[3],
            vert0t1[3], vert0t2[3], vert10[3], vert11[3], vert12[3],
            vert1tet[3], vert1t0[3], vert1t1[3], vert1t2[3], vert20[3],
            vert21[3], vert22[3], vert2tet[3], vert2t0[3], vert2t1[3],
            vert2t2[3], vert30[3], vert31[3], vert32[3], vert3tet[3],
            vert3t0[3], vert3t1[3], vert3t2[3];

          interpolatePoints(vPos[0], vPos[1], d0, vert00);
          interpolatePoints(vPos[0], vPos[2], d0, vert01);
          interpolatePoints(vPos[0], vPos[3], d0, vert02);
          interpolatePoints(vPos[0], tetCenter, dc, vert0tet);
          interpolatePoints(vPos[0], triCenter[0], dc, vert0t0);
          interpolatePoints(vPos[0], triCenter[1], dc, vert0t1);
          interpolatePoints(vPos[0], triCenter[2], dc, vert0t2);

          interpolatePoints(vPos[1], vPos[0], d0, vert10);
          interpolatePoints(vPos[1], vPos[2], d0, vert11);
          interpolatePoints(vPos[1], vPos[3], d0, vert12);
          interpolatePoints(vPos[1], tetCenter, dc, vert1tet);
          interpolatePoints(vPos[1], triCenter[0], dc, vert1t0);
          interpolatePoints(vPos[1], triCenter[1], dc, vert1t1);
          interpolatePoints(vPos[1], triCenter[3], dc, vert1t2);

          interpolatePoints(vPos[2], vPos[0], d0, vert20);
          interpolatePoints(vPos[2], vPos[1], d0, vert21);
          interpolatePoints(vPos[2], vPos[3], d0, vert22);
          interpolatePoints(vPos[2], tetCenter, dc, vert2tet);
          interpolatePoints(vPos[2], triCenter[0], dc, vert2t0);
          interpolatePoints(vPos[2], triCenter[2], dc, vert2t1);
          interpolatePoints(vPos[2], triCenter[3], dc, vert2t2);

          interpolatePoints(vPos[3], vPos[0], d0, vert30);
          interpolatePoints(vPos[3], vPos[1], d0, vert31);
          interpolatePoints(vPos[3], vPos[2], d0, vert32);
          interpolatePoints(vPos[3], tetCenter, dc, vert3tet);
          interpolatePoints(vPos[3], triCenter[1], dc, vert3t0);
          interpolatePoints(vPos[3], triCenter[2], dc, vert3t1);
          interpolatePoints(vPos[3], triCenter[3], dc, vert3t2);

          // Label Vert 0
          p[0] = vert00[0];
          p[1] = vert00[1];
          p[2] = vert00[2];
          p[3] = vert0t0[0];
          p[4] = vert0t0[1];
          p[5] = vert0t0[2];
          p[6] = vert0tet[0];
          p[7] = vert0tet[1];
          p[8] = vert0tet[2];
          p[9] = vert00[0];
          p[10] = vert00[1];
          p[11] = vert00[2];
          p[12] = vert0t1[0];
          p[13] = vert0t1[1];
          p[14] = vert0t1[2];
          p[15] = vert0tet[0];
          p[16] = vert0tet[1];
          p[17] = vert0tet[2];
          p[18] = vert01[0];
          p[19] = vert01[1];
          p[20] = vert01[2];
          p[21] = vert0t0[0];
          p[22] = vert0t0[1];
          p[23] = vert0t0[2];
          p[24] = vert0tet[0];
          p[25] = vert0tet[1];
          p[26] = vert0tet[2];
          p[27] = vert01[0];
          p[28] = vert01[1];
          p[29] = vert01[2];
          p[30] = vert0t2[0];
          p[31] = vert0t2[1];
          p[32] = vert0t2[2];
          p[33] = vert0tet[0];
          p[34] = vert0tet[1];
          p[35] = vert0tet[2];
          p[36] = vert02[0];
          p[37] = vert02[1];
          p[38] = vert02[2];
          p[39] = vert0t2[0];
          p[40] = vert0t2[1];
          p[41] = vert0t2[2];
          p[42] = vert0tet[0];
          p[43] = vert0tet[1];
          p[44] = vert0tet[2];
          p[45] = vert02[0];
          p[46] = vert02[1];
          p[47] = vert02[2];
          p[48] = vert0t1[0];
          p[49] = vert0t1[1];
          p[50] = vert0t1[2];
          p[51] = vert0tet[0];
          p[52] = vert0tet[1];
          p[53] = vert0tet[2];

          // Label Vert 1
          p[54] = vert10[0];
          p[55] = vert10[1];
          p[56] = vert10[2];
          p[57] = vert1t0[0];
          p[58] = vert1t0[1];
          p[59] = vert1t0[2];
          p[60] = vert1tet[0];
          p[61] = vert1tet[1];
          p[62] = vert1tet[2];
          p[63] = vert10[0];
          p[64] = vert10[1];
          p[65] = vert10[2];
          p[66] = vert1t1[0];
          p[67] = vert1t1[1];
          p[68] = vert1t1[2];
          p[69] = vert1tet[0];
          p[70] = vert1tet[1];
          p[71] = vert1tet[2];
          p[72] = vert11[0];
          p[73] = vert11[1];
          p[74] = vert11[2];
          p[75] = vert1t0[0];
          p[76] = vert1t0[1];
          p[77] = vert1t0[2];
          p[78] = vert1tet[0];
          p[79] = vert1tet[1];
          p[80] = vert1tet[2];
          p[81] = vert11[0];
          p[82] = vert11[1];
          p[83] = vert11[2];
          p[84] = vert1t2[0];
          p[85] = vert1t2[1];
          p[86] = vert1t2[2];
          p[87] = vert1tet[0];
          p[88] = vert1tet[1];
          p[89] = vert1tet[2];
          p[90] = vert12[0];
          p[91] = vert12[1];
          p[92] = vert12[2];
          p[93] = vert1t2[0];
          p[94] = vert1t2[1];
          p[95] = vert1t2[2];
          p[96] = vert1tet[0];
          p[97] = vert1tet[1];
          p[98] = vert1tet[2];
          p[99] = vert12[0];
          p[100] = vert12[1];
          p[101] = vert12[2];
          p[102] = vert1t1[0];
          p[103] = vert1t1[1];
          p[104] = vert1t1[2];
          p[105] = vert1tet[0];
          p[106] = vert1tet[1];
          p[107] = vert1tet[2];

          // Label Vert 2
          p[108] = vert20[0];
          p[109] = vert20[1];
          p[110] = vert20[2];
          p[111] = vert2t0[0];
          p[112] = vert2t0[1];
          p[113] = vert2t0[2];
          p[114] = vert2tet[0];
          p[115] = vert2tet[1];
          p[116] = vert2tet[2];
          p[117] = vert20[0];
          p[118] = vert20[1];
          p[119] = vert20[2];
          p[120] = vert2t1[0];
          p[121] = vert2t1[1];
          p[122] = vert2t1[2];
          p[123] = vert2tet[0];
          p[124] = vert2tet[1];
          p[125] = vert2tet[2];
          p[126] = vert21[0];
          p[127] = vert21[1];
          p[128] = vert21[2];
          p[129] = vert2t0[0];
          p[130] = vert2t0[1];
          p[131] = vert2t0[2];
          p[132] = vert2tet[0];
          p[133] = vert2tet[1];
          p[134] = vert2tet[2];
          p[135] = vert21[0];
          p[136] = vert21[1];
          p[137] = vert21[2];
          p[138] = vert2t2[0];
          p[139] = vert2t2[1];
          p[140] = vert2t2[2];
          p[141] = vert2tet[0];
          p[142] = vert2tet[1];
          p[143] = vert2tet[2];
          p[144] = vert22[0];
          p[145] = vert22[1];
          p[146] = vert22[2];
          p[147] = vert2t2[0];
          p[148] = vert2t2[1];
          p[149] = vert2t2[2];
          p[150] = vert2tet[0];
          p[151] = vert2tet[1];
          p[152] = vert2tet[2];
          p[153] = vert22[0];
          p[154] = vert22[1];
          p[155] = vert22[2];
          p[156] = vert2t1[0];
          p[157] = vert2t1[1];
          p[158] = vert2t1[2];
          p[159] = vert2tet[0];
          p[160] = vert2tet[1];
          p[161] = vert2tet[2];

          // Label Vert 3
          p[162] = vert30[0];
          p[163] = vert30[1];
          p[164] = vert30[2];
          p[165] = vert3t0[0];
          p[166] = vert3t0[1];
          p[167] = vert3t0[2];
          p[168] = vert3tet[0];
          p[169] = vert3tet[1];
          p[170] = vert3tet[2];
          p[171] = vert30[0];
          p[172] = vert30[1];
          p[173] = vert30[2];
          p[174] = vert3t1[0];
          p[175] = vert3t1[1];
          p[176] = vert3t1[2];
          p[177] = vert3tet[0];
          p[178] = vert3tet[1];
          p[179] = vert3tet[2];
          p[180] = vert31[0];
          p[181] = vert31[1];
          p[182] = vert31[2];
          p[183] = vert3t0[0];
          p[184] = vert3t0[1];
          p[185] = vert3t0[2];
          p[186] = vert3tet[0];
          p[187] = vert3tet[1];
          p[188] = vert3tet[2];
          p[189] = vert31[0];
          p[190] = vert31[1];
          p[191] = vert31[2];
          p[192] = vert3t2[0];
          p[193] = vert3t2[1];
          p[194] = vert3t2[2];
          p[195] = vert3tet[0];
          p[196] = vert3tet[1];
          p[197] = vert3tet[2];
          p[198] = vert32[0];
          p[199] = vert32[1];
          p[200] = vert32[2];
          p[201] = vert3t2[0];
          p[202] = vert3t2[1];
          p[203] = vert3t2[2];
          p[204] = vert3tet[0];
          p[205] = vert3tet[1];
          p[206] = vert3tet[2];
          p[207] = vert32[0];
          p[208] = vert32[1];
          p[209] = vert32[2];
          p[210] = vert3t1[0];
          p[211] = vert3t1[1];
          p[212] = vert3t1[2];
          p[213] = vert3tet[0];
          p[214] = vert3tet[1];
          p[215] = vert3tet[2];
          p += 216;

          c[0] = numThreadIndex + 0;
          c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2;
          c[3] = numThreadIndex + 3;
          c[4] = numThreadIndex + 4;
          c[5] = numThreadIndex + 5;
          c[6] = numThreadIndex + 6;
          c[7] = numThreadIndex + 7;
          c[8] = numThreadIndex + 8;
          c[9] = numThreadIndex + 9;
          c[10] = numThreadIndex + 10;
          c[11] = numThreadIndex + 11;
          c[12] = numThreadIndex + 12;
          c[13] = numThreadIndex + 13;
          c[14] = numThreadIndex + 14;
          c[15] = numThreadIndex + 15;
          c[16] = numThreadIndex + 16;
          c[17] = numThreadIndex + 17;
          c[18] = numThreadIndex + 18;
          c[19] = numThreadIndex + 19;
          c[20] = numThreadIndex + 20;
          c[21] = numThreadIndex + 21;
          c[22] = numThreadIndex + 22;
          c[23] = numThreadIndex + 23;
          c[24] = numThreadIndex + 24;
          c[25] = numThreadIndex + 25;
          c[26] = numThreadIndex + 26;
          c[27] = numThreadIndex + 27;
          c[28] = numThreadIndex + 28;
          c[29] = numThreadIndex + 29;
          c[30] = numThreadIndex + 30;
          c[31] = numThreadIndex + 31;
          c[32] = numThreadIndex + 32;
          c[33] = numThreadIndex + 33;
          c[34] = numThreadIndex + 34;
          c[35] = numThreadIndex + 35;
          c[36] = numThreadIndex + 36;
          c[37] = numThreadIndex + 37;
          c[38] = numThreadIndex + 38;
          c[39] = numThreadIndex + 39;
          c[40] = numThreadIndex + 40;
          c[41] = numThreadIndex + 41;
          c[42] = numThreadIndex + 42;
          c[43] = numThreadIndex + 43;
          c[44] = numThreadIndex + 44;
          c[45] = numThreadIndex + 45;
          c[46] = numThreadIndex + 46;
          c[47] = numThreadIndex + 47;
          c[48] = numThreadIndex + 48;
          c[49] = numThreadIndex + 49;
          c[50] = numThreadIndex + 50;
          c[51] = numThreadIndex + 51;
          c[52] = numThreadIndex + 52;
          c[53] = numThreadIndex + 53;
          c[54] = numThreadIndex + 54;
          c[55] = numThreadIndex + 55;
          c[56] = numThreadIndex + 56;
          c[57] = numThreadIndex + 57;
          c[58] = numThreadIndex + 58;
          c[59] = numThreadIndex + 59;
          c[60] = numThreadIndex + 60;
          c[61] = numThreadIndex + 61;
          c[62] = numThreadIndex + 62;
          c[63] = numThreadIndex + 63;
          c[64] = numThreadIndex + 64;
          c[65] = numThreadIndex + 65;
          c[66] = numThreadIndex + 66;
          c[67] = numThreadIndex + 67;
          c[68] = numThreadIndex + 68;
          c[69] = numThreadIndex + 69;
          c[70] = numThreadIndex + 70;
          c[71] = numThreadIndex + 71;
          c += 72;
          numThreadIndex += 72;

          m[0] = msm[0];
          m[1] = msm[0];
          m[2] = msm[0];
          m[3] = msm[0];
          m[4] = msm[0];
          m[5] = msm[0];
          m[6] = msm[1];
          m[7] = msm[1];
          m[8] = msm[1];
          m[9] = msm[1];
          m[10] = msm[1];
          m[11] = msm[1];
          m[12] = msm[2];
          m[13] = msm[2];
          m[14] = msm[2];
          m[15] = msm[2];
          m[16] = msm[2];
          m[17] = msm[2];
          m[18] = msm[3];
          m[19] = msm[3];
          m[20] = msm[3];
          m[21] = msm[3];
          m[22] = msm[3];
          m[23] = msm[3];
          m += 24;
        }
      }
    }
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Wrote detailed boundaries", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}