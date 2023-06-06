/// \ingroup base
/// \author Robin G. C. Maack <maack@rptu.de>
/// \date May 2023.
/// \brief Lookuptables used in MarchingTetrahedra algorithm

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
