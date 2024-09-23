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

#include <climits>
#include <queue>
#include <type_traits>

#include <MarchingTetrahedraLookupTables.inl>

using ttk::SimplexId;

namespace ttk {
  namespace mth {
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
             % ULLONG_MAX;
    }

    /**
     * @brief Get the center of two points
     *
     * @param[in] pos0 Position 0
     * @param[in] pos1 Position 1
     * @param[out] incenter Resulting position
     * @return int 0 on success
     */
    inline void getCenter(const std::array<float, 3> pos0,
                          const std::array<float, 3> pos1,
                          std::array<float, 3> &incenter) {
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
    inline void getCenter(const std::array<float, 3> pos0,
                          const std::array<float, 3> pos1,
                          const std::array<float, 3> pos2,
                          std::array<float, 3> &incenter) {
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
    inline void getCenter(const std::array<float, 3> pos0,
                          const std::array<float, 3> pos1,
                          const std::array<float, 3> pos2,
                          const std::array<float, 3> pos3,
                          std::array<float, 3> &incenter) {
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
    inline void interpolatePoints(const std::array<float, 3> pos0,
                                  const std::array<float, 3> pos1,
                                  const float lambda,
                                  std::array<float, 3> &result) {

      result[0] = lambda * pos0[0] + (1 - lambda) * pos1[0];
      result[1] = lambda * pos0[1] + (1 - lambda) * pos1[1];
      result[2] = lambda * pos0[2] + (1 - lambda) * pos1[2];
    }
  } // namespace mth

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
     * @tparam triangulationType Triangulationtype
     * @param[out] tetCases Binary codes
     * @param[out] numEdges Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangleCounter Table with binary code to number of triangles
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <typename triangulationType>
    int computeMarchingCases_2D(unsigned char *const tetCases,
                                size_t *const numEdges,
                                const unsigned long long *const scalars,
                                const size_t *const triangleCounter,
                                const triangulationType &triangulation) const;

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 2D basin separators.
     *
     * @tparam triangulationType Triangulationtype
     * @param[in] tetCases Binary codes
     * @param[in] numEdges Number of edges generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <typename triangulationType>
    int writeSeparators_2D(const unsigned char *const tetCases,
                           const size_t *numEdges,
                           const unsigned long long *const scalars,
                           const triangulationType &triangulation);

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 2D basin boundaries.
     *
     * @tparam triangulationType Triangulationtype
     * @param[in] tetCases Binary codes
     * @param[in] numEdges Number of edges generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <typename triangulationType>
    int writeBoundaries_2D(const unsigned char *const tetCases,
                           const size_t *numEdges,
                           const unsigned long long *const scalars,
                           const triangulationType &triangulation);

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 2D detailed basin boundaries.
     *
     * @tparam triangulationType Triangulationtype
     * @param[in] tetCases Binary codes
     * @param[in] numEdges Number of edges generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <typename triangulationType>
    int writeBoundariesDetailed_2D(const unsigned char *const tetCases,
                                   const size_t *numEdges,
                                   const unsigned long long *const scalars,
                                   const triangulationType &triangulation);

    /* 3D Datasets */

    /**
     * Computes a binary code for each tetrahedron and counts the number of
     * triangles generated, depending on the triangleCounter used.
     *
     * @tparam triangulationType Triangulation
     * @param[out] tetCases Binary codes
     * @param[out] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangleCounter Table with binary code to number of triangles
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <typename triangulationType>
    int computeMarchingCases_3D(unsigned char *const tetCases,
                                size_t *const numTriangles,
                                const unsigned long long *const scalars,
                                const size_t *const triangleCounter,
                                const triangulationType &triangulation) const;

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 3D basin separators.
     *
     * @tparam triangulationType Triangulationtype
     * @param[in] tetCases Binary codes
     * @param[in] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <typename triangulationType>
    int writeSeparators_3D(const unsigned char *const tetCases,
                           const size_t *numTriangles,
                           const unsigned long long *const scalars,
                           const triangulationType &triangulation);

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 3D basin boundaries.
     *
     * @tparam triangulationType Triangulationtype
     * @param[in] tetCases Binary codes
     * @param[in] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <typename triangulationType>
    int writeBoundaries_3D(const unsigned char *const tetCases,
                           const size_t *numTriangles,
                           const unsigned long long *const scalars,
                           const triangulationType &triangulation);

    /**
     * Writes the geometry of all Triangles to the ouputVariables set by
     * setOutput. 3D detailed basin boundaries.
     *
     * @tparam triangulationType Triangulationtype
     * @param[in] tetCases Binary codes
     * @param[in] numTriangles Number of triangles generated per thread
     * @param[in] scalars Scalars
     * @param[in] triangulation Triangulation
     * @return int
     */
    template <typename triangulationType>
    int writeBoundariesDetailed_3D(const unsigned char *const tetCases,
                                   const size_t *numTriangles,
                                   const unsigned long long *const scalars,
                                   const triangulationType &triangulation);

  protected:
    // Output options
    SURFACE_MODE SurfaceMode{SURFACE_MODE::SM_SEPARATORS};

    // Output data
    SimplexId output_numberOfPoints_{};
    SimplexId output_numberOfCells_{};
    std::vector<float> output_points_;
    std::vector<unsigned long long> output_cells_labels_;
    std::vector<SimplexId> output_cells_connectivity_;
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

  std::vector<unsigned char> tetCases;
  std::vector<size_t> numberOfTetCases;
  std::vector<unsigned long long> cScalars;

  cScalars.resize(nV);
  tetCases.resize(nC);

#ifdef TTK_ENABLE_OPENMP
  numberOfTetCases.resize(this->threadNumber_);
#else
  numberOfTetCases.resize(1);
#endif // TTK_ENABLE_OPENMP

  for(SimplexId vert = 0; vert < nV; vert++)
    std::memcpy(&cScalars[vert], &scalars[vert], sizeof(dataType));

  if(dim == 2) {
    if(SurfaceMode == SURFACE_MODE::SM_SEPARATORS) {
      computeMarchingCases_2D(&tetCases[0], &numberOfTetCases[0], &cScalars[0],
                              triangleNumberLookup, triangulation);
      writeSeparators_2D(
        &tetCases[0], &numberOfTetCases[0], &cScalars[0], triangulation);
    } else if(SurfaceMode == SURFACE_MODE::SM_BOUNDARIES) {
      computeMarchingCases_2D(&tetCases[0], &numberOfTetCases[0], &cScalars[0],
                              triangleNumberLookupBoundary, triangulation);
      writeBoundaries_2D(
        &tetCases[0], &numberOfTetCases[0], &cScalars[0], triangulation);
    } else if(SurfaceMode == SURFACE_MODE::SM_BOUNDARIES_DETAILED) {
      computeMarchingCases_2D(&tetCases[0], &numberOfTetCases[0], &cScalars[0],
                              triangleNumberLookupBoundaryDetailed,
                              triangulation);
      writeBoundariesDetailed_2D(
        &tetCases[0], &numberOfTetCases[0], &cScalars[0], triangulation);
    }
  } else if(dim == 3) {
    if(SurfaceMode == SURFACE_MODE::SM_SEPARATORS) {
      computeMarchingCases_3D(&tetCases[0], &numberOfTetCases[0], &cScalars[0],
                              tetLookupNumWallTriangles, triangulation);
      writeSeparators_3D(
        &tetCases[0], &numberOfTetCases[0], &cScalars[0], triangulation);
    } else if(SurfaceMode == SURFACE_MODE::SM_BOUNDARIES) {
      computeMarchingCases_3D(&tetCases[0], &numberOfTetCases[0], &cScalars[0],
                              tetLookupNumTrianglesBoundaries, triangulation);
      writeBoundaries_3D(
        &tetCases[0], &numberOfTetCases[0], &cScalars[0], triangulation);
    } else if(SurfaceMode == SURFACE_MODE::SM_BOUNDARIES_DETAILED) {
      computeMarchingCases_3D(&tetCases[0], &numberOfTetCases[0], &cScalars[0],
                              tetLookupNumTrianglesDetailedBoundary,
                              triangulation);
      writeBoundariesDetailed_3D(
        &tetCases[0], &numberOfTetCases[0], &cScalars[0], triangulation);
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

template <typename triangulationType>
int ttk::MarchingTetrahedra::computeMarchingCases_2D(
  unsigned char *const tetCases,
  size_t *const numEdges,
  const unsigned long long *const scalars,
  const size_t *const triangleCounter,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing separator cases", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  const SimplexId numTetra = triangulation.getNumberOfCells();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
  {
    SimplexId threadEdges = 0;
    const int tid = omp_get_thread_num();

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  SimplexId threadEdges = 0;
  const int tid = 0;
#endif // TTK_ENABLE_OPENMP
    for(SimplexId tet = 0; tet < numTetra; ++tet) {
      std::array<SimplexId, 3> vertices{};
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);

      const std::array<unsigned long long, 3> label
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

      // Set the third bit to 0 or 1
      const unsigned char index0 = (label[0] == label[1]) ? 0x00 : 0x04;

      // Set the first and second bit to 0, 1 or 2
      const unsigned char index1 = (label[0] == label[2])   ? 0x00
                                   : (label[1] == label[2]) ? 0x01
                                                            : 0x02;

      tetCases[tet] = index0 | index1;
      threadEdges += triangleCounter[tetCases[tet]];
    }
    numEdges[tid] = threadEdges;
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed separator cases", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MarchingTetrahedra::writeSeparators_2D(
  const unsigned char *const tetCases,
  const size_t *numEdges,
  const unsigned long long *const scalars,
  const triangulationType &triangulation) {

  ttk::Timer localTimer;

  this->printMsg("Writing Boundaries", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
  std::vector<size_t> edgeStartIndex(this->threadNumber_ + 1);

  edgeStartIndex[0] = 0;

  // Count triangle number and create iterator start indices
  for(int t = 1; t <= this->threadNumber_; ++t) {
    edgeStartIndex[t] = numEdges[t - 1] + edgeStartIndex[t - 1];
  }

  const size_t numTotalEdges = edgeStartIndex[this->threadNumber_];
#else // TTK_ENABLE_OPENMP
  const size_t numTotalEdges = numEdges[0];
#endif // TTK_ENABLE_OPENMP

  output_points_.resize(6 * numTotalEdges);
  output_cells_connectivity_.resize(2 * numTotalEdges);
  output_cells_labels_.resize(numTotalEdges);
  output_numberOfPoints_ = 2 * numTotalEdges;
  output_numberOfCells_ = numTotalEdges;

  const SimplexId numTets = triangulation.getNumberOfCells();
  float *p = output_points_.data();
  SimplexId *c = output_cells_connectivity_.data();
  unsigned long long *m = output_cells_labels_.data();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) firstprivate(p, c, m)
  {
    const int tid = omp_get_thread_num();
    size_t numThreadIndex = edgeStartIndex[tid];

    p += (numThreadIndex * 6);
    c += (numThreadIndex * 2);
    m += numThreadIndex;
    numThreadIndex = 2 * numThreadIndex;

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  size_t numThreadIndex = 0;
#endif // TTK_ENABLE_OPENMP
    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!triangleLookupIsMultiLabel[tetCases[tet]]) {
        continue;
      }

      const int *edgeVerts = triangleLookupEdgeVerts[tetCases[tet]];

      std::array<SimplexId, 3> vertices{};
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);

      std::array<std::array<float, 3>, 3> vertPos{};
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);

      const std::array<unsigned long long, 3> label
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

      if(triangleLookupIs2Label[tetCases[tet]]) {
        std::array<std::array<float, 3>, 2> eC{};
        mth::getCenter(vertPos[edgeVerts[0]], vertPos[edgeVerts[1]], eC[0]);
        mth::getCenter(vertPos[edgeVerts[2]], vertPos[edgeVerts[3]], eC[1]);

        // Create a hash from vertex label 0 and 1
        const unsigned long long sparseID
          = mth::getHash(label[edgeVerts[0]], label[edgeVerts[1]]);

        // Write the edge endpoints, cell indices, and label
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
        std::array<std::array<float, 3>, 4> eC{};
        mth::getCenter(vertPos[0], vertPos[1], eC[0]);
        mth::getCenter(vertPos[0], vertPos[2], eC[1]);
        mth::getCenter(vertPos[1], vertPos[2], eC[2]);
        mth::getCenter(vertPos[0], vertPos[1], vertPos[2], eC[3]);

        // Create a hash from all vertex label combinations
        const std::array<unsigned long long, 3> sparseID
          = {mth::getHash(label[0], label[1]), mth::getHash(label[0], label[2]),
             mth::getHash(label[1], label[2])};

        // Write all three lines, each connected to the triangle center
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
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg(
    "Wrote Boundaries", 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MarchingTetrahedra::writeBoundaries_2D(
  const unsigned char *const tetCases,
  const size_t *numEdges,
  const unsigned long long *const scalars,
  const triangulationType &triangulation) {

  ttk::Timer localTimer;

  this->printMsg("Writing Boundaries", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
  std::vector<size_t> edgeStartIndex(this->threadNumber_ + 1);

  edgeStartIndex[0] = 0;

  // Count triangle number and create iterator start indices
  for(int t = 1; t <= this->threadNumber_; ++t) {
    edgeStartIndex[t] = numEdges[t - 1] + edgeStartIndex[t - 1];
  }

  size_t const numTotalEdges = edgeStartIndex[this->threadNumber_];
#else // TTK_ENABLE_OPENMP
  size_t numTotalEdges = numEdges[0];
#endif // TTK_ENABLE_OPENMP

  output_points_.resize(6 * numTotalEdges);
  output_cells_connectivity_.resize(2 * numTotalEdges);
  output_cells_labels_.resize(numTotalEdges);
  output_numberOfPoints_ = 2 * numTotalEdges;
  output_numberOfCells_ = numTotalEdges;

  float *p = output_points_.data();
  SimplexId *c = output_cells_connectivity_.data();
  unsigned long long *m = output_cells_labels_.data();

  const SimplexId numTets = triangulation.getNumberOfCells();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) firstprivate(p, c, m)
  {
    const int tid = omp_get_thread_num();
    size_t numThreadIndex = edgeStartIndex[tid];

    p += (numThreadIndex * 6);
    c += (numThreadIndex * 2);
    m += numThreadIndex;

    numThreadIndex = 2 * numThreadIndex;

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  size_t numThreadIndex = 0;
#endif // TTK_ENABLE_OPENMP
    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!triangleLookupIsMultiLabel[tetCases[tet]]) {
        continue;
      }

      std::array<SimplexId, 3> vertices{};
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);

      std::array<std::array<float, 3>, 3> vertPos{};
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);

      const std::array<unsigned long long, 3> label
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

      if(triangleLookupIs2Label[tetCases[tet]]) {
        const int *edgeVerts = triangleLookupEdgeVerts[tetCases[tet]];

        // Write the edge endpoints, cell indices, and label
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

        m[0] = label[edgeVerts[0]];
        m += 1;
      }
    }
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg(
    "Wrote Boundaries", 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MarchingTetrahedra::writeBoundariesDetailed_2D(
  const unsigned char *const tetCases,
  const size_t *numEdges,
  const unsigned long long *const scalars,
  const triangulationType &triangulation) {

  ttk::Timer localTimer;

  this->printMsg("Writing Boundaries", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  constexpr float diff = 0.02;
  constexpr float d0 = 0.5 + diff;
  constexpr float d1 = 0.5 - diff;
  constexpr float dc = 2 * diff;

#ifdef TTK_ENABLE_OPENMP
  std::vector<size_t> edgeStartIndex(this->threadNumber_ + 1);

  edgeStartIndex[0] = 0;

  // Count triangle numbers and create iterator start indices
  for(int t = 1; t <= this->threadNumber_; ++t) {
    edgeStartIndex[t] = numEdges[t - 1] + edgeStartIndex[t - 1];
  }

  size_t const numTotalEdges = edgeStartIndex[this->threadNumber_];
#else // TTK_ENABLE_OPENMP
  size_t numTotalEdges = numEdges[0];
#endif // TTK_ENABLE_OPENMP

  output_points_.resize(6 * numTotalEdges);
  output_cells_connectivity_.resize(2 * numTotalEdges);
  output_cells_labels_.resize(numTotalEdges);
  output_numberOfPoints_ = 2 * numTotalEdges;
  output_numberOfCells_ = numTotalEdges;

  float *p = output_points_.data();
  SimplexId *c = output_cells_connectivity_.data();
  unsigned long long *m = output_cells_labels_.data();

  const SimplexId numTets = triangulation.getNumberOfCells();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) firstprivate(p, c, m)
  {
    const int tid = omp_get_thread_num();
    size_t numThreadIndex = edgeStartIndex[tid];

    p += (numThreadIndex * 6);
    c += (numThreadIndex * 2);
    m += numThreadIndex;

    numThreadIndex = 2 * numThreadIndex;

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  size_t numThreadIndex = 0;
#endif // TTK_ENABLE_OPENMP
    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!triangleLookupIsMultiLabel[tetCases[tet]]) {
        continue;
      }

      const int *vIds = triangleLookupEdgeVerts[tetCases[tet]];

      std::array<SimplexId, 3> vertices{};
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);

      std::array<std::array<float, 3>, 3> vPos{};
      triangulation.getVertexPoint(
        vertices[0], vPos[0][0], vPos[0][1], vPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vPos[1][0], vPos[1][1], vPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vPos[2][0], vPos[2][1], vPos[2][2]);

      const std::array<unsigned long long, 3> label
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]]};

      if(triangleLookupIs2Label[tetCases[tet]]) {
        std::array<float, 3> vert00{}, vert01{}, vert10{}, vert11{};

        mth::interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
        mth::interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
        mth::interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d1, vert10);
        mth::interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d1, vert11);

        // Write the edge endpoints, cell indices, and labels
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

        m[0] = label[vIds[0]];
        m[1] = label[vIds[1]];
        m += 2;

      } else {
        std::array<float, 3> vert00{}, vert01{}, vert10{}, vert11{}, vert20{},
          vert21{};
        mth::interpolatePoints(vPos[0], vPos[1], d0, vert00);
        mth::interpolatePoints(vPos[0], vPos[2], d0, vert01);
        mth::interpolatePoints(vPos[1], vPos[0], d0, vert10);
        mth::interpolatePoints(vPos[1], vPos[2], d0, vert11);
        mth::interpolatePoints(vPos[2], vPos[0], d0, vert20);
        mth::interpolatePoints(vPos[2], vPos[1], d0, vert21);

        std::array<float, 3> triCenter{};
        mth::getCenter(vPos[0], vPos[1], vPos[2], triCenter);

        std::array<float, 3> triS0{}, triS1{}, triS2{};
        mth::interpolatePoints(vPos[0], triCenter, dc, triS0);
        mth::interpolatePoints(vPos[1], triCenter, dc, triS1);
        mth::interpolatePoints(vPos[2], triCenter, dc, triS2);

        // Write the edge endpoints, cell indices, and labels
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

        m[0] = label[0];
        m[1] = label[0];
        m[2] = label[1];
        m[3] = label[1];
        m[4] = label[2];
        m[5] = label[2];
        m += 6;
      }
    }
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg(
    "Wrote Boundaries", 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MarchingTetrahedra::computeMarchingCases_3D(
  unsigned char *const tetCases,
  size_t *const numTriangles,
  const unsigned long long *const scalars,
  const size_t *const triangleCounter,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing separator cases", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  const SimplexId numTetra = triangulation.getNumberOfCells();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
  {
    SimplexId threadTriangles = 0;
    const int tid = omp_get_thread_num();

#pragma omp for schedule(static)
#else
  SimplexId threadTriangles = 0;
  const int tid = 0;
#endif
    for(SimplexId tet = 0; tet < numTetra; ++tet) {

      std::array<SimplexId, 4> vertices{};
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const std::array<unsigned long long, 4> label
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
           scalars[vertices[3]]};

      // Set the fifth bit to 0 or 1
      const unsigned char index1 = (label[0] == label[1]) ? 0x00 : 0x10;

      // Set the fourth and third bit to 0, 1 or 2
      const unsigned char index2 = (label[0] == label[2])   ? 0x00
                                   : (label[1] == label[2]) ? 0x04
                                                            : 0x08;

      // Set the first and second bit to 0, 1, 2 or 3
      const unsigned char index3 = (label[0] == label[3])   ? 0x00
                                   : (label[1] == label[3]) ? 0x01
                                   : (label[2] == label[3]) ? 0x02
                                                            : 0x03;

      tetCases[tet] = index1 | index2 | index3;
      threadTriangles += triangleCounter[tetCases[tet]];
    }
    numTriangles[tid] = threadTriangles;
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed separator cases", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MarchingTetrahedra::writeSeparators_3D(
  const unsigned char *const tetCases,
  const size_t *numTriangles,
  const unsigned long long *const scalars,
  const triangulationType &triangulation) {

  ttk::Timer localTimer;

  this->printMsg("Writing separators", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
  std::vector<size_t> triangleStartIndex(this->threadNumber_ + 1);

  triangleStartIndex[0] = 0;

  // Count triangle number and create iterator start indices
  for(int t = 1; t <= this->threadNumber_; ++t) {
    triangleStartIndex[t] = numTriangles[t - 1] + triangleStartIndex[t - 1];
  }

  size_t const numTotalTriangles = triangleStartIndex[this->threadNumber_];
#else // TTK_ENABLE_OPENMP
  size_t numTotalTriangles = numTriangles[0];
#endif // TTK_ENABLE_OPENMP

  output_points_.resize(9 * numTotalTriangles);
  output_cells_connectivity_.resize(3 * numTotalTriangles);
  output_cells_labels_.resize(numTotalTriangles);
  output_numberOfPoints_ = 3 * numTotalTriangles;
  output_numberOfCells_ = numTotalTriangles;

  float *p = output_points_.data();
  SimplexId *c = output_cells_connectivity_.data();
  unsigned long long *m = output_cells_labels_.data();

  const SimplexId numTets = triangulation.getNumberOfCells();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) firstprivate(p, c, m)
  {
    const int tid = omp_get_thread_num();
    size_t numThreadIndex = triangleStartIndex[tid];

    p += (numThreadIndex * 9);
    c += (numThreadIndex * 3);
    m += numThreadIndex;

    numThreadIndex = 3 * numThreadIndex;

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  size_t numThreadIndex = 0;
#endif // TTK_ENABLE_OPENMP
    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!tetLookupIsMultiLabel[tetCases[tet]]) {
        continue;
      }

      const int *tetEdgeIndices = tetLookupWall[tetCases[tet]];
      const int *tetVertLabel = tetLookupWallLabel[tetCases[tet]];

      std::array<SimplexId, 4> vertices{};
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const std::array<unsigned long long, 4> label
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
           scalars[vertices[3]]};

      std::array<std::array<float, 3>, 4> vertPos{};
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      std::array<std::array<float, 3>, 10> eC{};
      // 6 edge centers
      mth::getCenter(vertPos[0], vertPos[1], eC[0]);
      mth::getCenter(vertPos[0], vertPos[2], eC[1]);
      mth::getCenter(vertPos[0], vertPos[3], eC[2]);
      mth::getCenter(vertPos[1], vertPos[2], eC[3]);
      mth::getCenter(vertPos[1], vertPos[3], eC[4]);
      mth::getCenter(vertPos[2], vertPos[3], eC[5]);

      // 4 triangle centers
      mth::getCenter(vertPos[0], vertPos[1], vertPos[2], eC[6]);
      mth::getCenter(vertPos[0], vertPos[1], vertPos[3], eC[7]);
      mth::getCenter(vertPos[0], vertPos[2], vertPos[3], eC[8]);
      mth::getCenter(vertPos[1], vertPos[2], vertPos[3], eC[9]);

      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        std::array<float, 3> tetCenter{};
        mth::getCenter(
          vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

        // Create a hashes from all four label combinations
        unsigned long long const sparseMSIds[6] = {
          mth::getHash(label[0], label[1]), mth::getHash(label[0], label[2]),
          mth::getHash(label[0], label[3]), mth::getHash(label[1], label[2]),
          mth::getHash(label[1], label[3]), mth::getHash(label[2], label[3])};

        // Write the triangle endpoints, cell indices, and label
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

        for(size_t t = 0; t < numTris; ++t) {
          // Write the triangle endpoints, cell indices, and label
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

          // Create hash from the corresponsing vertex labels
          m[0] = mth::getHash(
            label[tetVertLabel[t * 2]], label[tetVertLabel[(t * 2) + 1]]);
          m += 1;
        }
      }
    }
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg(
    "Wrote separators", 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MarchingTetrahedra::writeBoundaries_3D(
  const unsigned char *const tetCases,
  const size_t *numTriangles,
  const unsigned long long *const scalars,
  const triangulationType &triangulation) {

  ttk::Timer localTimer;

  this->printMsg("Writing Boundaries", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
  std::vector<size_t> triangleStartIndex(this->threadNumber_ + 1);

  triangleStartIndex[0] = 0;

  // Count triangle number and create iterator start indices
  for(int t = 1; t <= this->threadNumber_; ++t) {
    triangleStartIndex[t] = numTriangles[t - 1] + triangleStartIndex[t - 1];
  }

  size_t const numTotalTriangles = triangleStartIndex[this->threadNumber_];
#else // TTK_ENABLE_OPENMP
  size_t numTotalTriangles = numTriangles[0];
#endif // TTK_ENABLE_OPENMP

  output_points_.resize(9 * numTotalTriangles);
  output_cells_connectivity_.resize(3 * numTotalTriangles);
  output_cells_labels_.resize(numTotalTriangles);
  output_numberOfPoints_ = 3 * numTotalTriangles;
  output_numberOfCells_ = numTotalTriangles;

  float *p = output_points_.data();
  SimplexId *c = output_cells_connectivity_.data();
  unsigned long long *m = output_cells_labels_.data();

  const SimplexId numTets = triangulation.getNumberOfCells();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) firstprivate(p, c, m)
  {
    const int tid = omp_get_thread_num();
    size_t numThreadIndex = triangleStartIndex[tid];

    p += (numThreadIndex * 9);
    c += (numThreadIndex * 3);
    m += numThreadIndex;

    numThreadIndex = 3 * numThreadIndex;

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  size_t numThreadIndex = 0;
#endif // TTK_ENABLE_OPENMP
    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!tetLookupFast[tetCases[tet]]) {
        continue;
      }

      std::array<SimplexId, 4> vertices{};
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const std::array<unsigned long long, 4> label
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
           scalars[vertices[3]]};

      // Get the tetrahedron local vertex ids that all three have the same
      // label
      const int id0 = tetLookupFastTri[tetLookupFastCase[tetCases[tet]]][0];
      const int id1 = tetLookupFastTri[tetLookupFastCase[tetCases[tet]]][1];
      const int id2 = tetLookupFastTri[tetLookupFastCase[tetCases[tet]]][2];

      std::array<std::array<float, 3>, 4> vertPos{};
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      // Write the triangle endpoints, cell indices, and label
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

      m[0] = label[id0];
      m += 1;
    }
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg(
    "Wrote Boundaries", 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MarchingTetrahedra::writeBoundariesDetailed_3D(
  const unsigned char *const tetCases,
  const size_t *numTriangles,
  const unsigned long long *const scalars,
  const triangulationType &triangulation) {

  ttk::Timer localTimer;

  this->printMsg("Writing detailed boundaries", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  constexpr float diff = 0.02;
  constexpr float d0 = 0.5 + diff;
  constexpr float d1 = 0.5 - diff;
  constexpr float dc = diff * 2;

#ifdef TTK_ENABLE_OPENMP
  std::vector<size_t> triangleStartIndex(this->threadNumber_ + 1);

  triangleStartIndex[0] = 0;

  // Count triangle number and create iterator start indices
  for(int t = 1; t <= this->threadNumber_; ++t) {
    triangleStartIndex[t] = numTriangles[t - 1] + triangleStartIndex[t - 1];
  }

  size_t const numTotalTriangles = triangleStartIndex[this->threadNumber_];
#else // TTK_ENABLE_OPENMP
  size_t numTotalTriangles = numTriangles[0];
#endif // TTK_ENABLE_OPENMP

  output_points_.resize(9 * numTotalTriangles);
  output_cells_connectivity_.resize(3 * numTotalTriangles);
  output_cells_labels_.resize(numTotalTriangles);
  output_numberOfPoints_ = 3 * numTotalTriangles;
  output_numberOfCells_ = numTotalTriangles;

  float *p = output_points_.data();
  SimplexId *c = output_cells_connectivity_.data();
  unsigned long long *m = output_cells_labels_.data();

  const SimplexId numTets = triangulation.getNumberOfCells();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) firstprivate(p, c, m)
  {
    const int tid = omp_get_thread_num();
    size_t numThreadIndex = triangleStartIndex[tid];

    p += (numThreadIndex * 9);
    c += (numThreadIndex * 3);
    m += numThreadIndex;

    numThreadIndex = 3 * numThreadIndex;

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  size_t numThreadIndex = 0;
#endif // TTK_ENABLE_OPENMP
    for(SimplexId tet = 0; tet < numTets; ++tet) {
      if(!tetLookupIsMultiLabel[tetCases[tet]]) {
        continue;
      }

      std::array<SimplexId, 4> vertices{};
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const std::array<unsigned long long, 4> label
        = {scalars[vertices[0]], scalars[vertices[1]], scalars[vertices[2]],
           scalars[vertices[3]]};

      std::array<std::array<float, 3>, 4> vPos{};
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

        std::array<float, 3> vert00{}, vert01{}, vert02{}, vert10{}, vert11{},
          vert12{};

        mth::interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
        mth::interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
        mth::interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d0, vert02);
        mth::interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d1, vert10);
        mth::interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d1, vert11);
        mth::interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d1, vert12);

        // Write the triangle endpoints, cell indices, and label
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

        m[0] = label[vIds[0]];
        m[1] = label[vIds[1]];
        m += 2;

        if(vIds[6] != -1) { // 2 vertices per label (e.g. AABB)
          std::array<float, 3> vert03{}, vert13{};

          mth::interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d0, vert03);
          mth::interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d1, vert13);

          // Write the triangle endpoints, cell indices, and label
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

          m[0] = label[vIds[0]];
          m[1] = label[vIds[1]];
          m += 2;
        }
      } else if(tetLookupIs3Label[tetCases[tet]]) {
        const int *vIds = tetLookupSplitBasisns3Label[tetCases[tet]];

        std::array<std::array<float, 3>, 4> triCenter{};
        mth::getCenter(vPos[0], vPos[1], vPos[2], triCenter[0]);
        mth::getCenter(vPos[0], vPos[1], vPos[3], triCenter[1]);
        mth::getCenter(vPos[0], vPos[2], vPos[3], triCenter[2]);
        mth::getCenter(vPos[1], vPos[2], vPos[3], triCenter[3]);

        std::array<std::array<float, 3>, 10> edgeCenters{};

        mth::getCenter(vPos[0], vPos[1], edgeCenters[0]);
        mth::getCenter(vPos[0], vPos[2], edgeCenters[1]);
        mth::getCenter(vPos[0], vPos[3], edgeCenters[2]);
        mth::getCenter(vPos[1], vPos[2], edgeCenters[3]);
        mth::getCenter(vPos[1], vPos[3], edgeCenters[4]);
        mth::getCenter(vPos[2], vPos[3], edgeCenters[5]);

        std::array<float, 3> edge00{}, edge01{}, edge02{}, edge03{}, tri00{},
          tri01{}, edge10{}, edge11{}, edge12{}, tri10{}, tri11{}, edge20{},
          edge21{}, edge22{}, tri20{}, tri21{};

        mth::interpolatePoints(vPos[vIds[0]], edgeCenters[vIds[1]], dc, edge00);
        mth::interpolatePoints(vPos[vIds[0]], edgeCenters[vIds[2]], dc, edge01);
        mth::interpolatePoints(vPos[vIds[0]], triCenter[vIds[3]], dc, tri00);
        mth::interpolatePoints(vPos[vIds[4]], edgeCenters[vIds[5]], dc, edge02);
        mth::interpolatePoints(vPos[vIds[4]], edgeCenters[vIds[6]], dc, edge03);
        mth::interpolatePoints(vPos[vIds[4]], triCenter[vIds[7]], dc, tri01);

        mth::interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[1]], dc, edge10);
        mth::interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[5]], dc, edge11);
        mth::interpolatePoints(
          vPos[vIds[8]], edgeCenters[vIds[10]], dc, edge12);
        mth::interpolatePoints(vPos[vIds[8]], triCenter[vIds[3]], dc, tri10);
        mth::interpolatePoints(vPos[vIds[8]], triCenter[vIds[7]], dc, tri11);

        mth::interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[2]], dc, edge20);
        mth::interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[6]], dc, edge21);
        mth::interpolatePoints(
          vPos[vIds[9]], edgeCenters[vIds[10]], dc, edge22);
        mth::interpolatePoints(vPos[vIds[9]], triCenter[vIds[3]], dc, tri20);
        mth::interpolatePoints(vPos[vIds[9]], triCenter[vIds[7]], dc, tri21);

        // Write the triangle endpoints, cell indices, and labels
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

        m[0] = label[vIds[0]];
        m[1] = label[vIds[0]];
        m[2] = label[vIds[0]];
        m[3] = label[vIds[0]];
        m[4] = label[vIds[8]];
        m[5] = label[vIds[8]];
        m[6] = label[vIds[8]];
        m[7] = label[vIds[9]];
        m[8] = label[vIds[9]];
        m[9] = label[vIds[9]];
        m += 10;

      } else { // 4 labels
        std::array<float, 3> tetCenter{};
        mth::getCenter(vPos[0], vPos[1], vPos[2], vPos[3], tetCenter);

        // the 4 triangle centers
        std::array<std::array<float, 3>, 4> triCenter{};
        mth::getCenter(vPos[0], vPos[1], vPos[2], triCenter[0]);
        mth::getCenter(vPos[0], vPos[1], vPos[3], triCenter[1]);
        mth::getCenter(vPos[0], vPos[2], vPos[3], triCenter[2]);
        mth::getCenter(vPos[1], vPos[2], vPos[3], triCenter[3]);

        std::array<float, 3> vert00{}, vert01{}, vert02{}, vert0tet{},
          vert0t0{}, vert0t1{}, vert0t2{}, vert10{}, vert11{}, vert12{},
          vert1tet{}, vert1t0{}, vert1t1{}, vert1t2{}, vert20{}, vert21{},
          vert22{}, vert2tet{}, vert2t0{}, vert2t1{}, vert2t2{}, vert30{},
          vert31{}, vert32{}, vert3tet{}, vert3t0{}, vert3t1{}, vert3t2{};

        mth::interpolatePoints(vPos[0], vPos[1], d0, vert00);
        mth::interpolatePoints(vPos[0], vPos[2], d0, vert01);
        mth::interpolatePoints(vPos[0], vPos[3], d0, vert02);
        mth::interpolatePoints(vPos[0], tetCenter, dc, vert0tet);
        mth::interpolatePoints(vPos[0], triCenter[0], dc, vert0t0);
        mth::interpolatePoints(vPos[0], triCenter[1], dc, vert0t1);
        mth::interpolatePoints(vPos[0], triCenter[2], dc, vert0t2);

        mth::interpolatePoints(vPos[1], vPos[0], d0, vert10);
        mth::interpolatePoints(vPos[1], vPos[2], d0, vert11);
        mth::interpolatePoints(vPos[1], vPos[3], d0, vert12);
        mth::interpolatePoints(vPos[1], tetCenter, dc, vert1tet);
        mth::interpolatePoints(vPos[1], triCenter[0], dc, vert1t0);
        mth::interpolatePoints(vPos[1], triCenter[1], dc, vert1t1);
        mth::interpolatePoints(vPos[1], triCenter[3], dc, vert1t2);

        mth::interpolatePoints(vPos[2], vPos[0], d0, vert20);
        mth::interpolatePoints(vPos[2], vPos[1], d0, vert21);
        mth::interpolatePoints(vPos[2], vPos[3], d0, vert22);
        mth::interpolatePoints(vPos[2], tetCenter, dc, vert2tet);
        mth::interpolatePoints(vPos[2], triCenter[0], dc, vert2t0);
        mth::interpolatePoints(vPos[2], triCenter[2], dc, vert2t1);
        mth::interpolatePoints(vPos[2], triCenter[3], dc, vert2t2);

        mth::interpolatePoints(vPos[3], vPos[0], d0, vert30);
        mth::interpolatePoints(vPos[3], vPos[1], d0, vert31);
        mth::interpolatePoints(vPos[3], vPos[2], d0, vert32);
        mth::interpolatePoints(vPos[3], tetCenter, dc, vert3tet);
        mth::interpolatePoints(vPos[3], triCenter[1], dc, vert3t0);
        mth::interpolatePoints(vPos[3], triCenter[2], dc, vert3t1);
        mth::interpolatePoints(vPos[3], triCenter[3], dc, vert3t2);

        // Write the triangle endpoints, cell indices, and labels
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

        m[0] = label[0];
        m[1] = label[0];
        m[2] = label[0];
        m[3] = label[0];
        m[4] = label[0];
        m[5] = label[0];
        m[6] = label[1];
        m[7] = label[1];
        m[8] = label[1];
        m[9] = label[1];
        m[10] = label[1];
        m[11] = label[1];
        m[12] = label[2];
        m[13] = label[2];
        m[14] = label[2];
        m[15] = label[2];
        m[16] = label[2];
        m[17] = label[2];
        m[18] = label[3];
        m[19] = label[3];
        m[20] = label[3];
        m[21] = label[3];
        m[22] = label[3];
        m[23] = label[3];
        m += 24;
      }
    }
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Wrote detailed boundaries", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}