/// \ingroup base
/// \class ttk::MeshGraph
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.12.2018
///
/// \brief TTK %meshGraph processing package.
///
/// %MeshGraph is a TTK processing package that generates for each one
/// dimensional cell (edge) a two dimensional cell by mapping a size value to
/// the width of the input cell. The output is a set of either quadratic cells
/// or linear polygons.
///
/// This filter supports two modes:
///
///     1) Each cell (a,b) is mapped to two quadric quads:
///
///     a-------------------b
///
///     a1--m1a1--m1--b1m1--b1
///     |         |         |
///     a         c         b
///     |         |         |
///     a0--a0m0--m0--m0b0--b0
///
///     2) Each cell (a,b) is subdivided into a linear polygon
///
///     a----------------------b
///
///     a1--s0u--s1u-- ...  b1
///     |                   |
///     a0--s0d--s1d-- ...  b0
///

#pragma once

// base code includes
#include <Wrapper.h>

namespace ttk {

  class MeshGraph : virtual public Debug {

  public:
    MeshGraph() {
      this->setDebugMsgPrefix("MeshGraph");
    }
    ~MeshGraph() = default;

    inline size_t computeNumberOfOutputPoints(const size_t &nInputPoints,
                                              const size_t &nInputCells,
                                              const bool &useQuadraticCells,
                                              const size_t &nSubdivisions
                                              = 0) const {
      return useQuadraticCells
               ? nInputPoints * 3
                   + nInputCells * 7 // 3 per input point (a -> a,a0,a1) + 7 per
                                     // input cell (a0m0,m0,m0b0,b1m1,m1,m1a1,c)
               : nInputPoints * 2
                   + nInputCells
                       * (nSubdivisions
                          * 2); // 2 per input point (a -> a0,a1) + 2 per cell
                                // subdivision (sIu+sId)
    }

    inline size_t
      computeNumberOfOutputCells(const size_t &nInputCells,
                                 const bool &useQuadraticCells) const {
      return useQuadraticCells
               ? nInputCells
                   * 2 // each cell gets converted into two quadratic quads
               : nInputCells; // each cell gets converted into a one cell with
                              // multiple points
    }

    inline size_t computeOutputCellSize(const size_t &nSubdivisions) const {
      return 4 + nSubdivisions * 2; // 4 corners + 2 for each subdivision
    }

    inline size_t
      computeOutputConnectivityArraySize(const size_t &nInputCells,
                                         const bool &useQuadraticCells,
                                         const size_t &nSubdivisions
                                         = 0) const {
      return useQuadraticCells
               ? nInputCells
                   * 16 // 8 corners (a0,m0,m1,a1, m0,b0,b1,m1) + 8
                        // mid-edge nodes (a0m0,c,m1a1,a, m0b0,b,b1m1,c)
               : nInputCells * this->computeOutputCellSize(nSubdivisions);
    }

    // Mesh graph with quadratic quads
    template <typename IT, typename CT, typename DT>
    int execute(
      // Output
      CT *outputPoints,
      IT *outputConnectivityArray,
      IT *outputOffsetArray,

      // Input
      const CT *inputPoints,
      const IT *inputConnectivityArray,
      const size_t &nInputPoints,
      const size_t &nInputCells,

      const DT *inputPointSizes,
      const CT &sizeScale,
      const size_t &sizeAxis) const;

    // Mesh graph with linear polygon
    template <typename IT, typename CT, typename DT>
    int execute2(
      // Output
      CT *outputPoints,
      IT *outputConnectivityArray,
      IT *outputOffsetArray,

      // Input
      const CT *inputPoints,
      const IT *inputConnectivityArray,
      const size_t nInputPoints,
      const size_t nInputCells,
      const size_t nSubdivisions,

      const DT *inputPointSizes,
      const CT sizeScale,
      const size_t sizeAxis) const;

    // Map input point data to output point data
    template <typename DT, typename IT>
    int mapInputPointDataToOutputPointData(DT *outputPointData,
                                           const size_t &nInputPoints,
                                           const size_t &nInputCells,
                                           const IT *inputConnectivityArray,
                                           const DT *inputPointData,
                                           const bool &useQuadraticCells,
                                           const size_t &nSubdivisions
                                           = 0) const;

    // Map input point data to output point data
    template <typename DT>
    int mapInputCellDataToOutputCellData(DT *outputCellData,
                                         const size_t &nInputCells,
                                         const DT *inputCellData,
                                         const bool &useQuadraticCells) const;
  };
} // namespace ttk

// =============================================================================
// Version 1: Mesh Graph with Quadratic Quads
// =============================================================================
template <typename IT, typename CT, typename DT>
int ttk::MeshGraph::execute(
  // Output
  CT *outputPoints,
  IT *outputConnectivityArray,
  IT *outputOffsetArray,

  // Input
  const CT *inputPoints,
  const IT *inputConnectivityArray,
  const size_t &nInputPoints,
  const size_t &nInputCells,
  const DT *inputPointSizes,
  const CT &sizeScale,
  const size_t &sizeAxis) const {

  // Print Input
  this->printMsg(debug::Separator::L1);
  this->printMsg({{"Mode", "Quadratic Quads"},
                  {"#Nodes", std::to_string(nInputPoints)},
                  {"#Edges", std::to_string(nInputCells)}});
  this->printMsg(debug::Separator::L2);

  const size_t edgePointOffset = nInputPoints * 3;
  // ---------------------------------------------------------------------------
  // Compute Output Point Locations
  // ---------------------------------------------------------------------------
  // outputPoints: [
  //     3 per input point (a,a0,a1,b,b0,b1,...),
  //     7 per input cell (m0,m1,a0m0,m0b0,b1m1,m1a1,c...)
  // ]
  {
    Timer t;
    this->printMsg("Computing node locations", 0, debug::LineMode::REPLACE);

// -----------------------------------------------------------------------------
// Compute points that result from input points
// -----------------------------------------------------------------------------
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputPoints; i++) {
      const CT *coord = &inputPoints[i * 3];

      size_t q = i * 9; // i*3*3
      // a
      outputPoints[q] = coord[0];
      outputPoints[q + 1] = coord[1];
      outputPoints[q + 2] = coord[2];

      // a0
      outputPoints[q + 3] = coord[0];
      outputPoints[q + 4] = coord[1];
      outputPoints[q + 5] = coord[2];

      // a1
      outputPoints[q + 6] = coord[0];
      outputPoints[q + 7] = coord[1];
      outputPoints[q + 8] = coord[2];

      const CT &size = ((CT)inputPointSizes[i]) * sizeScale;

      // Move a0 and a1 along size axis
      outputPoints[q + 3 + sizeAxis] += size / 2;
      outputPoints[q + 6 + sizeAxis] -= size / 2;
    }

    // -------------------------------------------------------------------------
    // Compute points that result from edges
    // -------------------------------------------------------------------------

    // Lambda function that linearly interpolates two point locations
    auto getMidPoint = [&](const size_t &m, const size_t &i, const size_t &j) {
      size_t mp = m * 3;
      size_t ip = i * 3;
      size_t jp = j * 3;
      for(size_t k = 0; k < 3; k++)
        outputPoints[mp + k]
          = (outputPoints[ip + k] + outputPoints[jp + k]) / 2;
    };

    // Lambda function that computes the output location of a mid point on a
    // bezier curve
    auto getMidPoint2
      = [&](const size_t &m, const size_t &p0, const size_t &p1) {
          auto bezierPoint = [](CT &n, const CT &q0, const CT &q2) {
            n = 0.5 * (0.5 * q0 + 0.5 * n) + 0.5 * (0.5 * n + 0.5 * q2);
          };

          size_t mi = m * 3;
          size_t p0i = p0 * 3; // first point
          size_t p1i = p1 * 3; // second point

          for(size_t i = 0; i < 3; i++)
            outputPoints[mi + i]
              = (outputPoints[p0i + i] + outputPoints[p1i + i]) / 2;
          outputPoints[mi + sizeAxis] = outputPoints[p0i + sizeAxis];

          for(size_t i = 0; i < 3; i++)
            bezierPoint(outputPoints[mi + i], outputPoints[p0i + i],
                        outputPoints[p1i + i]);
        };

// Iterate over input cells and generate new points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputCells; i++) {
      size_t temp = i * 2;
      size_t aInputIndex = (size_t)inputConnectivityArray[temp++];
      size_t bInputIndex = (size_t)inputConnectivityArray[temp];

      // already computed points
      size_t a = aInputIndex * 3;
      size_t a0 = a + 1;
      size_t a1 = a + 2;
      size_t b = bInputIndex * 3;
      size_t b0 = b + 1;
      size_t b1 = b + 2;

      // points to compute
      size_t offset = edgePointOffset + i * 7;
      size_t m0 = offset;
      size_t m1 = offset + 1;
      size_t a0m0 = offset + 2;
      size_t m0b0 = offset + 3;
      size_t b1m1 = offset + 4;
      size_t m1a1 = offset + 5;
      size_t c = offset + 6;

      getMidPoint(m0, a0, b0);
      getMidPoint(m1, a1, b1);

      getMidPoint(c, m0, m1);

      getMidPoint2(a0m0, a0, m0);
      getMidPoint2(m0b0, b0, m0);

      getMidPoint2(b1m1, b1, m1);
      getMidPoint2(m1a1, a1, m1);
    }

    // Print Status
    this->printMsg(
      "Computing mesh vertices", 1, t.getElapsedTime(), this->threadNumber_);
  }

  // ---------------------------------------------------------------------------
  // Compute Output Cells
  // ---------------------------------------------------------------------------
  {
    Timer t;
    this->printMsg("Computing mesh cells", 0, debug::LineMode::REPLACE);

    IT edgePointOffset_ = (IT)edgePointOffset;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputCells; i++) {
      size_t temp = i * 2;
      IT aInputIndex = inputConnectivityArray[temp++];
      IT bInputIndex = inputConnectivityArray[temp];

      // get point indicies
      IT a = aInputIndex * 3;
      IT a0 = a + 1;
      IT a1 = a + 2;
      IT b = bInputIndex * 3;
      IT b0 = b + 1;
      IT b1 = b + 2;

      IT i_ = (IT)i;
      IT offset = edgePointOffset_ + i_ * 7;
      IT m0 = offset;
      IT m1 = offset + 1;
      IT a0m0 = offset + 2;
      IT m0b0 = offset + 3;
      IT b1m1 = offset + 4;
      IT m1a1 = offset + 5;
      IT c = offset + 6;

      // output cell offset
      size_t q = i * 16;

      // first quadratic quad
      outputConnectivityArray[q++] = a0;
      outputConnectivityArray[q++] = m0;
      outputConnectivityArray[q++] = m1;
      outputConnectivityArray[q++] = a1;

      outputConnectivityArray[q++] = a0m0;
      outputConnectivityArray[q++] = c;
      outputConnectivityArray[q++] = m1a1;
      outputConnectivityArray[q++] = a;

      // second quadratic quad
      outputConnectivityArray[q++] = m0;
      outputConnectivityArray[q++] = b0;
      outputConnectivityArray[q++] = b1;
      outputConnectivityArray[q++] = m1;

      outputConnectivityArray[q++] = m0b0;
      outputConnectivityArray[q++] = b;
      outputConnectivityArray[q++] = b1m1;
      outputConnectivityArray[q++] = c;
    }

    for(size_t i = 0; i <= 2 * nInputCells; i++) {
      outputOffsetArray[i] = i * 8;
    }

    this->printMsg(
      "Computing mesh cells", 1, t.getElapsedTime(), this->threadNumber_);
  }

  return 1;
}

// =============================================================================
// Version 2: Mesh Graph with Linear Polygon
// =============================================================================
template <typename IT, typename CT, typename DT>
int ttk::MeshGraph::execute2(
  // Output
  CT *outputPoints,
  IT *outputConnectivityArray,
  IT *outputOffsetArray,

  // Input
  const CT *inputPoints,
  const IT *inputConnectivityArray,
  const size_t nInputPoints,
  const size_t nInputCells,
  const size_t nSubdivisions,
  const DT *inputPointSizes,
  const CT sizeScale,
  const size_t sizeAxis) const {

  this->printMsg(debug::Separator::L1);
  this->printMsg({{"Mode", "Linear Polygon"},
                  {"#Nodes", std::to_string(nInputPoints)},
                  {"#Edges", std::to_string(nInputCells)},
                  {"#Subdivisions", std::to_string(nSubdivisions)}});
  this->printMsg(debug::Separator::L2);

  size_t subdivisionOffset = nInputPoints * 2;
  size_t nSubdivisionPoints = nSubdivisions * 2;
  size_t outputPointsSubdivisonOffset = nSubdivisionPoints * 3;

  // ---------------------------------------------------------------------------
  // Compute Output Point Locations
  // ---------------------------------------------------------------------------
  // outputPoints: [
  //     corners: 2*inputPoints in inputPoints order;
  //     SubPoints: 2 per subdivison in cell order]
  // ]
  {
    Timer t;
    this->printMsg("Computing mesh vertices", 0, debug::LineMode::REPLACE);

// -----------------------------------------------------------------------------
// Compute Corners
// -----------------------------------------------------------------------------
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputPoints; i++) {
      const CT *coords = &inputPoints[i * 3];

      size_t q = i * 6;
      outputPoints[q] = coords[0];
      outputPoints[q + 1] = coords[1];
      outputPoints[q + 2] = coords[2];

      outputPoints[q + 3] = coords[0];
      outputPoints[q + 4] = coords[1];
      outputPoints[q + 5] = coords[2];

      const CT &size = ((CT)inputPointSizes[i]) * sizeScale;

      outputPoints[q + sizeAxis] += size / 2;
      outputPoints[q + 3 + sizeAxis] -= size / 2;
    }

    // -------------------------------------------------------------------------
    // Compute SubPoints
    // -------------------------------------------------------------------------
    size_t q = subdivisionOffset * 3;
    float nSubdivisionsP1 = nSubdivisions + 1;

    auto computeBezierPoint = [&](const size_t &no0, const size_t &no1,
                                  const size_t &subdOffset,
                                  const float lambda) {
      float lambdaI = 1 - lambda;

      float lambda_2 = lambda * lambda;
      float lambda_3 = lambda * lambda_2;

      float lambdaI_2 = lambdaI * lambdaI;
      float lambdaI_3 = lambdaI * lambdaI_2;

      float m0[3];
      float m1[3];
      for(size_t i = 0; i < 3; i++) {
        m0[i] = (outputPoints[no0 + i] + outputPoints[no1 + i]) / 2;
        m1[i] = m0[i];
      }
      m0[sizeAxis] = outputPoints[no0 + sizeAxis];
      m1[sizeAxis] = outputPoints[no1 + sizeAxis];

      for(size_t i = 0; i < 3; i++)
        outputPoints[subdOffset + i]
          = lambdaI_3 * outputPoints[no0 + i] + 3 * lambdaI_2 * lambda * m0[i]
            + 3 * lambdaI * lambda_2 * m1[i] + lambda_3 * outputPoints[no1 + i];
    };

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputCells; i++) {
      size_t temp = i * 2;
      IT n0 = inputConnectivityArray[temp++];
      IT n1 = inputConnectivityArray[temp];

      IT no0 = n0 * 6;
      IT no1 = n1 * 6;

      size_t q2 = q + i * outputPointsSubdivisonOffset;
      for(size_t j = 1; j <= nSubdivisions; j++) {
        computeBezierPoint(no0, no1, q2, j / nSubdivisionsP1);
        computeBezierPoint(no0 + 3, no1 + 3, q2 + 3, j / nSubdivisionsP1);

        q2 += 6;
      }
    }

    this->printMsg(
      "Computing mesh vertices", 1, t.getElapsedTime(), this->threadNumber_);
  }

  // ---------------------------------------------------------------------------
  // Compute Output Cells
  // ---------------------------------------------------------------------------
  {
    Timer t;
    this->printMsg("Computing mesh cells", 0, debug::LineMode::REPLACE);

    size_t cellSize = this->computeOutputCellSize(nSubdivisions);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputCells; i++) {
      size_t q = i * 2;
      IT in0 = inputConnectivityArray[q++] * 2;
      IT in1 = inputConnectivityArray[q] * 2;

      IT c0 = in0;
      IT c1 = in0 + 1;
      IT c2 = in1 + 1;
      IT c3 = in1;

      size_t q2 = cellSize * i;

      outputConnectivityArray[q2++] = c0;
      outputConnectivityArray[q2++] = c1;

      IT temp = subdivisionOffset + i * nSubdivisionPoints;

      for(size_t j = 0; j < nSubdivisions; j++)
        outputConnectivityArray[q2++] = temp + j * 2 + 1;

      outputConnectivityArray[q2++] = c2;
      outputConnectivityArray[q2++] = c3;

      for(int j = nSubdivisions - 1; j >= 0; j--)
        outputConnectivityArray[q2++] = temp + j * 2;
    }

    for(size_t i = 0; i <= nInputCells; i++)
      outputOffsetArray[i] = i * cellSize;

    this->printMsg(
      "Computing mesh cells", 1, t.getElapsedTime(), this->threadNumber_);
  }

  return 1;
}

// =============================================================================
// Map input point data to output point data
// =============================================================================
template <typename DT, typename IT>
int ttk::MeshGraph::mapInputPointDataToOutputPointData(
  DT *outputPointData,

  const size_t &nInputPoints,
  const size_t &nInputCells,
  const IT *inputConnectivityArray,
  const DT *inputPointData,
  const bool &useQuadraticCells,
  const size_t &nSubdivisions) const {

  if(useQuadraticCells) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputPoints; i++) {
      size_t q = i * 3;
      // a, a0, a1
      outputPointData[q] = inputPointData[i];
      outputPointData[q + 1] = inputPointData[i];
      outputPointData[q + 2] = inputPointData[i];
    }

    size_t edgePointOffset = nInputPoints * 3;

// Iterate over input cells and assign point data of intermediate points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputCells; i++) {
      size_t temp = i * 2;
      size_t aInputIndex = (size_t)inputConnectivityArray[temp++];
      size_t bInputIndex = (size_t)inputConnectivityArray[temp];

      size_t a = aInputIndex * 3;
      size_t b = bInputIndex * 3;

      size_t offset = edgePointOffset + i * 7;
      size_t m0 = offset;
      size_t m1 = offset + 1;
      size_t a0m0 = offset + 2;
      size_t m0b0 = offset + 3;
      size_t b1m1 = offset + 4;
      size_t m1a1 = offset + 5;
      size_t c = offset + 6;

      outputPointData[c] = (DT)((outputPointData[a] + outputPointData[b]) / 2);
      outputPointData[m0] = outputPointData[c];
      outputPointData[m1] = outputPointData[c];

      outputPointData[a0m0]
        = (DT)((outputPointData[a] + outputPointData[c]) / 2);
      outputPointData[m1a1] = outputPointData[a0m0];

      outputPointData[m0b0]
        = (DT)((outputPointData[c] + outputPointData[b]) / 2);
      outputPointData[b1m1] = outputPointData[m0b0];
    }
  } else {

// Corners
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputPoints; i++) {
      size_t offset = i * 2;
      auto &v = inputPointData[i];
      outputPointData[offset] = v;
      outputPointData[offset + 1] = v;
    }

    // Intermediate Points
    size_t subdivisionOffset = nInputPoints * 2;
    size_t nSubdivisionPoints = nSubdivisions * 2;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputCells; i++) {
      size_t q = i * 2;
      IT c0 = inputConnectivityArray[q];
      DT c0V = inputPointData[c0];

      size_t temp = subdivisionOffset + i * nSubdivisionPoints;
      for(size_t j = 0; j < nSubdivisions; j++) {
        size_t q2 = temp + j * 2;
        outputPointData[q2] = c0V;
        outputPointData[q2 + 1] = c0V;
      }
    }
  }

  return 1;
}

// =============================================================================
// Map input cell data to output cell data
// =============================================================================
template <typename DT>
int ttk::MeshGraph::mapInputCellDataToOutputCellData(
  DT *outputCellData,
  const size_t &nInputCells,
  const DT *inputCellData,
  const bool &useQuadraticCells) const {

  if(useQuadraticCells) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputCells; i++) {
      size_t offset = i * 2;
      outputCellData[offset] = inputCellData[i];
      outputCellData[offset + 1] = inputCellData[i];
    }
  } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputCells; i++) {
      outputCellData[i] = inputCellData[i];
    }
  }

  return 1;
}
