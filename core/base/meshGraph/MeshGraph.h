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

using namespace std;

namespace ttk {

  class MeshGraph : public Debug {

  public:
    MeshGraph(){};
    ~MeshGraph(){};

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
    };

    inline size_t
      computeNumberOfOutputCells(const size_t &nInputCells,
                                 const bool &useQuadraticCells) const {
      return useQuadraticCells
               ? nInputCells
                   * 2 // each cell gets converted into two quadratic quads
               : nInputCells; // each cell gets converted into a one cell with
                              // multiple points
    };

    inline size_t computeOutputCellSize(const size_t &nSubdivisions) const {
      return 5
             + nSubdivisions
                 * 2; // cellDim + 4 corners + 2 for each subdivision
    };

    inline size_t computeOutputTopologySize(const size_t &nInputCells,
                                            const bool &useQuadraticCells,
                                            const size_t &nSubdivisions
                                            = 0) const {
      return useQuadraticCells
               ? nInputCells
                   * 18 // 2*cellDim + 8 corners (a0,m0,m1,a1, m0,b0,b1,m1) + 8
                        // mid-edge nodes (a0m0,c,m1a1,a, m0b0,b,b1m1,c)
               : nInputCells * this->computeOutputCellSize(nSubdivisions);
    };

    // Mesh graph with quadratic quads
    template <typename topoType, typename sizeType>
    int execute(
      // Input
      const float *inputPoints,
      const topoType *inputTopology,
      const size_t &nInputPoints,
      const size_t &nInputCells,

      const sizeType *inputPointSizes,
      const float &sizeScale,
      const size_t &sizeAxis,

      // Output
      float *outputPoints,
      topoType *outputTopology) const;

    // Mesh graph with linear polygon
    template <typename topoType, typename sizeType>
    int execute2(
      // Input
      const float *inputPoints,
      const topoType *inputTopology,
      const size_t nInputPoints,
      const size_t nInputCells,
      const size_t nSubdivisions,

      const sizeType *inputPointSizes,
      const float sizeScale,
      const size_t sizeAxis,

      // Output
      float *outputPoints,
      topoType *outputTopology) const;

    // Map input point data to output point data
    template <typename topoType, typename dataType>
    int mapInputPointDataToOutputPointData(const topoType *inputTopology,
                                           const size_t &nInputPoints,
                                           const size_t &nInputCells,

                                           const dataType *inputPointData,
                                           dataType *outputPointData,

                                           const bool &useQuadraticCells,
                                           const size_t &nSubdivisions
                                           = 0) const;

    // Map input point data to output point data
    template <typename topoType, typename dataType>
    int mapInputCellDataToOutputCellData(const size_t &nInputCells,

                                         const dataType *inputCellData,
                                         dataType *outputCellData,

                                         const bool &useQuadraticCells,
                                         const size_t &nSubdivisions = 0) const;
  };
} // namespace ttk

// =============================================================================
// Version 1: Mesh Graph with Quadratic Quads
// =============================================================================
template <typename topoType, typename sizeType>
int ttk::MeshGraph::execute(
  // Input
  const float *inputPoints,
  const topoType *inputTopology,
  const size_t &nInputPoints,
  const size_t &nInputCells,

  const sizeType *inputPointSizes,
  const float &sizeScale,
  const size_t &sizeAxis,

  // Output
  float *outputPoints,
  topoType *outputTopology) const {

  Timer t;
  double t0 = 0;

  // Print Input
  {
    stringstream msg;
    msg << "[ttkMeshGraph] Computing quadratic quads for graph with" << endl
        << "[ttkMeshGraph]  - " << nInputPoints << " points" << endl
        << "[ttkMeshGraph]  - " << nInputCells << " edges" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  auto getInputPointData = [&](const size_t &pointIndex, float data[4]) {
    size_t i = pointIndex * 3;
    data[0] = inputPoints[i++];
    data[1] = inputPoints[i++];
    data[2] = inputPoints[i];

    data[3] = (inputPointSizes != nullptr)
                ? ((float)inputPointSizes[pointIndex]) * sizeScale
                : sizeScale;
  };

  // -------------------------------------------------------------------------
  // Compute Output Point Locations
  // -------------------------------------------------------------------------
  // outputPoints: [
  //     3 per input point (a,a0,a1,b,b0,b1,...),
  //     7 per input cell (m0,m1,a0m0,m0b0,b1m1,m1a1,c...)
  // ]
  size_t edgePointOffset = 3 * nInputPoints;
  {
    dMsg(cout, "[ttkMeshGraph] Computing output points ... ", timeMsg);
    t0 = t.getElapsedTime();

// ---------------------------------------------------------------------
// Compute points that result from input points
// ---------------------------------------------------------------------
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputPoints; i++) {
      float data[4];
      getInputPointData(i, data);

      size_t q = i * 9; // i*3*3
      // a
      outputPoints[q] = data[0];
      outputPoints[q + 1] = data[1];
      outputPoints[q + 2] = data[2];

      // a0
      outputPoints[q + 3] = data[0];
      outputPoints[q + 4] = data[1];
      outputPoints[q + 5] = data[2];

      // a1
      outputPoints[q + 6] = data[0];
      outputPoints[q + 7] = data[1];
      outputPoints[q + 8] = data[2];

      // Move a0 and a1 along size axis
      outputPoints[q + 3 + sizeAxis] += data[3] / 2;
      outputPoints[q + 6 + sizeAxis] -= data[3] / 2;
    }

    // ---------------------------------------------------------------------
    // Compute points that result from edges
    // ---------------------------------------------------------------------

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
          auto bezierPoint = [](float &n, const float &q0, const float &q2) {
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
      size_t temp = i * 3 + 1;
      size_t aInputIndex = (size_t)inputTopology[temp++];
      size_t bInputIndex = (size_t)inputTopology[temp];

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
    {
      stringstream msg;
      msg << "done (" << (t.getElapsedTime() - t0) << " s)." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  // -------------------------------------------------------------------------
  // Compute Output Cells
  // -------------------------------------------------------------------------
  {
    dMsg(cout, "[ttkMeshGraph] Computing output cells  ... ", timeMsg);
    t0 = t.getElapsedTime();

    topoType edgePointOffset_ = (topoType)edgePointOffset;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputCells; i++) {
      size_t temp = i * 3 + 1;
      topoType aInputIndex = inputTopology[temp++];
      topoType bInputIndex = inputTopology[temp];

      // get point indicies
      topoType a = aInputIndex * 3;
      topoType a0 = a + 1;
      topoType a1 = a + 2;
      topoType b = bInputIndex * 3;
      topoType b0 = b + 1;
      topoType b1 = b + 2;

      topoType i_ = (topoType)i;
      topoType offset = edgePointOffset_ + i_ * 7;
      topoType m0 = offset;
      topoType m1 = offset + 1;
      topoType a0m0 = offset + 2;
      topoType m0b0 = offset + 3;
      topoType b1m1 = offset + 4;
      topoType m1a1 = offset + 5;
      topoType c = offset + 6;

      // output cell offset
      size_t q = i * 18;

      // first quadratic quad
      outputTopology[q++] = 8;
      outputTopology[q++] = a0;
      outputTopology[q++] = m0;
      outputTopology[q++] = m1;
      outputTopology[q++] = a1;

      outputTopology[q++] = a0m0;
      outputTopology[q++] = c;
      outputTopology[q++] = m1a1;
      outputTopology[q++] = a;

      // second quadratic quad
      outputTopology[q++] = 8;
      outputTopology[q++] = m0;
      outputTopology[q++] = b0;
      outputTopology[q++] = b1;
      outputTopology[q++] = m1;

      outputTopology[q++] = m0b0;
      outputTopology[q++] = b;
      outputTopology[q++] = b1m1;
      outputTopology[q++] = c;
    }

    // Print Status
    {
      stringstream msg;
      msg << "done (" << (t.getElapsedTime() - t0) << " s)." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return 1;
}

// =============================================================================
// Version 2: Mesh Graph with Linear Polygon
// =============================================================================
template <typename topoType, typename sizeType>
int ttk::MeshGraph::execute2(
  // Input
  const float *inputPoints,
  const topoType *inputTopology,
  const size_t nInputPoints,
  const size_t nInputCells,
  const size_t nSubdivisions,

  const sizeType *inputPointSizes,
  const float sizeScale,
  const size_t sizeAxis,

  // Output
  float *outputPoints,
  topoType *outputTopology) const {

  Timer t;
  double t0 = 0;

  // Print Input
  {
    stringstream msg;
    msg << "[ttkMeshGraph] Computing linear polygon for graph with" << endl
        << "[ttkMeshGraph]  - " << nInputPoints << " points" << endl
        << "[ttkMeshGraph]  - " << nInputCells << " edges" << endl
        << "[ttkMeshGraph]  - " << nSubdivisions << " subdivisions" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  auto getInputPointData = [&](const size_t &pointIndex, float data[4]) {
    size_t i = pointIndex * 3;
    data[0] = inputPoints[i++];
    data[1] = inputPoints[i++];
    data[2] = inputPoints[i];

    data[3] = (inputPointSizes != nullptr)
                ? ((float)inputPointSizes[pointIndex]) * sizeScale
                : sizeScale;
  };

  size_t subdivisionOffset = nInputPoints * 2;
  size_t nSubdivisionPoints = nSubdivisions * 2;
  size_t outputPointsSubdivisonOffset = nSubdivisionPoints * 3;

  // -------------------------------------------------------------------------
  // Compute Output Point Locations
  // -------------------------------------------------------------------------
  // outputPoints: [
  //     corners: 2*inputPoints in inputPoints order;
  //     SubPoints: 2 per subdivison in cell order]
  // ]
  {
    dMsg(cout, "[ttkMeshGraph] Computing output points ... ", timeMsg);
    t0 = t.getElapsedTime();

// ---------------------------------------------------------------------
// Compute Corners
// ---------------------------------------------------------------------
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputPoints; i++) {
      float data[4];
      getInputPointData(i, data);

      size_t q = i * 6;
      outputPoints[q] = data[0];
      outputPoints[q + 1] = data[1];
      outputPoints[q + 2] = data[2];

      outputPoints[q + 3] = data[0];
      outputPoints[q + 4] = data[1];
      outputPoints[q + 5] = data[2];

      outputPoints[q + sizeAxis] += data[3] / 2;
      outputPoints[q + 3 + sizeAxis] -= data[3] / 2;
    }

    // ---------------------------------------------------------------------
    // Compute SubPoints
    // ---------------------------------------------------------------------
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
      size_t temp = i * 3 + 1;
      size_t n0 = (size_t)inputTopology[temp++];
      size_t n1 = (size_t)inputTopology[temp];

      size_t no0 = n0 * 6;
      size_t no1 = n1 * 6;

      size_t q2 = q + i * outputPointsSubdivisonOffset;
      for(float j = 1; j <= nSubdivisions; j++) {
        computeBezierPoint(no0, no1, q2, j / nSubdivisionsP1);
        computeBezierPoint(no0 + 3, no1 + 3, q2 + 3, j / nSubdivisionsP1);

        q2 += 6;
      }
    }

    // Print Status
    {
      stringstream msg;
      msg << "done (" << (t.getElapsedTime() - t0) << " s)." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  // -------------------------------------------------------------------------
  // Compute Output Cells
  // -------------------------------------------------------------------------
  {
    dMsg(cout, "[ttkMeshGraph] Computing output cells  ... ", timeMsg);
    t0 = t.getElapsedTime();

    size_t cellSize = this->computeOutputCellSize(nSubdivisions);
    topoType cellDim = ((topoType)cellSize) - 1;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < nInputCells; i++) {
      size_t q = i * 3 + 1;
      topoType in0 = inputTopology[q++] * 2;
      topoType in1 = inputTopology[q] * 2;

      topoType c0 = in0;
      topoType c1 = in0 + 1;
      topoType c2 = in1 + 1;
      topoType c3 = in1;

      size_t q2 = cellSize * i;
      outputTopology[q2++] = cellDim;

      outputTopology[q2++] = c0;
      outputTopology[q2++] = c1;

      size_t temp = subdivisionOffset + i * nSubdivisionPoints;

      for(size_t j = 0; j < nSubdivisions; j++)
        outputTopology[q2++] = (topoType)(temp + j * 2 + 1);

      outputTopology[q2++] = c2;
      outputTopology[q2++] = c3;

      for(int j = nSubdivisions - 1; j >= 0; j--)
        outputTopology[q2++] = (topoType)(temp + j * 2);
    }

    // Print Status
    {
      stringstream msg;
      msg << "done (" << (t.getElapsedTime() - t0) << " s)." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return 1;
}

// =============================================================================
// Map input point data to output point data
// =============================================================================
template <typename topoType, typename dataType>
int ttk::MeshGraph::mapInputPointDataToOutputPointData(
  const topoType *inputTopology,
  const size_t &nInputPoints,
  const size_t &nInputCells,

  const dataType *inputPointData,
  dataType *outputPointData,

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
      size_t temp = i * 3 + 1;
      size_t aInputIndex = (size_t)inputTopology[temp++];
      size_t bInputIndex = (size_t)inputTopology[temp];

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

      outputPointData[c]
        = (dataType)((outputPointData[a] + outputPointData[b]) / 2);
      outputPointData[m0] = outputPointData[c];
      outputPointData[m1] = outputPointData[c];

      outputPointData[a0m0]
        = (dataType)((outputPointData[a] + outputPointData[c]) / 2);
      outputPointData[m1a1] = outputPointData[a0m0];

      outputPointData[m0b0]
        = (dataType)((outputPointData[c] + outputPointData[b]) / 2);
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
      size_t q = i * 3 + 1;
      topoType c0 = inputTopology[q];
      dataType c0V = inputPointData[c0];

      size_t temp = subdivisionOffset + i * nSubdivisionPoints;
      for(size_t j = 0; j < nSubdivisions; j++) {
        size_t q2 = temp + j * 2;
        outputPointData[q2] = c0V;
        outputPointData[q2 + 1] = c0V;
      }
    }
  }

  return 1;
};

// =============================================================================
// Map input cell data to output cell data
// =============================================================================
template <typename topoType, typename dataType>
int ttk::MeshGraph::mapInputCellDataToOutputCellData(
  const size_t &nInputCells,

  const dataType *inputCellData,
  dataType *outputCellData,

  const bool &useQuadraticCells,
  const size_t &nSubdivisions) const {

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
};
