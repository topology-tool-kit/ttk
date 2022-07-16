/// \ingroup base
/// \class ttk::BarycentricSubdivision
/// \author Pierre Guillou (pierre.guillou@lip6.fr)
/// \date July 2019
///
/// \brief Subdivise a triangulation according to triangle barycenter
///
/// %BarycentricSubdivision generates a new, finer triangulation from
/// an input triangulation. Every triangle is divided in six new
/// triangles using the 3 edges middles and the triangle barycenter.
///
/// Scalar data on vertices (point data) with continuous values
/// (float/double) can be interpolated on the new
/// triangulation. Scalar data on input triangles can be replicated on
/// the new triangles.
///
/// \sa ttk::Triangulation
/// \sa ttkBarycentricSubdivision.cpp %for a usage example.

#pragma once

#include <array>
#include <type_traits>

// base code includes
#include <Geometry.h>
#include <Triangulation.h>

namespace ttk {

  class BarycentricSubdivision : virtual public Debug {

  public:
    BarycentricSubdivision() {
      this->setDebugMsgPrefix("BarycentricSubdivision");
    }

    inline void
      preconditionTriangulation(AbstractTriangulation *const triangulation) {
      if(triangulation == nullptr) {
        return;
      }
      triangulation->preconditionVertexNeighbors();
      triangulation->preconditionEdges();
      triangulation->preconditionTriangles();
      triangulation->preconditionTriangleEdges();
      nVertices_ = triangulation->getNumberOfVertices();
      nEdges_ = triangulation->getNumberOfEdges();
      nTriangles_ = triangulation->getNumberOfTriangles();
    }

    /** @brief Return the number of vertices in the output triangulation
     */
    inline SimplexId getNumberOfVertices() const {
      return nVertices_ + nEdges_ + nTriangles_;
    }

    /** @brief Return the number of triangles in the output triangulation
     */
    inline SimplexId getNumberOfTriangles() const {
      return nTriangles_ * 6;
    }

    template <typename triangulationType>
    int execute(const triangulationType &inputTriangl,
                ExplicitTriangulation &outputTriangl);

    /**
     * @brief Interpolate floating-point point data on subdivised triangulation
     *
     * Copy values on parent vertices, interpolate on edges and barycenters
     *
     * @param[in] data Pointer to input data on parent triangulation
     * @param[out] output Allocated buffer to be filled
     * @param[in] inputTriangl Input triangulation object
     *
     * @return 0 in case of success
     */
    template <typename T, typename triangulationType>
    int interpolateContinuousScalarField(
      const T *data, T *output, const triangulationType &inputTriangl) const {
      static_assert(
        std::is_floating_point<T>::value, "Floating point type required.");
      const auto nOutVerts = this->getNumberOfVertices();
      if(nOutVerts < 0 || nOutVerts != nVertices_ + nEdges_ + nTriangles_) {
        return 1;
      }

      // copy data on parent vertices
      std::copy(data, data + nVertices_, output);

      // interpolate on edges
      for(SimplexId i = 0; i < nEdges_; ++i) {
        SimplexId a{}, b{};
        inputTriangl.getEdgeVertex(i, 0, a);
        inputTriangl.getEdgeVertex(i, 0, b);
        output[nVertices_ + i] = (data[a] + data[b]) / T{2.0};
      }

      // interpolate on triangle barycenters
      for(SimplexId i = 0; i < nTriangles_; ++i) {
        SimplexId a{}, b{}, c{};
        inputTriangl.getTriangleVertex(i, 0, a);
        inputTriangl.getTriangleVertex(i, 1, b);
        inputTriangl.getTriangleVertex(i, 2, c);
        output[nVertices_ + nEdges_ + i]
          = (data[a] + data[b] + data[c]) / T{3.0};
      }
      return 0;
    }

    /**
     * @brief Interpolate integer point data on subdivised triangulation
     *
     * Copy values on parent vertices, put 0 elsewhere.
     *
     * @param[in] data Pointer to input data on parent triangulation
     * @param[out] output Allocated buffer to be filled
     *
     * @return 0 in case of success
     */
    template <typename T>
    int interpolateDiscreteScalarField(const T *data, T *output) const {

      static_assert(std::is_integral<T>::value, "Integral type required.");
      const auto nOutVerts = this->getNumberOfVertices();
      if(nOutVerts < 0 || nOutVerts < nVertices_) {
        return 1;
      }
      std::fill(output, output + nOutVerts, T{0});
      std::copy(data, data + nVertices_, output);
      return 0;
    }

    /**
     * @brief Interpolate cell data on subdivised triangulation
     *
     * Copy parent triangle value on the six children triangles.
     *
     * @param[in] data Pointer to input data on parent triangulation
     * @param[out] output Allocated buffer to be filled
     *
     * @return 0 in case of success
     */
    template <typename T>
    int interpolateCellDataField(const T *data, T *output) const {

      const size_t newTrianglesPerParent{6};
      const size_t nOutTriangles = this->getNumberOfTriangles();
      if(nOutTriangles != newTrianglesPerParent * nTriangles_) {
        return 1;
      }
      for(SimplexId i = 0; i < nTriangles_; ++i) {
        output[i * newTrianglesPerParent + 0] = data[i];
        output[i * newTrianglesPerParent + 1] = data[i];
        output[i * newTrianglesPerParent + 2] = data[i];
        output[i * newTrianglesPerParent + 3] = data[i];
        output[i * newTrianglesPerParent + 4] = data[i];
        output[i * newTrianglesPerParent + 5] = data[i];
      }
      return 0;
    }

  private:
    template <typename triangulationType>
    int subdiviseTriangulation(const triangulationType &inputTriangl);

    int buildOutputTriangulation(ExplicitTriangulation &outputTriangl);

    // input triangulation properties
    SimplexId nVertices_{};
    SimplexId nEdges_{};
    SimplexId nTriangles_{};

  protected:
    // output 3D coordinates of generated points: old points first, then edge
    // middles, then triangle barycenters
    std::vector<float> points_;
    // output triangles
    std::vector<LongSimplexId> cells_connectivity_;
    std::vector<LongSimplexId> cells_offsets_;
    // generated point cell id
    std::vector<SimplexId> pointId_;
    // generated points dimension: 0 vertex of parent triangulation, 1 edge
    // middle, 2 triangle barycenter
    std::vector<SimplexId> pointDim_;
  };
} // namespace ttk

template <typename triangulationType>
int ttk::BarycentricSubdivision::subdiviseTriangulation(
  const triangulationType &inputTriangl) {

  const SimplexId newPoints{nVertices_ + nEdges_ + nTriangles_};
  const size_t dataPerPoint{3};
  points_.clear();
  points_.resize(newPoints * dataPerPoint);

  const size_t dataPerCell{3};
  const size_t newTrianglesPerParent{6};
  cells_connectivity_.clear();
  cells_offsets_.clear();
  cells_connectivity_.resize(nTriangles_ * newTrianglesPerParent * dataPerCell);
  cells_offsets_.resize(nTriangles_ * newTrianglesPerParent + 1);

  pointId_.clear();
  pointId_.resize(newPoints);

  pointDim_.clear();
  pointDim_.resize(newPoints);

  // set input point coordinates and ids
  for(SimplexId i = 0; i < nVertices_; ++i) {
    inputTriangl.getVertexPoint(
      i, points_[3 * i + 0], points_[3 * i + 1], points_[3 * i + 2]);
    pointId_[i] = i;
  }

  // reserve memory for new points
  points_.reserve(newPoints * dataPerPoint);

  // generate edge middles
  for(SimplexId i = 0; i < nEdges_; ++i) {
    // edge vertices
    SimplexId a{}, b{};
    inputTriangl.getEdgeVertex(i, 0, a);
    inputTriangl.getEdgeVertex(i, 1, b);

    std::array<float, 3> pa{}, pb{}, mid{};
    inputTriangl.getVertexPoint(a, pa[0], pa[1], pa[2]);
    inputTriangl.getVertexPoint(b, pb[0], pb[1], pb[2]);

    mid[0] = (pa[0] + pb[0]) / 2.0F;
    mid[1] = (pa[1] + pb[1]) / 2.0F;
    mid[2] = (pa[2] + pb[2]) / 2.0F;

    const size_t offset = dataPerPoint * (nVertices_ + i);
    points_[offset + 0] = mid[0];
    points_[offset + 1] = mid[1];
    points_[offset + 2] = mid[2];
    pointId_[nVertices_ + i] = i;
    pointDim_[nVertices_ + i] = 1;
  }

  // generate triangle barycenters
  for(SimplexId i = 0; i < nTriangles_; ++i) {
    // triangle vertices
    SimplexId a{}, b{}, c{};
    inputTriangl.getTriangleVertex(i, 0, a);
    inputTriangl.getTriangleVertex(i, 1, b);
    inputTriangl.getTriangleVertex(i, 2, c);

    std::array<float, 3> pa{}, pb{}, pc{}, bary{};
    inputTriangl.getVertexPoint(a, pa[0], pa[1], pa[2]);
    inputTriangl.getVertexPoint(b, pb[0], pb[1], pb[2]);
    inputTriangl.getVertexPoint(c, pc[0], pc[1], pc[2]);

    bary[0] = (pa[0] + pb[0] + pc[0]) / 3.0F;
    bary[1] = (pa[1] + pb[1] + pc[1]) / 3.0F;
    bary[2] = (pa[2] + pb[2] + pc[2]) / 3.0F;

    const size_t offset = dataPerPoint * (nVertices_ + nEdges_ + i);
    points_[offset + 0] = bary[0];
    points_[offset + 1] = bary[1];
    points_[offset + 2] = bary[2];
    pointId_[nVertices_ + nEdges_ + i] = i;
    pointDim_[nVertices_ + nEdges_ + i] = 2;
  }

  LongSimplexId off_id = 0, co_id = 0;
  // subdivise every triangle
  for(SimplexId i = 0; i < nTriangles_; ++i) {
    // id of triangle barycenter
    SimplexId bary = nVertices_ + nEdges_ + i;

    for(SimplexId j = 0; j < inputTriangl.getTriangleEdgeNumber(i); ++j) {
      // edge id
      SimplexId e{};
      inputTriangl.getTriangleEdge(i, j, e);

      // id of middle of edge e
      SimplexId em = nVertices_ + e;

      // edge vertices
      SimplexId a{}, b{};
      inputTriangl.getEdgeVertex(e, 0, a);
      inputTriangl.getEdgeVertex(e, 1, b);

      // Ad triangle a - em - bary
      cells_offsets_[off_id] = co_id;
      off_id++;
      cells_connectivity_[co_id + 0] = a;
      cells_connectivity_[co_id + 1] = em;
      cells_connectivity_[co_id + 2] = bary;
      co_id += 3;

      // Ad triangle b - em - bary
      cells_offsets_[off_id] = co_id;
      off_id++;
      cells_connectivity_[co_id + 0] = b;
      cells_connectivity_[co_id + 1] = em;
      cells_connectivity_[co_id + 2] = bary;
      co_id += 3;
    }
  }
  // Last offset
  cells_offsets_[off_id] = co_id;

  return 0;
}

template <typename triangulationType>
int ttk::BarycentricSubdivision::execute(
  const triangulationType &inputTriangl,
  ttk::ExplicitTriangulation &outputTriangl) {

  // not implemented for dimension >= 3
  if(inputTriangl.getDimensionality() >= 3) {
    this->printErr("Not yet implemented for dimension 3 and above");
    return 1;
  }

  Timer tm;

  const SimplexId vertexNumber = inputTriangl.getNumberOfVertices();
  subdiviseTriangulation(inputTriangl);
  buildOutputTriangulation(outputTriangl);

  this->printMsg(
    "Data-set (" + std::to_string(vertexNumber) + " points) processed", 1.0,
    tm.getElapsedTime(), this->threadNumber_);

  return 0;
}
