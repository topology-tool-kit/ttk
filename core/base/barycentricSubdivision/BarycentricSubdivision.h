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

#include <type_traits>

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {

  class BarycentricSubdivision : public Debug {

  public:
    BarycentricSubdivision(std::vector<float> &points,
                           std::vector<LongSimplexId> &cells_co,
                           std::vector<LongSimplexId> &cells_off,
                           std::vector<SimplexId> &pointId,
                           std::vector<SimplexId> &pointDim)
      : points_{points}, cells_connectivity_{cells_co},
        cells_offsets_{cells_off}, pointId_{pointId}, pointDim_{pointDim} {
    }

    inline void setOutputTriangulation(Triangulation *const triangulation) {
      outputTriangl_ = triangulation;
    }
    inline void setInputPoints(const void *const addr) {
      inputPoints_ = static_cast<const float *>(addr);
    }
    inline void setupTriangulation(Triangulation *const triangulation) {
      inputTriangl_ = triangulation;
      if(inputTriangl_ != nullptr) {
        inputTriangl_->preconditionVertexNeighbors();
        inputTriangl_->preconditionEdges();
        inputTriangl_->preconditionTriangles();
        inputTriangl_->preconditionTriangleEdges();
      }
      nVertices_ = inputTriangl_->getNumberOfVertices();
      nEdges_ = inputTriangl_->getNumberOfEdges();
      nTriangles_ = inputTriangl_->getNumberOfTriangles();
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

    int execute();

    /**
     * @brief Interpolate floating-point point data on subdivised triangulation
     *
     * Copy values on parent vertices, interpolate on edges and barycenters
     *
     * @param[in] data Pointer to input data on parent triangulation
     * @param[out] output Allocated buffer to be filled
     *
     * @return 0 in case of success
     */
    template <typename T>
    int interpolateContinuousScalarField(const T *data, T *output) const {
      static_assert(
        std::is_floating_point<T>::value, "Floating point type required.");
      if(inputTriangl_ == nullptr || outputTriangl_ == nullptr) {
        return 1;
      }
      const auto nOutVerts = this->getNumberOfVertices();
      if(nOutVerts < 0 || nOutVerts != nVertices_ + nEdges_ + nTriangles_) {
        return 1;
      }

      // copy data on parent vertices
      std::copy(data, data + nVertices_, output);

      // interpolate on edges
      for(SimplexId i = 0; i < nEdges_; ++i) {
        SimplexId a{}, b{};
        inputTriangl_->getEdgeVertex(i, 0, a);
        inputTriangl_->getEdgeVertex(i, 0, b);
        output[nVertices_ + i] = (data[a] + data[b]) / T{2.0};
      }

      // interpolate on triangle barycenters
      for(SimplexId i = 0; i < nTriangles_; ++i) {
        SimplexId a{}, b{}, c{};
        inputTriangl_->getTriangleVertex(i, 0, a);
        inputTriangl_->getTriangleVertex(i, 1, b);
        inputTriangl_->getTriangleVertex(i, 2, c);
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
      if(inputTriangl_ == nullptr || outputTriangl_ == nullptr) {
        return 1;
      }
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
      if(inputTriangl_ == nullptr || outputTriangl_ == nullptr) {
        return 1;
      }
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
    int subdiviseTriangulation();
    int buildOutputTriangulation();

    // input triangulation properties
    SimplexId nVertices_{};
    SimplexId nEdges_{};
    SimplexId nTriangles_{};

    // input triangulation
    Triangulation *inputTriangl_{};
    // array of input points coordinates
    const float *inputPoints_{};

    // output 3D coordinates of generated points: old points first, then edge
    // middles, then triangle barycenters
    std::vector<float> &points_;
    // output triangles
    std::vector<LongSimplexId> &cells_connectivity_;
    std::vector<LongSimplexId> &cells_offsets_;
    // generated point cell id
    std::vector<SimplexId> &pointId_;
    // generated points dimension: 0 vertex of parent triangulation, 1 edge
    // middle, 2 triangle barycenter
    std::vector<SimplexId> &pointDim_;

    // output triangulation built on output points & output cells
    Triangulation *outputTriangl_{};
  };
} // namespace ttk
