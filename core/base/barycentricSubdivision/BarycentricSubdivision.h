/// \ingroup base
/// \class ttk::BarycentricSubdivision
/// \author Pierre Guillou (pierre.guillou@lip6.fr
/// \date July 2019
///
/// \brief TTK %barycentricSubdivision processing package.
///
/// %BarycentricSubdivision is a TTK processing package that takes a scalar
/// field on the input and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkBarycentricSubdivision.cpp %for a usage example.

#pragma once

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {

  class BarycentricSubdivision : public Debug {

  public:
    inline void setOutputTriangulation(Triangulation *const triangulation) {
      outputTriangl_ = triangulation;
    }
    inline void setInputPoints(const void *const addr) {
      inputPoints_ = static_cast<const float *>(addr);
    }
    inline void setupTriangulation(Triangulation *const triangulation) {
      inputTriangl_ = triangulation;
      if(inputTriangl_ != nullptr) {
        inputTriangl_->preprocessVertexNeighbors();
        inputTriangl_->preprocessEdges();
        inputTriangl_->preprocessTriangles();
      }
    }

    int execute();

  private:
    int subdiviseTriangulation();
    int buildOutputTriangulation();

    // input triangulation
    Triangulation *inputTriangl_{};
    // array of input points coordinates
    const float *inputPoints_{};
    // list of input point data
    std::vector<void *> pointData_{};
    // list of input cell data
    std::vector<void *> cellData_{};

  public:
    // output 3D coordinates of generated points: old points first, then edge
    // middles, then triangle barycenters
    std::vector<float> points_{};
    // output triangles
    std::vector<LongSimplexId> cells_{};
    // output triangulation built on output points & output cells
    Triangulation *outputTriangl_{};
    // generated points dimension: 0 vertex of parent triangulation, 1 edge
    // middle, 2 triangle barycenter
    std::vector<SimplexId> pointDim_{};
  };
} // namespace ttk
