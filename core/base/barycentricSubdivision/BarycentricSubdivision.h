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
    template <class dataType>
    int execute() const;

    inline void setOutputTriangulation(Triangulation *const triangulation) {
      outputTriangl_ = triangulation;
    }

    inline void setupTriangulation(Triangulation *const triangulation) {
      inputTriangl_ = triangulation;
      if(inputTriangl_ != nullptr) {
        inputTriangl_->preprocessVertexNeighbors();
        inputTriangl_->preprocessEdges();
        inputTriangl_->preprocessTriangles();
      }
    }

  private:
    int subdiviseTriangulation();
    int buildOutputTriangulation();

    // input triangulation
    Triangulation *inputTriangl_{};
    // array of input points coordinates
    const LongSimplexId *inputPoints_{};
    // list of input point data
    std::vector<void *> pointData_{};
    // list of input cell data
    std::vector<void *> cellData_{};

    // output 3D coordinates of generated points: old points first, then edge
    // middles, then triangle barycenters
    std::vector<float> points_{};
    // output triangles
    std::vector<LongSimplexId> cells_{};
    // output triangulation built on output points & output cells
    Triangulation *outputTriangl_{};
  };
} // namespace ttk

template <class dataType>
int ttk::BarycentricSubdivision::execute() const {

  Timer t;

  SimplexId vertexNumber = inputTriangl_->getNumberOfVertices();
  {
    std::stringstream msg;
    msg << "[BarycentricSubdivision] Data-set (" << vertexNumber
        << " points) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}
