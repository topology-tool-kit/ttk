/// \ingroup base
/// \class ttk::QuadrangulationSubdivision
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date March 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkQuadrangulationSubdivision.cpp % for a usage example.

#pragma once

// base code includes
#include <Geometry.h>
#include <MorseSmaleComplex.h>
#include <Triangulation.h>
#include <Wrapper.h>
#include <cmath>
#include <set>
#include <tuple>
#include <type_traits>

namespace ttk {

  class QuadrangulationSubdivision : public Debug {

  public:
    QuadrangulationSubdivision();

    // default destructor
    ~QuadrangulationSubdivision() override = default;
    // default copy constructor
    QuadrangulationSubdivision(const QuadrangulationSubdivision &) = default;
    // default move constructor
    QuadrangulationSubdivision(QuadrangulationSubdivision &&) = default;
    // default copy assignment operator
    QuadrangulationSubdivision &operator=(const QuadrangulationSubdivision &)
      = default;
    // default move assignment operator
    QuadrangulationSubdivision &operator=(QuadrangulationSubdivision &&)
      = default;

    inline void setSubdivisionLevel(const unsigned int value) {
      subdivisionLevel_ = value;
    }
    inline void setRelaxationIterations(const unsigned int value) {
      relaxationIterations_ = value;
    }
    inline void setInputQuadranglesNumber(const unsigned int value) {
      inputQuadVertexNumber_ = value;
    }
    inline void setInputQuadrangles(void *const address) {
      inputQuadrangles_ = static_cast<long long *>(address);
    }
    inline void setInputQuadIdentifiers(void *const address) {
      inputQuadIdentifiers_ = static_cast<SimplexId *>(address);
    }
    inline void setupTriangulation(Triangulation *const triangl) {
      triangulation_ = triangl;
      if(triangulation_ != nullptr) {
        vertexNumber_ = triangulation_->getNumberOfVertices();
        triangulation_->preprocessVertexNeighbors();
      }
    }
    inline void setOutputQuads(std::vector<long long> *const quads) {
      outputQuads_ = quads;
    }
    inline void setOutputPointNumber(unsigned int value) {
      outputPointNumber_ = value;
    }
    inline void setOutputPoints(void *const address) {
      outputPoints_ = static_cast<Points *>(address);
    }

    int execute() const;

  private:
    // vtkPoint instance with interleaved coordinates (AoS)
    struct Points {
      float x;
      float y;
      float z;
    };

  protected:
    // number of vertices in the mesh
    SimplexId vertexNumber_;

    // wanted number of subdivisions of the input quadrangles
    unsigned int subdivisionLevel_;
    // number of relaxation iterations
    unsigned int relaxationIterations_;

    // number of vertices in the input quadrangles
    unsigned int inputQuadVertexNumber_;
    // input quadrangles
    long long *inputQuadrangles_;
    // TTK identifiers of input quadrangles vertices
    SimplexId *inputQuadIdentifiers_;
    // input triangulation
    Triangulation *triangulation_;

    // array of output quadrangles
    std::vector<long long> *outputQuads_;
    unsigned int outputPointNumber_;
    // array of output quadrangle vertices
    Points *outputPoints_;
  };
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <QuadrangulationSubdivision.cpp>
