/// \ingroup base
/// \class ttk::SurfaceQuadrangulation
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date March 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkSurfaceQuadrangulation.cpp % for a usage example.

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

  class SurfaceQuadrangulation : public Debug {

  public:
    SurfaceQuadrangulation();

    // default destructor
    ~SurfaceQuadrangulation() override = default;
    // default copy constructor
    SurfaceQuadrangulation(const SurfaceQuadrangulation &) = default;
    // default move constructor
    SurfaceQuadrangulation(SurfaceQuadrangulation &&) = default;
    // default copy assignment operator
    SurfaceQuadrangulation &operator=(const SurfaceQuadrangulation &) = default;
    // default move assignment operator
    SurfaceQuadrangulation &operator=(SurfaceQuadrangulation &&) = default;

    inline int setVertexNumber(SimplexId vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }
    inline int setupTriangulation(Triangulation *triangulation) {
      triangulation_ = triangulation;
      if(triangulation_ != nullptr) {
        vertexNumber_ = triangulation_->getNumberOfVertices();
        triangulation_->preprocessVertexNeighbors();
      }
      return 0;
    }
    inline int setSubdivisionLevel(unsigned int value) {
      subdivisionLevel_ = value;
      return 0;
    }
    inline int setRelaxationIterations(unsigned int value) {
      relaxationIterations_ = value;
      return 0;
    }
    inline int setInputScalarFieldPointer(void *pointer) {
      inputScalarFieldPointer_ = pointer;
      return 0;
    }
    inline int setInputOffsetIdentifiersFieldPointer(void *pointer) {
      inputOffsetIdentifiersFieldPointer_ = pointer;
      return 0;
    }

    int execute() const;

  protected:
    // number of vertices in the mesh
    SimplexId vertexNumber_;
    // triangular input mesh
    Triangulation *triangulation_;

    // number of subdivisions of the Morse-Smale Complex cells
    unsigned int subdivisionLevel_;
    // number of relaxation iterations
    unsigned int relaxationIterations_;

    // input scalar field
    void *inputScalarFieldPointer_;
    // input offset field for identifiers
    void *inputOffsetIdentifiersFieldPointer_;
  };
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <SurfaceQuadrangulation.cpp>
