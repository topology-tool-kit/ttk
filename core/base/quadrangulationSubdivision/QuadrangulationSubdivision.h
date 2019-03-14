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

    inline void setVertexNumber(SimplexId vertexNumber) {
      vertexNumber_ = vertexNumber;
    }
    inline void setSubdivisionLevel(unsigned int value) {
      subdivisionLevel_ = value;
    }
    inline void setRelaxationIterations(unsigned int value) {
      relaxationIterations_ = value;
    }
    inline void setCriticalPointsNumber(unsigned int value) {
      criticalPointsNumber_ = value;
    }
    inline void setCriticalPoints(void *address) {
      criticalPoints_ = static_cast<Points *>(address);
    }
    inline void setCriticalPointsIdentifiers(void *address) {
      criticalPointsIdentifiers_ = static_cast<SimplexId *>(address);
    }
    inline void setCriticalPointsCellIds(void *address) {
      criticalPointsCellIds_ = static_cast<SimplexId *>(address);
    }
    inline void setCriticalPointsType(void *address) {
      criticalPointsType_ = static_cast<SimplexId *>(address);
    }
    inline void setSeparatriceNumber(unsigned int value) {
      separatriceNumber_ = value;
    }
    inline void setSepId(void *address) {
      sepId_ = static_cast<SimplexId *>(address);
    }
    inline void setSepSourceId(void *address) {
      sepSourceId_ = static_cast<SimplexId *>(address);
    }
    inline void setSepDestId(void *address) {
      sepDestId_ = static_cast<SimplexId *>(address);
    }
    inline void setOutputCells(std::vector<long long> *cells) {
      outputCells_ = cells;
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

    // number of subdivisions of the Morse-Smale Complex cells
    unsigned int subdivisionLevel_;
    // number of relaxation iterations
    unsigned int relaxationIterations_;

    // number of critical points from the Morse-Smale complex
    SimplexId criticalPointsNumber_;
    // interleaved array of critical points coordinates
    Points *criticalPoints_;
    // critical point identifiers in the source mesh
    SimplexId *criticalPointsIdentifiers_;
    SimplexId *criticalPointsCellIds_;
    SimplexId *criticalPointsType_;

    SimplexId separatriceNumber_;
    SimplexId *sepId_;
    SimplexId *sepSourceId_;
    SimplexId *sepDestId_;

    // array of output polygons
    std::vector<long long> *outputCells_;
  };
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <QuadrangulationSubdivision.cpp>
