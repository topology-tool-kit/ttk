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

    inline void setCriticalPointsNumber(unsigned int value) {
      criticalPointsNumber_ = value;
    }
    inline void setCriticalPointsCellIds(void *address) {
      criticalPointsCellIds_ = static_cast<SimplexId *>(address);
    }
    inline void setSeparatriceNumber(unsigned int value) {
      separatriceNumber_ = value;
    }
    inline void
      setSeparatrices(unsigned int number, void *sources, void *dests) {
      separatriceNumber_ = number;
      sepSourceId_ = static_cast<SimplexId *>(sources);
      sepDestId_ = static_cast<SimplexId *>(dests);
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
    // number of critical points from the Morse-Smale complex
    SimplexId criticalPointsNumber_;
    // mapping points id -> cells id
    SimplexId *criticalPointsCellIds_;

    // number of separatrices
    SimplexId separatriceNumber_;
    // unordered list of cells id sources
    SimplexId *sepSourceId_;
    // unordered list of cells id destinations
    SimplexId *sepDestId_;

    // array of output polygons
    std::vector<long long> *outputCells_;
  };
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <SurfaceQuadrangulation.cpp>
