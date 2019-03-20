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
    inline void setInputQuads(void *const address, unsigned int size) {
      inputQuads_ = reinterpret_cast<Quad *>(address);
      inputQuadNumber_ = size / 5; // 5 values per quad: 1 number and 4 indices
    }
    inline void setInputVertices(void *const address, unsigned int size) {
      inputVertices_ = reinterpret_cast<Point *>(address);
      inputVertexNumber_ = size;
    }
    inline void setInputVertexIdentifiers(void *const address) {
      inputVertexIdentifiers_ = static_cast<SimplexId *>(address);
    }
    inline void setupTriangulation(Triangulation *const triangl) {
      triangulation_ = triangl;
      if(triangulation_ != nullptr) {
        vertexNumber_ = triangulation_->getNumberOfVertices();
        triangulation_->preprocessVertexNeighbors();
      }
    }
    inline void setOutputQuads(std::vector<long long> *const quads) {
      outputQuads_ = reinterpret_cast<std::vector<Quad> *>(quads);
    }
    inline void setOutputPoints(std::vector<float> *const address) {
      outputPoints_ = reinterpret_cast<std::vector<Point> *>(address);
    }

    int execute();

  private:
    // vtkPoint instance with interleaved coordinates (AoS)
    struct Point {
      float x;
      float y;
      float z;
      Point operator+(Point other) {
        Point res;
        res.x = x + other.x;
        res.y = y + other.y;
        res.z = z + other.z;
        return res;
      }
      Point operator*(float scalar) {
        Point res;
        res.x = x * scalar;
        res.y = y * scalar;
        res.z = z * scalar;
        return res;
      }
      Point operator-(Point other) {
        return *this + other * (-1);
      }
    };
    // VTK_QUAD representation with vtkIdType
    struct Quad {
      long long n; // number of vertices, 4
      long long i; // index of first vertex
      long long j; // second vertex
      long long k; // third vertex
      long long l; // fourth vertex
    };

    int subdivise(std::vector<Quad> &currQuads,
                  const std::vector<Quad> &prevQuads);

    int project(const size_t firstPointIdx);

  protected:
    // number of vertices in the mesh
    SimplexId vertexNumber_;

    // wanted number of subdivisions of the input quadrangles
    unsigned int subdivisionLevel_;
    // number of relaxation iterations
    unsigned int relaxationIterations_;

    // number of input quadrangles
    unsigned int inputQuadNumber_;
    // input quadrangles
    Quad *inputQuads_;

    // number of input points (quad vertices)
    unsigned int inputVertexNumber_;
    // input quadrangle vertices (3D coordinates)
    Point *inputVertices_;
    // TTK identifiers of input quadrangles vertices
    SimplexId *inputVertexIdentifiers_;

    // input triangulation
    Triangulation *triangulation_;

    // array of output quadrangles
    std::vector<Quad> *outputQuads_;
    // array of output quadrangle vertices
    std::vector<Point> *outputPoints_;
  };
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <QuadrangulationSubdivision.cpp>
