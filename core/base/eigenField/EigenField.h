/// \ingroup base
/// \class ttk::EigenField
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkEigenField.cpp % for a usage example.

#pragma once

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {

  class EigenField : public Debug {

  public:
    // default constructor
    EigenField() = default;
    // default destructor
    ~EigenField() override = default;
    // default copy constructor
    EigenField(const EigenField &) = default;
    // default move constructor
    EigenField(EigenField &&) = default;
    // default copy assignment operator
    EigenField &operator=(const EigenField &) = default;
    // default move assignment operator
    EigenField &operator=(EigenField &&) = default;

    inline void setupTriangulation(Triangulation *triangulation) {
      triangulation_ = triangulation;
      if(triangulation_ != nullptr) {
        vertexNumber_ = triangulation_->getNumberOfVertices();
        triangulation_->preprocessVertexNeighbors();
        // cotan weights method needs more pre-processing
        triangulation_->preprocessEdgeTriangles();
      }
    }
    inline void setOutputScalarFieldPointer(void *data) {
      outputScalarFieldPointer_ = data;
    }
    inline void setEigenNumber(unsigned int eigenNumber) {
      eigenNumber_ = eigenNumber;
    }

    template <typename T>
    int execute() const;

  private:
    // the mesh
    Triangulation *triangulation_{};
    // number of vertices in the mesh
    SimplexId vertexNumber_{};
    // output eigenfunctions scalar fields
    void *outputScalarFieldPointer_{};
    // number of eigenfunctions to compute
    unsigned int eigenNumber_{20};
  };

} // namespace ttk
