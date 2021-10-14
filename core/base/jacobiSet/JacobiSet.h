/// \ingroup base
/// \class ttk::JacobiSet
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2015.
///
/// \brief TTK processing package for the computation of the Jacobi set of
/// bivariate volumetric data.
///
/// Given a bivariate scalar field defined on a PL 3-manifold, this package
/// produces the list of Jacobi edges (each entry is a pair given by the edge
/// identifier and the Jacobi edge type).
/// \param dataTypeU Data type of the input first component field (char, float,
/// etc.).
/// \param dataTypeV Data type of the input second component field (char, float,
/// etc.)
///
/// \b Related \b publication \n
/// "Jacobi sets of multiple Morse functions" \n
/// Herbert Edelsbrunner, John Harer \n
/// Foundations of Computational Mathematics. Cambridge University Press, 2002.
///
/// \sa ttkJacobiSet.cpp %for a usage example.

#pragma once

// base code includes
#include <Debug.h>
#include <ScalarFieldCriticalPoints.h>
#include <Triangulation.h>
#include <UnionFind.h>
#include <vector>

namespace ttk {

  class JacobiSet : virtual public Debug {
  public:
    JacobiSet();

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    int execute(std::vector<std::pair<SimplexId, char>> &jacobiSet,
                const dataTypeU *const uField,
                const dataTypeV *const vField,
                const triangulationType &triangulation,
                std::vector<char> *isPareto = NULL);

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    char getCriticalType(const SimplexId &edgeId,
                         const dataTypeU *const uField,
                         const dataTypeV *const vField,
                         const triangulationType &triangulation);

    template <class dataTypeU, class dataTypeV>
    int perturbate(const dataTypeU *const uField,
                   const dataTypeV *const vField,
                   const dataTypeU uEpsilon = Geometry::powIntTen(-DBL_DIG),
                   const dataTypeV vEpsilon
                   = Geometry::powIntTen(-DBL_DIG)) const;

    inline void
      setEdgeFans(const std::vector<std::vector<SimplexId>> *edgeFans) {
      edgeFans_ = edgeFans;
    }

    inline void setEdgeFanLinkEdgeList(
      const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
        *edgeFanLinkEdgeLists) {
      edgeFanLinkEdgeLists_ = edgeFanLinkEdgeLists;
    }

    inline void setEdgeList(
      const std::vector<std::pair<SimplexId, SimplexId>> *edgeList) {
      edgeList_ = edgeList;
    }

    inline void setSosOffsets(std::vector<SimplexId> *sosOffsets) {
      // legacy API
      setSosOffsetsU(sosOffsets->data());
    }

    /**
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p sosOffsets buffer prior to any
     * computation (the VTK wrapper already includes a mecanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    inline void setSosOffsetsU(const SimplexId *const sosOffsets) {
      sosOffsetsU_ = sosOffsets;
    }

    /**
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p sosOffsets buffer prior to any
     * computation (the VTK wrapper already includes a mecanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    inline void setSosOffsetsV(const SimplexId *const sosOffsets) {
      sosOffsetsV_ = sosOffsets;
    }

    // NOTE: here it's not clear how vtk builds vtkIdType
    // to check on bigger data-sets
    inline void setTetList(const SimplexId *tetList) {
      tetList_ = tetList;
    }

    inline void setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
    }

    inline void
      preconditionTriangulation(AbstractTriangulation *const triangulation) {
      if(triangulation) {
        triangulation->preconditionEdges();
        triangulation->preconditionEdgeStars();
      }
    }

  protected:
    template <class dataTypeU, class dataTypeV>
    int executeLegacy(std::vector<std::pair<SimplexId, char>> &jacobiSet,
                      const dataTypeU *const uField,
                      const dataTypeV *const vField);

    SimplexId vertexNumber_{};
    const SimplexId *tetList_{};
    const std::vector<std::pair<SimplexId, SimplexId>> *edgeList_{};
    // for each edge, one skeleton of its triangle fan
    const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
      *edgeFanLinkEdgeLists_{};
    // for each edge, the one skeleton of its triangle fan
    const std::vector<std::vector<SimplexId>> *edgeFans_{};
    const SimplexId *sosOffsetsU_{}, *sosOffsetsV_{};
  };
} // namespace ttk

// if the package is not a template, comment the following line
#include "JacobiSet_Template.h"
