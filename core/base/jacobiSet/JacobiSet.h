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

#ifndef _JACOBISET_H
#define _JACOBISET_H

// base code includes
#include <ScalarFieldCriticalPoints.h>
#include <Triangulation.h>
#include <UnionFind.h>
#include <Wrapper.h>

namespace ttk {

  template <class dataTypeU, class dataTypeV>
  class JacobiSet : public Debug {

  public:
    JacobiSet();

    ~JacobiSet();

    int connectivityPreprocessing(
      const std::vector<std::vector<SimplexId>> &edgeStarList,
      std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
        &edgeFanLinkEdgeLists,
      std::vector<std::vector<LongSimplexId>> &edgeFans,
      std::vector<SimplexId> &sosOffsets) const;

    int execute(std::vector<std::pair<SimplexId, char>> &jacobiSet);

    char getCriticalType(const SimplexId &edgeId);

    int perturbate(const dataTypeU &uEpsilon = pow(10, -DBL_DIG),
                   const dataTypeV &vEpsilon = pow(10, -DBL_DIG)) const;

    int setEdgeFans(const std::vector<std::vector<SimplexId>> *edgeFans) {
      edgeFans_ = edgeFans;
      return 0;
    }

    int setEdgeFanLinkEdgeList(
      const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
        *edgeFanLinkEdgeLists) {
      edgeFanLinkEdgeLists_ = edgeFanLinkEdgeLists;
      return 0;
    }

    int setEdgeList(
      const std::vector<std::pair<SimplexId, SimplexId>> *edgeList) {
      edgeList_ = edgeList;
      return 0;
    }

    int setInputField(const void *uField, const void *vField) {

      uField_ = uField;
      vField_ = vField;
      return 0;
    }

    int setSosOffsets(std::vector<SimplexId> *sosOffsets) {
      // legacy API
      return setSosOffsetsU(sosOffsets);
    }

    int setSosOffsetsU(std::vector<SimplexId> *sosOffsets) {
      sosOffsetsU_ = sosOffsets;
      return 0;
    }

    int setSosOffsetsV(std::vector<SimplexId> *sosOffsets) {
      sosOffsetsV_ = sosOffsets;
      return 0;
    }

    // NOTE: here it's not clear how vtk builds vtkIdType
    // to check on bigger data-sets
    int setTetList(const SimplexId *tetList) {
      tetList_ = tetList;
      return 0;
    }

    int setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    int setupTriangulation(Triangulation *triangulation) {

      triangulation_ = triangulation;

      // pre-condition functions
      if(triangulation_) {
        triangulation_->preprocessEdges();
        triangulation_->preprocessEdgeStars();
      }

      return 0;
    }

  protected:
    int executeLegacy(std::vector<std::pair<SimplexId, char>> &jacobiSet);

    SimplexId vertexNumber_;
    const SimplexId *tetList_;
    const void *uField_, *vField_;
    const std::vector<std::pair<SimplexId, SimplexId>> *edgeList_;
    // for each edge, one skeleton of its triangle fan
    const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
      *edgeFanLinkEdgeLists_;
    // for each edge, the one skeleton of its triangle fan
    const std::vector<std::vector<SimplexId>> *edgeFans_;
    std::vector<SimplexId> *sosOffsetsU_, *sosOffsetsV_;
    std::vector<SimplexId> localSosOffsetsU_, localSosOffsetsV_;
    Triangulation *triangulation_;
  };
} // namespace ttk

// if the package is not a template, comment the following line
#include <JacobiSet.inl>

#endif // JACOBISET_H
