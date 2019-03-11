/// \ingroup base
/// \class ttk::ScalarFieldCriticalPoints
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2015.
///
/// \brief TTK processing package for the computation of critical points in PL
/// scalar fields defined on PL manifolds.
///
/// This class computes the list of critical points of the input scalar field
/// and classify them according to their type.
///
/// \param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Critical points and curvature for embedded polyhedral surfaces" \n
/// Thomas Banchoff \n
/// American Mathematical Monthly, 1970.
///
/// \sa ttkScalarFieldCriticalPoints.cpp %for a usage example.

#ifndef _SCALARFIELDCRITICALPOINTS_H
#define _SCALARFIELDCRITICALPOINTS_H

#include <map>

// base code includes
#include <Triangulation.h>
#include <UnionFind.h>
#include <Wrapper.h>

namespace ttk {

  template <class dataType>
  class ScalarFieldCriticalPoints : public Debug {

  public:
    ScalarFieldCriticalPoints();

    ~ScalarFieldCriticalPoints();

    /// Execute the package.
    /// \param argment Dummy integer argument.
    /// \return Returns 0 upon success, negative values otherwise.
    int execute();

    std::pair<SimplexId, SimplexId>
      getNumberOfLowerUpperComponents(const SimplexId vertexId,
                                      Triangulation *triangulation) const;

    char getCriticalType(const SimplexId &vertexId) const {

      return getCriticalType(vertexId, triangulation_);
    }

    char getCriticalType(const SimplexId &vertexId,
                         Triangulation *triangulation) const;

    char getCriticalType(const SimplexId &vertexId,
                         const std::vector<std::pair<SimplexId, SimplexId>>
                           &vertexLinkEdgeList) const;

    static bool isSosHigherThan(const SimplexId &offset0,
                                const dataType &value0,
                                const SimplexId &offset1,
                                const dataType &value1) {

      return ((value0 > value1) || ((value0 == value1) && (offset0 > offset1)));
    }

    static bool isSosLowerThan(const SimplexId &offset0,
                               const dataType &value0,
                               const SimplexId &offset1,
                               const dataType &value1) {

      return ((value0 < value1) || ((value0 == value1) && (offset0 < offset1)));
    }

    int setDomainDimension(const int &dimension) {

      dimension_ = dimension;

      return 0;
    }

    int setOutput(std::vector<std::pair<SimplexId, char>> *criticalPoints) {

      criticalPoints_ = criticalPoints;

      return 0;
    }

    int setupTriangulation(Triangulation *triangulation) {

      triangulation_ = triangulation;

      // pre-condition functions
      if(triangulation_) {
        triangulation_->preprocessVertexNeighbors();
        triangulation_->preprocessVertexStars();
      }

      return 0;
    }

    int setScalarValues(const void *data) {

      scalarValues_ = (const dataType *)data;

      return 0;
    }

    int setSosOffsets(std::vector<SimplexId> *offsets) {

      sosOffsets_ = offsets;

      return 0;
    }

    int setVertexLinkEdgeLists(
      const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
        *edgeList) {

      vertexLinkEdgeLists_ = edgeList;

      return 0;
    }

    /// Set the number of vertices in the scalar field.
    /// \param vertexNumber Number of vertices in the data-set.
    /// \return Returns 0 upon success, negative values otherwise.
    int setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    void setNonManifold(const bool b) {
      forceNonManifoldCheck = b;
    }

  protected:
    int dimension_;
    SimplexId vertexNumber_;
    const dataType *scalarValues_;
    const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
      *vertexLinkEdgeLists_;
    std::vector<std::pair<SimplexId, char>> *criticalPoints_;
    std::vector<SimplexId> *sosOffsets_;
    std::vector<SimplexId> localSosOffSets_;
    Triangulation *triangulation_;

    bool forceNonManifoldCheck;
  };
} // namespace ttk

// if the package is not a template, comment the following line
#include <ScalarFieldCriticalPoints.inl>

#endif // SCALARFIELDCRITICALPOINTS_H
