/// \ingroup base
/// \class ttk::ExplicitTriangulation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date January 2016.
///
/// \brief ExplicitTriangulation is a class that provides time efficient
/// traversal methods on triangulations of piecewise linear manifolds.
/// \sa Triangulation

#ifndef _EXPLICITTRIANGULATION_H
#define _EXPLICITTRIANGULATION_H

// base code includes
#include <AbstractTriangulation.h>
#include <OneSkeleton.h>
#include <ThreeSkeleton.h>
#include <TwoSkeleton.h>
#include <ZeroSkeleton.h>

#include <memory>

namespace ttk {

  class ExplicitTriangulation final : public AbstractTriangulation {

  public:
    ExplicitTriangulation();

    virtual ~ExplicitTriangulation();

    int clear();

    size_t footprint(size_t size = 0) const;

    inline int getCellEdgeInternal(const SimplexId &cellId,
                                   const int &localEdgeId,
                                   SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)tetraEdgeList_.size()))
        return -1;
      if((localEdgeId < 0)
         || (localEdgeId >= (SimplexId)tetraEdgeList_[cellId].size()))
        return -2;
#endif
      edgeId = tetraEdgeList_[cellId][localEdgeId];
      return 0;
    }

    inline SimplexId
      getCellEdgeNumberInternal(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)tetraEdgeList_.size()))
        return -1;
#endif
      return tetraEdgeList_[cellId].size();
    }

    template <std::size_t N>
    inline void
      convertToVector(const std::vector<std::array<SimplexId, N>> &table,
                      std::vector<std::vector<SimplexId>> &vec) {
      for(size_t i = 0; i < table.size(); ++i) {
        vec[i] = {table[i].begin(), table[i].end()};
      }
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellEdgesInternal() override {

      convertToVector(tetraEdgeList_, cellEdgeVector_);
      return &cellEdgeVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
      const SimplexId &cellId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {

      neighborId = cellNeighborData_.get(cellId, localNeighborId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellNeighborNumber)(
      const SimplexId &cellId) const override {
      return cellNeighborData_.size(cellId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() override {
      cellNeighborData_.copyTo(cellNeighborList_);
      return &cellNeighborList_;
    }

    inline int getCellTriangleInternal(const SimplexId &cellId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)tetraTriangleList_.size()))
        return -1;
      if((localTriangleId < 0)
         || (localTriangleId >= (SimplexId)tetraTriangleList_[cellId].size()))
        return -2;
#endif
      triangleId = tetraTriangleList_[cellId][localTriangleId];

      return 0;
    }

    inline SimplexId
      getCellTriangleNumberInternal(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)tetraTriangleList_.size()))
        return -1;
#endif

      return tetraTriangleList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellTrianglesInternal() override {

      convertToVector(tetraTriangleList_, cellTriangleVector_);
      return &cellTriangleVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellVertex)(
      const SimplexId &cellId,
      const int &localVertexId,
      SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((!cellArray_) || (!cellNumber_))
        return -1;
#endif
      vertexId = cellArray_->getCellVertex(cellId, localVertexId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellVertexNumber)(
      const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((!cellArray_) || (!cellNumber_))
        return -1;
#endif
      return cellArray_->getCellVertexNumber(cellId);
    }

    int TTK_TRIANGULATION_INTERNAL(getDimensionality)() const override {
      return maxCellDim_;
    }

    inline const std::vector<std::array<SimplexId, 2>> *
      TTK_TRIANGULATION_INTERNAL(getEdges)() override {
      return &edgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeLink)(
      const SimplexId &edgeId,
      const int &localLinkId,
      SimplexId &linkId) const override {

      linkId = edgeLinkData_.get(edgeId, localLinkId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeLinkNumber)(
      const SimplexId &edgeId) const override {
      return edgeLinkData_.size(edgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() override {
      edgeLinkData_.copyTo(edgeLinkList_);
      return &edgeLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeStar)(
      const SimplexId &edgeId,
      const int &localStarId,
      SimplexId &starId) const override {

      starId = edgeStarData_.get(edgeId, localStarId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(
      const SimplexId &edgeId) const override {
      return edgeStarData_.size(edgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeStars)() override {
      edgeStarData_.copyTo(edgeStarList_);
      return &edgeStarList_;
    }

    inline int getEdgeTriangleInternal(const SimplexId &edgeId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {
      triangleId = edgeTriangleData_.get(edgeId, localTriangleId);
      return 0;
    }

    inline SimplexId
      getEdgeTriangleNumberInternal(const SimplexId &edgeId) const override {
      return edgeTriangleData_.size(edgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTrianglesInternal() override {
      edgeTriangleData_.copyTo(edgeTriangleList_);
      return &edgeTriangleList_;
    }

    inline int getEdgeVertexInternal(const SimplexId &edgeId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeList_.size()))
        return -1;
      if((localVertexId != 0) && (localVertexId != 1))
        return -2;
#endif
      if(!localVertexId)
        vertexId = edgeList_[edgeId][0];
      else
        vertexId = edgeList_[edgeId][1];
      return 0;
    }

    inline SimplexId
      TTK_TRIANGULATION_INTERNAL(getNumberOfCells)() const override {
      return cellNumber_;
    }

    inline SimplexId getNumberOfEdgesInternal() const override {
      return edgeList_.size();
    }

    inline SimplexId getNumberOfTrianglesInternal() const override {
      return triangleList_.size();
    }

    inline SimplexId
      TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() const override {
      return vertexNumber_;
    }

    inline const std::vector<std::array<SimplexId, 3>> *
      TTK_TRIANGULATION_INTERNAL(getTriangles)() override {
      return &triangleList_;
    }

    inline int getTriangleEdgeInternal(const SimplexId &triangleId,
                                       const int &localEdgeId,
                                       SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleEdgeList_.size()))
        return -1;
      if((localEdgeId < 0) || (localEdgeId > 2))
        return -2;
#endif

      edgeId = triangleEdgeList_[triangleId][localEdgeId];

      return 0;
    }

    inline SimplexId getTriangleEdgeNumberInternal(
      const SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleEdgeList_.size()))
        return -1;
#endif

      return triangleEdgeList_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getTriangleEdgesInternal() override {

      convertToVector(triangleEdgeList_, triangleEdgeVector_);
      return &triangleEdgeVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
      const SimplexId &triangleId,
      const int &localLinkId,
      SimplexId &linkId) const override {

      linkId = triangleLinkData_.get(triangleId, localLinkId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleLinkNumber)(
      const SimplexId &triangleId) const override {
      return triangleLinkData_.size(triangleId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() override {
      triangleLinkData_.copyTo(triangleLinkList_);
      return &triangleLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
      const SimplexId &triangleId,
      const int &localStarId,
      SimplexId &starId) const override {
      starId = triangleStarData_.get(triangleId, localStarId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(
      const SimplexId &triangleId) const override {
      return triangleStarData_.size(triangleId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleStars)() override {
      triangleStarData_.copyTo(triangleStarList_);
      return &triangleStarList_;
    }

    inline int getTriangleVertexInternal(const SimplexId &triangleId,
                                         const int &localVertexId,
                                         SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId >= (SimplexId)triangleList_.size()))
        return -1;
      if((localVertexId < 0)
         || (localVertexId >= (SimplexId)triangleList_[triangleId].size()))
        return -2;
#endif
      vertexId = triangleList_[triangleId][localVertexId];
      return 0;
    }

    inline int getVertexEdgeInternal(const SimplexId &vertexId,
                                     const int &localEdgeId,
                                     SimplexId &edgeId) const override {
      edgeId = vertexEdgeData_.get(vertexId, localEdgeId);
      return 0;
    }

    inline SimplexId
      getVertexEdgeNumberInternal(const SimplexId &vertexId) const override {
      return vertexEdgeData_.size(vertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexEdgesInternal() override {
      vertexEdgeData_.copyTo(vertexEdgeList_);
      return &vertexEdgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexLink)(
      const SimplexId &vertexId,
      const int &localLinkId,
      SimplexId &linkId) const override {

      linkId = vertexLinkData_.get(vertexId, localLinkId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexLinkNumber)(
      const SimplexId &vertexId) const override {
      return vertexLinkData_.size(vertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexLinks)() override {
      vertexLinkData_.copyTo(vertexLinkList_);
      return &vertexLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexNeighbor)(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {

      neighborId = vertexNeighborData_.get(vertexId, localNeighborId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexNeighborNumber)(
      const SimplexId &vertexId) const override {
      return vertexNeighborData_.size(vertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() override {
      vertexNeighborData_.copyTo(vertexNeighborList_);
      return &vertexNeighborList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexPoint)(
      const SimplexId &vertexId, float &x, float &y, float &z) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      if(doublePrecision_) {
        x = ((double *)pointSet_)[3 * vertexId];
        y = ((double *)pointSet_)[3 * vertexId + 1];
        z = ((double *)pointSet_)[3 * vertexId + 2];
      } else {
        x = ((float *)pointSet_)[3 * vertexId];
        y = ((float *)pointSet_)[3 * vertexId + 1];
        z = ((float *)pointSet_)[3 * vertexId + 2];
      }

      return 0;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexStar)(
      const SimplexId &vertexId,
      const int &localStarId,
      SimplexId &starId) const override {
      starId = vertexStarData_.get(vertexId, localStarId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexStarNumber)(
      const SimplexId &vertexId) const override {
      return vertexStarData_.size(vertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexStars)() override {
      vertexStarData_.copyTo(vertexStarList_);
      return &vertexStarList_;
    }

    inline int getVertexTriangleInternal(const SimplexId &vertexId,
                                         const int &localTriangleId,
                                         SimplexId &triangleId) const override {
      triangleId = vertexTriangleData_.get(vertexId, localTriangleId);
      return 0;
    }

    inline SimplexId getVertexTriangleNumberInternal(
      const SimplexId &vertexId) const override {
      return vertexTriangleData_.size(vertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() override {
      vertexTriangleData_.copyTo(vertexTriangleList_);
      return &vertexTriangleList_;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isEdgeOnBoundary)(
      const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)boundaryEdges_.size()))
        return false;
#endif
      return boundaryEdges_[edgeId];
    }

    inline bool isEmpty() const override {
      return !vertexNumber_;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isTriangleOnBoundary)(
      const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)boundaryTriangles_.size()))
        return false;
#endif
      return boundaryTriangles_[triangleId];
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isVertexOnBoundary)(
      const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)boundaryVertices_.size()))
        return false;
#endif
      return boundaryVertices_[vertexId];
    }

    inline int preconditionBoundaryEdgesInternal() override {

      if((!boundaryEdges_.empty())
         && (boundaryEdges_.size() == edgeList_.size())) {
        return 0;
      }

      preconditionEdgesInternal();
      boundaryEdges_.resize(edgeList_.size(), false);

      if(getDimensionality() == 2) {
        preconditionEdgeStarsInternal();
        for(SimplexId i = 0; i < (SimplexId)edgeStarData_.subvectorsNumber();
            i++) {
          if(edgeStarData_.size(i) == 1) {
            boundaryEdges_[i] = true;
          }
        }
      } else if(getDimensionality() == 3) {
        preconditionTriangleStarsInternal();
        preconditionTriangleEdgesInternal();

        for(size_t i = 0; i < triangleStarData_.subvectorsNumber(); i++) {
          if(triangleStarData_.size(i) == 1) {
            for(int j = 0; j < 3; j++) {
              boundaryEdges_[triangleEdgeList_[i][j]] = true;
            }
          }
        }
      } else {
        // unsupported dimension
        printErr("Unsupported dimension for boundary precondition");
        return -1;
      }

      return 0;
    }

    inline int preconditionBoundaryTrianglesInternal() override {

      if(getDimensionality() == 2)
        return 0;

      if((!boundaryTriangles_.empty())
         && (boundaryTriangles_.size() == triangleList_.size())) {
        return 0;
      }

      preconditionTrianglesInternal();
      boundaryTriangles_.resize(triangleList_.size(), false);

      if(getDimensionality() == 3) {
        preconditionTriangleStarsInternal();

        for(size_t i = 0; i < triangleStarData_.subvectorsNumber(); i++) {
          if(triangleStarData_.size(i) == 1) {
            boundaryTriangles_[i] = true;
          }
        }
      } else {
        // unsupported dimension
        printErr("Unsupported dimension for boundary precondition");
        return -1;
      }

      return 0;
    }

    inline int preconditionBoundaryVerticesInternal() override {

      if((!boundaryVertices_.empty())
         && ((SimplexId)boundaryVertices_.size() == vertexNumber_))
        return 0;

      boundaryVertices_.resize(vertexNumber_, false);

      // create the list of boundary elements
      // create their star
      // look for singletons
      if(getDimensionality() == 1) {
        preconditionVertexStarsInternal();
        for(size_t i = 0; i < vertexStarData_.subvectorsNumber(); i++) {
          if(vertexStarData_.size(i) == 1) {
            boundaryVertices_[i] = true;
          }
        }
      } else if(getDimensionality() == 2) {
        preconditionEdgesInternal();
        preconditionEdgeStarsInternal();

        for(SimplexId i = 0; i < (SimplexId)edgeStarData_.subvectorsNumber();
            i++) {
          if(edgeStarData_.size(i) == 1) {
            boundaryVertices_[edgeList_[i][0]] = true;
            boundaryVertices_[edgeList_[i][1]] = true;
          }
        }
      } else if(getDimensionality() == 3) {
        preconditionTrianglesInternal();
        preconditionTriangleStarsInternal();

        for(size_t i = 0; i < triangleStarData_.subvectorsNumber(); i++) {
          if(triangleStarData_.size(i) == 1) {
            boundaryVertices_[triangleList_[i][0]] = true;
            boundaryVertices_[triangleList_[i][1]] = true;
            boundaryVertices_[triangleList_[i][2]] = true;
          }
        }
      } else {
        // unsupported dimension
        printErr("Unsupported dimension for boundary precondition");
        return -1;
      }

      return 0;
    }

    inline int preconditionCellEdgesInternal() override {

      OneSkeleton os;
      os.setWrapper(this);

      if(tetraEdgeList_.empty() && getDimensionality() == 3) {
        os.buildEdgeList(
          vertexNumber_, *cellArray_, nullptr, nullptr, &tetraEdgeList_);
      } else if(triangleEdgeList_.empty() && getDimensionality() == 2) {
        os.buildEdgeList(
          vertexNumber_, *cellArray_, nullptr, nullptr, &triangleEdgeList_);
      }

      return 0;
    }

    inline int preconditionCellNeighborsInternal() override {

      if(cellNeighborData_.empty()) {
        if(getDimensionality() == 3) {
          ThreeSkeleton threeSkeleton;
          threeSkeleton.setWrapper(this);
          threeSkeleton.buildCellNeighborsFromTriangles(
            vertexNumber_, *cellArray_, cellNeighborData_, &triangleStarData_);
        } else if(getDimensionality() == 2) {
          TwoSkeleton twoSkeleton;
          twoSkeleton.setWrapper(this);
          twoSkeleton.buildCellNeighborsFromEdges(
            vertexNumber_, *cellArray_, cellNeighborData_, &edgeStarData_);
        }
      }

      return 0;
    }

    inline int preconditionCellTrianglesInternal() override {

      if(!tetraTriangleList_.size()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        if(triangleList_.size()) {
          // we already computed this guy, let's just get the cell triangles
          if(!triangleStarData_.empty()) {
            return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                                 nullptr, nullptr,
                                                 &tetraTriangleList_);
          } else {
            // let's compute the triangle star while we're at it...
            // it's just a tiny overhead.
            return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                                 nullptr, &triangleStarData_,
                                                 &tetraTriangleList_);
          }
        } else {
          // we have not computed this guy, let's do it while we're at it
          if(!triangleStarData_.empty()) {
            return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                                 &triangleList_, nullptr,
                                                 &tetraTriangleList_);
          } else {
            // let's compute the triangle star while we're at it...
            // it's just a tiny overhead.
            return twoSkeleton.buildTriangleList(
              vertexNumber_, *cellArray_, &triangleList_, &triangleStarData_,
              &tetraTriangleList_);
          }
        }
      }

      return 0;
    }

    inline int preconditionEdgesInternal() override {

      if(!edgeList_.size()) {
        OneSkeleton oneSkeleton;
        oneSkeleton.setWrapper(this);
        // also computes edgeStar and triangleEdge / tetraEdge lists for free...
        if(getDimensionality() == 2) {
          return oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_,
                                           &edgeList_, &edgeStarData_,
                                           &triangleEdgeList_);
        } else if(getDimensionality() == 3) {
          return oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_,
                                           &edgeList_, &edgeStarData_,
                                           &tetraEdgeList_);
        }
      }

      return 0;
    }

    inline int preconditionEdgeLinksInternal() override {

      if(edgeLinkData_.empty()) {

        if(getDimensionality() == 2) {
          preconditionEdgesInternal();
          preconditionEdgeStarsInternal();

          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          return oneSkeleton.buildEdgeLinks(
            edgeList_, edgeStarData_, *cellArray_, edgeLinkData_);
        } else if(getDimensionality() == 3) {
          preconditionEdgesInternal();
          preconditionEdgeStarsInternal();
          preconditionCellEdgesInternal();

          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          return oneSkeleton.buildEdgeLinks(
            edgeList_, edgeStarData_, tetraEdgeList_, edgeLinkData_);
        } else {
          // unsupported dimension
          printErr("Unsupported dimension for edge link precondition");
          return -1;
        }
      }

      return 0;
    }

    inline int preconditionEdgeStarsInternal() override {

      if(edgeStarData_.empty()) {
        OneSkeleton oneSkeleton;
        oneSkeleton.setWrapper(this);
        return oneSkeleton.buildEdgeList<3>(
          vertexNumber_, *cellArray_, nullptr, &edgeStarData_, nullptr);
      }
      return 0;
    }

    inline int preconditionEdgeTrianglesInternal() override {

      if(edgeTriangleData_.empty()) {
        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildEdgeTriangles(vertexNumber_, *cellArray_,
                                              edgeTriangleData_, &edgeList_,
                                              &triangleEdgeList_);
      }

      return 0;
    }

    inline int preconditionTrianglesInternal() override {

      if(!triangleList_.size()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                      &triangleList_, &triangleStarData_,
                                      &tetraTriangleList_);
      }

      return 0;
    }

    inline int preconditionTriangleEdgesInternal() override {

      if(!triangleEdgeList_.size()) {

        // WARNING
        // here triangleStarList and cellTriangleList will be computed (for
        // free) although they are not requireed to get the edgeTriangleList.
        // if memory usage is an issue, please change these pointers by nullptr.

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        return twoSkeleton.buildTriangleEdgeList(
          vertexNumber_, *cellArray_, triangleEdgeList_, &vertexEdgeData_,
          &edgeList_, &triangleList_, &triangleStarData_, &tetraTriangleList_);
      }

      return 0;
    }

    inline int preconditionTriangleLinksInternal() override {

      if(triangleLinkData_.empty()) {

        preconditionTriangleStarsInternal();

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildTriangleLinks(
          triangleList_, triangleStarData_, *cellArray_, triangleLinkData_);
      }

      return 0;
    }

    inline int preconditionTriangleStarsInternal() override {

      if(triangleStarData_.empty()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildTriangleList(
          vertexNumber_, *cellArray_, &triangleList_, &triangleStarData_);
      }

      return 0;
    }

    inline int preconditionVertexEdgesInternal() override {

      if((SimplexId)vertexEdgeData_.subvectorsNumber() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;

        if(!edgeList_.size()) {
          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_, &edgeList_);
        }

        zeroSkeleton.setWrapper(this);
        return zeroSkeleton.buildVertexEdges(
          vertexNumber_, edgeList_, vertexEdgeData_);
      }
      return 0;
    }

    inline int preconditionVertexLinksInternal() override {

      if((SimplexId)vertexLinkData_.subvectorsNumber() != vertexNumber_) {

        if(getDimensionality() == 2) {
          preconditionVertexStarsInternal();
          preconditionCellEdgesInternal();

          ZeroSkeleton zeroSkeleton;
          zeroSkeleton.setWrapper(this);
          return zeroSkeleton.buildVertexLinks(
            vertexStarData_, triangleEdgeList_, edgeList_, vertexLinkData_);
        } else if(getDimensionality() == 3) {
          preconditionVertexStarsInternal();
          preconditionCellTrianglesInternal();

          ZeroSkeleton zeroSkeleton;
          zeroSkeleton.setWrapper(this);
          return zeroSkeleton.buildVertexLinks(vertexStarData_,
                                               tetraTriangleList_,
                                               triangleList_, vertexLinkData_);
        } else {
          // unsupported dimension
          printErr("Unsupported dimension for vertex link precondition");
          return -1;
        }
      }
      return 0;
    }

    inline int preconditionVertexNeighborsInternal() override {

      if((SimplexId)vertexNeighborData_.subvectorsNumber() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;
        zeroSkeleton.setWrapper(this);
        return zeroSkeleton.buildVertexNeighbors(
          vertexNumber_, *cellArray_, vertexNeighborData_, &edgeList_);
      }
      return 0;
    }

    inline int preconditionVertexStarsInternal() override {

      if((SimplexId)vertexStarData_.subvectorsNumber() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;
        zeroSkeleton.setWrapper(this);

        return zeroSkeleton.buildVertexStars(
          vertexNumber_, *cellArray_, vertexStarData_);
      }
      return 0;
    }

    inline int preconditionVertexTrianglesInternal() override {

      if((SimplexId)vertexTriangleData_.subvectorsNumber() != vertexNumber_) {

        preconditionTrianglesInternal();

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        twoSkeleton.buildVertexTriangles(
          vertexNumber_, triangleList_, vertexTriangleData_);
      }

      return 0;
    }

#ifdef TTK_CELL_ARRAY_NEW
    // Layout with connectivity + offset array (new)
    inline int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *connectivity,
                             const LongSimplexId *offset) {

      // Cell Check
      {
        if(cellNumber > 0) {
          const auto &cellDimension = offset[1] - offset[0] - 1;

          if(cellDimension < 0 || cellDimension > 3) {
            this->printErr("Unable to create triangulation for cells of "
                           "dimension 4 or higher ("
                           + std::to_string(cellDimension) + ").");
            return -1;
          }

          bool error = false;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
          for(SimplexId i = 0; i < cellNumber; i++) {
            if(offset[i + 1] - offset[i] - 1 != cellDimension) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif // TTK_ENABLE_OPENMP
              error = true;
            }
          }

          if(error) {
            this->printErr("Unable to create triangulation for "
                           "inhomogeneous\ncell dimensions.");
            return -2;
          }
        }
      }

      if(cellNumber_)
        clear();

      cellNumber_ = cellNumber;

      cellArray_
        = std::make_shared<CellArray>(connectivity, offset, cellNumber);

      // TODO: ASSUME Regular Mesh Here to compute dimension!
      if(cellNumber) {
        if(cellArray_->getCellVertexNumber(0) == 3) {
          maxCellDim_ = 2;
        } else {
          maxCellDim_ = 3;
        }
      }
      return 0;
    }
#else
    // Flat layout with a single array (legacy & default one)
    inline int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *cellArray) {
      if(cellNumber_)
        clear();

      cellNumber_ = cellNumber;

      if(cellNumber) {
        // assume regular mesh here to compute dimension
        cellArray_ = std::make_shared<CellArray>(
          cellArray, cellNumber, cellArray[0] - 1);
        maxCellDim_ = cellArray[0] - 1;
      }
      return 0;
    }
#endif

    inline int setInputPoints(const SimplexId &pointNumber,
                              const void *pointSet,
                              const bool &doublePrecision = false) {

      if(vertexNumber_)
        clear();

      vertexNumber_ = pointNumber;
      pointSet_ = pointSet;
      doublePrecision_ = doublePrecision;
      return 0;
    }

  private:
    bool doublePrecision_;
    SimplexId cellNumber_, vertexNumber_;
    const void *pointSet_;
    int maxCellDim_;
    std::shared_ptr<CellArray> cellArray_;

    FlatJaggedArray vertexNeighborData_{};
    FlatJaggedArray cellNeighborData_{};
    FlatJaggedArray vertexEdgeData_{};
    FlatJaggedArray vertexTriangleData_{};
    FlatJaggedArray edgeTriangleData_{};
    FlatJaggedArray vertexStarData_{};
    FlatJaggedArray edgeStarData_{};
    FlatJaggedArray triangleStarData_{};
    FlatJaggedArray vertexLinkData_{};
    FlatJaggedArray edgeLinkData_{};
    FlatJaggedArray triangleLinkData_{};
  };
} // namespace ttk

#endif // _EXPLICITTRIANGULATION_H
