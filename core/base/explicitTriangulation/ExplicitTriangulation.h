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

    inline int getCellEdgeInternal(const SimplexId &cellId,
                                   const int &localEdgeId,
                                   SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellEdgeList_.size()))
        return -1;
      if((localEdgeId < 0)
         || (localEdgeId >= (SimplexId)cellEdgeList_[cellId].size()))
        return -2;
#endif
      edgeId = cellEdgeList_[cellId][localEdgeId];
      return 0;
    }

    inline SimplexId
      getCellEdgeNumberInternal(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellEdgeList_.size()))
        return -1;
#endif
      return cellEdgeList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellEdgesInternal() override {

      return &cellEdgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
      const SimplexId &cellId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellNeighborList_.size()))
        return -1;
      if((localNeighborId < 0)
         || (localNeighborId >= (SimplexId)cellNeighborList_[cellId].size()))
        return -2;
#endif
      neighborId = cellNeighborList_[cellId][localNeighborId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellNeighborNumber)(
      const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellNeighborList_.size()))
        return -1;
#endif
      return cellNeighborList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() override {
      return &cellNeighborList_;
    }

    inline int getCellTriangleInternal(const SimplexId &cellId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellTriangleList_.size()))
        return -1;
      if((localTriangleId < 0)
         || (localTriangleId >= (SimplexId)cellTriangleList_[cellId].size()))
        return -2;
#endif
      triangleId = cellTriangleList_[cellId][localTriangleId];

      return 0;
    }

    inline SimplexId
      getCellTriangleNumberInternal(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellTriangleList_.size()))
        return -1;
#endif

      return cellTriangleList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellTrianglesInternal() override {

      return &cellTriangleList_;
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

    inline const std::vector<std::pair<SimplexId, SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdges)() override {
      return &edgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeLink)(
      const SimplexId &edgeId,
      const int &localLinkId,
      SimplexId &linkId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeLinkList_.size()))
        return -1;
      if((localLinkId < 0)
         || (localLinkId >= (SimplexId)edgeLinkList_[edgeId].size()))
        return -2;
#endif
      linkId = edgeLinkList_[edgeId][localLinkId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeLinkNumber)(
      const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeLinkList_.size()))
        return -1;
#endif
      return edgeLinkList_[edgeId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() override {
      return &edgeLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeStar)(
      const SimplexId &edgeId,
      const int &localStarId,
      SimplexId &starId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeStarList_.size()))
        return -1;
      if((localStarId < 0)
         || (localStarId >= (SimplexId)edgeStarList_[edgeId].size()))
        return -2;
#endif
      starId = edgeStarList_[edgeId][localStarId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(
      const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeStarList_.size()))
        return -1;
#endif
      return edgeStarList_[edgeId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeStars)() override {
      return &edgeStarList_;
    }

    inline int getEdgeTriangleInternal(const SimplexId &edgeId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeTriangleList_.size()))
        return -1;
      if((localTriangleId < 0)
         || (localTriangleId >= (SimplexId)edgeTriangleList_[edgeId].size()))
        return -2;
#endif

      triangleId = edgeTriangleList_[edgeId][localTriangleId];

      return 0;
    }

    inline SimplexId
      getEdgeTriangleNumberInternal(const SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeTriangleList_.size()))
        return -1;
#endif

      return edgeTriangleList_[edgeId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTrianglesInternal() override {

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
        vertexId = edgeList_[edgeId].first;
      else
        vertexId = edgeList_[edgeId].second;
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

    inline const std::vector<std::vector<SimplexId>> *
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

      return &triangleEdgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
      const SimplexId &triangleId,
      const int &localLinkId,
      SimplexId &linkId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleLinkList_.size()))
        return -1;
      if((localLinkId < 0)
         || (localLinkId >= (SimplexId)triangleLinkList_[triangleId].size()))
        return -2;
#endif
      linkId = triangleLinkList_[triangleId][localLinkId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleLinkNumber)(
      const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleLinkList_.size()))
        return -1;
#endif
      return triangleLinkList_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() override {
      return &triangleLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
      const SimplexId &triangleId,
      const int &localStarId,
      SimplexId &starId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleStarList_.size()))
        return -1;
      if((localStarId < 0)
         || (localStarId >= (SimplexId)triangleStarList_[triangleId].size()))
        return -2;
#endif
      starId = triangleStarList_[triangleId][localStarId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(
      const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleStarList_.size()))
        return -1;
#endif
      return triangleStarList_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleStars)() override {
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
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexEdgeList_.size()))
        return -1;
      if((localEdgeId < 0)
         || (localEdgeId >= (SimplexId)vertexEdgeList_[vertexId].size()))
        return -2;
#endif
      edgeId = vertexEdgeList_[vertexId][localEdgeId];
      return 0;
    }

    inline SimplexId
      getVertexEdgeNumberInternal(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexEdgeList_.size()))
        return -1;
#endif
      return vertexEdgeList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexEdgesInternal() override {
      return &vertexEdgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexLink)(
      const SimplexId &vertexId,
      const int &localLinkId,
      SimplexId &linkId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexLinkList_.size()))
        return -1;
      if((localLinkId < 0)
         || (localLinkId >= (SimplexId)vertexLinkList_[vertexId].size()))
        return -2;
#endif
      linkId = vertexLinkList_[vertexId][localLinkId];

      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexLinkNumber)(
      const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexLinkList_.size()))
        return -1;
#endif
      return vertexLinkList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexLinks)() override {
      return &vertexLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexNeighbor)(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexNeighborList_.size()))
        return -1;
      if((localNeighborId < 0)
         || (localNeighborId
             >= (SimplexId)vertexNeighborList_[vertexId].size()))
        return -2;
#endif
      neighborId = vertexNeighborList_[vertexId][localNeighborId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexNeighborNumber)(
      const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif
      return vertexNeighborList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() override {
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
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexStarList_.size()))
        return -1;
      if((localStarId < 0)
         || (localStarId >= (SimplexId)vertexStarList_[vertexId].size()))
        return -2;
#endif
      starId = vertexStarList_[vertexId][localStarId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexStarNumber)(
      const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexStarList_.size()))
        return -1;
#endif
      return vertexStarList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexStars)() override {
      return &vertexStarList_;
    }

    inline int getVertexTriangleInternal(const SimplexId &vertexId,
                                         const int &localTriangleId,
                                         SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexTriangleList_.size()))
        return -1;
      if((localTriangleId < 0)
         || (localTriangleId
             >= (SimplexId)vertexTriangleList_[vertexId].size()))
        return -2;
#endif
      triangleId = vertexTriangleList_[vertexId][localTriangleId];
      return 0;
    }

    inline SimplexId getVertexTriangleNumberInternal(
      const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexTriangleList_.size()))
        return -1;
#endif
      return vertexTriangleList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() override {

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
        for(SimplexId i = 0; i < (SimplexId)edgeStarList_.size(); i++) {
          if(edgeStarList_[i].size() == 1) {
            boundaryEdges_[i] = true;
          }
        }
      } else if(getDimensionality() == 3) {
        preconditionTriangleStarsInternal();
        preconditionTriangleEdgesInternal();

        for(SimplexId i = 0; i < (SimplexId)triangleStarList_.size(); i++) {
          if(triangleStarList_[i].size() == 1) {
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

        for(SimplexId i = 0; i < (SimplexId)triangleStarList_.size(); i++) {
          if(triangleStarList_[i].size() == 1) {
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
        for(SimplexId i = 0; i < (SimplexId)vertexStarList_.size(); i++) {
          if(vertexStarList_[i].size() == 1) {
            boundaryVertices_[i] = true;
          }
        }
      } else if(getDimensionality() == 2) {
        preconditionEdgesInternal();
        preconditionEdgeStarsInternal();

        for(SimplexId i = 0; i < (SimplexId)edgeStarList_.size(); i++) {
          if(edgeStarList_[i].size() == 1) {
            boundaryVertices_[edgeList_[i].first] = true;
            boundaryVertices_[edgeList_[i].second] = true;
          }
        }
      } else if(getDimensionality() == 3) {
        preconditionTrianglesInternal();
        preconditionTriangleStarsInternal();

        for(SimplexId i = 0; i < (SimplexId)triangleStarList_.size(); i++) {
          if(triangleStarList_[i].size() == 1) {
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

      if(!cellEdgeList_.size()) {

        ThreeSkeleton threeSkeleton;
        threeSkeleton.setWrapper(this);

        threeSkeleton.buildCellEdges(vertexNumber_, *cellArray_, cellEdgeList_,
                                     &edgeList_, &vertexEdgeList_);
      }

      return 0;
    }

    inline int preconditionCellNeighborsInternal() override {

      if(!cellNeighborList_.size()) {
        ThreeSkeleton threeSkeleton;
        threeSkeleton.setWrapper(this);

        // choice here (for the more likely)
        threeSkeleton.buildCellNeighborsFromVertices(
          vertexNumber_, *cellArray_, cellNeighborList_, &vertexStarList_);
      }

      return 0;
    }

    inline int preconditionCellTrianglesInternal() override {

      if(!cellTriangleList_.size()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        if(triangleList_.size()) {
          // we already computed this guy, let's just get the cell triangles
          if(triangleStarList_.size()) {
            return twoSkeleton.buildTriangleList(
              vertexNumber_, *cellArray_, nullptr, nullptr, &cellTriangleList_);
          } else {
            // let's compute the triangle star while we're at it...
            // it's just a tiny overhead.
            return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                                 nullptr, &triangleStarList_,
                                                 &cellTriangleList_);
          }
        } else {
          // we have not computed this guy, let's do it while we're at it
          if(triangleStarList_.size()) {
            return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                                 &triangleList_, nullptr,
                                                 &cellTriangleList_);
          } else {
            // let's compute the triangle star while we're at it...
            // it's just a tiny overhead.
            return twoSkeleton.buildTriangleList(
              vertexNumber_, *cellArray_, &triangleList_, &triangleStarList_,
              &cellTriangleList_);
          }
        }
      }

      return 0;
    }

    inline int preconditionEdgesInternal() override {

      if(!edgeList_.size()) {
        OneSkeleton oneSkeleton;
        oneSkeleton.setWrapper(this);
        return oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_, edgeList_);
      }

      return 0;
    }

    inline int preconditionEdgeLinksInternal() override {

      if(!edgeLinkList_.size()) {

        if(getDimensionality() == 2) {
          preconditionEdgesInternal();
          preconditionEdgeStarsInternal();

          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          return oneSkeleton.buildEdgeLinks(
            edgeList_, edgeStarList_, *cellArray_, edgeLinkList_);
        } else if(getDimensionality() == 3) {
          preconditionEdgesInternal();
          preconditionEdgeStarsInternal();
          preconditionCellEdgesInternal();

          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          return oneSkeleton.buildEdgeLinks(
            edgeList_, edgeStarList_, cellEdgeList_, edgeLinkList_);
        } else {
          // unsupported dimension
          printErr("Unsupported dimension for edge link precondition");
          return -1;
        }
      }

      return 0;
    }

    inline int preconditionEdgeStarsInternal() override {

      if(!edgeStarList_.size()) {
        OneSkeleton oneSkeleton;
        oneSkeleton.setWrapper(this);
        return oneSkeleton.buildEdgeStars(vertexNumber_, *cellArray_,
                                          edgeStarList_, &edgeList_,
                                          &vertexStarList_);
      }
      return 0;
    }

    inline int preconditionEdgeTrianglesInternal() override {

      if(!edgeTriangleList_.size()) {

        // WARNING
        // here vertexStarList and triangleStarList will be computed (for
        // free) although they are not requireed to get the edgeTriangleList.
        // if memory usage is an issue, please change these pointers by nullptr.

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildEdgeTriangles(
          vertexNumber_, *cellArray_, edgeTriangleList_, &vertexStarList_,
          &edgeList_, &edgeStarList_, &triangleList_, &triangleStarList_,
          &cellTriangleList_);
      }

      return 0;
    }

    inline int preconditionTrianglesInternal() override {

      if(!triangleList_.size()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                      &triangleList_, &triangleStarList_,
                                      &cellTriangleList_);
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
          vertexNumber_, *cellArray_, triangleEdgeList_, &vertexEdgeList_,
          &edgeList_, &triangleList_, &triangleStarList_, &cellTriangleList_);
      }

      return 0;
    }

    inline int preconditionTriangleLinksInternal() override {

      if(!triangleLinkList_.size()) {

        preconditionTriangleStarsInternal();

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildTriangleLinks(
          triangleList_, triangleStarList_, *cellArray_, triangleLinkList_);
      }

      return 0;
    }

    inline int preconditionTriangleStarsInternal() override {

      if(!triangleStarList_.size()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildTriangleList(
          vertexNumber_, *cellArray_, &triangleList_, &triangleStarList_);
      }

      return 0;
    }

    inline int preconditionVertexEdgesInternal() override {

      if((SimplexId)vertexEdgeList_.size() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;

        if(!edgeList_.size()) {
          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_, edgeList_);
        }

        zeroSkeleton.setWrapper(this);
        return zeroSkeleton.buildVertexEdges(
          vertexNumber_, edgeList_, vertexEdgeList_);
      }
      return 0;
    }

    inline int preconditionVertexLinksInternal() override {

      if((SimplexId)vertexLinkList_.size() != vertexNumber_) {

        if(getDimensionality() == 2) {
          preconditionVertexStarsInternal();
          preconditionCellEdgesInternal();

          ZeroSkeleton zeroSkeleton;
          zeroSkeleton.setWrapper(this);
          return zeroSkeleton.buildVertexLinks(
            vertexStarList_, cellEdgeList_, edgeList_, vertexLinkList_);
        } else if(getDimensionality() == 3) {
          preconditionVertexStarsInternal();
          preconditionCellTrianglesInternal();

          ZeroSkeleton zeroSkeleton;
          zeroSkeleton.setWrapper(this);
          return zeroSkeleton.buildVertexLinks(
            vertexStarList_, cellTriangleList_, triangleList_, vertexLinkList_);
        } else {
          // unsupported dimension
          printErr("Unsupported dimension for vertex link precondition");
          return -1;
        }
      }
      return 0;
    }

    inline int preconditionVertexNeighborsInternal() override {

      if((SimplexId)vertexNeighborList_.size() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;
        zeroSkeleton.setWrapper(this);
        return zeroSkeleton.buildVertexNeighbors(
          vertexNumber_, *cellArray_, vertexNeighborList_, &edgeList_);
      }
      return 0;
    }

    inline int preconditionVertexStarsInternal() override {

      if((SimplexId)vertexStarList_.size() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;
        zeroSkeleton.setWrapper(this);

        return zeroSkeleton.buildVertexStars(
          vertexNumber_, *cellArray_, vertexStarList_);
      }
      return 0;
    }

    inline int preconditionVertexTrianglesInternal() override {

      if((SimplexId)vertexTriangleList_.size() != vertexNumber_) {

        preconditionTrianglesInternal();

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        twoSkeleton.buildVertexTriangles(
          vertexNumber_, triangleList_, vertexTriangleList_);
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
#pragma omp atomic write
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
  };
} // namespace ttk

#endif // _EXPLICITTRIANGULATION_H
