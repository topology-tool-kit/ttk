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

namespace ttk {

  class ExplicitTriangulation final : public AbstractTriangulation {

  public:
    ExplicitTriangulation();

    ~ExplicitTriangulation();

    int clear();

    inline int getCellEdgeInternal(const SimplexId &cellId,
                                   const int &localEdgeId,
                                   SimplexId &edgeId) const {

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

    inline SimplexId getCellEdgeNumberInternal(const SimplexId &cellId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellEdgeList_.size()))
        return -1;
#endif
      return cellEdgeList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *getCellEdgesInternal() {

      return &cellEdgeList_;
    }

    inline int getCellNeighborInternal(const SimplexId &cellId,
                                       const int &localNeighborId,
                                       SimplexId &neighborId) const {
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

    inline SimplexId
      getCellNeighborNumberInternal(const SimplexId &cellId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellNeighborList_.size()))
        return -1;
#endif
      return cellNeighborList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellNeighborsInternal() {
      return &cellNeighborList_;
    }

    inline int getCellTriangleInternal(const SimplexId &cellId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const {

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
      getCellTriangleNumberInternal(const SimplexId &cellId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellTriangleList_.size()))
        return -1;
#endif

      return cellTriangleList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellTrianglesInternal() {

      return &cellTriangleList_;
    }

    inline int getCellVertexInternal(const SimplexId &cellId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if((localVertexId < 0) || (localVertexId >= cellArray_[0]))
        return -2;
#endif
      vertexId = cellArray_[(cellArray_[0] + 1) * cellId + localVertexId + 1];
      return 0;
    }

    inline SimplexId
      getCellVertexNumberInternal(const SimplexId &cellId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if((!cellArray_) || (!cellNumber_))
        return -2;
#endif
      return cellArray_[0];
    }

    int getDimensionalityInternal() const {

      if((cellArray_) && (cellNumber_)) {
        return cellArray_[0] - 1;
      }

      return -1;
    }

    inline const std::vector<std::pair<SimplexId, SimplexId>> *
      getEdgesInternal() {
      return &edgeList_;
    }

    inline int getEdgeLinkInternal(const SimplexId &edgeId,
                                   const int &localLinkId,
                                   SimplexId &linkId) const {
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

    inline SimplexId getEdgeLinkNumberInternal(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeLinkList_.size()))
        return -1;
#endif
      return edgeLinkList_[edgeId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *getEdgeLinksInternal() {

      return &edgeLinkList_;
    }

    inline int getEdgeStarInternal(const SimplexId &edgeId,
                                   const int &localStarId,
                                   SimplexId &starId) const {
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

    inline SimplexId getEdgeStarNumberInternal(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeStarList_.size()))
        return -1;
#endif
      return edgeStarList_[edgeId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *getEdgeStarsInternal() {
      return &edgeStarList_;
    }

    inline int getEdgeTriangleInternal(const SimplexId &edgeId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const {

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
      getEdgeTriangleNumberInternal(const SimplexId &edgeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeTriangleList_.size()))
        return -1;
#endif

      return edgeTriangleList_[edgeId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTrianglesInternal() {

      return &edgeTriangleList_;
    }

    inline int getEdgeVertexInternal(const SimplexId &edgeId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const {
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

    inline SimplexId getNumberOfCellsInternal() const {
      return cellNumber_;
    }

    inline SimplexId getNumberOfEdgesInternal() const {
      return edgeList_.size();
    }

    inline SimplexId getNumberOfTrianglesInternal() const {
      return triangleList_.size();
    }

    inline SimplexId getNumberOfVerticesInternal() const {
      return vertexNumber_;
    }

    inline const std::vector<std::vector<SimplexId>> *getTrianglesInternal() {
      return &triangleList_;
    }

    inline int getTriangleEdgeInternal(const SimplexId &triangleId,
                                       const int &localEdgeId,
                                       SimplexId &edgeId) const {

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

    inline SimplexId
      getTriangleEdgeNumberInternal(const SimplexId &triangleId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleEdgeList_.size()))
        return -1;
#endif

      return triangleEdgeList_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getTriangleEdgesInternal() {

      return &triangleEdgeList_;
    }

    inline int getTriangleLinkInternal(const SimplexId &triangleId,
                                       const int &localLinkId,
                                       SimplexId &linkId) const {
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

    inline SimplexId
      getTriangleLinkNumberInternal(const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleLinkList_.size()))
        return -1;
#endif
      return triangleLinkList_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getTriangleLinksInternal() {
      return &triangleLinkList_;
    }

    inline int getTriangleStarInternal(const SimplexId &triangleId,
                                       const int &localStarId,
                                       SimplexId &starId) const {
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

    inline SimplexId
      getTriangleStarNumberInternal(const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleStarList_.size()))
        return -1;
#endif
      return triangleStarList_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getTriangleStarsInternal() {
      return &triangleStarList_;
    }

    inline int getTriangleVertexInternal(const SimplexId &triangleId,
                                         const int &localVertexId,
                                         SimplexId &vertexId) const {
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
                                     SimplexId &edgeId) const {
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
      getVertexEdgeNumberInternal(const SimplexId &vertexId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexEdgeList_.size()))
        return -1;
#endif
      return vertexEdgeList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *getVertexEdgesInternal() {
      return &vertexEdgeList_;
    }

    inline int getVertexLinkInternal(const SimplexId &vertexId,
                                     const int &localLinkId,
                                     SimplexId &linkId) const {

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

    inline SimplexId
      getVertexLinkNumberInternal(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexLinkList_.size()))
        return -1;
#endif
      return vertexLinkList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *getVertexLinksInternal() {
      return &vertexLinkList_;
    }

#ifdef TTK_ENABLE_KAMIKAZE
    inline int getVertexNeighbor(const SimplexId &vertexId,
                                 const int &localNeighborId,
                                 SimplexId &neighborId) const {
#else
    inline int getVertexNeighborInternal(const SimplexId &vertexId,
                                         const int &localNeighborId,
                                         SimplexId &neighborId) const {

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

#ifdef TTK_ENABLE_KAMIKAZE
    inline SimplexId getVertexNeighborNumber(const SimplexId &vertexId) const {
#else
    inline SimplexId
      getVertexNeighborNumberInternal(const SimplexId &vertexId) const {
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif
      return vertexNeighborList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexNeighborsInternal() {
      return &vertexNeighborList_;
    }

    inline int getVertexPointInternal(const SimplexId &vertexId,
                                      float &x,
                                      float &y,
                                      float &z) const {

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

    inline int getVertexStarInternal(const SimplexId &vertexId,
                                     const int &localStarId,
                                     SimplexId &starId) const {
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

    inline SimplexId
      getVertexStarNumberInternal(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexStarList_.size()))
        return -1;
#endif
      return vertexStarList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *getVertexStarsInternal() {
      return &vertexStarList_;
    }

    inline int getVertexTriangleInternal(const SimplexId &vertexId,
                                         const int &localTriangleId,
                                         SimplexId &triangleId) const {

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

    inline SimplexId
      getVertexTriangleNumberInternal(const SimplexId &vertexId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexTriangleList_.size()))
        return -1;
#endif
      return vertexTriangleList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() {

      return &vertexTriangleList_;
    }

    inline bool isEdgeOnBoundaryInternal(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)boundaryEdges_.size()))
        return false;
#endif
      return boundaryEdges_[edgeId];
    }

    inline bool isEmptyInternal() const {
      return !vertexNumber_;
    }

    inline bool
      isTriangleOnBoundaryInternal(const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)boundaryTriangles_.size()))
        return false;
#endif
      return boundaryTriangles_[triangleId];
    }

    inline bool isVertexOnBoundaryInternal(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)boundaryVertices_.size()))
        return false;
#endif
      return boundaryVertices_[vertexId];
    }

    inline int preconditionBoundaryEdgesInternal() {

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
        std::stringstream msg;
        msg << "[ExplicitTriangulation] Unsupported dimension for boundary "
            << "preprocessing." << std::endl;
        dMsg(std::cerr, msg.str(), infoMsg);
        return -1;
      }

      return 0;
    }

    inline int preconditionBoundaryTrianglesInternal() {

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
        std::stringstream msg;
        msg << "[ExplicitTriangulation] Unsupported dimension for boundary "
            << "preprocessing." << std::endl;
        dMsg(std::cerr, msg.str(), infoMsg);
        return -1;
      }

      return 0;
    }

    inline int preconditionBoundaryVerticesInternal() {

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
        std::stringstream msg;
        msg << "[ExplicitTriangulation] Unsupported dimension for boundary "
            << "preprocessing." << std::endl;
        dMsg(std::cerr, msg.str(), infoMsg);
        return -1;
      }

      return 0;
    }

    inline int preconditionCellEdgesInternal() {

      if(!cellEdgeList_.size()) {

        ThreeSkeleton threeSkeleton;
        threeSkeleton.setWrapper(this);

        threeSkeleton.buildCellEdges(vertexNumber_, cellNumber_, cellArray_,
                                     cellEdgeList_, &edgeList_,
                                     &vertexEdgeList_);
      }

      return 0;
    }

    inline int preconditionCellNeighborsInternal() {

      if(!cellNeighborList_.size()) {
        ThreeSkeleton threeSkeleton;
        threeSkeleton.setWrapper(this);

        // choice here (for the more likely)
        threeSkeleton.buildCellNeighborsFromVertices(
          vertexNumber_, cellNumber_, cellArray_, cellNeighborList_,
          &vertexStarList_);
      }

      return 0;
    }

    inline int preconditionCellTrianglesInternal() {

      if(!cellTriangleList_.size()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        if(triangleList_.size()) {
          // we already computed this guy, let's just get the cell triangles
          if(triangleStarList_.size()) {
            return twoSkeleton.buildTriangleList(vertexNumber_, cellNumber_,
                                                 cellArray_, NULL, NULL,
                                                 &cellTriangleList_);
          } else {
            // let's compute the triangle star while we're at it...
            // it's just a tiny overhead.
            return twoSkeleton.buildTriangleList(
              vertexNumber_, cellNumber_, cellArray_, NULL, &triangleStarList_,
              &cellTriangleList_);
          }
        } else {
          // we have not computed this guy, let's do it while we're at it
          if(triangleStarList_.size()) {
            return twoSkeleton.buildTriangleList(vertexNumber_, cellNumber_,
                                                 cellArray_, &triangleList_,
                                                 NULL, &cellTriangleList_);
          } else {
            // let's compute the triangle star while we're at it...
            // it's just a tiny overhead.
            return twoSkeleton.buildTriangleList(
              vertexNumber_, cellNumber_, cellArray_, &triangleList_,
              &triangleStarList_, &cellTriangleList_);
          }
        }
      }

      return 0;
    }

    inline int preconditionEdgesInternal() {

      if(!edgeList_.size()) {
        OneSkeleton oneSkeleton;
        oneSkeleton.setWrapper(this);
        return oneSkeleton.buildEdgeList(
          vertexNumber_, cellNumber_, cellArray_, edgeList_);
      }

      return 0;
    }

    inline int preconditionEdgeLinksInternal() {

      if(!edgeLinkList_.size()) {

        if(getDimensionality() == 2) {
          preconditionEdgesInternal();
          preconditionEdgeStarsInternal();

          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          return oneSkeleton.buildEdgeLinks(
            edgeList_, edgeStarList_, cellArray_, edgeLinkList_);
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
          std::stringstream msg;
          msg << "[ExplicitTriangulation] Unsupported dimension for edge link "
              << "preprocessing." << std::endl;
          dMsg(std::cerr, msg.str(), infoMsg);
          return -1;
        }
      }

      return 0;
    }

    inline int preconditionEdgeStarsInternal() {

      if(!edgeStarList_.size()) {
        OneSkeleton oneSkeleton;
        oneSkeleton.setWrapper(this);
        return oneSkeleton.buildEdgeStars(vertexNumber_, cellNumber_,
                                          cellArray_, edgeStarList_, &edgeList_,
                                          &vertexStarList_);
      }
      return 0;
    }

    inline int preconditionEdgeTrianglesInternal() {

      if(!edgeTriangleList_.size()) {

        // WARNING
        // here vertexStarList and triangleStarList will be computed (for
        // free) although they are not requireed to get the edgeTriangleList.
        // if memory usage is an issue, please change these pointers by NULL.

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildEdgeTriangles(
          vertexNumber_, cellNumber_, cellArray_, edgeTriangleList_,
          &vertexStarList_, &edgeList_, &edgeStarList_, &triangleList_,
          &triangleStarList_, &cellTriangleList_);
      }

      return 0;
    }

    inline int preconditionTrianglesInternal() {

      if(!triangleList_.size()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        twoSkeleton.buildTriangleList(vertexNumber_, cellNumber_, cellArray_,
                                      &triangleList_, &triangleStarList_,
                                      &cellTriangleList_);
      }

      return 0;
    }

    inline int preconditionTriangleEdgesInternal() {

      if(!triangleEdgeList_.size()) {

        // WARNING
        // here triangleStarList and cellTriangleList will be computed (for
        // free) although they are not requireed to get the edgeTriangleList.
        // if memory usage is an issue, please change these pointers by NULL.

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        return twoSkeleton.buildTriangleEdgeList(
          vertexNumber_, cellNumber_, cellArray_, triangleEdgeList_,
          &vertexEdgeList_, &edgeList_, &triangleList_, &triangleStarList_,
          &cellTriangleList_);
      }

      return 0;
    }

    inline int preconditionTriangleLinksInternal() {

      if(!triangleLinkList_.size()) {

        preconditionTriangleStarsInternal();

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildTriangleLinks(
          triangleList_, triangleStarList_, cellArray_, triangleLinkList_);
      }

      return 0;
    }

    inline int preconditionTriangleStarsInternal() {

      if(!triangleStarList_.size()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildTriangleList(vertexNumber_, cellNumber_,
                                             cellArray_, &triangleList_,
                                             &triangleStarList_);
      }

      return 0;
    }

    inline int preconditionVertexEdgesInternal() {

      if((SimplexId)vertexEdgeList_.size() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;

        if(!edgeList_.size()) {
          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          oneSkeleton.buildEdgeList(
            vertexNumber_, cellNumber_, cellArray_, edgeList_);
        }

        zeroSkeleton.setWrapper(this);
        return zeroSkeleton.buildVertexEdges(
          vertexNumber_, edgeList_, vertexEdgeList_);
      }
      return 0;
    }

    inline int preconditionVertexLinksInternal() {

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
          std::stringstream msg;
          msg << "[ExplicitTriangulation] Unsupported dimension for vertex"
              << " link preprocessing." << std::endl;
          dMsg(std::cerr, msg.str(), infoMsg);
          return -1;
        }
      }
      return 0;
    }

    inline int preconditionVertexNeighborsInternal() {

      if((SimplexId)vertexNeighborList_.size() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;
        zeroSkeleton.setWrapper(this);
        return zeroSkeleton.buildVertexNeighbors(
          vertexNumber_, cellNumber_, cellArray_, vertexNeighborList_,
          &edgeList_);
      }
      return 0;
    }

    inline int preconditionVertexStarsInternal() {

      if((SimplexId)vertexStarList_.size() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;
        zeroSkeleton.setWrapper(this);

        return zeroSkeleton.buildVertexStars(
          vertexNumber_, cellNumber_, cellArray_, vertexStarList_);
      }
      return 0;
    }

    inline int preconditionVertexTrianglesInternal() {

      if((SimplexId)vertexTriangleList_.size() != vertexNumber_) {

        preconditionTrianglesInternal();

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        twoSkeleton.buildVertexTriangles(
          vertexNumber_, triangleList_, vertexTriangleList_);
      }

      return 0;
    }

    inline int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *cellArray) {

      if(cellNumber_)
        clear();

      cellNumber_ = cellNumber;
      cellArray_ = cellArray;

      return 0;
    }

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
    const LongSimplexId *cellArray_;
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include                  <ExplicitTriangulation.cpp>

#endif // _EXPLICITTRIANGULATION_H
