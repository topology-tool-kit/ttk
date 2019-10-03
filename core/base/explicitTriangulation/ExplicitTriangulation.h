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

    inline int getCellEdge(const SimplexId &cellId,
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

    inline SimplexId getCellEdgeNumber(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellEdgeList_.size()))
        return -1;
#endif
      return cellEdgeList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *getCellEdges() override {

      return &cellEdgeList_;
    }

    inline int getCellNeighbor(const SimplexId &cellId,
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

    inline SimplexId
      getCellNeighborNumber(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellNeighborList_.size()))
        return -1;
#endif
      return cellNeighborList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellNeighbors() override {
      return &cellNeighborList_;
    }

    inline int getCellTriangle(const SimplexId &cellId,
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
      getCellTriangleNumber(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)cellTriangleList_.size()))
        return -1;
#endif

      return cellTriangleList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellTriangles() override {

      return &cellTriangleList_;
    }

    inline int getCellVertex(const SimplexId &cellId,
                             const int &localVertexId,
                             SimplexId &vertexId) const override {
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
      getCellVertexNumber(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= cellNumber_))
        return -1;
      if((!cellArray_) || (!cellNumber_))
        return -2;
#endif
      return cellArray_[0];
    }

    int getDimensionality() const override {

      if((cellArray_) && (cellNumber_)) {
        return cellArray_[0] - 1;
      }

      return -1;
    }

    inline const std::vector<std::pair<SimplexId, SimplexId>> *
      getEdges() override {
      return &edgeList_;
    }

    inline int getEdgeLink(const SimplexId &edgeId,
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

    inline SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeLinkList_.size()))
        return -1;
#endif
      return edgeLinkList_[edgeId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *getEdgeLinks() override {

      return &edgeLinkList_;
    }

    inline int getEdgeStar(const SimplexId &edgeId,
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

    inline SimplexId getEdgeStarNumber(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeStarList_.size()))
        return -1;
#endif
      return edgeStarList_[edgeId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *getEdgeStars() override {
      return &edgeStarList_;
    }

    inline int getEdgeTriangle(const SimplexId &edgeId,
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
      getEdgeTriangleNumber(const SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeTriangleList_.size()))
        return -1;
#endif

      return edgeTriangleList_[edgeId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTriangles() override {

      return &edgeTriangleList_;
    }

    inline int getEdgeVertex(const SimplexId &edgeId,
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

    inline SimplexId getNumberOfCells() const override {
      return cellNumber_;
    }

    inline SimplexId getNumberOfEdges() const override {
      return edgeList_.size();
    }

    inline SimplexId getNumberOfTriangles() const override {
      return triangleList_.size();
    }

    inline SimplexId getNumberOfVertices() const override {
      return vertexNumber_;
    }

    inline const std::vector<std::vector<SimplexId>> *getTriangles() override {
      return &triangleList_;
    }

    inline int getTriangleEdge(const SimplexId &triangleId,
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

    inline SimplexId
      getTriangleEdgeNumber(const SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleEdgeList_.size()))
        return -1;
#endif

      return triangleEdgeList_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getTriangleEdges() override {

      return &triangleEdgeList_;
    }

    inline int getTriangleLink(const SimplexId &triangleId,
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

    inline SimplexId
      getTriangleLinkNumber(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleLinkList_.size()))
        return -1;
#endif
      return triangleLinkList_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getTriangleLinks() override {
      return &triangleLinkList_;
    }

    inline int getTriangleStar(const SimplexId &triangleId,
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

    inline SimplexId
      getTriangleStarNumber(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleStarList_.size()))
        return -1;
#endif
      return triangleStarList_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getTriangleStars() override {
      return &triangleStarList_;
    }

    inline int getTriangleVertex(const SimplexId &triangleId,
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

    inline int getVertexEdge(const SimplexId &vertexId,
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
      getVertexEdgeNumber(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexEdgeList_.size()))
        return -1;
#endif
      return vertexEdgeList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexEdges() override {
      return &vertexEdgeList_;
    }

    inline int getVertexLink(const SimplexId &vertexId,
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

    inline SimplexId
      getVertexLinkNumber(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexLinkList_.size()))
        return -1;
#endif
      return vertexLinkList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexLinks() override {
      return &vertexLinkList_;
    }

    inline int getVertexNeighbor(const SimplexId &vertexId,
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

    inline SimplexId
      getVertexNeighborNumber(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif
      return vertexNeighborList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexNeighbors() override {
      return &vertexNeighborList_;
    }

    inline int getVertexPoint(const SimplexId &vertexId,
                              float &x,
                              float &y,
                              float &z) const override {

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

    inline int getVertexStar(const SimplexId &vertexId,
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

    inline SimplexId
      getVertexStarNumber(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexStarList_.size()))
        return -1;
#endif
      return vertexStarList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexStars() override {
      return &vertexStarList_;
    }

    inline int getVertexTriangle(const SimplexId &vertexId,
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

    inline SimplexId
      getVertexTriangleNumber(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)vertexTriangleList_.size()))
        return -1;
#endif
      return vertexTriangleList_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexTriangles() override {

      return &vertexTriangleList_;
    }

    inline bool hasPreprocessedBoundaryEdges() const override {
      if(getDimensionality() == 1)
        return true;
      return (boundaryEdges_.size() != 0);
    }

    inline bool hasPreprocessedBoundaryTriangles() const override {
      if(getDimensionality() == 2)
        return true;
      return (boundaryTriangles_.size() != 0);
    }

    inline bool hasPreprocessedBoundaryVertices() const override {
      return (boundaryVertices_.size() != 0);
    }

    inline bool hasPreprocessedCellEdges() const override {
      return (cellEdgeList_.size() != 0);
    }

    inline bool hasPreprocessedCellNeighbors() const override {
      return (cellNeighborList_.size() != 0);
    }

    inline bool hasPreprocessedCellTriangles() const override {
      return (cellTriangleList_.size() != 0);
    }

    inline bool hasPreprocessedEdges() const override {
      return (edgeList_.size() != 0);
    }

    inline bool hasPreprocessedEdgeLinks() const override {
      return (edgeLinkList_.size() != 0);
    }

    inline bool hasPreprocessedEdgeStars() const override {
      return (edgeStarList_.size() != 0);
    }

    inline bool hasPreprocessedEdgeTriangles() const override {
      return (edgeTriangleList_.size() != 0);
    }

    inline bool hasPreprocessedTriangles() const override {
      return (triangleList_.size() != 0);
    }

    inline bool hasPreprocessedTriangleEdges() const override {
      return (triangleEdgeList_.size() != 0);
    }

    inline bool hasPreprocessedTriangleLinks() const override {
      return (triangleLinkList_.size() != 0);
    }

    inline bool hasPreprocessedTriangleStars() const override {
      return (triangleStarList_.size() != 0);
    }

    inline bool hasPreprocessedVertexEdges() const override {
      return (vertexEdgeList_.size() != 0);
    }

    inline bool hasPreprocessedVertexLinks() const override {
      return (vertexLinkList_.size() != 0);
    }

    inline bool hasPreprocessedVertexNeighbors() const override {
      return (vertexNeighborList_.size() != 0);
    }

    inline bool hasPreprocessedVertexStars() const override {
      return (vertexStarList_.size() != 0);
    }

    inline bool hasPreprocessedVertexTriangles() const override {
      return (vertexTriangleList_.size() != 0);
    }

    inline bool isEdgeOnBoundary(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)boundaryEdges_.size()))
        return false;
#endif
      return boundaryEdges_[edgeId];
    }

    inline bool isEmpty() const override {
      return !vertexNumber_;
    }

    inline bool
      isTriangleOnBoundary(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)boundaryTriangles_.size()))
        return false;
#endif
      return boundaryTriangles_[triangleId];
    }

    inline bool isVertexOnBoundary(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)boundaryVertices_.size()))
        return false;
#endif
      return boundaryVertices_[vertexId];
    }

    inline int preprocessBoundaryEdges() override {

      if((!boundaryEdges_.empty())
         && (boundaryEdges_.size() == edgeList_.size())) {
        return 0;
      }

      preprocessEdges();
      boundaryEdges_.resize(edgeList_.size(), false);

      if(getDimensionality() == 2) {
        preprocessEdgeStars();
        for(SimplexId i = 0; i < (SimplexId)edgeStarList_.size(); i++) {
          if(edgeStarList_[i].size() == 1) {
            boundaryEdges_[i] = true;
          }
        }
      } else if(getDimensionality() == 3) {
        preprocessTriangleStars();
        preprocessTriangleEdges();

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

    inline int preprocessBoundaryTriangles() override {

      if(getDimensionality() == 2)
        return 0;

      if((!boundaryTriangles_.empty())
         && (boundaryTriangles_.size() == triangleList_.size())) {
        return 0;
      }

      preprocessTriangles();
      boundaryTriangles_.resize(triangleList_.size(), false);

      if(getDimensionality() == 3) {
        preprocessTriangleStars();

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

    inline int preprocessBoundaryVertices() override {

      if((!boundaryVertices_.empty())
         && ((SimplexId)boundaryVertices_.size() == vertexNumber_))
        return 0;

      boundaryVertices_.resize(vertexNumber_, false);

      // create the list of boundary elements
      // create their star
      // look for singletons
      if(getDimensionality() == 1) {
        preprocessVertexStars();
        for(SimplexId i = 0; i < (SimplexId)vertexStarList_.size(); i++) {
          if(vertexStarList_[i].size() == 1) {
            boundaryVertices_[i] = true;
          }
        }
      } else if(getDimensionality() == 2) {
        preprocessEdges();
        preprocessEdgeStars();

        for(SimplexId i = 0; i < (SimplexId)edgeStarList_.size(); i++) {
          if(edgeStarList_[i].size() == 1) {
            boundaryVertices_[edgeList_[i].first] = true;
            boundaryVertices_[edgeList_[i].second] = true;
          }
        }
      } else if(getDimensionality() == 3) {
        preprocessTriangles();
        preprocessTriangleStars();

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

    inline int preprocessCellEdges() override {

      if(!cellEdgeList_.size()) {

        ThreeSkeleton threeSkeleton;
        threeSkeleton.setWrapper(this);

        threeSkeleton.buildCellEdges(vertexNumber_, cellNumber_, cellArray_,
                                     cellEdgeList_, &edgeList_,
                                     &vertexEdgeList_);
      }

      return 0;
    }

    inline int preprocessCellNeighbors() override {

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

    inline int preprocessCellTriangles() override {

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

    inline int preprocessEdges() override {

      if(!edgeList_.size()) {
        OneSkeleton oneSkeleton;
        oneSkeleton.setWrapper(this);
        return oneSkeleton.buildEdgeList(
          vertexNumber_, cellNumber_, cellArray_, edgeList_);
      }

      return 0;
    }

    inline int preprocessEdgeLinks() override {

      if(!edgeLinkList_.size()) {

        if(getDimensionality() == 2) {
          preprocessEdges();
          preprocessEdgeStars();

          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          return oneSkeleton.buildEdgeLinks(
            edgeList_, edgeStarList_, cellArray_, edgeLinkList_);
        } else if(getDimensionality() == 3) {
          preprocessEdges();
          preprocessEdgeStars();
          preprocessCellEdges();

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

    inline int preprocessEdgeStars() override {

      if(!edgeStarList_.size()) {
        OneSkeleton oneSkeleton;
        oneSkeleton.setWrapper(this);
        return oneSkeleton.buildEdgeStars(vertexNumber_, cellNumber_,
                                          cellArray_, edgeStarList_, &edgeList_,
                                          &vertexStarList_);
      }
      return 0;
    }

    inline int preprocessEdgeTriangles() override {

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

    inline int preprocessTriangles() override {

      if(!triangleList_.size()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);

        twoSkeleton.buildTriangleList(vertexNumber_, cellNumber_, cellArray_,
                                      &triangleList_, &triangleStarList_,
                                      &cellTriangleList_);
      }

      return 0;
    }

    inline int preprocessTriangleEdges() override {

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

    inline int preprocessTriangleLinks() override {

      if(!triangleLinkList_.size()) {

        preprocessTriangleStars();

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildTriangleLinks(
          triangleList_, triangleStarList_, cellArray_, triangleLinkList_);
      }

      return 0;
    }

    inline int preprocessTriangleStars() override {

      if(!triangleStarList_.size()) {

        TwoSkeleton twoSkeleton;
        twoSkeleton.setWrapper(this);
        return twoSkeleton.buildTriangleList(vertexNumber_, cellNumber_,
                                             cellArray_, &triangleList_,
                                             &triangleStarList_);
      }

      return 0;
    }

    inline int preprocessVertexEdges() override {

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

    inline int preprocessVertexLinks() override {

      if((SimplexId)vertexLinkList_.size() != vertexNumber_) {

        if(getDimensionality() == 2) {
          preprocessVertexStars();
          preprocessCellEdges();

          ZeroSkeleton zeroSkeleton;
          zeroSkeleton.setWrapper(this);
          return zeroSkeleton.buildVertexLinks(
            vertexStarList_, cellEdgeList_, edgeList_, vertexLinkList_);
        } else if(getDimensionality() == 3) {
          preprocessVertexStars();
          preprocessCellTriangles();

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

    inline int preprocessVertexNeighbors() override {

      if((SimplexId)vertexNeighborList_.size() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;
        zeroSkeleton.setWrapper(this);
        return zeroSkeleton.buildVertexNeighbors(
          vertexNumber_, cellNumber_, cellArray_, vertexNeighborList_,
          &edgeList_);
      }
      return 0;
    }

    inline int preprocessVertexStars() override {

      if((SimplexId)vertexStarList_.size() != vertexNumber_) {
        ZeroSkeleton zeroSkeleton;
        zeroSkeleton.setWrapper(this);

        return zeroSkeleton.buildVertexStars(
          vertexNumber_, cellNumber_, cellArray_, vertexStarList_);
      }
      return 0;
    }

    inline int preprocessVertexTriangles() override {

      if((SimplexId)vertexTriangleList_.size() != vertexNumber_) {

        preprocessTriangles();

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

  protected:
    int clear();

    bool doublePrecision_;
    SimplexId cellNumber_, vertexNumber_;
    const void *pointSet_;
    const LongSimplexId *cellArray_;
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include                  <ExplicitTriangulation.cpp>

#endif // _EXPLICITTRIANGULATION_H
