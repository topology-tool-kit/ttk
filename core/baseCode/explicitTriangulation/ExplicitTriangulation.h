/// \ingroup baseCode
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
#include                  <AbstractTriangulation.h>
#include                  <OneSkeleton.h>
#include                  <ThreeSkeleton.h>
#include                  <TwoSkeleton.h>
#include                  <ZeroSkeleton.h>

namespace ttk{
  
  class ExplicitTriangulation : public AbstractTriangulation{

    public:
        
      ExplicitTriangulation();
      
      ~ExplicitTriangulation();
     
      inline int getCellEdge(const int &cellId, 
        const int &localEdgeId, int &edgeId) const{
       
#ifndef withKamikaze
        if((cellId < 0)||(cellId >= (int) cellEdgeList_.size()))
          return -1;
        if((localEdgeId < 0)
          ||(localEdgeId >= (int) cellEdgeList_[cellId].size()))
          return -2;
#endif
        edgeId = cellEdgeList_[cellId][localEdgeId];
        return 0;
      }
      
      inline int getCellEdgeNumber(const int &cellId) const {
#ifndef withKamikaze
        if((cellId < 0)||(cellId >= (int) cellEdgeList_.size()))
          return -1;
#endif
        return cellEdgeList_[cellId].size();
      }
      
      inline const vector<vector<int> > *getCellEdges(){
        
        return &cellEdgeList_;
      }
      
      inline int getCellNeighbor(const int &cellId,
        const int &localNeighborId, int &neighborId) const{
#ifndef withKamikaze
        if((cellId < 0)||(cellId >= (int) cellNeighborList_.size()))
          return -1;
        if((localNeighborId < 0)
          ||(localNeighborId >= (int) cellNeighborList_[cellId].size()))
          return -2;
#endif
        neighborId = cellNeighborList_[cellId][localNeighborId];
        return 0;
      }
      
      inline int getCellNeighborNumber(const int &cellId) const{
#ifndef withKamikaze
        if((cellId < 0)||(cellId >= (int) cellNeighborList_.size()))
          return -1;
#endif
        return cellNeighborList_[cellId].size();
      }
      
      inline const vector<vector<int> > *getCellNeighbors() { 
        return &cellNeighborList_;}
      
      inline int getCellTriangle(const int &cellId,
        const int &localTriangleId, int &triangleId) const{
        
#ifndef withKamikaze
        if((cellId < 0)||(cellId >= (int) cellTriangleList_.size()))
          return -1;
        if((localTriangleId < 0)
          ||(localTriangleId >= (int) cellTriangleList_[cellId].size()))
          return -2;
#endif
        triangleId = cellTriangleList_[cellId][localTriangleId];
        
        return 0;
      }
      
      inline int getCellTriangleNumber(const int &cellId) const{
       
#ifndef withKamikaze
        if((cellId < 0)||(cellId >= (int) cellTriangleList_.size()))
          return -1;
#endif
        
        return cellTriangleList_[cellId].size();
      }
      
      inline const vector<vector<int> > *getCellTriangles(){
        
        return &cellTriangleList_;
      }
      
      inline int getCellVertex(const int &cellId,
        const int &localVertexId, int &vertexId) const{
#ifndef withKamikaze
        if((cellId < 0)||(cellId >= cellNumber_))
          return -1;
        if((localVertexId < 0)
          ||(localVertexId >= cellArray_[0]))
          return -2;
#endif
        vertexId = cellArray_[(cellArray_[0] + 1)*cellId + localVertexId + 1];
        return 0;
      }
      
      inline int getCellVertexNumber(const int &cellId) const{
#ifndef withKamikaze
        if((cellId < 0)||(cellId >= cellNumber_))
          return -1;
        if((!cellArray_)||(!cellNumber_))
          return -2;
#endif
        return cellArray_[0];
      }
      
      int getDimensionality() const{
        
        if((cellArray_)&&(cellNumber_)){
          return cellArray_[0] - 1;
        }
        
        return -1;
      }
      
      inline const vector<pair<int, int> > *getEdges() {
        return &edgeList_;
      }
      
      inline int getEdgeLink(const int &edgeId, 
        const int &localLinkId, int &linkId) const{
#ifndef withKamikaze
        if((edgeId < 0)||(edgeId >= (int) edgeLinkList_.size()))
          return -1;
        if((localLinkId < 0)
          ||(localLinkId >= (int) edgeLinkList_[edgeId].size()))
          return -2;
#endif
        linkId = edgeLinkList_[edgeId][localLinkId];
        return 0;
      }
      
      inline int getEdgeLinkNumber(const int &edgeId) const{
#ifndef withKamikaze
        if((edgeId < 0)||(edgeId >= (int) edgeLinkList_.size()))
          return -1;
#endif
        return edgeLinkList_[edgeId].size();
      }
      
      inline const vector<vector<int> > *getEdgeLinks(){
        
        return &edgeLinkList_;
      }
      
      inline int getEdgeStar(const int &edgeId,
        const int &localStarId, int &starId) const{
#ifndef withKamikaze
        if((edgeId < 0)||(edgeId >= (int) edgeStarList_.size()))
          return -1;
        if((localStarId < 0)
          ||(localStarId >= (int) edgeStarList_[edgeId].size()))
          return -2;
#endif
        starId = edgeStarList_[edgeId][localStarId];
        return 0;
      }
      
      inline int getEdgeStarNumber(const int &edgeId) const{
#ifndef withKamikaze
        if((edgeId < 0)||(edgeId >= (int) edgeStarList_.size()))
          return -1;
#endif
        return edgeStarList_[edgeId].size();
      }
      
      inline const vector<vector<int> > *getEdgeStars(){
        return &edgeStarList_;
      }
      
      inline int getEdgeTriangle(const int &edgeId,
        const int &localTriangleId, int &triangleId) const{
          
#ifndef withKamikaze
        if((edgeId < 0)||(edgeId >= (int) edgeTriangleList_.size()))
          return -1;
        if((localTriangleId < 0)
          ||(localTriangleId >= (int) edgeTriangleList_[edgeId].size()))
          return -2;
#endif
        
        triangleId = edgeTriangleList_[edgeId][localTriangleId];
          
        return 0;
      }
      
      inline int getEdgeTriangleNumber(const int &edgeId) const{
        
#ifndef withKamikaze
        if((edgeId < 0)||(edgeId >= (int) edgeTriangleList_.size()))
          return -1;
#endif
        
        return edgeTriangleList_[edgeId].size();
      }
      
      inline const vector<vector<int> > *getEdgeTriangles(){
        
        return &edgeTriangleList_;
      }
      
      inline int getEdgeVertex(const int &edgeId,
        const int &localVertexId, int &vertexId) const{
#ifndef withKamikaze
        if((edgeId < 0)||(edgeId >= (int) edgeList_.size()))
          return -1;
        if((localVertexId != 0)&&(localVertexId != 1))
          return -2;
#endif
        if(!localVertexId)
          vertexId = edgeList_[edgeId].first;
        else
          vertexId = edgeList_[edgeId].second;
        return 0;
      }
      
      inline int getNumberOfCells() const { return cellNumber_;}
      
      inline int getNumberOfEdges() const{
        return edgeList_.size();
      }
      
      inline int getNumberOfTriangles() const{
        return triangleList_.size();
      }
      
      inline int getNumberOfVertices() const { return vertexNumber_;}
      
      inline const vector<vector<int> > *getTriangles(){
        return &triangleList_;
      }
      
      inline int getTriangleEdge(const int &triangleId, 
        const int &localEdgeId, int &edgeId) const{
          
#ifndef withKamikaze
        if((triangleId < 0)||(triangleId >= (int) triangleEdgeList_.size()))
          return -1;
        if((localEdgeId < 0)||(localEdgeId > 2))
          return -2;
#endif  
        
        edgeId = triangleEdgeList_[triangleId][localEdgeId];
          
        return 0;
      }
      
      inline int getTriangleEdgeNumber(const int &triangleId) const{
       
#ifndef withKamikaze
        if((triangleId < 0)||(triangleId >= (int) triangleEdgeList_.size()))
          return -1;
#endif 
        
        return triangleEdgeList_[triangleId].size();
      }
      
      inline const vector<vector<int> > *getTriangleEdges(){
        
        return &triangleEdgeList_;
      }
      
      inline int getTriangleLink(const int &triangleId, 
        const int &localLinkId, int &linkId) const{
#ifndef withKamikaze
        if((triangleId < 0)||(triangleId >= (int) triangleLinkList_.size()))
          return -1;
        if((localLinkId < 0)
          ||(localLinkId >= (int) triangleLinkList_[triangleId].size()))
          return -2;
#endif
        linkId = triangleLinkList_[triangleId][localLinkId];
        return 0;
      }
      
      inline int getTriangleLinkNumber(const int &triangleId) const{
#ifndef withKamikaze
        if((triangleId < 0)||(triangleId >= (int) triangleLinkList_.size()))
          return -1;
#endif
        return triangleLinkList_[triangleId].size();
      }
      
      inline const vector<vector<int> > *getTriangleLinks(){
        return &triangleLinkList_;
      }
      
      inline int getTriangleStar(const int &triangleId,
        const int &localStarId, int &starId) const{
#ifndef withKamikaze
        if((triangleId < 0)||(triangleId >= (int) triangleStarList_.size()))
          return -1;
        if((localStarId < 0)
          ||(localStarId >= (int) triangleStarList_[triangleId].size()))
          return -2;
#endif
        starId = triangleStarList_[triangleId][localStarId];
        return 0;
      }
      
      inline int getTriangleStarNumber(const int &triangleId) const{
#ifndef withKamikaze
        if((triangleId < 0)||(triangleId >= (int) triangleStarList_.size()))
          return -1;
#endif
        return triangleStarList_[triangleId].size();
      }
      
      inline const vector<vector<int> > *getTriangleStars(){
        return &triangleStarList_;
      }
      
      inline int getTriangleVertex(const int &triangleId,
        const int &localVertexId, int &vertexId) const{
#ifndef withKamikaze
        if((triangleId < 0)||(triangleId >= (int) triangleList_.size()))
          return -1;
        if((localVertexId < 0)
          ||(localVertexId >= (int) triangleList_[triangleId].size()))
          return -2;
#endif 
        vertexId = triangleList_[triangleId][localVertexId];
        return 0;
      }
      
      inline int getVertexEdge(const int &vertexId, 
        const int &localEdgeId, int &edgeId) const{
#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= (int) vertexEdgeList_.size()))
          return -1;
        if((localEdgeId < 0)
          ||(localEdgeId >= (int) vertexEdgeList_[vertexId].size()))
          return -2;
#endif
        edgeId = vertexEdgeList_[vertexId][localEdgeId];
        return 0;
      }
      
      inline int getVertexEdgeNumber(const int &vertexId) const{

#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= (int) vertexEdgeList_.size()))
          return -1;
#endif
        return vertexEdgeList_[vertexId].size();
      }
      
      inline const vector<vector<int> > *getVertexEdges(){
        return &vertexEdgeList_;
      }
      
      inline int getVertexLink(const int &vertexId, 
        const int &localLinkId, int &linkId) const{
         
        
#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= (int) vertexLinkList_.size()))
          return -1;
        if((localLinkId < 0)
          ||(localLinkId >= (int) vertexLinkList_[vertexId].size()))
          return -2;
#endif
        linkId = vertexLinkList_[vertexId][localLinkId];
        
        return 0;
      }

      inline int getVertexLinkNumber(const int &vertexId) const{
#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= (int) vertexLinkList_.size()))
          return -1;
#endif
        return vertexLinkList_[vertexId].size();
      }
      
      inline const vector<vector<int> > *getVertexLinks(){
        return &vertexLinkList_;
      }
      
      inline int getVertexNeighbor(const int &vertexId, 
        const int &localNeighborId, int &neighborId) const{
#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= (int) vertexNeighborList_.size()))
          return -1;
        if((localNeighborId < 0)
          ||(localNeighborId >= (int) vertexNeighborList_[vertexId].size()))
          return -2;
#endif
        neighborId = vertexNeighborList_[vertexId][localNeighborId];
        return 0;
      }
      
      inline int getVertexNeighborNumber(const int &vertexId) const{
#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= vertexNumber_))
          return -1;
#endif
        return vertexNeighborList_[vertexId].size();
      }
      
      inline const vector<vector<int> > *getVertexNeighbors(){
        return &vertexNeighborList_;
      }
      
      inline int getVertexPoint(const int &vertexId, 
        float &x, float &y, float &z) const{

#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= vertexNumber_))
          return -1;
#endif
        
        x = pointSet_[3*vertexId];
        y = pointSet_[3*vertexId + 1];
        z = pointSet_[3*vertexId + 2];
      
        return 0;
      }
      
      inline int getVertexStar(const int &vertexId, 
        const int &localStarId, int &starId) const{
#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= (int) vertexStarList_.size()))
          return -1;
        if((localStarId < 0)
          ||(localStarId >= (int) vertexStarList_[vertexId].size()))
          return -2;
#endif
        starId = vertexStarList_[vertexId][localStarId];
        return 0;
      }
      
      inline int getVertexStarNumber(const int &vertexId) const{
#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= (int) vertexStarList_.size()))
          return -1;
#endif
        return vertexStarList_[vertexId].size();
      }
      
      inline const vector<vector<int> > *getVertexStars(){
        return &vertexStarList_;
      }
      
      inline int getVertexTriangle(const int &vertexId,
        const int &localTriangleId, int &triangleId) const{
        
#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= (int) vertexTriangleList_.size()))
          return -1;
        if((localTriangleId < 0)
          ||(localTriangleId >= (int) vertexTriangleList_[vertexId].size()))
          return -2;
#endif
        triangleId = vertexTriangleList_[vertexId][localTriangleId];
        return 0;
      }
      
      inline int getVertexTriangleNumber(const int &vertexId) const{
        
#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= (int) vertexTriangleList_.size()))
          return -1;
#endif
        return vertexTriangleList_[vertexId].size();
      }
      
      inline const vector<vector<int> > *getVertexTriangles(){
        
        return &vertexTriangleList_;
      }
      
      inline bool hasPreprocessedBoundaryEdges() const{
        return (boundaryEdges_.size() != 0);
      }
      
      inline bool hasPreprocessedBoundaryTriangles() const{
        return (boundaryTriangles_.size() != 0);
      }
      
      inline bool hasPreprocessedBoundaryVertices() const{
        return (boundaryVertices_.size() != 0);
      }
      
      inline bool hasPreprocessedCellEdges() const{
        return (cellEdgeList_.size() != 0);
      }
      
      inline bool hasPreprocessedCellNeighbors() const{
        return (cellNeighborList_.size() != 0);
      }
      
      inline bool hasPreprocessedCellTriangles() const{
        return (cellTriangleList_.size() != 0);
      }
      
      inline bool hasPreprocessedEdges() const{
        return (edgeList_.size() != 0);
      }
      
      inline bool hasPreprocessedEdgeLinks() const{
        return (edgeLinkList_.size() != 0);
      }
      
      inline bool hasPreprocessedEdgeStars() const{
        return (edgeStarList_.size() != 0);
      }
      
      inline bool hasPreprocessedEdgeTriangles() const{
        return (edgeTriangleList_.size() != 0);
      }
      
      inline bool hasPreprocessedTriangles() const{
        return (triangleList_.size() != 0);
      }
      
      inline bool hasPreprocessedTriangleEdges() const{
        return (triangleEdgeList_.size() != 0);
      }
      
      inline bool hasPreprocessedTriangleLinks() const{
        return (triangleLinkList_.size() != 0);
      }
      
      inline bool hasPreprocessedTriangleStars() const{
        return (triangleStarList_.size() != 0);
      }
      
      inline bool hasPreprocessedVertexEdges() const{
        return (vertexEdgeList_.size() != 0);
      }
      
      inline bool hasPreprocessedVertexLinks() const{
        return (vertexLinkList_.size() != 0);
      }
      
      inline bool hasPreprocessedVertexNeighbors() const{
        return (vertexNeighborList_.size() != 0);
      }
      
      inline bool hasPreprocessedVertexStars() const{
        return (vertexStarList_.size() != 0);
      }
      
      inline bool hasPreprocessedVertexTriangles() const{
        return (vertexTriangleList_.size() != 0);
      }
      
      inline bool isEdgeOnBoundary(const int &edgeId) const{
#ifndef withKamikaze
        if((edgeId < 0)||(edgeId >= (int) boundaryEdges_.size()))
          return false;
#endif
        return boundaryEdges_[edgeId];
      }
      
      inline bool isEmpty() const { return !vertexNumber_;}
      
      inline bool isTriangleOnBoundary(const int &triangleId) const{
#ifndef withKamikaze
        if((triangleId < 0)||(triangleId >= (int) boundaryTriangles_.size()))
          return false;
#endif
        return boundaryTriangles_[triangleId];
      }
      
      inline bool isVertexOnBoundary(const int &vertexId) const{
#ifndef withKamikaze
        if((vertexId < 0)||(vertexId >= (int) boundaryVertices_.size()))
          return false;
#endif
        return boundaryVertices_[vertexId];
      }
      
      inline int preprocessBoundaryEdges(){
        
        if((!boundaryEdges_.empty())
          &&(boundaryEdges_.size() == edgeList_.size())){
          return 0; 
        }
        
        preprocessEdges();
        boundaryEdges_.resize(edgeList_.size(), false);
       
        if(getDimensionality() == 2){
          preprocessEdgeStars();
          for(int i = 0; i < (int) edgeStarList_.size(); i++){
            if(edgeStarList_[i].size() == 1){
              boundaryEdges_[i] = true;
            }
          }
        }
        else if(getDimensionality() == 3){
          preprocessTriangleStars();
          preprocessTriangleEdges();
          
          for(int i = 0; i < (int) triangleStarList_.size(); i++){
            if(triangleStarList_[i].size() == 1){
              for(int j = 0; j < 3; j++){
                boundaryEdges_[triangleEdgeList_[i][j]] = true;
              }
            }
          }
        }
        else{
          // unsupported dimension
          stringstream msg;
          msg << "[ExplicitTriangulation] Unsupported dimension for boundary "
            << "preprocessing."
            << endl;
          dMsg(cerr, msg.str(), infoMsg);
          return -1;
        }
        
        return 0;
      }
      
      inline int preprocessBoundaryTriangles(){
        
        if((!boundaryTriangles_.empty())
          &&(boundaryTriangles_.size() == triangleList_.size())){
          return 0; 
        }
        
        preprocessTriangles();
        boundaryTriangles_.resize(triangleList_.size(), false);
        
        if(getDimensionality() == 3){
          preprocessTriangleStars();
          
          for(int i = 0; i < (int) triangleStarList_.size(); i++){
            if(triangleStarList_[i].size() == 1){
              boundaryTriangles_[i] = true;
            }
          }
        }
        else{
          // unsupported dimension
          stringstream msg;
          msg << "[ExplicitTriangulation] Unsupported dimension for boundary "
            << "preprocessing."
            << endl;
          dMsg(cerr, msg.str(), infoMsg);
          return -1;
        }
       
        return 0;
      }
      
      inline int preprocessBoundaryVertices(){
        
        if((!boundaryVertices_.empty())
          &&((int) boundaryVertices_.size() == vertexNumber_))
          return 0;
        
        boundaryVertices_.resize(vertexNumber_, false);
        
        // create the list of boundary elements
        // create their star
        // look for singletons
        if(getDimensionality() == 2){
          preprocessEdges();
          preprocessEdgeStars();
          
          for(int i = 0; i < (int) edgeStarList_.size(); i++){
            if(edgeStarList_[i].size() == 1){
              boundaryVertices_[edgeList_[i].first] = true;
              boundaryVertices_[edgeList_[i].second] = true;
            }
          }
        }
        else if(getDimensionality() == 3){
          preprocessTriangles();
          preprocessTriangleStars();
          
          for(int i = 0; i < (int) triangleStarList_.size(); i++){
            if(triangleStarList_[i].size() == 1){
              boundaryVertices_[triangleList_[i][0]] = true;
              boundaryVertices_[triangleList_[i][1]] = true;
              boundaryVertices_[triangleList_[i][2]] = true;
            }
          }
        }
        else{
          // unsupported dimension
          stringstream msg;
          msg << "[ExplicitTriangulation] Unsupported dimension for boundary "
            << "preprocessing."
            << endl;
          dMsg(cerr, msg.str(), infoMsg);
          return -1;
        }
        
        return 0;
      }
      
      inline int preprocessCellEdges(){
        
        if(!cellEdgeList_.size()){
          
          ThreeSkeleton threeSkeleton;
          threeSkeleton.setWrapper(this);
          
          threeSkeleton.buildCellEdges(vertexNumber_,
            cellNumber_, cellArray_, cellEdgeList_,
            &edgeList_, &vertexEdgeList_);
        }
        
        return 0;
      }
      
      inline int preprocessCellNeighbors(){
        
        if(!cellNeighborList_.size()){
          ThreeSkeleton threeSkeleton;
          threeSkeleton.setWrapper(this);
          
          // choice here (for the more likely)
          threeSkeleton.buildCellNeighborsFromVertices(vertexNumber_,
            cellNumber_, cellArray_, cellNeighborList_, &vertexStarList_);
        }
        
        return 0;
      }
      
      inline int preprocessCellTriangles(){
        
        if(!cellTriangleList_.size()){
          
          TwoSkeleton twoSkeleton;
          twoSkeleton.setWrapper(this);
          
          if(triangleList_.size()){
            // we already computed this guy, let's just get the cell triangles
            if(triangleStarList_.size()){
              return twoSkeleton.buildTriangleList(
                vertexNumber_, cellNumber_, cellArray_,
                NULL, NULL, &cellTriangleList_);
            }
            else{
              // let's compute the triangle star while we're at it...
              // it's just a tiny overhead.
              return twoSkeleton.buildTriangleList(
                vertexNumber_, cellNumber_, cellArray_,
                NULL, &triangleStarList_, &cellTriangleList_);
            }
          }
          else{
            // we have not computed this guy, let's do it while we're at it
            if(triangleStarList_.size()){
              return twoSkeleton.buildTriangleList(
                vertexNumber_, cellNumber_, cellArray_,
                &triangleList_, NULL, &cellTriangleList_);
            }
            else{
              // let's compute the triangle star while we're at it...
              // it's just a tiny overhead.
              return twoSkeleton.buildTriangleList(
                vertexNumber_, cellNumber_, cellArray_,
                &triangleList_, &triangleStarList_, &cellTriangleList_);
            }
          }
        }
        
        return 0;
      }
      
      inline int preprocessEdges(){
        
        if(!edgeList_.size()){
          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          return oneSkeleton.buildEdgeList(vertexNumber_, cellNumber_,
            cellArray_, edgeList_);
        }
        
        return 0;
      }
      
      inline int preprocessEdgeLinks(){
        
        if(!edgeLinkList_.size()){
          
          if(getDimensionality() == 2){
            preprocessEdges();
            preprocessEdgeStars();
            
            
            OneSkeleton oneSkeleton;
            oneSkeleton.setWrapper(this);
            return oneSkeleton.buildEdgeLinks(
              edgeList_, edgeStarList_, 
              cellArray_, edgeLinkList_);
          }
          else if(getDimensionality() == 3){
            preprocessEdges();
            preprocessEdgeStars();
            preprocessCellEdges();
            
            OneSkeleton oneSkeleton;
            oneSkeleton.setWrapper(this);
            return oneSkeleton.buildEdgeLinks(
              edgeList_, edgeStarList_, cellEdgeList_, edgeLinkList_);
          }
          else{
            // unsupported dimension
            stringstream msg;
            msg 
              << "[ExplicitTriangulation] Unsupported dimension for edge link "
              << "preprocessing."
              << endl;
            dMsg(cerr, msg.str(), infoMsg);
            return -1;
          }
        }
        
        return 0;
      }
      
      inline int preprocessEdgeStars(){
        
        if(!edgeStarList_.size()){
          OneSkeleton oneSkeleton;
          oneSkeleton.setWrapper(this);
          return oneSkeleton.buildEdgeStars(vertexNumber_, cellNumber_,
            cellArray_, edgeStarList_, &edgeList_, &vertexStarList_);
        }
        return 0;
      }
      
      inline int preprocessEdgeTriangles(){
        
        if(!edgeTriangleList_.size()){
          
          // WARNING
          // here vertexStarList and triangleStarList will be computed (for 
          // free) although they are not requireed to get the edgeTriangleList.
          // if memory usage is an issue, please change these pointers by NULL.
          
          TwoSkeleton twoSkeleton;
          twoSkeleton.setWrapper(this);
          return twoSkeleton.buildEdgeTriangles(
            vertexNumber_, cellNumber_, cellArray_,
            edgeTriangleList_,
            &vertexStarList_,
            &edgeList_, &edgeStarList_,
            &triangleList_, &triangleStarList_, &cellTriangleList_);
        }
        
        return 0;
      }
      
      inline int preprocessTriangles(){

        if(!triangleList_.size()){

          TwoSkeleton twoSkeleton;
          twoSkeleton.setWrapper(this);

          twoSkeleton.buildTriangleList(
            vertexNumber_, cellNumber_, cellArray_,
            &triangleList_, &triangleStarList_);
        }

        return 0;
      }
      
      inline int preprocessTriangleEdges(){
        
        if(!triangleEdgeList_.size()){
        
          // WARNING
          // here triangleStarList and cellTriangleList will be computed (for 
          // free) although they are not requireed to get the edgeTriangleList.
          // if memory usage is an issue, please change these pointers by NULL.
          
          TwoSkeleton twoSkeleton;
          twoSkeleton.setWrapper(this);
          
          return twoSkeleton.buildTriangleEdgeList(
            vertexNumber_, cellNumber_, cellArray_,
            triangleEdgeList_, 
            &vertexEdgeList_,
            &edgeList_,
            &triangleList_, &triangleStarList_,
            &cellTriangleList_);
        }
        
        return 0;
      }
      
      inline int preprocessTriangleLinks(){
        
        if(!triangleLinkList_.size()){
          
          preprocessTriangleStars();
          
          TwoSkeleton twoSkeleton;
          twoSkeleton.setWrapper(this);
          return twoSkeleton.buildTriangleLinks(
            triangleList_, triangleStarList_, cellArray_, triangleLinkList_);
        }
        
        return 0;
      }
      
      inline int preprocessTriangleStars(){
        
        if(!triangleStarList_.size()){
          
          TwoSkeleton twoSkeleton;
          twoSkeleton.setWrapper(this);
          return twoSkeleton.buildTriangleList(vertexNumber_, cellNumber_,
            cellArray_, &triangleList_, &triangleStarList_);
        }
        
        return 0;
      }
      
      inline int preprocessVertexEdges() {
        
        if((int) vertexEdgeList_.size() != vertexNumber_){
          ZeroSkeleton zeroSkeleton;
          
          if(!edgeList_.size()){
            OneSkeleton oneSkeleton;
            oneSkeleton.setWrapper(this);
            oneSkeleton.buildEdgeList(vertexNumber_, cellNumber_,
              cellArray_, edgeList_);
          }
          
          zeroSkeleton.setWrapper(this);
          return zeroSkeleton.buildVertexEdges(vertexNumber_,
            edgeList_, vertexEdgeList_);
        }
        return 0;
      }
      
      inline int preprocessVertexLinks(){
       
        if((int) vertexLinkList_.size() != vertexNumber_){
          
          if(getDimensionality() == 2){
            preprocessVertexStars();
            preprocessCellEdges();
            
            ZeroSkeleton zeroSkeleton;
            zeroSkeleton.setWrapper(this);
            return zeroSkeleton.buildVertexLinks(
              vertexStarList_, cellEdgeList_, edgeList_, vertexLinkList_);
          }
          else if(getDimensionality() == 3){
            preprocessVertexStars();
            preprocessCellTriangles();
            
            ZeroSkeleton zeroSkeleton;
            zeroSkeleton.setWrapper(this);
            return zeroSkeleton.buildVertexLinks(
              vertexStarList_, cellTriangleList_, triangleList_, 
              vertexLinkList_);
          }
          else{
            // unsupported dimension
            stringstream msg;
            msg 
              << "[ExplicitTriangulation] Unsupported dimension for vertex"
              << " link preprocessing."
              << endl;
            dMsg(cerr, msg.str(), infoMsg);
            return -1;
          }
        }
        return 0;
      }
      
      inline int preprocessVertexNeighbors(){
        
        if((int) vertexNeighborList_.size() != vertexNumber_){
          ZeroSkeleton zeroSkeleton;
          zeroSkeleton.setWrapper(this);
          return zeroSkeleton.buildVertexNeighbors(vertexNumber_, cellNumber_,
            cellArray_, vertexNeighborList_, &edgeList_);
        }
        return 0;
      }
      
      inline int preprocessVertexStars(){
       
        if((int) vertexStarList_.size() != vertexNumber_){
          ZeroSkeleton zeroSkeleton;
          zeroSkeleton.setWrapper(this);
     
          return zeroSkeleton.buildVertexStars(vertexNumber_, cellNumber_,
            cellArray_, vertexStarList_);
        }
        return 0;
      }
     
      inline int preprocessVertexTriangles(){
        
        if((int) vertexTriangleList_.size() != vertexNumber_){
       
          preprocessTriangles();

          TwoSkeleton twoSkeleton;
          twoSkeleton.setWrapper(this);

          twoSkeleton.buildVertexTriangles(vertexNumber_, triangleList_, 
            vertexTriangleList_);
        }

        return 0;
      }

      inline int setInputCells(const int &cellNumber, 
        const long long int *cellArray){
        
        if(cellNumber_)
          clear();
        
        cellNumber_ = cellNumber;
        cellArray_ = cellArray;
        
        return 0;
      }
      
      inline int setInputPoints(const int &pointNumber, const float *pointSet){
        
        if(vertexNumber_)
          clear();
        
        vertexNumber_ = pointNumber;
        pointSet_ = pointSet;
        return 0;
      }

      
      
    protected:
    
      int clear();
      
      int                 cellNumber_, vertexNumber_;
      const float         *pointSet_;
      const long long int *cellArray_;
      
  };
}

// if the package is not a template, comment the following line
// #include                  <ExplicitTriangulation.cpp>

#endif // _EXPLICITTRIANGULATION_H
