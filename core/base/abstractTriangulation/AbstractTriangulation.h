/// \ingroup base
/// \class ttk::AbstractTriangulation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date January 2016.
/// 
/// \brief AbstractTriangulation is a virtual class that defines an interface 
/// for efficient traversal methods on triangulations of piecewise linear 
/// manifolds. 
///
/// \sa Triangulation

#ifndef _ABSTRACTTRIANGULATION_H
#define _ABSTRACTTRIANGULATION_H

// base code includes
#include                  <Wrapper.h>

namespace ttk{
  
  class AbstractTriangulation : public Wrapper{

    
    public:
        
      AbstractTriangulation();
      
      ~AbstractTriangulation();

      virtual int clear();
      
      virtual size_t footprint() const;
      
      virtual int getCellEdge(const CellId &cellId, 
        const int &localEdgeId, EdgeId &edgeId) const = 0;
        
      virtual int getCellEdgeNumber(const CellId &cellId) const = 0;
      
      virtual const std::vector<std::vector<EdgeId> > *getCellEdges() = 0;
      
      virtual int getCellNeighbor(const CellId &cellId,
        const int &localNeighborId, CellId &neighborId) const = 0;
        
      virtual int getCellNeighborNumber(const CellId &cellId) const = 0;
      
      virtual const std::vector<std::vector<CellId> > *getCellNeighbors() = 0;
      
      virtual int getCellTriangle(const CellId &cellId, 
        const int &localTriangleId, TriangleId &triangleId) const = 0;
        
      virtual int getCellTriangleNumber(const CellId &cellId) const = 0;
        
      virtual const std::vector<std::vector<TriangleId> > *getCellTriangles() = 0;
      
      virtual int getCellVertex(const CellId &cellId,
        const int &localVertexId, VertexId &vertexId) const = 0;
    
      virtual int getCellVertexNumber(const CellId &cellId) const = 0;
        
      virtual int getDimensionality() const = 0;
      
      virtual const std::vector<std::pair<VertexId, VertexId> > *getEdges() = 0;
        
      virtual int getEdgeLink(const EdgeId &edgeId, 
        const int &localLinkId, LinkId &linkId) const = 0;
        
      virtual int getEdgeLinkNumber(const EdgeId &edgeId) const = 0;
      
      virtual const std::vector<std::vector<LinkId> > *getEdgeLinks() = 0;
      
      virtual int getEdgeStar(const EdgeId &edgeId, 
        const int &localStarId, CellId &starId) const = 0;
        
      virtual int getEdgeStarNumber(const EdgeId &edgeId) const = 0;
      
      virtual const std::vector<std::vector<CellId> > *getEdgeStars() = 0;
     
      virtual int getEdgeTriangle(const EdgeId &edgeId, 
        const int &localTriangleId, TriangleId &triangleId) const = 0;
        
      virtual int getEdgeTriangleNumber(const EdgeId &edgeId) const = 0;
        
      virtual const std::vector<std::vector<TriangleId> > *getEdgeTriangles() = 0;
      
      virtual int getEdgeVertex(const EdgeId &edgeId, 
        const int &localVertexId, VertexId &vertexId) const = 0;
      
      virtual CellId getNumberOfCells() const = 0;
      
      virtual EdgeId getNumberOfEdges() const = 0;
      
      virtual TriangleId getNumberOfTriangles() const = 0;
      
      virtual VertexId getNumberOfVertices() const = 0;
      
      virtual const std::vector<std::vector<VertexId> > *getTriangles() = 0;
      
      virtual int getTriangleEdge(const TriangleId &triangleId,
        const int &localEdgeId, EdgeId &edgeId) const = 0;
      
      virtual int getTriangleEdgeNumber(const TriangleId &triangleId) const = 0;
      
      virtual const std::vector<std::vector<EdgeId> > *getTriangleEdges() = 0;
      
      virtual int getTriangleLink(const TriangleId &triangleId, 
        const int &localLinkId, LinkId &linkId) const = 0;
        
      virtual int getTriangleLinkNumber(const TriangleId &triangleId) const = 0;
      
      virtual const std::vector<std::vector<LinkId> > *getTriangleLinks() = 0;
      
      virtual int getTriangleStar(const TriangleId &triangleId,
        const int &localStarId, CellId &starId) const = 0;  
        
      virtual int getTriangleStarNumber(const TriangleId &triangleId) const = 0;
      
      virtual const std::vector<std::vector<CellId> > *getTriangleStars() = 0;
      
      virtual int getTriangleVertex(const TriangleId &triangleId,
        const int &localVertexId, VertexId &vertexId) const = 0;
      
      virtual int getVertexEdge(const VertexId &vertexId, 
        const int &localEdgeId, EdgeId &edgeId) const = 0;
        
      virtual int getVertexEdgeNumber(const VertexId &vertexId) const = 0;
      
      virtual const std::vector<std::vector<EdgeId> > *getVertexEdges() = 0;
      
      virtual int getVertexLink(const VertexId &vertexId, 
        const int &localLinkId, LinkId &linkId) const = 0;
        
      virtual int getVertexLinkNumber(const VertexId &vertexId) const = 0;
      
      virtual const std::vector<std::vector<LinkId> > *getVertexLinks() = 0;
      
      virtual int getVertexNeighbor(const VertexId &vertexId, 
        const int &localNeighborId, VertexId &neighborId) const = 0;
        
      virtual int getVertexNeighborNumber(const VertexId &vertexId) const = 0;
      
      virtual const std::vector<std::vector<VertexId> > *getVertexNeighbors() = 0;
      
      virtual int getVertexPoint(const VertexId &vertexId,
        float &x, float &y, float &z) const = 0;
        
      virtual int getVertexStar(const VertexId &vertexId, const int &localStarId,
        CellId &starId) const = 0;
        
      virtual int getVertexStarNumber(const VertexId &vertexId) const = 0;
        
      virtual const std::vector<std::vector<CellId> > *getVertexStars() = 0;
      
      virtual int getVertexTriangle(const VertexId &vertexId, 
        const int &localTriangleId, TriangleId &triangleId) const = 0;
        
      virtual int getVertexTriangleNumber(const VertexId &vertexId) const = 0;
        
      virtual const std::vector<std::vector<TriangleId> > *getVertexTriangles() = 0;
        
      virtual inline bool hasPreprocessedBoundaryEdges() const{
        return hasPreprocessedBoundaryEdges_;
      }
      
      virtual inline bool hasPreprocessedBoundaryTriangles() const{
        return hasPreprocessedBoundaryTriangles_;
      }
       
      virtual inline bool hasPreprocessedBoundaryVertices() const{
        return hasPreprocessedBoundaryVertices_;
      }
     
      virtual inline bool hasPreprocessedCellEdges() const{
        return hasPreprocessedCellEdges_;
      }
      
      virtual inline bool hasPreprocessedCellNeighbors() const{
        return hasPreprocessedCellNeighbors_;
      }
      
      virtual inline bool hasPreprocessedCellTriangles() const{
        return hasPreprocessedCellTriangles_;
      }
       
      virtual inline bool hasPreprocessedEdgeLinks() const{
        return hasPreprocessedEdgeLinks_;
      }
       
      virtual inline bool hasPreprocessedEdgeStars() const{
        return hasPreprocessedEdgeStars_;
      }
      
      virtual inline bool hasPreprocessedEdgeTriangles() const{
        return hasPreprocessedEdgeTriangles_;
      }
        
      virtual inline bool hasPreprocessedEdges() const{
        return hasPreprocessedEdges_;
      }
      
      virtual inline bool hasPreprocessedTriangles() const{
        return hasPreprocessedTriangles_;
      }
      
      virtual inline bool hasPreprocessedTriangleEdges() const{
        return hasPreprocessedTriangleEdges_;
      }
      
      virtual inline bool hasPreprocessedTriangleLinks() const{
        return hasPreprocessedTriangleLinks_;
      }
      
      virtual inline bool hasPreprocessedTriangleStars() const{
        return hasPreprocessedTriangleStars_;
      }
        
      virtual inline bool hasPreprocessedVertexEdges() const{
        return hasPreprocessedVertexEdges_;
      }
      
      virtual inline bool hasPreprocessedVertexLinks() const{
        return hasPreprocessedVertexLinks_;
      }
      
      virtual inline bool hasPreprocessedVertexNeighbors() const{
        return hasPreprocessedVertexNeighbors_;
      }
      
      virtual inline bool hasPreprocessedVertexStars() const{
        return hasPreprocessedVertexStars_;
      }
      
      virtual inline bool hasPreprocessedVertexTriangles() const{
        return hasPreprocessedVertexTriangles_;
      }
      
      virtual bool isEdgeOnBoundary(const EdgeId &edgeId) const = 0;
        
      virtual bool isEmpty() const = 0;
      
      virtual bool isTriangleOnBoundary(const TriangleId &triangleId) const = 0;
      
      virtual bool isVertexOnBoundary(const VertexId &vertexId) const = 0;

      virtual int preprocessBoundaryEdges(){
        preprocessEdges();
        hasPreprocessedBoundaryEdges_ = true;
        return 0;
      }
      
      virtual int preprocessBoundaryTriangles(){
        preprocessTriangles();
        hasPreprocessedBoundaryTriangles_ = true;
        return 0;
      }
      
      virtual int preprocessBoundaryVertices(){
        hasPreprocessedBoundaryVertices_ = true;
        return 0;
      }
      
      virtual int preprocessCellEdges(){
        preprocessEdges();
        hasPreprocessedCellEdges_ = true;
        return 0;
      }
      
      virtual int preprocessCellNeighbors(){
        hasPreprocessedCellNeighbors_ = true;
        return 0;
      }
      
      virtual int preprocessCellTriangles(){
        preprocessTriangles();
        hasPreprocessedCellTriangles_ = true;
        return 0;
      }
      
      virtual int preprocessEdges(){
        hasPreprocessedEdges_ = true;
        return 0;
      }
      
      virtual int preprocessEdgeLinks(){
        preprocessEdges();
        hasPreprocessedEdgeLinks_ = true;
        return 0;
      }
      
      virtual int preprocessEdgeStars(){
        preprocessEdges();
        hasPreprocessedEdgeStars_ = true;
        return 0;
      }
      
      virtual int preprocessEdgeTriangles(){
        preprocessEdges();
        preprocessTriangles();
        hasPreprocessedEdgeTriangles_ = true;
        return 0;
      }
      
      virtual int preprocessTriangles(){
        hasPreprocessedTriangles_ = true;
        return 0;
      }
      
      virtual int preprocessTriangleEdges(){
        preprocessEdges();
        preprocessTriangles();
        hasPreprocessedTriangleEdges_ = true;
        return 0;
      }
      
      virtual int preprocessTriangleLinks(){
        preprocessTriangles();
        hasPreprocessedTriangleLinks_ = true;
        return 0;
      }
      
      virtual int preprocessTriangleStars(){
        preprocessTriangles();
        hasPreprocessedTriangleStars_ = true;
        return 0;
      }
      
      virtual int preprocessVertexEdges() { 
        preprocessEdges();
        hasPreprocessedVertexEdges_ = true;
        return 0;
      }
      
      virtual int preprocessVertexLinks(){
        hasPreprocessedVertexLinks_ = true;
        return 0;
      }
      
      virtual int preprocessVertexNeighbors(){
        hasPreprocessedVertexNeighbors_ = true;
        return 0;
      }
      
      virtual int preprocessVertexStars(){
        hasPreprocessedVertexStars_ = true;
        return 0;
      }
      
      virtual int preprocessVertexTriangles(){
        preprocessTriangles();
        hasPreprocessedVertexTriangles_ = true;
        return 0;
      }

    protected:
      
      // empty wrapping to VTK for now
      bool needsToAbort(){ return false;};
      
      template <class itemType>
        int tableFootprint(const std::vector<itemType> &table,
          const std::string tableName = "", 
          std::stringstream *msg = NULL) const{
        
        if((table.size())&&(tableName.length())&&(msg)){
          (*msg) << "[AbstractTriangulation] " << tableName << ": "
            << table.size()*sizeof(itemType) << " bytes" << std::endl;
        }
            
        return table.size()*sizeof(itemType);
      }
      
      template <class itemType>
        size_t tableTableFootprint(
          const std::vector<std::vector<itemType> > &table,
          const std::string tableName = "", 
          std::stringstream *msg = NULL) const;
      
      int updateProgress(const float &progress) {return 0;};
      
      bool                hasPreprocessedBoundaryEdges_,
                          hasPreprocessedBoundaryTriangles_,
                          hasPreprocessedBoundaryVertices_,
                          hasPreprocessedCellEdges_,
                          hasPreprocessedCellNeighbors_,
                          hasPreprocessedCellTriangles_,
                          hasPreprocessedEdges_,
                          hasPreprocessedEdgeLinks_,
                          hasPreprocessedEdgeStars_,
                          hasPreprocessedEdgeTriangles_,
                          hasPreprocessedTriangles_,
                          hasPreprocessedTriangleEdges_,
                          hasPreprocessedTriangleLinks_,
                          hasPreprocessedTriangleStars_,
                          hasPreprocessedVertexEdges_,
                          hasPreprocessedVertexLinks_,
                          hasPreprocessedVertexNeighbors_,
                          hasPreprocessedVertexStars_,
                          hasPreprocessedVertexTriangles_;
      
      std::vector<bool>        boundaryEdges_,
                          boundaryTriangles_,
                          boundaryVertices_;
      std::vector<std::vector<EdgeId> > 
                          cellEdgeList_;
      std::vector<std::vector<CellId> >
                          cellNeighborList_;
      std::vector<std::vector<TriangleId> > 
                          cellTriangleList_;
      std::vector<std::vector<LinkId> >
                          edgeLinkList_;
      std::vector<std::pair<VertexId, VertexId> >
                          edgeList_;
      std::vector<std::vector<CellId> >
                          edgeStarList_;
      std::vector<std::vector<Triangle> >
                          edgeTriangleList_;
      std::vector<std::vector<VertexId> >
                          triangleList_;
      std::vector<std::vector<EdgeId> > 
                          triangleEdgeList_;
      std::vector<std::vector<LinkId> >
                          triangleLinkList_;
      std::vector<std::vector<CellId> >
                          triangleStarList_;
      std::vector<std::vector<EdgeId> > 
                          vertexEdgeList_;
      std::vector<std::vector<Link> >
                          vertexLinkList_;
      std::vector<std::vector<VertexId> > 
                          vertexNeighborList_;
      std::vector<std::vector<CellId> >
                          vertexStarList_;
      std::vector<std::vector<TriangleId> >
                          vertexTriangleList_;
  };
}

// if the package is not a template, comment the following line
// #include                  <AbstractTriangulation.cpp>

#endif // _ABSTRACTTRIANGULATION_H
