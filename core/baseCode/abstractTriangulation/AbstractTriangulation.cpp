#include                  <AbstractTriangulation.h>

AbstractTriangulation::AbstractTriangulation(){
  
  clear();
}

AbstractTriangulation::~AbstractTriangulation(){

}

int AbstractTriangulation::clear(){
  
  hasPreprocessedBoundaryEdges_ = false;
  hasPreprocessedBoundaryTriangles_ = false;
  hasPreprocessedBoundaryVertices_ = false;
  hasPreprocessedCellEdges_ = false;
  hasPreprocessedCellNeighbors_ = false;
  hasPreprocessedCellTriangles_ = false;
  hasPreprocessedEdges_ = false;
  hasPreprocessedEdgeLinks_ = false;
  hasPreprocessedEdgeStars_ = false;
  hasPreprocessedEdgeTriangles_ = false;
  hasPreprocessedTriangles_ = false;
  hasPreprocessedTriangleEdges_ = false;
  hasPreprocessedTriangleLinks_ = false;
  hasPreprocessedTriangleStars_ = false;
  hasPreprocessedVertexEdges_ = false;
  hasPreprocessedVertexLinks_ = false;
  hasPreprocessedVertexNeighbors_ = false;
  hasPreprocessedVertexStars_ = false;
  hasPreprocessedVertexTriangles_ = false;
  
  boundaryEdges_.clear();
  boundaryTriangles_.clear();
  boundaryVertices_.clear();
  
  cellEdgeList_.clear();
  cellNeighborList_.clear();
  cellTriangleList_.clear();
  
  edgeLinkList_.clear();
  edgeList_.clear();
  edgeStarList_.clear();
  edgeTriangleList_.clear();
  
  triangleList_.clear();
  triangleEdgeList_.clear();
  triangleLinkList_.clear();
  triangleStarList_.clear();
  
  vertexEdgeList_.clear();
  vertexLinkList_.clear();
  vertexNeighborList_.clear();
  vertexStarList_.clear();
  vertexTriangleList_.clear();
  
  return 0;
}