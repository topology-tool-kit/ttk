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

template <class itemType> int AbstractTriangulation::tableTableFootprint(
  const vector<vector<itemType> > &table,
  const string tableName,
  stringstream *msg) const{
  
  int localByteNumber = 0;
  
  for(int i = 0; i < (int) table.size(); i++){
    localByteNumber += table[i].size()*sizeof(itemType);
  }

  if((localByteNumber)&&(tableName.length())&&(msg)){
    (*msg) << "[AbstractTriangulation] " << tableName << ": "
      << localByteNumber << " bytes" << endl;
  }
  
  return localByteNumber;
}

int AbstractTriangulation::footprint() const{

  int size = sizeof(*this);
  stringstream msg;
  
  size += tableFootprint<bool>(boundaryEdges_, "boundaryEdges_", &msg);
 
  size += tableFootprint<bool>(boundaryTriangles_, "boundaryTriangles_", &msg);
  
  size += tableFootprint<bool>(boundaryVertices_, "boundaryVertices_", &msg);
  
  size += tableTableFootprint<int>(cellEdgeList_, "cellEdgeList_", &msg);
  
  size += 
    tableTableFootprint<int>(cellNeighborList_, "cellNeighborList_", &msg);
    
  size +=
    tableTableFootprint<int>(cellTriangleList_, "cellTriangleList_", &msg);
 
  size +=
    tableTableFootprint<int>(edgeLinkList_, "edgeLinkList_", &msg);
 
  size +=
    tableFootprint<pair<int, int>>(edgeList_, "edgeList_", &msg);
    
  size += 
    tableTableFootprint<int>(edgeStarList_, "edgeStarList_", &msg);
    
  size +=
    tableTableFootprint<int>(edgeTriangleList_, "edgeTriangleList_", &msg);
    
  size +=
    tableTableFootprint<int>(triangleList_, "triangleList_", &msg);
    
  size +=
    tableTableFootprint<int>(triangleEdgeList_, "triangleEdgeList_", &msg);
  
  size +=
    tableTableFootprint<int>(triangleLinkList_, "triangleLinkList_", &msg);
    
  size +=
    tableTableFootprint<int>(triangleStarList_, "triangleStarList_", &msg);
    
  size +=
    tableTableFootprint<int>(vertexEdgeList_, "vertexEdgeList_", &msg);
  
  size +=
    tableTableFootprint<int>(vertexLinkList_, "vertexLinkList_", &msg);
    
  size +=
    tableTableFootprint<int>(vertexNeighborList_, "vertexNeighborList_", &msg);
    
  size +=
    tableTableFootprint<int>(vertexStarList_, "vertexStarList_", &msg);
    
  size +=
    tableTableFootprint<int>(vertexTriangleList_, "vertexTriangleList_", &msg);
    
  msg << "[AbstractTriangulation] Total footprint: "
    << (size/1024)/1024 << " MB." << endl;
 
  dMsg(cout, msg.str(), memoryMsg);
  
  return size;
}
