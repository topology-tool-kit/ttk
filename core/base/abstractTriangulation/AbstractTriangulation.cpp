#include <AbstractTriangulation.h>

using namespace std;
using namespace ttk;

AbstractTriangulation::AbstractTriangulation() {

  clear();
}

AbstractTriangulation::~AbstractTriangulation() {
}

int AbstractTriangulation::clear() {

  hasPreconditionedBoundaryEdges_ = false;
  hasPreconditionedBoundaryTriangles_ = false;
  hasPreconditionedBoundaryVertices_ = false;
  hasPreconditionedCellEdges_ = false;
  hasPreconditionedCellNeighbors_ = false;
  hasPreconditionedCellTriangles_ = false;
  hasPreconditionedEdges_ = false;
  hasPreconditionedEdgeLinks_ = false;
  hasPreconditionedEdgeStars_ = false;
  hasPreconditionedEdgeTriangles_ = false;
  hasPreconditionedTriangles_ = false;
  hasPreconditionedTriangleEdges_ = false;
  hasPreconditionedTriangleLinks_ = false;
  hasPreconditionedTriangleStars_ = false;
  hasPreconditionedVertexEdges_ = false;
  hasPreconditionedVertexLinks_ = false;
  hasPreconditionedVertexNeighbors_ = false;
  hasPreconditionedVertexStars_ = false;
  hasPreconditionedVertexTriangles_ = false;

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

template <class itemType>
size_t AbstractTriangulation::tableTableFootprint(
  const vector<vector<itemType>> &table,
  const string tableName,
  stringstream *msg) const {

  size_t localByteNumber = 0;

  for(size_t i = 0; i < table.size(); i++) {
    localByteNumber += table[i].size() * sizeof(itemType);
  }

  if((localByteNumber) && (tableName.length()) && (msg)) {
    (*msg) << "[AbstractTriangulation] " << tableName << ": " << localByteNumber
           << " bytes" << endl;
  }

  return localByteNumber;
}

size_t AbstractTriangulation::footprint() const {

  size_t size = sizeof(*this);
  stringstream msg;

  size += tableFootprint<bool>(boundaryEdges_, "boundaryEdges_", &msg);

  size += tableFootprint<bool>(boundaryTriangles_, "boundaryTriangles_", &msg);

  size += tableFootprint<bool>(boundaryVertices_, "boundaryVertices_", &msg);

  size += tableTableFootprint<SimplexId>(cellEdgeList_, "cellEdgeList_", &msg);

  size += tableTableFootprint<SimplexId>(
    cellNeighborList_, "cellNeighborList_", &msg);

  size += tableTableFootprint<SimplexId>(
    cellTriangleList_, "cellTriangleList_", &msg);

  size += tableTableFootprint<SimplexId>(edgeLinkList_, "edgeLinkList_", &msg);

  size
    += tableFootprint<pair<SimplexId, SimplexId>>(edgeList_, "edgeList_", &msg);

  size += tableTableFootprint<SimplexId>(edgeStarList_, "edgeStarList_", &msg);

  size += tableTableFootprint<SimplexId>(
    edgeTriangleList_, "edgeTriangleList_", &msg);

  size += tableTableFootprint<SimplexId>(triangleList_, "triangleList_", &msg);

  size += tableTableFootprint<SimplexId>(
    triangleEdgeList_, "triangleEdgeList_", &msg);

  size += tableTableFootprint<SimplexId>(
    triangleLinkList_, "triangleLinkList_", &msg);

  size += tableTableFootprint<SimplexId>(
    triangleStarList_, "triangleStarList_", &msg);

  size
    += tableTableFootprint<SimplexId>(vertexEdgeList_, "vertexEdgeList_", &msg);

  size
    += tableTableFootprint<SimplexId>(vertexLinkList_, "vertexLinkList_", &msg);

  size += tableTableFootprint<SimplexId>(
    vertexNeighborList_, "vertexNeighborList_", &msg);

  size
    += tableTableFootprint<SimplexId>(vertexStarList_, "vertexStarList_", &msg);

  size += tableTableFootprint<SimplexId>(
    vertexTriangleList_, "vertexTriangleList_", &msg);

  msg << "[AbstractTriangulation] Total footprint: " << (size / 1024) / 1024
      << " MB." << endl;

  dMsg(cout, msg.str(), memoryMsg);

  return size;
}
