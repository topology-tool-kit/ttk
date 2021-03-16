#include <AbstractTriangulation.h>

using namespace std;
using namespace ttk;

AbstractTriangulation::AbstractTriangulation() {

  setDebugMsgPrefix("AbstractTriangulation");

  clear();
}

AbstractTriangulation::~AbstractTriangulation() {
}

int AbstractTriangulation::clear() {

  hasPeriodicBoundaries_ = false;
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

  tetraEdgeList_.clear();
  cellNeighborList_.clear();
  tetraTriangleList_.clear();

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
  const string &tableName,
  ostream &stream) const {

  size_t localByteNumber = 0;
  stringstream msg;

  for(size_t i = 0; i < table.size(); i++) {
    localByteNumber += table[i].size() * sizeof(itemType);
  }

  if((localByteNumber) && (tableName.length()) && (msg)) {
    msg << tableName << ": " << localByteNumber << " bytes";
    printMsg(msg.str(), debug::Priority::INFO, debug::LineMode::NEW, stream);
  }

  return localByteNumber;
}

size_t AbstractTriangulation::footprint(size_t size) const {

  size += sizeof(*this);
  stringstream msg;

  size += tableFootprint<bool>(boundaryEdges_, "boundaryEdges_");

  size += tableFootprint<bool>(boundaryTriangles_, "boundaryTriangles_");

  size += tableFootprint<bool>(boundaryVertices_, "boundaryVertices_");

  size += tableFootprint(tetraEdgeList_, "tetraEdgeList_");

  size
    += tableTableFootprint<SimplexId>(cellNeighborList_, "cellNeighborList_");

  size += tableFootprint(tetraTriangleList_, "tetraTriangleList_");

  size += tableTableFootprint<SimplexId>(edgeLinkList_, "edgeLinkList_");

  size += tableFootprint(edgeList_, "edgeList_");

  size += tableTableFootprint<SimplexId>(edgeStarList_, "edgeStarList_");

  size
    += tableTableFootprint<SimplexId>(edgeTriangleList_, "edgeTriangleList_");

  size += tableFootprint(triangleList_, "triangleList_");

  size += tableFootprint(triangleEdgeList_, "triangleEdgeList_");

  size
    += tableTableFootprint<SimplexId>(triangleLinkList_, "triangleLinkList_");

  size
    += tableTableFootprint<SimplexId>(triangleStarList_, "triangleStarList_");

  size += tableTableFootprint<SimplexId>(vertexEdgeList_, "vertexEdgeList_");

  size += tableTableFootprint<SimplexId>(vertexLinkList_, "vertexLinkList_");

  size += tableTableFootprint<SimplexId>(
    vertexNeighborList_, "vertexNeighborList_");

  size += tableTableFootprint<SimplexId>(vertexStarList_, "vertexStarList_");

  size += tableTableFootprint<SimplexId>(
    vertexTriangleList_, "vertexTriangleList_");

  size += tableTableFootprint(cellEdgeVector_, "cellEdgeVector_");
  size += tableTableFootprint(cellTriangleVector_, "cellTriangleVector_");
  size += tableTableFootprint(triangleEdgeVector_, "triangleEdgeVector_");

  msg << "Total footprint: " << (size / 1024) / 1024 << " MB.";
  printMsg(msg.str());

  return size;
}
