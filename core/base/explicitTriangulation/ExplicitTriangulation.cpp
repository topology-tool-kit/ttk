#include <ExplicitTriangulation.h>

using namespace std;
using namespace ttk;

ExplicitTriangulation::ExplicitTriangulation() {

  setDebugMsgPrefix("ExplicitTriangulation");

  clear();
}

ExplicitTriangulation::~ExplicitTriangulation() {
}

int ExplicitTriangulation::clear() {
  vertexNumber_ = 0;
  cellNumber_ = 0;
  doublePrecision_ = false;

  printMsg("Triangulation cleared.", debug::Priority::DETAIL);

  return AbstractTriangulation::clear();
}

size_t ExplicitTriangulation::footprint(size_t size) const {

  const auto printArrayFootprint
    = [this](const FlatJaggedArray &array, const std::string &name) {
        if(!array.empty() && !name.empty()) {
          this->printMsg(name + std::string{": "}
                         + std::to_string(array.footprint()) + " bytes");
        }
        return array.footprint();
      };

  size += printArrayFootprint(vertexNeighborData_, "vertexNeighborData_");
  size += printArrayFootprint(cellNeighborData_, "cellNeighborData_");
  size += printArrayFootprint(vertexEdgeData_, "vertexEdgeData_");
  size += printArrayFootprint(vertexTriangleData_, "vertexTriangleData_");
  size += printArrayFootprint(edgeTriangleData_, "edgeTriangleData_");
  size += printArrayFootprint(vertexStarData_, "vertexStarData_");
  size += printArrayFootprint(edgeStarData_, "edgeStarData_");
  size += printArrayFootprint(triangleStarData_, "triangleStarData_");
  size += printArrayFootprint(vertexLinkData_, "vertexLinkData_");
  size += printArrayFootprint(edgeLinkData_, "edgeLinkData_");
  size += printArrayFootprint(triangleLinkData_, "triangleLinkData_");

  return AbstractTriangulation::footprint(size);
}
