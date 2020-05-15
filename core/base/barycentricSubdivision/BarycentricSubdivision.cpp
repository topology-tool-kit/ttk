#include "BarycentricSubdivision.h"
#include "Geometry.h"

#include <array>

#define MODULE_S "[BarycentricSubdivision] "

int ttk::BarycentricSubdivision::subdiviseTriangulation() {

  // not implemented for dimension >= 3
  if(inputTriangl_->getDimensionality() >= 3) {
    std::stringstream msg;
    msg << MODULE_S "Not yet implemented for dimension 3 and above"
        << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
    return 1;
  }

  const SimplexId newPoints{nVertices_ + nEdges_ + nTriangles_};
  const size_t dataPerPoint{3};
  points_.clear();
  points_.resize(newPoints * dataPerPoint);

  const size_t dataPerCell{3};
  const size_t newTrianglesPerParent{6};
  cells_connectivity_.clear();
  cells_offsets_.clear();
  cells_connectivity_.resize(nTriangles_ * newTrianglesPerParent * dataPerCell);
  cells_offsets_.resize(nTriangles_ * newTrianglesPerParent + 1);

  pointId_.clear();
  pointId_.resize(newPoints);

  pointDim_.clear();
  pointDim_.resize(newPoints);

  // copy input points
  std::copy(
    inputPoints_, inputPoints_ + nVertices_ * dataPerPoint, points_.begin());

  // set input point ids
  for(SimplexId i = 0; i < nVertices_; ++i) {
    pointId_[i] = i;
  }

  // reserve memory for new points
  points_.reserve(newPoints * dataPerPoint);

  // generate edge middles
  for(SimplexId i = 0; i < nEdges_; ++i) {
    // edge vertices
    SimplexId a{}, b{};
    inputTriangl_->getEdgeVertex(i, 0, a);
    inputTriangl_->getEdgeVertex(i, 1, b);

    std::array<float, 3> pa{}, pb{}, mid{};
    inputTriangl_->getVertexPoint(a, pa[0], pa[1], pa[2]);
    inputTriangl_->getVertexPoint(b, pb[0], pb[1], pb[2]);

    mid[0] = (pa[0] + pb[0]) / 2.0F;
    mid[1] = (pa[1] + pb[1]) / 2.0F;
    mid[2] = (pa[2] + pb[2]) / 2.0F;

    const size_t offset = dataPerPoint * (nVertices_ + i);
    points_[offset + 0] = mid[0];
    points_[offset + 1] = mid[1];
    points_[offset + 2] = mid[2];
    pointId_[nVertices_ + i] = i;
    pointDim_[nVertices_ + i] = 1;
  }

  // generate triangle barycenters
  for(SimplexId i = 0; i < nTriangles_; ++i) {
    // triangle vertices
    SimplexId a{}, b{}, c{};
    inputTriangl_->getTriangleVertex(i, 0, a);
    inputTriangl_->getTriangleVertex(i, 1, b);
    inputTriangl_->getTriangleVertex(i, 2, c);

    std::array<float, 3> pa{}, pb{}, pc{}, bary{};
    inputTriangl_->getVertexPoint(a, pa[0], pa[1], pa[2]);
    inputTriangl_->getVertexPoint(b, pb[0], pb[1], pb[2]);
    inputTriangl_->getVertexPoint(c, pc[0], pc[1], pc[2]);

    bary[0] = (pa[0] + pb[0] + pc[0]) / 3.0F;
    bary[1] = (pa[1] + pb[1] + pc[1]) / 3.0F;
    bary[2] = (pa[2] + pb[2] + pc[2]) / 3.0F;

    const size_t offset = dataPerPoint * (nVertices_ + nEdges_ + i);
    points_[offset + 0] = bary[0];
    points_[offset + 1] = bary[1];
    points_[offset + 2] = bary[2];
    pointId_[nVertices_ + nEdges_ + i] = i;
    pointDim_[nVertices_ + nEdges_ + i] = 2;
  }

  LongSimplexId off_id = 0, co_id = 0;
  // subdivise every triangle
  for(SimplexId i = 0; i < nTriangles_; ++i) {
    // id of triangle barycenter
    SimplexId bary = nVertices_ + nEdges_ + i;

    for(SimplexId j = 0; j < inputTriangl_->getTriangleEdgeNumber(i); ++j) {
      // edge id
      SimplexId e{};
      inputTriangl_->getTriangleEdge(i, j, e);

      // id of middle of edge e
      SimplexId em = nVertices_ + e;

      // edge vertices
      SimplexId a{}, b{};
      inputTriangl_->getEdgeVertex(e, 0, a);
      inputTriangl_->getEdgeVertex(e, 1, b);

      // Ad triangle a - em - bary
      cells_offsets_[off_id] = co_id;
      off_id++;
      cells_connectivity_[co_id + 0] = a;
      cells_connectivity_[co_id + 1] = em;
      cells_connectivity_[co_id + 2] = bary;
      co_id += 3;

      // Ad triangle b - em - bary
      cells_offsets_[off_id] = co_id;
      off_id++;
      cells_connectivity_[co_id + 0] = b;
      cells_connectivity_[co_id + 1] = em;
      cells_connectivity_[co_id + 2] = bary;
      co_id += 3;
    }
  }
  // Last offset
  cells_offsets_[off_id] = co_id;

  return 0;
}

int ttk::BarycentricSubdivision::buildOutputTriangulation() {
  // ensure subdivision is already performed
  if(points_.empty() || cells_connectivity_.empty()) {
    return 1;
  }

  // ensure output triangulation allocated by caller
  if(outputTriangl_ == nullptr) {
    return 2;
  }

  outputTriangl_->setInputPoints(points_.size() / 3, points_.data());
#ifdef CELL_ARRAY_NEW
  outputTriangl_->setInputCells(cells_offsets_.size() - 1,
                                cells_connectivity_.data(),
                                cells_offsets_.data());
#else
  // TODO use revert translator
  LongSimplexId* cells = nullptr;
  CellArray::TranslateToFlatLayout(cells_connectivity_, cells_offsets_, cells);
  outputTriangl_->setInputCells(cells_offsets_.size() - 1, cells);
#endif

  return 0;
}

int ttk::BarycentricSubdivision::execute() {

  Timer t;

  SimplexId vertexNumber = inputTriangl_->getNumberOfVertices();
  subdiviseTriangulation();
  buildOutputTriangulation();

  {
    std::stringstream msg;
    msg << MODULE_S "Data-set (" << vertexNumber << " points) processed in "
        << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}
