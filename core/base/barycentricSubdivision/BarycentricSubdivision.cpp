#include "BarycentricSubdivision.h"
#include "Geometry.h"
#include <array>

int ttk::BarycentricSubdivision::subdiviseTriangulation() {

  // not implemented for dimension >= 3
  if(inputTriangl_->getDimensionality() >= 3) {
    std::stringstream msg;
    msg
      << "Barycentric subdivision not yet implemented for dimension 3 and above"
      << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
    return 1;
  }

  const SimplexId newPoints{inputTriangl_->getNumberOfVertices()
                            + inputTriangl_->getNumberOfEdges()
                            + inputTriangl_->getNumberOfTriangles()};
  const size_t dataPerPoint{3};
  points_.clear();
  points_.resize(inputTriangl_->getNumberOfVertices() * dataPerPoint);

  const size_t dataPerCell{4};
  const size_t newTrianglesPerParent{6};
  cells_.clear();
  cells_.reserve(inputTriangl_->getNumberOfTriangles() * dataPerCell
                 * newTrianglesPerParent);

  // copy input points
  std::copy(inputPoints_,
            inputPoints_ + inputTriangl_->getNumberOfVertices() * dataPerPoint,
            points_.begin());

  // reserve memory for new points
  points_.reserve(newPoints * dataPerPoint);

  // generate edge middles
  for(SimplexId i = 0; i < inputTriangl_->getNumberOfEdges(); ++i) {
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

    points_.emplace_back(mid[0]);
    points_.emplace_back(mid[1]);
    points_.emplace_back(mid[2]);
  }

  // generate triangle barycenters
  for(SimplexId i = 0; i < inputTriangl_->getNumberOfTriangles(); ++i) {
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

    points_.emplace_back(bary[0]);
    points_.emplace_back(bary[1]);
    points_.emplace_back(bary[2]);
  }

  // subdivise every triangle
  for(SimplexId i = 0; i < inputTriangl_->getNumberOfTriangles(); ++i) {
    // id of triangle barycenter
    SimplexId bary = inputTriangl_->getNumberOfVertices()
                     + inputTriangl_->getNumberOfEdges() + i;

    for(SimplexId j = 0; j < inputTriangl_->getTriangleEdgeNumber(i); ++j) {
      // edge id
      SimplexId e{};
      inputTriangl_->getTriangleEdge(i, j, e);

      // id of middle of edge e
      SimplexId em = inputTriangl_->getNumberOfVertices() + e;

      // edge vertices
      SimplexId a{}, b{};
      inputTriangl_->getEdgeVertex(e, 0, a);
      inputTriangl_->getEdgeVertex(e, 1, b);

      // new triangles: a, em, bary & b, em, bary
      cells_.emplace_back(3);
      cells_.emplace_back(a);
      cells_.emplace_back(em);
      cells_.emplace_back(bary);
      cells_.emplace_back(3);
      cells_.emplace_back(b);
      cells_.emplace_back(em);
      cells_.emplace_back(bary);
    }
  }

  return 0;
}

int ttk::BarycentricSubdivision::buildOutputTriangulation() {
  // ensure subdivision is already performed
  if(points_.empty() || cells_.empty()) {
    return 1;
  }

  // ensure output triangulation allocated by caller
  if(outputTriangl_ == nullptr) {
    return 2;
  }

  outputTriangl_->setInputPoints(points_.size(), points_.data());
  outputTriangl_->setInputCells(cells_.size(), cells_.data());

  return 0;
}

int ttk::BarycentricSubdivision::execute() {

  Timer t;

  SimplexId vertexNumber = inputTriangl_->getNumberOfVertices();
  subdiviseTriangulation();
  buildOutputTriangulation();

  {
    std::stringstream msg;
    msg << "[BarycentricSubdivision] Data-set (" << vertexNumber
        << " points) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}
