#include <ExplicitTriangulation.h>
#include <OneSkeleton.h>
#include <ThreeSkeleton.h>
#include <TwoSkeleton.h>
#include <ZeroSkeleton.h>

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

int ExplicitTriangulation::preconditionBoundaryEdgesInternal() {
  if((!boundaryEdges_.empty()) && (boundaryEdges_.size() == edgeList_.size())) {
    return 0;
  }

  Timer tm{};

  preconditionEdgesInternal();
  boundaryEdges_.resize(edgeList_.size(), false);

  if(getDimensionality() == 2) {
    preconditionEdgeStarsInternal();
    for(SimplexId i = 0; i < (SimplexId)edgeStarData_.subvectorsNumber(); i++) {
      if(edgeStarData_.size(i) == 1) {
        boundaryEdges_[i] = true;
      }
    }
  } else if(getDimensionality() == 3) {
    preconditionTriangleStarsInternal();
    preconditionTriangleEdgesInternal();

    for(size_t i = 0; i < triangleStarData_.subvectorsNumber(); i++) {
      if(triangleStarData_.size(i) == 1) {
        for(int j = 0; j < 3; j++) {
          boundaryEdges_[triangleEdgeList_[i][j]] = true;
        }
      }
    }
  } else {
    // unsupported dimension
    printErr("Unsupported dimension for boundary precondition");
    return -1;
  }

  this->printMsg("Extracted boundary edges", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionBoundaryTrianglesInternal() {
  if(getDimensionality() == 2)
    return 0;

  Timer tm{};

  if((!boundaryTriangles_.empty())
     && (boundaryTriangles_.size() == triangleList_.size())) {
    return 0;
  }

  preconditionTrianglesInternal();
  boundaryTriangles_.resize(triangleList_.size(), false);

  if(getDimensionality() == 3) {
    preconditionTriangleStarsInternal();

    for(size_t i = 0; i < triangleStarData_.subvectorsNumber(); i++) {
      if(triangleStarData_.size(i) == 1) {
        boundaryTriangles_[i] = true;
      }
    }
  } else {
    // unsupported dimension
    printErr("Unsupported dimension for boundary precondition");
    return -1;
  }

  this->printMsg("Extracted boundary triangles", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionBoundaryVerticesInternal() {
  if((!boundaryVertices_.empty())
     && ((SimplexId)boundaryVertices_.size() == vertexNumber_))
    return 0;

  Timer tm{};

  boundaryVertices_.resize(vertexNumber_, false);

  // create the list of boundary elements
  // create their star
  // look for singletons
  if(getDimensionality() == 1) {
    preconditionVertexStarsInternal();
    for(size_t i = 0; i < vertexStarData_.subvectorsNumber(); i++) {
      if(vertexStarData_.size(i) == 1) {
        boundaryVertices_[i] = true;
      }
    }
  } else if(getDimensionality() == 2) {
    preconditionEdgesInternal();
    preconditionEdgeStarsInternal();

    for(SimplexId i = 0; i < (SimplexId)edgeStarData_.subvectorsNumber(); i++) {
      if(edgeStarData_.size(i) == 1) {
        boundaryVertices_[edgeList_[i][0]] = true;
        boundaryVertices_[edgeList_[i][1]] = true;
      }
    }
  } else if(getDimensionality() == 3) {
    preconditionTrianglesInternal();
    preconditionTriangleStarsInternal();

    for(size_t i = 0; i < triangleStarData_.subvectorsNumber(); i++) {
      if(triangleStarData_.size(i) == 1) {
        boundaryVertices_[triangleList_[i][0]] = true;
        boundaryVertices_[triangleList_[i][1]] = true;
        boundaryVertices_[triangleList_[i][2]] = true;
      }
    }
  } else {
    // unsupported dimension
    printErr("Unsupported dimension for boundary precondition");
    return -1;
  }

  this->printMsg("Extracted boundary vertices", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionCellEdgesInternal() {
  OneSkeleton os;
  os.setWrapper(this);

  if(tetraEdgeList_.empty() && getDimensionality() == 3) {
    os.buildEdgeList(
      vertexNumber_, *cellArray_, nullptr, nullptr, &tetraEdgeList_);
  } else if(triangleEdgeList_.empty() && getDimensionality() == 2) {
    os.buildEdgeList(
      vertexNumber_, *cellArray_, nullptr, nullptr, &triangleEdgeList_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionCellNeighborsInternal() {

  if(cellNeighborData_.empty()) {
    if(getDimensionality() == 3) {
      ThreeSkeleton threeSkeleton;
      threeSkeleton.setWrapper(this);
      threeSkeleton.buildCellNeighborsFromTriangles(
        vertexNumber_, *cellArray_, cellNeighborData_, &triangleStarData_);
    } else if(getDimensionality() == 2) {
      TwoSkeleton twoSkeleton;
      twoSkeleton.setWrapper(this);
      twoSkeleton.buildCellNeighborsFromEdges(
        vertexNumber_, *cellArray_, cellNeighborData_, &edgeStarData_);
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionCellTrianglesInternal() {

  if(!tetraTriangleList_.size()) {

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    if(triangleList_.size()) {
      // we already computed this guy, let's just get the cell triangles
      if(!triangleStarData_.empty()) {
        return twoSkeleton.buildTriangleList(
          vertexNumber_, *cellArray_, nullptr, nullptr, &tetraTriangleList_);
      } else {
        // let's compute the triangle star while we're at it...
        // it's just a tiny overhead.
        return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                             nullptr, &triangleStarData_,
                                             &tetraTriangleList_);
      }
    } else {
      // we have not computed this guy, let's do it while we're at it
      if(!triangleStarData_.empty()) {
        return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                             &triangleList_, nullptr,
                                             &tetraTriangleList_);
      } else {
        // let's compute the triangle star while we're at it...
        // it's just a tiny overhead.
        return twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_,
                                             &triangleList_, &triangleStarData_,
                                             &tetraTriangleList_);
      }
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionEdgesInternal() {

  if(!edgeList_.size()) {
    OneSkeleton oneSkeleton;
    oneSkeleton.setWrapper(this);
    // also computes edgeStar and triangleEdge / tetraEdge lists for free...
    if(getDimensionality() == 2) {
      return oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_, &edgeList_,
                                       &edgeStarData_, &triangleEdgeList_);
    } else if(getDimensionality() == 3) {
      return oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_, &edgeList_,
                                       &edgeStarData_, &tetraEdgeList_);
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionEdgeLinksInternal() {

  if(edgeLinkData_.empty()) {

    if(getDimensionality() == 2) {
      preconditionEdgesInternal();
      preconditionEdgeStarsInternal();

      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      return oneSkeleton.buildEdgeLinks(
        edgeList_, edgeStarData_, *cellArray_, edgeLinkData_);
    } else if(getDimensionality() == 3) {
      preconditionEdgesInternal();
      preconditionEdgeStarsInternal();
      preconditionCellEdgesInternal();

      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      return oneSkeleton.buildEdgeLinks(
        edgeList_, edgeStarData_, tetraEdgeList_, edgeLinkData_);
    } else {
      // unsupported dimension
      printErr("Unsupported dimension for edge link precondition");
      return -1;
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionEdgeStarsInternal() {

  if(edgeStarData_.empty()) {
    OneSkeleton oneSkeleton;
    oneSkeleton.setWrapper(this);
    return oneSkeleton.buildEdgeList<3>(
      vertexNumber_, *cellArray_, nullptr, &edgeStarData_, nullptr);
  }
  return 0;
}

int ExplicitTriangulation::preconditionEdgeTrianglesInternal() {

  if(edgeTriangleData_.empty()) {
    preconditionTriangleEdgesInternal();

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);
    return twoSkeleton.buildEdgeTriangles(vertexNumber_, *cellArray_,
                                          edgeTriangleData_, &edgeList_,
                                          &triangleEdgeList_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionTrianglesInternal() {

  if(!triangleList_.size()) {

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_, &triangleList_,
                                  &triangleStarData_, &tetraTriangleList_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionTriangleEdgesInternal() {

  if(!triangleEdgeList_.size()) {

    // WARNING
    // here triangleStarList and cellTriangleList will be computed (for
    // free) although they are not requireed to get the edgeTriangleList.
    // if memory usage is an issue, please change these pointers by nullptr.

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    return twoSkeleton.buildTriangleEdgeList(
      vertexNumber_, *cellArray_, triangleEdgeList_, &vertexEdgeData_,
      &edgeList_, &triangleList_, &triangleStarData_, &tetraTriangleList_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionTriangleLinksInternal() {

  if(triangleLinkData_.empty()) {

    preconditionTriangleStarsInternal();

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);
    return twoSkeleton.buildTriangleLinks(
      triangleList_, triangleStarData_, *cellArray_, triangleLinkData_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionTriangleStarsInternal() {

  if(triangleStarData_.empty()) {

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);
    return twoSkeleton.buildTriangleList(
      vertexNumber_, *cellArray_, &triangleList_, &triangleStarData_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionVertexEdgesInternal() {

  if((SimplexId)vertexEdgeData_.subvectorsNumber() != vertexNumber_) {
    ZeroSkeleton zeroSkeleton;

    if(!edgeList_.size()) {
      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_, &edgeList_);
    }

    zeroSkeleton.setWrapper(this);
    return zeroSkeleton.buildVertexEdges(
      vertexNumber_, edgeList_, vertexEdgeData_);
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexLinksInternal() {

  if((SimplexId)vertexLinkData_.subvectorsNumber() != vertexNumber_) {

    if(getDimensionality() == 2) {
      preconditionEdgesInternal();
      preconditionVertexStarsInternal();

      ZeroSkeleton zeroSkeleton;
      zeroSkeleton.setWrapper(this);
      return zeroSkeleton.buildVertexLinks(
        vertexStarData_, triangleEdgeList_, edgeList_, vertexLinkData_);
    } else if(getDimensionality() == 3) {
      preconditionTrianglesInternal();
      preconditionVertexStarsInternal();

      ZeroSkeleton zeroSkeleton;
      zeroSkeleton.setWrapper(this);
      return zeroSkeleton.buildVertexLinks(
        vertexStarData_, tetraTriangleList_, triangleList_, vertexLinkData_);
    } else {
      // unsupported dimension
      printErr("Unsupported dimension for vertex link precondition");
      return -1;
    }
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexNeighborsInternal() {

  if((SimplexId)vertexNeighborData_.subvectorsNumber() != vertexNumber_) {
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setWrapper(this);
    return zeroSkeleton.buildVertexNeighbors(
      vertexNumber_, *cellArray_, vertexNeighborData_, &edgeList_);
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexStarsInternal() {

  if((SimplexId)vertexStarData_.subvectorsNumber() != vertexNumber_) {
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setWrapper(this);

    return zeroSkeleton.buildVertexStars(
      vertexNumber_, *cellArray_, vertexStarData_);
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexTrianglesInternal() {

  if((SimplexId)vertexTriangleData_.subvectorsNumber() != vertexNumber_) {

    preconditionTrianglesInternal();

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    twoSkeleton.buildVertexTriangles(
      vertexNumber_, triangleList_, vertexTriangleData_);
  }

  return 0;
}
