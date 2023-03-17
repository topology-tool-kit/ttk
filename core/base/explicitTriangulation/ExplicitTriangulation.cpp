#include <ExplicitTriangulation.h>

#include <OneSkeleton.h>
#include <ThreeSkeleton.h>
#include <TwoSkeleton.h>
#include <ZeroSkeleton.h>

#include <cstring>
#include <numeric>

using namespace ttk;

ExplicitTriangulation::ExplicitTriangulation() {

  setDebugMsgPrefix("ExplicitTriangulation");

  clear();
}

ExplicitTriangulation::~ExplicitTriangulation() = default;

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

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((!boundaryEdges_.empty()) && (boundaryEdges_.size() == edgeList_.size())) {
    return 0;
  }

  Timer tm{};

  preconditionEdgesInternal();
  boundaryEdges_.resize(edgeList_.size(), false);

  if(getDimensionality() == 2) {
    preconditionEdgeStarsInternal();
    for(SimplexId i = 0; i < (SimplexId)edgeStarData_.size(); i++) {
      if(edgeStarData_[i].size() == 1) {
        boundaryEdges_[i] = true;
      }
    }
  } else if(getDimensionality() == 3) {
    preconditionTriangleStarsInternal();
    preconditionTriangleEdgesInternal();

    for(size_t i = 0; i < triangleStarData_.size(); i++) {
      if(triangleStarData_[i].size() == 1) {
        for(int j = 0; j < 3; j++) {
          boundaryEdges_[triangleEdgeList_[i][j]] = true;
        }
      }
    }
  } else {
    // unsupported dimension
    printErr("Unsupported dimension for edge boundary precondition");
    return -1;
  }

#ifdef TTK_ENABLE_MPI

  if(ttk::isRunningWithMPI()) {
    ttk::SimplexId edgeNumber = edgeList_.size();
    std::vector<unsigned char> charBoundary(edgeNumber, false);
    for(int i = 0; i < edgeNumber; ++i) {
      charBoundary[i] = boundaryEdges_[i] ? '1' : '0';
    }
    ttk::exchangeGhostDataWithoutTriangulation(
      charBoundary.data(),
      [this](const SimplexId a) { return this->edgeRankArray_[a]; },
      [this](const SimplexId a) { return this->edgeLidToGid_[a]; },
      [this](const SimplexId a) { return this->edgeGidToLid_[a]; }, edgeNumber,
      ttk::MPIcomm_, this->getNeighborRanks());
    for(int i = 0; i < edgeNumber; ++i) {
      boundaryEdges_[i] = (charBoundary[i] == '1');
    }
  }

#endif // TTK_ENABLE_MPI

  this->printMsg("Extracted boundary edges", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

#ifdef TTK_ENABLE_MPI

int ExplicitTriangulation::preconditionVertexRankArray() {
  if(vertexRankArray_.size() == 0) {
    this->vertexRankArray_.resize(this->vertexNumber_, 0);
    if(ttk::isRunningWithMPI()) {
      ttk::produceRankArray(this->vertexRankArray_, this->vertGid_,
                            this->vertexGhost_, this->vertexNumber_,
                            this->boundingBox_.data(), this->neighborRanks_);
      ttk::preconditionNeighborsUsingRankArray<ttk::SimplexId>(
        this->neighborRanks_,
        [this](const ttk::SimplexId a) { return this->getVertexRank(a); },
        this->vertexNumber_, ttk::MPIcomm_);
    }
  }
  return 0;
}

int ExplicitTriangulation::preconditionCellRankArray() {
  if(cellRankArray_.size() == 0) {
    this->cellRankArray_.resize(this->cellNumber_, 0);
    if(ttk::isRunningWithMPI()) {
      ttk::produceRankArray(this->cellRankArray_, this->cellGid_,
                            this->cellGhost_, this->cellNumber_,
                            this->boundingBox_.data(), this->neighborRanks_);
    }
  }
  return 0;
}

int ExplicitTriangulation::preconditionEdgeRankArray() {
  ttk::SimplexId edgeNumber = this->getNumberOfEdgesInternal();
  edgeRankArray_.resize(edgeNumber, 0);
  if(ttk::isRunningWithMPI()) {
    ttk::SimplexId min_id;
    for(ttk::SimplexId id = 0; id < edgeNumber; id++) {
      this->TTK_TRIANGULATION_INTERNAL(getEdgeStar)(id, 0, min_id);
      const auto nStar{this->TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(id)};
      for(SimplexId i = 1; i < nStar; ++i) {
        SimplexId sid{-1};
        this->TTK_TRIANGULATION_INTERNAL(getEdgeStar)(id, i, sid);
        // rule: an edge is owned by the cell in its star with the
        // lowest global id
        if(this->cellGid_[sid] < this->cellGid_[min_id]) {
          min_id = sid;
        }
      }
      edgeRankArray_[id] = cellRankArray_[min_id];
    }
  }
  return 0;
}

int ExplicitTriangulation::preconditionTriangleRankArray() {
  ttk::SimplexId triangleNumber = this->getNumberOfTrianglesInternal();
  triangleRankArray_.resize(triangleNumber, 0);
  if(ttk::isRunningWithMPI()) {
    ttk::SimplexId min_id;
    for(ttk::SimplexId id = 0; id < triangleNumber; id++) {
      this->TTK_TRIANGULATION_INTERNAL(getTriangleStar)(id, 0, min_id);
      const auto nStar{
        this->TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(id)};
      for(SimplexId i = 1; i < nStar; ++i) {
        SimplexId sid{-1};
        this->TTK_TRIANGULATION_INTERNAL(getTriangleStar)(id, i, sid);
        // rule: an triangle is owned by the cell in its star with the
        // lowest global id
        if(this->cellGid_[sid] < this->cellGid_[min_id]) {
          min_id = sid;
        }
      }
      triangleRankArray_[id] = cellRankArray_[min_id];
    }
  }
  return 0;
}

#endif // TTK_ENABLE_MPI

int ExplicitTriangulation::preconditionBoundaryTrianglesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

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

    for(size_t i = 0; i < triangleStarData_.size(); i++) {
      if(triangleStarData_[i].size() == 1) {
        boundaryTriangles_[i] = true;
      }
    }
  } else {
    // unsupported dimension
    printErr("Unsupported dimension for triangle boundary precondition");
    return -1;
  }

#ifdef TTK_ENABLE_MPI

  if(ttk::isRunningWithMPI()) {
    ttk::SimplexId triangleNumber = triangleList_.size();
    std::vector<unsigned char> charBoundary(triangleNumber, false);
    for(int i = 0; i < triangleNumber; ++i) {
      charBoundary[i] = boundaryTriangles_[i] ? '1' : '0';
    }
    ttk::exchangeGhostDataWithoutTriangulation(
      charBoundary.data(),
      [this](const SimplexId a) { return this->triangleRankArray_[a]; },
      [this](const SimplexId a) { return this->triangleLidToGid_[a]; },
      [this](const SimplexId a) { return this->triangleGidToLid_[a]; },
      triangleNumber, ttk::MPIcomm_, this->getNeighborRanks());
    for(int i = 0; i < triangleNumber; ++i) {
      boundaryTriangles_[i] = (charBoundary[i] == '1');
    }
  }

#endif // TTK_ENABLE_MPI

  this->printMsg("Extracted boundary triangles", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionBoundaryVerticesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

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
    for(size_t i = 0; i < vertexStarData_.size(); i++) {
      if(vertexStarData_[i].size() == 1) {
        boundaryVertices_[i] = true;
      }
    }
  } else if(getDimensionality() == 2) {
    preconditionEdgesInternal();
    preconditionEdgeStarsInternal();

    for(SimplexId i = 0; i < (SimplexId)edgeStarData_.size(); i++) {
      if(edgeStarData_[i].size() == 1) {
        boundaryVertices_[edgeList_[i][0]] = true;
        boundaryVertices_[edgeList_[i][1]] = true;
      }
    }
  } else if(getDimensionality() == 3) {
    preconditionTrianglesInternal();
    preconditionTriangleStarsInternal();

    for(size_t i = 0; i < triangleStarData_.size(); i++) {
      if(triangleStarData_[i].size() == 1) {
        boundaryVertices_[triangleList_[i][0]] = true;
        boundaryVertices_[triangleList_[i][1]] = true;
        boundaryVertices_[triangleList_[i][2]] = true;
      }
    }
  } else {
    // unsupported dimension
    printErr("Unsupported dimension for vertex boundary precondition");
    return -1;
  }

#ifdef TTK_ENABLE_MPI

  if(ttk::isRunningWithMPI()) {
    this->preconditionDistributedVertices();
    std::vector<unsigned char> charBoundary(vertexNumber_, false);
    for(int i = 0; i < vertexNumber_; ++i) {
      charBoundary[i] = boundaryVertices_[i] ? '1' : '0';
    }
    ttk::exchangeGhostVertices<unsigned char, ExplicitTriangulation>(
      charBoundary.data(), this, ttk::MPIcomm_);

    for(int i = 0; i < vertexNumber_; ++i) {
      if(vertexRankArray_[i] != ttk::MPIrank_) {
        boundaryVertices_[i] = (charBoundary[i] == '1');
      }
    }
  }

#endif // TTK_ENABLE_MPI

  this->printMsg("Extracted boundary vertices", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionCellEdgesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((tetraEdgeList_.empty() && getDimensionality() == 3)
     || (triangleEdgeList_.empty() && getDimensionality() == 2)) {
    this->preconditionEdgesInternal();
  }

  return 0;
}

int ExplicitTriangulation::preconditionCellNeighborsInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if(cellNeighborData_.empty()) {
    if(getDimensionality() == 3) {
      ThreeSkeleton threeSkeleton;
      threeSkeleton.setWrapper(this);
      threeSkeleton.buildCellNeighborsFromTriangles(
        vertexNumber_, *cellArray_, cellNeighborData_, &triangleStarData_);
    } else if(getDimensionality() == 2) {
      this->preconditionEdgeStarsInternal();
      TwoSkeleton twoSkeleton;
      twoSkeleton.setWrapper(this);
      twoSkeleton.buildCellNeighborsFromEdges(
        *cellArray_, cellNeighborData_, edgeStarData_);
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionCellTrianglesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if(tetraTriangleList_.empty()) {

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

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if(edgeList_.empty()) {
    OneSkeleton oneSkeleton;
    oneSkeleton.setWrapper(this);
    // also computes edgeStar and triangleEdge / tetraEdge lists for free...
    int ret{};
    if(getDimensionality() == 1) {
      std::vector<std::array<SimplexId, 1>> tmp{};
      return oneSkeleton.buildEdgeList<1>(
        vertexNumber_, *cellArray_, edgeList_, edgeStarData_, tmp);
    } else if(getDimensionality() == 2) {
      ret = oneSkeleton.buildEdgeList(vertexNumber_, *cellArray_, edgeList_,
                                      edgeStarData_, triangleEdgeList_);
    } else if(getDimensionality() == 3) {
      ret = oneSkeleton.buildEdgeList(
        vertexNumber_, *cellArray_, edgeList_, edgeStarData_, tetraEdgeList_);
    }

    if(ret != 0) {
      return ret;
    }

#ifdef TTK_ENABLE_MPI
    if(this->getDimensionality() == 2 || this->getDimensionality() == 3) {
      return this->preconditionDistributedEdges();
    }
#endif // TTK_ENABLE_MPI
  }

  return 0;
}

int ExplicitTriangulation::preconditionEdgeLinksInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

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

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if(edgeStarData_.empty()) {
    this->preconditionEdgesInternal();
  }
  return 0;
}

int ExplicitTriangulation::preconditionEdgeTrianglesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if(edgeTriangleData_.empty()) {
    this->preconditionEdgesInternal();
    this->preconditionTriangleEdgesInternal();

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);
    return twoSkeleton.buildEdgeTriangles(vertexNumber_, *cellArray_,
                                          edgeTriangleData_, edgeList_,
                                          &triangleEdgeList_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionTrianglesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if(triangleList_.empty()) {

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    twoSkeleton.buildTriangleList(vertexNumber_, *cellArray_, &triangleList_,
                                  &triangleStarData_, &tetraTriangleList_);

#ifdef TTK_ENABLE_MPI
    this->preconditionDistributedTriangles();
#endif // TTK_ENABLE_MPI
  }

  return 0;
}

int ExplicitTriangulation::preconditionTriangleEdgesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if(triangleEdgeList_.empty()) {
    this->preconditionEdgesInternal();

    // WARNING
    // here triangleStarList and cellTriangleList will be computed (for
    // free) although they are not requireed to get the edgeTriangleList.
    // if memory usage is an issue, please change these pointers by nullptr.

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    return twoSkeleton.buildTriangleEdgeList(
      vertexNumber_, *cellArray_, triangleEdgeList_, edgeList_,
      &vertexEdgeData_, &triangleList_, &triangleStarData_,
      &tetraTriangleList_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionTriangleLinksInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

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

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if(triangleStarData_.empty()) {

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);
    return twoSkeleton.buildTriangleList(
      vertexNumber_, *cellArray_, &triangleList_, &triangleStarData_);
  }

  return 0;
}

int ExplicitTriangulation::preconditionVertexEdgesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((SimplexId)vertexEdgeData_.size() != vertexNumber_) {
    ZeroSkeleton zeroSkeleton;

    if(edgeList_.empty()) {
      this->preconditionEdgesInternal();
    }

    zeroSkeleton.setWrapper(this);
    return zeroSkeleton.buildVertexEdges(
      vertexNumber_, edgeList_, vertexEdgeData_);
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexLinksInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((SimplexId)vertexLinkData_.size() != vertexNumber_) {

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

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((SimplexId)vertexNeighborData_.size() != vertexNumber_) {
    this->preconditionEdgesInternal();
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setWrapper(this);
    return zeroSkeleton.buildVertexNeighbors(
      vertexNumber_, vertexNeighborData_, edgeList_);
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexStarsInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((SimplexId)vertexStarData_.size() != vertexNumber_) {
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setWrapper(this);

    return zeroSkeleton.buildVertexStars(
      vertexNumber_, *cellArray_, vertexStarData_);
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexTrianglesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
#ifdef TTK_ENABLE_MPI
    if(!(ttk::isRunningWithMPI()))
#endif
      this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((SimplexId)vertexTriangleData_.size() != vertexNumber_) {

    preconditionTrianglesInternal();

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    twoSkeleton.buildVertexTriangles(
      vertexNumber_, triangleList_, vertexTriangleData_);
  }

  return 0;
}

int ttk::ExplicitTriangulation::preconditionManifoldInternal() {

  // quick check by numbering (d-1)-simplices star
  FlatJaggedArray *simplexStar{};

  if(this->getDimensionality() == 3) {
    this->preconditionTriangleStarsInternal();
    simplexStar = &this->triangleStarData_;
  } else if(this->getDimensionality() == 2) {
    this->preconditionEdgeStarsInternal();
    simplexStar = &this->edgeStarData_;
  } else if(this->getDimensionality() == 1) {
    this->preconditionVertexStarsInternal();
    simplexStar = &this->vertexStarData_;
  }

  if(simplexStar == nullptr) {
    return 0;
  }

  for(const auto star : *simplexStar) {
    if(star.size() < 1 || star.size() > 2) {
      this->isManifold_ = false;
      this->printWrn("Non manifold data-set detected");
      break;
    }
  }

  return 0;
}

#ifdef TTK_ENABLE_MPI

int ExplicitTriangulation::preconditionDistributedCells() {
  if(this->hasPreconditionedDistributedCells_) {
    return 0;
  }
  if(!ttk::hasInitializedMPI()) {
    return -1;
  }
  if(this->cellGid_ == nullptr) {
    this->printErr("Missing global cell identifiers array!");
    return -2;
  }

  Timer tm{};

  this->preconditionCellRankArray();

  // number of local cells (with ghost cells...)
  const auto nLocCells{this->getNumberOfCells()};

  // global cell id -> local cell id (reverse of this->cellGid_)
  this->cellGidToLid_.reserve(nLocCells);
  for(LongSimplexId lcid = 0; lcid < nLocCells; ++lcid) {
    this->cellGidToLid_[this->cellGid_[lcid]] = lcid;
  }

  this->ghostCellsPerOwner_.resize(ttk::MPIsize_);

  for(LongSimplexId lcid = 0; lcid < nLocCells; ++lcid) {
    if(this->cellRankArray_[lcid] != ttk::MPIrank_) {
      // store ghost cell global ids (per rank)
      this->ghostCellsPerOwner_[this->cellRankArray_[lcid]].emplace_back(
        this->cellGid_[lcid]);
    }
  }

  // for each rank, store the global id of local cells that are ghost cells of
  // other ranks.
  const auto MIT{ttk::getMPIType(ttk::SimplexId{})};
  this->remoteGhostCells_.resize(ttk::MPIsize_);
  // number of owned cells that are ghost cells of other ranks
  std::vector<size_t> nOwnedGhostCellsPerRank(ttk::MPIsize_);

  for(const auto neigh : this->neighborRanks_) {
    // 1. send to neigh number of ghost cells owned by neigh
    const auto nCells{this->ghostCellsPerOwner_[neigh].size()};
    MPI_Sendrecv(&nCells, 1, ttk::getMPIType(nCells), neigh, ttk::MPIrank_,
                 &nOwnedGhostCellsPerRank[neigh], 1, ttk::getMPIType(nCells),
                 neigh, neigh, ttk::MPIcomm_, MPI_STATUS_IGNORE);
    this->remoteGhostCells_[neigh].resize(nOwnedGhostCellsPerRank[neigh]);

    // 2. send to neigh list of ghost cells owned by neigh
    MPI_Sendrecv(this->ghostCellsPerOwner_[neigh].data(),
                 this->ghostCellsPerOwner_[neigh].size(), MIT, neigh,
                 ttk::MPIrank_, this->remoteGhostCells_[neigh].data(),
                 this->remoteGhostCells_[neigh].size(), MIT, neigh, neigh,
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);
  }

  this->preconditionDistributedCellRanges();

  this->hasPreconditionedDistributedCells_ = true;

  if(ttk::MPIrank_ == 0) {
    this->printMsg("Domain contains "
                     + std::to_string(this->gatheredCellRanges_.back().end + 1)
                     + " cells",
                   1.0, tm.getElapsedTime(), this->threadNumber_);
  }

  return 0;
}

int ttk::ExplicitTriangulation::preconditionDistributedCellRanges() {

  // 1. store all local cells owned by current rank by global id

  std::vector<SimplexId> localCellIds{};
  localCellIds.reserve(this->getNumberOfCells());
  for(SimplexId i = 0; i < this->getNumberOfCells(); ++i) {
    if(this->cellRankArray_[i] == ttk::MPIrank_) {
      localCellIds.emplace_back(i);
    }
  }

  TTK_PSORT(this->threadNumber_, localCellIds.begin(), localCellIds.end(),
            [this](const SimplexId a, const SimplexId b) {
              return this->cellGid_[a] < this->cellGid_[b];
            });

  // 2. determine ranges of contiguous cell global ids

  size_t begRange{};
  while(begRange < localCellIds.size()) {
    size_t endRange{begRange + 1};

    if(begRange < localCellIds.size() - 1) {
      for(size_t j = begRange + 1; j < localCellIds.size(); ++j) {
        if(this->cellGid_[localCellIds[j]]
           > this->cellGid_[localCellIds[j - 1]] + 1) {
          endRange = j;
          break;
        }
      }
      if(endRange == begRange + 1
         && this->cellGid_[localCellIds[endRange]]
              == this->cellGid_[localCellIds[endRange - 1]] + 1) {
        endRange = localCellIds.size();
      }
    }

    const size_t gbeg = this->cellGid_[localCellIds[begRange]];
    const size_t gend = this->cellGid_[localCellIds[endRange - 1]];
    const auto nRanges{this->localCellRanges_.size()};

    // inclusive range
    this->localCellRanges_.emplace_back(
      CellRange{nRanges, gbeg, gend, static_cast<size_t>(ttk::MPIrank_)});

    begRange = endRange;
  }

  // 3. send to rank 0 the vector of ranges so it can compute range offsets

  if(ttk::MPIrank_ == 0) {
    this->nRangesPerRank_.resize(ttk::MPIsize_);
  }

  const int rangeSize = this->localCellRanges_.size();
  MPI_Gather(&rangeSize, 1, ttk::getMPIType(rangeSize),
             this->nRangesPerRank_.data(), 1, ttk::getMPIType(rangeSize), 0,
             ttk::MPIcomm_);

  std::vector<int> displacements{};

  if(ttk::MPIrank_ == 0) {
    const auto nRanges{std::accumulate(
      this->nRangesPerRank_.begin(), this->nRangesPerRank_.end(), 0)};
    this->gatheredCellRanges_.resize(nRanges);
    displacements.resize(this->nRangesPerRank_.size());

    for(size_t i = 0; i < this->nRangesPerRank_.size() - 1; ++i) {
      displacements[i + 1] = displacements[i] + this->nRangesPerRank_[i];
    }
  }

  auto cellRangeDT{CellRange::getMPIType()};
  MPI_Type_commit(&cellRangeDT);

  MPI_Gatherv(this->localCellRanges_.data(), this->localCellRanges_.size(),
              cellRangeDT, this->gatheredCellRanges_.data(),
              this->nRangesPerRank_.data(), displacements.data(), cellRangeDT,
              0, ttk::MPIcomm_);

  MPI_Type_free(&cellRangeDT);

  // 4. sort range vector on rank 0

  if(ttk::MPIrank_ == 0) {
    TTK_PSORT(
      this->threadNumber_, this->gatheredCellRanges_.begin(),
      this->gatheredCellRanges_.end(),
      [](const CellRange &a, const CellRange &b) { return a.begin < b.begin; });
  }

  return 0;
}

size_t ttk::ExplicitTriangulation::computeCellRangeOffsets(
  std::vector<size_t> &nSimplicesPerRange) const {

  // 1. send to rank 0 number of edges per cell range

  std::vector<std::vector<size_t>> nSimplicesPerRangePerRank{};

  if(ttk::MPIrank_ == 0) {
    nSimplicesPerRangePerRank.resize(this->nRangesPerRank_.size());
    for(int i = 0; i < ttk::MPIsize_; ++i) {
      nSimplicesPerRangePerRank[i].resize(this->nRangesPerRank_[i]);
    }
    for(int i = 0; i < ttk::MPIsize_; ++i) {
      if(i == 0) {
        continue;
      }
      // receive src content from other ranks
      MPI_Recv(nSimplicesPerRangePerRank[i].data(),
               nSimplicesPerRangePerRank[i].size(), ttk::getMPIType(size_t{}),
               i, MPI_ANY_TAG, ttk::MPIcomm_, MPI_STATUS_IGNORE);
    }
    std::swap(nSimplicesPerRangePerRank[0], nSimplicesPerRange);
  } else {
    MPI_Send(nSimplicesPerRange.data(), nSimplicesPerRange.size(),
             ttk::getMPIType(size_t{}), 0, 0, ttk::MPIcomm_);
  }

  // 2. compute range offsets on rank 0

  size_t nSimplices{};
  if(ttk::MPIrank_ == 0) {

    for(const auto &range : this->gatheredCellRanges_) {
      auto &pSum{nSimplicesPerRangePerRank[range.rank][range.id]};
      std::swap(pSum, nSimplices);
      nSimplices += pSum;
    }
  }

  // 3. send back range offsets to other ranks

  if(ttk::MPIrank_ == 0) {
    for(int i = 1; i < ttk::MPIsize_; ++i) {
      MPI_Send(nSimplicesPerRangePerRank[i].data(),
               nSimplicesPerRangePerRank[i].size(), ttk::getMPIType(size_t{}),
               i, 0, ttk::MPIcomm_);
    }
    std::swap(nSimplicesPerRange, nSimplicesPerRangePerRank[0]);
  } else {
    MPI_Recv(nSimplicesPerRange.data(), nSimplicesPerRange.size(),
             ttk::getMPIType(size_t{}), MPI_ANY_TAG, 0, ttk::MPIcomm_,
             MPI_STATUS_IGNORE);
  }

  return nSimplices;
}

template <typename Func0, typename Func1, typename Func2>
int ttk::ExplicitTriangulation::exchangeDistributedInternal(
  const Func0 &getGlobalSimplexId,
  const Func1 &storeGlobalSimplexId,
  const Func2 &iterCond,
  const int nSimplicesPerCell) {

  // per neighbor, owned ghost cell simplex global ids to transfer back
  std::vector<std::vector<SimplexId>> globalIdPerOwnedGhostCell(ttk::MPIsize_);
  // per neighbor, non-owned ghost cell simplex global ids to transfer back
  std::vector<std::vector<SimplexId>> globalIdPerLocalGhostCell(ttk::MPIsize_);

  const auto MIT{ttk::getMPIType(ttk::SimplexId{})};

  // make sure that all simplices are correctly labelled: for a given
  // rank, a simplex can be owned by a ghost cell from a neighboring
  // rank but in reality can be owned by another ghost cell in a third
  // rank
  bool doIter{true};

  while(doIter) {

    doIter = false;

    // 3. for each list of ghost cell, accumulate the global simplex id
    for(const auto neigh : this->neighborRanks_) {
      // sending side
      globalIdPerOwnedGhostCell[neigh].resize(
        nSimplicesPerCell * this->remoteGhostCells_[neigh].size());
      for(size_t i = 0; i < this->remoteGhostCells_[neigh].size(); ++i) {
        const auto lcid{this->cellGidToLid_[this->remoteGhostCells_[neigh][i]]};
        for(int j = 0; j < nSimplicesPerCell; ++j) {
          globalIdPerOwnedGhostCell[neigh][nSimplicesPerCell * i + j]
            = getGlobalSimplexId(lcid, j);
        }
      }
      // receiving side
      globalIdPerLocalGhostCell[neigh].resize(
        nSimplicesPerCell * this->ghostCellsPerOwner_[neigh].size());

      // 4. transfer back global simplex ids
      MPI_Sendrecv(globalIdPerOwnedGhostCell[neigh].data(),
                   globalIdPerOwnedGhostCell[neigh].size(), MIT, neigh,
                   ttk::MPIrank_, globalIdPerLocalGhostCell[neigh].data(),
                   globalIdPerLocalGhostCell[neigh].size(), MIT, neigh, neigh,
                   ttk::MPIcomm_, MPI_STATUS_IGNORE);
    }

    // 5. extend local <-> global simplex ids mappings
    for(const auto neigh : this->neighborRanks_) {
      for(size_t i = 0; i < this->ghostCellsPerOwner_[neigh].size(); ++i) {
        const auto gcid{this->ghostCellsPerOwner_[neigh][i]};
        const auto lcid{this->cellGidToLid_[gcid]};
        for(int j = 0; j < nSimplicesPerCell; ++j) {
          const auto geid{
            globalIdPerLocalGhostCell[neigh][nSimplicesPerCell * i + j]};
          storeGlobalSimplexId(lcid, geid, j);
        }
      }
    }

    // do an additional transmission if there still is some locally
    // non-labelled simplices
    int doNextIter{0};
    if(iterCond()) {
      doNextIter = 1;
      doIter = true;
    }
    for(int i = 0; i < ttk::MPIsize_; ++i) {
      if(doIter) {
        // reset doNextIter (might have been erased by the MPI_Bcast)
        doNextIter = 1;
      }
      MPI_Bcast(&doNextIter, 1, ttk::getMPIType(doNextIter), i, ttk::MPIcomm_);
      doIter |= (doNextIter == 1);
    }

    if(doIter && ttk::MPIrank_ == 0) {
      this->printMsg("Re-sending global ids to neighbors...");
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionDistributedEdges() {
  if(this->hasPreconditionedDistributedEdges_) {
    return 0;
  }
  if(!ttk::hasInitializedMPI()) {
    return -1;
  }
  if(this->cellGid_ == nullptr) {
    this->printErr("Missing global cell identifiers array!");
    return -2;
  }

  if(this->getDimensionality() != 2 && this->getDimensionality() != 3) {
    return -3;
  }

  if(this->getDimensionality() == 2) {
    this->preconditionTriangleEdges();
  }

  Timer tm{};

  this->preconditionDistributedCells();

  // allocate memory
  this->edgeLidToGid_.resize(this->getNumberOfEdgesInternal(), -1);
  this->edgeGidToLid_.reserve(this->getNumberOfEdgesInternal());

  // 1. for every range of local cells, number the edges locally

  std::vector<SimplexId> edgeLidToRangeId(this->getNumberOfEdgesInternal(), -1);
  std::vector<size_t> nEdgesPerRange(this->localCellRanges_.size());

  const auto edgeAlreadyProcessed = [this](const SimplexId leid,
                                           const SimplexId lcid) {
    const auto nStar{this->TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(leid)};
    for(SimplexId i = 0; i < nStar; ++i) {
      SimplexId sid{-1};
      this->TTK_TRIANGULATION_INTERNAL(getEdgeStar)(leid, i, sid);
      if(sid == -1 || sid == lcid) {
        continue;
      }
      // rule: an edge is owned by the cell in its star with the
      // lowest global id
      if(this->cellGid_[sid] < this->cellGid_[lcid]) {
        return true;
        break;
      }
    }
    return false;
  };

  const auto countCellEdges
    = [this, &edgeAlreadyProcessed](const SimplexId lcid,
                                    std::vector<SimplexId> &edgeGid,
                                    std::vector<SimplexId> &edgeRangeId,
                                    const size_t rangeId, size_t &edgeCount) {
        SimplexId nEdges{};
        if(this->maxCellDim_ == 3) {
          nEdges = this->getCellEdgeNumberInternal(lcid);
        } else if(this->maxCellDim_ == 2) {
          nEdges = this->getTriangleEdgeNumberInternal(lcid);
        }
        for(SimplexId k = 0; k < nEdges; ++k) {
          SimplexId leid{-1};
          if(this->maxCellDim_ == 3) {
            this->getCellEdgeInternal(lcid, k, leid);
          } else if(this->maxCellDim_ == 2) {
            this->getTriangleEdge(lcid, k, leid);
          }
          const auto alreadyProcessed = edgeAlreadyProcessed(leid, lcid);
          if(!alreadyProcessed) {
            edgeGid[leid] = edgeCount;
            edgeRangeId[leid] = rangeId;
            edgeCount++;
          }
        }
      };

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < this->localCellRanges_.size(); ++i) {
    auto &range{this->localCellRanges_[i]};
    range.id = i;
    for(size_t j = range.begin; j <= range.end; ++j) {
      // local cell id
      const auto lcid{this->cellGidToLid_[j]};
      countCellEdges(lcid, this->edgeLidToGid_, edgeLidToRangeId, range.id,
                     nEdgesPerRange[i]);
    }
  }

  // 2. compute range offset on rank 0
  const auto nEdges = this->computeCellRangeOffsets(nEdgesPerRange);

  // 3. locally edit the edge global id with range offsets

  for(SimplexId leid = 0; leid < this->getNumberOfEdgesInternal(); ++leid) {
    if(this->edgeLidToGid_[leid] == -1) {
      // not owned by a cell of this rank
      continue;
    }
    const auto geid{this->edgeLidToGid_[leid]
                    + nEdgesPerRange[edgeLidToRangeId[leid]]};
    this->edgeLidToGid_[leid] = geid;
    this->edgeGidToLid_[geid] = leid;
  }

  // 4. exchange global ids between ghost cells

  const auto nEdgesPerCell{this->getDimensionality() == 3 ? 6 : 3};
  this->exchangeDistributedInternal(
    [this](const SimplexId lcid, const int j) {
      SimplexId leid{};
      this->getCellEdgeInternal(lcid, j, leid);
      return this->edgeLidToGid_[leid];
    },
    [this](const SimplexId lcid, const SimplexId geid, const int j) {
      SimplexId leid{};
      this->getCellEdgeInternal(lcid, j, leid);
      if(this->edgeLidToGid_[leid] == -1 && geid != -1) {
        this->edgeLidToGid_[leid] = geid;
        this->edgeGidToLid_[geid] = leid;
      }
    },
    [this]() {
      return std::count(
               this->edgeLidToGid_.begin(), this->edgeLidToGid_.end(), -1)
             > 0;
    },
    nEdgesPerCell);

  this->preconditionEdgeRankArray();

  if(MPIrank_ == 0) {
    this->printMsg("Domain contains " + std::to_string(nEdges) + " edges", 1.0,
                   tm.getElapsedTime(), 1);
  }

  this->hasPreconditionedDistributedEdges_ = true;

  return 0;
}

int ExplicitTriangulation::preconditionDistributedTriangles() {
  if(this->hasPreconditionedDistributedTriangles_) {
    return 0;
  }
  if(!ttk::hasInitializedMPI()) {
    return -1;
  }
  if(this->cellGid_ == nullptr) {
    this->printErr("Missing global cell identifiers array!");
    return -2;
  }

  if(this->getDimensionality() != 3) {
    return -3;
  }

  Timer tm{};

  this->preconditionDistributedCells();

  // allocate memory
  this->triangleLidToGid_.resize(this->getNumberOfTrianglesInternal(), -1);
  this->triangleGidToLid_.reserve(this->getNumberOfTrianglesInternal());

  // 1. for every range of local cells, number the edges locally

  std::vector<SimplexId> triangleLidToRangeId(
    this->getNumberOfTrianglesInternal(), -1);
  std::vector<size_t> nTrianglesPerRange(this->localCellRanges_.size());

  const auto triangleAlreadyProcessed
    = [this](const SimplexId leid, const SimplexId lcid) {
        const auto nStar{
          this->TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(leid)};
        for(SimplexId i = 0; i < nStar; ++i) {
          SimplexId sid{-1};
          this->TTK_TRIANGULATION_INTERNAL(getTriangleStar)(leid, i, sid);
          if(sid == -1 || sid == lcid) {
            continue;
          }
          // rule: an triangle is owned by the cell in its star with the
          // lowest global id
          if(this->cellGid_[sid] < this->cellGid_[lcid]) {
            return true;
            break;
          }
        }
        return false;
      };

  const auto countCellTriangles
    = [this, &triangleAlreadyProcessed](
        const SimplexId lcid, std::vector<SimplexId> &triangleGid,
        std::vector<SimplexId> &triangleRangeId, const size_t rangeId,
        size_t &triangleCount) {
        const auto nTriangles{this->getCellTriangleNumberInternal(lcid)};
        for(SimplexId k = 0; k < nTriangles; ++k) {
          SimplexId leid{-1};
          this->getCellTriangleInternal(lcid, k, leid);
          const auto alreadyProcessed = triangleAlreadyProcessed(leid, lcid);
          if(!alreadyProcessed) {
            triangleGid[leid] = triangleCount;
            triangleRangeId[leid] = rangeId;
            triangleCount++;
          }
        }
      };

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < this->localCellRanges_.size(); ++i) {
    auto &range{this->localCellRanges_[i]};
    range.id = i;
    for(size_t j = range.begin; j <= range.end; ++j) {
      // local cell id
      const auto lcid{this->cellGidToLid_[j]};
      countCellTriangles(lcid, this->triangleLidToGid_, triangleLidToRangeId,
                         range.id, nTrianglesPerRange[i]);
    }
  }

  // 2. compute range offset on rank 0
  const auto nTriangles = this->computeCellRangeOffsets(nTrianglesPerRange);

  // 3. locally edit the triangle global id with range offsets

  for(SimplexId leid = 0; leid < this->getNumberOfTrianglesInternal(); ++leid) {
    if(this->triangleLidToGid_[leid] == -1) {
      // not owned by a cell of this rank
      continue;
    }
    const auto geid{this->triangleLidToGid_[leid]
                    + nTrianglesPerRange[triangleLidToRangeId[leid]]};
    this->triangleLidToGid_[leid] = geid;
    this->triangleGidToLid_[geid] = leid;
  }

  // 4. exchange global ids between ghost cells

  const auto nTrianglesPerCell{4};
  this->exchangeDistributedInternal(
    [this](const SimplexId lcid, const int j) {
      SimplexId ltid{};
      this->getCellTriangleInternal(lcid, j, ltid);
      return this->triangleLidToGid_[ltid];
    },
    [this](const SimplexId lcid, const SimplexId gtid, const int j) {
      SimplexId ltid{};
      this->getCellTriangleInternal(lcid, j, ltid);
      if(this->triangleLidToGid_[ltid] == -1 && gtid != -1) {
        this->triangleLidToGid_[ltid] = gtid;
        this->triangleGidToLid_[gtid] = ltid;
      }
    },
    [this]() {
      return std::count(this->triangleLidToGid_.begin(),
                        this->triangleLidToGid_.end(), -1)
             > 0;
    },
    nTrianglesPerCell);

  this->preconditionTriangleRankArray();

  if(MPIrank_ == 0) {
    this->printMsg(
      "Domain contains " + std::to_string(nTriangles) + " triangles", 1.0,
      tm.getElapsedTime(), 1);
  }

  this->hasPreconditionedDistributedTriangles_ = true;

  return 0;
}

int ExplicitTriangulation::preconditionDistributedVertices() {
  if(this->hasPreconditionedDistributedVertices_) {
    return 0;
  }
  if(!hasInitializedMPI()) {
    return -1;
  }
  if(this->vertGid_ == nullptr) {
    this->printErr("Missing global vertex identifiers array!");
    return -2;
  }

  this->hasPreconditionedDistributedVertices_ = true;
  this->preconditionVertexRankArray();

  // number of local vertices (with ghost vertices...)
  const auto nLocVertices{this->getNumberOfVertices()};

  // global vertex id -> local vertex id (reverse of this->vertGid_)
  this->vertexGidToLid_.reserve(nLocVertices);
  for(LongSimplexId lvid = 0; lvid < nLocVertices; ++lvid) {
    this->vertexGidToLid_[this->vertGid_[lvid]] = lvid;
  }
  this->ghostVerticesPerOwner_.resize(ttk::MPIsize_);

  for(LongSimplexId lvid = 0; lvid < nLocVertices; ++lvid) {
    if(this->vertexRankArray_[lvid] != ttk::MPIrank_) {
      // store ghost cell global ids (per rank)
      this->ghostVerticesPerOwner_[this->vertexRankArray_[lvid]].emplace_back(
        this->vertGid_[lvid]);
    }
  }

  // for each rank, store the global id of local cells that are ghost cells of
  // other ranks.
  const auto MIT{ttk::getMPIType(ttk::SimplexId{})};
  this->remoteGhostVertices_.resize(ttk::MPIsize_);
  // number of owned cells that are ghost cells of other ranks
  std::vector<size_t> nOwnedGhostVerticesPerRank(ttk::MPIsize_);

  for(const auto neigh : this->neighborRanks_) {
    // 1. send to neigh number of ghost cells owned by neigh
    const auto nVerts{this->ghostVerticesPerOwner_[neigh].size()};
    MPI_Sendrecv(&nVerts, 1, ttk::getMPIType(nVerts), neigh, ttk::MPIrank_,
                 &nOwnedGhostVerticesPerRank[neigh], 1, ttk::getMPIType(nVerts),
                 neigh, neigh, ttk::MPIcomm_, MPI_STATUS_IGNORE);
    this->remoteGhostVertices_[neigh].resize(nOwnedGhostVerticesPerRank[neigh]);

    // 2. send to neigh list of ghost cells owned by neigh
    MPI_Sendrecv(this->ghostVerticesPerOwner_[neigh].data(),
                 this->ghostVerticesPerOwner_[neigh].size(), MIT, neigh,
                 ttk::MPIrank_, this->remoteGhostVertices_[neigh].data(),
                 this->remoteGhostVertices_[neigh].size(), MIT, neigh, neigh,
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);
  }

  return 0;
}

#endif // TTK_ENABLE_MPI

template <typename T>
void writeBin(std::ofstream &stream, const T var) {
  stream.write(reinterpret_cast<const char *>(&var), sizeof(var));
}

template <typename T>
void writeBinArray(std::ofstream &stream,
                   const T *const buff,
                   const size_t size) {
  stream.write(reinterpret_cast<const char *>(buff), size * sizeof(T));
}

// initialize static member variables
const char *ExplicitTriangulation::magicBytes_ = "TTKTriangulationFileFormat";
const unsigned long ExplicitTriangulation::formatVersion_ = 1;

int ExplicitTriangulation::writeToFile(std::ofstream &stream) const {

  // 1. magic bytes (char *)
  stream.write(ttk::ExplicitTriangulation::magicBytes_,
               std::strlen(ttk::ExplicitTriangulation::magicBytes_));
  // 2. format version (unsigned long)
  writeBin(stream, ttk::ExplicitTriangulation::formatVersion_);
  // 3. dimensionality (int)
  const auto dim = this->getDimensionality();
  writeBin(stream, dim);
  // 4. number of vertices (SimplexId)
  const auto nVerts = this->getNumberOfVertices();
  writeBin(stream, nVerts);
  // 5. number of edges (SimplexId)
  const auto edgesNumber = [this, dim]() -> SimplexId {
    if(dim == 1) {
      return this->getNumberOfCells();
    } else if(dim > 1) {
      return this->edgeList_.size();
    }
    return 0;
  };
  const auto nEdges = edgesNumber();
  writeBin(stream, nEdges);
  // 6. number of triangles (SimplexId, 0 in 1D)
  const auto trianglesNumber = [this, dim]() -> SimplexId {
    if(dim == 2) {
      return this->getNumberOfCells();
    } else if(dim == 3) {
      return this->triangleList_.size();
    }
    return 0;
  };
  const auto nTriangles = trianglesNumber();
  writeBin(stream, nTriangles);
  // 7. number of tetrahedron (SimplexId, 0 in 2D)
  const auto nTetras = dim > 2 ? this->getNumberOfCells() : 0;
  writeBin(stream, nTetras);

  // only write buffers oustside this->cellArray_ (cellVertex, vertexCoords),
  // those ones will be provided by VTK

  // fixed-size arrays (in AbstractTriangulation.h)

#define WRITE_FIXED(ARRAY)                             \
  if(ARRAY.empty()) {                                  \
    writeBin(stream, char{0});                         \
  } else {                                             \
    writeBin(stream, char{1});                         \
    writeBinArray(stream, ARRAY.data(), ARRAY.size()); \
  }

  // 8. edgeList (SimplexId array)
  WRITE_FIXED(this->edgeList_);
  // 9. triangleList (SimplexId array)
  WRITE_FIXED(this->triangleList_);
  // 10. triangleEdgeList (SimplexId array)
  WRITE_FIXED(this->triangleEdgeList_);
  // 11. tetraEdgeList (SimplexId array)
  WRITE_FIXED(this->tetraEdgeList_);
  // 12. tetraTriangleList (SimplexId array)
  WRITE_FIXED(this->tetraTriangleList_);

  // variable-size arrays (FlatJaggedArray in ExplicitTriangulation.h)

#define WRITE_GUARD(ARRAY)     \
  if(ARRAY.empty()) {          \
    writeBin(stream, char{0}); \
    return;                    \
  } else {                     \
    writeBin(stream, char{1}); \
  }

  const auto write_variable = [&stream](const FlatJaggedArray &arr) {
    // empty array guard
    WRITE_GUARD(arr);
    writeBinArray(stream, arr.offset_ptr(), arr.size() + 1);
    writeBinArray(stream, arr.get_ptr(0, 0), arr.dataSize());
  };

  // 13. vertexNeighbors (SimplexId array, offsets then data)
  write_variable(this->vertexNeighborData_);
  // 14. cellNeighbors (SimplexId array, offsets then data)
  write_variable(this->cellNeighborData_);
  // 15. vertexEdges (SimplexId array, offsets then data)
  write_variable(this->vertexEdgeData_);
  // 16. vertexTriangles (SimplexId array, offsets then data)
  write_variable(this->vertexTriangleData_);
  // 17. edgeTriangles (SimplexId array, offsets then data)
  write_variable(this->edgeTriangleData_);
  // 18. vertexStars (SimplexId array, offsets then data)
  write_variable(this->vertexStarData_);
  // 19. edgeStars (SimplexId array, offsets then data)
  write_variable(this->edgeStarData_);
  // 20. triangleStars (SimplexId array, offsets then data)
  write_variable(this->triangleStarData_);
  // 21. vertexLinks (SimplexId array, offsets then data)
  write_variable(this->vertexLinkData_);
  // 22. edgeLinks (SimplexId array, offsets then data)
  write_variable(this->edgeLinkData_);
  // 23. triangleLinks (SimplexId array, offsets then data)
  write_variable(this->triangleLinkData_);

  const auto write_bool = [&stream](const std::vector<bool> &arr) {
    // empty array guard
    WRITE_GUARD(arr);
    for(size_t i = 0; i < arr.size(); ++i) {
      const auto b = static_cast<char>(arr[i]);
      writeBin(stream, b);
    }
  };

  // 24. boundary vertices (bool array)
  write_bool(this->boundaryVertices_);
  // 25. boundary edges (bool array)
  write_bool(this->boundaryEdges_);
  // 26. boundary triangles (bool array)
  write_bool(this->boundaryTriangles_);

  return 0;
}

int ExplicitTriangulation::writeToFileASCII(std::ofstream &stream) const {
  // 1. magic bytes
  stream << ttk::ExplicitTriangulation::magicBytes_ << '\n';
  // 2. format version
  stream << ttk::ExplicitTriangulation::formatVersion_ + 1 << '\n';
  // 3. dimensionality
  const auto dim = this->getDimensionality();
  stream << "dim " << dim << '\n';
  // 4. -> 7. number of simplices
  const auto nVerts = this->getNumberOfVertices();
  const auto edgesNumber = [this, dim]() -> SimplexId {
    if(dim == 1) {
      return this->getNumberOfCells();
    } else if(dim > 1) {
      return this->edgeList_.size();
    }
    return 0;
  };
  const auto nEdges = edgesNumber();
  const auto trianglesNumber = [this, dim]() -> SimplexId {
    if(dim == 2) {
      return this->getNumberOfCells();
    } else if(dim == 3) {
      return this->triangleList_.size();
    }
    return 0;
  };
  const auto nTriangles = trianglesNumber();
  const auto nTetras = dim > 2 ? this->getNumberOfCells() : 0;
  stream << nVerts << ' ' << nEdges << ' ' << nTriangles << ' ' << nTetras
         << '\n';

#define WRITE_ASCII(ARRAY)                     \
  stream << #ARRAY << '\n';                    \
  for(const auto &slice : ARRAY) {             \
    for(size_t i = 0; i < slice.size(); ++i) { \
      if(i > 0) {                              \
        stream << ' ';                         \
      }                                        \
      stream << slice[i];                      \
    }                                          \
    stream << '\n';                            \
  }

  const auto writeCellArray = [this, &stream]() {
    stream << "this->cellArray_\n";
    for(SimplexId i = 0; i < this->cellNumber_; ++i) {
      for(SimplexId j = 0; j < this->cellArray_->getCellVertexNumber(i); ++j) {
        stream << this->cellArray_->getCellVertex(i, j) << ' ';
      }
      stream << '\n';
    }
  };

  // ?. cellArray_
  writeCellArray();

  // 8. -> 12. fixed-size arrays
  WRITE_ASCII(this->edgeList_);
  WRITE_ASCII(this->triangleList_);
  WRITE_ASCII(this->triangleEdgeList_);
  WRITE_ASCII(this->tetraEdgeList_);

  // 13. -> 23. variable-size arrays
  WRITE_ASCII(this->vertexNeighborData_);
  WRITE_ASCII(this->cellNeighborData_);
  WRITE_ASCII(this->vertexEdgeData_);
  WRITE_ASCII(this->vertexTriangleData_);
  WRITE_ASCII(this->edgeTriangleData_);
  WRITE_ASCII(this->vertexStarData_);
  WRITE_ASCII(this->edgeStarData_);
  WRITE_ASCII(this->triangleStarData_);
  WRITE_ASCII(this->vertexLinkData_);
  WRITE_ASCII(this->edgeLinkData_);
  WRITE_ASCII(this->triangleLinkData_);

#define WRITE_BOOL(ARRAY)      \
  stream << #ARRAY << '\n';    \
  for(const auto el : ARRAY) { \
    stream << el << '\n';      \
  }

  // 24. -> 26. boolean arrays
  WRITE_BOOL(this->boundaryVertices_);
  WRITE_BOOL(this->boundaryEdges_);
  WRITE_BOOL(this->boundaryTriangles_);

  return 0;
}

template <typename T>
void readBin(std::ifstream &stream, T &res) {
  stream.read(reinterpret_cast<char *>(&res), sizeof(res));
}

template <typename T>
void readBinArray(std::ifstream &stream, T *const res, const size_t size) {
  stream.read(reinterpret_cast<char *>(res), size * sizeof(T));
}

int ExplicitTriangulation::readFromFile(std::ifstream &stream) {

  // 1. magic bytes (char *)
  const auto magicBytesLen
    = std::strlen(ttk::ExplicitTriangulation::magicBytes_);
  std::vector<char> mBytes(magicBytesLen + 1);
  stream.read(mBytes.data(), magicBytesLen);
  const auto hasMagicBytes
    = std::strcmp(mBytes.data(), ttk::ExplicitTriangulation::magicBytes_) == 0;
  if(!hasMagicBytes) {
    this->printErr("Could not find magic bytes in input files!");
    this->printErr("Aborting...");
    return 0;
  }
  // 2. format version (unsigned long)
  unsigned long version{};
  readBin(stream, version);
  if(version != ttk::ExplicitTriangulation::formatVersion_) {
    this->printWrn("File format version (" + std::to_string(version)
                   + ") and software version ("
                   + std::to_string(ttk::ExplicitTriangulation::formatVersion_)
                   + ") are different!");
  }

  int dim{};
  SimplexId nVerts{}, nEdges{}, nTriangles{}, nTetras{};

  // 3. dimensionality (int)
  readBin(stream, dim);
  // 4. number of vertices (SimplexId)
  readBin(stream, nVerts);
  // 5. number of edges (SimplexId)
  readBin(stream, nEdges);
  // 6. number of triangles (SimplexId, 0 in 1D)
  readBin(stream, nTriangles);
  // 7. number of tetrahedron (SimplexId, 0 in 2D)
  readBin(stream, nTetras);

  if(dim != this->getDimensionality()) {
    this->printErr("Incorrect dimension!");
    return 0;
  }
  if(nVerts != this->getNumberOfVertices()) {
    this->printErr("Incorrect number of vertices!");
    return 0;
  }
  if((dim == 2 && nTriangles != this->getNumberOfCells())
     || (dim == 3 && nTetras != this->getNumberOfCells())) {
    this->printErr("Incorrect number of cells!");
    return 0;
  }

  // fixed-size arrays (in AbstractTriangulation.h)

  const auto read_guard = [&stream]() {
    char g{};
    readBin(stream, g);
    return g == 0;
  };

#define READ_FIXED(ARRAY, N_ITEMS)               \
  if(!read_guard()) {                            \
    ARRAY.resize(N_ITEMS);                       \
    readBinArray(stream, ARRAY.data(), N_ITEMS); \
  }

  // 8. edgeList (SimplexId array)
  READ_FIXED(this->edgeList_, nEdges);
  // 9. triangleList (SimplexId array)
  READ_FIXED(this->triangleList_, nTriangles);
  // 10. triangleEdgeList (SimplexId array)
  READ_FIXED(this->triangleEdgeList_, nTriangles);
  // 11. tetraEdgeList (SimplexId array)
  READ_FIXED(this->tetraEdgeList_, nTetras);
  // 12. tetraTriangleList (SimplexId array)
  READ_FIXED(this->tetraTriangleList_, nTetras);

  // variable-size arrays (FlagJaggedArrays in ExplicitTriangulation.h)

  const auto read_variable
    = [&stream, &read_guard](FlatJaggedArray &arr, const SimplexId n_items) {
        // empty array guard
        if(read_guard()) {
          return;
        }
        std::vector<SimplexId> offsets{}, data{};
        offsets.resize(n_items + 1);
        readBinArray(stream, offsets.data(), offsets.size());
        data.resize(offsets.back());
        readBinArray(stream, data.data(), data.size());
        arr.setData(std::move(data), std::move(offsets));
      };

  // 13. vertexNeighbors (SimplexId array, offsets then data)
  read_variable(this->vertexNeighborData_, nVerts);
  // 14. cellNeighbors (SimplexId array, offsets then data)
  read_variable(this->cellNeighborData_, this->getNumberOfCells());
  // 15. vertexEdges (SimplexId array, offsets then data)
  read_variable(this->vertexEdgeData_, nVerts);
  // 16. vertexTriangles (SimplexId array, offsets then data)
  read_variable(this->vertexTriangleData_, nVerts);
  // 17. edgeTriangles (SimplexId array, offsets then data)
  read_variable(this->edgeTriangleData_, nEdges);
  // 18. vertexStars (SimplexId array, offsets then data)
  read_variable(this->vertexStarData_, nVerts);
  // 19. edgeStars (SimplexId array, offsets then data)
  read_variable(this->edgeStarData_, nEdges);
  // 20. triangleStars (SimplexId array, offsets then data)
  read_variable(this->triangleStarData_, nTriangles);
  // 21. vertexLinks (SimplexId array, offsets then data)
  read_variable(this->vertexLinkData_, nVerts);
  // 22. edgeLinks (SimplexId array, offsets then data)
  read_variable(this->edgeLinkData_, nEdges);
  // 23. triangleLinks (SimplexId array, offsets then data)
  read_variable(this->triangleLinkData_, nTriangles);

  const auto read_bool
    = [&stream, &read_guard](std::vector<bool> &arr, const SimplexId n_items) {
        // empty array guard
        if(read_guard()) {
          return;
        }
        // resize vector
        arr.resize(n_items);
        for(SimplexId i = 0; i < n_items; ++i) {
          char b{};
          stream.read(&b, sizeof(b));
          arr[i] = static_cast<bool>(b);
        }
      };

  // 24. boundary vertices (bool array)
  read_bool(this->boundaryVertices_, nVerts);
  // 25. boundary edges (bool array)
  read_bool(this->boundaryEdges_, nEdges);
  // 26. boundary triangles (bool array)
  read_bool(this->boundaryTriangles_, nTriangles);

  return 0;
}
