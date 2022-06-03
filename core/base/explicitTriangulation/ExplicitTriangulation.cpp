#include <ExplicitTriangulation.h>
#include <OneSkeleton.h>
#include <ThreeSkeleton.h>
#include <TwoSkeleton.h>
#include <ZeroSkeleton.h>
#include <cstring>

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
    printErr("Unsupported dimension for edge boundary precondition");
    return -1;
  }

  this->printMsg("Extracted boundary edges", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionBoundaryTrianglesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
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

    for(size_t i = 0; i < triangleStarData_.subvectorsNumber(); i++) {
      if(triangleStarData_.size(i) == 1) {
        boundaryTriangles_[i] = true;
      }
    }
  } else {
    // unsupported dimension
    printErr("Unsupported dimension for triangle boundary precondition");
    return -1;
  }

  this->printMsg("Extracted boundary triangles", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionBoundaryVerticesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
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
    printErr("Unsupported dimension for vertex boundary precondition");
    return -1;
  }

  this->printMsg("Extracted boundary vertices", 1.0, tm.getElapsedTime(), 1);

  return 0;
}

int ExplicitTriangulation::preconditionCellEdgesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
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
    this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((SimplexId)vertexEdgeData_.subvectorsNumber() != vertexNumber_) {
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
    this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

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

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
    this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((SimplexId)vertexNeighborData_.subvectorsNumber() != vertexNumber_) {
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
    this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((SimplexId)vertexStarData_.subvectorsNumber() != vertexNumber_) {
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setWrapper(this);

    return zeroSkeleton.buildVertexStars(
      vertexNumber_, *cellArray_, vertexStarData_);
  }
  return 0;
}

int ExplicitTriangulation::preconditionVertexTrianglesInternal() {

  if(this->cellArray_ == nullptr || this->vertexNumber_ == 0) {
    this->printErr("Empty dataset, precondition skipped");
    return 1;
  }

  if((SimplexId)vertexTriangleData_.subvectorsNumber() != vertexNumber_) {

    preconditionTrianglesInternal();

    TwoSkeleton twoSkeleton;
    twoSkeleton.setWrapper(this);

    twoSkeleton.buildVertexTriangles(
      vertexNumber_, triangleList_, vertexTriangleData_);
  }

  return 0;
}

#ifdef TTK_ENABLE_MPI

int ExplicitTriangulation::preconditionDistributedCells() {
  if(this->hasPreconditionedDistributedCells_) {
    return 0;
  }
  if(!ttk::isRunningWithMPI()) {
    return -1;
  }
  if(this->cellGid_ == nullptr) {
    this->printWrn("Missing global identifiers on cells");
    return -2;
  }
  if(this->ghostCellMask_ == nullptr) {
    this->printWrn("Missing ghost mask on cells");
    return -3;
  }

  // number of local cells (with ghost cells...)
  const auto nLocCells{this->getNumberOfCells()};

  // global cell id -> local cell id (reverse of this->cellGid_)
  this->cellGidToLid_.reserve(nLocCells);
  for(LongSimplexId lcid = 0; lcid < nLocCells; ++lcid) {
    this->cellGidToLid_[this->cellGid_[lcid]] = lcid;
  }

  // global cell ids for cell owned by current rank
  std::vector<LongSimplexId> globalCellsInRank{};
  globalCellsInRank.reserve(nLocCells);

  for(LongSimplexId lcid = 0; lcid < nLocCells; ++lcid) {
    if(this->ghostCellMask_[lcid] == 0) {
      globalCellsInRank.emplace_back(this->cellGid_[lcid]);
    }
  }

  // global cell ids per rank (used by rank 0)
  std::vector<std::vector<LongSimplexId>> cellsRankPerProc{};

  gatherVectors(cellsRankPerProc, globalCellsInRank, 0);

  if(ttk::MPIrank_ == 0) {
    // global number of cells
    auto nGlobalCells = cellsRankPerProc[0].size();
    for(int i = 1; i < MPIsize_; ++i) {
      nGlobalCells += cellsRankPerProc[i].size();
    }

    this->cellGidToRank_.resize(nGlobalCells, -1);
    for(int i = 0; i < MPIsize_; ++i) {
      for(const auto gcid : cellsRankPerProc[i]) {
        this->cellGidToRank_[gcid] = i;
      }
    }
  }

  this->hasPreconditionedDistributedCells_ = true;

  return 0;
}

template <typename Func>
int preconditionDistributedIntermediate(size_t &globalCount,
                                        const std::vector<int> &cellGidToRank,
                                        const Func &processCells) {

  if(!ttk::isRunningWithMPI()) {
    return -1;
  }

  const int COMPUTE = 1, STOP = 2;

  if(ttk::MPIrank_ == 0) {
    size_t gcid{};
    while(gcid < cellGidToRank.size()) {
      const auto currRank{cellGidToRank[gcid]};
      size_t endCurrRank{gcid + 1};
      if(gcid < cellGidToRank.size() - 1) {
        for(size_t j = gcid + 1; j < cellGidToRank.size(); ++j) {
          if(cellGidToRank[j] != currRank) {
            endCurrRank = j;
            break;
          }
        }
        if(endCurrRank == gcid + 1 && cellGidToRank[endCurrRank] == currRank) {
          endCurrRank = cellGidToRank.size();
        }
      }

      if(currRank == 0) {
        // process local edges
        processCells(gcid, endCurrRank, globalCount);
      } else {
        // send to relevant process
        std::array<unsigned long, 3> toSend{globalCount, gcid, endCurrRank};
        MPI_Send(toSend.data(), 3, MPI_UNSIGNED_LONG, currRank, COMPUTE,
                 MPI_COMM_WORLD);
        // receive updated edgeCount
        MPI_Recv(&globalCount, 1, MPI_UNSIGNED_LONG, currRank, MPI_ANY_TAG,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      gcid = endCurrRank;
    }

    // send STOP signal to break infinite loop
    std::array<unsigned long, 3> dummy{};
    for(int i = 1; i < ttk::MPIsize_; ++i) {
      MPI_Send(dummy.data(), 3, MPI_UNSIGNED_LONG, i, STOP, MPI_COMM_WORLD);
    }

  } else {
    while(true) {
      // receive data
      std::array<unsigned long, 3> toRecv{};
      MPI_Status status;
      MPI_Recv(toRecv.data(), 3, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG,
               MPI_COMM_WORLD, &status);
      if(status.MPI_TAG == STOP) {
        break;
      }
      // process local edges
      globalCount = toRecv[0];
      processCells(toRecv[1], toRecv[2], globalCount);
      // send back updated edgeCount to rank 0
      MPI_Send(&globalCount, 1, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
    }
  }

  return 0;
}

int ExplicitTriangulation::preconditionDistributedEdges() {
  if(this->hasPreconditionedDistributedEdges_) {
    return 0;
  }
  if(!ttk::isRunningWithMPI()) {
    return -1;
  }
  if(this->cellGid_ == nullptr) {
    this->printWrn("Missing global identifiers on cells");
    return -2;
  }

  if(this->getDimensionality() != 2 && this->getDimensionality() != 3) {
    return -3;
  }

  this->preconditionDistributedCells();

  // allocate memory
  this->edgeLidToGid_.resize(this->getNumberOfEdgesInternal(), -1);
  this->edgeGidToLid_.reserve(this->getNumberOfEdgesInternal());

  const auto processEdge = [this](const SimplexId leid, const SimplexId lcid,
                                  size_t &edgeCount) {
    bool alreadyProcessed{false};
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
        alreadyProcessed = true;
        break;
      }
    }
    if(!alreadyProcessed) {
      this->edgeLidToGid_[leid] = edgeCount;
      this->edgeGidToLid_[edgeCount] = leid;
      edgeCount++;
    }
  };

  const auto processCellsEdges
    = [this, &processEdge](
        const size_t gcid, const size_t endCurrRank, size_t &edgeCount) {
        for(size_t j = gcid; j < endCurrRank; ++j) {
          // local cell id
          const auto lcid{this->cellGidToLid_[j]};
          for(SimplexId i = 0; i < this->getCellEdgeNumberInternal(lcid); ++i) {
            SimplexId leid{-1};
            this->getCellEdgeInternal(lcid, i, leid);
            processEdge(leid, lcid, edgeCount);
          }
        }
      };

  size_t edgeCount{};
  preconditionDistributedIntermediate(
    edgeCount, this->cellGidToRank_, processCellsEdges);

  if(MPIrank_ == 0) {
    this->printMsg("Domain contains " + std::to_string(edgeCount) + " edges");
  }

  this->hasPreconditionedDistributedEdges_ = true;

  return 0;
}

int ExplicitTriangulation::preconditionDistributedTriangles() {
  if(this->hasPreconditionedDistributedTriangles_) {
    return 0;
  }
  if(!ttk::isRunningWithMPI()) {
    return -1;
  }
  if(this->cellGid_ == nullptr) {
    this->printWrn("Missing global identifiers on cells");
    return -2;
  }

  if(this->getDimensionality() != 3) {
    return -3;
  }

  this->preconditionDistributedCells();

  // allocate memory
  this->triangleLidToGid_.resize(this->getNumberOfTrianglesInternal(), -1);
  this->triangleGidToLid_.reserve(this->getNumberOfTrianglesInternal());

  const auto processTriangle =
    [this](const SimplexId ltid, const SimplexId lcid, size_t &triangleCount) {
      bool alreadyProcessed{false};
      const auto nStar{
        this->TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(ltid)};
      for(SimplexId i = 0; i < nStar; ++i) {
        SimplexId sid{-1};
        this->TTK_TRIANGULATION_INTERNAL(getTriangleStar)(ltid, i, sid);
        if(sid == -1 || sid == lcid) {
          continue;
        }
        // rule: an edge is owned by the cell in its star with the
        // lowest global id
        if(this->cellGid_[sid] < this->cellGid_[lcid]) {
          alreadyProcessed = true;
          break;
        }
      }
      if(!alreadyProcessed) {
        this->triangleLidToGid_[ltid] = triangleCount;
        this->triangleGidToLid_[triangleCount] = ltid;
        triangleCount++;
      }
    };

  const auto processCellsTriangles
    = [this, &processTriangle](
        const size_t gcid, const size_t endCurrRank, size_t &triangleCount) {
        for(size_t j = gcid; j < endCurrRank; ++j) {
          // local cell id
          const auto lcid{this->cellGidToLid_[j]};
          for(SimplexId i = 0; i < this->getCellTriangleNumberInternal(lcid);
              ++i) {
            SimplexId ltid{-1};
            this->getCellTriangleInternal(lcid, i, ltid);
            processTriangle(ltid, lcid, triangleCount);
          }
        }
      };

  size_t triangleCount{};
  preconditionDistributedIntermediate(
    triangleCount, this->cellGidToRank_, processCellsTriangles);

  if(MPIrank_ == 0) {
    this->printMsg("Domain contains " + std::to_string(triangleCount)
                   + " triangles");
  }

  this->hasPreconditionedDistributedTriangles_ = true;

  return 0;
}

int ExplicitTriangulation::preconditionDistributedVertices() {
  if(this->hasPreconditionedDistributedVertices_) {
    return 0;
  }
  if(!isRunningWithMPI()) {
    return -1;
  }
  if(this->globalIdsArray_ == nullptr) {
    this->printWrn("Missing global identifiers array!");
    return -2;
  }

  // allocate memory
  this->vertexLidToGid_.resize(this->vertexNumber_, -1);
  this->vertexGidToLid_.reserve(this->vertexNumber_);

  for(SimplexId i = 0; i < this->vertexNumber_; ++i) {
    this->vertexLidToGid_[i] = this->globalIdsArray_[i];
    this->vertexGidToLid_[this->globalIdsArray_[i]] = i;
  }

  if(MPIrank_ == 0) {
    this->printMsg("Domain contains "
                   + std::to_string(this->getNumberOfVerticesInternal())
                   + " vertices");
  }

  this->hasPreconditionedDistributedVertices_ = true;

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
  stream.write(this->magicBytes_, std::strlen(this->magicBytes_));
  // 2. format version (unsigned long)
  writeBin(stream, this->formatVersion_);
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
    writeBinArray(stream, arr.offset_ptr(), arr.subvectorsNumber() + 1);
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
  const auto magicBytesLen = std::strlen(this->magicBytes_);
  std::vector<char> mBytes(magicBytesLen + 1);
  stream.read(mBytes.data(), magicBytesLen);
  const auto hasMagicBytes = std::strcmp(mBytes.data(), this->magicBytes_) == 0;
  if(!hasMagicBytes) {
    this->printErr("Could not find magic bytes in input files!");
    this->printErr("Aborting...");
    return 0;
  }
  // 2. format version (unsigned long)
  unsigned long version{};
  readBin(stream, version);
  if(version != this->formatVersion_) {
    this->printWrn("File format version (" + std::to_string(version)
                   + ") and software version ("
                   + std::to_string(this->formatVersion_) + ") are different!");
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
    if(g == 0) {
      return true;
    }
    return false;
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
