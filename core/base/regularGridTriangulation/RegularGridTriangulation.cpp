#include <RegularGridTriangulation.h>
#include <numeric>
#include <string>

ttk::RegularGridTriangulation::RegularGridTriangulation()
  : dimensionality_{-1} {
  setDebugMsgPrefix("RegularGridTriangulation");
}

ttk::SimplexId ttk::RegularGridTriangulation::findEdgeFromVertices(
  const SimplexId v0, const SimplexId v1) const {
  // loop over v0 edges to find the one between v0 and v1
  const auto nEdges = this->getVertexEdgeNumberInternal(v0);
  for(SimplexId i = 0; i < nEdges; ++i) {
    SimplexId e{};
    std::array<SimplexId, 2> eVerts{};
    this->getVertexEdgeInternal(v0, i, e);
    this->getEdgeVertexInternal(e, 0, eVerts[0]);
    this->getEdgeVertexInternal(e, 1, eVerts[1]);
    if((v0 == eVerts[0] && v1 == eVerts[1])
       || (v0 == eVerts[1] && v1 == eVerts[0])) {
      return e;
    }
  }

  return -1;
}

ttk::SimplexId ttk::RegularGridTriangulation::findTriangleFromVertices(
  std::array<SimplexId, 3> &verts) const {
  std::sort(verts.begin(), verts.end());

  // loop over verts[0] triangles to find the one shared by all 3
  const auto nTriangles = this->getVertexTriangleNumberInternal(verts[0]);
  for(SimplexId i = 0; i < nTriangles; ++i) {
    SimplexId t{};
    std::array<SimplexId, 3> tVerts{};
    this->getVertexTriangleInternal(verts[0], i, t);
    this->getTriangleVertexInternal(t, 0, tVerts[0]);
    this->getTriangleVertexInternal(t, 1, tVerts[1]);
    this->getTriangleVertexInternal(t, 2, tVerts[2]);
    std::sort(tVerts.begin(), tVerts.end());
    if(tVerts == verts) {
      return t;
    }
  }

  return -1;
}

#ifdef TTK_ENABLE_MPI

int ttk::RegularGridTriangulation::getVertexRankInternal(
  const SimplexId lvid) const {

  if(this->vertexGhost_[lvid] == 0) {
    return ttk::MPIrank_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(this->neighborRanks_.empty()) {
    this->printErr("Empty neighborsRanks_!");
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE
  const auto p{this->getVertGlobalCoords(lvid)};
  for(const auto neigh : this->neighborRanks_) {
    const auto &bbox{this->neighborVertexBBoxes_[neigh]};
    if(p[0] >= bbox[0] && p[0] <= bbox[1] && p[1] >= bbox[2] && p[1] <= bbox[3]
       && p[2] >= bbox[4] && p[2] <= bbox[5]) {
      return neigh;
    }
  }
  return -1;
}

ttk::SimplexId ttk::RegularGridTriangulation::getVertexGlobalIdInternal(
  const SimplexId lvid) const {
  if(!ttk::isRunningWithMPI()) {
    return lvid;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(lvid > this->TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() - 1
     || lvid < 0) {
    return -1;
  }
  if(this->metaGrid_ == nullptr) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  const auto p{this->getVertGlobalCoords(lvid)};
  const auto &dims{this->metaGrid_->getGridDimensions()};

  // global coordinates to identifier (inverse of vertexToPosition)
  return p[0] + p[1] * dims[0] + p[2] * dims[0] * dims[1];
}

ttk::SimplexId ttk::RegularGridTriangulation::getVertexLocalIdInternal(
  const SimplexId gvid) const {
  if(!ttk::isRunningWithMPI()) {
    return gvid;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(this->metaGrid_ == nullptr) {
    return -1;
  }
  if(gvid
       > this->metaGrid_->TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() - 1
     || gvid < 0) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  const auto p{this->getVertLocalCoords(gvid)};
  const auto &dims{this->getGridDimensions()};

  if(p[0] < 0 || p[1] < 0 || p[2] < 0 || p[0] > dims[0] - 1
     || p[1] > dims[1] - 1 || p[2] > dims[2] - 1) {
    return -1;
  }

  // local coordinates to identifier (inverse of vertexToPosition)
  return p[0] + p[1] * dims[0] + p[2] * dims[0] * dims[1];
}

ttk::SimplexId ttk::RegularGridTriangulation::getCellGlobalIdInternal(
  const SimplexId lcid) const {
  if(!ttk::isRunningWithMPI()) {
    return lcid;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(lcid > this->TTK_TRIANGULATION_INTERNAL(getNumberOfCells)() - 1
     || lcid < 0) {
    return -1;
  }
  if(this->metaGrid_ == nullptr) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  // local cube coordinates
  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {
    this->tetrahedronToPosition(lcid, p.data());
  } else if(this->dimensionality_ == 2) {
    this->triangleToPosition2d(lcid, p.data());
  }

  // global cube coordinates
  p[0] += this->localGridOffset_[0];
  p[1] += this->localGridOffset_[1];
  p[2] += this->localGridOffset_[2];

  const auto &dims{this->metaGrid_->getGridDimensions()};

  // global coordinates to identifier (inverse of tetrahedronToPosition)
  const auto globCubeId{p[0] + p[1] * (dims[0] - 1)
                        + p[2] * (dims[0] - 1) * (dims[1] - 1)};

  const auto nCellsPerCube{this->dimensionality_ == 3 ? 6 : 2};
  return globCubeId * nCellsPerCube + lcid % nCellsPerCube;
}

ttk::SimplexId ttk::RegularGridTriangulation::getCellLocalIdInternal(
  const SimplexId gcid) const {
  if(!ttk::isRunningWithMPI()) {
    return gcid;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(this->metaGrid_ == nullptr) {
    return -1;
  }
  if(gcid > this->metaGrid_->TTK_TRIANGULATION_INTERNAL(getNumberOfCells)() - 1
     || gcid < 0) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  // global cube coordinates
  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {
    this->metaGrid_->tetrahedronToPosition(gcid, p.data());
  } else if(this->dimensionality_ == 2) {
    this->metaGrid_->triangleToPosition2d(gcid, p.data());
  }

  // local cube coordinates
  p[0] -= this->localGridOffset_[0];
  p[1] -= this->localGridOffset_[1];
  p[2] -= this->localGridOffset_[2];

  const auto &dims{this->getGridDimensions()};

  // local coordinates to identifier (inverse of tetrahedronToPosition)
  const auto locCubeId{p[0] + p[1] * (dims[0] - 1)
                       + p[2] * (dims[0] - 1) * (dims[1] - 1)};

  const auto nCellsPerCube{this->dimensionality_ == 3 ? 6 : 2};
  return locCubeId * nCellsPerCube + gcid % nCellsPerCube;
}

ttk::SimplexId ttk::RegularGridTriangulation::getEdgeGlobalIdInternal(
  const SimplexId leid) const {
  if(!ttk::isRunningWithMPI()) {
    return leid;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(leid > this->getNumberOfEdgesInternal() - 1 || leid < 0) {
    return -1;
  }
  if(this->metaGrid_ == nullptr) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  if(this->dimensionality_ == 1) {
    return this->getCellGlobalIdInternal(leid);
  }

  // local vertices ids
  SimplexId lv0{}, lv1{};
  this->getEdgeVertexInternal(leid, 0, lv0);
  this->getEdgeVertexInternal(leid, 1, lv1);

  // global vertices ids
  const auto gv0 = this->getVertexGlobalId(lv0);
  const auto gv1 = this->getVertexGlobalId(lv1);
  if(gv0 == -1 || gv1 == -1) {
    return -1;
  }

  return this->metaGrid_->findEdgeFromVertices(gv0, gv1);
}

ttk::SimplexId ttk::RegularGridTriangulation::getEdgeLocalIdInternal(
  const SimplexId geid) const {
  if(!ttk::isRunningWithMPI()) {
    return geid;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(this->metaGrid_ == nullptr) {
    return -1;
  }
  if(geid > this->metaGrid_->getNumberOfEdgesInternal() - 1 || geid < 0) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  if(this->dimensionality_ == 1) {
    return this->getCellLocalIdInternal(geid);
  }

  // global vertices ids
  SimplexId gv0{}, gv1{};
  this->metaGrid_->getEdgeVertexInternal(geid, 0, gv0);
  this->metaGrid_->getEdgeVertexInternal(geid, 1, gv1);

  // local vertices ids
  const auto lv0 = this->getVertexLocalId(gv0);
  const auto lv1 = this->getVertexLocalId(gv1);
  if(lv0 == -1 || lv1 == -1) {
    return -1;
  }

  return this->findEdgeFromVertices(lv0, lv1);
}

ttk::SimplexId ttk::RegularGridTriangulation::getTriangleGlobalIdInternal(
  const SimplexId ltid) const {
  if(!ttk::isRunningWithMPI()) {
    return ltid;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(ltid > this->getNumberOfTrianglesInternal() - 1 || ltid < 0) {
    return -1;
  }
  if(this->metaGrid_ == nullptr) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  if(this->dimensionality_ == 2) {
    return this->getCellGlobalIdInternal(ltid);
  }

  // local vertices ids
  SimplexId lv0{}, lv1{}, lv2{};
  this->getTriangleVertexInternal(ltid, 0, lv0);
  this->getTriangleVertexInternal(ltid, 1, lv1);
  this->getTriangleVertexInternal(ltid, 2, lv2);

  // global vertices ids
  std::array<SimplexId, 3> globVerts{
    this->getVertexGlobalId(lv0),
    this->getVertexGlobalId(lv1),
    this->getVertexGlobalId(lv2),
  };
  for(const auto gv : globVerts) {
    if(gv == -1) {
      return -1;
    }
  }

  return this->metaGrid_->findTriangleFromVertices(globVerts);
}

ttk::SimplexId ttk::RegularGridTriangulation::getTriangleLocalIdInternal(
  const SimplexId gtid) const {
  if(!ttk::isRunningWithMPI()) {
    return gtid;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(this->metaGrid_ == nullptr) {
    return -1;
  }
  if(gtid > this->metaGrid_->getNumberOfTrianglesInternal() - 1 || gtid < 0) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  if(this->dimensionality_ == 2) {
    return this->getCellGlobalIdInternal(gtid);
  }

  // local vertices ids
  SimplexId gv0{}, gv1{}, gv2{};
  this->metaGrid_->getTriangleVertexInternal(gtid, 0, gv0);
  this->metaGrid_->getTriangleVertexInternal(gtid, 1, gv1);
  this->metaGrid_->getTriangleVertexInternal(gtid, 2, gv2);

  // global vertices ids
  std::array<SimplexId, 3> locVerts{
    this->getVertexLocalId(gv0),
    this->getVertexLocalId(gv1),
    this->getVertexLocalId(gv2),
  };
  for(const auto lv : locVerts) {
    if(lv == -1) {
      return -1;
    }
  }

  return this->findTriangleFromVertices(locVerts);
}

int ttk::RegularGridTriangulation::preconditionDistributedVertices() {
  if(this->hasPreconditionedDistributedVertices_) {
    return 0;
  }
  if(!isRunningWithMPI()) {
    return -1;
  }
  if(this->vertexGhost_ == nullptr) {
    if(ttk::isRunningWithMPI()) {
      this->printErr("Missing vertex ghost array!");
    }
    return -3;
  }
  // number of local vertices (with ghost vertices...)
  const auto nLocVertices{this->getNumberOfVertices()};
  this->neighborVertexBBoxes_.resize(ttk::MPIsize_);
  auto &localBBox{this->neighborVertexBBoxes_[ttk::MPIrank_]};
  // "good" starting values?
  ttk::SimplexId localBBox_x_min{this->localGridOffset_[0]
                                 + this->dimensions_[0]},
    localBBox_y_min{this->localGridOffset_[1] + this->dimensions_[1]},
    localBBox_z_min{this->localGridOffset_[2] + this->dimensions_[2]};
  ttk::SimplexId localBBox_x_max{this->localGridOffset_[0]},
    localBBox_y_max{this->localGridOffset_[1]},
    localBBox_z_max{this->localGridOffset_[2]};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for reduction(                    \
  min                                                  \
  : localBBox_x_min, localBBox_y_min, localBBox_z_min) \
  reduction(max                                        \
            : localBBox_x_max, localBBox_y_max, localBBox_z_max)
#endif
  for(SimplexId lvid = 0; lvid < nLocVertices; ++lvid) {
    // only keep non-ghost vertices
    if(this->vertexGhost_[lvid] != 0) {
      continue;
    }
    // local vertex coordinates
    std::array<SimplexId, 3> p{};
    p = this->getVertGlobalCoords(lvid);

    if(p[0] < localBBox_x_min) {
      localBBox_x_min = p[0];
    }
    if(p[0] > localBBox_x_max) {
      localBBox_x_max = p[0];
    }
    if(p[1] < localBBox_y_min) {
      localBBox_y_min = p[1];
    }
    if(p[1] > localBBox_y_max) {
      localBBox_y_max = p[1];
    }
    if(p[2] < localBBox_z_min) {
      localBBox_z_min = p[2];
    }
    if(p[2] > localBBox_z_max) {
      localBBox_z_max = p[2];
    }
  }

  localBBox = {
    localBBox_x_min, localBBox_x_max, localBBox_y_min,
    localBBox_y_max, localBBox_z_min, localBBox_z_max,
  };

  for(size_t i = 0; i < this->neighborRanks_.size(); ++i) {
    const auto neigh{this->neighborRanks_[i]};
    MPI_Sendrecv(this->neighborVertexBBoxes_[ttk::MPIrank_].data(), 6,
                 ttk::getMPIType(SimplexId{}), neigh, ttk::MPIrank_,
                 this->neighborVertexBBoxes_[neigh].data(), 6,
                 ttk::getMPIType(SimplexId{}), neigh, neigh, ttk::MPIcomm_,
                 MPI_STATUS_IGNORE);
  }

  this->hasPreconditionedDistributedVertices_ = true;

  return 0;
}

int ttk::RegularGridTriangulation::preconditionExchangeGhostCells() {

  if(this->hasPreconditionedExchangeGhostCells_) {
    return 0;
  }
  if(!ttk::hasInitializedMPI()) {
    return -1;
  }

  // number of local cells (with ghost cells...)
  const auto nLocCells{this->getNumberOfCells()};

  int cellRank = 0;
  for(LongSimplexId lcid = 0; lcid < nLocCells; ++lcid) {
    cellRank = this->getCellRankInternal(lcid);
    if(cellRank != ttk::MPIrank_) {
      // store ghost cell global ids (per rank)
      this->ghostCellsPerOwner_[cellRank].emplace_back(
        this->getCellGlobalIdInternal(lcid));
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

  this->hasPreconditionedExchangeGhostCells_ = true;
  return 0;
}

int ttk::RegularGridTriangulation::preconditionExchangeGhostVertices() {

  if(this->hasPreconditionedExchangeGhostVertices_) {
    return 0;
  }
  if((!ttk::hasInitializedMPI()) || (!ttk::isRunningWithMPI())) {
    return -1;
  }

  this->ghostVerticesPerOwner_.resize(ttk::MPIsize_);

  // number of local vertices (with ghost vertices...)
  const auto nLocVertices{this->getNumberOfVertices()};

  for(LongSimplexId lvid = 0; lvid < nLocVertices; ++lvid) {
    if(this->getVertexRankInternal(lvid) != ttk::MPIrank_) {
      // store ghost cell global ids (per rank)
      this->ghostVerticesPerOwner_[this->getVertexRankInternal(lvid)]
        .emplace_back(this->getVertexGlobalIdInternal(lvid));
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

  this->hasPreconditionedExchangeGhostVertices_ = true;
  return 0;
}
#endif
