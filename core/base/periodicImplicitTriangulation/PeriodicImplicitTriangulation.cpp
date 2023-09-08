#include <PeriodicImplicitTriangulation.h>
using namespace std;
using namespace ttk;

PeriodicImplicitTriangulation::PeriodicImplicitTriangulation()
  : dimensionality_{-1}, cellNumber_{}, vertexNumber_{}, edgeNumber_{},
    triangleNumber_{}, tetrahedronNumber_{}, isAccelerated_{} {
  setDebugMsgPrefix("PeriodicImplicitTriangulation");
  hasPeriodicBoundaries_ = true;
}

PeriodicImplicitTriangulation::~PeriodicImplicitTriangulation() = default;

int PeriodicImplicitTriangulation::setInputGrid(const float &xOrigin,
                                                const float &yOrigin,
                                                const float &zOrigin,
                                                const float &xSpacing,
                                                const float &ySpacing,
                                                const float &zSpacing,
                                                const SimplexId &xDim,
                                                const SimplexId &yDim,
                                                const SimplexId &zDim) {

  // Dimensionality //
  if(xDim < 1 or yDim < 1 or zDim < 1)
    dimensionality_ = -1;
  else if(xDim > 1 and yDim > 1 and zDim > 1)
    dimensionality_ = 3;
  else if((xDim > 1 and yDim > 1) or (yDim > 1 and zDim > 1)
          or (xDim > 1 and zDim > 1))
    dimensionality_ = 2;
  else if(xDim > 1 or yDim > 1 or zDim > 1)
    dimensionality_ = 1;
  else
    dimensionality_ = 0;

  // Essentials //
  origin_[0] = xOrigin;
  origin_[1] = yOrigin;
  origin_[2] = zOrigin;
  spacing_[0] = xSpacing;
  spacing_[1] = ySpacing;
  spacing_[2] = zSpacing;
  dimensions_[0] = xDim;
  dimensions_[1] = yDim;
  dimensions_[2] = zDim;
  nbvoxels_[0] = xDim - 1;
  nbvoxels_[1] = yDim - 1;
  nbvoxels_[2] = zDim - 1;

  if(dimensionality_ == 3) {
    // VertexShift
    vshift_[0] = xDim;
    vshift_[1] = xDim * yDim;
    // EdgeSetDimensions
    esetdims_[0] = xDim * yDim * zDim;
    esetdims_[1] = xDim * yDim * zDim;
    esetdims_[2] = xDim * yDim * zDim;
    esetdims_[3] = xDim * yDim * zDim;
    esetdims_[4] = xDim * yDim * zDim;
    esetdims_[5] = xDim * yDim * zDim;
    esetdims_[6] = xDim * yDim * zDim;
    // EdgeSetShift
    esetshift_[0] = esetdims_[0];
    for(int k = 1; k < 7; ++k)
      esetshift_[k] = esetshift_[k - 1] + esetdims_[k];
    // EdgeShift
    eshift_[0] = xDim;
    eshift_[1] = xDim * yDim;
    eshift_[2] = xDim;
    eshift_[3] = xDim * yDim;
    eshift_[4] = xDim;
    eshift_[5] = xDim * yDim;
    eshift_[6] = xDim;
    eshift_[7] = xDim * yDim;
    eshift_[8] = xDim;
    eshift_[9] = xDim * yDim;
    eshift_[10] = xDim;
    eshift_[11] = xDim * yDim;
    eshift_[12] = xDim;
    eshift_[13] = xDim * yDim;
    // TriangleSetDimensions
    tsetdims_[0] = xDim * yDim * zDim * 2;
    tsetdims_[1] = xDim * yDim * zDim * 2;
    tsetdims_[2] = xDim * yDim * zDim * 2;
    tsetdims_[3] = xDim * yDim * zDim * 2;
    tsetdims_[4] = xDim * yDim * zDim * 2;
    tsetdims_[5] = xDim * yDim * zDim * 2;
    // TriangleSetShift
    tsetshift_[0] = tsetdims_[0];
    for(int k = 1; k < 6; ++k)
      tsetshift_[k] = tsetshift_[k - 1] + tsetdims_[k];
    // TriangleShift
    tshift_[0] = xDim * 2;
    tshift_[1] = xDim * yDim * 2;
    tshift_[2] = xDim * 2;
    tshift_[3] = xDim * yDim * 2;
    tshift_[4] = xDim * 2;
    tshift_[5] = xDim * yDim * 2;
    tshift_[6] = xDim * 2;
    tshift_[7] = xDim * yDim * 2;
    tshift_[8] = xDim * 2;
    tshift_[9] = xDim * yDim * 2;
    tshift_[10] = xDim * 2;
    tshift_[11] = xDim * yDim * 2;
    // TetrahedronShift
    tetshift_[0] = xDim * 6;
    tetshift_[1] = xDim * yDim * 6;

    wrap_[0] = xDim;
    wrap_[1] = xDim * yDim;
    wrap_[2] = xDim * yDim * zDim;

    // Numbers
    vertexNumber_ = xDim * yDim * zDim;
    edgeNumber_ = 0;
    for(int k = 0; k < 7; ++k)
      edgeNumber_ += esetdims_[k];
    triangleNumber_ = 0;
    for(int k = 0; k < 6; ++k)
      triangleNumber_ += tsetdims_[k];
    tetrahedronNumber_ = xDim * yDim * zDim * 6;
    cellNumber_ = tetrahedronNumber_;

    checkAcceleration();
  } else if(dimensionality_ == 2) {
    // dimensions selectors
    if(xDim == 1) {
      Di_ = 1;
      Dj_ = 2;
    } else if(yDim == 1) {
      Di_ = 0;
      Dj_ = 2;
    } else {
      Di_ = 0;
      Dj_ = 1;
    }
    // VertexShift
    vshift_[0] = dimensions_[Di_];
    // EdgeSetDimensions
    esetdims_[0] = dimensions_[Di_] * dimensions_[Dj_];
    esetdims_[1] = dimensions_[Di_] * dimensions_[Dj_];
    esetdims_[2] = dimensions_[Di_] * dimensions_[Dj_];
    // EdgeSetShift
    esetshift_[0] = esetdims_[0];
    for(int k = 1; k < 3; ++k)
      esetshift_[k] = esetshift_[k - 1] + esetdims_[k];
    // EdgeShift
    eshift_[0] = dimensions_[Di_];
    eshift_[2] = dimensions_[Di_];
    eshift_[4] = dimensions_[Di_];
    // TriangleShift
    tshift_[0] = dimensions_[Di_] * 2;

    wrap_[0] = dimensions_[Di_];
    wrap_[1] = dimensions_[Di_] * dimensions_[Dj_];

    // Numbers
    vertexNumber_ = dimensions_[Di_] * dimensions_[Dj_];
    edgeNumber_ = 0;
    for(int k = 0; k < 3; ++k)
      edgeNumber_ += esetdims_[k];
    triangleNumber_ = dimensions_[Di_] * dimensions_[Dj_] * 2;
    cellNumber_ = triangleNumber_;

    checkAcceleration();
  } else if(dimensionality_ == 1) {
    // dimensions selectors
    for(int k = 0; k < 3; ++k) {
      if(dimensions_[k] > 1) {
        Di_ = k;
        break;
      }
    }

    // Numbers
    vertexNumber_ = dimensions_[Di_];
    edgeNumber_ = vertexNumber_;
    cellNumber_ = edgeNumber_;
  }

  return 0;
}

int PeriodicImplicitTriangulation::checkAcceleration() {
  isAccelerated_ = false;

  unsigned long long int msb[3];
  if(dimensionality_ == 3) {
    bool allDimensionsArePowerOfTwo = true;
    for(int k = 0; k < 3; ++k)
      if(!isPowerOfTwo(dimensions_[k], msb[k]))
        allDimensionsArePowerOfTwo = false;

    if(allDimensionsArePowerOfTwo) {
      mod_[0] = dimensions_[0] - 1;
      mod_[1] = dimensions_[0] * dimensions_[1] - 1;
      div_[0] = msb[0];
      div_[1] = msb[0] + msb[1];
      isAccelerated_ = true;
    }
  } else if(dimensionality_ == 2) {
    bool const isDi = isPowerOfTwo(dimensions_[Di_], msb[Di_]);
    bool const isDj = isPowerOfTwo(dimensions_[Dj_], msb[Dj_]);
    bool const allDimensionsArePowerOfTwo = (isDi and isDj);

    if(allDimensionsArePowerOfTwo) {
      mod_[0] = dimensions_[Di_] - 1;
      div_[0] = msb[Di_];
      isAccelerated_ = true;
    }
  }

  if(isAccelerated_) {
    printMsg("Accelerated getVertex*() requests.", debug::Priority::INFO);
  }

  return 0;
}

bool PeriodicImplicitTriangulation::isPowerOfTwo(unsigned long long int v,
                                                 unsigned long long int &r) {
  if(v && !(v & (v - 1))) {
    r = 0;
    while(v >>= 1)
      r++;
    return true;
  }
  return false;
}

bool PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  isVertexOnBoundary)(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#else
  TTK_FORCE_USE(vertexId);
#endif

  return false;
}

bool PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  isEdgeOnBoundary)(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#else
  TTK_FORCE_USE(edgeId);
#endif

  return false;
}

bool PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  isTriangleOnBoundary)(const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#else
  TTK_FORCE_USE(triangleId);
#endif

  return false;
}

SimplexId PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getVertexNeighborNumber)(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#else
  TTK_FORCE_USE(vertexId);
#endif

  if(dimensionality_ == 3) {
    return 14; // abcdefgh
  } else if(dimensionality_ == 2) {
    return 6; // abcd
  } else if(dimensionality_ == 1) {
    return 2; // ab
  }

  return -1;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexNeighbor)(const SimplexId &vertexId,
                     const int &localNeighborId,
                     SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= getVertexNeighborNumber(vertexId))
    return -1;
#endif

  neighborId = -1;
  const auto &p = this->underlying().getVertexCoords(vertexId);

  if(dimensionality_ == 3) {
    neighborId = getVertexNeighbor3d(p.data(), vertexId, localNeighborId);
  } else if(dimensionality_ == 2) {
    neighborId = getVertexNeighbor2d(p.data(), vertexId, localNeighborId);
  } else if(dimensionality_ == 1) {
    // ab
    if(vertexId > 0 and vertexId < nbvoxels_[Di_]) {
      if(localNeighborId == 0)
        neighborId = vertexId + 1;
      else
        neighborId = vertexId - 1;
    } else if(vertexId == 0) {
      if(localNeighborId == 0)
        neighborId = vertexId + 1; // a
      else
        neighborId = nbvoxels_[Di_];
    } else {
      if(localNeighborId == 0)
        neighborId = 0; // a
      else
        neighborId = vertexId - 1; // b
    }
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
    getVertexNeighbors)() {
  if(vertexNeighborList_.empty()) {
    Timer t;
    vertexNeighborList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexNeighborList_[i].resize(getVertexNeighborNumber(i));
      for(SimplexId j = 0; j < (SimplexId)vertexNeighborList_[i].size(); ++j)
        getVertexNeighbor(i, j, vertexNeighborList_[i][j]);
    }

    printMsg("Built " + to_string(vertexNumber_) + " vertex neighbors.", 1,
             t.getElapsedTime(), 1);
  }

  return &vertexNeighborList_;
}

SimplexId PeriodicImplicitTriangulation::getVertexEdgeNumberInternal(
  const SimplexId &vertexId) const {
  return getVertexNeighborNumber(vertexId);
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getVertexEdgeInternal(
  const SimplexId &vertexId, const int &localEdgeId, SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localEdgeId < 0 or localEdgeId >= getVertexEdgeNumberInternal(vertexId))
    return -1;
#endif

  edgeId = -1;
  const auto &p = this->underlying().getVertexCoords(vertexId);

  if(dimensionality_ == 3) {
    edgeId = getVertexEdge3d(p.data(), localEdgeId);
  } else if(dimensionality_ == 2) {
    edgeId = getVertexEdge2d(p.data(), localEdgeId);
  } else if(dimensionality_ == 1) {
    if(vertexId > 0 and vertexId < nbvoxels_[Di_]) {
      // ab
      edgeId = localEdgeId == 0 ? vertexId : vertexId - 1;
    } else if(vertexId == 0) {
      // a
      edgeId = localEdgeId == 0 ? vertexId : 0;
    } else {
      // b
      edgeId = localEdgeId == 0 ? 0 : vertexId - 1;
    }
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getVertexEdgesInternal() {
  if(vertexEdgeList_.empty()) {
    Timer t;

    vertexEdgeList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexEdgeList_[i].resize(getVertexEdgeNumberInternal(i));
      for(SimplexId j = 0; j < (SimplexId)vertexEdgeList_[i].size(); ++j)
        getVertexEdgeInternal(i, j, vertexEdgeList_[i][j]);
    }

    printMsg("Built " + to_string(vertexNumber_) + " vertex edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &vertexEdgeList_;
}

SimplexId PeriodicImplicitTriangulation::getVertexTriangleNumberInternal(
  const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#else
  TTK_FORCE_USE(vertexId);
#endif

  if(dimensionality_ == 3) {
    return 36; // abcdefgh
  }

  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getVertexTriangleInternal(
  const SimplexId &vertexId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= getVertexTriangleNumberInternal(vertexId))
    return -1;
#endif
  triangleId = -1;

  const auto &p = this->underlying().getVertexCoords(vertexId);

  if(dimensionality_ == 3) {
    triangleId = getVertexTriangle3d(p.data(), localTriangleId);
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getVertexTrianglesInternal() {
  if(vertexTriangleList_.empty()) {
    Timer t;

    vertexTriangleList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexTriangleList_[i].resize(getVertexTriangleNumberInternal(i));
      for(SimplexId j = 0; j < (SimplexId)vertexTriangleList_[i].size(); ++j)
        getVertexTriangleInternal(i, j, vertexTriangleList_[i][j]);
    }

    printMsg("Built " + to_string(vertexNumber_) + " vertex triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &vertexTriangleList_;
}

SimplexId PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getVertexLinkNumber)(const SimplexId &vertexId) const {
  return getVertexStarNumber(vertexId);
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexLink)(const SimplexId &vertexId,
                 const int &localLinkId,
                 SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getVertexLinkNumber(vertexId))
    return -1;
#endif

  linkId = -1;
  const auto &p = this->underlying().getVertexCoords(vertexId);

  if(dimensionality_ == 3) {
    linkId = getVertexLink3d(p.data(), localLinkId);
  } else if(dimensionality_ == 2) {
    linkId = getVertexLink2d(p.data(), localLinkId); // abcd
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexLinks)() {
  if(vertexLinkList_.empty()) {
    Timer t;

    vertexLinkList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexLinkList_[i].resize(getVertexLinkNumber(i));
      for(SimplexId j = 0; j < (SimplexId)vertexLinkList_[i].size(); ++j)
        getVertexLink(i, j, vertexLinkList_[i][j]);
    }

    printMsg("Built " + to_string(vertexNumber_) + " vertex links.", 1,
             t.getElapsedTime(), 1);
  }

  return &vertexLinkList_;
}

SimplexId PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getVertexStarNumber)(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#else
  TTK_FORCE_USE(vertexId);
#endif

  if(dimensionality_ == 3) {
    return 24; // abcdefgh
  } else if(dimensionality_ == 2) {
    return 6; // abcd
  }

  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexStar)(const SimplexId &vertexId,
                 const int &localStarId,
                 SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getVertexStarNumber(vertexId))
    return -1;
#endif

  starId = -1;
  const auto &p = this->underlying().getVertexCoords(vertexId);

  if(dimensionality_ == 3) {
    starId = getVertexStar3d(p.data(), localStarId);
  } else if(dimensionality_ == 2) {
    starId = getVertexStar2d(p.data(), localStarId);
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexStars)() {
  if(vertexStarList_.empty()) {
    Timer t;
    vertexStarList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexStarList_[i].resize(getVertexStarNumber(i));
      for(SimplexId j = 0; j < (SimplexId)vertexStarList_[i].size(); ++j)
        getVertexStar(i, j, vertexStarList_[i][j]);
    }

    printMsg("Built " + to_string(vertexNumber_) + " vertex stars.", 1,
             t.getElapsedTime(), 1);
  }

  return &vertexStarList_;
}

#ifdef TTK_ENABLE_MPI

void PeriodicImplicitTriangulation::setIsBoundaryPeriodic(
  std::array<unsigned char, 6> boundary) {
  this->isBoundaryPeriodic = boundary;
}

void PeriodicImplicitTriangulation::createMetaGrid(const double *const bounds) {
  // only works with 2 processes or more
  if(!ttk::isRunningWithMPI()) {
    return;
  }

  // no need to create it anew?
  if(this->metaGrid_ != nullptr) {
    return;
  }

  // Reorganize bounds to only execute Allreduce twice
  std::array<double, 6> tempBounds = {
    bounds[0], bounds[2], bounds[4], bounds[1], bounds[3], bounds[5],
  };

  for(int i = 0; i < 3; i++) {
    if(dimensionality_ > i) {
      if(this->isBoundaryPeriodic[2 * i] == 1) {
        tempBounds[i] += spacing_[i];
      }
      if(this->isBoundaryPeriodic[2 * i + 1] == 1) {
        tempBounds[3 + i] -= spacing_[i];
      }
    }
  }

  std::array<double, 6> tempGlobalBounds{};
  // Compute and send to all processes the lower bounds of the data set
  MPI_Allreduce(tempBounds.data(), tempGlobalBounds.data(), 3, MPI_DOUBLE,
                MPI_MIN, ttk::MPIcomm_);
  // Compute and send to all processes the higher bounds of the data set
  MPI_Allreduce(&tempBounds[3], &tempGlobalBounds[3], 3, MPI_DOUBLE, MPI_MAX,
                ttk::MPIcomm_);

  // re-order tempGlobalBounds
  std::array<double, 6> globalBounds{
    tempGlobalBounds[0], tempGlobalBounds[3], tempGlobalBounds[1],
    tempGlobalBounds[4], tempGlobalBounds[2], tempGlobalBounds[5],
  };

  const std::array<ttk::SimplexId, 3> dimensions = {
    static_cast<ttk::SimplexId>(
      std::round((globalBounds[1] - globalBounds[0]) / this->spacing_[0]))
      + 1,
    static_cast<ttk::SimplexId>(
      std::round((globalBounds[3] - globalBounds[2]) / this->spacing_[1]))
      + 1,
    static_cast<ttk::SimplexId>(
      std::round((globalBounds[5] - globalBounds[4]) / this->spacing_[2]))
      + 1,
  };

  this->localGridOffset_ = {
    static_cast<SimplexId>(
      std::round((this->origin_[0] - globalBounds[0]) / this->spacing_[0])),
    static_cast<SimplexId>(
      std::round((this->origin_[1] - globalBounds[2]) / this->spacing_[1])),
    static_cast<SimplexId>(
      std::round((this->origin_[2] - globalBounds[4]) / this->spacing_[2])),
  };

  this->metaGrid_ = std::make_shared<PeriodicNoPreconditions>();
  this->metaGrid_->setInputGrid(globalBounds[0], globalBounds[1],
                                globalBounds[2], this->spacing_[0],
                                this->spacing_[1], this->spacing_[2],
                                dimensions[0], dimensions[1], dimensions[2]);
}

int PeriodicImplicitTriangulation::preconditionDistributedCells() {
  if(this->hasPreconditionedDistributedCells_) {
    return 0;
  }
  if(!ttk::hasInitializedMPI()) {
    return -1;
  }
  if(this->metaGrid_ == nullptr) {
    return 0;
  }
  if(this->cellGhost_ == nullptr) {
    if(ttk::isRunningWithMPI()) {
      this->printErr("Missing cell ghost array!");
    }
    return -3;
  }

  Timer tm{};

  // there are 6 tetrahedra per cubic cell (and 2 triangles per square)
  const int nTetraPerCube{this->dimensionality_ == 3 ? 6 : 2};

  // number of local cells (with ghost cells but without the additional periodic
  // cells)
  const auto nLocCells{(this->dimensions_[0] - 1) * (this->dimensions_[1] - 1)
                       * (this->dimensions_[2] - 1) * nTetraPerCube};

  std::vector<unsigned char> fillCells(nLocCells / nTetraPerCube);

  this->neighborCellBBoxes_.resize(ttk::MPIsize_);
  auto &localBBox{this->neighborCellBBoxes_[ttk::MPIrank_]};
  // "good" starting values?
  ttk::SimplexId localBBox_x_min{this->localGridOffset_[0]
                                 + this->dimensions_[0]},
    localBBox_y_min{this->localGridOffset_[1] + this->dimensions_[1]},
    localBBox_z_min{this->localGridOffset_[2] + this->dimensions_[2]};
  ttk::SimplexId localBBox_x_max{this->localGridOffset_[0]},
    localBBox_y_max{this->localGridOffset_[1]},
    localBBox_z_max{this->localGridOffset_[2]};
  const auto &dims{this->metaGrid_->getGridDimensions()};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for reduction(                    \
  min                                                  \
  : localBBox_x_min, localBBox_y_min, localBBox_z_min) \
  reduction(max                                        \
            : localBBox_x_max, localBBox_y_max, localBBox_z_max)
#endif
  for(SimplexId lcid = 0; lcid < nLocCells; ++lcid) {
    // only keep non-ghost cells
    if(this->cellGhost_[lcid / nTetraPerCube] != 0) {
      continue;
    }
    // local vertex coordinates
    std::array<SimplexId, 3> p{};
    if(this->dimensionality_ == 3) {
      this->tetrahedronToPosition(lcid, p.data());
    } else if(this->dimensionality_ == 2) {
      this->triangleToPosition2d(lcid, p.data());
      // compatibility with tetrahedronToPosition; fix a bounding box
      // error in the first axis
      p[0] /= 2;
    }

    // global vertex coordinates
    p[0] += this->localGridOffset_[0];
    p[1] += this->localGridOffset_[1];
    p[2] += this->localGridOffset_[2];

    if(p[0] < localBBox_x_min) {
      localBBox_x_min = max(p[0], static_cast<ttk::SimplexId>(0));
    }
    if(p[0] > localBBox_x_max) {
      localBBox_x_max = min(p[0], dims[0]);
    }
    if(p[1] < localBBox_y_min) {
      localBBox_y_min = max(p[1], static_cast<ttk::SimplexId>(0));
    }
    if(p[1] > localBBox_y_max) {
      localBBox_y_max = min(p[1], dims[1]);
    }
    if(p[2] < localBBox_z_min) {
      localBBox_z_min = max(p[2], static_cast<ttk::SimplexId>(0));
    }
    if(p[2] > localBBox_z_max) {
      localBBox_z_max = min(p[2], dims[2]);
    }
  }
  localBBox_x_min -= isBoundaryPeriodic[0];
  if(dimensionality_ > 1) {
    localBBox_y_min -= isBoundaryPeriodic[2];
    if(dimensionality_ > 2)
      localBBox_z_min -= isBoundaryPeriodic[4];
  }

  localBBox_x_max++;
  localBBox_y_max++;
  localBBox_z_max++;

  localBBox = {
    localBBox_x_min, localBBox_x_max, localBBox_y_min,
    localBBox_y_max, localBBox_z_min, localBBox_z_max,
  };

  for(size_t i = 0; i < this->neighborRanks_.size(); ++i) {
    const auto neigh{this->neighborRanks_[i]};
    MPI_Sendrecv(this->neighborCellBBoxes_[ttk::MPIrank_].data(), 6,
                 ttk::getMPIType(SimplexId{}), neigh, ttk::MPIrank_,
                 this->neighborCellBBoxes_[neigh].data(), 6,
                 ttk::getMPIType(SimplexId{}), neigh, neigh, ttk::MPIcomm_,
                 MPI_STATUS_IGNORE);
  }

  this->hasPreconditionedDistributedCells_ = true;

  return 0;
}

std::array<SimplexId, 3> PeriodicImplicitTriangulation::getVertGlobalCoords(
  const SimplexId lvid) const {
  // local vertex coordinates
  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {
    this->vertexToPosition(lvid, p.data());
  } else if(this->dimensionality_ == 2) {
    this->vertexToPosition2d(lvid, p.data());
  }
  // global vertex coordinates
  p[0] += this->localGridOffset_[0];
  p[1] += this->localGridOffset_[1];
  p[2] += this->localGridOffset_[2];

  const auto &dims{this->metaGrid_->getGridDimensions()};

  p[0] = (p[0] + dims[0]) % dims[0];
  if(dimensionality_ > 1) {
    p[1] = (p[1] + dims[1]) % dims[1];
    if(dimensionality_ > 2)
      p[2] = (p[2] + dims[2]) % dims[2];
  }

  return p;
}

std::array<SimplexId, 3> PeriodicImplicitTriangulation::getVertLocalCoords(
  const SimplexId gvid) const {
  // global vertex coordinates
  std::array<SimplexId, 3> pGlobal{};
  if(this->dimensionality_ == 3) {
    this->metaGrid_->vertexToPosition(gvid, pGlobal.data());
  } else if(this->dimensionality_ == 2) {
    this->metaGrid_->vertexToPosition2d(gvid, pGlobal.data());
  }
  std::array<SimplexId, 3> p{pGlobal};
  // local vertex coordinates
  p[0] -= this->localGridOffset_[0];
  p[1] -= this->localGridOffset_[1];
  p[2] -= this->localGridOffset_[2];

  const auto &dims{this->getGridDimensions()};

  if(p[0] >= 0 && p[1] >= 0 && p[2] >= 0 && p[0] <= dims[0] - 1
     && p[1] <= dims[1] - 1 && p[2] <= dims[2] - 1) {
    return p;
  }
  for(int i = 0; i < 3; i++) {
    if((p[i] < 0 || p[i] > dims[i] - 1) && pGlobal[i] == 0) {
      p[i] = dims[i] - 1;
    }
    if((p[i] < 0 || p[i] > dims[i] - 1)
       && pGlobal[i] == this->metaGrid_->dimensions_[i] - 1) {
      p[i] = 0;
    }
  }

  if(p[0] >= 0 && p[1] >= 0 && p[2] >= 0 && p[0] <= dims[0] - 1
     && p[1] <= dims[1] - 1 && p[2] <= dims[2] - 1) {
    if(this->vertexGhost_[p[0] + p[1] * dims[0] + p[2] * dims[0] * dims[1]]
       != 0) {
      return p;
    }
  }
  return std::array<SimplexId, 3>{-1, -1, -1};
}

int ttk::PeriodicImplicitTriangulation::getCellRankInternal(
  const SimplexId lcid) const {

  const int nTetraPerCube{this->dimensionality_ == 3 ? 6 : 2};
  const auto locCubeId{lcid / nTetraPerCube};

  if(this->cellGhost_[locCubeId] == 0) {
    return ttk::MPIrank_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(this->neighborRanks_.empty()) {
    this->printErr("Empty neighborsRanks_!");
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  const auto nVertsCell{this->getCellVertexNumber(lcid)};
  std::vector<bool> inRank(nVertsCell);
  std::map<int, int> neighborOccurrences;
  for(const auto neigh : this->neighborRanks_) {
    std::fill(inRank.begin(), inRank.end(), false);
    const auto &bbox{this->neighborCellBBoxes_[neigh]};
    for(SimplexId i = 0; i < nVertsCell; ++i) {
      SimplexId v{};
      this->getCellVertex(lcid, i, v);
      if(this->vertexGhost_[v] == 0) {
        inRank[i] = true;
      } else {
        const auto p{this->getVertGlobalCoords(v)};
        if(p[0] >= bbox[0] && p[0] <= bbox[1] && p[1] >= bbox[2]
           && p[1] <= bbox[3] && p[2] >= bbox[4] && p[2] <= bbox[5]) {
          inRank[i] = true;
        }
      }
    }
    if(std::all_of(
         inRank.begin(), inRank.end(), [](const bool v) { return v; })) {
      return neigh;
    }
    neighborOccurrences[neigh]
      = std::accumulate(inRank.begin(), inRank.end(), 0);
  }

  auto pr = std::max_element(
    std::begin(neighborOccurrences), std::end(neighborOccurrences),
    [](const std::pair<int, int> &p1, const std::pair<int, int> &p2) {
      return p1.second < p2.second;
    });
  return pr->first;
}

#endif // TTK_ENABLE_MPI

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexPoint)(const SimplexId &vertexId,
                  float &x,
                  float &y,
                  float &z) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    const auto &p = this->underlying().getVertexCoords(vertexId);

    x = origin_[0] + spacing_[0] * p[0];
    y = origin_[1] + spacing_[1] * p[1];
    z = origin_[2] + spacing_[2] * p[2];
  } else if(dimensionality_ == 2) {
    const auto &p = this->underlying().getVertexCoords(vertexId);

    if(dimensions_[0] > 1 and dimensions_[1] > 1) {
      x = origin_[0] + spacing_[0] * p[0];
      y = origin_[1] + spacing_[1] * p[1];
      z = origin_[2];
    } else if(dimensions_[1] > 1 and dimensions_[2] > 1) {
      x = origin_[0];
      y = origin_[1] + spacing_[1] * p[0];
      z = origin_[2] + spacing_[2] * p[1];
    } else if(dimensions_[0] > 1 and dimensions_[2] > 1) {
      x = origin_[0] + spacing_[0] * p[0];
      y = origin_[1];
      z = origin_[2] + spacing_[2] * p[1];
    }
  } else if(dimensionality_ == 1) {
    if(dimensions_[0] > 1) {
      x = origin_[0] + spacing_[0] * vertexId;
      y = origin_[1];
      z = origin_[2];
    } else if(dimensions_[1] > 1) {
      x = origin_[0];
      y = origin_[1] + spacing_[1] * vertexId;
      z = origin_[2];
    } else if(dimensions_[2] > 1) {
      x = origin_[0];
      y = origin_[1];
      z = origin_[2] + spacing_[2] * vertexId;
    }
  }

  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getEdgeVertexInternal(
  const SimplexId &edgeId,
  const int &localVertexId,
  SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 2)
    return -2;
#endif

  vertexId = -1;
  const auto &p = this->underlying().getEdgeCoords(edgeId);
  const SimplexId wrapXRight = (p[0] == nbvoxels_[0] ? -wrap_[0] : 0);
  const SimplexId wrapYBottom = (p[1] == nbvoxels_[1] ? -wrap_[1] : 0);
  const SimplexId wrapZFront = (p[2] == nbvoxels_[2] ? -wrap_[2] : 0);
  const auto a = p[0] + this->underlying().getEdgeVertexAccelerated(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
    case EdgePosition::L_3D:
      vertexId = a + (localVertexId == 0 ? 0 : (1 + wrapXRight));
      break;
    case EdgePosition::H_3D:
      vertexId = a + (localVertexId == 0 ? 0 : (vshift_[0] + wrapYBottom));
      break;
    case EdgePosition::P_3D:
      vertexId = a + (localVertexId == 0 ? 0 : (vshift_[1] + wrapZFront));
      break;
    case EdgePosition::D1_3D:
      vertexId = a
                 + (localVertexId == 0 ? (1 + wrapXRight)
                                       : (vshift_[0] + wrapYBottom));
      break;
    case EdgePosition::D2_3D:
      vertexId = a
                 + (localVertexId == 0
                      ? 0
                      : (vshift_[0] + wrapYBottom + vshift_[1] + wrapZFront));
      break;
    case EdgePosition::D3_3D:
      vertexId
        = a
          + (localVertexId == 0 ? (1 + wrapXRight) : (vshift_[1] + wrapZFront));
      break;
    case EdgePosition::D4_3D:
      vertexId = a
                 + (localVertexId == 0
                      ? (1 + wrapXRight)
                      : (vshift_[0] + wrapYBottom + vshift_[1] + wrapZFront));
      break;
    case EdgePosition::L_2D:
      vertexId = a + (localVertexId == 0 ? 0 : (1 + wrapXRight));
      break;
    case EdgePosition::H_2D:
      vertexId = a + (localVertexId == 0 ? 0 : (vshift_[0] + wrapYBottom));
      break;
    case EdgePosition::D1_2D:
      vertexId = a
                 + (localVertexId == 0 ? (1 + wrapXRight)
                                       : (vshift_[0] + wrapYBottom));
      break;
    case EdgePosition::FIRST_EDGE_1D:
      vertexId = localVertexId == 0 ? 0 : 1;
      break;
    case EdgePosition::LAST_EDGE_1D:
      vertexId = localVertexId == 0 ? edgeId : 0;
      break;
    case EdgePosition::CENTER_1D:
      vertexId = localVertexId == 0 ? edgeId : edgeId + 1;
      break;
    default:
      break;
  }

  return 0;
}

const vector<std::array<SimplexId, 2>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdges)() {

  if(edgeList_.empty()) {
    Timer t;

    edgeList_.resize(edgeNumber_);
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      SimplexId id0, id1;
      getEdgeVertexInternal(i, 0, id0);
      getEdgeVertexInternal(i, 1, id1);
      edgeList_[i] = {id0, id1};
    }

    printMsg(
      "Built " + to_string(edgeNumber_) + " edges.", 1, t.getElapsedTime(), 1);
  }

  return &edgeList_;
}

template <typename Derived>
SimplexId
  PeriodicImplicitTriangulationCRTP<Derived>::getEdgeTriangleNumberInternal(
    const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  switch(this->underlying().getEdgePosition(edgeId)) {
    case EdgePosition::L_3D:
    case EdgePosition::H_3D:
    case EdgePosition::P_3D:
    case EdgePosition::D4_3D:
      return 6;
    case EdgePosition::D1_3D:
    case EdgePosition::D2_3D:
    case EdgePosition::D3_3D:
      return 4;
    case EdgePosition::L_2D:
    case EdgePosition::H_2D:
    case EdgePosition::D1_2D:
      return 2;
    default:
      return 0;
  }
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getEdgeTriangleInternal(
  const SimplexId &edgeId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= getEdgeTriangleNumberInternal(edgeId))
    return -1;
#endif

  triangleId = -1;
  const auto &p = this->underlying().getEdgeCoords(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
    case EdgePosition::L_3D:
      triangleId = getEdgeTriangle3dL(p.data(), localTriangleId);
      break;
    case EdgePosition::H_3D:
      triangleId = getEdgeTriangle3dH(p.data(), localTriangleId);
      break;
    case EdgePosition::P_3D:
      triangleId = getEdgeTriangle3dP(p.data(), localTriangleId);
      break;
    case EdgePosition::D1_3D:
      triangleId = getEdgeTriangle3dD1(p.data(), localTriangleId);
      break;
    case EdgePosition::D2_3D:
      triangleId = getEdgeTriangle3dD2(p.data(), localTriangleId);
      break;
    case EdgePosition::D3_3D:
      triangleId = getEdgeTriangle3dD3(p.data(), localTriangleId);
      break;
    case EdgePosition::D4_3D:
      triangleId = getEdgeTriangle3dD4(p.data(), localTriangleId);
      break;
    case EdgePosition::L_2D:
      triangleId = getEdgeTriangle2dL(p.data(), localTriangleId);
      break;
    case EdgePosition::H_2D:
      triangleId = getEdgeTriangle2dH(p.data(), localTriangleId);
      break;
    case EdgePosition::D1_2D:
      triangleId = getEdgeTriangle2dD1(p.data(), localTriangleId);
      break;
    default:
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getEdgeTrianglesInternal() {
  if(edgeTriangleList_.empty()) {
    Timer t;

    edgeTriangleList_.resize(edgeNumber_);
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      edgeTriangleList_[i].resize(getEdgeTriangleNumberInternal(i));
      for(SimplexId j = 0; j < (SimplexId)edgeTriangleList_[i].size(); ++j)
        getEdgeTriangleInternal(i, j, edgeTriangleList_[i][j]);
    }

    printMsg("Built " + to_string(edgeNumber_) + " edge triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &edgeTriangleList_;
}

SimplexId PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getEdgeLinkNumber)(const SimplexId &edgeId) const {

  return getEdgeStarNumber(edgeId);
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::TTK_TRIANGULATION_INTERNAL(
  getEdgeLink)(const SimplexId &edgeId,
               const int &localLinkId,
               SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getEdgeLinkNumber(edgeId))
    return -1;
#endif

  linkId = -1;
  const auto &p = this->underlying().getEdgeCoords(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
    case EdgePosition::L_3D:
      linkId = getEdgeLinkL(p.data(), localLinkId);
      break;
    case EdgePosition::H_3D:
      linkId = getEdgeLinkH(p.data(), localLinkId);
      break;
    case EdgePosition::P_3D:
      linkId = getEdgeLinkP(p.data(), localLinkId);
      break;
    case EdgePosition::D1_3D:
      linkId = getEdgeLinkD1(p.data(), localLinkId);
      break;
    case EdgePosition::D2_3D:
      linkId = getEdgeLinkD2(p.data(), localLinkId);
      break;
    case EdgePosition::D3_3D:
      linkId = getEdgeLinkD3(p.data(), localLinkId);
      break;
    case EdgePosition::D4_3D:
      linkId = getEdgeLinkD4(p.data(), localLinkId);
      break;
    case EdgePosition::L_2D:
      linkId = getEdgeLink2dL(p.data(), localLinkId);
      break;
    case EdgePosition::H_2D:
      linkId = getEdgeLink2dH(p.data(), localLinkId);
      break;
    case EdgePosition::D1_2D:
      linkId = getEdgeLink2dD1(p.data(), localLinkId);
      break;
    default:
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() {

  if(edgeLinkList_.empty()) {
    Timer t;

    edgeLinkList_.resize(edgeNumber_);
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      edgeLinkList_[i].resize(getEdgeLinkNumber(i));
      for(SimplexId j = 0; j < (SimplexId)edgeLinkList_[i].size(); ++j)
        getEdgeLink(i, j, edgeLinkList_[i][j]);
    }

    printMsg("Built " + to_string(edgeNumber_) + " edge links.", 1,
             t.getElapsedTime(), 1);
  }

  return &edgeLinkList_;
}

template <typename Derived>
SimplexId
  PeriodicImplicitTriangulationCRTP<Derived>::TTK_TRIANGULATION_INTERNAL(
    getEdgeStarNumber)(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  switch(this->underlying().getEdgePosition(edgeId)) {
    case EdgePosition::L_3D:
    case EdgePosition::H_3D:
    case EdgePosition::P_3D:
    case EdgePosition::D4_3D:
      return 6;
    case EdgePosition::D1_3D:
    case EdgePosition::D2_3D:
    case EdgePosition::D3_3D:
      return 4;
    case EdgePosition::L_2D:
    case EdgePosition::H_2D:
    case EdgePosition::D1_2D:
      return 2;
    default:
      return 0;
  }

  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::TTK_TRIANGULATION_INTERNAL(
  getEdgeStar)(const SimplexId &edgeId,
               const int &localStarId,
               SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getEdgeStarNumber(edgeId))
    return -1;
#endif

  starId = -1;
  const auto &p = this->underlying().getEdgeCoords(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
    case EdgePosition::L_3D:
      starId = getEdgeStarL(p.data(), localStarId);
      break;
    case EdgePosition::H_3D:
      starId = getEdgeStarH(p.data(), localStarId);
      break;
    case EdgePosition::P_3D:
      starId = getEdgeStarP(p.data(), localStarId);
      break;
    case EdgePosition::D1_3D:
      starId = getEdgeStarD1(p.data(), localStarId);
      break;
    case EdgePosition::D2_3D:
      starId = getEdgeStarD2(p.data(), localStarId);
      break;
    case EdgePosition::D3_3D:
      starId = getEdgeStarD3(p.data(), localStarId);
      break;
    case EdgePosition::D4_3D:
      starId
        = p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + localStarId;
      break;
    case EdgePosition::L_2D:
      starId = getEdgeStar2dL(p.data(), localStarId);
      break;
    case EdgePosition::H_2D:
      starId = getEdgeStar2dH(p.data(), localStarId);
      break;
    case EdgePosition::D1_2D:
      starId = p[0] * 2 + p[1] * tshift_[0] + localStarId;
      break;
    default:
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdgeStars)() {

  if(edgeStarList_.empty()) {
    Timer t;

    edgeStarList_.resize(edgeNumber_);
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      edgeStarList_[i].resize(getEdgeStarNumber(i));
      for(SimplexId j = 0; j < (SimplexId)edgeStarList_[i].size(); ++j)
        getEdgeStar(i, j, edgeStarList_[i][j]);
    }

    printMsg("Built " + to_string(edgeNumber_) + " edge stars.", 1,
             t.getElapsedTime(), 1);
  }

  return &edgeStarList_;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getTriangleVertexInternal(
  const SimplexId &triangleId,
  const int &localVertexId,
  SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 3)
    return -2;
#endif

  vertexId = -1;
  const auto &p = this->underlying().getTriangleCoords(triangleId);
  const SimplexId wrapXRight = (p[0] / 2 == nbvoxels_[Di_]) ? -wrap_[0] : 0;
  const SimplexId wrapYBottom = (p[1] == nbvoxels_[Dj_]) ? -wrap_[1] : 0;

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case TrianglePosition::F_3D:
      vertexId = getTriangleVertexF(p.data(), localVertexId);
      break;
    case TrianglePosition::H_3D:
      vertexId = getTriangleVertexH(p.data(), localVertexId);
      break;
    case TrianglePosition::C_3D:
      vertexId = getTriangleVertexC(p.data(), localVertexId);
      break;
    case TrianglePosition::D1_3D:
      vertexId = getTriangleVertexD1(p.data(), localVertexId);
      break;
    case TrianglePosition::D2_3D:
      vertexId = getTriangleVertexD2(p.data(), localVertexId);
      break;
    case TrianglePosition::D3_3D:
      vertexId = getTriangleVertexD3(p.data(), localVertexId);
      break;
    case TrianglePosition::TOP_2D:
      if(localVertexId == 0) {
        vertexId = p[0] / 2 + p[1] * vshift_[0];
      } else if(localVertexId == 1) {
        vertexId = p[0] / 2 + p[1] * vshift_[0] + 1 + wrapXRight;
      } else if(localVertexId == 2) {
        vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0] + wrapYBottom;
      }
      break;
    case TrianglePosition::BOTTOM_2D:
      if(localVertexId == 0) {
        vertexId = p[0] / 2 + p[1] * vshift_[0] + 1 + wrapXRight;
      } else if(localVertexId == 1) {
        vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0] + 1 + wrapXRight
                   + wrapYBottom;
      } else if(localVertexId == 2) {
        vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0] + wrapYBottom;
      }
      break;
    default:
      break;
  }

  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getTriangleEdgeInternal(
  const SimplexId &triangleId,
  const int &localEdgeId,
  SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
  if(localEdgeId < 0 or localEdgeId >= 3)
    return -2;
#endif

  edgeId = -1;
  const auto &p = this->underlying().getTriangleCoords(triangleId);
  const SimplexId wrapXRight = (p[0] / 2 == nbvoxels_[Di_]) ? -wrap_[0] : 0;
  const SimplexId wrapYBottom = (p[1] == nbvoxels_[Dj_]) ? -wrap_[1] : 0;
  const SimplexId id = triangleId % 2;

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case TrianglePosition::F_3D:
      edgeId = (id == 1) ? getTriangleEdgeF_1(p.data(), localEdgeId)
                         : getTriangleEdgeF_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::H_3D:
      edgeId = (id == 1) ? getTriangleEdgeH_1(p.data(), localEdgeId)
                         : getTriangleEdgeH_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::C_3D:
      edgeId = (id == 1) ? getTriangleEdgeC_1(p.data(), localEdgeId)
                         : getTriangleEdgeC_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::D1_3D:
      edgeId = (id == 1) ? getTriangleEdgeD1_1(p.data(), localEdgeId)
                         : getTriangleEdgeD1_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::D2_3D:
      edgeId = (id == 1) ? getTriangleEdgeD2_1(p.data(), localEdgeId)
                         : getTriangleEdgeD2_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::D3_3D:
      edgeId = (id == 1) ? getTriangleEdgeD3_1(p.data(), localEdgeId)
                         : getTriangleEdgeD3_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::TOP_2D:
      if(localEdgeId == 0) {
        edgeId = p[0] / 2 + p[1] * eshift_[0];
      } else if(localEdgeId == 1) {
        edgeId = esetshift_[0] + p[0] / 2 + p[1] * eshift_[2];
      } else if(localEdgeId == 2) {
        edgeId = esetshift_[1] + p[0] / 2 + p[1] * eshift_[4];
      }
      break;
    case TrianglePosition::BOTTOM_2D:
      if(localEdgeId == 0) {
        edgeId = p[0] / 2 + (p[1] + 1) * eshift_[0] + wrapYBottom;
      } else if(localEdgeId == 1) {
        edgeId
          = esetshift_[0] + (p[0] + 1) / 2 + p[1] * eshift_[2] + wrapXRight;
      } else if(localEdgeId == 2) {
        edgeId = esetshift_[1] + p[0] / 2 + p[1] * eshift_[4];
      }
      break;
    default:
      break;
  }

  return 0;
}

int PeriodicImplicitTriangulation::getTriangleEdgesInternal(
  vector<vector<SimplexId>> &edges) const {
  edges.resize(triangleNumber_);
  for(SimplexId i = 0; i < triangleNumber_; ++i) {
    edges[i].resize(3);
    for(int j = 0; j < 3; ++j)
      getTriangleEdgeInternal(i, j, edges[i][j]);
  }
  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getTriangleEdgesInternal() {
  if(triangleEdgeVector_.empty()) {
    Timer t;

    getTriangleEdgesInternal(triangleEdgeVector_);

    printMsg("Built " + to_string(triangleNumber_) + " triangle edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &triangleEdgeVector_;
}

const vector<std::array<SimplexId, 3>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangles)() {

  if(triangleList_.empty()) {
    Timer t;

    triangleList_.resize(triangleNumber_);
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      for(int j = 0; j < 3; ++j)
        getTriangleVertexInternal(i, j, triangleList_[i][j]);
    }

    printMsg("Built " + to_string(triangleNumber_) + " triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &triangleList_;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::TTK_TRIANGULATION_INTERNAL(
  getTriangleLink)(const SimplexId &triangleId,
                   const int &localLinkId,
                   SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getTriangleLinkNumber(triangleId))
    return -1;
#endif

  linkId = -1;
  const auto &p = this->underlying().getTriangleCoords(triangleId);

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case TrianglePosition::F_3D:
      linkId = getTriangleLinkF(p.data(), localLinkId);
      break;
    case TrianglePosition::H_3D:
      linkId = getTriangleLinkH(p.data(), localLinkId);
      break;
    case TrianglePosition::C_3D:
      linkId = getTriangleLinkC(p.data(), localLinkId);
      break;
    case TrianglePosition::D1_3D:
      linkId = getTriangleLinkD1(p.data(), localLinkId);
      break;
    case TrianglePosition::D2_3D:
      linkId = getTriangleLinkD2(p.data(), localLinkId);
      break;
    case TrianglePosition::D3_3D:
      linkId = getTriangleLinkD3(p.data(), localLinkId);
      break;
    default:
      break;
  }

  return 0;
}

SimplexId PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getTriangleLinkNumber)(const SimplexId &triangleId) const {

  return getTriangleStarNumber(triangleId);
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
    getTriangleLinks)() {

  if(triangleLinkList_.empty()) {
    Timer t;

    triangleLinkList_.resize(triangleNumber_);
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      triangleLinkList_[i].resize(getTriangleLinkNumber(i));
      for(SimplexId j = 0; j < (SimplexId)triangleLinkList_[i].size(); ++j)
        getTriangleLink(i, j, triangleLinkList_[i][j]);
    }

    printMsg("Built " + to_string(triangleNumber_) + " triangle links.", 1,
             t.getElapsedTime(), 1);
  }
  return &triangleLinkList_;
}

template <typename Derived>
SimplexId
  PeriodicImplicitTriangulationCRTP<Derived>::TTK_TRIANGULATION_INTERNAL(
    getTriangleStarNumber)(const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#else
  TTK_FORCE_USE(triangleId);
#endif

  if(dimensionality_ == 3) {
    return 2;
  }
  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::TTK_TRIANGULATION_INTERNAL(
  getTriangleStar)(const SimplexId &triangleId,
                   const int &localStarId,
                   SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getTriangleStarNumber(triangleId))
    return -1;
#endif

  starId = -1;
  const auto &p = this->underlying().getTriangleCoords(triangleId);

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case TrianglePosition::F_3D:
      starId = getTriangleStarF(p.data(), localStarId);
      break;
    case TrianglePosition::H_3D:
      starId = getTriangleStarH(p.data(), localStarId);
      break;
    case TrianglePosition::C_3D:
      starId = getTriangleStarC(p.data(), localStarId);
      break;
    case TrianglePosition::D1_3D:
      starId = getTriangleStarD1(p.data(), localStarId);
      break;
    case TrianglePosition::D2_3D:
      starId = getTriangleStarD2(p.data(), localStarId);
      break;
    case TrianglePosition::D3_3D:
      starId = getTriangleStarD3(p.data(), localStarId);
      break;
    default:
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
    getTriangleStars)() {

  if(triangleStarList_.empty()) {
    Timer t;

    triangleStarList_.resize(triangleNumber_);
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      triangleStarList_[i].resize(getTriangleStarNumber(i));
      for(SimplexId j = 0; j < (SimplexId)triangleStarList_[i].size(); ++j)
        getTriangleStar(i, j, triangleStarList_[i][j]);
    }

    printMsg("Built " + to_string(triangleNumber_) + " triangle stars.", 1,
             t.getElapsedTime(), 1);
  }
  return &triangleStarList_;
}

SimplexId PeriodicImplicitTriangulation::getTriangleNeighborNumber(
  const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  TTK_FORCE_USE(triangleId);
  if(dimensionality_ == 2) {
    return 3;
  }

  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getTriangleNeighbor(
  const SimplexId &triangleId,
  const int &localNeighborId,
  SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= getTriangleNeighborNumber(triangleId))
    return -1;
#endif

  neighborId = -1;
  const auto &p = this->underlying().getTriangleCoords(triangleId);

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case TrianglePosition::BOTTOM_2D:

      if(p[0] / 2 == nbvoxels_[Di_] and p[1] == nbvoxels_[Dj_]) {
        if(localNeighborId == 0) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId + 1 - wrap_[0] * 2;
        } else if(localNeighborId == 2) {
          neighborId = triangleId + tshift_[0] - 1 - wrap_[1] * 2;
        }

      } else if(p[0] / 2 == nbvoxels_[Di_]) {
        if(localNeighborId == 0) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId + 1 - wrap_[0] * 2;
        } else if(localNeighborId == 2) {
          neighborId = triangleId + tshift_[0] - 1;
        }

      } else if(p[1] == nbvoxels_[Dj_]) {
        if(localNeighborId == 0) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 2) {
          neighborId = triangleId + tshift_[0] - 1 - wrap_[1] * 2;
        }

      } else {
        if(localNeighborId == 0) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 2) {
          neighborId = triangleId + tshift_[0] - 1;
        }
      }
      break;

    case TrianglePosition::TOP_2D:

      if(p[0] / 2 == 0 and p[1] == 0) {
        if(localNeighborId == 0) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId - 1 + wrap_[0] * 2;
        } else if(localNeighborId == 2) {
          neighborId = triangleId - tshift_[0] + 1 + wrap_[1] * 2;
        }

      } else if(p[0] / 2 == 0) {
        if(localNeighborId == 0) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId - 1 + wrap_[0] * 2;
        } else if(localNeighborId == 2) {
          neighborId = triangleId - tshift_[0] + 1;
        }

      } else if(p[1] == 0) {
        if(localNeighborId == 0) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 2) {
          neighborId = triangleId - tshift_[0] + 1 + wrap_[1] * 2;
        }

      } else {
        if(localNeighborId == 0) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 2) {
          neighborId = triangleId - tshift_[0] + 1;
        }
      }
      break;
    default:
      break;
  }

  return 0;
}

int PeriodicImplicitTriangulation::getTriangleNeighbors(
  vector<vector<SimplexId>> &neighbors) {
  neighbors.resize(triangleNumber_);
  for(SimplexId i = 0; i < triangleNumber_; ++i) {
    neighbors[i].resize(getTriangleNeighborNumber(i));
    for(SimplexId j = 0; j < (SimplexId)neighbors[i].size(); ++j)
      getTriangleNeighbor(i, j, neighbors[i][j]);
  }
  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getTetrahedronVertex(
  const SimplexId &tetId, const int &localVertexId, SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 4)
    return -2;
#endif

  vertexId = -1;

  if(dimensionality_ == 3) {
    const auto &p = this->underlying().getTetrahedronCoords(tetId);
    const SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        vertexId = getTetrahedronVertexABCG(p.data(), localVertexId);
        break;
      case 1:
        vertexId = getTetrahedronVertexBCDG(p.data(), localVertexId);
        break;
      case 2:
        vertexId = getTetrahedronVertexABEG(p.data(), localVertexId);
        break;
      case 3:
        vertexId = getTetrahedronVertexBEFG(p.data(), localVertexId);
        break;
      case 4:
        vertexId = getTetrahedronVertexBFGH(p.data(), localVertexId);
        break;
      case 5:
        vertexId = getTetrahedronVertexBDGH(p.data(), localVertexId);
        break;
    }
  }
  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getTetrahedronEdge(
  const SimplexId &tetId, const int &localEdgeId, SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localEdgeId < 0 or localEdgeId >= 6)
    return -2;
#endif

  edgeId = -1;

  if(dimensionality_ == 3) {
    const auto &p = this->underlying().getTetrahedronCoords(tetId);
    const SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        edgeId = getTetrahedronEdgeABCG(p.data(), localEdgeId);
        break;
      case 1:
        edgeId = getTetrahedronEdgeBCDG(p.data(), localEdgeId);
        break;
      case 2:
        edgeId = getTetrahedronEdgeABEG(p.data(), localEdgeId);
        break;
      case 3:
        edgeId = getTetrahedronEdgeBEFG(p.data(), localEdgeId);
        break;
      case 4:
        edgeId = getTetrahedronEdgeBFGH(p.data(), localEdgeId);
        break;
      case 5:
        edgeId = getTetrahedronEdgeBDGH(p.data(), localEdgeId);
        break;
    }
  }

  return 0;
}

int PeriodicImplicitTriangulation::getTetrahedronEdges(
  vector<vector<SimplexId>> &edges) const {
  edges.resize(tetrahedronNumber_);
  for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    edges[i].resize(6);
    for(int j = 0; j < 6; ++j)
      getTetrahedronEdge(i, j, edges[i][j]);
  }

  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getTetrahedronTriangle(
  const SimplexId &tetId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localTriangleId < 0 or localTriangleId >= 4)
    return -2;
#endif

  triangleId = -1;

  if(dimensionality_ == 3) {
    const auto &p = this->underlying().getTetrahedronCoords(tetId);
    const SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        triangleId = getTetrahedronTriangleABCG(p.data(), localTriangleId);
        break;
      case 1:
        triangleId = getTetrahedronTriangleBCDG(p.data(), localTriangleId);
        break;
      case 2:
        triangleId = getTetrahedronTriangleABEG(p.data(), localTriangleId);
        break;
      case 3:
        triangleId = getTetrahedronTriangleBEFG(p.data(), localTriangleId);
        break;
      case 4:
        triangleId = getTetrahedronTriangleBFGH(p.data(), localTriangleId);
        break;
      case 5:
        triangleId = getTetrahedronTriangleBDGH(p.data(), localTriangleId);
        break;
    }
  }

  return 0;
}

int PeriodicImplicitTriangulation::getTetrahedronTriangles(
  vector<vector<SimplexId>> &triangles) const {
  triangles.resize(tetrahedronNumber_);
  for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    triangles[i].resize(4);
    for(int j = 0; j < 4; ++j)
      getTetrahedronTriangle(i, j, triangles[i][j]);
  }

  return 0;
}

SimplexId PeriodicImplicitTriangulation::getTetrahedronNeighborNumber(
  const SimplexId &tetId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
#endif

  TTK_FORCE_USE(tetId);
  if(dimensionality_ == 3) {
    return 4;
  }

  return 0;
}

template <typename Derived>
int PeriodicImplicitTriangulationCRTP<Derived>::getTetrahedronNeighbor(
  const SimplexId &tetId,
  const int &localNeighborId,
  SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= getTetrahedronNeighborNumber(tetId))
    return -1;
#endif

  neighborId = -1;

  if(dimensionality_ == 3) {
    const auto &p = this->underlying().getTetrahedronCoords(tetId);
    const SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        neighborId
          = getTetrahedronNeighborABCG(tetId, p.data(), localNeighborId);
        break;
      case 1:
        neighborId
          = getTetrahedronNeighborBCDG(tetId, p.data(), localNeighborId);
        break;
      case 2:
        neighborId
          = getTetrahedronNeighborABEG(tetId, p.data(), localNeighborId);
        break;
      case 3:
        neighborId
          = getTetrahedronNeighborBEFG(tetId, p.data(), localNeighborId);
        break;
      case 4:
        neighborId
          = getTetrahedronNeighborBFGH(tetId, p.data(), localNeighborId);
        break;
      case 5:
        neighborId
          = getTetrahedronNeighborBDGH(tetId, p.data(), localNeighborId);
        break;
    }
  }

  return 0;
}

int PeriodicImplicitTriangulation::getTetrahedronNeighbors(
  vector<vector<SimplexId>> &neighbors) {
  neighbors.resize(tetrahedronNumber_);
  for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    neighbors[i].resize(getTetrahedronNeighborNumber(i));
    for(SimplexId j = 0; j < (SimplexId)neighbors[i].size(); ++j)
      getTetrahedronNeighbor(i, j, neighbors[i][j]);
  }

  return 0;
}

SimplexId PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getCellVertexNumber)(const SimplexId &ttkNotUsed(cellId)) const {

  return dimensionality_ + 1;
}

int PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getCellVertex)(
  const SimplexId &cellId,
  const int &localVertexId,
  SimplexId &vertexId) const {

  if(dimensionality_ == 3)
    getTetrahedronVertex(cellId, localVertexId, vertexId);
  else if(dimensionality_ == 2)
    getTriangleVertexInternal(cellId, localVertexId, vertexId);
  else if(dimensionality_ == 1)
    getEdgeVertexInternal(cellId, localVertexId, vertexId);

  return 0;
}

SimplexId PeriodicImplicitTriangulation::getCellEdgeNumberInternal(
  const SimplexId &ttkNotUsed(cellId)) const {
  if(dimensionality_ == 3)
    return 6;
  else if(dimensionality_ == 2)
    return 3;

  return 0;
}

int PeriodicImplicitTriangulation::getCellEdgeInternal(
  const SimplexId &cellId, const int &localEdgeId, SimplexId &edgeId) const {
  if(dimensionality_ == 3)
    getTetrahedronEdge(cellId, localEdgeId, edgeId);
  else if(dimensionality_ == 2)
    getTriangleEdgeInternal(cellId, localEdgeId, edgeId);
  else if(dimensionality_ == 1)
    getCellNeighbor(cellId, localEdgeId, edgeId);

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getCellEdgesInternal() {
  if(cellEdgeVector_.empty()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronEdges(cellEdgeVector_);
    else if(dimensionality_ == 2)
      getTriangleEdgesInternal(cellEdgeVector_);

    printMsg("Built " + to_string(cellNumber_) + " cell edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &cellEdgeVector_;
}

int PeriodicImplicitTriangulation::getCellTriangleInternal(
  const SimplexId &cellId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
  if(dimensionality_ == 3)
    getTetrahedronTriangle(cellId, localTriangleId, triangleId);

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getCellTrianglesInternal() {
  if(cellTriangleVector_.empty()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronTriangles(cellTriangleVector_);

    printMsg("Built " + to_string(cellNumber_) + " cell triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &cellTriangleVector_;
}

SimplexId PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getCellNeighborNumber)(const SimplexId &cellId) const {

  if(dimensionality_ == 3)
    return getTetrahedronNeighborNumber(cellId);
  else if(dimensionality_ == 2)
    return getTriangleNeighborNumber(cellId);
  else if(dimensionality_ == 1) {
    printErr("getCellNeighborNumber() not implemented in 1D! (TODO)");
    return -1;
  }

  return 0;
}

int PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
  const SimplexId &cellId,
  const int &localNeighborId,
  SimplexId &neighborId) const {

  if(dimensionality_ == 3)
    getTetrahedronNeighbor(cellId, localNeighborId, neighborId);
  else if(dimensionality_ == 2)
    getTriangleNeighbor(cellId, localNeighborId, neighborId);
  else if(dimensionality_ == 1) {
    printErr("getCellNeighbor() not implemented in 1D! (TODO)");
    return -1;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
    getCellNeighbors)() {

  if(cellNeighborList_.empty()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronNeighbors(cellNeighborList_);
    else if(dimensionality_ == 2)
      getTriangleNeighbors(cellNeighborList_);
    else if(dimensionality_ == 1) {
      printErr("getCellNeighbors() not implemented in 1D! (TODO)");
      return nullptr;
    }

    printMsg("Built " + to_string(cellNumber_) + " cell neighbors.", 1,
             t.getElapsedTime(), 1);
  }

  return &cellNeighborList_;
}

// explicit instantiations
template class ttk::PeriodicImplicitTriangulationCRTP<
  ttk::PeriodicWithPreconditions>;
template class ttk::PeriodicImplicitTriangulationCRTP<
  ttk::PeriodicNoPreconditions>;
