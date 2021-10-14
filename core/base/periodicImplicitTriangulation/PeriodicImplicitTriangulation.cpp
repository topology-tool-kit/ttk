#include <PeriodicImplicitTriangulation.h>

using namespace std;
using namespace ttk;

PeriodicImplicitTriangulation::PeriodicImplicitTriangulation()
  : dimensionality_{-1}, cellNumber_{}, vertexNumber_{}, edgeNumber_{},
    triangleNumber_{}, tetrahedronNumber_{}, isAccelerated_{} {
  setDebugMsgPrefix("PeriodicImplicitTriangulation");
  hasPeriodicBoundaries_ = true;
}

PeriodicImplicitTriangulation::~PeriodicImplicitTriangulation() {
}

int PeriodicImplicitTriangulation::setInputGrid(const float &xOrigin,
                                                const float &yOrigin,
                                                const float &zOrigin,
                                                const float &xSpacing,
                                                const float &ySpacing,
                                                const float &zSpacing,
                                                const int &xDim,
                                                const int &yDim,
                                                const int &zDim) {

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

  // ensure preconditionned vertices and cells
  this->preconditionVerticesInternal();
  this->preconditionCellsInternal();

  return 0;
}

int PeriodicImplicitTriangulation::preconditionVerticesInternal() {
  vertexCoords_.resize(vertexNumber_);

  if(dimensionality_ == 1) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexCoords_[i][0] = i;
    }
  } else if(dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      auto &p = vertexCoords_[i];
      vertexToPosition2d(i, p.data());
    }
  } else if(dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      auto &p = vertexCoords_[i];
      vertexToPosition(i, p.data());
    }
  }
  return 0;
}

int PeriodicImplicitTriangulation::preconditionEdgesInternal() {
  edgePositions_.resize(edgeNumber_);
  edgeCoords_.resize(edgeNumber_);

  if(dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      auto &p = edgeCoords_[i];
      if(i < esetshift_[0]) {
        edgeToPosition(i, 0, p.data());
        edgePositions_[i] = EdgePosition::L_3D;
      } else if(i < esetshift_[1]) {
        edgeToPosition(i, 1, p.data());
        edgePositions_[i] = EdgePosition::H_3D;
      } else if(i < esetshift_[2]) {
        edgeToPosition(i, 2, p.data());
        edgePositions_[i] = EdgePosition::P_3D;
      } else if(i < esetshift_[3]) {
        edgeToPosition(i, 3, p.data());
        edgePositions_[i] = EdgePosition::D1_3D;
      } else if(i < esetshift_[4]) {
        edgeToPosition(i, 4, p.data());
        edgePositions_[i] = EdgePosition::D2_3D;
      } else if(i < esetshift_[5]) {
        edgeToPosition(i, 5, p.data());
        edgePositions_[i] = EdgePosition::D3_3D;
      } else if(i < esetshift_[6]) {
        edgeToPosition(i, 6, p.data());
        edgePositions_[i] = EdgePosition::D4_3D;
      }
    }

  } else if(this->dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      auto &p = edgeCoords_[i];
      if(i < esetshift_[0]) {
        edgeToPosition2d(i, 0, p.data());
        edgePositions_[i] = EdgePosition::L_2D;
      } else if(i < esetshift_[1]) {
        edgeToPosition2d(i, 1, p.data());
        edgePositions_[i] = EdgePosition::H_2D;
      } else if(i < esetshift_[2]) {
        edgeToPosition2d(i, 2, p.data());
        edgePositions_[i] = EdgePosition::D1_2D;
      }
    }
  }

  else if(this->dimensionality_ == 1) {
    edgePositions_[0] = EdgePosition::FIRST_EDGE_1D;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 1; i < edgeNumber_ - 1; ++i) {
      edgePositions_[i] = EdgePosition::CENTER_1D;
    }
    edgePositions_[edgeNumber_ - 1] = EdgePosition::LAST_EDGE_1D;
  }

  edgeVertexAccelerated_.resize(edgeNumber_);

  if(isAccelerated_) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      const auto &p = edgeCoords_[i];
      edgeVertexAccelerated_[i] = (p[1] << div_[0]) + (p[2] << div_[1]);
    }
  } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      const auto &p = edgeCoords_[i];
      edgeVertexAccelerated_[i] = p[1] * vshift_[0] + p[2] * vshift_[1];
    }
  }

  return 0;
}

int PeriodicImplicitTriangulation::preconditionTrianglesInternal() {
  if(this->dimensionality_ != 3 && this->dimensionality_ != 2) {
    return 1;
  }

  trianglePositions_.resize(triangleNumber_);
  triangleCoords_.resize(triangleNumber_);

  if(dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      if(i < tsetshift_[0]) {
        triangleToPosition(i, 0, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::F_3D;
      } else if(i < tsetshift_[1]) {
        triangleToPosition(i, 1, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::H_3D;
      } else if(i < tsetshift_[2]) {
        triangleToPosition(i, 2, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::C_3D;
      } else if(i < tsetshift_[3]) {
        triangleToPosition(i, 3, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D1_3D;
      } else if(i < tsetshift_[4]) {
        triangleToPosition(i, 4, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D2_3D;
      } else if(i < tsetshift_[5]) {
        triangleToPosition(i, 5, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D3_3D;
      }
    }

  } else if(dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      triangleToPosition2d(i, triangleCoords_[i].data());
      if(i % 2 == 0) {
        trianglePositions_[i] = TrianglePosition::TOP_2D;
      } else {
        trianglePositions_[i] = TrianglePosition::BOTTOM_2D;
      }
    }
  }

  return 0;
}

int PeriodicImplicitTriangulation::preconditionTetrahedronsInternal() {
  if(this->dimensionality_ != 3) {
    return 1;
  }

  tetrahedronCoords_.resize(tetrahedronNumber_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    tetrahedronToPosition(i, tetrahedronCoords_[i].data());
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
    bool isDi = isPowerOfTwo(dimensions_[Di_], msb[Di_]);
    bool isDj = isPowerOfTwo(dimensions_[Dj_], msb[Dj_]);
    bool allDimensionsArePowerOfTwo = (isDi and isDj);

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

int PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getVertexNeighbor)(const SimplexId &vertexId,
                     const int &localNeighborId,
                     SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= getVertexNeighborNumber(vertexId))
    return -1;
#endif

  neighborId = -1;
  const auto &p = vertexCoords_[vertexId];

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
  if(!vertexNeighborList_.size()) {
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

int PeriodicImplicitTriangulation::getVertexEdgeInternal(
  const SimplexId &vertexId, const int &localEdgeId, SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localEdgeId < 0 or localEdgeId >= getVertexEdgeNumberInternal(vertexId))
    return -1;
#endif

  edgeId = -1;
  const auto &p = vertexCoords_[vertexId];

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
  if(!vertexEdgeList_.size()) {
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

int PeriodicImplicitTriangulation::getVertexTriangleInternal(
  const SimplexId &vertexId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= getVertexTriangleNumberInternal(vertexId))
    return -1;
#endif
  triangleId = -1;

  const auto &p = vertexCoords_[vertexId];

  if(dimensionality_ == 3) {
    triangleId = getVertexTriangle3d(p.data(), localTriangleId);
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getVertexTrianglesInternal() {
  if(!vertexTriangleList_.size()) {
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

int PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexLink)(
  const SimplexId &vertexId, const int &localLinkId, SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getVertexLinkNumber(vertexId))
    return -1;
#endif

  linkId = -1;
  const auto &p = vertexCoords_[vertexId];

  if(dimensionality_ == 3) {
    linkId = getVertexLink3d(p.data(), localLinkId);
  } else if(dimensionality_ == 2) {
    linkId = getVertexLink2d(p.data(), localLinkId); // abcd
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexLinks)() {
  if(!vertexLinkList_.size()) {
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

int PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexStar)(
  const SimplexId &vertexId, const int &localStarId, SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getVertexStarNumber(vertexId))
    return -1;
#endif

  starId = -1;
  const auto &p = vertexCoords_[vertexId];

  if(dimensionality_ == 3) {
    starId = getVertexStar3d(p.data(), localStarId);
  } else if(dimensionality_ == 2) {
    starId = getVertexStar2d(p.data(), localStarId);
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexStars)() {
  if(!vertexStarList_.size()) {
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

int PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexPoint)(
  const SimplexId &vertexId, float &x, float &y, float &z) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    const auto &p = vertexCoords_[vertexId];

    x = origin_[0] + spacing_[0] * p[0];
    y = origin_[1] + spacing_[1] * p[1];
    z = origin_[2] + spacing_[2] * p[2];
  } else if(dimensionality_ == 2) {
    const auto &p = vertexCoords_[vertexId];

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

int PeriodicImplicitTriangulation::getEdgeVertexInternal(
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
  const auto &p = edgeCoords_[edgeId];
  const SimplexId wrapXRight = (p[0] == nbvoxels_[0] ? -wrap_[0] : 0);
  const SimplexId wrapYBottom = (p[1] == nbvoxels_[1] ? -wrap_[1] : 0);
  const SimplexId wrapZFront = (p[2] == nbvoxels_[2] ? -wrap_[2] : 0);
  const auto a = p[0] + edgeVertexAccelerated_[edgeId];

  switch(edgePositions_[edgeId]) {
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

  if(!edgeList_.size()) {
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

SimplexId PeriodicImplicitTriangulation::getEdgeTriangleNumberInternal(
  const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  switch(edgePositions_[edgeId]) {
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

int PeriodicImplicitTriangulation::getEdgeTriangleInternal(
  const SimplexId &edgeId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= getEdgeTriangleNumberInternal(edgeId))
    return -1;
#endif

  triangleId = -1;
  const auto &p = edgeCoords_[edgeId];

  switch(edgePositions_[edgeId]) {
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
  if(!edgeTriangleList_.size()) {
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

int PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdgeLink)(
  const SimplexId &edgeId, const int &localLinkId, SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getEdgeLinkNumber(edgeId))
    return -1;
#endif

  linkId = -1;
  const auto &p = edgeCoords_[edgeId];

  switch(edgePositions_[edgeId]) {
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

  if(!edgeLinkList_.size()) {
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

SimplexId PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getEdgeStarNumber)(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  switch(edgePositions_[edgeId]) {
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

int PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdgeStar)(
  const SimplexId &edgeId, const int &localStarId, SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getEdgeStarNumber(edgeId))
    return -1;
#endif

  starId = -1;
  const auto &p = edgeCoords_[edgeId];

  switch(edgePositions_[edgeId]) {
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

  if(!edgeStarList_.size()) {
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

int PeriodicImplicitTriangulation::getTriangleVertexInternal(
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
  const auto &p = triangleCoords_[triangleId];
  const SimplexId wrapXRight = (p[0] / 2 == nbvoxels_[Di_]) ? -wrap_[0] : 0;
  const SimplexId wrapYBottom = (p[1] == nbvoxels_[Dj_]) ? -wrap_[1] : 0;

  switch(trianglePositions_[triangleId]) {
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

int PeriodicImplicitTriangulation::getTriangleEdgeInternal(
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
  const auto &p = triangleCoords_[triangleId];
  const SimplexId wrapXRight = (p[0] / 2 == nbvoxels_[Di_]) ? -wrap_[0] : 0;
  const SimplexId wrapYBottom = (p[1] == nbvoxels_[Dj_]) ? -wrap_[1] : 0;
  const SimplexId id = triangleId % 2;

  switch(trianglePositions_[triangleId]) {
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

  if(!triangleList_.size()) {
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

int PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
  const SimplexId &triangleId,
  const int &localLinkId,
  SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getTriangleLinkNumber(triangleId))
    return -1;
#endif

  linkId = -1;
  const auto &p = triangleCoords_[triangleId];

  switch(trianglePositions_[triangleId]) {
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

  if(!triangleLinkList_.size()) {
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

SimplexId PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
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

int PeriodicImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
  const SimplexId &triangleId,
  const int &localStarId,
  SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getTriangleStarNumber(triangleId))
    return -1;
#endif

  starId = -1;
  const auto &p = triangleCoords_[triangleId];

  switch(trianglePositions_[triangleId]) {
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

  if(!triangleStarList_.size()) {
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

int PeriodicImplicitTriangulation::getTriangleNeighbor(
  const SimplexId &triangleId,
  const int &localNeighborId,
  SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= getTriangleNeighborNumber(triangleId))
    return -1;
#endif

  neighborId = -1;
  const auto &p = triangleCoords_[triangleId];

  switch(trianglePositions_[triangleId]) {
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

int PeriodicImplicitTriangulation::getTetrahedronVertex(
  const SimplexId &tetId, const int &localVertexId, SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 4)
    return -2;
#endif

  vertexId = -1;

  if(dimensionality_ == 3) {
    const auto &p = tetrahedronCoords_[tetId];
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

int PeriodicImplicitTriangulation::getTetrahedronEdge(const SimplexId &tetId,
                                                      const int &localEdgeId,
                                                      SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localEdgeId < 0 or localEdgeId >= 6)
    return -2;
#endif

  edgeId = -1;

  if(dimensionality_ == 3) {
    const auto &p = tetrahedronCoords_[tetId];
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

int PeriodicImplicitTriangulation::getTetrahedronTriangle(
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
    const auto &p = tetrahedronCoords_[tetId];
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

int PeriodicImplicitTriangulation::getTetrahedronNeighbor(
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
    const auto &p = tetrahedronCoords_[tetId];
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

  if(!cellNeighborList_.size()) {
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
