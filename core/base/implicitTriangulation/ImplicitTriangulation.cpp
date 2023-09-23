#include <ImplicitTriangulation.h>

#include <numeric>

using namespace std;
using namespace ttk;

ImplicitTriangulation::ImplicitTriangulation()
  : cellNumber_{}, vertexNumber_{}, edgeNumber_{}, triangleNumber_{},
    tetrahedronNumber_{}, isAccelerated_{} {
  setDebugMsgPrefix("ImplicitTriangulation");
#ifdef TTK_ENABLE_MPI
  this->hasPreconditionedDistributedEdges_ = true;
  this->hasPreconditionedDistributedTriangles_ = true;
#endif // TTK_ENABLE_MPI
}

ImplicitTriangulation::~ImplicitTriangulation() = default;

int ImplicitTriangulation::setInputGrid(const float &xOrigin,
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
    esetdims_[0] = (xDim - 1) * yDim * zDim;
    esetdims_[1] = xDim * (yDim - 1) * zDim;
    esetdims_[2] = xDim * yDim * (zDim - 1);
    esetdims_[3] = (xDim - 1) * (yDim - 1) * zDim;
    esetdims_[4] = xDim * (yDim - 1) * (zDim - 1);
    esetdims_[5] = (xDim - 1) * yDim * (zDim - 1);
    esetdims_[6] = (xDim - 1) * (yDim - 1) * (zDim - 1);
    // EdgeSetShift
    esetshift_[0] = esetdims_[0];
    for(int k = 1; k < 7; ++k)
      esetshift_[k] = esetshift_[k - 1] + esetdims_[k];
    // EdgeShift
    eshift_[0] = xDim - 1;
    eshift_[1] = (xDim - 1) * yDim;
    eshift_[2] = xDim;
    eshift_[3] = xDim * (yDim - 1);
    eshift_[4] = xDim;
    eshift_[5] = xDim * yDim;
    eshift_[6] = xDim - 1;
    eshift_[7] = (xDim - 1) * (yDim - 1);
    eshift_[8] = xDim;
    eshift_[9] = xDim * (yDim - 1);
    eshift_[10] = xDim - 1;
    eshift_[11] = (xDim - 1) * yDim;
    eshift_[12] = xDim - 1;
    eshift_[13] = (xDim - 1) * (yDim - 1);
    // TriangleSetDimensions
    tsetdims_[0] = (xDim - 1) * (yDim - 1) * zDim * 2;
    tsetdims_[1] = (xDim - 1) * yDim * (zDim - 1) * 2;
    tsetdims_[2] = xDim * (yDim - 1) * (zDim - 1) * 2;
    tsetdims_[3] = (xDim - 1) * (yDim - 1) * (zDim - 1) * 2;
    tsetdims_[4] = (xDim - 1) * (yDim - 1) * (zDim - 1) * 2;
    tsetdims_[5] = (xDim - 1) * (yDim - 1) * (zDim - 1) * 2;
    // TriangleSetShift
    tsetshift_[0] = tsetdims_[0];
    for(int k = 1; k < 6; ++k)
      tsetshift_[k] = tsetshift_[k - 1] + tsetdims_[k];
    // TriangleShift
    tshift_[0] = (xDim - 1) * 2;
    tshift_[1] = (xDim - 1) * (yDim - 1) * 2;
    tshift_[2] = (xDim - 1) * 2;
    tshift_[3] = (xDim - 1) * yDim * 2;
    tshift_[4] = xDim * 2;
    tshift_[5] = xDim * (yDim - 1) * 2;
    tshift_[6] = (xDim - 1) * 2;
    tshift_[7] = (xDim - 1) * (yDim - 1) * 2;
    tshift_[8] = (xDim - 1) * 2;
    tshift_[9] = (xDim - 1) * (yDim - 1) * 2;
    tshift_[10] = (xDim - 1) * 2;
    tshift_[11] = (xDim - 1) * (yDim - 1) * 2;
    // TetrahedronShift
    tetshift_[0] = (xDim - 1) * 6;
    tetshift_[1] = (xDim - 1) * (yDim - 1) * 6;

    // Numbers
    vertexNumber_ = xDim * yDim * zDim;
    edgeNumber_ = 0;
    for(int k = 0; k < 7; ++k)
      edgeNumber_ += esetdims_[k];
    triangleNumber_ = 0;
    for(int k = 0; k < 6; ++k)
      triangleNumber_ += tsetdims_[k];
    tetrahedronNumber_ = (xDim - 1) * (yDim - 1) * (zDim - 1) * 6;
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
    esetdims_[0] = (dimensions_[Di_] - 1) * dimensions_[Dj_];
    esetdims_[1] = dimensions_[Di_] * (dimensions_[Dj_] - 1);
    esetdims_[2] = (dimensions_[Di_] - 1) * (dimensions_[Dj_] - 1);
    // EdgeSetShift
    esetshift_[0] = esetdims_[0];
    for(int k = 1; k < 3; ++k)
      esetshift_[k] = esetshift_[k - 1] + esetdims_[k];
    // EdgeShift
    eshift_[0] = dimensions_[Di_] - 1;
    eshift_[2] = dimensions_[Di_];
    eshift_[4] = dimensions_[Di_] - 1;
    // TriangleShift
    tshift_[0] = (dimensions_[Di_] - 1) * 2;

    // Numbers
    vertexNumber_ = dimensions_[Di_] * dimensions_[Dj_];
    edgeNumber_ = 0;
    for(int k = 0; k < 3; ++k)
      edgeNumber_ += esetdims_[k];
    triangleNumber_ = (dimensions_[Di_] - 1) * (dimensions_[Dj_] - 1) * 2;
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
    edgeNumber_ = vertexNumber_ - 1;
    cellNumber_ = edgeNumber_;
  }

  return 0;
}

int ImplicitTriangulation::checkAcceleration() {
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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() {
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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::getVertexEdgesInternal() {
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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::getVertexTrianglesInternal() {
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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexLinks)() {
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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexStars)() {

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

const vector<std::array<SimplexId, 2>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdges)() {

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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::getEdgeTrianglesInternal() {
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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() {

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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdgeStars)() {

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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::getTriangleEdgesInternal() {
  if(triangleEdgeVector_.empty()) {
    Timer t;

    getTriangleEdgesInternal(triangleEdgeVector_);

    printMsg("Built " + to_string(triangleNumber_) + " triangle edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &triangleEdgeVector_;
}

const vector<std::array<SimplexId, 3>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangles)() {

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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() {
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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangleStars)() {

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

const vector<vector<SimplexId>> *ImplicitTriangulation::getCellEdgesInternal() {
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

const vector<vector<SimplexId>> *
  ImplicitTriangulation::getCellTrianglesInternal() {
  if(cellTriangleVector_.empty()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronTriangles(cellTriangleVector_);

    printMsg("Built " + to_string(cellNumber_) + " cell triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &cellTriangleVector_;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() {
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

int ImplicitTriangulation::preconditionVertexNeighborsInternal() {
  // V(abcdefgh)=V(g)+V(d)::{g,h}+V(h)::{g}+V(b)::{c,d,g,h}
  this->vertexNeighborABCDEFGH_ = {
    -vshift_[0] - vshift_[1], // a
    1 - vshift_[0] - vshift_[1], // b
    -vshift_[1], // c
    1 - vshift_[1], // d
    -vshift_[0], // e
    1 - vshift_[0], // f
    1, // h
    -1 + vshift_[1], // V(d)::{g}
    vshift_[1], // V(d)::{h}
    -1, // V(h)::{g}
    -1 + vshift_[0], // V(b)::{c}
    vshift_[0], // V(b)::{d}
    -1 + vshift_[0] + vshift_[1], // V(b)::{g}
    vshift_[0] + vshift_[1] // V(b)::{h}
  };

  // V(abdc)=V(b)+V(d)::{b}+V(c)::{b}+V(a)::{b}
  this->vertexNeighborABCD_ = {
    -1, // a
    -1 + vshift_[0], // c
    vshift_[0], // d
    -1 + vshift_[1], // e
    vshift_[1], // f
    -1 + vshift_[0] + vshift_[1], // g
    vshift_[0] + vshift_[1], // h
    -vshift_[0], // V(d)::{b}
    1 - vshift_[0], // V(c)::{b}
    1, // V(a)::{b}
  };
  // V(efhg)=V(g)+V(h)::{g}+V(f)::{g,h}
  this->vertexNeighborEFGH_ = {
    -vshift_[0] - vshift_[1], // a
    1 - vshift_[0] - vshift_[1], // b
    -vshift_[1], // c
    1 - vshift_[1], // d
    -vshift_[0], // e
    1 - vshift_[0], // f
    1, // h
    -1, // V(h)::{g}
    -1 + vshift_[0], // V(f)::{g}
    vshift_[0], // V(f)::{h}
  };
  // V(aefb)=V(b)+V(a)::{b}+V(e)::{b}+V(f)::{b}
  this->vertexNeighborAEFB_ = {
    -1, // a
    -1 + vshift_[0], // c
    vshift_[0], // d
    -1 + vshift_[1], // e
    vshift_[1], // f
    -1 + vshift_[0] + vshift_[1], // g
    vshift_[0] + vshift_[1], // h
    1, // V(a)::{b}
    1 - vshift_[1], // V(e)::{b}
    -vshift_[1], // V(f)::{b}
  };
  // V(ghdc)=V(g)+V(h)::{g}+V(d)::{g,h}
  this->vertexNeighborGHDC_ = {
    -vshift_[0] - vshift_[1], // a
    1 - vshift_[0] - vshift_[1], // b
    -vshift_[1], // c
    1 - vshift_[1], // d
    -vshift_[0], // e
    1 - vshift_[0], // f
    1, // h
    -1, // V(h)::{g}
    -1 + vshift_[1], // V(d)::{g}
    vshift_[1] // V(d)::{h}
  };
  // V(aegc)=V(g)+V(a)::{c,g}+V(c)::{g}
  this->vertexNeighborAEGC_ = {
    -vshift_[0] - vshift_[1], // a
    1 - vshift_[0] - vshift_[1], // b
    -vshift_[1], // c
    1 - vshift_[1], // d
    -vshift_[0], // e
    1 - vshift_[0], // f
    1, // h
    vshift_[0], // V(a)::{c}
    vshift_[0] + vshift_[1], // V(a)::{g}
    vshift_[1], // V(c)::{g}
  };
  // V(bfhd)=V(b)+V(f)::{b}+V(h)::{b}+V(d)::{b}
  this->vertexNeighborBFHD_ = {
    -1, // a
    -1 + vshift_[0], // c
    vshift_[0], // d
    -1 + vshift_[1], // e
    vshift_[1], // f
    -1 + vshift_[0] + vshift_[1], // g
    vshift_[0] + vshift_[1], // h
    -vshift_[1], // V(f)::{b}
    -vshift_[0] - vshift_[1], // V(h)::{b}
    -vshift_[0], // V(d)::{b}
  };

  // V(ab)=V(b)+V(a)::{b}
  this->vertexNeighborAB_ = {
    -1, // a
    -1 + vshift_[0], // c
    vshift_[0], // d
    -1 + vshift_[1], // e
    vshift_[1], // f
    -1 + vshift_[0] + vshift_[1], // g
    vshift_[0] + vshift_[1], // h
    1, // V(a)::{b}
  };
  // V(bd)=V(b)+V(d)::{b}
  this->vertexNeighborBD_ = {
    -1, // a
    -1 + vshift_[0], // c
    vshift_[0], // d
    -1 + vshift_[1], // e
    vshift_[1], // f
    -1 + vshift_[0] + vshift_[1], // g
    vshift_[0] + vshift_[1], // h
    -vshift_[0], // V(d)::{b}
  };
  // V(gh)=V(g)+V(h)::{g}
  this->vertexNeighborGH_ = {
    -vshift_[0] - vshift_[1], // a
    1 - vshift_[0] - vshift_[1], // b
    -vshift_[1], // c
    1 - vshift_[1], // d
    -vshift_[0], // e
    1 - vshift_[0], // f
    1, // h
    -1, // V(h)::{g}
  };
  // V(eg)=V(g)+V(e)::{g}
  this->vertexNeighborEG_ = {
    -vshift_[0] - vshift_[1], // a
    1 - vshift_[0] - vshift_[1], // b
    -vshift_[1], // c
    1 - vshift_[1], // d
    -vshift_[0], // e
    1 - vshift_[0], // f
    1, // h
    vshift_[0], // V(e)::{g}
  };
  // V(cg)=V(g)+V(c)::{g}
  this->vertexNeighborCG_ = {
    -vshift_[0] - vshift_[1], // a
    1 - vshift_[0] - vshift_[1], // b
    -vshift_[1], // c
    1 - vshift_[1], // d
    -vshift_[0], // e
    1 - vshift_[0], // f
    1, // h
    vshift_[1], // V(c)::{g}
  };
  // V(bf)=V(b)+V(f)::{b}
  this->vertexNeighborBF_ = {
    -1, // a
    -1 + vshift_[0], // c
    vshift_[0], // d
    -1 + vshift_[1], // e
    vshift_[1], // f
    -1 + vshift_[0] + vshift_[1], // g
    vshift_[0] + vshift_[1], // h
    -vshift_[1], // V(f)::{b}
  };

  // V(b)={a,c,d,e,f,g,h}
  this->vertexNeighborB_ = {
    -1, // a
    -1 + vshift_[0], // c
    vshift_[0], // d
    -1 + vshift_[1], // e
    vshift_[1], // f
    -1 + vshift_[0] + vshift_[1], // g
    vshift_[0] + vshift_[1], // h
  };
  // V(g)={a,b,c,d,e,f,h}
  this->vertexNeighborG_ = {
    -vshift_[0] - vshift_[1], // a
    1 - vshift_[0] - vshift_[1], // b
    -vshift_[1], // c
    1 - vshift_[1], // d
    -vshift_[0], // e
    1 - vshift_[0], // f
    1, // h
  };

  // V(ef)=V(f)+V(e)::{b,f}
  this->vertexNeighborEF_ = {
    -vshift_[1], // b
    -1, // e
    -1 + vshift_[0], // g
    vshift_[0], // h
    1 - vshift_[1], // V(e)::{b}
    1, // V(e)::{f}
  };
  // V(cd)=V(d)+V(c)::{b,d}
  this->vertexNeighborCD_ = {
    -vshift_[0], // b
    -1, // c
    -1 + vshift_[1], // g
    vshift_[1], // h
    1 - vshift_[0], // V(c)::{b}
    1, // V(c)::{d}
  };
  // V(ac)=V(c)+V(a)::{c,g}
  this->vertexNeighborAC_ = {
    -vshift_[0], // a
    1 - vshift_[0], // b
    1, // d
    vshift_[1], // g
    vshift_[0], // V(a)::{c}
    vshift_[0] + vshift_[1], // V(a)::{c}
  };
  // V(ae)=V(a)+V(e)::{a,b}
  this->vertexNeighborAE_ = {
    1, // b
    vshift_[0], // c
    vshift_[1], // e
    vshift_[0] + vshift_[1], // g
    -vshift_[1], // V(e)::{a}
    1 - vshift_[1], // V(e)::{b}
  };
  // V(fh)=V(f)+V(h)::{b,f}
  this->vertexNeighborFH_ = {
    -vshift_[1], // b
    -1, // e
    -1 + vshift_[0], // g
    vshift_[0], // h
    -vshift_[0] - vshift_[1], // V(h)::{b}
    -vshift_[0], // V(h)::{f}
  };
  // V(dh)=V(d)+V(h)::{b,d}
  this->vertexNeighborDH_ = {
    -vshift_[0], // b
    -1, // c
    -1 + vshift_[1], // g
    vshift_[1], // h
    -vshift_[0] - vshift_[1], // V(h)::{b}
    -vshift_[1], // V(h)::{d}
  };

  // V(a)={b,c,e,g}
  this->vertexNeighborA_ = {
    1, // b
    vshift_[0], // c
    vshift_[1], // e
    vshift_[0] + vshift_[1], // g
  };
  // V(c)={a,b,d,g}
  this->vertexNeighborC_ = {
    -vshift_[0], // a
    1 - vshift_[0], // b
    1, // d
    +vshift_[1], // g
  };
  // V(d)={b,c,g,h}
  this->vertexNeighborD_ = {
    -vshift_[0], // b
    -1, // c
    -1 + vshift_[1], // g
    vshift_[1], // h
  };
  // V(e)={a,b,f,g}
  this->vertexNeighborE_ = {
    -vshift_[1], // a
    1 - vshift_[1], // b
    1, // f
    +vshift_[0], // g
  };
  // V(f)={b,e,g,h}
  this->vertexNeighborF_ = {
    -vshift_[1], // b
    -1, // e
    -1 + vshift_[0], // g
    vshift_[0], // h
  };
  // V(h)={b,d,f,g}
  this->vertexNeighborH_ = {
    -vshift_[0] - vshift_[1], // b
    -vshift_[1], // d
    -vshift_[0], // f
    -1, // g
  };

  // V(abcd)=V(d)::{b,c}+V(c)::{b,d}+V(a)::{b}+V(b)::{c}
  this->vertexNeighbor2dABCD_ = {
    -1, -vshift_[0], -vshift_[0] + 1, 1, vshift_[0], vshift_[0] - 1,
  };
  // V(ab)=V(b)::{a,c,d}+V(a)::{b}
  this->vertexNeighbor2dAB_ = {
    -1, // V(b)::a
    vshift_[0] - 1, // V(b)::c
    vshift_[0], // V(b)::d
    +1, // V(a)::b
  };
  // V(cd)=V(c)::{a,b,d}+V(d)::{c}
  this->vertexNeighbor2dCD_ = {
    -1, // V(d)::c
    -vshift_[0], // V(c)::a
    -vshift_[0] + 1, // V(c)::b
    1, // V(c)::d
  };
  // V(ac)=V(c)::{a,b,d}+V(a)::{c}
  this->vertexNeighbor2dAC_ = {
    -vshift_[0], // V(c)::{a}
    -vshift_[0] + 1, // V(c)::{b}
    1, // V(c)::{d}
    vshift_[0], // V(a)::{c}
  };
  // V(bd)=V(b)::{c,d}+V(d)::{b,c}
  this->vertexNeighbor2dBD_ = {
    vshift_[0] - 1, // V(b)::{c}
    vshift_[0], // V(b)::{d}
    -vshift_[0], // V(d)::{b}
    -1, // V(d)::{c}
  };
  // V(b)={a,c,d}
  this->vertexNeighbor2dB_ = {
    -1, // a
    vshift_[0], // d
    vshift_[0] - 1, // c
  };
  // V(c)={a,b,d}
  this->vertexNeighbor2dC_ = {
    1, // d
    -vshift_[0], // a
    -vshift_[0] + 1, // b
  };
  this->vertexNeighbor2dA_ = {};
  // V(a)={b,c}
  this->vertexNeighbor2dA_ = {
    1, // b
    vshift_[0] // c
  };
  // V(d)={c,b}
  this->vertexNeighbor2dD_ = {
    -1, // c
    -vshift_[0], // b

  };

  return 0;
}

#ifdef TTK_ENABLE_MPI

int ttk::ImplicitTriangulation::preconditionDistributedCells() {
  if(this->hasPreconditionedDistributedCells_) {
    return 0;
  }
  if(!ttk::hasInitializedMPI()) {
    return -1;
  }
  if(this->cellGhost_ == nullptr) {
    if(ttk::isRunningWithMPI()) {
      this->printErr("Missing cell ghost array!");
    }
    return -3;
  }

  Timer tm{};

  // number of local cells (with ghost cells...)
  const auto nLocCells{this->getNumberOfCells()};

  // there are 6 tetrahedra per cubic cell (and 2 triangles per square)
  const int nTetraPerCube{this->dimensionality_ == 3 ? 6 : 2};
  std::vector<unsigned char> fillCells(nLocCells / nTetraPerCube);

  this->ghostCellsPerOwner_.resize(ttk::MPIsize_);

  this->neighborCellBBoxes_.resize(ttk::MPIsize_);
  auto &localBBox{this->neighborCellBBoxes_[ttk::MPIrank_]};

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
  for(SimplexId lcid = 0; lcid < nLocCells; ++lcid) {
    // only keep non-ghost cells
    if(this->cellGhost_[lcid / nTetraPerCube] == 1) {
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

void ttk::ImplicitTriangulation::createMetaGrid(const double *const bounds) {

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

  const std::array<int, 3> dimensions = {
    static_cast<int>(
      std::round((globalBounds[1] - globalBounds[0]) / this->spacing_[0]))
      + 1,
    static_cast<int>(
      std::round((globalBounds[3] - globalBounds[2]) / this->spacing_[1]))
      + 1,
    static_cast<int>(
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

  this->metaGrid_ = std::make_shared<ImplicitNoPreconditions>();
  this->metaGrid_->setInputGrid(globalBounds[0], globalBounds[1],
                                globalBounds[2], this->spacing_[0],
                                this->spacing_[1], this->spacing_[2],
                                dimensions[0], dimensions[1], dimensions[2]);
  this->metaGrid_->preconditionBoundaryVertices();
  this->metaGrid_->preconditionBoundaryEdges();
  this->metaGrid_->preconditionBoundaryTriangles();
}

#endif // TTK_ENABLE_MPI

// explicit instantiations
template class ttk::ImplicitTriangulationCRTP<ttk::ImplicitWithPreconditions>;
template class ttk::ImplicitTriangulationCRTP<ttk::ImplicitNoPreconditions>;
