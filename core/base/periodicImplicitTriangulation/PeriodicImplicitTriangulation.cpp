#include <PeriodicImplicitTriangulation.h>

using namespace std;
using namespace ttk;

PeriodicImplicitTriangulation::PeriodicImplicitTriangulation()
  : dimensionality_{-1}, cellNumber_{}, vertexNumber_{}, edgeNumber_{},
    triangleNumber_{}, tetrahedronNumber_{}, isAccelerated_{} {
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
    stringstream msg;
    msg << "[PeriodicImplicitTriangulation] The getVertex*() requests are "
           "accelerated."
        << endl;
    dMsg(cout, msg.str(), infoMsg);
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

bool PeriodicImplicitTriangulation::isVertexOnBoundary(
  const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  return false;
}

bool PeriodicImplicitTriangulation::isEdgeOnBoundary(
  const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  return false;
}

bool PeriodicImplicitTriangulation::isTriangleOnBoundary(
  const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  return false;
}

inline SimplexId PeriodicImplicitTriangulation::getVertexNeighborNumber(
  const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
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

int PeriodicImplicitTriangulation::getVertexNeighbor(
  const SimplexId &vertexId,
  const int &localNeighborId,
  SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= getVertexNeighborNumber(vertexId))
    return -1;
#endif

  neighborId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];
    vertexToPosition(vertexId, p);
    neighborId = getVertexNeighbor3d(p, vertexId, localNeighborId);
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);
    neighborId = getVertexNeighbor2d(p, vertexId, localNeighborId);
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
  PeriodicImplicitTriangulation::getVertexNeighbors() {
  if(!vertexNeighborList_.size()) {
    Timer t;
    vertexNeighborList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexNeighborList_[i].resize(getVertexNeighborNumber(i));
      for(SimplexId j = 0; j < (SimplexId)vertexNeighborList_[i].size(); ++j)
        getVertexNeighbor(i, j, vertexNeighborList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Vertex neighbors built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &vertexNeighborList_;
}

SimplexId PeriodicImplicitTriangulation::getVertexEdgeNumber(
  const SimplexId &vertexId) const {
  return getVertexNeighborNumber(vertexId);
}

int PeriodicImplicitTriangulation::getVertexEdge(const SimplexId &vertexId,
                                                 const int &localEdgeId,
                                                 SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localEdgeId < 0 or localEdgeId >= getVertexEdgeNumber(vertexId))
    return -1;
#endif
  //    e--------f
  //   /|       /|
  //  / |      / |
  // a--g-----b--h
  // | /      | /
  // |/       |/
  // c--------d
  //
  // Classement des "Edges" et dans cet ordre:
  // L: largeur (type ab)
  // H: hauteur (type ac)
  // P: profondeur (type ae)
  // D1: diagonale1 (type bc)
  // D2: diagonale2 (type ag)
  // D3: diagonale3 (type be)
  // D4: diagonale4 (type bg)

  edgeId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];
    vertexToPosition(vertexId, p);
    edgeId = getVertexEdge3d(p, localEdgeId);
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);
    edgeId = getVertexEdge2d(p, localEdgeId);
  } else if(dimensionality_ == 1) {
    // ab
    if(vertexId > 0 and vertexId < nbvoxels_[Di_]) {
      if(localEdgeId == 0)
        edgeId = vertexId;
      else
        edgeId = vertexId - 1;
    } else if(vertexId == 0) {
      if(localEdgeId == 0)
        edgeId = vertexId; // a
      else
        edgeId = 0;
    } else {
      if(localEdgeId == 0)
        edgeId = 0;
      else
        edgeId = vertexId - 1; // b
    }
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getVertexEdges() {
  if(!vertexEdgeList_.size()) {
    Timer t;

    vertexEdgeList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexEdgeList_[i].resize(getVertexEdgeNumber(i));
      for(SimplexId j = 0; j < (SimplexId)vertexEdgeList_[i].size(); ++j)
        getVertexEdge(i, j, vertexEdgeList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Vertex edges built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &vertexEdgeList_;
}

inline SimplexId PeriodicImplicitTriangulation::getVertexTriangleNumber(
  const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    return 36; // abcdefgh
  }

  return 0;
}

int PeriodicImplicitTriangulation::getVertexTriangle(
  const SimplexId &vertexId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= getVertexTriangleNumber(vertexId))
    return -1;
#endif
  triangleId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];
    vertexToPosition(vertexId, p);
    triangleId = getVertexTriangle3d(p, localTriangleId);
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getVertexTriangles() {
  if(!vertexTriangleList_.size()) {
    Timer t;

    vertexTriangleList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexTriangleList_[i].resize(getVertexTriangleNumber(i));
      for(SimplexId j = 0; j < (SimplexId)vertexTriangleList_[i].size(); ++j)
        getVertexTriangle(i, j, vertexTriangleList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Vertex triangles built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &vertexTriangleList_;
}

SimplexId PeriodicImplicitTriangulation::getVertexLinkNumber(
  const SimplexId &vertexId) const {
  return getVertexStarNumber(vertexId);
}

int PeriodicImplicitTriangulation::getVertexLink(const SimplexId &vertexId,
                                                 const int &localLinkId,
                                                 SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getVertexLinkNumber(vertexId))
    return -1;
#endif

  linkId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];
    vertexToPosition(vertexId, p);
    linkId = getVertexLink3d(p, localLinkId);
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);
    linkId = getVertexLink2d(p, localLinkId); // abcd
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getVertexLinks() {
  if(!vertexLinkList_.size()) {
    Timer t;

    vertexLinkList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexLinkList_[i].resize(getVertexLinkNumber(i));
      for(SimplexId j = 0; j < (SimplexId)vertexLinkList_[i].size(); ++j)
        getVertexLink(i, j, vertexLinkList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Vertex links built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &vertexLinkList_;
}

inline SimplexId PeriodicImplicitTriangulation::getVertexStarNumber(
  const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    return 24; // abcdefgh
  } else if(dimensionality_ == 2) {
    return 6; // abcd
  }

  return 0;
}

int PeriodicImplicitTriangulation::getVertexStar(const SimplexId &vertexId,
                                                 const int &localStarId,
                                                 SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getVertexStarNumber(vertexId))
    return -1;
#endif

  starId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];
    vertexToPosition(vertexId, p);
    starId = getVertexStar3d(p, localStarId);
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);
    starId = getVertexStar2d(p, localStarId);
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getVertexStars() {
  if(!vertexStarList_.size()) {
    Timer t;
    vertexStarList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      vertexStarList_[i].resize(getVertexStarNumber(i));
      for(SimplexId j = 0; j < (SimplexId)vertexStarList_[i].size(); ++j)
        getVertexStar(i, j, vertexStarList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Vertex stars built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &vertexStarList_;
}

int PeriodicImplicitTriangulation::getVertexPoint(const SimplexId &vertexId,
                                                  float &x,
                                                  float &y,
                                                  float &z) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    SimplexId p[3];
    vertexToPosition(vertexId, p);

    x = origin_[0] + spacing_[0] * p[0];
    y = origin_[1] + spacing_[1] * p[1];
    z = origin_[2] + spacing_[2] * p[2];
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);

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

int PeriodicImplicitTriangulation::getEdgeVertex(const SimplexId &edgeId,
                                                 const int &localVertexId,
                                                 SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 2)
    return -2;
#endif

  vertexId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];
    SimplexId wrapXRight = 0;
    SimplexId wrapYBottom = 0;
    SimplexId wrapZFront = 0;
    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition(edgeId, 0, p);
      if(p[0] == nbvoxels_[0])
        wrapXRight = -wrap_[0];
      if(p[1] == nbvoxels_[1])
        wrapYBottom = -wrap_[1];
      if(p[2] == nbvoxels_[2])
        wrapZFront = -wrap_[2];
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]);
        else
          vertexId
            = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + 1 + wrapXRight;
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
        else
          vertexId
            = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
      }
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition(edgeId, 1, p);
      if(p[0] == nbvoxels_[0])
        wrapXRight = -wrap_[0];
      if(p[1] == nbvoxels_[1])
        wrapYBottom = -wrap_[1];
      if(p[2] == nbvoxels_[2])
        wrapZFront = -wrap_[2];
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[0]
                     + wrapYBottom;
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
                     + wrapYBottom;
      }
    }
    // P
    else if(edgeId < esetshift_[2]) {
      edgeToPosition(edgeId, 2, p);
      if(p[0] == nbvoxels_[0])
        wrapXRight = -wrap_[0];
      if(p[1] == nbvoxels_[1])
        wrapYBottom = -wrap_[1];
      if(p[2] == nbvoxels_[2])
        wrapZFront = -wrap_[2];
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[1]
                     + wrapZFront;
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
                     + wrapZFront;
      }
    }
    // D1
    else if(edgeId < esetshift_[3]) {
      edgeToPosition(edgeId, 3, p);
      if(p[0] == nbvoxels_[0])
        wrapXRight = -wrap_[0];
      if(p[1] == nbvoxels_[1])
        wrapYBottom = -wrap_[1];
      if(p[2] == nbvoxels_[2])
        wrapZFront = -wrap_[2];
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId
            = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + 1 + wrapXRight;
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[0]
                     + wrapYBottom;
      } else {
        if(localVertexId == 0)
          vertexId
            = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
                     + wrapYBottom;
      }
    }
    // D2
    else if(edgeId < esetshift_[4]) {
      edgeToPosition(edgeId, 4, p);
      if(p[0] == nbvoxels_[0])
        wrapXRight = -wrap_[0];
      if(p[1] == nbvoxels_[1])
        wrapYBottom = -wrap_[1];
      if(p[2] == nbvoxels_[2])
        wrapZFront = -wrap_[2];
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[0]
                     + vshift_[1] + wrapYBottom + wrapZFront;
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
                     + vshift_[1] + wrapYBottom + wrapZFront;
      }
    }
    // D3
    else if(edgeId < esetshift_[5]) {
      edgeToPosition(edgeId, 5, p);
      if(p[0] == nbvoxels_[0])
        wrapXRight = -wrap_[0];
      if(p[1] == nbvoxels_[1])
        wrapYBottom = -wrap_[1];
      if(p[2] == nbvoxels_[2])
        wrapZFront = -wrap_[2];
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId
            = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + 1 + wrapXRight;
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[1]
                     + wrapZFront;

      } else {
        if(localVertexId == 0)
          vertexId
            = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
                     + wrapZFront;
      }
    }
    // D4
    else if(edgeId < esetshift_[6]) {
      edgeToPosition(edgeId, 6, p);
      if(p[0] == nbvoxels_[0])
        wrapXRight = -wrap_[0];
      if(p[1] == nbvoxels_[1])
        wrapYBottom = -wrap_[1];
      if(p[2] == nbvoxels_[2])
        wrapZFront = -wrap_[2];
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId
            = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + 1 + wrapXRight;
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[0]
                     + vshift_[1] + wrapYBottom + wrapZFront;

      } else {
        if(localVertexId == 0)
          vertexId
            = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
                     + vshift_[1] + wrapYBottom + wrapZFront;
      }
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    SimplexId wrapXRight = 0;
    SimplexId wrapYBottom = 0;
    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition2d(edgeId, 0, p);
      if(p[0] == nbvoxels_[Di_])
        wrapXRight = -wrap_[0];
      if(p[1] == nbvoxels_[Dj_])
        wrapYBottom = -wrap_[1];
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + 1 + wrapXRight;
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0];
        else
          vertexId = p[0] + p[1] * vshift_[0] + 1 + wrapXRight;
      }
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition2d(edgeId, 1, p);
      if(p[0] == nbvoxels_[Di_])
        wrapXRight = -wrap_[0];
      if(p[1] == nbvoxels_[Dj_])
        wrapYBottom = -wrap_[1];
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + vshift_[0] + wrapYBottom;
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0];
        else
          vertexId = p[0] + p[1] * vshift_[0] + vshift_[0] + wrapYBottom;
      }
    }
    // D1
    else if(edgeId < esetshift_[2]) {
      edgeToPosition2d(edgeId, 2, p);
      if(p[0] == nbvoxels_[Di_])
        wrapXRight = -wrap_[0];
      if(p[1] == nbvoxels_[Dj_])
        wrapYBottom = -wrap_[1];
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + 1 + wrapXRight;
        else
          vertexId = p[0] + (p[1] << div_[0]) + vshift_[0] + wrapYBottom;
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + 1 + wrapXRight;
        else
          vertexId = p[0] + p[1] * vshift_[0] + vshift_[0] + wrapYBottom;
      }
    }
  } else if(dimensionality_ == 1) {
    if(edgeId > 0 and edgeId < (edgeNumber_)) {
      if(localVertexId == 0)
        vertexId = edgeId;
      else
        vertexId = edgeId + 1;
    } else if(edgeId == 0) {
      if(localVertexId == 0)
        vertexId = 0;
      else
        vertexId = 1;
    } else {
      if(localVertexId == 0)
        vertexId = edgeId;
      else
        vertexId = 0;
    }
  }

  return 0;
}

const vector<pair<SimplexId, SimplexId>> *
  PeriodicImplicitTriangulation::getEdges() {
  if(!edgeList_.size()) {
    Timer t;

    edgeList_.resize(edgeNumber_);
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      SimplexId id0, id1;
      getEdgeVertex(i, 0, id0);
      getEdgeVertex(i, 1, id1);
      edgeList_[i].first = id0;
      edgeList_[i].second = id1;
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Edge-list built in "
          << t.getElapsedTime() << " s. (" << edgeList_.size() << " edges, ("
          << 1 << " thread(s))" << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &edgeList_;
}

inline SimplexId PeriodicImplicitTriangulation::getEdgeTriangleNumber(
  const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    // L
    if(edgeId < esetshift_[0]) {
      return 6;
    }
    // H
    else if(edgeId < esetshift_[1]) {
      return 6;
    }
    // P
    else if(edgeId < esetshift_[2]) {
      return 6;
    }
    // D1
    else if(edgeId < esetshift_[3]) {
      return 4;
    }
    // D2
    else if(edgeId < esetshift_[4]) {
      return 4;
    }
    // D3
    else if(edgeId < esetshift_[5]) {
      return 4;
    }
    // D4
    else if(edgeId < esetshift_[6])
      return 6;
  } else if(dimensionality_ == 2) {
    return 2;
  }

  return 0;
}

int PeriodicImplicitTriangulation::getEdgeTriangle(
  const SimplexId &edgeId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0 or localTriangleId >= getEdgeTriangleNumber(edgeId))
    return -1;
#endif

  triangleId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition(edgeId, 0, p);
      triangleId = getEdgeTriangle3dL(p, localTriangleId);
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition(edgeId, 1, p);
      triangleId = getEdgeTriangle3dH(p, localTriangleId);
    }
    // P
    else if(edgeId < esetshift_[2]) {
      edgeToPosition(edgeId, 2, p);
      triangleId = getEdgeTriangle3dP(p, localTriangleId);
    }
    // D1
    else if(edgeId < esetshift_[3]) {
      edgeToPosition(edgeId, 3, p);
      triangleId = getEdgeTriangle3dD1(p, localTriangleId);
    }
    // D2
    else if(edgeId < esetshift_[4]) {
      edgeToPosition(edgeId, 4, p);
      triangleId = getEdgeTriangle3dD2(p, localTriangleId);
    }
    // D3
    else if(edgeId < esetshift_[5]) {
      edgeToPosition(edgeId, 5, p);
      triangleId = getEdgeTriangle3dD3(p, localTriangleId);
    }
    // D4
    else if(edgeId < esetshift_[6]) {
      edgeToPosition(edgeId, 6, p);
      triangleId = getEdgeTriangle3dD4(p, localTriangleId);
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition2d(edgeId, 0, p);
      triangleId = getEdgeTriangle2dL(p, localTriangleId);
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition2d(edgeId, 1, p);
      triangleId = getEdgeTriangle2dH(p, localTriangleId);
    }
    // D1
    else if(edgeId < esetshift_[2]) {
      edgeToPosition2d(edgeId, 2, p);
      triangleId = getEdgeTriangle2dD1(p, localTriangleId);
    }
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getEdgeTriangles() {
  if(!edgeTriangleList_.size()) {
    Timer t;

    edgeTriangleList_.resize(edgeNumber_);
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      edgeTriangleList_[i].resize(getEdgeTriangleNumber(i));
      for(SimplexId j = 0; j < (SimplexId)edgeTriangleList_[i].size(); ++j)
        getEdgeTriangle(i, j, edgeTriangleList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Triangle edges built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &edgeTriangleList_;
}

inline SimplexId PeriodicImplicitTriangulation::getEdgeLinkNumber(
  const SimplexId &edgeId) const {
  return getEdgeStarNumber(edgeId);
}

int PeriodicImplicitTriangulation::getEdgeLink(const SimplexId &edgeId,
                                               const int &localLinkId,
                                               SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getEdgeLinkNumber(edgeId))
    return -1;
#endif

  linkId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];

    if(edgeId < esetshift_[0]) {
      edgeToPosition(edgeId, 0, p);
      linkId = getEdgeLinkL(p, localLinkId); // L
    } else if(edgeId < esetshift_[1]) {
      edgeToPosition(edgeId, 1, p);
      linkId = getEdgeLinkH(p, localLinkId); // H
    } else if(edgeId < esetshift_[2]) {
      edgeToPosition(edgeId, 2, p);
      linkId = getEdgeLinkP(p, localLinkId); // P
    } else if(edgeId < esetshift_[3]) {
      edgeToPosition(edgeId, 3, p);
      linkId = getEdgeLinkD1(p, localLinkId); // D1
    } else if(edgeId < esetshift_[4]) {
      edgeToPosition(edgeId, 4, p);
      linkId = getEdgeLinkD2(p, localLinkId); // D2
    } else if(edgeId < esetshift_[5]) {
      edgeToPosition(edgeId, 5, p);
      linkId = getEdgeLinkD3(p, localLinkId); // D3
    } else if(edgeId < esetshift_[6]) {
      edgeToPosition(edgeId, 6, p);
      linkId = getEdgeLinkD4(p, localLinkId); // D4
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition2d(edgeId, 0, p);
      linkId = getEdgeLink2dL(p, localLinkId);
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition2d(edgeId, 1, p);
      linkId = getEdgeLink2dH(p, localLinkId);
    }
    // D1
    else if(edgeId < esetshift_[2]) {
      edgeToPosition2d(edgeId, 2, p);
      linkId = getEdgeLink2dD1(p, localLinkId);
    }
  }

  return 0;
}

const vector<vector<SimplexId>> *PeriodicImplicitTriangulation::getEdgeLinks() {
  if(!edgeLinkList_.size()) {
    Timer t;

    edgeLinkList_.resize(edgeNumber_);
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      edgeLinkList_[i].resize(getEdgeLinkNumber(i));
      for(SimplexId j = 0; j < (SimplexId)edgeLinkList_[i].size(); ++j)
        getEdgeLink(i, j, edgeLinkList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] List of edge links built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &edgeLinkList_;
}

inline SimplexId PeriodicImplicitTriangulation::getEdgeStarNumber(
  const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    // L
    if(edgeId < esetshift_[0]) {
      return 6; // ABCG,ABEG,BCDG,BEFG,BFGH,BDGH
    }
    // H
    else if(edgeId < esetshift_[1]) {
      return 6; // ABCG,ABEG,BEFG,BFGH,BCDG,BDGH
    }
    // P
    else if(edgeId < esetshift_[2]) {
      return 6; // BDGH,ABCG,BCDG,ABEG,BEFG,BFGH
    }
    // D1
    else if(edgeId < esetshift_[3]) {
      return 4; // ABCG,BCDG,BEFG,BFGH
    }
    // D2
    else if(edgeId < esetshift_[4]) {
      return 4; // ABCG,ABEG,BDGH,BFGH
    }
    // D3
    else if(edgeId < esetshift_[5]) {
      return 4; // ABEG,BEFG,BCDG,BDGH
    }
    // D4
    else if(edgeId < esetshift_[6])
      return 6;
  } else if(dimensionality_ == 2) {
    return 2;
  }

  return 0;
}

int PeriodicImplicitTriangulation::getEdgeStar(const SimplexId &edgeId,
                                               const int &localStarId,
                                               SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getEdgeStarNumber(edgeId))
    return -1;
#endif

  starId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];

    if(edgeId < esetshift_[0]) {
      edgeToPosition(edgeId, 0, p);
      starId = getEdgeStarL(p, localStarId); // L
    } else if(edgeId < esetshift_[1]) {
      edgeToPosition(edgeId, 1, p);
      starId = getEdgeStarH(p, localStarId); // H
    } else if(edgeId < esetshift_[2]) {
      edgeToPosition(edgeId, 2, p);
      starId = getEdgeStarP(p, localStarId); // P
    } else if(edgeId < esetshift_[3]) {
      edgeToPosition(edgeId, 3, p);
      starId = getEdgeStarD1(p, localStarId); // D1
    } else if(edgeId < esetshift_[4]) {
      edgeToPosition(edgeId, 4, p);
      starId = getEdgeStarD2(p, localStarId); // D2
    } else if(edgeId < esetshift_[5]) {
      edgeToPosition(edgeId, 5, p);
      starId = getEdgeStarD3(p, localStarId); // D3
    } else if(edgeId < esetshift_[6]) {
      edgeToPosition(edgeId, 6, p);
      starId = p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + localStarId; // D4
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition2d(edgeId, 0, p);
      starId = getEdgeStar2dL(p, localStarId); // L
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition2d(edgeId, 1, p);
      starId = getEdgeStar2dH(p, localStarId); // L
    }
    // D1
    else if(edgeId < esetshift_[2]) {
      edgeToPosition2d(edgeId, 2, p);
      starId = p[0] * 2 + p[1] * tshift_[0] + localStarId; // D1
    }
  }

  return 0;
}

const vector<vector<SimplexId>> *PeriodicImplicitTriangulation::getEdgeStars() {
  if(!edgeStarList_.size()) {
    Timer t;

    edgeStarList_.resize(edgeNumber_);
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      edgeStarList_[i].resize(getEdgeStarNumber(i));
      for(SimplexId j = 0; j < (SimplexId)edgeStarList_[i].size(); ++j)
        getEdgeStar(i, j, edgeStarList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] List of edge stars built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &edgeStarList_;
}

int PeriodicImplicitTriangulation::getTriangleVertex(
  const SimplexId &triangleId,
  const int &localVertexId,
  SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 3)
    return -2;
#endif

  //    e--------f
  //   /|       /|
  //  / |      / |
  // a--g-----b--h
  // | /      | /
  // |/       |/
  // c--------d
  //
  // Classement des "Triangles" et dans cet ordre:
  // F: face (type abc/bcd)
  // C: cote (type abe/bef)
  // H: haut (type acg/aeg)
  // D1: diagonale1 (type bdg/beg)
  // D2: diagonale2 (type abg/bgh)
  // D3: diagonale3 (type bcg/bfg)

  vertexId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];

    // F
    if(triangleId < tsetshift_[0]) {
      triangleToPosition(triangleId, 0, p);
      vertexId = getTriangleVertexF(p, localVertexId);
    }
    // H
    else if(triangleId < tsetshift_[1]) {
      triangleToPosition(triangleId, 1, p);
      vertexId = getTriangleVertexH(p, localVertexId);
    }
    // C
    else if(triangleId < tsetshift_[2]) {
      triangleToPosition(triangleId, 2, p);
      vertexId = getTriangleVertexC(p, localVertexId);
    }
    // D1
    else if(triangleId < tsetshift_[3]) {
      triangleToPosition(triangleId, 3, p);
      vertexId = getTriangleVertexD1(p, localVertexId);
    }
    // D2
    else if(triangleId < tsetshift_[4]) {
      triangleToPosition(triangleId, 4, p);
      vertexId = getTriangleVertexD2(p, localVertexId);
    }
    // D3
    else if(triangleId < tsetshift_[5]) {
      triangleToPosition(triangleId, 5, p);
      vertexId = getTriangleVertexD3(p, localVertexId);
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    triangleToPosition2d(triangleId, p);
    const SimplexId id = triangleId % 2;

    SimplexId wrapXRight = 0;
    SimplexId wrapYBottom = 0;
    if(p[0] / 2 == nbvoxels_[Di_])
      wrapXRight = -wrap_[0];
    if(p[1] == nbvoxels_[Dj_])
      wrapYBottom = -wrap_[1];
    if(id == 0) {
      switch(localVertexId) {
        case 0:
          vertexId = p[0] / 2 + p[1] * vshift_[0];
          break;
        case 1:
          vertexId = p[0] / 2 + p[1] * vshift_[0] + 1 + wrapXRight;
          break;
        case 2:
          vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0] + wrapYBottom;
          break;
      }
    } else {
      switch(localVertexId) {
        case 0:
          vertexId = p[0] / 2 + p[1] * vshift_[0] + 1 + wrapXRight;
          break;
        case 1:
          vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0] + 1 + wrapXRight
                     + wrapYBottom;
          break;
        case 2:
          vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0] + wrapYBottom;
          break;
      }
    }
  }

  return 0;
}

int PeriodicImplicitTriangulation::getTriangleEdge(const SimplexId &triangleId,
                                                   const int &localEdgeId,
                                                   SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
  if(localEdgeId < 0 or localEdgeId >= 3)
    return -2;
#endif

  edgeId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];
    const SimplexId id = triangleId % 2;

    // F
    if(triangleId < tsetshift_[0]) {
      triangleToPosition(triangleId, 0, p);

      if(id)
        edgeId = getTriangleEdgeF_1(p, localEdgeId);
      else
        edgeId = getTriangleEdgeF_0(p, localEdgeId);
    }
    // H
    else if(triangleId < tsetshift_[1]) {
      triangleToPosition(triangleId, 1, p);

      if(id)
        edgeId = getTriangleEdgeH_1(p, localEdgeId);
      else
        edgeId = getTriangleEdgeH_0(p, localEdgeId);
    }
    // C
    else if(triangleId < tsetshift_[2]) {
      triangleToPosition(triangleId, 2, p);

      if(id)
        edgeId = getTriangleEdgeC_1(p, localEdgeId);
      else
        edgeId = getTriangleEdgeC_0(p, localEdgeId);
    }
    // D1
    else if(triangleId < tsetshift_[3]) {
      triangleToPosition(triangleId, 3, p);

      if(id)
        edgeId = getTriangleEdgeD1_1(p, localEdgeId);
      else
        edgeId = getTriangleEdgeD1_0(p, localEdgeId);
    }
    // D2
    else if(triangleId < tsetshift_[4]) {
      triangleToPosition(triangleId, 4, p);

      if(id)
        edgeId = getTriangleEdgeD2_1(p, localEdgeId);
      else
        edgeId = getTriangleEdgeD2_0(p, localEdgeId);
    }
    // D3
    else if(triangleId < tsetshift_[5]) {
      triangleToPosition(triangleId, 5, p);

      if(id)
        edgeId = getTriangleEdgeD3_1(p, localEdgeId);
      else
        edgeId = getTriangleEdgeD3_0(p, localEdgeId);
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    const SimplexId id = triangleId % 2;
    triangleToPosition2d(triangleId, p);

    SimplexId wrapXRight = 0;
    SimplexId wrapYBottom = 0;
    if(p[0] / 2 == nbvoxels_[Di_])
      wrapXRight = -wrap_[0];
    if(p[1] == nbvoxels_[Dj_])
      wrapYBottom = -wrap_[1];
    if(id == 0) {
      switch(localEdgeId) {
        case 0:
          edgeId = p[0] / 2 + p[1] * eshift_[0];
          break;
        case 1:
          edgeId = esetshift_[0] + p[0] / 2 + p[1] * eshift_[2];
          break;
        case 2:
          edgeId = esetshift_[1] + p[0] / 2 + p[1] * eshift_[4];
          break;
      }
    } else {
      switch(localEdgeId) {
        case 0:
          edgeId = p[0] / 2 + (p[1] + 1) * eshift_[0] + wrapYBottom;
          break;
        case 1:
          edgeId
            = esetshift_[0] + (p[0] + 1) / 2 + p[1] * eshift_[2] + wrapXRight;
          break;
        case 2:
          edgeId = esetshift_[1] + p[0] / 2 + p[1] * eshift_[4];
          break;
      }
    }
  }

  return 0;
}

int PeriodicImplicitTriangulation::getTriangleEdges(
  vector<vector<SimplexId>> &edges) const {
  edges.resize(triangleNumber_);
  for(SimplexId i = 0; i < triangleNumber_; ++i) {
    edges[i].resize(3);
    for(int j = 0; j < 3; ++j)
      getTriangleEdge(i, j, edges[i][j]);
  }
  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getTriangleEdges() {
  if(!triangleEdgeList_.size()) {
    Timer t;

    getTriangleEdges(triangleEdgeList_);

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Triangle edges ("
          << triangleNumber_ << " triangle(s), " << edgeNumber_
          << " edge(s)) computed in " << t.getElapsedTime() << " s. (" << 1
          << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &triangleEdgeList_;
}

const vector<vector<SimplexId>> *PeriodicImplicitTriangulation::getTriangles() {
  if(!triangleList_.size()) {
    Timer t;

    triangleList_.resize(triangleNumber_);
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      triangleList_[i].resize(3);
      for(int j = 0; j < 3; ++j)
        getTriangleVertex(i, j, triangleList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Triangle list ("
          << triangleNumber_ << " triangles) computed in " << t.getElapsedTime()
          << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &triangleList_;
}

int PeriodicImplicitTriangulation::getTriangleLink(const SimplexId &triangleId,
                                                   const int &localLinkId,
                                                   SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getTriangleLinkNumber(triangleId))
    return -1;
#endif

  linkId = -1;

  if(dimensionality_ == 3) {
    SimplexId p[3];

    // F
    if(triangleId < tsetshift_[0]) {
      triangleToPosition(triangleId, 0, p);
      linkId = getTriangleLinkF(p, localLinkId);
    }
    // H
    else if(triangleId < tsetshift_[1]) {
      triangleToPosition(triangleId, 1, p);
      linkId = getTriangleLinkH(p, localLinkId);
    }
    // C
    else if(triangleId < tsetshift_[2]) {
      triangleToPosition(triangleId, 2, p);
      linkId = getTriangleLinkC(p, localLinkId);
    }
    // D1
    else if(triangleId < tsetshift_[3]) {
      triangleToPosition(triangleId, 3, p);
      linkId = getTriangleLinkD1(p, localLinkId);
    }
    // D2
    else if(triangleId < tsetshift_[4]) {
      triangleToPosition(triangleId, 4, p);
      linkId = getTriangleLinkD2(p, localLinkId);
    }
    // D3
    else if(triangleId < tsetshift_[5]) {
      triangleToPosition(triangleId, 5, p);
      linkId = getTriangleLinkD3(p, localLinkId);
    }
  }

  return 0;
}

inline SimplexId PeriodicImplicitTriangulation::getTriangleLinkNumber(
  const SimplexId &triangleId) const {
  return getTriangleStarNumber(triangleId);
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getTriangleLinks() {
  if(!triangleLinkList_.size()) {
    Timer t;

    triangleLinkList_.resize(triangleNumber_);
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      triangleLinkList_[i].resize(getTriangleLinkNumber(i));
      for(SimplexId j = 0; j < (SimplexId)triangleLinkList_[i].size(); ++j)
        getTriangleLink(i, j, triangleLinkList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[TriangulationVTI] Triangle links built in " << t.getElapsedTime()
          << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }
  return &triangleLinkList_;
}

inline SimplexId PeriodicImplicitTriangulation::getTriangleStarNumber(
  const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    return 2;
  }
  return 0;
}

int PeriodicImplicitTriangulation::getTriangleStar(const SimplexId &triangleId,
                                                   const int &localStarId,
                                                   SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getTriangleStarNumber(triangleId))
    return -1;
#endif

  starId = -1;
  if(dimensionality_ == 3) {
    SimplexId p[3];

    // F
    if(triangleId < tsetshift_[0]) {
      triangleToPosition(triangleId, 0, p);
      starId = getTriangleStarF(p, localStarId);
    }
    // H
    else if(triangleId < tsetshift_[1]) {
      triangleToPosition(triangleId, 1, p);
      starId = getTriangleStarH(p, localStarId);
    }
    // C
    else if(triangleId < tsetshift_[2]) {
      triangleToPosition(triangleId, 2, p);
      starId = getTriangleStarC(p, localStarId);
    }
    // D1
    else if(triangleId < tsetshift_[3]) {
      triangleToPosition(triangleId, 3, p);
      starId = getTriangleStarD1(p, localStarId);
    }
    // D2
    else if(triangleId < tsetshift_[4]) {
      triangleToPosition(triangleId, 4, p);
      starId = getTriangleStarD2(p, localStarId);
    }
    // D3
    else if(triangleId < tsetshift_[5]) {
      triangleToPosition(triangleId, 5, p);
      starId = getTriangleStarD3(p, localStarId);
    }
  }
  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getTriangleStars() {
  if(!triangleStarList_.size()) {
    Timer t;

    triangleStarList_.resize(triangleNumber_);
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      triangleStarList_[i].resize(getTriangleStarNumber(i));
      for(SimplexId j = 0; j < (SimplexId)triangleStarList_[i].size(); ++j)
        getTriangleStar(i, j, triangleStarList_[i][j]);
    }

    {
      stringstream msg;
      msg << "[TriangulationVTI] Triangle stars built in " << t.getElapsedTime()
          << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }
  return &triangleStarList_;
}

inline SimplexId PeriodicImplicitTriangulation::getTriangleNeighborNumber(
  const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

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

  if(dimensionality_ == 2) {
    SimplexId p[2];
    triangleToPosition2d(triangleId, p);
    const SimplexId id = triangleId % 2;

    if(id) {
      if(p[0] / 2 == nbvoxels_[Di_] and p[1] == nbvoxels_[Dj_]) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId - 1;
            break;
          case 1:
            neighborId = triangleId + 1 - wrap_[0] * 2;
            break;
          case 2:
            neighborId = triangleId + tshift_[0] - 1 - wrap_[1] * 2;
            break;
        }
      } else if(p[0] / 2 == nbvoxels_[Di_]) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId - 1;
            break;
          case 1:
            neighborId = triangleId + 1 - wrap_[0] * 2;
            break;
          case 2:
            neighborId = triangleId + tshift_[0] - 1;
            break;
        }
      } else if(p[1] == nbvoxels_[Dj_]) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId - 1;
            break;
          case 1:
            neighborId = triangleId + 1;
            break;
          case 2:
            neighborId = triangleId + tshift_[0] - 1 - wrap_[1] * 2;
            break;
        }
      } else {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId - 1;
            break;
          case 1:
            neighborId = triangleId + 1;
            break;
          case 2:
            neighborId = triangleId + tshift_[0] - 1;
            break;
        }
      }
    } else {
      if(p[0] / 2 == 0 and p[1] == 0) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId + 1;
            break;
          case 1:
            neighborId = triangleId - 1 + wrap_[0] * 2;
            break;
          case 2:
            neighborId = triangleId - tshift_[0] + 1 + wrap_[1] * 2;
            break;
        }
      } else if(p[0] / 2 == 0) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId + 1;
            break;
          case 1:
            neighborId = triangleId - 1 + wrap_[0] * 2;
            break;
          case 2:
            neighborId = triangleId - tshift_[0] + 1;
            break;
        }
      } else if(p[1] == 0) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId + 1;
            break;
          case 1:
            neighborId = triangleId - 1;
            break;
          case 2:
            neighborId = triangleId - tshift_[0] + 1 + wrap_[1] * 2;
            break;
        }
      } else {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId + 1;
            break;
          case 1:
            neighborId = triangleId - 1;
            break;
          case 2:
            neighborId = triangleId - tshift_[0] + 1;
            break;
        }
      }
    }
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
    SimplexId p[3];
    tetrahedronToPosition(tetId, p);
    const SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        vertexId = getTetrahedronVertexABCG(p, localVertexId);
        break;
      case 1:
        vertexId = getTetrahedronVertexBCDG(p, localVertexId);
        break;
      case 2:
        vertexId = getTetrahedronVertexABEG(p, localVertexId);
        break;
      case 3:
        vertexId = getTetrahedronVertexBEFG(p, localVertexId);
        break;
      case 4:
        vertexId = getTetrahedronVertexBFGH(p, localVertexId);
        break;
      case 5:
        vertexId = getTetrahedronVertexBDGH(p, localVertexId);
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
    SimplexId p[3];
    tetrahedronToPosition(tetId, p);
    const SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        edgeId = getTetrahedronEdgeABCG(p, localEdgeId);
        break;
      case 1:
        edgeId = getTetrahedronEdgeBCDG(p, localEdgeId);
        break;
      case 2:
        edgeId = getTetrahedronEdgeABEG(p, localEdgeId);
        break;
      case 3:
        edgeId = getTetrahedronEdgeBEFG(p, localEdgeId);
        break;
      case 4:
        edgeId = getTetrahedronEdgeBFGH(p, localEdgeId);
        break;
      case 5:
        edgeId = getTetrahedronEdgeBDGH(p, localEdgeId);
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
    SimplexId p[3];
    tetrahedronToPosition(tetId, p);
    const SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        triangleId = getTetrahedronTriangleABCG(p, localTriangleId);
        break;
      case 1:
        triangleId = getTetrahedronTriangleBCDG(p, localTriangleId);
        break;
      case 2:
        triangleId = getTetrahedronTriangleABEG(p, localTriangleId);
        break;
      case 3:
        triangleId = getTetrahedronTriangleBEFG(p, localTriangleId);
        break;
      case 4:
        triangleId = getTetrahedronTriangleBFGH(p, localTriangleId);
        break;
      case 5:
        triangleId = getTetrahedronTriangleBDGH(p, localTriangleId);
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
    const SimplexId id = tetId % 6;
    SimplexId p[3];
    tetrahedronToPosition(tetId, p);

    switch(id) {
      case 0:
        neighborId = getTetrahedronNeighborABCG(tetId, p, localNeighborId);
        break;
      case 1:
        neighborId = getTetrahedronNeighborBCDG(tetId, p, localNeighborId);
        break;
      case 2:
        neighborId = getTetrahedronNeighborABEG(tetId, p, localNeighborId);
        break;
      case 3:
        neighborId = getTetrahedronNeighborBEFG(tetId, p, localNeighborId);
        break;
      case 4:
        neighborId = getTetrahedronNeighborBFGH(tetId, p, localNeighborId);
        break;
      case 5:
        neighborId = getTetrahedronNeighborBDGH(tetId, p, localNeighborId);
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

SimplexId PeriodicImplicitTriangulation::getCellVertexNumber(
  const SimplexId &cellId) const {
  return dimensionality_ + 1;
}

int PeriodicImplicitTriangulation::getCellVertex(const SimplexId &cellId,
                                                 const int &localVertexId,
                                                 SimplexId &vertexId) const {
  if(dimensionality_ == 3)
    getTetrahedronVertex(cellId, localVertexId, vertexId);
  else if(dimensionality_ == 2)
    getTriangleVertex(cellId, localVertexId, vertexId);
  else if(dimensionality_ == 1)
    getEdgeVertex(cellId, localVertexId, vertexId);

  return 0;
}

SimplexId PeriodicImplicitTriangulation::getCellEdgeNumber(
  const SimplexId &cellId) const {
  if(dimensionality_ == 3)
    return 6;
  else if(dimensionality_ == 2)
    return 3;

  return 0;
}

int PeriodicImplicitTriangulation::getCellEdge(const SimplexId &cellId,
                                               const int &localEdgeId,
                                               SimplexId &edgeId) const {
  if(dimensionality_ == 3)
    getTetrahedronEdge(cellId, localEdgeId, edgeId);
  else if(dimensionality_ == 2)
    getTriangleEdge(cellId, localEdgeId, edgeId);

  return 0;
}

const vector<vector<SimplexId>> *PeriodicImplicitTriangulation::getCellEdges() {
  if(!cellEdgeList_.size()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronEdges(cellEdgeList_);
    else if(dimensionality_ == 2)
      getTriangleEdges(cellEdgeList_);

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Cell edges ("
          << getNumberOfCells() << " cell(s), " << edgeNumber_
          << "edge(s)) computed in " << t.getElapsedTime() << " s. (" << 1
          << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &cellEdgeList_;
}

int PeriodicImplicitTriangulation::getCellTriangle(
  const SimplexId &cellId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
  if(dimensionality_ == 3)
    getTetrahedronTriangle(cellId, localTriangleId, triangleId);

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getCellTriangles() {
  if(!cellTriangleList_.size()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronTriangles(cellTriangleList_);

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Cell triangles (" << cellNumber_
          << " cell(s), " << triangleNumber_ << "edge(s)) computed in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &cellTriangleList_;
}

SimplexId PeriodicImplicitTriangulation::getCellNeighborNumber(
  const SimplexId &cellId) const {
  if(dimensionality_ == 3)
    return getTetrahedronNeighborNumber(cellId);
  else if(dimensionality_ == 2)
    return getTriangleNeighborNumber(cellId);
  else if(dimensionality_ == 1) {
    stringstream msg;
    msg << "[PeriodicImplicitTriangulation] getCellNeighborNumber() in 1D:"
        << endl;
    msg << "[PeriodicImplicitTriangulation] Not implemented! TODO!" << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
    return -1;
  }

  return 0;
}

int PeriodicImplicitTriangulation::getCellNeighbor(
  const SimplexId &cellId,
  const int &localNeighborId,
  SimplexId &neighborId) const {
  if(dimensionality_ == 3)
    getTetrahedronNeighbor(cellId, localNeighborId, neighborId);
  else if(dimensionality_ == 2)
    getTriangleNeighbor(cellId, localNeighborId, neighborId);
  else if(dimensionality_ == 1) {
    stringstream msg;
    msg << "[PeriodicImplicitTriangulation] getCellNeighbor() in 1D:" << endl;
    msg << "[PeriodicImplicitTriangulation] Not implemented! TODO!" << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
    return -1;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  PeriodicImplicitTriangulation::getCellNeighbors() {
  if(!cellNeighborList_.size()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronNeighbors(cellNeighborList_);
    else if(dimensionality_ == 2)
      getTriangleNeighbors(cellNeighborList_);
    else if(dimensionality_ == 1) {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] getCellNeighbors() in 1D:"
          << endl;
      msg << "[PeriodicImplicitTriangulation] Not implemented! TODO!" << endl;
      dMsg(cerr, msg.str(), Debug::fatalMsg);
      return nullptr;
    }

    {
      stringstream msg;
      msg << "[PeriodicImplicitTriangulation] Cell neighbors ("
          << getNumberOfCells() << " cells) computed in " << t.getElapsedTime()
          << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &cellNeighborList_;
}
