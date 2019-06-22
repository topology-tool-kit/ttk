#include <ImplicitTriangulation.h>

using namespace std;
using namespace ttk;

ImplicitTriangulation::ImplicitTriangulation()
  : dimensionality_{-1}, cellNumber_{}, vertexNumber_{}, edgeNumber_{},
    triangleNumber_{}, tetrahedronNumber_{}, isAccelerated_{} {
}

ImplicitTriangulation::~ImplicitTriangulation() {
}

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
    msg << "[ImplicitTriangulation] The getVertex*() requests are accelerated."
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}

bool ImplicitTriangulation::isPowerOfTwo(unsigned long long int v,
                                         unsigned long long int &r) {
  if(v && !(v & (v - 1))) {
    r = 0;
    while(v >>= 1)
      r++;
    return true;
  }
  return false;
}

bool ImplicitTriangulation::isVertexOnBoundary(
  const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    SimplexId p[3];
    vertexToPosition(vertexId, p);

    return (p[0] == 0 or p[1] == 0 or p[2] == 0 or p[0] == nbvoxels_[0]
            or p[1] == nbvoxels_[1] or p[2] == nbvoxels_[2]);
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);

    return (p[0] == 0 or p[1] == 0 or p[0] == nbvoxels_[Di_]
            or p[1] == nbvoxels_[Dj_]);
  }

  return false;
}

bool ImplicitTriangulation::isEdgeOnBoundary(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    SimplexId p[3];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition(edgeId, 0, p);
      return (p[1] == 0 or p[1] == nbvoxels_[1] or p[2] == 0
              or p[2] == nbvoxels_[2]);
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition(edgeId, 1, p);
      return (p[0] == 0 or p[0] == nbvoxels_[0] or p[2] == 0
              or p[2] == nbvoxels_[2]);
    }
    // P
    else if(edgeId < esetshift_[2]) {
      edgeToPosition(edgeId, 2, p);
      return (p[0] == 0 or p[1] == 0 or p[0] == nbvoxels_[0]
              or p[1] == nbvoxels_[1]);
    }
    // D1
    else if(edgeId < esetshift_[3]) {
      edgeToPosition(edgeId, 3, p);
      return (p[2] == 0 or p[2] == nbvoxels_[2]);
    }
    // D2
    else if(edgeId < esetshift_[4]) {
      edgeToPosition(edgeId, 4, p);
      return (p[0] == 0 or p[0] == nbvoxels_[0]);
    }
    // D3
    else if(edgeId < esetshift_[5]) {
      edgeToPosition(edgeId, 5, p);
      return (p[1] == 0 or p[1] == nbvoxels_[1]);
    } else
      return false;
  } else if(dimensionality_ == 2)
    return (getEdgeStarNumber(edgeId) == 1);

  return false;
}

bool ImplicitTriangulation::isTriangleOnBoundary(
  const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  if(dimensionality_ == 3)
    return (getTriangleStarNumber(triangleId) == 1);

  return false;
}

inline SimplexId ImplicitTriangulation::getVertexNeighborNumber(
  const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    SimplexId p[3];
    vertexToPosition(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[0]) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 14; // abcdefgh
        else
          return 10; // abdc ou efhg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 10; // aefb
        else if(p[2] == 0)
          return 8; // ab
        else
          return 6; // ef
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 10; // ghdc
        else if(p[2] == 0)
          return 6; // cd
        else
          return 8; // gh
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 10; // aegc
        else if(p[2] == 0)
          return 6; // ac
        else
          return 8; // eg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 6; // ae
        else
          return 4; // a ou e
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 8; // cg
        else if(p[2] == 0)
          return 4; // c
        else
          return 7; // g
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 10; // bfhd
        else if(p[2] == 0)
          return 8; // bd
        else
          return 6; // fh
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 8; // bf
        else if(p[2] == 0)
          return 7; // b
        else
          return 4; // f
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 6; // dh
        else
          return 4; // d ou h
      }
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[Di_]) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        return 6; // abcd
      else if(p[1] == 0)
        return 4; // ab
      else
        return 4; // cd
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        return 4; // ac
      else if(p[1] == 0)
        return 2; // a
      else
        return 3; // c
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        return 4; // bd
      else if(p[1] == 0)
        return 3; // b
      else
        return 2; // d
    }
  } else if(dimensionality_ == 1) {
    if(vertexId > 0 and vertexId < nbvoxels_[Di_])
      return 2; // ab
    else
      return 1; // a ou b
  }

  return -1;
}

int ImplicitTriangulation::getVertexNeighbor(const SimplexId &vertexId,
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

    if(0 < p[0] and p[0] < nbvoxels_[0]) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          neighborId
            = getVertexNeighborABCDEFGH(vertexId, localNeighborId); // abcdefgh
        else if(p[2] == 0)
          neighborId = getVertexNeighborABDC(vertexId, localNeighborId); // abdc
        else
          neighborId = getVertexNeighborEFHG(vertexId, localNeighborId); // efhg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          neighborId = getVertexNeighborAEFB(vertexId, localNeighborId); // aefb
        else if(p[2] == 0)
          neighborId = getVertexNeighborAB(vertexId, localNeighborId); // ab
        else
          neighborId = getVertexNeighborEF(vertexId, localNeighborId); // ef
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          neighborId = getVertexNeighborGHDC(vertexId, localNeighborId); // ghdc
        else if(p[2] == 0)
          neighborId = getVertexNeighborCD(vertexId, localNeighborId); // cd
        else
          neighborId = getVertexNeighborGH(vertexId, localNeighborId); // gh
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          neighborId = getVertexNeighborAEGC(vertexId, localNeighborId); // aegc
        else if(p[2] == 0)
          neighborId = getVertexNeighborAC(vertexId, localNeighborId); // ac
        else
          neighborId = getVertexNeighborEG(vertexId, localNeighborId); // eg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          neighborId = getVertexNeighborAE(vertexId, localNeighborId); // ae
        else if(p[2] == 0)
          neighborId = getVertexNeighborA(vertexId, localNeighborId); // a
        else
          neighborId = getVertexNeighborE(vertexId, localNeighborId); // e
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          neighborId = getVertexNeighborCG(vertexId, localNeighborId); // cg
        else if(p[2] == 0)
          neighborId = getVertexNeighborC(vertexId, localNeighborId); // c
        else
          neighborId = getVertexNeighborG(vertexId, localNeighborId); // g
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          neighborId = getVertexNeighborBFHD(vertexId, localNeighborId); // bfhd
        else if(p[2] == 0)
          neighborId = getVertexNeighborBD(vertexId, localNeighborId); // bd
        else
          neighborId = getVertexNeighborFH(vertexId, localNeighborId); // fh
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          neighborId = getVertexNeighborBF(vertexId, localNeighborId); // bf
        else if(p[2] == 0)
          neighborId = getVertexNeighborB(vertexId, localNeighborId); // b
        else
          neighborId = getVertexNeighborF(vertexId, localNeighborId); // f
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          neighborId = getVertexNeighborDH(vertexId, localNeighborId); // dh
        else if(p[2] == 0)
          neighborId = getVertexNeighborD(vertexId, localNeighborId); // d
        else
          neighborId = getVertexNeighborH(vertexId, localNeighborId); // h
      }
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[Di_]) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        neighborId = getVertexNeighbor2dABCD(vertexId, localNeighborId); // abcd
      else if(p[1] == 0)
        neighborId = getVertexNeighbor2dAB(vertexId, localNeighborId); // ab
      else
        neighborId = getVertexNeighbor2dCD(vertexId, localNeighborId); // cd
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        neighborId = getVertexNeighbor2dAC(vertexId, localNeighborId); // ac
      else if(p[1] == 0)
        neighborId = getVertexNeighbor2dA(vertexId, localNeighborId); // a
      else
        neighborId = getVertexNeighbor2dC(vertexId, localNeighborId); // c
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        neighborId = getVertexNeighbor2dBD(vertexId, localNeighborId); // bd
      else if(p[1] == 0)
        neighborId = getVertexNeighbor2dB(vertexId, localNeighborId); // b
      else
        neighborId = getVertexNeighbor2dD(vertexId, localNeighborId); // d
    }
  } else if(dimensionality_ == 1) {
    // ab
    if(vertexId > 0 and vertexId < nbvoxels_[Di_]) {
      if(localNeighborId == 0)
        neighborId = vertexId + 1;
      else
        neighborId = vertexId - 1;
    } else if(vertexId == 0)
      neighborId = vertexId + 1; // a
    else
      neighborId = vertexId - 1; // b
  }

  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getVertexNeighbors() {
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
      msg << "[ImplicitTriangulation] Vertex neighbors built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &vertexNeighborList_;
}

SimplexId
  ImplicitTriangulation::getVertexEdgeNumber(const SimplexId &vertexId) const {
  return getVertexNeighborNumber(vertexId);
}

int ImplicitTriangulation::getVertexEdge(const SimplexId &vertexId,
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

    if(0 < p[0] and p[0] < nbvoxels_[0]) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          edgeId = getVertexEdgeABCDEFGH(p, localEdgeId); // abcdefgh
        else if(p[2] == 0)
          edgeId = getVertexEdgeABDC(p, localEdgeId); // abdc
        else
          edgeId = getVertexEdgeEFHG(p, localEdgeId); // efhg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          edgeId = getVertexEdgeAEFB(p, localEdgeId); // aefb
        else if(p[2] == 0)
          edgeId = getVertexEdgeAB(p, localEdgeId); // ab
        else
          edgeId = getVertexEdgeEF(p, localEdgeId); // ef
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          edgeId = getVertexEdgeGHDC(p, localEdgeId); // ghdc
        else if(p[2] == 0)
          edgeId = getVertexEdgeCD(p, localEdgeId); // cd
        else
          edgeId = getVertexEdgeGH(p, localEdgeId); // gh
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          edgeId = getVertexEdgeAEGC(p, localEdgeId); // aegc
        else if(p[2] == 0)
          edgeId = getVertexEdgeAC(p, localEdgeId); // ac
        else
          edgeId = getVertexEdgeEG(p, localEdgeId); // eg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          edgeId = getVertexEdgeAE(p, localEdgeId); // ae
        else if(p[2] == 0)
          edgeId = getVertexEdgeA(p, localEdgeId); // a
        else
          edgeId = getVertexEdgeE(p, localEdgeId); // e
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          edgeId = getVertexEdgeCG(p, localEdgeId); // cg
        else if(p[2] == 0)
          edgeId = getVertexEdgeC(p, localEdgeId); // c
        else
          edgeId = getVertexEdgeG(p, localEdgeId); // g
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          edgeId = getVertexEdgeBFHD(p, localEdgeId); // bfhd
        else if(p[2] == 0)
          edgeId = getVertexEdgeBD(p, localEdgeId); // bd
        else
          edgeId = getVertexEdgeFH(p, localEdgeId); // fh
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          edgeId = getVertexEdgeBF(p, localEdgeId); // bf
        else if(p[2] == 0)
          edgeId = getVertexEdgeB(p, localEdgeId); // b
        else
          edgeId = getVertexEdgeF(p, localEdgeId); // f
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          edgeId = getVertexEdgeDH(p, localEdgeId); // dh
        else if(p[2] == 0)
          edgeId = getVertexEdgeD(p, localEdgeId); // d
        else
          edgeId = getVertexEdgeH(p, localEdgeId); // h
      }
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[Di_]) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        edgeId = getVertexEdge2dABCD(p, localEdgeId); // abcd
      else if(p[1] == 0)
        edgeId = getVertexEdge2dAB(p, localEdgeId); // ab
      else
        edgeId = getVertexEdge2dCD(p, localEdgeId); // cd
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        edgeId = getVertexEdge2dAC(p, localEdgeId); // ac
      else if(p[1] == 0)
        edgeId = getVertexEdge2dA(p, localEdgeId); // a
      else
        edgeId = getVertexEdge2dC(p, localEdgeId); // c
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        edgeId = getVertexEdge2dBD(p, localEdgeId); // bd
      else if(p[1] == 0)
        edgeId = getVertexEdge2dB(p, localEdgeId); // b
      else
        edgeId = getVertexEdge2dD(p, localEdgeId); // d
    }
  } else if(dimensionality_ == 1) {
    // ab
    if(vertexId > 0 and vertexId < nbvoxels_[Di_]) {
      if(localEdgeId == 0)
        edgeId = vertexId;
      else
        edgeId = vertexId - 1;
    } else if(vertexId == 0)
      edgeId = vertexId; // a
    else
      edgeId = vertexId - 1; // b
  }

  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getVertexEdges() {
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
      msg << "[ImplicitTriangulation] Vertex edges built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &vertexEdgeList_;
}

inline SimplexId ImplicitTriangulation::getVertexTriangleNumber(
  const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    SimplexId p[3];
    vertexToPosition(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[0]) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 36; // abcdefgh
        else
          return 21; // abdc ou efhg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 21; // aefb
        else if(p[2] == 0)
          return 15; // ab
        else
          return 9; // ef
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 21; // ghdc
        else if(p[2] == 0)
          return 9; // cd
        else
          return 15; // gh
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 21; // aegc
        else if(p[2] == 0)
          return 9; // ac
        else
          return 15; // eg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 9; // ae
        else
          return 5; // a ou e
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 15; // cg
        else if(p[2] == 0)
          return 5; // c
        else
          return 12; // g
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 21; // bfhd
        else if(p[2] == 0)
          return 15; // bd
        else
          return 9; // fh
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 15; // bf
        else if(p[2] == 0)
          return 12; // b
        else
          return 5; // f
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 9; // dh
        else
          return 5; // d ou h
      }
    }
  }

  return 0;
}

int ImplicitTriangulation::getVertexTriangle(const SimplexId &vertexId,
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

    if(0 < p[0] and p[0] < nbvoxels_[0]) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          triangleId
            = getVertexTriangleABCDEFGH(p, localTriangleId); // abcdefgh
        else if(p[2] == 0)
          triangleId = getVertexTriangleABDC(p, localTriangleId); // abdc
        else
          triangleId = getVertexTriangleEFHG(p, localTriangleId); // efhg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          triangleId = getVertexTriangleAEFB(p, localTriangleId); // aefb
        else if(p[2] == 0)
          triangleId = getVertexTriangleAB(p, localTriangleId); // ab
        else
          triangleId = getVertexTriangleEF(p, localTriangleId); // ef
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          triangleId = getVertexTriangleGHDC(p, localTriangleId); // ghdc
        else if(p[2] == 0)
          triangleId = getVertexTriangleCD(p, localTriangleId); // cd
        else
          triangleId = getVertexTriangleGH(p, localTriangleId); // gh
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          triangleId = getVertexTriangleAEGC(p, localTriangleId); // aegc
        else if(p[2] == 0)
          triangleId = getVertexTriangleAC(p, localTriangleId); // ac
        else
          triangleId = getVertexTriangleEG(p, localTriangleId); // eg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          triangleId = getVertexTriangleAE(p, localTriangleId); // ae
        else if(p[2] == 0)
          triangleId = getVertexTriangleA(p, localTriangleId); // a
        else
          triangleId = getVertexTriangleE(p, localTriangleId); // e
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          triangleId = getVertexTriangleCG(p, localTriangleId); // cg
        else if(p[2] == 0)
          triangleId = getVertexTriangleC(p, localTriangleId); // c
        else
          triangleId = getVertexTriangleG(p, localTriangleId); // g
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          triangleId = getVertexTriangleBFHD(p, localTriangleId); // bfhd
        else if(p[2] == 0)
          triangleId = getVertexTriangleBD(p, localTriangleId); // bd
        else
          triangleId = getVertexTriangleFH(p, localTriangleId); // fh
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          triangleId = getVertexTriangleBF(p, localTriangleId); // bf
        else if(p[2] == 0)
          triangleId = getVertexTriangleB(p, localTriangleId); // b
        else
          triangleId = getVertexTriangleF(p, localTriangleId); // f
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          triangleId = getVertexTriangleDH(p, localTriangleId); // dh
        else if(p[2] == 0)
          triangleId = getVertexTriangleD(p, localTriangleId); // d
        else
          triangleId = getVertexTriangleH(p, localTriangleId); // h
      }
    }
  }

  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getVertexTriangles() {
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
      msg << "[ImplicitTriangulation] Vertex triangles built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &vertexTriangleList_;
}

SimplexId
  ImplicitTriangulation::getVertexLinkNumber(const SimplexId &vertexId) const {
  return getVertexStarNumber(vertexId);
}

int ImplicitTriangulation::getVertexLink(const SimplexId &vertexId,
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

    if(0 < p[0] and p[0] < nbvoxels_[0]) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          linkId = getVertexLinkABCDEFGH(p, localLinkId); // abcdefgh
        else if(p[2] == 0)
          linkId = getVertexLinkABDC(p, localLinkId); // abdc
        else
          linkId = getVertexLinkEFHG(p, localLinkId); // efhg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          linkId = getVertexLinkAEFB(p, localLinkId); // aefb
        else if(p[2] == 0)
          linkId = getVertexLinkAB(p, localLinkId); // ab
        else
          linkId = getVertexLinkEF(p, localLinkId); // ef
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          linkId = getVertexLinkGHDC(p, localLinkId); // ghdc
        else if(p[2] == 0)
          linkId = getVertexLinkCD(p, localLinkId); // cd
        else
          linkId = getVertexLinkGH(p, localLinkId); // gh
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          linkId = getVertexLinkAEGC(p, localLinkId); // aegc
        else if(p[2] == 0)
          linkId = getVertexLinkAC(p, localLinkId); // ac
        else
          linkId = getVertexLinkEG(p, localLinkId); // eg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          linkId = getVertexLinkAE(p, localLinkId); // ae
        else if(p[2] == 0)
          linkId = getVertexLinkA(p, localLinkId); // a
        else
          linkId = getVertexLinkE(p, localLinkId); // e
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          linkId = getVertexLinkCG(p, localLinkId); // cg
        else if(p[2] == 0)
          linkId = getVertexLinkC(p, localLinkId); // c
        else
          linkId = getVertexLinkG(p, localLinkId); // g
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          linkId = getVertexLinkBFHD(p, localLinkId); // bfhd
        else if(p[2] == 0)
          linkId = getVertexLinkBD(p, localLinkId); // bd
        else
          linkId = getVertexLinkFH(p, localLinkId); // fh
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          linkId = getVertexLinkBF(p, localLinkId); // bf
        else if(p[2] == 0)
          linkId = getVertexLinkB(p, localLinkId); // b
        else
          linkId = getVertexLinkF(p, localLinkId); // f
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          linkId = getVertexLinkDH(p, localLinkId); // dh
        else if(p[2] == 0)
          linkId = getVertexLinkD(p, localLinkId); // d
        else
          linkId = getVertexLinkH(p, localLinkId); // h
      }
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[Di_]) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        linkId = getVertexLink2dABCD(p, localLinkId); // abcd
      else if(p[1] == 0)
        linkId = getVertexLink2dAB(p, localLinkId); // ab
      else
        linkId = getVertexLink2dCD(p, localLinkId); // cd
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        linkId = getVertexLink2dAC(p, localLinkId); // ac
      else if(p[1] == 0)
        linkId = getVertexLink2dA(p, localLinkId); // a
      else
        linkId = getVertexLink2dC(p, localLinkId); // c
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        linkId = getVertexLink2dBD(p, localLinkId); // bd
      else if(p[1] == 0)
        linkId = getVertexLink2dB(p, localLinkId); // b
      else
        linkId = getVertexLink2dD(p, localLinkId); // d
    }
  }

  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getVertexLinks() {
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
      msg << "[ImplicitTriangulation] Vertex links built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &vertexLinkList_;
}

inline SimplexId
  ImplicitTriangulation::getVertexStarNumber(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    SimplexId p[3];
    vertexToPosition(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[0]) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 24; // abcdefgh
        else
          return 12; // abdc ou efhg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 12; // aefb
        else if(p[2] == 0)
          return 8; // ab
        else
          return 4; // ef
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 12; // ghdc
        else if(p[2] == 0)
          return 4; // cd
        else
          return 8; // gh
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 12; // aegc
        else if(p[2] == 0)
          return 4; // ac
        else
          return 8; // eg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 4; // ae
        else
          return 2; // a ou e
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 8; // cg
        else if(p[2] == 0)
          return 2; // c
        else
          return 6; // g
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 12; // bfhd
        else if(p[2] == 0)
          return 8; // bd
        else
          return 4; // fh
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 8; // bf
        else if(p[2] == 0)
          return 6; // b
        else
          return 2; // f
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 4; // dh
        else
          return 2; // d ou h
      }
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[Di_]) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        return 6; // abcd
      else if(p[1] == 0)
        return 3; // ab
      else
        return 3; // cd
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        return 3; // ac
      else if(p[1] == 0)
        return 1; // a
      else
        return 2; // c
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        return 3; // bd
      else if(p[1] == 0)
        return 2; // b
      else
        return 1; // d
    }
  }

  return 0;
}

int ImplicitTriangulation::getVertexStar(const SimplexId &vertexId,
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

    if(0 < p[0] and p[0] < nbvoxels_[0]) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          starId = getVertexStarABCDEFGH(p, localStarId); // abcdefgh
        else if(p[2] == 0)
          starId = getVertexStarABDC(p, localStarId); // abdc
        else
          starId = getVertexStarEFHG(p, localStarId); // efhg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          starId = getVertexStarAEFB(p, localStarId); // aefb
        else if(p[2] == 0)
          starId = getVertexStarAB(p, localStarId); // ab
        else
          starId = getVertexStarEF(p, localStarId); // ef
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          starId = getVertexStarGHDC(p, localStarId); // ghdc
        else if(p[2] == 0)
          starId = getVertexStarCD(p, localStarId); // cd
        else
          starId = getVertexStarGH(p, localStarId); // gh
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          starId = getVertexStarAEGC(p, localStarId); // aegc
        else if(p[2] == 0)
          starId = getVertexStarAC(p, localStarId); // ac
        else
          starId = getVertexStarEG(p, localStarId); // eg
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          starId = getVertexStarAE(p, localStarId); // ae
        else if(p[2] == 0)
          starId = getVertexStarA(p, localStarId); // a
        else
          starId = getVertexStarE(p, localStarId); // e
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          starId = getVertexStarCG(p, localStarId); // cg
        else if(p[2] == 0)
          starId = getVertexStarC(p, localStarId); // c
        else
          starId = getVertexStarG(p, localStarId); // g
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          starId = getVertexStarBFHD(p, localStarId); // bfhd
        else if(p[2] == 0)
          starId = getVertexStarBD(p, localStarId); // bd
        else
          starId = getVertexStarFH(p, localStarId); // fh
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          starId = getVertexStarBF(p, localStarId); // bf
        else if(p[2] == 0)
          starId = getVertexStarB(p, localStarId); // b
        else
          starId = getVertexStarF(p, localStarId); // f
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          starId = getVertexStarDH(p, localStarId); // dh
        else if(p[2] == 0)
          starId = getVertexStarD(p, localStarId); // d
        else
          starId = getVertexStarH(p, localStarId); // h
      }
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];
    vertexToPosition2d(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[Di_]) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        starId = getVertexStar2dABCD(p, localStarId); // abcd
      else if(p[1] == 0)
        starId = getVertexStar2dAB(p, localStarId); // ab
      else
        starId = getVertexStar2dCD(p, localStarId); // cd
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        starId = getVertexStar2dAC(p, localStarId); // ac
      else if(p[1] == 0)
        starId = getVertexStar2dA(p, localStarId); // a
      else
        starId = getVertexStar2dC(p, localStarId); // c
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_])
        starId = getVertexStar2dBD(p, localStarId); // bd
      else if(p[1] == 0)
        starId = getVertexStar2dB(p, localStarId); // b
      else
        starId = getVertexStar2dD(p, localStarId); // d
    }
  }

  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getVertexStars() {
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
      msg << "[ImplicitTriangulation] Vertex stars built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &vertexStarList_;
}

int ImplicitTriangulation::getVertexPoint(const SimplexId &vertexId,
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

int ImplicitTriangulation::getEdgeVertex(const SimplexId &edgeId,
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

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition(edgeId, 0, p);
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + 1;
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
      }
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition(edgeId, 1, p);
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[0];
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0];
      }
    }
    // P
    else if(edgeId < esetshift_[2]) {
      edgeToPosition(edgeId, 2, p);
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[1];
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1];
      }
    }
    // D1
    else if(edgeId < esetshift_[3]) {
      edgeToPosition(edgeId, 3, p);
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + 1;
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[0];
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0];
      }
    }
    // D2
    else if(edgeId < esetshift_[4]) {
      edgeToPosition(edgeId, 4, p);
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[0]
                     + vshift_[1];
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
                     + vshift_[1];
      }
    }
    // D3
    else if(edgeId < esetshift_[5]) {
      edgeToPosition(edgeId, 5, p);
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + 1;
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[1];

      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1];
      }
    }
    // D4
    else if(edgeId < esetshift_[6]) {
      edgeToPosition(edgeId, 6, p);
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + 1;
        else
          vertexId = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]) + vshift_[0]
                     + vshift_[1];

      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
        else
          vertexId = p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
                     + vshift_[1];
      }
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition2d(edgeId, 0, p);
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + 1;
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0];
        else
          vertexId = p[0] + p[1] * vshift_[0] + 1;
      }
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition2d(edgeId, 1, p);
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]);
        else
          vertexId = p[0] + (p[1] << div_[0]) + vshift_[0];
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0];
        else
          vertexId = p[0] + p[1] * vshift_[0] + vshift_[0];
      }
    }
    // D1
    else if(edgeId < esetshift_[2]) {
      edgeToPosition2d(edgeId, 2, p);
      if(isAccelerated_) {
        if(localVertexId == 0)
          vertexId = p[0] + (p[1] << div_[0]) + 1;
        else
          vertexId = p[0] + (p[1] << div_[0]) + vshift_[0];
      } else {
        if(localVertexId == 0)
          vertexId = p[0] + p[1] * vshift_[0] + 1;
        else
          vertexId = p[0] + p[1] * vshift_[0] + vshift_[0];
      }
    }
  } else if(dimensionality_ == 1) {
    if(edgeId > 0 and edgeId < (edgeNumber_ - 1)) {
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
        vertexId = edgeNumber_ - 1;
      else
        vertexId = edgeNumber_;
    }
  }

  return 0;
}

const vector<pair<SimplexId, SimplexId>> *ImplicitTriangulation::getEdges() {
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
      msg << "[ImplicitTriangulation] Edge-list built in " << t.getElapsedTime()
          << " s. (" << edgeList_.size() << " edges, (" << 1 << " thread(s))"
          << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &edgeList_;
}

inline SimplexId
  ImplicitTriangulation::getEdgeTriangleNumber(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    SimplexId p[3];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition(edgeId, 0, p);

      if(p[2] > 0 and p[2] < nbvoxels_[2]) {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 6;
        else if(p[1] == 0)
          return 4;
        else
          return 4;
      } else if(p[2] == 0) {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 4;
        else if(p[1] == 0)
          return 3;
        else
          return 2;
      } else {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 4;
        else if(p[1] == 0)
          return 2;
        else
          return 3;
      }
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition(edgeId, 1, p);

      if(p[0] > 0 and p[0] < nbvoxels_[0]) {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          return 6;
        else if(p[2] == 0)
          return 4;
        else
          return 4;
      } else if(p[0] == 0) {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          return 4;
        else if(p[2] == 0)
          return 2;
        else
          return 3;
      } else {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          return 4;
        else if(p[2] == 0)
          return 3;
        else
          return 2;
      }
    }
    // P
    else if(edgeId < esetshift_[2]) {
      edgeToPosition(edgeId, 2, p);

      if(p[0] > 0 and p[0] < nbvoxels_[0]) {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 6;
        else if(p[1] == 0)
          return 4;
        else
          return 4;
      } else if(p[0] == 0) {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 4;
        else if(p[1] == 0)
          return 2;
        else
          return 3;
      } else {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 4;
        else if(p[1] == 0)
          return 3;
        else
          return 2;
      }
    }
    // D1
    else if(edgeId < esetshift_[3]) {
      edgeToPosition(edgeId, 3, p);

      if(p[2] > 0 and p[2] < nbvoxels_[2])
        return 4;
      else
        return 3;
    }
    // D2
    else if(edgeId < esetshift_[4]) {
      edgeToPosition(edgeId, 4, p);

      if(p[0] > 0 and p[0] < nbvoxels_[0])
        return 4;
      else
        return 3;
    }
    // D3
    else if(edgeId < esetshift_[5]) {
      edgeToPosition(edgeId, 5, p);

      if(p[1] > 0 and p[1] < nbvoxels_[1])
        return 4;
      else
        return 3;
    }
    // D4
    else if(edgeId < esetshift_[6])
      return 6;
  } else if(dimensionality_ == 2) {
    SimplexId p[2];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition2d(edgeId, 0, p);

      if(p[1] > 0 and p[1] < nbvoxels_[Dj_])
        return 2;
      else
        return 1;
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition2d(edgeId, 1, p);

      if(p[0] > 0 and p[0] < nbvoxels_[Di_])
        return 2;
      else
        return 1;
    }
    // D1
    else if(edgeId < esetshift_[2])
      return 2;
  }

  return 0;
}

int ImplicitTriangulation::getEdgeTriangle(const SimplexId &edgeId,
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

      if(p[1] > 0 and p[1] < nbvoxels_[1]) {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          triangleId = getEdgeTriangleL_xnn(p, localTriangleId);
        else if(p[2] == 0)
          triangleId = getEdgeTriangleL_xn0(p, localTriangleId);
        else
          triangleId = getEdgeTriangleL_xnN(p, localTriangleId);
      } else if(p[1] == 0) {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          triangleId = getEdgeTriangleL_x0n(p, localTriangleId);
        else if(p[2] == 0)
          triangleId = getEdgeTriangleL_x00(p, localTriangleId);
        else
          triangleId = getEdgeTriangleL_x0N(p, localTriangleId);
      } else {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          triangleId = getEdgeTriangleL_xNn(p, localTriangleId);
        else if(p[2] == 0)
          triangleId = getEdgeTriangleL_xN0(p, localTriangleId);
        else
          triangleId = getEdgeTriangleL_xNN(p, localTriangleId);
      }
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition(edgeId, 1, p);

      if(p[0] > 0 and p[0] < nbvoxels_[0]) {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          triangleId = getEdgeTriangleH_nyn(p, localTriangleId);
        else if(p[2] == 0)
          triangleId = getEdgeTriangleH_ny0(p, localTriangleId);
        else
          triangleId = getEdgeTriangleH_nyN(p, localTriangleId);
      } else if(p[0] == 0) {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          triangleId = getEdgeTriangleH_0yn(p, localTriangleId);
        else if(p[2] == 0)
          triangleId = getEdgeTriangleH_0y0(p, localTriangleId);
        else
          triangleId = getEdgeTriangleH_0yN(p, localTriangleId);
      } else {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          triangleId = getEdgeTriangleH_Nyn(p, localTriangleId);
        else if(p[2] == 0)
          triangleId = getEdgeTriangleH_Ny0(p, localTriangleId);
        else
          triangleId = getEdgeTriangleH_NyN(p, localTriangleId);
      }
    }
    // P
    else if(edgeId < esetshift_[2]) {
      edgeToPosition(edgeId, 2, p);

      if(p[0] > 0 and p[0] < nbvoxels_[0]) {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          triangleId = getEdgeTriangleP_nnz(p, localTriangleId);
        else if(p[1] == 0)
          triangleId = getEdgeTriangleP_n0z(p, localTriangleId);
        else
          triangleId = getEdgeTriangleP_nNz(p, localTriangleId);
      } else if(p[0] == 0) {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          triangleId = getEdgeTriangleP_0nz(p, localTriangleId);
        else if(p[1] == 0)
          triangleId = getEdgeTriangleP_00z(p, localTriangleId);
        else
          triangleId = getEdgeTriangleP_0Nz(p, localTriangleId);
      } else {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          triangleId = getEdgeTriangleP_Nnz(p, localTriangleId);
        else if(p[1] == 0)
          triangleId = getEdgeTriangleP_N0z(p, localTriangleId);
        else
          triangleId = getEdgeTriangleP_NNz(p, localTriangleId);
      }
    }
    // D1
    else if(edgeId < esetshift_[3]) {
      edgeToPosition(edgeId, 3, p);

      if(p[2] > 0 and p[2] < nbvoxels_[2])
        triangleId = getEdgeTriangleD1_xyn(p, localTriangleId);
      else if(p[2] == 0)
        triangleId = getEdgeTriangleD1_xy0(p, localTriangleId);
      else
        triangleId = getEdgeTriangleD1_xyN(p, localTriangleId);
    }
    // D2
    else if(edgeId < esetshift_[4]) {
      edgeToPosition(edgeId, 4, p);

      if(p[0] > 0 and p[0] < nbvoxels_[0])
        triangleId = getEdgeTriangleD2_nyz(p, localTriangleId);
      else if(p[0] == 0)
        triangleId = getEdgeTriangleD2_0yz(p, localTriangleId);
      else
        triangleId = getEdgeTriangleD2_Nyz(p, localTriangleId);
    }
    // D3
    else if(edgeId < esetshift_[5]) {
      edgeToPosition(edgeId, 5, p);

      if(p[1] > 0 and p[1] < nbvoxels_[1])
        triangleId = getEdgeTriangleD3_xnz(p, localTriangleId);
      else if(p[1] == 0)
        triangleId = getEdgeTriangleD3_x0z(p, localTriangleId);
      else
        triangleId = getEdgeTriangleD3_xNz(p, localTriangleId);
    }
    // D4
    else if(edgeId < esetshift_[6]) {
      edgeToPosition(edgeId, 6, p);

      triangleId = getEdgeTriangleD4_xyz(p, localTriangleId);
    }
  } else if(dimensionality_ == 2) {
    SimplexId p[2];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition2d(edgeId, 0, p);

      if(p[1] > 0 and p[1] < nbvoxels_[Dj_])
        triangleId = getEdgeTriangleL_xn(p, localTriangleId);
      else if(p[1] == 0)
        triangleId = getEdgeTriangleL_x0(p, localTriangleId);
      else
        triangleId = getEdgeTriangleL_xN(p, localTriangleId);
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition2d(edgeId, 1, p);

      if(p[0] > 0 and p[0] < nbvoxels_[Di_])
        triangleId = getEdgeTriangleH_ny(p, localTriangleId);
      else if(p[0] == 0)
        triangleId = getEdgeTriangleH_0y(p, localTriangleId);
      else
        triangleId = getEdgeTriangleH_Ny(p, localTriangleId);
    }
    // D1
    else if(edgeId < esetshift_[2]) {
      edgeToPosition2d(edgeId, 2, p);

      triangleId = getEdgeTriangleD1_xy(p, localTriangleId);
    }
  }

  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getEdgeTriangles() {
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
      msg << "[ImplicitTriangulation] Triangle edges built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &edgeTriangleList_;
}

inline SimplexId
  ImplicitTriangulation::getEdgeLinkNumber(const SimplexId &edgeId) const {
  return getEdgeStarNumber(edgeId);
}

int ImplicitTriangulation::getEdgeLink(const SimplexId &edgeId,
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

const vector<vector<SimplexId>> *ImplicitTriangulation::getEdgeLinks() {
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
      msg << "[ImplicitTriangulation] List of edge links built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &edgeLinkList_;
}

inline SimplexId
  ImplicitTriangulation::getEdgeStarNumber(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    SimplexId p[3];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition(edgeId, 0, p);

      if(p[2] > 0 and p[2] < nbvoxels_[2]) {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 6; // ABCG,ABEG,BCDG,BEFG,BFGH,BDGH
        else
          return 3; // BCDG,BFGH,BDGH or ABCG,ABEG,BEFG
      } else if(p[2] == 0) {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 3; // ABCG,ABEG,BCDG
        else if(p[1] == 0)
          return 2; // ABCG,ABEG
        else
          return 1; // BCDG
      } else {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 3; // BEFG,BFGH,BDGH
        else if(p[1] == 0)
          return 1; // BEFG
        else
          return 2; // BFGH,BDGH
      }
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition(edgeId, 1, p);

      if(p[0] > 0 and p[0] < nbvoxels_[0]) {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          return 6; // ABCG,ABEG,BEFG,BFGH,BCDG,BDGH
        else if(p[2] == 0)
          return 3; // BCDG,BDGH,ABCG
        else
          return 3; // ABEG,BEFG,BFGH
      } else if(p[0] == 0) {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          return 3; // ABEG,BEFG,ABCG
        else if(p[2] == 0)
          return 1; // ABCG
        else
          return 2; // ABEG,BEFG
      } else {
        if(p[2] > 0 and p[2] < nbvoxels_[2])
          return 3; // BCDG,BDGH,BFGH
        else if(p[2] == 0)
          return 2; // BCDG,BDGH
        else
          return 1; // BFGH
      }
    }
    // P
    else if(edgeId < esetshift_[2]) {
      edgeToPosition(edgeId, 2, p);

      if(p[0] > 0 and p[0] < nbvoxels_[0]) {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 6; // BDGH,ABCG,BCDG,ABEG,BEFG,BFGH
        else if(p[1] == 0)
          return 3; // BEFG,BFGH,ABEG
        else
          return 3; // ABCG,BCDG,BDGH
      } else if(p[0] == 0) {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 3; // ABCG,BCDG,ABEG
        else if(p[1] == 0)
          return 1; // ABEG
        else
          return 2; // ABCG,BCDG
      } else {
        if(p[1] > 0 and p[1] < nbvoxels_[1])
          return 3; // BEFG,BFGH,BDGH
        else if(p[1] == 0)
          return 2; // BEFG,BFGH
        else
          return 1; // BDGH
      }
    }
    // D1
    else if(edgeId < esetshift_[3]) {
      edgeToPosition(edgeId, 3, p);

      if(p[2] > 0 and p[2] < nbvoxels_[2])
        return 4; // ABCG,BCDG,BEFG,BFGH
      else
        return 2; // ABCG,BCDG ou BEFG,BFGH
    }
    // D2
    else if(edgeId < esetshift_[4]) {
      edgeToPosition(edgeId, 4, p);

      if(p[0] > 0 and p[0] < nbvoxels_[0])
        return 4; // ABCG,ABEG,BDGH,BFGH
      else
        return 2; // ABCG,ABEG ou BDGH,BFGH
    }
    // D3
    else if(edgeId < esetshift_[5]) {
      edgeToPosition(edgeId, 5, p);

      if(p[1] > 0 and p[1] < nbvoxels_[1])
        return 4; // ABEG,BEFG,BCDG,BDGH
      else
        return 2; // ABEG,BEFG ou BCDG,BDGH
    }
    // D4
    else if(edgeId < esetshift_[6])
      return 6;
  } else if(dimensionality_ == 2) {
    SimplexId p[2];

    // L
    if(edgeId < esetshift_[0]) {
      edgeToPosition2d(edgeId, 0, p);

      if(p[1] > 0 and p[1] < nbvoxels_[Dj_])
        return 2;
      else
        return 1;
    }
    // H
    else if(edgeId < esetshift_[1]) {
      edgeToPosition2d(edgeId, 1, p);

      if(p[0] > 0 and p[0] < nbvoxels_[Di_])
        return 2;
      else
        return 1;
    }
    // D1
    else if(edgeId < esetshift_[2])
      return 2;
  }

  return 0;
}

int ImplicitTriangulation::getEdgeStar(const SimplexId &edgeId,
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

const vector<vector<SimplexId>> *ImplicitTriangulation::getEdgeStars() {
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
      msg << "[ImplicitTriangulation] List of edge stars built in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &edgeStarList_;
}

int ImplicitTriangulation::getTriangleVertex(const SimplexId &triangleId,
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

    if(id == 0) {
      switch(localVertexId) {
        case 0:
          vertexId = p[0] / 2 + p[1] * vshift_[0];
          break;
        case 1:
          vertexId = p[0] / 2 + p[1] * vshift_[0] + 1;
          break;
        case 2:
          vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0];
          break;
      }
    } else {
      switch(localVertexId) {
        case 0:
          vertexId = p[0] / 2 + p[1] * vshift_[0] + 1;
          break;
        case 1:
          vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0] + 1;
          break;
        case 2:
          vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0];
          break;
      }
    }
  }

  return 0;
}

int ImplicitTriangulation::getTriangleEdge(const SimplexId &triangleId,
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
          edgeId = p[0] / 2 + (p[1] + 1) * eshift_[0];
          break;
        case 1:
          edgeId = esetshift_[0] + (p[0] + 1) / 2 + p[1] * eshift_[2];
          break;
        case 2:
          edgeId = esetshift_[1] + p[0] / 2 + p[1] * eshift_[4];
          break;
      }
    }
  }

  return 0;
}

int ImplicitTriangulation::getTriangleEdges(
  vector<vector<SimplexId>> &edges) const {
  edges.resize(triangleNumber_);
  for(SimplexId i = 0; i < triangleNumber_; ++i) {
    edges[i].resize(3);
    for(int j = 0; j < 3; ++j)
      getTriangleEdge(i, j, edges[i][j]);
  }
  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getTriangleEdges() {
  if(!triangleEdgeList_.size()) {
    Timer t;

    getTriangleEdges(triangleEdgeList_);

    {
      stringstream msg;
      msg << "[ImplicitTriangulation] Triangle edges (" << triangleNumber_
          << " triangle(s), " << edgeNumber_ << " edge(s)) computed in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &triangleEdgeList_;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getTriangles() {
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
      msg << "[ImplicitTriangulation] Triangle list (" << triangleNumber_
          << " triangles) computed in " << t.getElapsedTime() << " s. (" << 1
          << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &triangleList_;
}

int ImplicitTriangulation::getTriangleLink(const SimplexId &triangleId,
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

inline SimplexId ImplicitTriangulation::getTriangleLinkNumber(
  const SimplexId &triangleId) const {
  return getTriangleStarNumber(triangleId);
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getTriangleLinks() {
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

inline SimplexId ImplicitTriangulation::getTriangleStarNumber(
  const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    SimplexId p[3];

    // F
    if(triangleId < tsetshift_[0]) {
      triangleToPosition(triangleId, 0, p);
      if(p[2] > 0 and p[2] < nbvoxels_[2])
        return 2;
      else
        return 1;
    }
    // H
    else if(triangleId < tsetshift_[1]) {
      triangleToPosition(triangleId, 1, p);
      if(p[1] > 0 and p[1] < nbvoxels_[1])
        return 2;
      else
        return 1;
    }
    // C
    else if(triangleId < tsetshift_[2]) {
      triangleToPosition(triangleId, 2, p);
      if(p[0] < 2 or p[0] >= (dimensions_[0] * 2 - 2))
        return 1;
      else
        return 2;
    }
    // D1
    else if(triangleId < tsetshift_[3]) {
      return 2;
    }
    // D2
    else if(triangleId < tsetshift_[4]) {
      triangleToPosition(triangleId, 4, p);
      return 2;
    }
    // D3
    else if(triangleId < tsetshift_[5]) {
      triangleToPosition(triangleId, 5, p);
      return 2;
    }
  }
  return 0;
}

int ImplicitTriangulation::getTriangleStar(const SimplexId &triangleId,
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

const vector<vector<SimplexId>> *ImplicitTriangulation::getTriangleStars() {
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

inline SimplexId ImplicitTriangulation::getTriangleNeighborNumber(
  const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  if(dimensionality_ == 2) {
    SimplexId p[2];
    triangleToPosition2d(triangleId, p);
    const SimplexId id = triangleId % 2;

    if(id) {
      if(p[0] / 2 == nbvoxels_[Di_] - 1 and p[1] == nbvoxels_[Dj_] - 1)
        return 1;
      else if(p[0] / 2 == nbvoxels_[Di_] - 1 or p[1] == nbvoxels_[Dj_] - 1)
        return 2;
      else
        return 3;
    } else {
      if(p[0] == 0 and p[1] == 0)
        return 1;
      else if(p[0] == 0 or p[1] == 0)
        return 2;
      else
        return 3;
    }
  }

  return 0;
}

int ImplicitTriangulation::getTriangleNeighbor(const SimplexId &triangleId,
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
      if(p[0] / 2 == nbvoxels_[Di_] - 1 and p[1] == nbvoxels_[Dj_] - 1)
        neighborId = triangleId - 1;
      else if(p[0] / 2 == nbvoxels_[Di_] - 1) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId - 1;
            break;
          case 1:
            neighborId = triangleId + tshift_[0] - 1;
            break;
        }
      } else if(p[1] == nbvoxels_[Dj_] - 1) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId - 1;
            break;
          case 1:
            neighborId = triangleId + 1;
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
      if(p[0] == 0 and p[1] == 0)
        neighborId = triangleId + 1;
      else if(p[0] == 0) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId + 1;
            break;
          case 1:
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

int ImplicitTriangulation::getTriangleNeighbors(
  vector<vector<SimplexId>> &neighbors) {
  neighbors.resize(triangleNumber_);
  for(SimplexId i = 0; i < triangleNumber_; ++i) {
    neighbors[i].resize(getTriangleNeighborNumber(i));
    for(SimplexId j = 0; j < (SimplexId)neighbors[i].size(); ++j)
      getTriangleNeighbor(i, j, neighbors[i][j]);
  }
  return 0;
}

int ImplicitTriangulation::getTetrahedronVertex(const SimplexId &tetId,
                                                const int &localVertexId,
                                                SimplexId &vertexId) const {
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

int ImplicitTriangulation::getTetrahedronEdge(const SimplexId &tetId,
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

int ImplicitTriangulation::getTetrahedronEdges(
  vector<vector<SimplexId>> &edges) const {
  edges.resize(tetrahedronNumber_);
  for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    edges[i].resize(6);
    for(int j = 0; j < 6; ++j)
      getTetrahedronEdge(i, j, edges[i][j]);
  }

  return 0;
}

int ImplicitTriangulation::getTetrahedronTriangle(const SimplexId &tetId,
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

int ImplicitTriangulation::getTetrahedronTriangles(
  vector<vector<SimplexId>> &triangles) const {
  triangles.resize(tetrahedronNumber_);
  for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    triangles[i].resize(4);
    for(int j = 0; j < 4; ++j)
      getTetrahedronTriangle(i, j, triangles[i][j]);
  }

  return 0;
}

SimplexId ImplicitTriangulation::getTetrahedronNeighborNumber(
  const SimplexId &tetId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
#endif

  if(dimensionality_ == 3) {
    const SimplexId id = tetId % 6;
    SimplexId p[3];
    tetrahedronToPosition(tetId, p);

    switch(id) {
      case 0: // ABCG
        if(p[0] == 0 and p[2] == 0)
          return 2;
        else if(p[0] == 0 or p[2] == 0)
          return 3;
        else
          return 4;
        break;
      case 1: // BCDG
        if(p[1] == nbvoxels_[1] - 1 and p[2] == 0)
          return 2;
        else if(p[1] == nbvoxels_[1] - 1 or p[2] == 0)
          return 3;
        else
          return 4;
        break;
      case 2: // ABEG
        if(p[0] == 0 and p[1] == 0)
          return 2;
        else if(p[0] == 0 or p[1] == 0)
          return 3;
        else
          return 4;
        break;
      case 3: // BEFG
        if(p[1] == 0 and p[2] == nbvoxels_[2] - 1)
          return 2;
        else if(p[1] == 0 or p[2] == nbvoxels_[2] - 1)
          return 3;
        else
          return 4;
        break;
      case 4: // BFGH
        if(p[0] == nbvoxels_[0] - 1 and p[2] == nbvoxels_[2] - 1)
          return 2;
        else if(p[0] == nbvoxels_[0] - 1 or p[2] == nbvoxels_[2] - 1)
          return 3;
        else
          return 4;
        break;
      case 5: // BDGH
        if(p[0] == nbvoxels_[0] - 1 and p[1] == nbvoxels_[1] - 1)
          return 2;
        else if(p[0] == nbvoxels_[0] - 1 or p[1] == nbvoxels_[1] - 1)
          return 3;
        else
          return 4;
        break;
    }
  }

  return 0;
}

int ImplicitTriangulation::getTetrahedronNeighbor(const SimplexId &tetId,
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

int ImplicitTriangulation::getTetrahedronNeighbors(
  vector<vector<SimplexId>> &neighbors) {
  neighbors.resize(tetrahedronNumber_);
  for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    neighbors[i].resize(getTetrahedronNeighborNumber(i));
    for(SimplexId j = 0; j < (SimplexId)neighbors[i].size(); ++j)
      getTetrahedronNeighbor(i, j, neighbors[i][j]);
  }

  return 0;
}

SimplexId
  ImplicitTriangulation::getCellVertexNumber(const SimplexId &cellId) const {
  return dimensionality_ + 1;
}

int ImplicitTriangulation::getCellVertex(const SimplexId &cellId,
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

SimplexId
  ImplicitTriangulation::getCellEdgeNumber(const SimplexId &cellId) const {
  if(dimensionality_ == 3)
    return 6;
  else if(dimensionality_ == 2)
    return 3;

  return 0;
}

int ImplicitTriangulation::getCellEdge(const SimplexId &cellId,
                                       const int &localEdgeId,
                                       SimplexId &edgeId) const {
  if(dimensionality_ == 3)
    getTetrahedronEdge(cellId, localEdgeId, edgeId);
  else if(dimensionality_ == 2)
    getTriangleEdge(cellId, localEdgeId, edgeId);

  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getCellEdges() {
  if(!cellEdgeList_.size()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronEdges(cellEdgeList_);
    else if(dimensionality_ == 2)
      getTriangleEdges(cellEdgeList_);

    {
      stringstream msg;
      msg << "[ImplicitTriangulation] Cell edges (" << getNumberOfCells()
          << " cell(s), " << edgeNumber_ << "edge(s)) computed in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &cellEdgeList_;
}

int ImplicitTriangulation::getCellTriangle(const SimplexId &cellId,
                                           const int &localTriangleId,
                                           SimplexId &triangleId) const {
  if(dimensionality_ == 3)
    getTetrahedronTriangle(cellId, localTriangleId, triangleId);

  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getCellTriangles() {
  if(!cellTriangleList_.size()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronTriangles(cellTriangleList_);

    {
      stringstream msg;
      msg << "[ImplicitTriangulation] Cell triangles (" << cellNumber_
          << " cell(s), " << triangleNumber_ << "edge(s)) computed in "
          << t.getElapsedTime() << " s. (" << 1 << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &cellTriangleList_;
}

SimplexId
  ImplicitTriangulation::getCellNeighborNumber(const SimplexId &cellId) const {
  if(dimensionality_ == 3)
    return getTetrahedronNeighborNumber(cellId);
  else if(dimensionality_ == 2)
    return getTriangleNeighborNumber(cellId);
  else if(dimensionality_ == 1) {
    stringstream msg;
    msg << "[ImplicitTriangulation] getCellNeighborNumber() in 1D:" << endl;
    msg << "[ImplicitTriangulation] Not implemented! TODO!" << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
    return -1;
  }

  return 0;
}

int ImplicitTriangulation::getCellNeighbor(const SimplexId &cellId,
                                           const int &localNeighborId,
                                           SimplexId &neighborId) const {
  if(dimensionality_ == 3)
    getTetrahedronNeighbor(cellId, localNeighborId, neighborId);
  else if(dimensionality_ == 2)
    getTriangleNeighbor(cellId, localNeighborId, neighborId);
  else if(dimensionality_ == 1) {
    stringstream msg;
    msg << "[ImplicitTriangulation] getCellNeighbor() in 1D:" << endl;
    msg << "[ImplicitTriangulation] Not implemented! TODO!" << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
    return -1;
  }

  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getCellNeighbors() {
  if(!cellNeighborList_.size()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronNeighbors(cellNeighborList_);
    else if(dimensionality_ == 2)
      getTriangleNeighbors(cellNeighborList_);
    else if(dimensionality_ == 1) {
      stringstream msg;
      msg << "[ImplicitTriangulation] getCellNeighbors() in 1D:" << endl;
      msg << "[ImplicitTriangulation] Not implemented! TODO!" << endl;
      dMsg(cerr, msg.str(), Debug::fatalMsg);
      return nullptr;
    }

    {
      stringstream msg;
      msg << "[ImplicitTriangulation] Cell neighbors (" << getNumberOfCells()
          << " cells) computed in " << t.getElapsedTime() << " s. (" << 1
          << " thread(s))." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  return &cellNeighborList_;
}
