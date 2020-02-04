#include <ImplicitTriangulation.h>

using namespace std;
using namespace ttk;

#define CASE_EDGE_POSITION_L_3D \
  case EdgePosition::L_xnn_3D:  \
  case EdgePosition::L_xn0_3D:  \
  case EdgePosition::L_xnN_3D:  \
  case EdgePosition::L_x0n_3D:  \
  case EdgePosition::L_x00_3D:  \
  case EdgePosition::L_x0N_3D:  \
  case EdgePosition::L_xNn_3D:  \
  case EdgePosition::L_xN0_3D:  \
  case EdgePosition::L_xNN_3D
#define CASE_EDGE_POSITION_H_3D \
  case EdgePosition::H_nyn_3D:  \
  case EdgePosition::H_ny0_3D:  \
  case EdgePosition::H_nyN_3D:  \
  case EdgePosition::H_0yn_3D:  \
  case EdgePosition::H_0y0_3D:  \
  case EdgePosition::H_0yN_3D:  \
  case EdgePosition::H_Nyn_3D:  \
  case EdgePosition::H_Ny0_3D:  \
  case EdgePosition::H_NyN_3D
#define CASE_EDGE_POSITION_P_3D \
  case EdgePosition::P_nnz_3D:  \
  case EdgePosition::P_n0z_3D:  \
  case EdgePosition::P_nNz_3D:  \
  case EdgePosition::P_0nz_3D:  \
  case EdgePosition::P_00z_3D:  \
  case EdgePosition::P_0Nz_3D:  \
  case EdgePosition::P_Nnz_3D:  \
  case EdgePosition::P_N0z_3D:  \
  case EdgePosition::P_NNz_3D
#define CASE_EDGE_POSITION_D1_3D \
  case EdgePosition::D1_xyn_3D:  \
  case EdgePosition::D1_xy0_3D:  \
  case EdgePosition::D1_xyN_3D
#define CASE_EDGE_POSITION_D2_3D \
  case EdgePosition::D2_nyz_3D:  \
  case EdgePosition::D2_0yz_3D:  \
  case EdgePosition::D2_Nyz_3D
#define CASE_EDGE_POSITION_D3_3D \
  case EdgePosition::D3_xnz_3D:  \
  case EdgePosition::D3_x0z_3D:  \
  case EdgePosition::D3_xNz_3D
#define CASE_EDGE_POSITION_L_2D \
  case EdgePosition::L_xn_2D:   \
  case EdgePosition::L_x0_2D:   \
  case EdgePosition::L_xN_2D
#define CASE_EDGE_POSITION_H_2D \
  case EdgePosition::H_ny_2D:   \
  case EdgePosition::H_0y_2D:   \
  case EdgePosition::H_Ny_2D

ImplicitTriangulation::ImplicitTriangulation()
  : dimensionality_{-1}, cellNumber_{}, vertexNumber_{}, edgeNumber_{},
    triangleNumber_{}, tetrahedronNumber_{}, isAccelerated_{} {
  setDebugMsgPrefix("ImplicitTriangulation");
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
    printMsg("Accelerated getVertex*() requests.", debug::Priority::INFO);
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

bool ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(isVertexOnBoundary)(
  const SimplexId &vertexId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return false;
#endif // !TTK_ENABLE_KAMIKAZE

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
    case VertexPosition::CENTER_2D:
    case VertexPosition::CENTER_1D:
      return false;
    default:
      return true;
  }
}

bool ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(isEdgeOnBoundary)(
  const SimplexId &edgeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return false;
#endif // !TTK_ENABLE_KAMIKAZE

  switch(edgePositions_[edgeId]) {
    case EdgePosition::L_xnn_3D:
    case EdgePosition::H_nyn_3D:
    case EdgePosition::P_nnz_3D:
    case EdgePosition::D1_xyn_3D:
    case EdgePosition::D2_nyz_3D:
    case EdgePosition::D3_xnz_3D:
    case EdgePosition::D4_3D:
    case EdgePosition::L_xn_2D:
    case EdgePosition::H_ny_2D:
    case EdgePosition::D1_2D:
      return false;
    default:
      break;
  }
  return true;
}

bool ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(isTriangleOnBoundary)(
  const SimplexId &triangleId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return false;
#endif // !TTK_ENABLE_KAMIKAZE

  if(dimensionality_ == 3)
    return (getTriangleStarNumber(triangleId) == 1);

  return false;
}

int ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexNeighbor)(
  const SimplexId &vertexId,
  const int &localNeighborId,
  SimplexId &neighborId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= getVertexNeighborNumber(vertexId))
    return -1;
#endif // !TTK_ENABLE_KAMIKAZE

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      neighborId = getVertexNeighborABCDEFGH(vertexId, localNeighborId);
      break;
    case VertexPosition::FRONT_FACE_3D:
      neighborId = getVertexNeighborABDC(vertexId, localNeighborId);
      break;
    case VertexPosition::BACK_FACE_3D:
      neighborId = getVertexNeighborEFHG(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_FACE_3D:
      neighborId = getVertexNeighborAEFB(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_FACE_3D:
      neighborId = getVertexNeighborGHDC(vertexId, localNeighborId);
      break;
    case VertexPosition::LEFT_FACE_3D:
      neighborId = getVertexNeighborAEGC(vertexId, localNeighborId);
      break;
    case VertexPosition::RIGHT_FACE_3D:
      neighborId = getVertexNeighborBFHD(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_FRONT_EDGE_3D: // ab
      neighborId = getVertexNeighborAB(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      neighborId = getVertexNeighborCD(vertexId, localNeighborId);
      break;
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      neighborId = getVertexNeighborAC(vertexId, localNeighborId);
      break;
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      neighborId = getVertexNeighborBD(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_BACK_EDGE_3D: // ef
      neighborId = getVertexNeighborEF(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      neighborId = getVertexNeighborGH(vertexId, localNeighborId);
      break;
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
      neighborId = getVertexNeighborEG(vertexId, localNeighborId);
      break;
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      neighborId = getVertexNeighborFH(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
      neighborId = getVertexNeighborAE(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      neighborId = getVertexNeighborBF(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      neighborId = getVertexNeighborCG(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      neighborId = getVertexNeighborDH(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      neighborId = getVertexNeighborA(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      neighborId = getVertexNeighborB(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      neighborId = getVertexNeighborC(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      neighborId = getVertexNeighborD(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      neighborId = getVertexNeighborE(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      neighborId = getVertexNeighborF(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      neighborId = getVertexNeighborG(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      neighborId = getVertexNeighborH(vertexId, localNeighborId);
      break;
    case VertexPosition::CENTER_2D:
      neighborId = getVertexNeighbor2dABCD(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_EDGE_2D:
      neighborId = getVertexNeighbor2dAB(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_EDGE_2D:
      neighborId = getVertexNeighbor2dCD(vertexId, localNeighborId);
      break;
    case VertexPosition::LEFT_EDGE_2D:
      neighborId = getVertexNeighbor2dAC(vertexId, localNeighborId);
      break;
    case VertexPosition::RIGHT_EDGE_2D:
      neighborId = getVertexNeighbor2dBD(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_LEFT_CORNER_2D: // a
      neighborId = getVertexNeighbor2dA(vertexId, localNeighborId);
      break;
    case VertexPosition::TOP_RIGHT_CORNER_2D: // b
      neighborId = getVertexNeighbor2dB(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      neighborId = getVertexNeighbor2dC(vertexId, localNeighborId);
      break;
    case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      neighborId = getVertexNeighbor2dD(vertexId, localNeighborId);
      break;
    case VertexPosition::CENTER_1D:
      neighborId = (localNeighborId == 0 ? vertexId + 1 : vertexId - 1);
      break;
    case VertexPosition::LEFT_CORNER_1D:
      neighborId = vertexId + 1;
      break;
    case VertexPosition::RIGHT_CORNER_1D:
      neighborId = vertexId - 1;
      break;
    default:
      neighborId = -1;
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() {
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

SimplexId ImplicitTriangulation::getVertexEdgeNumberInternal(
  const SimplexId &vertexId) const {
  return getVertexNeighborNumber(vertexId);
}

int ImplicitTriangulation::getVertexEdgeInternal(const SimplexId &vertexId,
                                                 const int &localEdgeId,
                                                 SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localEdgeId < 0 or localEdgeId >= getVertexEdgeNumberInternal(vertexId))
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

  const auto &p = vertexCoords_[vertexId];

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      edgeId = getVertexEdgeABCDEFGH(p.data(), localEdgeId);
      break;
    case VertexPosition::FRONT_FACE_3D:
      edgeId = getVertexEdgeABDC(p.data(), localEdgeId);
      break;
    case VertexPosition::BACK_FACE_3D:
      edgeId = getVertexEdgeEFHG(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_FACE_3D:
      edgeId = getVertexEdgeAEFB(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_FACE_3D:
      edgeId = getVertexEdgeGHDC(p.data(), localEdgeId);
      break;
    case VertexPosition::LEFT_FACE_3D:
      edgeId = getVertexEdgeAEGC(p.data(), localEdgeId);
      break;
    case VertexPosition::RIGHT_FACE_3D:
      edgeId = getVertexEdgeBFHD(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_FRONT_EDGE_3D: // ab
      edgeId = getVertexEdgeAB(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      edgeId = getVertexEdgeCD(p.data(), localEdgeId);
      break;
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      edgeId = getVertexEdgeAC(p.data(), localEdgeId);
      break;
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      edgeId = getVertexEdgeBD(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_BACK_EDGE_3D: // ef
      edgeId = getVertexEdgeEF(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      edgeId = getVertexEdgeGH(p.data(), localEdgeId);
      break;
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
      edgeId = getVertexEdgeEG(p.data(), localEdgeId);
      break;
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      edgeId = getVertexEdgeFH(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
      edgeId = getVertexEdgeAE(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      edgeId = getVertexEdgeBF(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      edgeId = getVertexEdgeCG(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      edgeId = getVertexEdgeDH(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      edgeId = getVertexEdgeA(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      edgeId = getVertexEdgeB(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      edgeId = getVertexEdgeC(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      edgeId = getVertexEdgeD(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      edgeId = getVertexEdgeE(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      edgeId = getVertexEdgeF(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      edgeId = getVertexEdgeG(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      edgeId = getVertexEdgeH(p.data(), localEdgeId);
      break;
    case VertexPosition::CENTER_2D:
      edgeId = getVertexEdge2dABCD(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_EDGE_2D:
      edgeId = getVertexEdge2dAB(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_EDGE_2D:
      edgeId = getVertexEdge2dCD(p.data(), localEdgeId);
      break;
    case VertexPosition::LEFT_EDGE_2D:
      edgeId = getVertexEdge2dAC(p.data(), localEdgeId);
      break;
    case VertexPosition::RIGHT_EDGE_2D:
      edgeId = getVertexEdge2dBD(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_LEFT_CORNER_2D: // a
      edgeId = getVertexEdge2dA(p.data(), localEdgeId);
      break;
    case VertexPosition::TOP_RIGHT_CORNER_2D: // b
      edgeId = getVertexEdge2dB(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      edgeId = getVertexEdge2dC(p.data(), localEdgeId);
      break;
    case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      edgeId = getVertexEdge2dD(p.data(), localEdgeId);
      break;
    case VertexPosition::CENTER_1D:
      edgeId = (localEdgeId == 0 ? vertexId : vertexId - 1);
      break;
    case VertexPosition::LEFT_CORNER_1D:
      edgeId = vertexId;
      break;
    case VertexPosition::RIGHT_CORNER_1D:
      edgeId = vertexId - 1;
      break;
    default:
      edgeId = -1;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::getVertexEdgesInternal() {
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

inline SimplexId ImplicitTriangulation::getVertexTriangleNumberInternal(
  const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      return 36;
    case VertexPosition::FRONT_FACE_3D:
    case VertexPosition::BACK_FACE_3D:
    case VertexPosition::TOP_FACE_3D:
    case VertexPosition::BOTTOM_FACE_3D:
    case VertexPosition::LEFT_FACE_3D:
    case VertexPosition::RIGHT_FACE_3D:
      return 21;
    case VertexPosition::TOP_FRONT_EDGE_3D: // ab
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      return 15;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      return 12;
    case VertexPosition::TOP_BACK_EDGE_3D: // ef
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      return 9;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      return 5;
    default: // 1D + 2D
      break;
  }

  return 0;
}

int ImplicitTriangulation::getVertexTriangleInternal(
  const SimplexId &vertexId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= getVertexTriangleNumberInternal(vertexId))
    return -1;
#endif

  const auto &p = vertexCoords_[vertexId];

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      triangleId = getVertexTriangleABCDEFGH(p.data(), localTriangleId);
      break;
    case VertexPosition::FRONT_FACE_3D:
      triangleId = getVertexTriangleABDC(p.data(), localTriangleId);
      break;
    case VertexPosition::BACK_FACE_3D:
      triangleId = getVertexTriangleEFHG(p.data(), localTriangleId);
      break;
    case VertexPosition::TOP_FACE_3D:
      triangleId = getVertexTriangleAEFB(p.data(), localTriangleId);
      break;
    case VertexPosition::BOTTOM_FACE_3D:
      triangleId = getVertexTriangleGHDC(p.data(), localTriangleId);
      break;
    case VertexPosition::LEFT_FACE_3D:
      triangleId = getVertexTriangleAEGC(p.data(), localTriangleId);
      break;
    case VertexPosition::RIGHT_FACE_3D:
      triangleId = getVertexTriangleBFHD(p.data(), localTriangleId);
      break;
    case VertexPosition::TOP_FRONT_EDGE_3D: // ab
      triangleId = getVertexTriangleAB(p.data(), localTriangleId);
      break;
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      triangleId = getVertexTriangleCD(p.data(), localTriangleId);
      break;
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      triangleId = getVertexTriangleAC(p.data(), localTriangleId);
      break;
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      triangleId = getVertexTriangleBD(p.data(), localTriangleId);
      break;
    case VertexPosition::TOP_BACK_EDGE_3D: // ef
      triangleId = getVertexTriangleEF(p.data(), localTriangleId);
      break;
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      triangleId = getVertexTriangleGH(p.data(), localTriangleId);
      break;
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
      triangleId = getVertexTriangleEG(p.data(), localTriangleId);
      break;
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      triangleId = getVertexTriangleFH(p.data(), localTriangleId);
      break;
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
      triangleId = getVertexTriangleAE(p.data(), localTriangleId);
      break;
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      triangleId = getVertexTriangleBF(p.data(), localTriangleId);
      break;
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      triangleId = getVertexTriangleCG(p.data(), localTriangleId);
      break;
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      triangleId = getVertexTriangleDH(p.data(), localTriangleId);
      break;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      triangleId = getVertexTriangleA(p.data(), localTriangleId);
      break;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      triangleId = getVertexTriangleB(p.data(), localTriangleId);
      break;
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      triangleId = getVertexTriangleC(p.data(), localTriangleId);
      break;
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      triangleId = getVertexTriangleD(p.data(), localTriangleId);
      break;
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      triangleId = getVertexTriangleE(p.data(), localTriangleId);
      break;
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      triangleId = getVertexTriangleF(p.data(), localTriangleId);
      break;
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      triangleId = getVertexTriangleG(p.data(), localTriangleId);
      break;
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      triangleId = getVertexTriangleH(p.data(), localTriangleId);
      break;
    default: // 1D + 2D
      triangleId = -1;
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::getVertexTrianglesInternal() {
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

SimplexId ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getVertexLinkNumber)(const SimplexId &vertexId) const {
  return getVertexStarNumber(vertexId);
}

int ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexLink)(
  const SimplexId &vertexId, const int &localLinkId, SimplexId &linkId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getVertexLinkNumber(vertexId))
    return -1;
#endif // !TTK_ENABLE_KAMIKAZE

  const auto &p = vertexCoords_[vertexId];

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      linkId = getVertexLinkABCDEFGH(p.data(), localLinkId);
      break;
    case VertexPosition::FRONT_FACE_3D:
      linkId = getVertexLinkABDC(p.data(), localLinkId);
      break;
    case VertexPosition::BACK_FACE_3D:
      linkId = getVertexLinkEFHG(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_FACE_3D:
      linkId = getVertexLinkAEFB(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_FACE_3D:
      linkId = getVertexLinkGHDC(p.data(), localLinkId);
      break;
    case VertexPosition::LEFT_FACE_3D:
      linkId = getVertexLinkAEGC(p.data(), localLinkId);
      break;
    case VertexPosition::RIGHT_FACE_3D:
      linkId = getVertexLinkBFHD(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_FRONT_EDGE_3D: // ab
      linkId = getVertexLinkAB(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      linkId = getVertexLinkCD(p.data(), localLinkId);
      break;
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      linkId = getVertexLinkAC(p.data(), localLinkId);
      break;
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      linkId = getVertexLinkBD(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_BACK_EDGE_3D: // ef
      linkId = getVertexLinkEF(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      linkId = getVertexLinkGH(p.data(), localLinkId);
      break;
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
      linkId = getVertexLinkEG(p.data(), localLinkId);
      break;
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      linkId = getVertexLinkFH(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
      linkId = getVertexLinkAE(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      linkId = getVertexLinkBF(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      linkId = getVertexLinkCG(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      linkId = getVertexLinkDH(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      linkId = getVertexLinkA(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      linkId = getVertexLinkB(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      linkId = getVertexLinkC(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      linkId = getVertexLinkD(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      linkId = getVertexLinkE(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      linkId = getVertexLinkF(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      linkId = getVertexLinkG(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      linkId = getVertexLinkH(p.data(), localLinkId);
      break;
    case VertexPosition::CENTER_2D:
      linkId = getVertexLink2dABCD(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_EDGE_2D:
      linkId = getVertexLink2dAB(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_EDGE_2D:
      linkId = getVertexLink2dCD(p.data(), localLinkId);
      break;
    case VertexPosition::LEFT_EDGE_2D:
      linkId = getVertexLink2dAC(p.data(), localLinkId);
      break;
    case VertexPosition::RIGHT_EDGE_2D:
      linkId = getVertexLink2dBD(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_LEFT_CORNER_2D: // a
      linkId = getVertexLink2dA(p.data(), localLinkId);
      break;
    case VertexPosition::TOP_RIGHT_CORNER_2D: // b
      linkId = getVertexLink2dB(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      linkId = getVertexLink2dC(p.data(), localLinkId);
      break;
    case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      linkId = getVertexLink2dD(p.data(), localLinkId);
      break;
    default: // 1D
      linkId = -1;
      break;
  };

  return 0;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexLinks)() {
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

inline SimplexId ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getVertexStarNumber)(const SimplexId &vertexId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif // !TTK_ENABLE_KAMIKAZE

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      return 24;
    case VertexPosition::FRONT_FACE_3D:
    case VertexPosition::BACK_FACE_3D:
    case VertexPosition::TOP_FACE_3D:
    case VertexPosition::BOTTOM_FACE_3D:
    case VertexPosition::LEFT_FACE_3D:
    case VertexPosition::RIGHT_FACE_3D:
      return 12;
    case VertexPosition::TOP_FRONT_EDGE_3D: // ab
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      return 8;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
    case VertexPosition::CENTER_2D:
      return 6;
    case VertexPosition::TOP_BACK_EDGE_3D: // ef
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      return 4;
    case VertexPosition::TOP_EDGE_2D: // ab
    case VertexPosition::BOTTOM_EDGE_2D: // cd
    case VertexPosition::LEFT_EDGE_2D: // ac
    case VertexPosition::RIGHT_EDGE_2D: // bd
      return 3;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
    case VertexPosition::TOP_RIGHT_CORNER_2D: // b
    case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      return 2;
    case VertexPosition::TOP_LEFT_CORNER_2D: // a
    case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      return 1;
    default: // 1D
      break;
  }

  return 0;
}

int ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexStar)(
  const SimplexId &vertexId, const int &localStarId, SimplexId &starId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getVertexStarNumber(vertexId))
    return -1;
#endif // !TTK_ENABLE_KAMIKAZE

  const auto &p = vertexCoords_[vertexId];

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      starId = getVertexStarABCDEFGH(p.data(), localStarId);
      break;
    case VertexPosition::FRONT_FACE_3D:
      starId = getVertexStarABDC(p.data(), localStarId);
      break;
    case VertexPosition::BACK_FACE_3D:
      starId = getVertexStarEFHG(p.data(), localStarId);
      break;
    case VertexPosition::TOP_FACE_3D:
      starId = getVertexStarAEFB(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_FACE_3D:
      starId = getVertexStarGHDC(p.data(), localStarId);
      break;
    case VertexPosition::LEFT_FACE_3D:
      starId = getVertexStarAEGC(p.data(), localStarId);
      break;
    case VertexPosition::RIGHT_FACE_3D:
      starId = getVertexStarBFHD(p.data(), localStarId);
      break;
    case VertexPosition::TOP_FRONT_EDGE_3D: // ab
      starId = getVertexStarAB(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      starId = getVertexStarCD(p.data(), localStarId);
      break;
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      starId = getVertexStarAC(p.data(), localStarId);
      break;
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      starId = getVertexStarBD(p.data(), localStarId);
      break;
    case VertexPosition::TOP_BACK_EDGE_3D: // ef
      starId = getVertexStarEF(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      starId = getVertexStarGH(p.data(), localStarId);
      break;
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
      starId = getVertexStarEG(p.data(), localStarId);
      break;
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      starId = getVertexStarFH(p.data(), localStarId);
      break;
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
      starId = getVertexStarAE(p.data(), localStarId);
      break;
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      starId = getVertexStarBF(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      starId = getVertexStarCG(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      starId = getVertexStarDH(p.data(), localStarId);
      break;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      starId = getVertexStarA(p.data(), localStarId);
      break;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      starId = getVertexStarB(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      starId = getVertexStarC(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      starId = getVertexStarD(p.data(), localStarId);
      break;
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      starId = getVertexStarE(p.data(), localStarId);
      break;
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      starId = getVertexStarF(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      starId = getVertexStarG(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      starId = getVertexStarH(p.data(), localStarId);
      break;
    case VertexPosition::CENTER_2D:
      starId = getVertexStar2dABCD(p.data(), localStarId);
      break;
    case VertexPosition::TOP_EDGE_2D:
      starId = getVertexStar2dAB(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_EDGE_2D:
      starId = getVertexStar2dCD(p.data(), localStarId);
      break;
    case VertexPosition::LEFT_EDGE_2D:
      starId = getVertexStar2dAC(p.data(), localStarId);
      break;
    case VertexPosition::RIGHT_EDGE_2D:
      starId = getVertexStar2dBD(p.data(), localStarId);
      break;
    case VertexPosition::TOP_LEFT_CORNER_2D: // a
      starId = getVertexStar2dA(p.data(), localStarId);
      break;
    case VertexPosition::TOP_RIGHT_CORNER_2D: // b
      starId = getVertexStar2dB(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      starId = getVertexStar2dC(p.data(), localStarId);
      break;
    case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      starId = getVertexStar2dD(p.data(), localStarId);
      break;
    default: // 1D
      starId = -1;
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexStars)() {

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

int ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getVertexPoint)(
  const SimplexId &vertexId, float &x, float &y, float &z) const {

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

int ImplicitTriangulation::getEdgeVertexInternal(const SimplexId &edgeId,
                                                 const int &localVertexId,
                                                 SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 2)
    return -2;
#endif

  const auto &p = edgeCoords_[edgeId];

  const auto helper3d = [&](const SimplexId a, const SimplexId b) -> SimplexId {
    if(isAccelerated_) {
      const auto tmp = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]);
      return (localVertexId == 0) ? tmp + a : tmp + b;
    } else {
      const auto tmp = p[0] + (p[1] * vshift_[0]) + (p[2] * vshift_[1]);
      return (localVertexId == 0) ? tmp + a : tmp + b;
    }
  };

  const auto helper2d = [&](const SimplexId a, const SimplexId b) -> SimplexId {
    if(isAccelerated_) {
      const auto tmp = p[0] + (p[1] << div_[0]);
      return localVertexId == 0 ? tmp + a : tmp + b;
    } else {
      const auto tmp = p[0] + (p[1] * vshift_[0]);
      return localVertexId == 0 ? tmp + a : tmp + b;
    }
  };

  switch(edgePositions_[edgeId]) {
  CASE_EDGE_POSITION_L_3D:
    vertexId = helper3d(0, 1);
    break;
  CASE_EDGE_POSITION_H_3D:
    vertexId = helper3d(0, vshift_[0]);
    break;
  CASE_EDGE_POSITION_P_3D:
    vertexId = helper3d(0, vshift_[1]);
    break;
  CASE_EDGE_POSITION_D1_3D:
    vertexId = helper3d(1, vshift_[0]);
    break;
  CASE_EDGE_POSITION_D2_3D:
    vertexId = helper3d(0, vshift_[0] + vshift_[1]);
    break;
  CASE_EDGE_POSITION_D3_3D:
    vertexId = helper3d(1, vshift_[1]);
    break;
    case EdgePosition::D4_3D:
      vertexId = helper3d(1, vshift_[0] + vshift_[1]);
      break;

    CASE_EDGE_POSITION_L_2D:
      vertexId = helper2d(0, 1);
      break;
    CASE_EDGE_POSITION_H_2D:
      vertexId = helper2d(0, vshift_[0]);
      break;
    case EdgePosition::D1_2D:
      vertexId = helper2d(1, vshift_[0]);
      break;

    case EdgePosition::FIRST_EDGE_1D:
      vertexId = localVertexId == 0 ? 0 : 1;
      break;
    case EdgePosition::LAST_EDGE_1D:
      vertexId = localVertexId == 0 ? edgeNumber_ - 1 : edgeNumber_;
      break;
    case EdgePosition::CENTER_1D:
      vertexId = localVertexId == 0 ? edgeId : edgeId + 1;
      break;
  }

  return 0;
}

const vector<pair<SimplexId, SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdges)() {

  if(!edgeList_.size()) {
    Timer t;

    edgeList_.resize(edgeNumber_);
    for(SimplexId i = 0; i < edgeNumber_; ++i) {
      SimplexId id0, id1;
      getEdgeVertexInternal(i, 0, id0);
      getEdgeVertexInternal(i, 1, id1);
      edgeList_[i].first = id0;
      edgeList_[i].second = id1;
    }

    printMsg(
      "Built " + to_string(edgeNumber_) + " edges.", 1, t.getElapsedTime(), 1);
  }

  return &edgeList_;
}

inline SimplexId ImplicitTriangulation::getEdgeTriangleNumberInternal(
  const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  switch(edgePositions_[edgeId]) {
    case EdgePosition::L_xnn_3D:
    case EdgePosition::H_nyn_3D:
    case EdgePosition::P_nnz_3D:
    case EdgePosition::D4_3D:
      return 6;
    case EdgePosition::L_x0n_3D:
    case EdgePosition::L_xNn_3D:
    case EdgePosition::L_xn0_3D:
    case EdgePosition::L_xnN_3D:
    case EdgePosition::H_ny0_3D:
    case EdgePosition::H_nyN_3D:
    case EdgePosition::H_0yn_3D:
    case EdgePosition::H_Nyn_3D:
    case EdgePosition::P_n0z_3D:
    case EdgePosition::P_nNz_3D:
    case EdgePosition::P_0nz_3D:
    case EdgePosition::P_Nnz_3D:
    case EdgePosition::D1_xyn_3D:
    case EdgePosition::D2_nyz_3D:
    case EdgePosition::D3_xnz_3D:
      return 4;
    case EdgePosition::L_x00_3D:
    case EdgePosition::L_xNN_3D:
    case EdgePosition::H_0yN_3D:
    case EdgePosition::H_Ny0_3D:
    case EdgePosition::P_0Nz_3D:
    case EdgePosition::P_N0z_3D:
    case EdgePosition::D1_xy0_3D:
    case EdgePosition::D1_xyN_3D:
    case EdgePosition::D2_0yz_3D:
    case EdgePosition::D2_Nyz_3D:
    case EdgePosition::D3_x0z_3D:
    case EdgePosition::D3_xNz_3D:
      return 3;
    case EdgePosition::L_xN0_3D:
    case EdgePosition::L_x0N_3D:
    case EdgePosition::H_0y0_3D:
    case EdgePosition::H_NyN_3D:
    case EdgePosition::P_00z_3D:
    case EdgePosition::P_NNz_3D:
    case EdgePosition::L_xn_2D:
    case EdgePosition::H_ny_2D:
    case EdgePosition::D1_2D:
      return 2;
    case EdgePosition::L_x0_2D:
    case EdgePosition::L_xN_2D:
    case EdgePosition::H_0y_2D:
    case EdgePosition::H_Ny_2D:
      return 1;

    default: // 1D
      break;
  }

  return 0;
}

int ImplicitTriangulation::getEdgeTriangleInternal(
  const SimplexId &edgeId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= getEdgeTriangleNumberInternal(edgeId))
    return -1;
#endif

  const auto &p = edgeCoords_[edgeId];

  switch(edgePositions_[edgeId]) {
    case EdgePosition::L_xnn_3D:
      triangleId = getEdgeTriangleL_xnn(p.data(), localTriangleId);
      break;
    case EdgePosition::L_xn0_3D:
      triangleId = getEdgeTriangleL_xn0(p.data(), localTriangleId);
      break;
    case EdgePosition::L_xnN_3D:
      triangleId = getEdgeTriangleL_xnN(p.data(), localTriangleId);
      break;
    case EdgePosition::L_x0n_3D:
      triangleId = getEdgeTriangleL_x0n(p.data(), localTriangleId);
      break;
    case EdgePosition::L_x00_3D:
      triangleId = getEdgeTriangleL_x00(p.data(), localTriangleId);
      break;
    case EdgePosition::L_x0N_3D:
      triangleId = getEdgeTriangleL_x0N(p.data(), localTriangleId);
      break;
    case EdgePosition::L_xNn_3D:
      triangleId = getEdgeTriangleL_xNn(p.data(), localTriangleId);
      break;
    case EdgePosition::L_xN0_3D:
      triangleId = getEdgeTriangleL_xN0(p.data(), localTriangleId);
      break;
    case EdgePosition::L_xNN_3D:
      triangleId = getEdgeTriangleL_xNN(p.data(), localTriangleId);
      break;
    case EdgePosition::H_nyn_3D:
      triangleId = getEdgeTriangleH_nyn(p.data(), localTriangleId);
      break;
    case EdgePosition::H_ny0_3D:
      triangleId = getEdgeTriangleH_ny0(p.data(), localTriangleId);
      break;
    case EdgePosition::H_nyN_3D:
      triangleId = getEdgeTriangleH_nyN(p.data(), localTriangleId);
      break;
    case EdgePosition::H_0yn_3D:
      triangleId = getEdgeTriangleH_0yn(p.data(), localTriangleId);
      break;
    case EdgePosition::H_0y0_3D:
      triangleId = getEdgeTriangleH_0y0(p.data(), localTriangleId);
      break;
    case EdgePosition::H_0yN_3D:
      triangleId = getEdgeTriangleH_0yN(p.data(), localTriangleId);
      break;
    case EdgePosition::H_Nyn_3D:
      triangleId = getEdgeTriangleH_Nyn(p.data(), localTriangleId);
      break;
    case EdgePosition::H_Ny0_3D:
      triangleId = getEdgeTriangleH_Ny0(p.data(), localTriangleId);
      break;
    case EdgePosition::H_NyN_3D:
      triangleId = getEdgeTriangleH_NyN(p.data(), localTriangleId);
      break;
    case EdgePosition::P_nnz_3D:
      triangleId = getEdgeTriangleP_nnz(p.data(), localTriangleId);
      break;
    case EdgePosition::P_n0z_3D:
      triangleId = getEdgeTriangleP_n0z(p.data(), localTriangleId);
      break;
    case EdgePosition::P_nNz_3D:
      triangleId = getEdgeTriangleP_nNz(p.data(), localTriangleId);
      break;
    case EdgePosition::P_0nz_3D:
      triangleId = getEdgeTriangleP_0nz(p.data(), localTriangleId);
      break;
    case EdgePosition::P_00z_3D:
      triangleId = getEdgeTriangleP_00z(p.data(), localTriangleId);
      break;
    case EdgePosition::P_0Nz_3D:
      triangleId = getEdgeTriangleP_0Nz(p.data(), localTriangleId);
      break;
    case EdgePosition::P_Nnz_3D:
      triangleId = getEdgeTriangleP_Nnz(p.data(), localTriangleId);
      break;
    case EdgePosition::P_N0z_3D:
      triangleId = getEdgeTriangleP_N0z(p.data(), localTriangleId);
      break;
    case EdgePosition::P_NNz_3D:
      triangleId = getEdgeTriangleP_NNz(p.data(), localTriangleId);
      break;
    case EdgePosition::D1_xyn_3D:
      triangleId = getEdgeTriangleD1_xyn(p.data(), localTriangleId);
      break;
    case EdgePosition::D1_xy0_3D:
      triangleId = getEdgeTriangleD1_xy0(p.data(), localTriangleId);
      break;
    case EdgePosition::D1_xyN_3D:
      triangleId = getEdgeTriangleD1_xyN(p.data(), localTriangleId);
      break;
    case EdgePosition::D2_nyz_3D:
      triangleId = getEdgeTriangleD2_nyz(p.data(), localTriangleId);
      break;
    case EdgePosition::D2_0yz_3D:
      triangleId = getEdgeTriangleD2_0yz(p.data(), localTriangleId);
      break;
    case EdgePosition::D2_Nyz_3D:
      triangleId = getEdgeTriangleD2_Nyz(p.data(), localTriangleId);
      break;
    case EdgePosition::D3_xnz_3D:
      triangleId = getEdgeTriangleD3_xnz(p.data(), localTriangleId);
      break;
    case EdgePosition::D3_x0z_3D:
      triangleId = getEdgeTriangleD3_x0z(p.data(), localTriangleId);
      break;
    case EdgePosition::D3_xNz_3D:
      triangleId = getEdgeTriangleD3_xNz(p.data(), localTriangleId);
      break;
    case EdgePosition::D4_3D:
      triangleId = getEdgeTriangleD4_xyz(p.data(), localTriangleId);
      break;

    case EdgePosition::L_xn_2D:
      triangleId = getEdgeTriangleL_xn(p.data(), localTriangleId);
      break;
    case EdgePosition::L_x0_2D:
      triangleId = getEdgeTriangleL_x0(p.data(), localTriangleId);
      break;
    case EdgePosition::L_xN_2D:
      triangleId = getEdgeTriangleL_xN(p.data(), localTriangleId);
      break;
    case EdgePosition::H_ny_2D:
      triangleId = getEdgeTriangleH_ny(p.data(), localTriangleId);
      break;
    case EdgePosition::H_0y_2D:
      triangleId = getEdgeTriangleH_0y(p.data(), localTriangleId);
      break;
    case EdgePosition::H_Ny_2D:
      triangleId = getEdgeTriangleH_Ny(p.data(), localTriangleId);
      break;
    case EdgePosition::D1_2D:
      triangleId = getEdgeTriangleD1_xy(p.data(), localTriangleId);
      break;

    default: // 1D
      triangleId = -1;
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::getEdgeTrianglesInternal() {
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

inline SimplexId ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getEdgeLinkNumber)(const SimplexId &edgeId) const {
  return getEdgeStarNumber(edgeId);
}

int ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdgeLink)(
  const SimplexId &edgeId, const int &localLinkId, SimplexId &linkId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getEdgeLinkNumber(edgeId))
    return -1;
#endif

  const auto &p = edgeCoords_[edgeId];

  switch(edgePositions_[edgeId]) {
  CASE_EDGE_POSITION_L_3D:
    linkId = getEdgeLinkL(p.data(), localLinkId);
    break;
  CASE_EDGE_POSITION_H_3D:
    linkId = getEdgeLinkH(p.data(), localLinkId);
    break;
  CASE_EDGE_POSITION_P_3D:
    linkId = getEdgeLinkP(p.data(), localLinkId);
    break;
  CASE_EDGE_POSITION_D1_3D:
    linkId = getEdgeLinkD1(p.data(), localLinkId);
    break;
  CASE_EDGE_POSITION_D2_3D:
    linkId = getEdgeLinkD2(p.data(), localLinkId);
    break;
  CASE_EDGE_POSITION_D3_3D:
    linkId = getEdgeLinkD3(p.data(), localLinkId);
    break;
    case EdgePosition::D4_3D:
      linkId = getEdgeLinkD4(p.data(), localLinkId);
      break;

    CASE_EDGE_POSITION_L_2D:
      linkId = getEdgeLink2dL(p.data(), localLinkId);
      break;
    CASE_EDGE_POSITION_H_2D:
      linkId = getEdgeLink2dH(p.data(), localLinkId);
      break;
    case EdgePosition::D1_2D:
      linkId = getEdgeLink2dD1(p.data(), localLinkId);
      break;

    default: // 1D
      linkId = -1;
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() {

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

inline SimplexId ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getEdgeStarNumber)(const SimplexId &edgeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= edgeNumber_)
    return -1;
#endif

  switch(edgePositions_[edgeId]) {
    case EdgePosition::L_xnn_3D:
    case EdgePosition::H_nyn_3D:
    case EdgePosition::P_nnz_3D:
    case EdgePosition::D4_3D:
      return 6;
    case EdgePosition::D1_xyn_3D:
    case EdgePosition::D2_nyz_3D:
    case EdgePosition::D3_xnz_3D:
      return 4;
    case EdgePosition::L_x0n_3D:
    case EdgePosition::L_xNn_3D:
    case EdgePosition::L_xn0_3D:
    case EdgePosition::L_xnN_3D:
    case EdgePosition::H_ny0_3D:
    case EdgePosition::H_nyN_3D:
    case EdgePosition::H_0yn_3D:
    case EdgePosition::H_Nyn_3D:
    case EdgePosition::P_n0z_3D:
    case EdgePosition::P_nNz_3D:
    case EdgePosition::P_0nz_3D:
    case EdgePosition::P_Nnz_3D:
      return 3;
    case EdgePosition::L_x00_3D:
    case EdgePosition::L_xNN_3D:
    case EdgePosition::H_0yN_3D:
    case EdgePosition::H_Ny0_3D:
    case EdgePosition::P_0Nz_3D:
    case EdgePosition::P_N0z_3D:
    case EdgePosition::D1_xy0_3D:
    case EdgePosition::D1_xyN_3D:
    case EdgePosition::D2_0yz_3D:
    case EdgePosition::D2_Nyz_3D:
    case EdgePosition::D3_x0z_3D:
    case EdgePosition::D3_xNz_3D:
    case EdgePosition::L_xn_2D:
    case EdgePosition::H_ny_2D:
    case EdgePosition::D1_2D:
      return 2;
    case EdgePosition::L_xN0_3D:
    case EdgePosition::L_x0N_3D:
    case EdgePosition::H_0y0_3D:
    case EdgePosition::H_NyN_3D:
    case EdgePosition::P_00z_3D:
    case EdgePosition::P_NNz_3D:
    case EdgePosition::L_x0_2D:
    case EdgePosition::L_xN_2D:
    case EdgePosition::H_0y_2D:
    case EdgePosition::H_Ny_2D:
      return 1;

    default: // 1D
      break;
  }

  return 0;
}

int ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdgeStar)(
  const SimplexId &edgeId, const int &localStarId, SimplexId &starId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getEdgeStarNumber(edgeId))
    return -1;
#endif

  const auto &p = edgeCoords_[edgeId];

  switch(edgePositions_[edgeId]) {
  CASE_EDGE_POSITION_L_3D:
    starId = getEdgeStarL(p.data(), localStarId);
    break;
  CASE_EDGE_POSITION_H_3D:
    starId = getEdgeStarH(p.data(), localStarId);
    break;
  CASE_EDGE_POSITION_P_3D:
    starId = getEdgeStarP(p.data(), localStarId);
    break;
  CASE_EDGE_POSITION_D1_3D:
    starId = getEdgeStarD1(p.data(), localStarId);
    break;
  CASE_EDGE_POSITION_D2_3D:
    starId = getEdgeStarD2(p.data(), localStarId);
    break;
  CASE_EDGE_POSITION_D3_3D:
    starId = getEdgeStarD3(p.data(), localStarId);
    break;
    case EdgePosition::D4_3D:
      starId
        = p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + localStarId;
      break;

    CASE_EDGE_POSITION_L_2D:
      starId = getEdgeStar2dL(p.data(), localStarId);
      break;
    CASE_EDGE_POSITION_H_2D:
      starId = getEdgeStar2dH(p.data(), localStarId);
      break;
    case EdgePosition::D1_2D:
      starId = p[0] * 2 + p[1] * tshift_[0] + localStarId;
      break;

    default: // 1D
      starId = -1;
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getEdgeStars)() {

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

int ImplicitTriangulation::getTriangleVertexInternal(
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

  const auto &p = triangleCoords_[triangleId];
  vertexId = -1;

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
      switch(localVertexId) {
        break;
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
      break;
    case TrianglePosition::BOTTOM_2D:
      switch(localVertexId) {
        break;
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

  return 0;
}

int ImplicitTriangulation::getTriangleEdgeInternal(const SimplexId &triangleId,
                                                   const int &localEdgeId,
                                                   SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
  if(localEdgeId < 0 or localEdgeId >= 3)
    return -2;
#endif

  const auto &p = triangleCoords_[triangleId];
  const auto par = triangleId % 2;
  edgeId = -1;

  switch(trianglePositions_[triangleId]) {
    case TrianglePosition::F_3D:
      edgeId = (par == 1) ? getTriangleEdgeF_1(p.data(), localEdgeId)
                          : getTriangleEdgeF_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::H_3D:
      edgeId = (par == 1) ? getTriangleEdgeH_1(p.data(), localEdgeId)
                          : getTriangleEdgeH_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::C_3D:
      edgeId = (par == 1) ? getTriangleEdgeC_1(p.data(), localEdgeId)
                          : getTriangleEdgeC_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::D1_3D:
      edgeId = (par == 1) ? getTriangleEdgeD1_1(p.data(), localEdgeId)
                          : getTriangleEdgeD1_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::D2_3D:
      edgeId = (par == 1) ? getTriangleEdgeD2_1(p.data(), localEdgeId)
                          : getTriangleEdgeD2_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::D3_3D:
      edgeId = (par == 1) ? getTriangleEdgeD3_1(p.data(), localEdgeId)
                          : getTriangleEdgeD3_0(p.data(), localEdgeId);
      break;
    case TrianglePosition::TOP_2D:
      switch(localEdgeId) {
        break;
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
      break;
    case TrianglePosition::BOTTOM_2D:
      switch(localEdgeId) {
        break;
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

  return 0;
}

int ImplicitTriangulation::getTriangleEdgesInternal(
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
  ImplicitTriangulation::getTriangleEdgesInternal() {
  if(!triangleEdgeList_.size()) {
    Timer t;

    getTriangleEdgesInternal(triangleEdgeList_);

    printMsg("Built " + to_string(triangleNumber_) + " triangle edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &triangleEdgeList_;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangles)() {

  if(!triangleList_.size()) {
    Timer t;

    triangleList_.resize(triangleNumber_);
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      triangleList_[i].resize(3);
      for(int j = 0; j < 3; ++j)
        getTriangleVertexInternal(i, j, triangleList_[i][j]);
    }

    printMsg("Built " + to_string(triangleNumber_) + " triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &triangleList_;
}

int ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
  const SimplexId &triangleId,
  const int &localLinkId,
  SimplexId &linkId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getTriangleLinkNumber(triangleId))
    return -1;
#endif

  const auto p = triangleCoords_[triangleId];

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
    default: // 2D
      linkId = -1;
      break;
  }

  return 0;
}

inline SimplexId ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getTriangleLinkNumber)(const SimplexId &triangleId) const {
  return getTriangleStarNumber(triangleId);
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() {
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

inline SimplexId ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getTriangleStarNumber)(const SimplexId &triangleId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  const auto p = triangleCoords_[triangleId];

  switch(trianglePositions_[triangleId]) {
    case TrianglePosition::F_3D:
      return (p[2] > 0 and p[2] < nbvoxels_[2]) ? 2 : 1;
    case TrianglePosition::H_3D:
      return (p[1] > 0 and p[1] < nbvoxels_[1]) ? 2 : 1;
    case TrianglePosition::C_3D:
      return (p[0] < 2 or p[0] >= (dimensions_[0] * 2 - 2)) ? 1 : 2;

    case TrianglePosition::D1_3D:
    case TrianglePosition::D2_3D:
    case TrianglePosition::D3_3D:
      return 2;
    default: // 2D
      break;
  }
  return 0;
}

int ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
  const SimplexId &triangleId,
  const int &localStarId,
  SimplexId &starId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getTriangleStarNumber(triangleId))
    return -1;
#endif

  const auto p = triangleCoords_[triangleId];

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
    default: // 2D
      starId = -1;
      break;
  }

  return 0;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getTriangleStars)() {

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

inline SimplexId ImplicitTriangulation::getTriangleNeighborNumber(
  const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  if(dimensionality_ == 2) {
    const auto p = triangleCoords_[triangleId];
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
    const auto p = triangleCoords_[triangleId].data();
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
    const SimplexId id = tetId % 6;
    const auto p = tetrahedronCoords_[tetId].data();

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
    const SimplexId id = tetId % 6;
    const auto p = tetrahedronCoords_[tetId].data();

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
    const SimplexId id = tetId % 6;
    const auto p = tetrahedronCoords_[tetId].data();

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
    const auto p = tetrahedronCoords_[tetId].data();

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
    const auto p = tetrahedronCoords_[tetId].data();

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

SimplexId ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
  getCellVertexNumber)(const SimplexId &cellId) const {
  return dimensionality_ + 1;
}

int ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getCellVertex)(
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

SimplexId ImplicitTriangulation::getCellEdgeNumberInternal(
  const SimplexId &cellId) const {
  if(dimensionality_ == 3)
    return 6;
  else if(dimensionality_ == 2)
    return 3;

  return 0;
}

int ImplicitTriangulation::getCellEdgeInternal(const SimplexId &cellId,
                                               const int &localEdgeId,
                                               SimplexId &edgeId) const {
  if(dimensionality_ == 3)
    getTetrahedronEdge(cellId, localEdgeId, edgeId);
  else if(dimensionality_ == 2)
    getTriangleEdgeInternal(cellId, localEdgeId, edgeId);
  else if(dimensionality_ == 1)
    getCellNeighbor(cellId, localEdgeId, edgeId);

  return 0;
}

const vector<vector<SimplexId>> *ImplicitTriangulation::getCellEdgesInternal() {
  if(!cellEdgeList_.size()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronEdges(cellEdgeList_);
    else if(dimensionality_ == 2)
      getTriangleEdgesInternal(cellEdgeList_);

    printMsg("Built " + to_string(cellNumber_) + " cell edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &cellEdgeList_;
}

int ImplicitTriangulation::getCellTriangleInternal(
  const SimplexId &cellId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
  if(dimensionality_ == 3)
    getTetrahedronTriangle(cellId, localTriangleId, triangleId);

  return 0;
}

const vector<vector<SimplexId>> *
  ImplicitTriangulation::getCellTrianglesInternal() {
  if(!cellTriangleList_.size()) {
    Timer t;

    if(dimensionality_ == 3)
      getTetrahedronTriangles(cellTriangleList_);

    printMsg("Built " + to_string(cellNumber_) + " cell triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &cellTriangleList_;
}

SimplexId ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(
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

int ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
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
  ImplicitTriangulation::TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() {
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
