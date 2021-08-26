#include <MultiresTriangulation.h>

using ttk::MultiresTriangulation;
using ttk::SimplexId;

SimplexId MultiresTriangulation::getVertexNeighborAtDecimation(
  const SimplexId &vertexId,
  const int &localNeighborId,
  SimplexId &neighborId,
  int decimation) const {

  const auto &p = this->vertexCoords_[vertexId];

  SimplexId shiftX = decimation;
  SimplexId shiftY = decimation;
  SimplexId shiftZ = decimation;

  if(dimensionality_ == 2) {
    if((nbvoxels_[Di_] % decimation) and (p[0] + decimation > nbvoxels_[Di_])) {
      shiftX = nbvoxels_[0] % decimation;
    }
    if((nbvoxels_[Dj_] % decimation) and (p[1] + decimation > nbvoxels_[Dj_])) {
      shiftY = nbvoxels_[1] % decimation;
    }
  } else if(dimensionality_ == 3) {
    if((nbvoxels_[0] % decimation) and (p[0] + decimation > nbvoxels_[0])) {
      shiftX = nbvoxels_[0] % decimation;
    }
    if((nbvoxels_[1] % decimation) and (p[1] + decimation > nbvoxels_[1])) {
      shiftY = nbvoxels_[1] % decimation;
    }
    if((nbvoxels_[2] % decimation) and (p[2] + decimation > nbvoxels_[2])) {
      shiftZ = nbvoxels_[2] % decimation;
    }
  }

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      neighborId = getVertexNeighborAtDecimationABCDEFGH(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::FRONT_FACE_3D:
      neighborId = getVertexNeighborAtDecimationABDC(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::BACK_FACE_3D:
      neighborId = getVertexNeighborAtDecimationEFHG(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::TOP_FACE_3D:
      neighborId = getVertexNeighborAtDecimationAEFB(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::BOTTOM_FACE_3D:
      neighborId = getVertexNeighborAtDecimationGHDC(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::LEFT_FACE_3D:
      neighborId = getVertexNeighborAtDecimationAEGC(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::RIGHT_FACE_3D:
      neighborId = getVertexNeighborAtDecimationBFHD(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::TOP_FRONT_EDGE_3D: // ab
      neighborId = getVertexNeighborAtDecimationAB(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      neighborId = getVertexNeighborAtDecimationCD(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      neighborId = getVertexNeighborAtDecimationAC(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      neighborId = getVertexNeighborAtDecimationBD(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::TOP_BACK_EDGE_3D: // ef
      neighborId = getVertexNeighborAtDecimationEF(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      neighborId = getVertexNeighborAtDecimationGH(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
      neighborId = getVertexNeighborAtDecimationEG(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      neighborId = getVertexNeighborAtDecimationFH(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
      neighborId = getVertexNeighborAtDecimationAE(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      neighborId = getVertexNeighborAtDecimationBF(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      neighborId = getVertexNeighborAtDecimationCG(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      neighborId = getVertexNeighborAtDecimationDH(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      neighborId = getVertexNeighborAtDecimationA(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      neighborId = getVertexNeighborAtDecimationB(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      neighborId = getVertexNeighborAtDecimationC(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      neighborId = getVertexNeighborAtDecimationD(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      neighborId = getVertexNeighborAtDecimationE(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      neighborId = getVertexNeighborAtDecimationF(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      neighborId = getVertexNeighborAtDecimationG(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      neighborId = getVertexNeighborAtDecimationH(
        vertexId, localNeighborId, shiftX, shiftY, shiftZ, decimation);
      break;
    case VertexPosition::CENTER_2D:
      neighborId = getVertexNeighborAtDecimation2dABCD(
        vertexId, localNeighborId, shiftX, shiftY, decimation);
      break;
    case VertexPosition::TOP_EDGE_2D:
      neighborId = getVertexNeighborAtDecimation2dAB(
        vertexId, localNeighborId, shiftX, shiftY, decimation);
      break;
    case VertexPosition::BOTTOM_EDGE_2D:
      neighborId = getVertexNeighborAtDecimation2dCD(
        vertexId, localNeighborId, shiftX, shiftY, decimation);
      break;
    case VertexPosition::LEFT_EDGE_2D:
      neighborId = getVertexNeighborAtDecimation2dAC(
        vertexId, localNeighborId, shiftX, shiftY, decimation);
      break;
    case VertexPosition::RIGHT_EDGE_2D:
      neighborId = getVertexNeighborAtDecimation2dBD(
        vertexId, localNeighborId, shiftX, shiftY, decimation);
      break;
    case VertexPosition::TOP_LEFT_CORNER_2D: // a
      neighborId = getVertexNeighborAtDecimation2dA(
        vertexId, localNeighborId, shiftX, shiftY, decimation);
      break;
    case VertexPosition::TOP_RIGHT_CORNER_2D: // b
      neighborId = getVertexNeighborAtDecimation2dB(
        vertexId, localNeighborId, shiftX, shiftY, decimation);
      break;
    case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      neighborId = getVertexNeighborAtDecimation2dC(
        vertexId, localNeighborId, shiftX, shiftY, decimation);
      break;
    case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      neighborId = getVertexNeighborAtDecimation2dD(
        vertexId, localNeighborId, shiftX, shiftY, decimation);
      break;
    case VertexPosition::CENTER_1D:
      neighborId = (localNeighborId == 0 ? vertexId + decimation
                                         : vertexId - decimation);
      break;
    case VertexPosition::LEFT_CORNER_1D:
      neighborId = vertexId + decimation;
      break;
    case VertexPosition::RIGHT_CORNER_1D:
      neighborId = vertexId - decimation;
      break;
    default:
      neighborId = -1;
      break;
  }

  return 0;
}

SimplexId MultiresTriangulation::getVertexNeighborAtDecimation2dA(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const int ttkNotUsed(decimation)) const {
  // V(a)={b,c}
  switch(id) {
    case 0:
      return v + shiftX; // b
    case 1:
      return v + gridDimensions_[Di_] * shiftY; // c
  }
  return -1;
}
SimplexId MultiresTriangulation::getVertexNeighborAtDecimation2dB(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const int ttkNotUsed(decimation)) const {
  // V(b)={a,c,d}
  switch(id) {
    case 0:
      return v - shiftX; // a
    case 1:
      return v + gridDimensions_[Di_] * shiftY; // d
    case 2:
      return v + (gridDimensions_[Di_] * shiftY - shiftX); // c
  }
  return -1;
}

SimplexId MultiresTriangulation::getVertexNeighborAtDecimation2dC(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const int ttkNotUsed(decimation)) const {
  // V(c)={a,b,d}
  switch(id) {
    case 0:
      return v + shiftX;
    case 1:
      return v - gridDimensions_[Di_] * shiftY; // a
    case 2:
      return v + (shiftX - gridDimensions_[Di_] * shiftY); // b
  }
  return -1;
}
SimplexId MultiresTriangulation::getVertexNeighborAtDecimation2dD(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const int ttkNotUsed(decimation)) const {
  // V(d)={c,b}
  switch(id) {
    case 0:
      return v - shiftX; // c
    case 1:
      return v - gridDimensions_[Di_] * shiftY; // b
  }
  return -1;
}
SimplexId MultiresTriangulation::getVertexNeighborAtDecimation2dAB(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const int decimation) const {
  // V(ab)=V(b)::{a,c,d}+V(a)::{b}
  switch(id) {
    case 0:
      return v - decimation; // V(b)::a
    case 1:
      return v + (gridDimensions_[Di_] * shiftY - decimation); // V(b)::c
    case 2:
      return v + gridDimensions_[Di_] * shiftY; // V(b)::d
    case 3:
      return v + shiftX; // V(a)::b
  }
  return -1;
}
SimplexId MultiresTriangulation::getVertexNeighborAtDecimation2dCD(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const int decimation) const {
  // V(cd)=V(c)::{a,b,d}+V(d)::{c}
  switch(id) {
    case 0:
      return v - decimation; // V(d)::c
    case 1:
      return v - gridDimensions_[Di_] * shiftY; // V(c)::a
    case 2:
      return v + (shiftX - gridDimensions_[Di_] * shiftY); // V(c)::b
    case 3:
      return v + shiftX; // V(c)::d
  }
  return -1;
}
SimplexId MultiresTriangulation::getVertexNeighborAtDecimation2dAC(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const int decimation) const {
  // V(ac)=V(c)::{a,b,d}+V(a)::{c}
  switch(id) {
    case 0:
      return v - gridDimensions_[Di_] * decimation; // V(c)::{a}
    case 1:
      return v + (shiftX - gridDimensions_[Di_] * decimation); // V(c)::{b}
    case 2:
      return v + shiftX; // V(c)::{d}
    case 3:
      return v + gridDimensions_[Di_] * shiftY; // V(a)::{c}
  }
  return -1;
}
SimplexId MultiresTriangulation::getVertexNeighborAtDecimation2dBD(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const int decimation) const {
  // V(bd)=V(b)::{c,d}+V(d)::{b,c}
  switch(id) {
    case 0:
      return v + (gridDimensions_[Di_] * shiftY - shiftX); // V(b)::{c}
    case 1:
      return v + gridDimensions_[Di_] * shiftY; // V(b)::{d}
    case 2:
      return v - gridDimensions_[Di_] * decimation; // V(d)::{b}
    case 3:
      return v - shiftX; // V(d)::{c}
  }
  return -1;
}

SimplexId MultiresTriangulation::getVertexNeighborAtDecimation2dABCD(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const int decimation) const {
  // V(abcd)=V(d)::{b,c}+V(c)::{b,d}+V(a)::{b}+V(b)::{c}
  switch(id) {
    case 0:
      return v - decimation;
    case 1:
      return v - gridDimensions_[Di_] * decimation;
    case 2:
      return v + (shiftX - gridDimensions_[Di_] * decimation);
    case 3:
      return v + shiftX;
    case 4:
      return v + gridDimensions_[Di_] * shiftY;
    case 5:
      return v + (gridDimensions_[Di_] * shiftY - decimation);
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationA(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int ttkNotUsed(decimation)) const {
  // V(a)={b,c,e,g}
  switch(id) {
    case 0:
      return v + shiftX; // b
    case 1:
      return v + vshift_[0] * shiftY; // c
    case 2:
      return v + vshift_[1] * shiftZ; // e
    case 3:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // g
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationB(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int ttkNotUsed(decimation)) const {
  // V(b)={a,c,d,e,f,g,h}
  switch(id) {
    case 0:
      return v - shiftX; // a
    case 1:
      return v + (vshift_[0] * shiftY - shiftX); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - shiftX); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
  }
  return -1;
}
inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationC(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int ttkNotUsed(decimation)) const {
  // V(c)={a,b,d,g}
  switch(id) {
    case 0:
      return v - vshift_[0] * shiftY; // a
    case 1:
      return v + (shiftX - vshift_[0] * shiftY); // b
    case 2:
      return v + shiftX; // d
    case 3:
      return v + vshift_[1] * shiftZ; // g
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationD(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int ttkNotUsed(decimation)) const {
  // V(d)={b,c,g,h}
  switch(id) {
    case 0:
      return v - vshift_[0] * shiftY; // b
    case 1:
      return v - shiftX; // c
    case 2:
      return v + (-shiftX + vshift_[1] * shiftZ); // g
    case 3:
      return v + vshift_[1] * shiftZ; // h
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationE(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int ttkNotUsed(decimation)) const {
  // V(e)={a,b,f,g}
  switch(id) {
    case 0:
      return v - vshift_[1] * shiftZ; // a
    case 1:
      return v + (shiftX - vshift_[1] * shiftZ); // b
    case 2:
      return v + shiftX; // f
    case 3:
      return v + vshift_[0] * shiftY; // g
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationF(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int ttkNotUsed(decimation)) const {
  // V(f)={b,e,g,h}
  switch(id) {
    case 0:
      return v - vshift_[1] * shiftZ; // b
    case 1:
      return v - shiftX; // e
    case 2:
      return v + (vshift_[0] * shiftY - shiftX); // g
    case 3:
      return v + vshift_[0] * shiftY; // h
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationG(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int ttkNotUsed(decimation)) const {
  // V(g)={a,b,c,d,e,f,h}
  switch(id) {
    case 0:
      return v - (vshift_[0] * shiftY + vshift_[1] * shiftZ); // a
    case 1:
      return v + (shiftX - vshift_[0] * shiftY - vshift_[1] * shiftZ); // b
    case 2:
      return v - vshift_[1] * shiftZ; // c
    case 3:
      return v + (shiftX - vshift_[1] * shiftZ); // d
    case 4:
      return v - vshift_[0] * shiftY; // e
    case 5:
      return v + (shiftX - vshift_[0] * shiftY); // f
    case 6:
      return v + shiftX; // h
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationH(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int ttkNotUsed(decimation)) const {
  // V(h)={b,d,f,g}
  switch(id) {
    case 0:
      return v - (vshift_[0] * shiftY + vshift_[1] * shiftZ); // b
    case 1:
      return v - vshift_[1] * shiftZ; // d
    case 2:
      return v - vshift_[0] * shiftY; // f
    case 3:
      return v - shiftX; // g
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationAB(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(ab)=V(b)+V(a)::{b}
  switch(id) {
    case 0:
      return v - decimation; // a
    case 1:
      return v + (vshift_[0] * shiftY - decimation); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - decimation); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - decimation); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v + shiftX; // V(a)::{b}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationCD(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(cd)=V(d)+V(c)::{b,d}
  switch(id) {
    case 0:
      return v - vshift_[0] * shiftY; // b
    case 1:
      return v - decimation; // c
    case 2:
      return v + (vshift_[1] * shiftZ - decimation); // g
    case 3:
      return v + vshift_[1] * shiftZ; // h
    case 4:
      return v + (shiftX - vshift_[0] * shiftY); // V(c)::{b}
    case 5:
      return v + shiftX; // V(c)::{d}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationEF(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(ef)=V(f)+V(e)::{b,f}
  switch(id) {
    case 0:
      return v - vshift_[1] * shiftZ; // b
    case 1:
      return v - decimation; // e
    case 2:
      return v + (vshift_[0] * shiftY - decimation); // g
    case 3:
      return v + vshift_[0] * shiftY; // h
    case 4:
      return v + (shiftX - vshift_[1] * shiftZ); // V(e)::{b}
    case 5:
      return v + shiftX; // V(e)::{f}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationGH(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(gh)=V(g)+V(h)::{g}
  switch(id) {
    case 0:
      return v - (vshift_[0] * shiftY + vshift_[1] * shiftZ); // a
    case 1:
      return v + (shiftX - vshift_[0] * shiftY - vshift_[1] * shiftZ); // b
    case 2:
      return v - vshift_[1] * shiftZ; // c
    case 3:
      return v + (shiftX - vshift_[1] * shiftZ); // d
    case 4:
      return v - vshift_[0] * shiftY; // e
    case 5:
      return v + (shiftX - vshift_[0] * shiftY); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v - decimation; // V(h)::{g}
  }

  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationAC(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(ac)=V(c)+V(a)::{c,g}
  switch(id) {
    case 0:
      return v - vshift_[0] * decimation; // a
    case 1:
      return v + (shiftX - vshift_[0] * decimation); // b
    case 2:
      return v + shiftX; // d
    case 3:
      return v + vshift_[1] * shiftZ; // g
    case 4:
      return v + vshift_[0] * shiftY; // V(a)::{c}
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // V(a)::{c}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationBD(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(bd)=V(b)+V(d)::{b}
  switch(id) {
    case 0:
      return v - shiftX; // a
    case 1:
      return v + (vshift_[0] * shiftY - shiftX); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - shiftX); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v - vshift_[0] * decimation; // V(d)::{b}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationEG(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(eg)=V(g)+V(e)::{g}
  switch(id) {
    case 0:
      return v - (vshift_[0] * decimation + vshift_[1] * shiftZ); // a
    case 1:
      return v + (shiftX - vshift_[0] * decimation - vshift_[1] * shiftZ); // b
    case 2:
      return v - vshift_[1] * shiftZ; // c
    case 3:
      return v + (shiftX - vshift_[1] * shiftZ); // d
    case 4:
      return v - vshift_[0] * decimation; // e
    case 5:
      return v + (shiftX - vshift_[0] * decimation); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v + vshift_[0] * shiftY; // V(e)::{g}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationFH(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(fh)=V(f)+V(h)::{b,f}
  switch(id) {
    case 0:
      return v - vshift_[1] * shiftZ; // b
    case 1:
      return v - shiftX; // e
    case 2:
      return v + (vshift_[0] * shiftY - shiftX); // g
    case 3:
      return v + vshift_[0] * shiftY; // h
    case 4:
      return v - (vshift_[0] * decimation + vshift_[1] * shiftZ); // V(h)::{b}
    case 5:
      return v - vshift_[0] * decimation; // V(h)::{f}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationAE(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(ae)=V(a)+V(e)::{a,b}
  switch(id) {
    case 0:
      return v + shiftX; // b
    case 1:
      return v + vshift_[0] * shiftY; // c
    case 2:
      return v + vshift_[1] * shiftZ; // e
    case 3:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // g
    case 4:
      return v - vshift_[1] * decimation; // V(e)::{a}
    case 5:
      return v + (shiftX - vshift_[1] * decimation); // V(e)::{b}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationBF(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(bf)=V(b)+V(f)::{b}
  switch(id) {
    case 0:
      return v - shiftX; // a
    case 1:
      return v + (vshift_[0] * shiftY - shiftX); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - shiftX); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v - vshift_[1] * decimation; // V(f)::{b}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationCG(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(cg)=V(g)+V(c)::{g}
  switch(id) {
    case 0:
      return v - (vshift_[0] * shiftY + vshift_[1] * decimation); // a
    case 1:
      return v + (shiftX - vshift_[0] * shiftY - vshift_[1] * decimation); // b
    case 2:
      return v - vshift_[1] * decimation; // c
    case 3:
      return v + (shiftX - vshift_[1] * decimation); // d
    case 4:
      return v - vshift_[0] * shiftY; // e
    case 5:
      return v + (shiftX - vshift_[0] * shiftY); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v + vshift_[1] * shiftZ; // V(c)::{g}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationDH(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(dh)=V(d)+V(h)::{b,d}
  switch(id) {
    case 0:
      return v - vshift_[0] * shiftY; // b
    case 1:
      return v - shiftX; // c
    case 2:
      return v + (vshift_[1] * shiftZ - shiftX); // g
    case 3:
      return v + vshift_[1] * shiftZ; // h
    case 4:
      return v - (vshift_[0] * shiftY + vshift_[1] * decimation); // V(h)::{b}
    case 5:
      return v - vshift_[1] * decimation; // V(h)::{d}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationABDC(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(abdc)=V(b)+V(d)::{b}+V(c)::{b}+V(a)::{b}
  switch(id) {
    case 0:
      return v - decimation; // a
    case 1:
      return v + (vshift_[0] * shiftY - decimation); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - decimation); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - decimation); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v - vshift_[0] * decimation; // V(d)::{b}
    case 8:
      return v + (shiftX - vshift_[0] * decimation); // V(c)::{b}
    case 9:
      return v + shiftX; // V(a)::{b}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationEFHG(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(efhg)=V(g)+V(h)::{g}+V(f)::{g,h}
  switch(id) {
    case 0:
      return v - (vshift_[0] * decimation + vshift_[1] * shiftZ); // a
    case 1:
      return v + (shiftX - vshift_[0] * decimation - vshift_[1] * shiftZ); // b
    case 2:
      return v - vshift_[1] * shiftZ; // c
    case 3:
      return v + (shiftX - vshift_[1] * shiftZ); // d
    case 4:
      return v - vshift_[0] * decimation; // e
    case 5:
      return v + (shiftX - vshift_[0] * decimation); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v - decimation; // V(h)::{g}
    case 8:
      return v + (vshift_[0] * shiftY - decimation); // V(f)::{g}
    case 9:
      return v + vshift_[0] * shiftY; // V(f)::{h}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationAEGC(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(aegc)=V(g)+V(a)::{c,g}+V(c)::{g}
  switch(id) {
    case 0:
      return v - (vshift_[0] * decimation + vshift_[1] * decimation); // a
    case 1:
      return v
             + (shiftX - vshift_[0] * decimation
                - vshift_[1] * decimation); // b
    case 2:
      return v - vshift_[1] * decimation; // c
    case 3:
      return v + (shiftX - vshift_[1] * decimation); // d
    case 4:
      return v - vshift_[0] * decimation; // e
    case 5:
      return v + (shiftX - vshift_[0] * decimation); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v + vshift_[0] * shiftY; // V(a)::{c}
    case 8:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // V(a)::{g}
    case 9:
      return v + vshift_[1] * shiftZ; // V(c)::{g}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationBFHD(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(bfhd)=V(b)+V(f)::{b}+V(h)::{b}+V(d)::{b}
  switch(id) {
    case 0:
      return v - shiftX; // a
    case 1:
      return v + (vshift_[0] * shiftY - shiftX); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - shiftX); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v - vshift_[1] * decimation; // V(f)::{b}
    case 8:
      return v
             - (vshift_[0] * decimation + vshift_[1] * decimation); // V(h)::{b}
    case 9:
      return v - vshift_[0] * decimation; // V(d)::{b}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationAEFB(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(aefb)=V(b)+V(a)::{b}+V(e)::{b}+V(f)::{b}
  switch(id) {
    case 0:
      return v - decimation; // a
    case 1:
      return v + (vshift_[0] * shiftY - decimation); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - decimation); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - decimation); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v + shiftX; // V(a)::{b}
    case 8:
      return v + (shiftX - vshift_[1] * decimation); // V(e)::{b}
    case 9:
      return v - vshift_[1] * decimation; // V(f)::{b}
  }
  return -1;
}

inline ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationGHDC(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(ghdc)=V(g)+V(h)::{g}+V(d)::{g,h}
  switch(id) {
    case 0:
      return v - (vshift_[0] * shiftY + vshift_[1] * decimation); // a
    case 1:
      return v + (shiftX - vshift_[0] * shiftY - vshift_[1] * decimation); // b
    case 2:
      return v - vshift_[1] * decimation; // c
    case 3:
      return v + (shiftX - vshift_[1] * decimation); // d
    case 4:
      return v - vshift_[0] * shiftY; // e
    case 5:
      return v + (shiftX - vshift_[0] * shiftY); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v - decimation; // V(h)::{g}
    case 8:
      return v + (vshift_[1] * shiftZ - decimation); // V(d)::{g}
    case 9:
      return v + vshift_[1] * shiftZ; // V(d)::{h}
  }
  return -1;
}

ttk::SimplexId MultiresTriangulation::getVertexNeighborAtDecimationABCDEFGH(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  const int decimation) const {
  // V(abcdefgh)=V(g)+V(d)::{g,h}+V(h)::{g}+V(b)::{c,d,g,h}
  switch(id) {
    case 0:
      return v - (vshift_[0] * decimation + vshift_[1] * decimation); // a
    case 1:
      return v
             + (shiftX - vshift_[0] * decimation
                - vshift_[1] * decimation); // b
    case 2:
      return v - vshift_[1] * decimation; // c
    case 3:
      return v + (shiftX - vshift_[1] * decimation); // d
    case 4:
      return v - vshift_[0] * decimation; // e
    case 5:
      return v + (shiftX - vshift_[0] * decimation); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v + (vshift_[1] * shiftZ - decimation); // V(d)::{g}
    case 8:
      return v + vshift_[1] * shiftZ; // V(d)::{h}
    case 9:
      return v - decimation; // V(h)::{g}
    case 10:
      return v + (vshift_[0] * shiftY - decimation); // V(b)::{c}
    case 11:
      return v + vshift_[0] * shiftY; // V(b)::{d}
    case 12:
      return v
             + (vshift_[0] * shiftY + vshift_[1] * shiftZ
                - decimation); // V(b)::{g}
    case 13:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // V(b)::{h}
  }
  return -1;
}
