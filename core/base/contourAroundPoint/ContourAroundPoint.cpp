#include "ContourAroundPoint.hpp"


int ttk::ContourAroundPoint::setupDomain(Triangulation *triangulation,
                                         void *scalars) {
  msg(std::string(60, '-').c_str());

  _inpFieldTriangulation = triangulation;
  _inpFieldScalars = scalars;
  _inpDimMax = triangulation->getDimensionality();
  if(!triangulation)
    return -1;
  if(!scalars)
    return -2;
  if(_inpDimMax < 2 || _inpDimMax > 3)
    return -3;

  // Call all the required precondition functions here!
  
  //for getVertexEdgeNumber, getVertexEdge
  triangulation->preconditionVertexEdges();
  // for getVertexNeighbor
  triangulation->preconditionVertexNeighbors();
  // for getEdgeVertex
  triangulation->preconditionEdges();
  // for getEdgeTriangleNumber, getEdgeTriangle
  triangulation->preconditionEdgeTriangles();
  // for getTriangleEdge
  triangulation->preconditionTriangleEdges();
  
  return 0;
}

//----------------------------------------------------------------------------//

int ttk::ContourAroundPoint::setupConstraints(float *coords, float *isovalues,
                                              int *flags, std::size_t np) {
  _inpPointCoords = coords;
  _inpPointIsovals = isovalues;
  _inpPointFlags = flags;
  _np = np;
  if(!coords)
    return -1;
  if(!isovalues)
    return -2;
  if(!flags)
    return -3;
  return 0;
}

//----------------------------------------------------------------------------//

ttk::SimplexId ttk::ContourAroundPoint::findInpVert(SimplexId p) const {
  // This implementation is based on a naive nearest neighbor search
  SimplexId minv = 0;
  float mind = compDist2(minv, p);
  const auto nv = _inpFieldTriangulation->getNumberOfVertices();
  for(SimplexId v = 1; v < nv; ++v) {
    const auto d = compDist2(v, p);
    if(d < mind) {
      minv = v;
      mind = d;
    }
  }
  return minv;
}

//----------------------------------------------------------------------------//

float ttk::ContourAroundPoint::compDist2(SimplexId v, SimplexId p) const {
  float vx, vy, vz;
  _inpFieldTriangulation->getVertexPoint(v, vx, vy, vz);
  const auto pCoords = &(_inpPointCoords[p * 3]);
  const float dx = pCoords[0] - vx;
  const float dy = pCoords[1] - vy;
  const float dz = pCoords[2] - vz;
  return dx * dx + dy * dy + dz * dz;
}

//----------------------------------------------------------------------------//

void ttk::ContourAroundPoint::getOutputField(
    SimplexId* &cinfos, SimplexId &nc,
    float* &coords, float* &scalars, int* &flags, SimplexId &nv) const {
  
  nc = _outNc;
  
  cinfos = new SimplexId[_outFieldCinfos.size()];
  std::copy(_outFieldCinfos.begin(), _outFieldCinfos.end(), cinfos);
  
  nv = _outFieldScalars.size();
  
  coords = new float[nv * 3];
  std::copy(_outFieldCoords.begin(), _outFieldCoords.end(), coords);
  scalars = new float[nv];
  std::copy(_outFieldScalars.begin(), _outFieldScalars.end(), scalars);
  flags = new int[nv];
  std::copy(_outFieldFlags.begin(), _outFieldFlags.end(), flags);
}

//----------------------------------------------------------------------------//

void ttk::ContourAroundPoint::getOutputPoints(
    float* &coords, float* &scalars, SimplexId &nv) const {
  
  nv = _outPointScalars.size();
  
  coords = new float[nv * 3];
  std::copy(_outPointCoords.begin(), _outPointCoords.end(), coords);
  scalars = new float[nv];
  std::copy(_outPointScalars.begin(), _outPointScalars.end(), scalars);
}
