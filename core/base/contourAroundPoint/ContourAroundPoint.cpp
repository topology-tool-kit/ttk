#include "ContourAroundPoint.hpp"

#include <limits> // for nan

using module = ttk::ContourAroundPoint;


int module::setInputField(Triangulation *triangulation, void *scalars,
                          double sizeFilter, double radius) {
  msg(std::string(60, '-').c_str());
  
  if(!triangulation) return -1;
  if(!scalars) return -2;
  _inpFldTriang = triangulation;
  _inpFldScalars = scalars;
  
  _inpDimMax = triangulation->getDimensionality();
  if(_inpDimMax < 2 || _inpDimMax > 3)
    return -3;
  
  _sizeMin = _inpFldTriang->getNumberOfVertices() * sizeFilter / 10000. + 1;
  
  if(radius == -1.)
    compRadius();
  else
    _radius = radius;
  
  // Call all the required precondition functions here!
  
  // for getVertexEdgeNumber, getVertexEdge
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

int module::setInputPoints(float *coords, float *scalars, float *isovals,
                           int *flags, std::size_t np) {
  _inpPtsCoords = coords;
  _inpPtsScalars = scalars;
  _inpPtsIsovals = isovals;
  _inpPtsFlags = flags;
  _inpPtsNum = np;
  if(!coords) return -1;
  if(!scalars) return -2;
  if(!isovals) return -3;
  if(!flags) return -4;
  return 0;
}

//----------------------------------------------------------------------------//

ttk::SimplexId module::findInpFldVert(SimplexId p) const {
  // This implementation is based on a naive nearest neighbor search
  SimplexId minv = 0;
  float mind = compDist2(minv, p);
  const auto nv = _inpFldTriang->getNumberOfVertices();
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

float module::compDist2(SimplexId v, SimplexId p) const {
  float vx, vy, vz;
  _inpFldTriang->getVertexPoint(v, vx, vy, vz);
  const auto pCoords = &(_inpPtsCoords[p * 3]);
  const float dx = pCoords[0] - vx;
  const float dy = pCoords[1] - vy;
  const float dz = pCoords[2] - vz;
  return dx * dx + dy * dy + dz * dz;
}

//----------------------------------------------------------------------------//

void module::getOutputContours(
    SimplexId *&cinfos, SimplexId &nc,
    float *&coords, float *&scalars, int *&flags, SimplexId &nv) const {

  nc = _outContoursNc;

  cinfos = new SimplexId[_outContoursCinfos.size()];
  std::copy(_outContoursCinfos.begin(), _outContoursCinfos.end(), cinfos);

  nv = _outContoursScalars.size();

  coords = new float[nv * 3];
  std::copy(_outContoursCoords.begin(), _outContoursCoords.end(), coords);
  scalars = new float[nv];
  std::copy(_outContoursScalars.begin(), _outContoursScalars.end(), scalars);
  flags = new int[nv];
  std::copy(_outContoursFlags.begin(), _outContoursFlags.end(), flags);
}

//----------------------------------------------------------------------------//

void module::getOutputCentroids(
    float *&coords, float *&scalars, int* &flags, SimplexId &nv) const {

  nv = _outCentroidsScalars.size();

  coords = new float[nv * 3];
  std::copy(_outCentroidsCoords.begin(), _outCentroidsCoords.end(), coords);
  scalars = new float[nv];
  std::copy(_outCentroidsScalars.begin(), _outCentroidsScalars.end(), scalars);
  flags = new int[nv];
  std::copy(_outCentroidsFlags.begin(), _outCentroidsFlags.end(), flags);
}

//----------------------------------------------------------------------------//

double module::compRadius() {
  if(_inpFldTriang) {
    float x, y, z;
    _inpFldTriang->getVertexPoint(0, x, y, z);
    _radius = std::sqrt(x*x + y*y + z*z);
  } else
    _radius = std::numeric_limits<double>::signaling_NaN();
  return _radius;
}
