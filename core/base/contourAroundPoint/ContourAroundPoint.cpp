#include "ContourAroundPoint.hpp"

using module = ttk::ContourAroundPoint;
using ttkIdx = ttk::SimplexId; // for Windows

//----------------------------------------------------------------------------//

int module::setInputPoints(
  float *coords, float *scalars, float *isovals, int *flags, std::size_t np) {
  _inpPtsCoords = coords;
  _inpPtsScalars = scalars;
  _inpPtsIsovals = isovals;
  _inpPtsFlags = flags;
  _inpPtsNum = np;
  if(!coords)
    return -1;
  if(!scalars)
    return -2;
  if(!isovals)
    return -3;
  if(!flags)
    return -4;
  return 0;
}

//----------------------------------------------------------------------------//

void module::getOutputContours(ttkIdx *&cinfos,
                               ttkIdx &nc,
                               float *&coords,
                               float *&scalars,
                               int *&flags,
                               ttkIdx &nv) const {

  nc = _outContoursNc;

  cinfos = new ttkIdx[_outContoursCinfos.size()];
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

void module::getOutputCentroids(float *&coords,
                                float *&scalars,
                                int *&flags,
                                ttkIdx &nv) const {

  nv = _outCentroidsScalars.size();

  coords = new float[nv * 3];
  std::copy(_outCentroidsCoords.begin(), _outCentroidsCoords.end(), coords);
  scalars = new float[nv];
  std::copy(_outCentroidsScalars.begin(), _outCentroidsScalars.end(), scalars);
  flags = new int[nv];
  std::copy(_outCentroidsFlags.begin(), _outCentroidsFlags.end(), flags);
}
