#include "ContourAroundPoint.hpp"

//----------------------------------------------------------------------------//

int ttk::ContourAroundPoint::setInputPoints(
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
