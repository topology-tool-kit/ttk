#include "ContourAroundPoint.hpp"

//#include <Triangulation.h>

int ttk::ContourAroundPoint::setupDomain(Triangulation *triangulation,
                                         void *scalars) {
  msg(std::string(60, '-').c_str());

  _inpFieldTriangulation = triangulation;
  _inpFieldScalars = scalars;
  if(!triangulation)
    return -1;
  if(!scalars)
    return -2;
  const SimplexId maxNumDim = triangulation->getDimensionality();
  _inpNvPerC = maxNumDim + 1;
  // There is no "getMinNumDim" method :-(
  const SimplexId nc = triangulation->getNumberOfCells();

#ifndef NDEBUG
  for(SimplexId c = 0; c < nc; ++c) {
    const SimplexId lnv = triangulation->getCellVertexNumber(c);
    if(lnv == _inpNvPerC)
      continue;
    std::ostringstream stream;
    stream << "Cell " << c << " has " << lnv << " vertices but it should be "
           << _inpNvPerC;
    err(stream.str().c_str());
    return -3;
  }
#endif

  {
    std::ostringstream stream;
    stream << "Number of input cells: " << nc
           << "  number of vertices per cell: " << _inpNvPerC;
    msg(stream.str().c_str());
  }

  // Call all the required precondition functions here!
  triangulation->preprocessVertexStars(); // for findCell --> getVertexStar
  triangulation
    ->preprocessCellNeighbors(); // for compOneContour --> getCellNeighbor
  triangulation->preprocessCellEdges(); // for addOutput --> getCellEdge
  triangulation->preprocessEdges(); //  i.a. for addOutput --> getEdgeVertex
  return 0;
}

//------------------------------------------------------------------------------------------------//

int ttk::ContourAroundPoint::setupConstraints(float *coords,
                                              float *isovalues,
                                              std::size_t np,
                                              int *flags) {
  _inpPointCoords = coords;
  _inpPointIsovals = isovalues;
  _np = np;
  _inpPointFlags = flags;
  if(!coords)
    return -1;
  if(!isovalues)
    return -2;
  if(!flags)
    return -3;
  return 0;
}

//------------------------------------------------------------------------------------------------//

ttk::SimplexId ttk::ContourAroundPoint::findCell(std::size_t p) const {
  // NOTE This whole method is just a hack. Eventually, the best solution would
  // probably be to use a `vtkAbstractCellLocator` in the wrapped algorithm and
  // pass an array with the cell index for each point to this module. This
  // implementation is based on a naive nearest neighbor search.
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

  SimplexId c;
#ifndef NDEBUG
  const auto errCode =
#endif
    _inpFieldTriangulation->getVertexStar(minv, 0, c);
  assert(errCode == 0);
  return c;
}

//------------------------------------------------------------------------------------------------//

float ttk::ContourAroundPoint::compDist2(SimplexId v, std::size_t p) const {
  float vx, vy, vz;
#ifndef NDEBUG
  const auto errCode =
#endif
    _inpFieldTriangulation->getVertexPoint(v, vx, vy, vz);
  assert(errCode == 0);
  const auto pCoords = &(_inpPointCoords[p * 3]);
  const float dx = pCoords[0] - vx;
  const float dy = pCoords[1] - vy;
  const float dz = pCoords[2] - vz;
  return dx * dx + dy * dy + dz * dz;
}

//------------------------------------------------------------------------------------------------//

void ttk::ContourAroundPoint::enqueueNeighbors(
  SimplexId c,
  std::stack<SimplexId> &q,
  const std::set<SimplexId> &visited) const {
  const SimplexId numNei = _inpFieldTriangulation->getCellNeighborNumber(c);
  for(SimplexId i = 0; i < numNei; ++i) {
    SimplexId nei;
#ifndef NDEBUG
    const auto errCode =
#endif
      _inpFieldTriangulation->getCellNeighbor(c, i, nei);
    assert(errCode == 0);
    if(visited.count(nei) == 0)
      q.push(nei);
  }
}

//------------------------------------------------------------------------------------------------//

void ttk::ContourAroundPoint::getOutputField(float *&coords,
                                             SimplexId &nv,
                                             SimplexId *&cinfos,
                                             SimplexId &nc,
                                             float *&scalars,
                                             int *&flags) {
  nv = _outFieldScalars.size();
  assert(SimplexId(_outFieldCoords.size()) == nv * 3);

  coords = new float[nv * 3];
  std::copy(_outFieldCoords.begin(), _outFieldCoords.end(), coords);
  scalars = new float[nv];
  std::copy(_outFieldScalars.begin(), _outFieldScalars.end(), scalars);

  const auto nvPerC = _outNvPerC;
  nc = _outFieldCinfos.size() / nvPerC;
  {
    std::ostringstream stream;
    stream << "Number of output cells: " << nc
           << "  number of vertices per cell: " << nvPerC;
    msg(stream.str().c_str());
  }
  cinfos = new SimplexId[nc * (nvPerC + 1)];
  for(SimplexId c = 0, tgt = 0, src = 0; c < nc; ++c) {
    cinfos[tgt++] = nvPerC;
    for(SimplexId i = 0; i < nvPerC; ++i)
      cinfos[tgt++] = _outFieldCinfos[src++];
  }

  flags = new int[nv];
  std::copy(_outFieldFlags.begin(), _outFieldFlags.end(), flags);
}
