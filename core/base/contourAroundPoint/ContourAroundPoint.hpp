/** \ingroup base
 * \class ttk::ContourAroundPoint
 * \author Christopher Kappe <kappe@cs.uni-kl.de>
 * \date January 25, 2019
 *
 * \brief TTK %contourAroundPoint processing package.
 *
 * %ContourAroundPoint is a TTK processing package that takes a scalar field on the input
 * and produces a scalar field on the output.
 *
 * \sa ttk::Triangulation
 * \sa ttkContourAroundPoint.cpp %for a usage example.
 */
#pragma once

#include <Triangulation.h> // needed in template method
//#include <Wrapper.h>

#include <cassert>
#include <set>
#include <stack>
#include <vector>

namespace ttk {
class ContourAroundPoint : public Debug {

public:
  
  ContourAroundPoint() {}
  
  ~ContourAroundPoint() {}

  /**
   * Pass the point input data (e.g. from the wrapped algorithm).
   * @param coords 3D point coordinates in an interleaved array.
   * @param isovalues Isovalue for each contour.
   * @return 0 upon success, negative values otherwise.
   */
  int setupConstraints(float* coords, float* isovalues, std::size_t np, int* flags);

  /**
   * Setup a (valid) triangulation object for this TTK base object.
   *
   * \pre This function should be called prior to any usage of this TTK
   * object, in a clearly distinct pre-processing step that involves no
   * traversal or computation at all. An error will be returned otherwise.
   *
   * \note It is recommended to exclude this pre-processing function from
   * any time performance measurement. Therefore, it is recommended to
   * call this function ONLY in the pre-processing steps of your program.
   * Note however, that your triangulation object must be valid when
   * calling this function (i.e. you should have filled it at this point,
   * see the setInput*() functions of ttk::Triangulation). See ttkContourAroundPoint
   * for further examples.
   *
   * \param triangulation Pointer to a valid triangulation.
   * \param scalars Scalar for each vertex.
   * \return 0 upon success, negative values otherwise.
   * \sa ttk::Triangulation
   */
  int setupDomain(Triangulation *triangulation, void *scalars);

  /**
   * Execute the package.
   * @pre setupDomain.
   * @return 0 upon success, negative values otherwise.
   */
  template<class scalarT>
  int execute() const;

  /**
   * Get the output field data (e.g. for the wrapped algorithm).
   * To be called after a successful `execute`.
   * The ownership of the pointers moves to the caller.
   * @param coords 3D Vertex coordinates as interleaved array.
   * @param nv Number of vertices.
   * @param cinfos Sequence of cell infos like `n v0 ... vn-1`.
   * @param nc Number of cells.
   * @param scalars Scalar value for each vertex.
   */
  void getOutputField(float*& coords, SimplexId& nv,
                      SimplexId*& cinfos, SimplexId& nc, float*& scalars, int*& flags);

  // NOTE code-clone from vtk/ttkContourAroundPoint
  /// Override this method in order to always prepend a class- and debugLevel-specific prefix
  /// and include a line end after the message.
  virtual int dMsg(std::ostream &stream, std::string msg, const int &debugLevel = infoMsg) const
  {
    if(debugLevel > debugLevel_ && debugLevel > ttk::globalDebugLevel_)
      return 0;
    stream << "[ttkContourAroundPoint] ";
    switch(debugLevel)
    {
    case fatalMsg: stream << "Error: "; break; // something went wrong
    case timeMsg: stream << "Time consumption: "; break; // x.yyy s
    case memoryMsg: stream << "Memory usage: "; break; // x.yyy MB
    }
    const auto res = ttk::Debug::dMsg(stream, msg, debugLevel);
    stream << std::endl;
    return res;
  }
  
  
protected:

  /**
   * Given one of the input points, find the containing cell in the input field.
   * The point is given by its index, and the cell is also returned by its index.
   * If `p` happens to be one of the vertices, the first cell of its star is returned.
   * If `p` is outside of the domain of the field, the nearest cell is returned.
   */
  SimplexId findCell(std::size_t p) const;

  /// Compute the squared distance between input field vertex v and input point p.
  float compDist2(SimplexId v, std::size_t p) const;

  /// Compute the contour for one point (given by its cell) passing the isovalue and isBelow flag.
  template<typename scalarT>
  int compOneContour(SimplexId cCenter, float isovalue, int flag) const;

  /// Does cell `c` contain the `isovalue`?
  template<typename scalarT>
  bool checkContains(SimplexId c, float isovalue) const;

  /// Enqueue the neighbors of cell `c` into the queue `q` that have not yet been visited.
  void enqueueNeighbors(SimplexId c, std::stack<SimplexId>& q,
                        const std::set<SimplexId>& visited) const;

  /// For testing: Add the cells that contain the isovalue for one point.
  template<typename scalarT>
  void addRoughOutput(const std::vector<SimplexId>& cellIds) const;

  /// Exctract the contour from the cells.
  template<typename scalarT>
  void addOutput(const std::vector<SimplexId>& cellIds, float isoval, int flag) const;

  float* _inpPointCoords  = nullptr; // 3D position of the points (as interleaved array)
  float* _inpPointIsovals = nullptr; // scalar value of the points
  std::size_t _np = 0;
  int* _inpPointFlags = nullptr; // 0 if the isovalue for a point is above the points value

  Triangulation *_inpFieldTriangulation = nullptr; // scalar field domain
  void *_inpFieldScalars = nullptr; // scalar field image
  SimplexId _inpNvPerC; // number of dimensions of the domain + 1

  // (Re)computed in every call to `execute`.
  mutable std::vector<float> _outFieldCoords;
  mutable std::vector<float> _outFieldScalars;
  mutable std::vector<SimplexId> _outFieldCinfos; // cell vertex indices (without nvPerC info)
  mutable SimplexId _outNvPerC; // _inpNvPerC - 1 if not debug output
  mutable std::vector<int> _outFieldFlags;
  // Because global vertex indices should be a gapless sequence of integers and the output field
  // will have only a subset of the input vertices, and vertices may appear in several cells,
  // we need to establish independent output vertex indices and a respective mapping.
  mutable std::map<SimplexId, SimplexId> _inpVert2outVert; // NOTE Only needed for addRoughOutput
};
}

//------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------//
template<class scalarT>
int ttk::ContourAroundPoint::execute() const
{
  _outFieldCoords.resize(0);
  _outFieldCinfos.resize(0);
  _outFieldScalars.resize(0);
  _outFieldFlags.resize(0);
  _inpVert2outVert.clear();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!_inpPointCoords) return -1;
  if(!_inpPointIsovals) return -2;
  if(!_inpFieldTriangulation) return -3;
  if(!_inpFieldScalars) return -4;
#endif

  Timer timUseObj;

  // The following open-mp processing is only relevant for embarrassingly parallel algorithms.
//#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for num_threads(threadNumber_)
//#endif
  for(std::size_t p = 0; p < _np; ++p)
  {
    const auto pStr = "p" + std::to_string(p);
    const auto isoval = _inpPointIsovals[p];
    std::ostringstream msgStream;
    msgStream << pStr << " isoval="<<std::setw(7)<<isoval;
    msg(msgStream.str().c_str(), detailedInfoMsg);

    switch(compOneContour<scalarT>(findCell(p), isoval, _inpPointFlags[p]))
    {
    case -1:
      msg(("Warning: no isocontour found around " + pStr).c_str()); break;
    case -2:
      msg(("Warning: isocontour around " + pStr + " goes only through one cell").c_str()); break;
    }
  }
   
  std::ostringstream timUseStream;
  timUseStream << std::fixed<< std::setprecision(3)<<timUseObj.getElapsedTime() << " s  ("
               << _np << " points using " << threadNumber_ << " threads)";
  dMsg(std::cout, timUseStream.str(), timeMsg);
  return 0;
}

//------------------------------------------------------------------------------------------------//
template<typename scalarT>
int ttk::ContourAroundPoint::compOneContour(SimplexId cCenter, float isovalue, int flag) const
{
  const SimplexId nc = _inpFieldTriangulation->getNumberOfCells();
  // With a DFS starting at cell `cCenter` we find the first cell that contains the isovalue
  // (one vertex is below, one above the isovalue), call it `cStart`.
  // From `cStart` we extend the isocontour by searching for all neighboring cells that also contain
  // the isovalue. The difference to the first search is that this time we don't enqueue the
  // neighbors of cells that don't contain the isovalue.
  // (But we can still use the same `visited` set.)
  // NOTE This approach assumes that the isocontour is a single connected component.
  std::set<SimplexId> visited;

  // Find an initial cell which contais the isovalue.
  SimplexId cStart = nc;
  auto q = std::stack<SimplexId>({cCenter}); // LIFO queue for DFS
  while(!q.empty())
  {
    const auto c = q.top(); q.pop();
    visited.insert(c);
    if(checkContains<scalarT>(c, isovalue)) {
      cStart = c;
      break; // We don't have to search further, we have a cell that contains the isovalue.
    }
    enqueueNeighbors(c, q, visited);
  }
  if(cStart == nc)
    return -1;

  // Compute the cells contributing to the contour starting at the initial cell.
  auto isocells = std::vector<SimplexId>{cStart};
  q = std::stack<SimplexId>(); // reset q
  enqueueNeighbors(cStart, q, visited);
  while(!q.empty())
  {
    const auto c = q.top(); q.pop();
    visited.insert(c);
    if(checkContains<scalarT>(c, isovalue)) {
      isocells.push_back(c);
      enqueueNeighbors(c, q, visited);
    }
  }
#if 0
  addRoughOutput<scalarT>(isocells);
#else
  addOutput<scalarT>(isocells, isovalue, flag);
#endif
  if(isocells.size() == 1)
    return -2;

  return 0;
}

//------------------------------------------------------------------------------------------------//
template<typename scalarT>
bool ttk::ContourAroundPoint::checkContains(SimplexId c, float isovalue) const
{
  auto inpFieldScalars = reinterpret_cast<const scalarT*>(_inpFieldScalars);
  // BEGIN debug
//  std::ostringstream stream;
//  stream << "c=" << c << " iso=" << isovalue << ": ";
//  for(SimplexId i = 0; i < _inpFieldTriangulation->getCellVertexNumber(c); ++i)
//  {
//    SimplexId v; _inpFieldTriangulation->getCellVertex(c, i, v);
//    stream << inpFieldScalars[v] << " | ";
//  }
//  msg(stream.str().c_str());
  // END debug
  SimplexId v0; _inpFieldTriangulation->getCellVertex(c, 0, v0);
  // If true, imagine a "-" on the vertex v0.
  const bool hasBelow = inpFieldScalars[v0] < isovalue;
  const SimplexId lnv = _inpFieldTriangulation->getCellVertexNumber(c);
  for(SimplexId i = 1; i < lnv; ++i)
  {
    SimplexId v;
#ifndef NDEBUG
    const auto errCode =
#endif
    _inpFieldTriangulation->getCellVertex(c, i, v);
    assert(errCode == 0);
    // If true, imagine the opoosite sign on the vertex v.
    if((inpFieldScalars[v] < isovalue) != hasBelow) {
      return true; // We don't have to look at any other vertices of `c`, we have a "-+ edge".
    }
  }
  return false;
}

//------------------------------------------------------------------------------------------------//
template<typename scalarT>
void ttk::ContourAroundPoint::addRoughOutput(const std::vector<SimplexId>& cellIds) const
{
  _outNvPerC = _inpNvPerC;
  auto vinp2vout = _inpVert2outVert;
  auto inpScalars = reinterpret_cast<const scalarT*>(_inpFieldScalars);
  float vx, vy, vz;

  for(const auto c : cellIds)
  {
    const auto lnv = _inpFieldTriangulation->getCellVertexNumber(c);
    for(SimplexId i = 0; i < lnv; ++i)
    {
      SimplexId vinp; _inpFieldTriangulation->getCellVertex(c, i, vinp);
      if(vinp2vout.count(vinp) == 0) // new output vertex
      {
        const SimplexId vout = _outFieldScalars.size();
        vinp2vout[vinp] = vout;
        _outFieldCinfos.push_back(vout);
        _outFieldScalars.push_back(inpScalars[vinp]);
        _inpFieldTriangulation->getVertexPoint(vinp, vx, vy, vz);
        _outFieldCoords.push_back(vx);
        _outFieldCoords.push_back(vy);
        _outFieldCoords.push_back(vz);
      }
      else // scalars and coords for this vertex are already set
      {
        _outFieldCinfos.push_back(vinp2vout[vinp]);
      }
    }
  }
}

//------------------------------------------------------------------------------------------------//
template<typename scalarT>
void ttk::ContourAroundPoint::addOutput(const std::vector<SimplexId>& cellIds, float isoval,
                                        int flag) const
{
  _outNvPerC = _inpNvPerC - 1;
  const SimplexId inpNePerC = _inpNvPerC == 3 ? 3 : 6;
  std::map<SimplexId, SimplexId> e2v;
  // NOTE We simply create line segments when it would be better to create line strips.

  auto inpScalars = reinterpret_cast<const scalarT*>(_inpFieldScalars);

  for(const auto c : cellIds)
  {
    for(SimplexId i = 0; i < inpNePerC; ++i)
    {
      SimplexId e; _inpFieldTriangulation->getCellEdge(c, i, e);
      if(e2v.count(e) != 0) {
        _outFieldCinfos.push_back(e2v[e]);
        continue;
      }

      SimplexId p; _inpFieldTriangulation->getEdgeVertex(e, 0, p);
      SimplexId q; _inpFieldTriangulation->getEdgeVertex(e, 1, q);
      const float pVal = static_cast<float>(inpScalars[p]);
      const float qVal = static_cast<float>(inpScalars[q]);
      if((pVal < isoval && qVal > isoval) || (pVal > isoval && qVal < isoval))
      {
        const SimplexId v = _outFieldScalars.size();
        e2v[e] = v;
        float px, py, pz; _inpFieldTriangulation->getVertexPoint(p, px, py, pz);
        float qx, qy, qz; _inpFieldTriangulation->getVertexPoint(q, qx, qy, qz);
        const double pFac = (qVal - isoval) / (qVal - pVal);
        const double qFac = 1 - pFac;
        const float x = px * pFac + qx * qFac;
        const float y = py * pFac + qy * qFac;
        const float z = pz * pFac + qz * qFac;
        _outFieldCoords.push_back(x);
        _outFieldCoords.push_back(y);
        _outFieldCoords.push_back(z);
        _outFieldCinfos.push_back(v);
        _outFieldScalars.push_back(isoval);
        _outFieldFlags.push_back(flag);
      }
    }
  }
}
