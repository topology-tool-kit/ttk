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

#include <Triangulation.h>

#include <set>
#include <vector>

namespace ttk {
class ContourAroundPoint : public Debug {

public:
  
  ContourAroundPoint() {}
  
  ~ContourAroundPoint() {}
  
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
   * see the setInput*() functions of ttk::Triangulation).
   * See ttkContourAroundPoint for further examples.
   *
   * \param triangulation Pointer to a valid triangulation.
   * \param scalars Scalar for each vertex.
   * \return 0 upon success, negative values otherwise.
   * \sa ttk::Triangulation
   */
  int setupDomain(Triangulation *triangulation, void *scalars);

  /**
   * Pass the point input data (e.g. from the wrapped algorithm).
   * @param coords 3D point coordinates in an interleaved array.
   * @param isovalues Isovalue for each contour.
   * @return 0 upon success, negative values otherwise.
   */
  int setupConstraints(float *coords, float *isovalues, int *flags,
                       std::size_t np);

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
  void getOutputField(
      SimplexId* &cinfos, SimplexId &nc,
      float* &coords, float* &scalars, int* &flags, SimplexId &nv);

  // NOTE code-clone from vtk/ttkContourAroundPoint
  /// Override this method in order to always prepend a class- and
  /// debugLevel-specific prefix and include a line end after the message.
  virtual int dMsg(std::ostream &stream, std::string msg,
                   const int &debugLevel = infoMsg) const
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

//  /**
//   * Given one of the input points, find the containing cell in the input field.
//   * The point is given by its index, and the cell is also returned by
//   * its index. If `p` happens to be one of the vertices, the first cell of
//   * its star is returned. If `p` is outside of the domain of the field,
//   * the nearest cell is returned.
//   */
//  SimplexId findCell(std::size_t p) const;
  SimplexId findInpVert(SimplexId p) const;

  /// Compute the squared distance between input field vertex v and
  /// input point p.
  float compDist2(SimplexId v, SimplexId p) const;

  template<typename scalarT>
  int handleOneInpPt(SimplexId cCenter, float isovalue, int flag) const;

//  /// Does cell `c` contain the `isovalue`?
//  template<typename scalarT>
//  bool checkContains(SimplexId c, float isovalue) const;

//  /// Enqueue the neighbors of cell `c` into the queue `q` that have not yet
//  /// been visited.
//  void enqueueNeighbors(SimplexId c, std::stack<SimplexId>& q,
//                        const std::set<SimplexId>& visited) const;

//  /// For testing: Add the cells that contain the isovalue for one point.
//  template<typename scalarT>
//  void addRoughOutput(const std::vector<SimplexId>& cellIds) const;

  /// Exctract the contour from the input edges that intersect it.
  template<typename scalarT>
  void addOutput(const std::set<SimplexId>& inpEdges,
                 float isoval, int flag) const;

  /// 3D positions of the points (as interleaved array)
  float *_inpPointCoords  = nullptr;
  /// scalar value of the points
  float *_inpPointIsovals = nullptr;
  /// 0 if the isovalue for a point is above the point's value;
  /// "point is a minimum is default"
  int *_inpPointFlags = nullptr;
  /// number of input points
  std::size_t _np = 0;

  Triangulation *_inpFieldTriangulation = nullptr; // scalar field domain
  void *_inpFieldScalars = nullptr; // scalar field image
  // up to _inpDimMax-dimensional cells appear in _inpFieldTriangulation;
  // lower-dimensional cells may appear too.
  SimplexId _inpDimMax = 0;

  // (Re)computed in every call to `execute`
  mutable std::vector<SimplexId> _outFieldCinfos;
  mutable SimplexId _outNc = 0;
  mutable std::vector<float> _outFieldCoords;
  mutable std::vector<float> _outFieldScalars;
  mutable std::vector<int> _outFieldFlags;
};
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
template<class scalarT>
int ttk::ContourAroundPoint::execute() const
{
  _outFieldCinfos.resize(0);
  _outNc = 0;
  _outFieldCoords.resize(0);
  _outFieldScalars.resize(0);
  _outFieldFlags.resize(0);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!_inpPointCoords) return -1;
  if(!_inpPointIsovals) return -2;
  if(!_inpFieldTriangulation) return -3;
  if(!_inpFieldScalars) return -4;
#endif

  Timer timUseObj;

  // The following open-mp processing is only relevant for
  // embarrassingly parallel algorithms.
//#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for num_threads(threadNumber_)
//#endif
  for(SimplexId p = 0; p < _np; ++p)
  {
    const auto pStr = "p" + std::to_string(p);
    const auto isoval = _inpPointIsovals[p];
    std::ostringstream msgStream;
    msgStream << pStr << " isoval="<<std::setw(7)<<isoval;
    msg(msgStream.str().c_str(), advancedInfoMsg);

    handleOneInpPt<scalarT>(findInpVert(p), isoval, _inpPointFlags[p]);
    // TODO check ret val
  }
   
  std::ostringstream timUseStream;
  timUseStream<< std::fixed<< std::setprecision(3)<< timUseObj.getElapsedTime()
              << " s  ("<< _np<< " points using "<< threadNumber_<< " threads)";
  dMsg(std::cout, timUseStream.str(), timeMsg);
  return 0;
}

//----------------------------------------------------------------------------//
template<typename scalarT>
int ttk::ContourAroundPoint::handleOneInpPt(SimplexId vBeg,
                                            float isovalue, int flag) const
{
  auto triangu = _inpFieldTriangulation; // just a short alias
  auto inpScalars = reinterpret_cast<const scalarT*>(_inpFieldScalars);
  const bool vBegIsMin = flag == 0;
  
  // The general idea of the algorithm is to start a (breadth-first) search
  // on the triangulation starting at the input vertex.
  // Until the isocontour is met, all vertices are collected to approximate
  // the (weighted) centroid of the enclosed area/volume.
  // The edges leading outside of the region are also collected, to compute
  // the actual isocontour (that intersects these edges).
  
  auto innerVerts = std::vector<SimplexId>{vBeg};
  std::set<SimplexId> xEdges; // set data structure is needed in addOutput
  
  struct VertEdge {
    SimplexId v; // where am i
    SimplexId e; // how did i get here
  };
  std::vector<VertEdge> q;
  std::set<SimplexId> handledEdges;
  
  auto enqueueNeighbors = [triangu, &q, &handledEdges](SimplexId vSrc) {
    const auto n = triangu->getVertexEdgeNumber(vSrc);
    for(SimplexId i = 0; i < n; ++i) {
      SimplexId e; triangu->getVertexEdge(vSrc, i, e);
      if(handledEdges.count(e))
        continue;
      handledEdges.insert(e);
      // We assume the local indexing of vertex neighbors is the same
      // as for the edges (edge 0 leads to neighbor 0, etc.).
      SimplexId vTgt; triangu->getVertexNeighbor(vSrc, i, vTgt);
      q.push_back({vTgt, e});
    }
  };
  
  enqueueNeighbors(vBeg);
  while(!q.empty())
  {
    const auto ve = q.back(); q.pop_back();
    const auto v = ve.v;
    const bool vIsAboveIso = inpScalars[v] > isovalue;
    if(vBegIsMin == vIsAboveIso) { // crossed the contour
      xEdges.insert(ve.e);
     } else {
      innerVerts.push_back(v);
      enqueueNeighbors(v);
    }
  }
  
  addOutput<scalarT>(xEdges, isovalue, flag);
  // TODO make second output
//  std::ostringstream out;
//  out << "num. inner vertices: " << innerVerts.size();
//  msg(out.str().c_str());
  return 0;
}

//----------------------------------------------------------------------------//
template<typename scalarT>
void ttk::ContourAroundPoint::addOutput(const std::set<SimplexId> &inpEdges,
                                        float isoval, int flag) const
{
  auto triangu = _inpFieldTriangulation; // just a short alias
  auto inpScalars = reinterpret_cast<const scalarT*>(_inpFieldScalars);
  
  // Every edge e of the given edges maps to one output vertex v;
  // v is connected to a number of other such vertices.
  // Zero if e is part of an "antenna", one if e is part of the boundary
  // of the grid, two in the default 2D case, possibly more in 3D.
  // These neighbors of v are reached by looking at the triangles that e is a
  // part of.
  
  std::map<SimplexId, SimplexId> inpe2outv;
  
  auto getOrMakeOutVert = [&](SimplexId e) {
    if(inpe2outv.count(e))
      return inpe2outv[e];
      
    // geometry of the vertex (linear interpolation on e)
    const SimplexId v = _outFieldScalars.size();
    inpe2outv[e] = v;
    
    SimplexId p; triangu->getEdgeVertex(e, 0, p);
    SimplexId q; triangu->getEdgeVertex(e, 1, q);
    
    const float pVal = static_cast<float>(inpScalars[p]);
    const float qVal = static_cast<float>(inpScalars[q]);
    const double pFac = (qVal - isoval) / (qVal - pVal);
    const double qFac = 1 - pFac;
    
    float px, py, pz; triangu->getVertexPoint(p, px, py, pz);
    float qx, qy, qz; triangu->getVertexPoint(q, qx, qy, qz);
    const float x = px * pFac + qx * qFac;
    const float y = py * pFac + qy * qFac;
    const float z = pz * pFac + qz * qFac;
    _outFieldCoords.push_back(x);
    _outFieldCoords.push_back(y);
    _outFieldCoords.push_back(z);
    
    _outFieldScalars.push_back(isoval);
    _outFieldFlags.push_back(flag);
    return v;
  };
  
  // Each input triangle can contribute to at most one output edge
  std::set<SimplexId> handledTris;
  
  if(_inpDimMax == 3) {
    auto msg = "! Currently only the *edges* of contour surfaces are computed";
    dMsg(std::cout, msg);
  }
  
  for(const auto e : inpEdges) {
    const auto nTriLoc = triangu->getEdgeTriangleNumber(e);
    if(nTriLoc == 0)
      continue; // we do not output "contour points" (0D cells)
    
    const auto v = getOrMakeOutVert(e);
    
    for(SimplexId tLoc = 0; tLoc < nTriLoc; ++tLoc) {
      SimplexId tGlo; triangu->getEdgeTriangle(e, tLoc, tGlo);
      if(handledTris.count(tGlo))
        continue;
      
      handledTris.insert(tGlo);
      // Find the other edge that contains the isovalue
      for(SimplexId eLoc = 0; eLoc < 3; ++eLoc) {
        SimplexId eGlo; triangu->getTriangleEdge(tGlo, eLoc, eGlo);
        if(eGlo == e || !inpEdges.count(eGlo))
          continue;
        
        // new output edge from v(e) to v(eGlo)
        _outFieldCinfos.push_back(2);
        _outFieldCinfos.push_back(v);
        _outFieldCinfos.push_back(getOrMakeOutVert(eGlo));
        ++_outNc;
      } // END target edge
    } // END triangle
  } // END source edge
}
#if 0 // TODO remove old code
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
#endif // old
