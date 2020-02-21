/** \ingroup base
 * \class ttk::ContourAroundPoint
 * \author Christopher Kappe <kappe@cs.uni-kl.de>
 * \date January 25, 2019
 *
 * \brief TTK %contourAroundPoint processing package.
 *
 * %ContourAroundPoint is a TTK processing package that takes a scalar field
 * on the input and produces a scalar field on the output.
 *
 * \sa ttk::Triangulation
 * \sa ttkContourAroundPoint.cpp %for a usage example.
 */
#pragma once

#include <Triangulation.h>

#include <cassert>
#include <cmath>
#include <functional> // see comment below
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
   * \param sizeFilter 0 --> all pass, 10000 none pass.
   * \return 0 upon success, negative values otherwise.
   * \sa ttk::Triangulation
   */
  int setInputField(Triangulation *triangulation, void *scalars,
                    double sizeFilter);

  /**
   * Input the point data (e.g. from the wrapped algorithm).
   * @param coords 3D point coordinates in an interleaved array.
   * @param scalars Scalar value for each point.
   * @param isovals Isovalue corresponding to each point.
   * @param flags isMax-flag for each point.
   * @return 0 upon success, negative values otherwise.
   */
  int setInputPoints(float *coords, float *scalars, float *isovals, int *flags,
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
   * @param cinfos Sequence of cell infos like `n v0 ... vn-1`.
   * @param nc Number of cells.
   * @param coords 3D Vertex coordinates as interleaved array.
   * @param flags isMax-flag for each vertex.
   * @param scalars Scalar value for each vertex.
   * @param nv Number of vertices.
   */
  void getOutputContours(
      SimplexId* &cinfos, SimplexId &nc,
      float* &coords, float* &scalars, int* &flags, SimplexId &nv) const;
  
  /**
   * Get the output point data (e.g. for the wrapped algorithm).
   * To be called after a successful `execute`.
   * The ownership of the pointers moves to the caller.
   * @param coords 3D point coordinates as interleaved array.
   * @param scalars Scalar value for each point.
   * @param flags isMax-flag for each point.
   * @param nv Number of vertices.
   */
  void getOutputCentroids(
      float* &coords, float* &scalars, int* &flags, SimplexId &nv) const;

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
  
  // This static function might as well go into some kind of utility file.
  // In C++>11 auto return type is possible (for a function returning a lambda).
  /**
   * @brief Get a sigmoid function.
   * @param lox Small x-value (e.g. the expected minimum input); default 0.
   * @param hix Large x-value (e.g. the expected maximum input); default 1.
   * @param slopeFac Factor to multiply with the default slope
   *        (hiy-loy)/(hix-lox); a negative value creates a "zee" function that
   *        goes from hiy to loy.
   * @param loy Infimum  y-value; default 0.
   * @param hiy Supremum y-value; default 1.
   */
  static std::function<double(double)> getSigmoidFunc(
      double lox=0., double hix=1., double slopeFac=1.,
      double loy=0., double hiy=1.) {
    const double ry = hiy - loy;
    // inflection position where y = loy + ry/2
    const double midx = (lox + hix)/2.;
    const double slope = slopeFac * ry / (hix - lox);
    const double mulWithX = -slope;
    const double addToScaledX = slope * midx;
    return [mulWithX, addToScaledX, ry, loy](double x) {
      return ry / (1. + std::exp(mulWithX*x + addToScaledX)) + loy;
    };
  }
  
  
protected:

  /// Given one of the input points, find the nearest vertex in the input field.
  /// N.B.: Typically, the points are actually vertices of the input field.
  SimplexId findInpFldVert(SimplexId p) const;

  /// Compute the squared distance between input field vertex v and
  /// input point p.
  float compDist2(SimplexId v, SimplexId p) const;

  template<typename scalarT>
  void handleOneInpPt(
      SimplexId vBeg, float isoval, int flag, float scalar) const;

  /// Exctract the contour from the input edges that intersect it.
  /// `inpEdges` may be empty.
  template<typename scalarT>
  void extendOutFld(const std::set<SimplexId>& inpEdges,
                    float isoval, int flag) const;
  
  /// Compute one point based on the vertices within one contoured region.
  /// `vertices` must contain at least one vertex.
  template<typename scalarT>
  void extendOutPts(const std::vector<SimplexId> &vertices,
                    float isoval, int flag, float scalar) const;
  
  
  /* Input data (field and points) */
  
  Triangulation *_inpFldTriang = nullptr;
  void *_inpFldScalars = nullptr; // scalar field image
  // up to _inpDimMax-dimensional cells appear in _inpFldTriang;
  // lower-dimensional cells may appear too.
  SimplexId _inpDimMax = 0;
  // minimum required output region size,
  // in number of contained input field vertices
  std::size_t _sizeMin = 0;
  
  float *_inpPtsCoords  = nullptr;
  float *_inpPtsScalars = nullptr;
  float *_inpPtsIsovals = nullptr;
  // 0 if the isovalue for a point is above the point's value;
  // "point is a minimum is default"
  int *_inpPtsFlags = nullptr;
  std::size_t _inpPtsNum = 0;
  
  
  /* Output data (contours and centroids) */
  
  mutable std::vector<SimplexId> _outContoursCinfos;
  mutable SimplexId _outContoursNc = 0;
  mutable std::vector<float> _outContoursCoords;
  mutable std::vector<float> _outContoursScalars;
  mutable std::vector<int> _outContoursFlags;
  
  mutable std::vector<float> _outCentroidsCoords;
  mutable std::vector<float> _outCentroidsScalars;
  mutable std::vector<int> _outCentroidsFlags;
};
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
template<class scalarT>
int ttk::ContourAroundPoint::execute() const
{
  _outContoursCinfos.resize(0);
  _outContoursNc = 0;
  _outContoursCoords.resize(0);
  _outContoursScalars.resize(0);
  _outContoursFlags.resize(0);
  
  _outCentroidsCoords.resize(0);
  _outCentroidsScalars.resize(0);
  _outCentroidsFlags.resize(0);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!_inpFldTriang) return -1;
  if(!_inpFldScalars) return -2;
  if(!_inpPtsCoords) return -3;
  if(!_inpPtsScalars) return -4;
  if(!_inpPtsIsovals) return -5;
  if(!_inpPtsFlags) return -6;
#endif

  Timer timUseObj;

  // The following open-mp processing is only relevant for
  // embarrassingly parallel algorithms.
//#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for num_threads(threadNumber_)
//#endif
  for(SimplexId p = 0; p < _inpPtsNum; ++p)
  {
    handleOneInpPt<scalarT>(findInpFldVert(p), _inpPtsIsovals[p],
                            _inpPtsFlags[p], _inpPtsScalars[p]);
  }
   
  std::ostringstream timUseStream;
  timUseStream<< std::fixed<< std::setprecision(3)<< timUseObj.getElapsedTime()
              << " s  ("<< _inpPtsNum<< " points using "<< threadNumber_<< " threads)";
  dMsg(std::cout, timUseStream.str(), timeMsg);
  return 0;
}

//----------------------------------------------------------------------------//
template<typename scalarT>
void ttk::ContourAroundPoint::handleOneInpPt(
    SimplexId vBeg, float isoval, int flag, float scalar) const
{
  auto triangu = _inpFldTriang; // just a short alias
  auto inpScalars = reinterpret_cast<const scalarT*>(_inpFldScalars);
  const bool vBegIsMin = flag == 0;
  
  // The general idea of the algorithm is to start a (breadth-first) search
  // on the triangulation starting at the input vertex.
  // Until the isocontour is met, all vertices are collected to approximate
  // the (weighted) centroid of the enclosed area/volume.
  // The edges leading outside of the region are also collected, to compute
  // the actual isocontour (that intersects these edges).
  
  auto innerVerts = std::vector<SimplexId>{vBeg};
  std::set<SimplexId> xEdges; // set data-structure is needed in addOutput
  
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
    const bool vIsAboveIso = inpScalars[v] > isoval;
    if(vBegIsMin == vIsAboveIso) { // crossed the contour
      xEdges.insert(ve.e);
     } else {
      innerVerts.push_back(v);
      enqueueNeighbors(v);
    }
  }
  
  if(innerVerts.size() < _sizeMin)
    return;
  extendOutFld<scalarT>(xEdges, isoval, flag);
  extendOutPts<scalarT>(innerVerts, isoval, flag, scalar);
}

//----------------------------------------------------------------------------//
template<typename scalarT>
void ttk::ContourAroundPoint::extendOutFld(
    const std::set<SimplexId> &inpEdges, float isoval, int flag) const
{
  auto triangu = _inpFldTriang; // just a short alias
  auto inpScalars = reinterpret_cast<const scalarT*>(_inpFldScalars);
  
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
    const SimplexId v = _outContoursScalars.size();
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
    _outContoursCoords.push_back(x);
    _outContoursCoords.push_back(y);
    _outContoursCoords.push_back(z);
    
    _outContoursScalars.push_back(isoval);
    _outContoursFlags.push_back(flag);
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
        _outContoursCinfos.push_back(2);
        _outContoursCinfos.push_back(v);
        _outContoursCinfos.push_back(getOrMakeOutVert(eGlo));
        ++_outContoursNc;
      } // END target edge
    } // END triangle
  } // END source edge
}

//----------------------------------------------------------------------------//
template<typename scalarT>
void ttk::ContourAroundPoint::extendOutPts(
    const std::vector<SimplexId> &vertices,
    float isoval, int flag, float scalar) const {
  
  auto inpScalars = reinterpret_cast<const scalarT*>(_inpFldScalars);
  // Vertices that are close to the isovalue have little weight.
  using ScaPair = std::pair<float,float>;
  auto scaMinMax = flag ? ScaPair{isoval, scalar} : ScaPair{scalar, isoval};
  auto wOfSca = getSigmoidFunc(scaMinMax.first, scaMinMax.second, -1.);
  
  // do the computation in double precision
  double wSum = 0.;
  double outX = 0.;
  double outY = 0.;
  double outZ = 0.;
  double outSca = 0.;
  // There are several imaginable ways to compute the scalar value at the
  // output point. We choose the weighted mean, where the weights are the same
  // as for the output coordinates.
  
  for(const auto v : vertices) {
    const scalarT sca = inpScalars[v];
//    const double w = flag ? sca - isoval : isoval - sca;
//    const double diffSca = sca - isoval;
//    const double w = diffSca * diffSca;
    const double w = wOfSca(sca);
    wSum += w;
    
    float x, y, z;
    _inpFldTriang->getVertexPoint(v, x, y, z);
    outX += x * w;
    outY += y * w;
    outZ += z * w;
    
    outSca += sca * w;
  }
  assert(wSum != 0.);
  outX /= wSum;
  outY /= wSum;
  outZ /= wSum;
  outSca /= wSum;
  
  _outCentroidsCoords.push_back(static_cast<float>(outX));
  _outCentroidsCoords.push_back(static_cast<float>(outY));
  _outCentroidsCoords.push_back(static_cast<float>(outZ));
  _outCentroidsScalars.push_back(static_cast<float>(outSca));
  _outCentroidsFlags.push_back(flag);
}
