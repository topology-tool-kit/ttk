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

#include <AbstractTriangulation.h>

#include <cassert>
#include <cmath>
#include <limits> // for nan
#include <map>
#include <set>
#include <vector>

namespace ttk {
  class ContourAroundPoint : virtual public Debug {

  public:
    ContourAroundPoint() {
      this->setDebugMsgPrefix("ContourAroundPoint");
    }

    ~ContourAroundPoint() {
    }

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
     * \param radius If all vertices lie on a shpere and the output is supposed
     *        to do so as well, pass the radius of the sphere here or -1 to have
     *        it computed. The default 0 signals that the data is not spherical.
     * \return 0 upon success, negative values otherwise.
     * \sa ttk::Triangulation
     */
    template <class triangulationType = AbstractTriangulation>
    int setInputField(triangulationType *triangulation,
                      void *scalars,
                      double sizeFilter,
                      double radius = 0.);

    /**
     * Input the point data (e.g. from the wrapped algorithm).
     * @param coords 3D point coordinates in an interleaved array.
     * @param scalars Scalar value for each point.
     * @param isovals Isovalue corresponding to each point.
     * @param flags isMax-flag for each point.
     * @param np Number of points.
     * @return 0 upon success, negative values otherwise.
     */
    int setInputPoints(float *coords,
                       float *scalars,
                       float *isovals,
                       int *flags,
                       std::size_t np);

    /**
     * Execute the package.
     * @pre setupDomain.
     * @return 0 upon success, negative values otherwise.
     */
    template <class scalarT>
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
    void getOutputContours(SimplexId *&cinfos,
                           SimplexId &nc,
                           float *&coords,
                           float *&scalars,
                           int *&flags,
                           SimplexId &nv) const;

    /**
     * Get the output point data (e.g. for the wrapped algorithm).
     * To be called after a successful `execute`.
     * The ownership of the pointers moves to the caller.
     * @param coords 3D point coordinates as interleaved array.
     * @param scalars Scalar value for each point.
     * @param flags isMax-flag for each point.
     * @param nv Number of vertices.
     */
    void getOutputCentroids(float *&coords,
                            float *&scalars,
                            int *&flags,
                            SimplexId &nv) const;

  protected:
    /// Given one of the input points, find the nearest vertex in the input
    /// field. N.B.: Typically, the points are actually vertices of the input
    /// field.
    template <class triangulationType = AbstractTriangulation>
    SimplexId findInpFldVert(SimplexId p) const;

    /// Compute the squared distance between input field vertex v and
    /// input point p.
    template <class triangulationType = AbstractTriangulation>
    float compDist2(SimplexId v, SimplexId p) const;

    template <typename scalarT, class triangulationType = AbstractTriangulation>
    void handleOneInpPt(SimplexId vBeg,
                        float isoval,
                        int flag,
                        float scalar) const;

    /// Exctract the contour from the input edges that intersect it.
    /// `inpEdges` may be empty.
    template <typename scalarT, class triangulationType = AbstractTriangulation>
    void extendOutFld(const std::set<SimplexId> &inpEdges,
                      float isoval,
                      int flag) const;

    /// Compute one point based on the vertices within one contoured region.
    /// `vertices` must contain at least one vertex.
    template <typename scalarT, class triangulationType = AbstractTriangulation>
    void extendOutPts(const std::vector<SimplexId> &vertices,
                      float isoval,
                      int flag,
                      float scalar) const;

    /// Set and return _radius
    template <class triangulationType = AbstractTriangulation>
    double compRadius() {
      if(_inpFldTriang) {
        float x, y, z;
        reinterpret_cast<triangulationType *>(_inpFldTriang)
          ->getVertexPoint(0, x, y, z);
        _radius = std::sqrt(x * x + y * y + z * z);
      } else
        _radius = std::numeric_limits<double>::signaling_NaN();
      return _radius;
    }

    /// If greater than 0, means the coordinates are supposed to lie on a sphere
    /// with fixed radius.
    double _radius = 0.;

    /* Input data (field and points) */

    void *_inpFldTriang = nullptr; // scalar field domain
    void *_inpFldScalars = nullptr; // scalar field image
    // up to _inpDimMax-dimensional cells appear in _inpFldTriang;
    // lower-dimensional cells may appear too.
    SimplexId _inpDimMax = 0;
    // minimum required output region size,
    // in number of contained input field vertices
    std::size_t _sizeMin = 0;

    float *_inpPtsCoords = nullptr;
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
} // namespace ttk

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

template <class Triang>
int ttk::ContourAroundPoint::setInputField(Triang *triangulation,
                                           void *scalars,
                                           double sizeFilter,
                                           double radius) {

  if(!triangulation)
    return -1;
  if(!scalars)
    return -2;
  _inpFldTriang = triangulation;
  _inpFldScalars = scalars;

  _inpDimMax = triangulation->getDimensionality();
  if(_inpDimMax < 2 || _inpDimMax > 3)
    return -3;

  _sizeMin = triangulation->getNumberOfVertices() * sizeFilter / 10000. + 1;

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

template <class scalarT>
int ttk::ContourAroundPoint::execute() const {
  _outContoursCinfos.resize(0);
  _outContoursNc = 0;
  _outContoursCoords.resize(0);
  _outContoursScalars.resize(0);
  _outContoursFlags.resize(0);

  _outCentroidsCoords.resize(0);
  _outCentroidsScalars.resize(0);
  _outCentroidsFlags.resize(0);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!_inpFldTriang)
    return -1;
  if(!_inpFldScalars)
    return -2;
  if(!_inpPtsCoords)
    return -3;
  if(!_inpPtsScalars)
    return -4;
  if(!_inpPtsIsovals)
    return -5;
  if(!_inpPtsFlags)
    return -6;
#endif

  Timer timUseObj;

  // The following open-mp processing is only relevant for
  // embarrassingly parallel algorithms.
  //#ifdef TTK_ENABLE_OPENMP
  //#pragma omp parallel for num_threads(threadNumber_)
  //#endif
  for(size_t p = 0; p < _inpPtsNum; ++p) {
    handleOneInpPt<scalarT>(
      findInpFldVert(p), _inpPtsIsovals[p], _inpPtsFlags[p], _inpPtsScalars[p]);
  }

  printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
  printMsg("Complete", 1, timUseObj.getElapsedTime());
  printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
  return 1;
}

//----------------------------------------------------------------------------//
template <class Triang>
float ttk::ContourAroundPoint::compDist2(SimplexId v, SimplexId p) const {
  float vx, vy, vz;
  reinterpret_cast<Triang *>(_inpFldTriang)->getVertexPoint(v, vx, vy, vz);
  const auto pCoords = &(_inpPtsCoords[p * 3]);
  const float dx = pCoords[0] - vx;
  const float dy = pCoords[1] - vy;
  const float dz = pCoords[2] - vz;
  return dx * dx + dy * dy + dz * dz;
}

//----------------------------------------------------------------------------//
template <class Triang>
ttk::SimplexId ttk::ContourAroundPoint::findInpFldVert(SimplexId p) const {
  // This implementation is based on a naive nearest neighbor search
  SimplexId minv = 0;
  float mind = compDist2(minv, p);
  const auto nv
    = reinterpret_cast<Triang *>(_inpFldTriang)->getNumberOfVertices();
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
template <typename scalarT, class Triang>
void ttk::ContourAroundPoint::handleOneInpPt(SimplexId vBeg,
                                             float isoval,
                                             int flag,
                                             float scalar) const {
  auto triang = reinterpret_cast<Triang *>(_inpFldTriang);
  auto inpScalars = reinterpret_cast<const scalarT *>(_inpFldScalars);
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

  auto enqueueNeighbors = [triang, &q, &handledEdges](SimplexId vSrc) {
    const auto n = triang->getVertexEdgeNumber(vSrc);
    for(SimplexId i = 0; i < n; ++i) {
      SimplexId e;
      triang->getVertexEdge(vSrc, i, e);
      if(handledEdges.count(e))
        continue;
      handledEdges.insert(e);
      // We assume the local indexing of vertex neighbors is the same
      // as for the edges (edge 0 leads to neighbor 0, etc.).
      SimplexId vTgt;
      triang->getVertexNeighbor(vSrc, i, vTgt);
      q.push_back({vTgt, e});
    }
  };

  enqueueNeighbors(vBeg);
  while(!q.empty()) {
    const auto ve = q.back();
    q.pop_back();
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
template <typename scalarT, class Triang>
void ttk::ContourAroundPoint::extendOutFld(const std::set<SimplexId> &inpEdges,
                                           float isoval,
                                           int flag) const {
  auto triang = reinterpret_cast<Triang *>(_inpFldTriang);
  auto inpScalars = reinterpret_cast<const scalarT *>(_inpFldScalars);

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

    SimplexId p;
    triang->getEdgeVertex(e, 0, p);
    SimplexId q;
    triang->getEdgeVertex(e, 1, q);

    const float pVal = static_cast<float>(inpScalars[p]);
    const float qVal = static_cast<float>(inpScalars[q]);
    const double pFac = (qVal - isoval) / (qVal - pVal);
    const double qFac = 1 - pFac;

    float px, py, pz;
    triang->getVertexPoint(p, px, py, pz);
    float qx, qy, qz;
    triang->getVertexPoint(q, qx, qy, qz);
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

  if(_inpDimMax == 3)
    printWrn("Currently only the *edges* of contour surfaces are computed");

  for(const auto e : inpEdges) {
    const auto nTriLoc = triang->getEdgeTriangleNumber(e);
    if(nTriLoc == 0)
      continue; // we do not output "contour points" (0D cells)

    const auto v = getOrMakeOutVert(e);

    for(SimplexId tLoc = 0; tLoc < nTriLoc; ++tLoc) {
      SimplexId tGlo;
      triang->getEdgeTriangle(e, tLoc, tGlo);
      if(handledTris.count(tGlo))
        continue;

      handledTris.insert(tGlo);
      // Find the other edge that contains the isovalue
      for(SimplexId eLoc = 0; eLoc < 3; ++eLoc) {
        SimplexId eGlo;
        triang->getTriangleEdge(tGlo, eLoc, eGlo);
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
template <typename scalarT, class Triang>
void ttk::ContourAroundPoint::extendOutPts(
  const std::vector<SimplexId> &vertices,
  float isoval,
  int flag,
  float extremeVal) const {

  auto triang = reinterpret_cast<Triang *>(_inpFldTriang);
  auto inpScalars = reinterpret_cast<const scalarT *>(_inpFldScalars);

  // We weight the vertices based on the difference to the extreme value.
  // Vertices that are close to the extreme value get much weight,
  // vertices that are close to the isovalue get little weight.
  const double dMax = std::abs(isoval - extremeVal);
#ifndef NDEBUG
  const double dMaxRelaxed = dMax * 1.001;
#endif

  // For now, we use a "full cosine" kernel.
  const double kernelInputScaler = M_PI / dMax;
  // d=0 ... cos(0)=1 ... w=1, d=dMax ... cos(pi)=-1 ... w = 0
  auto k_func = [kernelInputScaler](double d) {
    return (std::cos(d * kernelInputScaler) + 1.) / 2.;
  };

  // do the computation in double precision
  double wSum = 0.;
  double outX = 0.;
  double outY = 0.;
  double outZ = 0.;
  double outSca = 0.;

  for(const auto v : vertices) {
    const scalarT vSca = inpScalars[v];
    // Because we only use symmetric kernels for the weighting, the sign of
    // the difference does not matter.
    const double d = vSca - extremeVal;
#ifndef NDEBUG
    if(std::abs(d) > dMaxRelaxed)
      printWrn("d: " + std::to_string(d) + "  dMax: " + std::to_string(dMax));
#endif
    const double w = k_func(d);
    wSum += w;

    float x, y, z;
    triang->getVertexPoint(v, x, y, z);
    outX += x * w;
    outY += y * w;
    outZ += z * w;
    outSca += vSca * w;
  }
  assert(wSum != 0.);
  outX /= wSum;
  outY /= wSum;
  outZ /= wSum;
  outSca /= wSum;

  if(_radius > 0.) {
    const double radiusCur = std::sqrt(outX * outX + outY * outY + outZ * outZ);
    const double radiusScaler = _radius / radiusCur;
    outX *= radiusScaler;
    outY *= radiusScaler;
    outZ *= radiusScaler;
  }

  _outCentroidsCoords.push_back(static_cast<float>(outX));
  _outCentroidsCoords.push_back(static_cast<float>(outY));
  _outCentroidsCoords.push_back(static_cast<float>(outZ));
  _outCentroidsScalars.push_back(static_cast<float>(outSca));
  _outCentroidsFlags.push_back(flag);
}
