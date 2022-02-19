/// \ingroup base
/// \class ttk::FiberSurface
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2015.
///
/// \brief TTK processing package that computes fiber surfaces.
///
/// Fiber surfaces are defined as the pre-images of curves drawn in the range
/// of bivariate volumetric functions, typically on top of the continuous
/// scatterplot. Fiber surfaces generalize the segmentation features of
/// isosurfaces to bivariate data.
/// This TTK processing package implements an exact, parallel and fast
/// algorithm for fiber surface computation on (explicit or implicit)
/// tetrahedral meshes.
///
/// \b Related \b publication \n
/// "Fast and Exact Fiber Surface Extraction for Tetrahedral Meshes" \n
/// Pavol Klacansky, Julien Tierny, Hamish Carr, Zhao Geng \n
/// IEEE Transactions on Visualization and Computer Graphics, 2016.
///
/// \param dataTypeU Data type of the input first component field (char, float,
/// etc.)
/// \param dataTypeV Data type of the input second component field (char, float,
/// etc.)
///
/// \sa ttkFiberSurface.cpp %for a usage example.

#pragma once

// standard includes
#include <array>
#include <queue>

// base code includes
#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
#include <RangeDrivenOctree.h>
#endif

#include <Debug.h>
#include <Geometry.h>
#include <Triangulation.h>

namespace ttk {

  class FiberSurface : virtual public Debug {

  public:
    struct Vertex {
      bool isBasePoint_{}, isIntersectionPoint_{};
      SimplexId localId_{}, globalId_{}, polygonEdgeId_{};
      // TODO also encode the vertex ids of the triangle of the input mesh
      // where this point has been computed (for constrained triangulation)
      std::pair<SimplexId, SimplexId> meshEdge_{};
      std::array<double, 3> p_{};
      double t_{};
      std::pair<double, double> uv_{};
    };

    struct Triangle {
      std::array<SimplexId, 3> vertexIds_{};
      SimplexId tetId_{}, caseId_{}, polygonEdgeId_{};
    };

    FiberSurface();

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int buildOctree(const triangulationType *const triangulation);
#endif

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int computeContour(const std::pair<double, double> &rangePoint0,
                              const std::pair<double, double> &rangePoint1,
                              const std::vector<SimplexId> &seedTetList,
                              const triangulationType *const triangulation,
                              const SimplexId &polygonEdgeId = 0) const;

    template <class dataTypeU, class dataTypeV>
    inline int computeContour(
      const std::vector<std::pair<std::pair<double, double>,
                                  std::pair<double, double>>> &edgeList,
      const std::vector<SimplexId> &seedTetList,
      const std::vector<SimplexId> *edgeIdList = NULL) const;

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int computeSurface(const std::pair<double, double> &rangePoint0,
                              const std::pair<double, double> &rangePoint1,
                              const triangulationType *const triangulation,
                              const SimplexId &polygonEdgeId = 0) const;

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int computeSurface(const triangulationType *const triangulation);

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int
      computeSurfaceWithOctree(const std::pair<double, double> &rangePoint0,
                               const std::pair<double, double> &rangePoint1,
                               const triangulationType *const triangulation,
                               const SimplexId &polygonEdgeId) const;
#endif

    template <class dataTypeU, class dataTypeV>
    inline int finalize(const bool &mergeDuplicatedVertices = false,
                        const bool &removeSmallEdges = false,
                        const bool &edgeFlips = false,
                        const bool &intersectionRemesh = false);

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
    inline void flushOctree() {
      octree_.flush();
    }
#endif

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int processTetrahedron(const SimplexId &tetId,
                                  const std::pair<double, double> &rangePoint0,
                                  const std::pair<double, double> &rangePoint1,
                                  const triangulationType *const triangulation,
                                  const SimplexId &polygonEdgeId = 0) const;

    inline int setGlobalVertexList(std::vector<Vertex> *globalList) {
      globalVertexList_ = globalList;
      return 0;
    }

    inline int setInputField(const void *uField, const void *vField) {

      uField_ = uField;
      vField_ = vField;

      return 0;
    }

    inline int setPointMerging(const bool &onOff) {
      pointSnapping_ = onOff;
      return 0;
    }

    inline int setPointMergingThreshold(const double &threshold) {
      pointSnappingThreshold_ = threshold;
      return 0;
    }

    inline int setPointNumber(const SimplexId &number) {
      pointNumber_ = number;
      return 0;
    }

    inline int setPointSet(const float *pointSet) {
      pointSet_ = pointSet;
      return 0;
    }

    inline int setPolygon(
      const std::vector<std::pair<std::pair<double, double>,
                                  std::pair<double, double>>> *polygon) {
      polygon_ = polygon;
      return 0;
    }

    inline int setPolygonEdgeNumber(const SimplexId &polygonEdgeNumber) {
      polygonEdgeNumber_ = polygonEdgeNumber;
      polygonEdgeVertexLists_.resize(polygonEdgeNumber, NULL);
      polygonEdgeTriangleLists_.resize(polygonEdgeNumber, NULL);
      return 0;
    }

    inline int setTetList(const SimplexId *tetList) {
      tetList_ = tetList;
      return 0;
    }

    inline int
      setTetNeighbors(const std::vector<std::vector<SimplexId>> *tetNeighbors) {
      tetNeighbors_ = tetNeighbors;
      return 0;
    }

    inline int setTetNumber(const SimplexId &tetNumber) {
      tetNumber_ = tetNumber;
      return 0;
    }

    inline int setTriangleList(const SimplexId &polygonEdgeId,
                               std::vector<Triangle> *triangleList) {

#ifndef TTK_ENABLE_KAMIKAZE
      if((polygonEdgeId >= 0)
         && (polygonEdgeId < (SimplexId)polygonEdgeTriangleLists_.size()))
#endif
        polygonEdgeTriangleLists_[polygonEdgeId] = triangleList;

      return 0;
    }

    inline void
      preconditionTriangulation(AbstractTriangulation *triangulation) {
      // for breadth-first search traversals
      if(triangulation != nullptr) {
        triangulation->preconditionCellNeighbors();
      }
    }

    inline int setVertexList(const SimplexId &polygonEdgeId,
                             std::vector<Vertex> *vertexList) {

#ifndef TTK_ENABLE_KAMIKAZE
      if((polygonEdgeId >= 0)
         && (polygonEdgeId < (SimplexId)polygonEdgeVertexLists_.size()))
#endif
        polygonEdgeVertexLists_[polygonEdgeId] = vertexList;

      return 0;
    }

  protected:
    struct IntersectionTriangle {
      SimplexId caseId_{};
      // use negative values for new triangles
      SimplexId triangleId_{};
      SimplexId polygonEdgeId_{};
      // use negative values for new vertices
      std::array<SimplexId, 3> vertexIds_{};
      std::array<std::pair<double, double>, 3> uv_{};
      std::array<double, 3> t_{};
      std::array<std::array<double, 3>, 3> p_{};
      std::pair<double, double> intersection_{};
    };

    template <class dataTypeU, class dataTypeV, typename triangulationType>
    inline int computeBaseTriangle(
      const SimplexId &tetId,
      const SimplexId &localEdgeId0,
      const double &t0,
      const double &u0,
      const double &v0,
      const SimplexId &localEdgeId1,
      const double &t1,
      const double &u1,
      const double &v1,
      const SimplexId &localEdgeId2,
      const double &t2,
      const double &u2,
      const double &v2,
      std::vector<std::vector<double>> &basePoints,
      std::vector<std::pair<double, double>> &basePointProections,
      std::vector<double> &basePointParameterization,
      std::vector<std::pair<SimplexId, SimplexId>> &baseEdges,
      const triangulationType *const triangulation) const;

    template <class dataTypeU, class dataTYpeV, typename triangulationType>
    inline int computeCase0(const SimplexId &polygonEdgeId,
                            const SimplexId &tetId,
                            const SimplexId &localEdgeId0,
                            const double &t0,
                            const double &u0,
                            const double &v0,
                            const SimplexId &localEdgeId1,
                            const double &t1,
                            const double &u1,
                            const double &v1,
                            const SimplexId &localEdgeId2,
                            const double &t2,
                            const double &u2,
                            const double &v2,
                            const triangulationType *const triangulation) const;

    template <class dataTypeU, class dataTYpeV, typename triangulationType>
    inline int computeCase1(const SimplexId &polygonEdgeId,
                            const SimplexId &tetId,
                            const SimplexId &localEdgeId0,
                            const double &t0,
                            const double &u0,
                            const double &v0,
                            const SimplexId &localEdgeId1,
                            const double &t1,
                            const double &u1,
                            const double &v1,
                            const SimplexId &localEdgeId2,
                            const double &t2,
                            const double &u2,
                            const double &v2,
                            const triangulationType *const triangulation) const;

    template <class dataTypeU, class dataTYpeV, typename triangulationType>
    inline int computeCase2(const SimplexId &polygonEdgeId,
                            const SimplexId &tetId,
                            const SimplexId &localEdgeId0,
                            const double &t0,
                            const double &u0,
                            const double &v0,
                            const SimplexId &localEdgeId1,
                            const double &t1,
                            const double &u1,
                            const double &v1,
                            const SimplexId &localEdgeId2,
                            const double &t2,
                            const double &u2,
                            const double &v2,
                            const triangulationType *const triangulation) const;

    template <class dataTypeU, class dataTYpeV, typename triangulationType>
    inline int computeCase3(const SimplexId &polygonEdgeId,
                            const SimplexId &tetId,
                            const SimplexId &localEdgeId0,
                            const double &t0,
                            const double &u0,
                            const double &v0,
                            const SimplexId &localEdgeId1,
                            const double &t1,
                            const double &u1,
                            const double &v1,
                            const SimplexId &localEdgeId2,
                            const double &t2,
                            const double &u2,
                            const double &v2,
                            const triangulationType *const triangulation) const;

    template <class dataTypeU, class dataTYpeV, typename triangulationType>
    inline int computeCase4(const SimplexId &polygonEdgeId,
                            const SimplexId &tetId,
                            const SimplexId &localEdgeId0,
                            const double &t0,
                            const double &u0,
                            const double &v0,
                            const SimplexId &localEdgeId1,
                            const double &t1,
                            const double &u1,
                            const double &v1,
                            const SimplexId &localEdgeId2,
                            const double &t2,
                            const double &u2,
                            const double &v2,
                            const triangulationType *const triangulation) const;

    int computeTriangleFiber(
      const SimplexId &tetId,
      const SimplexId &triangleId,
      const std::pair<double, double> &intersection,
      const std::vector<std::vector<IntersectionTriangle>> &tetIntersections,
      std::vector<double> &pA,
      std::vector<double> &pB,
      SimplexId &pivotVertexId,
      bool &edgeFiber) const;

    int computeTriangleIntersection(
      const SimplexId &tetId,
      const SimplexId &triangleId0,
      const SimplexId &triangleId1,
      const SimplexId &polygonEdgeId0,
      const SimplexId &polygonEdgeId1,
      const std::pair<double, double> &intersection,
      SimplexId &newVertexNumber,
      SimplexId &newTriangleNumber,
      std::vector<std::vector<IntersectionTriangle>> &tetIntersections,
      std::vector<std::vector<Vertex>> &tetNewVertices) const;

    int computeTriangleIntersection(
      const SimplexId &tetId,
      const SimplexId &triangleId,
      const SimplexId &polygonEdgeId,
      const std::pair<double, double> &intersection,
      const std::vector<double> &pA,
      const std::vector<double> &pB,
      const SimplexId &pivotVertexId,
      SimplexId &newVertexNumber,
      SimplexId &newTriangleNumber,
      std::vector<std::vector<IntersectionTriangle>> &tetIntersections,
      std::vector<std::vector<Vertex>> &tetNewVertices) const;

    int createNewIntersectionTriangle(
      const SimplexId &tetId,
      const SimplexId &triangleId,
      const SimplexId &vertexId0,
      const SimplexId &vertexId1,
      const SimplexId &vertexId2,
      const std::vector<std::vector<Vertex>> &tetNewVertices,
      SimplexId &newTriangleNumber,
      std::vector<std::vector<IntersectionTriangle>> &tetIntersections,
      const std::pair<double, double> *intersection = NULL) const;

    int flipEdges() const;

    int
      flipEdges(std::vector<std::pair<SimplexId, SimplexId>> &triangles) const;

    int getNumberOfCommonVertices(
      const SimplexId &tetId,
      const SimplexId &triangleId0,
      const SimplexId &triangleId1,
      const std::vector<std::vector<IntersectionTriangle>> &tetIntersections)
      const;

    int getTriangleRangeExtremities(
      const SimplexId &tetId,
      const SimplexId &triangleId,
      const std::vector<std::vector<IntersectionTriangle>> &tetIntersections,
      std::pair<double, double> &extremity0,
      std::pair<double, double> &extremity1) const;

    bool hasDuplicatedVertices(const double *p0,
                               const double *p1,
                               const double *p2) const;

    int interpolateBasePoints(const std::vector<double> &p0,
                              const std::pair<double, double> &uv0,
                              const double &t0,
                              const std::vector<double> &p1,
                              const std::pair<double, double> &uv1,
                              const double &t1,
                              const double &t,
                              Vertex &v) const;

    bool isEdgeAngleCollapsible(
      const SimplexId &source,
      const SimplexId &destination,
      const SimplexId &pivotVertexId,
      const std::vector<std::pair<SimplexId, SimplexId>> &starNeighbors) const;

    bool isEdgeFlippable(const SimplexId &edgeVertexId0,
                         const SimplexId &edgeVertexId1,
                         const SimplexId &otherVertexId0,
                         const SimplexId &otherVertexId1) const;

    inline bool isIntersectionTriangleColinear(
      const SimplexId &tetId,
      const SimplexId &triangleId,
      const std::vector<std::vector<IntersectionTriangle>> &tetIntersections,
      const std::vector<std::vector<Vertex>> &tetNewVertices,
      const SimplexId &vertexId0,
      const SimplexId &vertexId1,
      const SimplexId &vertexId2) const {

      std::vector<std::vector<double>> points(3);
      for(int i = 0; i < 3; i++) {
        SimplexId vertexId = vertexId0;
        if(i == 1)
          vertexId = vertexId1;
        if(i == 2)
          vertexId = vertexId2;

        points[i].resize(3);
        if(vertexId >= 0) {
          for(int j = 0; j < 3; j++) {
            points[i][j] = tetIntersections[tetId][triangleId].p_[vertexId][j];
          }
        } else {
          for(int j = 0; j < 3; j++) {
            points[i][j] = tetNewVertices[tetId][(-vertexId) - 1].p_[j];
          }
        }
      }

      return Geometry::isTriangleColinear(
        points[0].data(), points[1].data(), points[2].data());
    }

    int mergeEdges(const double &distanceThreshold) const;

    int mergeVertices(const double &distanceThreshold) const;

    template <class dataTypeU, class dataTypeV>
    inline int remeshIntersections() const;

    int snapToBasePoint(const std::vector<std::vector<double>> &basePoints,
                        const std::vector<std::pair<double, double>> &uv,
                        const std::vector<double> &t,
                        Vertex &v) const;

    int snapVertexBarycentrics() const;

    int snapVertexBarycentrics(
      const SimplexId &tetId,
      const std::vector<std::pair<SimplexId, SimplexId>> &triangles) const;

    bool pointSnapping_{false};

    SimplexId pointNumber_{}, tetNumber_{}, polygonEdgeNumber_{};
    const void *uField_{}, *vField_{};
    const float *pointSet_{};
    const SimplexId *tetList_{};
    const std::vector<std::vector<SimplexId>> *tetNeighbors_{};
    std::array<SimplexId, 12> edgeImplicitEncoding_{
      0, 1, 0, 2, 0, 3, 3, 1, 2, 1, 2, 3};

    double edgeCollapseThreshold_{Geometry::powIntTen(-FLT_DIG + 2)},
      pointSnappingThreshold_{Geometry::powIntTen(-FLT_DIG + 1)};

    const std::vector<std::pair<std::pair<double, double>,
                                std::pair<double, double>>> *polygon_{};

    std::vector<Vertex> *globalVertexList_{};
    std::vector<std::vector<Vertex> *> polygonEdgeVertexLists_{};
    std::vector<std::vector<Triangle> *> polygonEdgeTriangleLists_{};

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
    RangeDrivenOctree octree_{};
#endif
  };
} // namespace ttk

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int
  ttk::FiberSurface::buildOctree(const triangulationType *const triangulation) {

  if(!uField_)
    return -1;
  if(!vField_)
    return -2;

  if(octree_.empty()) {

    octree_.setDebugLevel(debugLevel_);
    octree_.setThreadNumber(threadNumber_);
    if(!triangulation) {
      octree_.setCellList(tetList_);
      octree_.setCellNumber(tetNumber_);
      octree_.setPointList(pointSet_);
      octree_.setVertexNumber(pointNumber_);
    }
    octree_.setRange(uField_, vField_);

    octree_.build<dataTypeU, dataTypeV>(triangulation);
  }

  return 0;
}
#endif

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::computeBaseTriangle(
  const SimplexId &tetId,
  const SimplexId &localEdgeId0,
  const double &t0,
  const double &u0,
  const double &v0,
  const SimplexId &localEdgeId1,
  const double &t1,
  const double &u1,
  const double &v1,
  const SimplexId &localEdgeId2,
  const double &t2,
  const double &u2,
  const double &v2,
  std::vector<std::vector<double>> &basePoints,
  std::vector<std::pair<double, double>> &basePointProjections,
  std::vector<double> &basePointParameterization,
  std::vector<std::pair<SimplexId, SimplexId>> &baseEdges,
  const triangulationType *const triangulation) const {

  basePoints.resize(3);
  basePointProjections.resize(3);
  basePointParameterization.resize(3);
  baseEdges.resize(3);

  for(int i = 0; i < 3; i++) {

    SimplexId vertexId0 = 0, vertexId1 = 0;

    switch(i) {

      case 0:
        if(!triangulation) {
          vertexId0
            = tetList_[5 * tetId + 1 + edgeImplicitEncoding_[2 * localEdgeId0]];
          vertexId1 = tetList_[5 * tetId + 1
                               + edgeImplicitEncoding_[2 * localEdgeId0 + 1]];
        } else {
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId0], vertexId0);
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId0 + 1], vertexId1);
        }
        basePointProjections[i].first = u0;
        basePointProjections[i].second = v0;
        basePointParameterization[i] = t0;
        break;

      case 1:
        if(!triangulation) {
          vertexId0
            = tetList_[5 * tetId + 1 + edgeImplicitEncoding_[2 * localEdgeId1]];
          vertexId1 = tetList_[5 * tetId + 1
                               + edgeImplicitEncoding_[2 * localEdgeId1 + 1]];
        } else {
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId1], vertexId0);
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId1 + 1], vertexId1);
        }
        basePointProjections[i].first = u1;
        basePointProjections[i].second = v1;
        basePointParameterization[i] = t1;
        break;

      case 2:
        if(!triangulation) {
          vertexId0
            = tetList_[5 * tetId + 1 + edgeImplicitEncoding_[2 * localEdgeId2]];
          vertexId1 = tetList_[5 * tetId + 1
                               + edgeImplicitEncoding_[2 * localEdgeId2 + 1]];
        } else {
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId2], vertexId0);
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId2 + 1], vertexId1);
        }
        basePointProjections[i].first = u2;
        basePointProjections[i].second = v2;
        basePointParameterization[i] = t2;
        break;
    }

    std::vector<double> baryCentrics;
    std::vector<double> p0(2), p1(2), p(2);
    p0[0] = ((dataTypeU *)uField_)[vertexId0];
    p0[1] = ((dataTypeV *)vField_)[vertexId0];
    p1[0] = ((dataTypeU *)uField_)[vertexId1];
    p1[1] = ((dataTypeV *)vField_)[vertexId1];
    p[0] = basePointProjections[i].first;
    p[1] = basePointProjections[i].second;
    Geometry::computeBarycentricCoordinates(
      p0.data(), p1.data(), p.data(), baryCentrics, 2);

    basePoints[i].resize(3);

    float pA[3], pB[3];
    if(triangulation) {
      triangulation->getVertexPoint(vertexId0, pA[0], pA[1], pA[2]);
      triangulation->getVertexPoint(vertexId1, pB[0], pB[1], pB[2]);
    }

    for(int j = 0; j < 3; j++) {

      double c0, c1;
      if(!triangulation) {
        c0 = pointSet_[3 * vertexId0 + j];
        c1 = pointSet_[3 * vertexId1 + j];
      } else {
        c0 = pA[j];
        c1 = pB[j];
      }

      basePoints[i][j] = baryCentrics[0] * c0 + baryCentrics[1] * c1;
    }

    if(vertexId0 < vertexId1) {
      baseEdges[i] = std::pair<SimplexId, SimplexId>(vertexId0, vertexId1);
    } else {
      baseEdges[i] = std::pair<SimplexId, SimplexId>(vertexId1, vertexId0);
    }
  }

  return 0;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::computeCase0(
  const SimplexId &polygonEdgeId,
  const SimplexId &tetId,
  const SimplexId &localEdgeId0,
  const double &t0,
  const double &u0,
  const double &v0,
  const SimplexId &localEdgeId1,
  const double &t1,
  const double &u1,
  const double &v1,
  const SimplexId &localEdgeId2,
  const double &t2,
  const double &u2,
  const double &v2,
  const triangulationType *const triangulation) const {

  // that one's easy, make just one triangle
  SimplexId vertexId = (*polygonEdgeVertexLists_[polygonEdgeId]).size();

  // alloc 1 more triangle
  (*polygonEdgeTriangleLists_[polygonEdgeId])
    .resize((*polygonEdgeTriangleLists_[polygonEdgeId]).size() + 1);
  (*polygonEdgeTriangleLists_[polygonEdgeId]).back().tetId_ = tetId;
  (*polygonEdgeTriangleLists_[polygonEdgeId]).back().caseId_ = 0;
  (*polygonEdgeTriangleLists_[polygonEdgeId]).back().polygonEdgeId_
    = polygonEdgeId;

  (*polygonEdgeTriangleLists_[polygonEdgeId]).back().vertexIds_[0] = vertexId;
  (*polygonEdgeTriangleLists_[polygonEdgeId]).back().vertexIds_[1]
    = vertexId + 1;
  (*polygonEdgeTriangleLists_[polygonEdgeId]).back().vertexIds_[2]
    = vertexId + 2;

  // alloc 3 more vertices
  (*polygonEdgeVertexLists_[polygonEdgeId]).resize(vertexId + 3);
  for(int i = 0; i < 3; i++) {
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].isBasePoint_ = true;
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].isIntersectionPoint_
      = false;
  }

  // get the vertex coordinates
  for(int i = 0; i < 3; i++) {

    SimplexId vertexId0 = 0, vertexId1 = 0;

    switch(i) {
      case 0:
        if(!triangulation) {
          vertexId0
            = tetList_[5 * tetId + 1 + edgeImplicitEncoding_[2 * localEdgeId0]];
          vertexId1 = tetList_[5 * tetId + 1
                               + edgeImplicitEncoding_[2 * localEdgeId0 + 1]];
        } else {
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId0], vertexId0);
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId0 + 1], vertexId1);
        }
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].uv_.first = u0;
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].uv_.second = v0;
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].t_ = t0;
        break;

      case 1:
        if(!triangulation) {
          vertexId0
            = tetList_[5 * tetId + 1 + edgeImplicitEncoding_[2 * localEdgeId1]];
          vertexId1 = tetList_[5 * tetId + 1
                               + edgeImplicitEncoding_[2 * localEdgeId1 + 1]];
        } else {
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId1], vertexId0);
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId1 + 1], vertexId1);
        }
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].uv_.first = u1;
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].uv_.second = v1;
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].t_ = t1;
        break;

      case 2:
        if(!triangulation) {
          vertexId0
            = tetList_[5 * tetId + 1 + edgeImplicitEncoding_[2 * localEdgeId2]];
          vertexId1 = tetList_[5 * tetId + 1
                               + edgeImplicitEncoding_[2 * localEdgeId2 + 1]];
        } else {
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId2], vertexId0);
          triangulation->getCellVertex(
            tetId, edgeImplicitEncoding_[2 * localEdgeId2 + 1], vertexId1);
        }
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].uv_.first = u2;
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].uv_.second = v2;
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].t_ = t2;
        break;
    }

    std::vector<double> baryCentrics;
    std::vector<double> p0(2), p1(2), p(2);
    p0[0] = ((dataTypeU *)uField_)[vertexId0];
    p0[1] = ((dataTypeV *)vField_)[vertexId0];
    p1[0] = ((dataTypeU *)uField_)[vertexId1];
    p1[1] = ((dataTypeV *)vField_)[vertexId1];
    p[0] = (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].uv_.first;
    p[1] = (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].uv_.second;
    Geometry::computeBarycentricCoordinates(
      p0.data(), p1.data(), p.data(), baryCentrics, 2);

    float pA[3], pB[3];
    if(triangulation) {
      triangulation->getVertexPoint(vertexId0, pA[0], pA[1], pA[2]);
      triangulation->getVertexPoint(vertexId1, pB[0], pB[1], pB[2]);
    }

    for(int j = 0; j < 3; j++) {

      double c0, c1;
      if(!triangulation) {
        c0 = pointSet_[3 * vertexId0 + j];
        c1 = pointSet_[3 * vertexId1 + j];
      } else {
        c0 = pA[j];
        c1 = pB[j];
      }

      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].p_[j]
        = baryCentrics[0] * c0 + baryCentrics[1] * c1;
    }

    if(vertexId0 < vertexId1)
      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].meshEdge_
        = std::pair<SimplexId, SimplexId>(vertexId0, vertexId1);
    else
      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].meshEdge_
        = std::pair<SimplexId, SimplexId>(vertexId1, vertexId0);
  }

  // return the number of created vertices
  return 3;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::computeCase1(
  const SimplexId &polygonEdgeId,
  const SimplexId &tetId,
  const SimplexId &localEdgeId0,
  const double &t0,
  const double &u0,
  const double &v0,
  const SimplexId &localEdgeId1,
  const double &t1,
  const double &u1,
  const double &v1,
  const SimplexId &localEdgeId2,
  const double &t2,
  const double &u2,
  const double &v2,
  const triangulationType *const triangulation) const {

  SimplexId vertexId = (*polygonEdgeVertexLists_[polygonEdgeId]).size();

  // alloc 5 more vertices
  (*polygonEdgeVertexLists_[polygonEdgeId]).resize(vertexId + 5);
  for(int i = 0; i < 5; i++) {
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].isBasePoint_ = true;
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].isIntersectionPoint_
      = false;
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].meshEdge_
      = std::pair<SimplexId, SimplexId>(-1, -1);
  }

  // alloc 3 more triangles
  SimplexId triangleId = (*polygonEdgeTriangleLists_[polygonEdgeId]).size();
  (*polygonEdgeTriangleLists_[polygonEdgeId]).resize(triangleId + 3);

  for(int i = 0; i < 3; i++) {

    (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].tetId_ = tetId;
    (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].caseId_ = 1;
    (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].polygonEdgeId_
      = polygonEdgeId;

    switch(i) {
      case 0:
        (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i]
          .vertexIds_[0]
          = vertexId;
        (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i]
          .vertexIds_[1]
          = vertexId + 1;
        (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i]
          .vertexIds_[2]
          = vertexId + 2;
        break;
      case 1:
        (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i]
          .vertexIds_[0]
          = vertexId + 1;
        (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i]
          .vertexIds_[1]
          = vertexId + 2;
        (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i]
          .vertexIds_[2]
          = vertexId + 3;
        break;
      case 2:
        (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i]
          .vertexIds_[0]
          = vertexId + 2;
        (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i]
          .vertexIds_[1]
          = vertexId + 3;
        (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i]
          .vertexIds_[2]
          = vertexId + 4;
        break;
    }
  }

  // compute the base triangle vertices like in case 1
  std::vector<std::vector<double>> basePoints(3);
  std::vector<std::pair<double, double>> basePointProjections(3);
  std::vector<double> basePointParameterization(3);
  std::vector<std::pair<SimplexId, SimplexId>> baseEdges(3);

  computeBaseTriangle<dataTypeU, dataTypeV>(
    tetId, localEdgeId0, t0, u0, v0, localEdgeId1, t1, u1, v1, localEdgeId2, t2,
    u2, v2, basePoints, basePointProjections, basePointParameterization,
    baseEdges, triangulation);

  // find the pivot vertex for this case
  SimplexId pivotVertexId = -1;

  if((t0 >= 0) && (t0 <= 1)) {
    pivotVertexId = 0;
  }
  if((t1 >= 0) && (t1 <= 1)) {
    pivotVertexId = 1;
  }
  if((t2 >= 0) && (t2 <= 1)) {
    pivotVertexId = 2;
  }

  // now get the vertex coordinates
  for(int i = 0; i < 5; i++) {

    SimplexId vertexId0 = -1, vertexId1 = -1;
    double t{};

    if(!i) {
      // just take the pivot vertex
      for(int j = 0; j < 3; j++) {
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId].p_[j]
          = basePoints[pivotVertexId][j];
      }

      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId].t_
        = basePointParameterization[pivotVertexId];
      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId].uv_
        = basePointProjections[pivotVertexId];
      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId].meshEdge_
        = baseEdges[pivotVertexId];
    } else {

      switch(i) {

        case 1:
          // interpolation between pivotVertexId and (pivotVertexId-1)%3
          // if pivot is positive: interpolate at 1
          // if not interpolate at 0
          vertexId0 = pivotVertexId;
          vertexId1 = (pivotVertexId + 2) % 3;
          if(basePointParameterization[(pivotVertexId + 2) % 3] > 1) {
            t = 1;
          } else {
            t = 0;
          }
          break;

        case 2:
          // interpolation between pivotVertexId and (pivotVertexId+1)%3
          // if pivot is positive: interpolate at 1
          // if not interpolate at 0
          vertexId0 = pivotVertexId;
          vertexId1 = (pivotVertexId + 1) % 3;
          if(basePointParameterization[(pivotVertexId + 1) % 3] > 1) {
            t = 1;
          } else {
            t = 0;
          }
          break;

        case 3:
          vertexId0 = (pivotVertexId + 2) % 3;
          vertexId1 = (pivotVertexId + 1) % 3;
          if(basePointParameterization[(pivotVertexId + 2) % 3] < 0) {
            t = 0;
          } else {
            t = 1;
          }
          break;

        case 4:
          vertexId0 = (pivotVertexId + 2) % 3;
          vertexId1 = (pivotVertexId + 1) % 3;
          if(basePointParameterization[(pivotVertexId + 2) % 3] < 0) {
            t = 1;
          } else {
            t = 0;
          }
          break;
      }

      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].t_ = t;

      interpolateBasePoints(
        basePoints[vertexId0], basePointProjections[vertexId0],
        basePointParameterization[vertexId0], basePoints[vertexId1],
        basePointProjections[vertexId1], basePointParameterization[vertexId1],
        t, (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i]);
      //       snapToBasePoint(
      //         basePoints, basePointProjections, basePointParameterization,
      //         (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i]);
    }
  }

  // return the number of created vertices
  return 5;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::computeCase2(
  const SimplexId &polygonEdgeId,
  const SimplexId &tetId,
  const SimplexId &localEdgeId0,
  const double &t0,
  const double &u0,
  const double &v0,
  const SimplexId &localEdgeId1,
  const double &t1,
  const double &u1,
  const double &v1,
  const SimplexId &localEdgeId2,
  const double &t2,
  const double &u2,
  const double &v2,
  const triangulationType *const triangulation) const {

  SimplexId vertexId = (*polygonEdgeVertexLists_[polygonEdgeId]).size();

  // alloc 4 more vertices
  (*polygonEdgeVertexLists_[polygonEdgeId]).resize(vertexId + 4);
  for(int i = 0; i < 4; i++) {
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].isBasePoint_ = true;
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].isIntersectionPoint_
      = false;
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].meshEdge_
      = std::pair<SimplexId, SimplexId>(-1, -1);
  }

  // alloc 2 more triangles
  SimplexId triangleId = (*polygonEdgeTriangleLists_[polygonEdgeId]).size();
  (*polygonEdgeTriangleLists_[polygonEdgeId]).resize(triangleId + 2);

  for(int i = 0; i < 2; i++) {

    (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].tetId_ = tetId;
    (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].caseId_ = 2;
    (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].polygonEdgeId_
      = polygonEdgeId;

    if(!i) {
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[0]
        = vertexId;
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[1]
        = vertexId + 1;
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[2]
        = vertexId + 2;
    } else {
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[0]
        = vertexId + 1;
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[1]
        = vertexId + 3;
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[2]
        = vertexId + 2;
    }
  }

  // compute the base triangle vertices like in case 1
  std::vector<std::vector<double>> basePoints(3);
  std::vector<std::pair<double, double>> basePointProjections(3);
  std::vector<double> basePointParameterization(3);
  std::vector<std::pair<SimplexId, SimplexId>> baseEdges(3);

  computeBaseTriangle<dataTypeU, dataTypeV>(
    tetId, localEdgeId0, t0, u0, v0, localEdgeId1, t1, u1, v1, localEdgeId2, t2,
    u2, v2, basePoints, basePointProjections, basePointParameterization,
    baseEdges, triangulation);

  // find the pivot for this case
  bool isPivotPositive = false;
  SimplexId pivotVertexId = -1;

  if(((t0 < 0) && ((t1 < 0) || (t2 < 0)))
     || ((t1 < 0) && ((t0 < 0) || (t2 < 0)))
     || ((t2 < 0) && ((t1 < 0) || (t0 < 0)))) {
    isPivotPositive = true;
  }
  if(isPivotPositive) {
    if(t0 >= 1)
      pivotVertexId = 0;
  } else {
    if(t0 <= 0)
      pivotVertexId = 0;
  }
  if(isPivotPositive) {
    if(t1 >= 1)
      pivotVertexId = 1;
  } else {
    if(t1 <= 0)
      pivotVertexId = 1;
  }
  if(isPivotPositive) {
    if(t2 >= 1)
      pivotVertexId = 2;
  } else {
    if(t2 <= 0)
      pivotVertexId = 2;
  }

  // now get the vertex coordinates
  for(int i = 0; i < 4; i++) {

    SimplexId vertexId0 = -1, vertexId1 = -1;
    double t{};

    switch(i) {

      case 0:
        // interpolation between pivotVertexId and (pivotVertexId-1)%3
        // if pivot is positive: interpolate at 1
        // if not interpolate at 0
        vertexId0 = pivotVertexId;
        vertexId1 = (pivotVertexId + 2) % 3;
        if(isPivotPositive)
          t = 1;
        else
          t = 0;
        break;

      case 1:
        // interpolation between pivotVertexId and (pivotVertexId+1)%3
        // if pivot is positive: interpolate at 1
        // if not interpolate at 0
        vertexId0 = pivotVertexId;
        vertexId1 = (pivotVertexId + 1) % 3;
        if(isPivotPositive)
          t = 1;
        else
          t = 0;
        break;

      case 2:
        // interpolation between pivotVertexId and (pivotVertexId-1)%3
        // if pivot is positive: interpolate at 0
        // if not interpolate at 1
        vertexId0 = pivotVertexId;
        vertexId1 = (pivotVertexId + 2) % 3;
        if(isPivotPositive)
          t = 0;
        else
          t = 1;
        break;

      case 3:
        // interpolation between pivotVertexId and (pivotVertexId+1)%3
        // if pivot is positive: interpolate at 0
        // if not interpolate at 1
        vertexId0 = pivotVertexId;
        vertexId1 = (pivotVertexId + 1) % 3;
        if(isPivotPositive)
          t = 0;
        else
          t = 1;
        break;
    }

    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].t_ = t;

    interpolateBasePoints(
      basePoints[vertexId0], basePointProjections[vertexId0],
      basePointParameterization[vertexId0], basePoints[vertexId1],
      basePointProjections[vertexId1], basePointParameterization[vertexId1], t,
      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i]);
    //     snapToBasePoint(
    //       basePoints, basePointProjections, basePointParameterization,
    //       (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i]);
  }

  // return the number of created vertices
  return 4;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::computeCase3(
  const SimplexId &polygonEdgeId,
  const SimplexId &tetId,
  const SimplexId &localEdgeId0,
  const double &t0,
  const double &u0,
  const double &v0,
  const SimplexId &localEdgeId1,
  const double &t1,
  const double &u1,
  const double &v1,
  const SimplexId &localEdgeId2,
  const double &t2,
  const double &u2,
  const double &v2,
  const triangulationType *const triangulation) const {

  SimplexId vertexId = (*polygonEdgeVertexLists_[polygonEdgeId]).size();

  // alloc 3 more vertices
  (*polygonEdgeVertexLists_[polygonEdgeId]).resize(vertexId + 3);
  for(int i = 0; i < 3; i++) {
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].isBasePoint_ = true;
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].isIntersectionPoint_
      = false;
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].meshEdge_
      = std::pair<SimplexId, SimplexId>(-1, -1);
  }

  // alloc 1 more triangle
  SimplexId triangleId = (*polygonEdgeTriangleLists_[polygonEdgeId]).size();
  (*polygonEdgeTriangleLists_[polygonEdgeId]).resize(triangleId + 1);

  (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId].tetId_ = tetId;
  (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId].caseId_ = 3;
  (*polygonEdgeTriangleLists_[polygonEdgeId]).back().polygonEdgeId_
    = polygonEdgeId;

  (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId].vertexIds_[0]
    = vertexId;
  (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId].vertexIds_[1]
    = vertexId + 1;
  (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId].vertexIds_[2]
    = vertexId + 2;

  // compute the base triangle vertices like in case 1
  std::vector<std::vector<double>> basePoints(3);
  std::vector<std::pair<double, double>> basePointProjections(3);
  std::vector<double> basePointParameterization(3);
  std::vector<std::pair<SimplexId, SimplexId>> baseEdges(3);

  computeBaseTriangle<dataTypeU, dataTypeV>(
    tetId, localEdgeId0, t0, u0, v0, localEdgeId1, t1, u1, v1, localEdgeId2, t2,
    u2, v2, basePoints, basePointProjections, basePointParameterization,
    baseEdges, triangulation);

  // now find the pivot
  bool isPivotPositive = false;
  SimplexId pivotVertexId = -1;

  if((t0 <= 1) && (t0 >= 0)) {
    pivotVertexId = 0;
  } else {
    if(t0 < 0)
      isPivotPositive = true;
  }
  if((t1 <= 1) && (t1 >= 0)) {
    pivotVertexId = 1;
  } else {
    if(t1 < 0)
      isPivotPositive = true;
  }
  if((t2 <= 1) && (t2 >= 0)) {
    pivotVertexId = 2;
  } else {
    if(t2 < 0)
      isPivotPositive = true;
  }

  // now get the vertex coordinates
  for(int i = 0; i < 3; i++) {

    SimplexId vertexId0 = -1, vertexId1 = -1;
    double t{};

    if(!i) {
      // special case of the pivot vertex
      for(int j = 0; j < 3; j++) {
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId].p_[j]
          = basePoints[pivotVertexId][j];
      }

      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId].t_
        = basePointParameterization[pivotVertexId];
      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId].uv_
        = basePointProjections[pivotVertexId];
      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId].meshEdge_
        = baseEdges[pivotVertexId];
    } else {
      if(i == 1) {
        // interpolation between pivotVertexId and pivotVertexId+1
        // if pivot is positive, interpolate at value 0
        // else interpolate at value 1
        vertexId0 = pivotVertexId;
        vertexId1 = (pivotVertexId + 1) % 3;
        if(isPivotPositive)
          t = 0;
        else
          t = 1;
      }
      if(i == 2) {
        // interpolation between pivotVertexId and pivotVertexId-1
        // if pivot is positive, interpolate at value 0
        // else interpolate at value 1
        vertexId0 = pivotVertexId;
        vertexId1 = (pivotVertexId + 2) % 3;
        if(isPivotPositive)
          t = 0;
        else
          t = 1;
      }

      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].t_ = t;

      interpolateBasePoints(
        basePoints[vertexId0], basePointProjections[vertexId0],
        basePointParameterization[vertexId0], basePoints[vertexId1],
        basePointProjections[vertexId1], basePointParameterization[vertexId1],
        t, (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i]);
      //       snapToBasePoint(
      //         basePoints, basePointProjections, basePointParameterization,
      //         (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i]);
    }
  }

  // return the number of created vertices
  return 3;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::computeCase4(
  const SimplexId &polygonEdgeId,
  const SimplexId &tetId,
  const SimplexId &localEdgeId0,
  const double &t0,
  const double &u0,
  const double &v0,
  const SimplexId &localEdgeId1,
  const double &t1,
  const double &u1,
  const double &v1,
  const SimplexId &localEdgeId2,
  const double &t2,
  const double &u2,
  const double &v2,
  const triangulationType *const triangulation) const {

  SimplexId vertexId = (*polygonEdgeVertexLists_[polygonEdgeId]).size();

  // alloc 4 more vertices
  (*polygonEdgeVertexLists_[polygonEdgeId]).resize(vertexId + 4);
  for(int i = 0; i < 4; i++) {
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].isBasePoint_ = true;
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].isIntersectionPoint_
      = false;
    (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].meshEdge_
      = std::pair<SimplexId, SimplexId>(-1, -1);
  }

  // alloc 2 more triangles
  SimplexId triangleId = (*polygonEdgeTriangleLists_[polygonEdgeId]).size();
  (*polygonEdgeTriangleLists_[polygonEdgeId]).resize(triangleId + 2);

  for(int i = 0; i < 2; i++) {

    (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].tetId_ = tetId;
    (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].caseId_ = 4;
    (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].polygonEdgeId_
      = polygonEdgeId;

    if(!i) {
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[0]
        = vertexId;
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[1]
        = vertexId + 1;
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[2]
        = vertexId + 2;
    } else {
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[0]
        = vertexId + 1;
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[1]
        = vertexId + 3;
      (*polygonEdgeTriangleLists_[polygonEdgeId])[triangleId + i].vertexIds_[2]
        = vertexId + 2;
    }
  }

  // compute the base triangle vertices like in case 1
  std::vector<std::vector<double>> basePoints(3);
  std::vector<std::pair<double, double>> basePointProjections(3);
  std::vector<double> basePointParameterization(3);
  std::vector<std::pair<SimplexId, SimplexId>> baseEdges(3);

  computeBaseTriangle<dataTypeU, dataTypeV>(
    tetId, localEdgeId0, t0, u0, v0, localEdgeId1, t1, u1, v1, localEdgeId2, t2,
    u2, v2, basePoints, basePointProjections, basePointParameterization,
    baseEdges, triangulation);

  // find the pivot vertex for this case
  bool isPivotPositive = false;
  SimplexId pivotVertexId = -1;

  if(t0 > 1) {
    pivotVertexId = 0;
    isPivotPositive = true;
  } else if(t0 < 0) {
    pivotVertexId = 0;
    isPivotPositive = false;
  }

  if(t1 > 1) {
    pivotVertexId = 1;
    isPivotPositive = true;
  } else if(t1 < 0) {
    pivotVertexId = 1;
    isPivotPositive = false;
  }
  if(t2 > 1) {
    pivotVertexId = 2;
    isPivotPositive = true;
  } else if(t2 < 0) {
    pivotVertexId = 2;
    isPivotPositive = false;
  }

  // now get the vertex coordinates depending on the case
  for(int i = 0; i < 4; i++) {

    SimplexId vertexId0 = -1, vertexId1 = -1;
    double t{};

    if(i < 2) {
      if(!i) {
        // interpolation between pivotVertexId and (pivotVertexId-1)%3
        // if pivot is positive: interpolate at 1
        // if not interpolate at 0
        vertexId0 = pivotVertexId;
        vertexId1 = (pivotVertexId + 2) % 3;
        if(isPivotPositive)
          t = 1;
        else
          t = 0;
      } else {
        // interpolation between pivotVertexId and (pivotVertexId+1)%3
        // if pivot is positive: interpolate at 1
        // if not interpolate at 0
        vertexId0 = pivotVertexId;
        vertexId1 = (pivotVertexId + 1) % 3;
        if(isPivotPositive)
          t = 1;
        else
          t = 0;
      }

      (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].t_ = t;

      interpolateBasePoints(
        basePoints[vertexId0], basePointProjections[vertexId0],
        basePointParameterization[vertexId0], basePoints[vertexId1],
        basePointProjections[vertexId1], basePointParameterization[vertexId1],
        t, (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i]);
      //       snapToBasePoint(
      //         basePoints, basePointProjections, basePointParameterization,
      //         (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i]);

    } else {
      if(i == 2) {
        // take (pivotVertexId-1)%3
        for(int j = 0; j < 3; j++) {
          (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].p_[j]
            = basePoints[(pivotVertexId + 2) % 3][j];
        }

        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].t_
          = basePointParameterization[(pivotVertexId + 2) % 3];
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].uv_
          = basePointProjections[(pivotVertexId + 2) % 3];
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].meshEdge_
          = baseEdges[(pivotVertexId + 2) % 3];
      } else {
        // take (pivtoVertexId+1)%3
        for(int j = 0; j < 3; j++) {
          (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].p_[j]
            = basePoints[(pivotVertexId + 1) % 3][j];
        }
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].t_
          = basePointParameterization[(pivotVertexId + 1) % 3];
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].uv_
          = basePointProjections[(pivotVertexId + 1) % 3];
        (*polygonEdgeVertexLists_[polygonEdgeId])[vertexId + i].meshEdge_
          = baseEdges[(pivotVertexId + 1) % 3];
      }
    }
  }

  // return the number of created vertices
  return 4;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::computeContour(
  const std::pair<double, double> &rangePoint0,
  const std::pair<double, double> &rangePoint1,
  const std::vector<SimplexId> &seedTetList,
  const triangulationType *const triangulation,
  const SimplexId &polygonEdgeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(!uField_)
    return -4;
  if(!vField_)
    return -5;
  if(!triangulation)
    return -6;
  if(!polygonEdgeNumber_)
    return -7;
  if(!globalVertexList_)
    return -8;
#endif

  std::vector<bool> visitedTets(triangulation->getNumberOfCells(), false);
  std::queue<SimplexId> tetQueue;

  // init the queue
  for(SimplexId i = 0; i < (SimplexId)seedTetList.size(); i++) {
    tetQueue.push(seedTetList[i]);
  }

  SimplexId createdVertices = 0;

  do {

    SimplexId tetId = tetQueue.front();
    tetQueue.pop();

    if(!visitedTets[tetId]) {

      createdVertices = processTetrahedron<dataTypeU, dataTypeV>(
        tetId, rangePoint0, rangePoint1, triangulation, polygonEdgeId);

      if(createdVertices) {
        // only propagate if we created a triangle
        SimplexId tetNeighborNumber
          = triangulation->getCellNeighborNumber(tetId);

        for(SimplexId i = 0; i < tetNeighborNumber; i++) {

          SimplexId neighborId = -1;
          triangulation->getCellNeighbor(tetId, i, neighborId);

          if(!visitedTets[neighborId])
            tetQueue.push(neighborId);
        }
      }

      visitedTets[tetId] = true;
    }

  } while(tetQueue.size());

  return 0;
}

template <class dataTypeU, class dataTypeV>
int ttk::FiberSurface::computeContour(
  const std::vector<
    std::pair<std::pair<double, double>, std::pair<double, double>>> &edgeList,
  const std::vector<SimplexId> &seedTetList,
  const std::vector<SimplexId> *edgeIdList) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(!tetNeighbors_)
    return -1;
  if(!tetNumber_)
    return -2;
  if(!tetList_)
    return -3;
  if(!uField_)
    return -4;
  if(!vField_)
    return -5;
  if(!pointSet_)
    return -6;
  if(!polygonEdgeNumber_)
    return -7;
  if(!globalVertexList_)
    return -8;
#endif

  std::vector<bool> visitedTets(tetNumber_, false);
  std::queue<SimplexId> tetQueue;

  // init the queue
  for(SimplexId i = 0; i < (SimplexId)seedTetList.size(); i++) {
    tetQueue.push(seedTetList[i]);
  }

  SimplexId createdVertices = 0;

  do {

    SimplexId tetId = tetQueue.front();
    tetQueue.pop();

    if(!visitedTets[tetId]) {

      std::vector<std::vector<SimplexId>> threadedTetQueue(edgeList.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(SimplexId i = 0; i < (SimplexId)edgeList.size(); i++) {

        SimplexId polygonEdgeId = 0;

        if(edgeIdList) {
          polygonEdgeId = (*edgeIdList)[i];
        }

        createdVertices = processTetrahedron<dataTypeU, dataTypeV>(
          tetId, edgeList[i].first, edgeList[i].second, polygonEdgeId);

        if(createdVertices) {
          // only propagate if we created a triangle
          for(SimplexId j = 0; j < (SimplexId)(*tetNeighbors_)[tetId].size();
              j++) {
            if(!visitedTets[(*tetNeighbors_)[tetId][j]]) {
              threadedTetQueue[i].push_back((*tetNeighbors_)[tetId][j]);
            }
          }
        }
      }

      visitedTets[tetId] = true;

      for(SimplexId i = 0; i < (SimplexId)threadedTetQueue.size(); i++) {
        for(SimplexId j = 0; j < (SimplexId)threadedTetQueue[i].size(); j++) {
          tetQueue.push(threadedTetQueue[i][j]);
        }
      }
    }

  } while(tetQueue.size());

  return 0;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::computeSurface(
  const std::pair<double, double> &rangePoint0,
  const std::pair<double, double> &rangePoint1,
  const triangulationType *const triangulation,
  const SimplexId &polygonEdgeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if((!tetNumber_) && (!triangulation))
    return -1;
  if((!tetList_) && (!triangulation))
    return -2;
  if(!uField_)
    return -3;
  if(!vField_)
    return -4;
  if((!pointSet_) && (!triangulation))
    return -5;
  if(!polygonEdgeNumber_)
    return -6;
  if(!globalVertexList_)
    return -7;
#endif

  SimplexId tetNumber = tetNumber_;

  if(triangulation) {
    tetNumber = triangulation->getNumberOfCells();
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < tetNumber; i++) {

    processTetrahedron<dataTypeU, dataTypeV>(
      i, rangePoint0, rangePoint1, triangulation, polygonEdgeId);
  }

  return 0;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::computeSurface(
  const triangulationType *const triangulation) {

#ifndef TTK_ENABLE_KAMIKAZE
  if((!tetNumber_) && (!triangulation))
    return -1;
  if((!tetList_) && (!triangulation))
    return -2;
  if(!uField_)
    return -3;
  if(!vField_)
    return -4;
  if((!pointSet_) && (!triangulation))
    return -5;
  if(!polygon_)
    return -6;
  if(polygonEdgeNumber_ != (SimplexId)polygon_->size())
    return -7;
  if(!globalVertexList_)
    return -8;
#endif

  Timer t;

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
  if(!octree_.empty()) {

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < polygonEdgeNumber_; i++) {

      computeSurfaceWithOctree<dataTypeU, dataTypeV>(
        (*polygon_)[i].first, (*polygon_)[i].second, triangulation, i);
    }
  } else {
    // regular extraction (the octree has not been computed)
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < polygonEdgeNumber_; i++) {
      computeSurface<dataTypeU, dataTypeV>(
        (*polygon_)[i].first, (*polygon_)[i].second, triangulation, i);
    }
  }

#else
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < polygonEdgeNumber_; i++) {

    computeSurface<dataTypeU, dataTypeV>(
      (*polygon_)[i].first, (*polygon_)[i].second, i);
  }
#endif

  finalize<dataTypeU, dataTypeV>(pointSnapping_, false, false, false);

  this->printMsg("Extracted", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::computeSurfaceWithOctree(
  const std::pair<double, double> &rangePoint0,
  const std::pair<double, double> &rangePoint1,
  const triangulationType *const triangulation,
  const SimplexId &polygonEdgeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if((!tetNumber_) && (!triangulation))
    return -1;
  if((!tetList_) && (!triangulation))
    return -2;
  if(!uField_)
    return -3;
  if(!vField_)
    return -4;
  if((!pointSet_) && (!triangulation))
    return -5;
  if(!polygonEdgeNumber_)
    return -6;
  if(!globalVertexList_)
    return -7;
#endif

  std::vector<SimplexId> tetList;
  octree_.rangeSegmentQuery(rangePoint0, rangePoint1, tetList);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)tetList.size(); i++) {
    processTetrahedron<dataTypeU, dataTypeV>(
      tetList[i], rangePoint0, rangePoint1, triangulation, polygonEdgeId);
  }

  return 0;
}
#endif

template <class dataTypeU, class dataTypeV>
int ttk::FiberSurface::finalize(const bool &mergeDuplicatedVertices,
                                const bool &removeSmallEdges,
                                const bool &edgeFlips,
                                const bool &intersectionRemesh) {

  // make only one vertex list
  SimplexId fiberSurfaceVertexNumber = 0;
  for(SimplexId i = 0; i < (SimplexId)polygonEdgeVertexLists_.size(); i++) {
    fiberSurfaceVertexNumber += (*polygonEdgeVertexLists_[i]).size();
  }

  (*globalVertexList_).resize(fiberSurfaceVertexNumber);
  fiberSurfaceVertexNumber = 0;
  for(SimplexId i = 0; i < (SimplexId)polygonEdgeVertexLists_.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)polygonEdgeVertexLists_[i]->size();
        j++) {
      (*polygonEdgeVertexLists_[i])[j].polygonEdgeId_ = i;
      (*polygonEdgeVertexLists_[i])[j].localId_ = j;
      (*polygonEdgeVertexLists_[i])[j].globalId_ = fiberSurfaceVertexNumber;
      (*globalVertexList_)[fiberSurfaceVertexNumber]
        = (*polygonEdgeVertexLists_[i])[j];
      fiberSurfaceVertexNumber++;
    }
  }
  for(SimplexId i = 0; i < (SimplexId)polygonEdgeTriangleLists_.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)polygonEdgeTriangleLists_[i]->size();
        j++) {
      for(int k = 0; k < 3; k++) {
        (*polygonEdgeTriangleLists_[i])[j].vertexIds_[k]
          = (*polygonEdgeVertexLists_[i])[(*polygonEdgeTriangleLists_[i])[j]
                                            .vertexIds_[k]]
              .globalId_;
      }
    }
  }

  // NOTE:
  // 1) in a first step, make everything work in a POST process.
  // this is what really matters for tet gen.
  // 2) probably come back and apply in a pre-process (more degenerate
  // intersection cases).
  if(intersectionRemesh) {
    remeshIntersections<dataTypeU, dataTypeV>();
  }

  if((mergeDuplicatedVertices) || (removeSmallEdges)) {
    //     we need to have a complex to perform the edge collapses
    mergeVertices(pointSnappingThreshold_);
  }

  if(edgeFlips)
    flipEdges();

  if(removeSmallEdges)
    mergeEdges(edgeCollapseThreshold_);

  // now we can release the memory for the threaded vertices
  for(SimplexId i = 0; i < (SimplexId)polygonEdgeVertexLists_.size(); i++) {
    polygonEdgeVertexLists_[i]->clear();
  }

  return 0;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
inline int ttk::FiberSurface::processTetrahedron(
  const SimplexId &tetId,
  const std::pair<double, double> &rangePoint0,
  const std::pair<double, double> &rangePoint1,
  const triangulationType *const triangulation,
  const SimplexId &polygonEdgeId) const {

  double rangeEdge[2];
  rangeEdge[0] = rangePoint0.first - rangePoint1.first;
  rangeEdge[1] = rangePoint0.second - rangePoint1.second;

  double rangeNormal[2];
  rangeNormal[0] = -rangeEdge[1];
  rangeNormal[1] = rangeEdge[0];

  const double prec_dbl = Geometry::powInt(10.0, -DBL_DIG);

  // 1. compute the distance to the range line carrying the saddleEdge
  SimplexId upperNumber = 0;
  SimplexId lowerNumber = 0;
  SimplexId equalVertexLocalId = -1;
  double d[4];
  for(int i = 0; i < 4; i++) {

    SimplexId vertexId = 0;
    if(!triangulation) {
      vertexId = tetList_[5 * tetId + 1 + i];
    } else {
      triangulation->getCellVertex(tetId, i, vertexId);
    }

    double projectedVertex[2];
    projectedVertex[0] = ((dataTypeU *)uField_)[vertexId];
    projectedVertex[1] = ((dataTypeV *)vField_)[vertexId];

    double vertexRangeEdge[2];
    vertexRangeEdge[0] = projectedVertex[0] - rangePoint0.first;
    vertexRangeEdge[1] = projectedVertex[1] - rangePoint0.second;

    d[i] = vertexRangeEdge[0] * rangeNormal[0]
           + vertexRangeEdge[1] * rangeNormal[1];

    if(fabs(d[i]) < prec_dbl)
      d[i] = 0;

    if(d[i] > 0)
      upperNumber++;
    if(d[i] < 0)
      lowerNumber++;

    if(d[i] == 0) {
      equalVertexLocalId = i;
    }
  }

  // 2. compute the base triangle(s)
  if(!((upperNumber == 0) || (lowerNumber == 0))) {

    // the fiber surface is passing through this tetrahedron.
    std::vector<bool> lonelyVertex(4, false);
    std::vector<SimplexId> triangleEdgeNumbers(2, 0);
    std::vector<std::vector<SimplexId>> triangleEdges(2);
    triangleEdges[0].resize(3, -1);
    triangleEdges[1].resize(3, -1);

    // implicit edge encoding
    // 0: O-1
    // 1: 0-2
    // 2: 0-3
    // 3: 3-1 [order!]
    // 4: 2-1 [order!]
    // 5: 2-3
    SimplexId edgeCounter = 0;
    for(int i = 0; i < 4; i++) {

      SimplexId jStart = i + 1;
      SimplexId jEnd = 4;
      SimplexId jStep = 1;

      if(i == 1) {
        // special ordering here
        // (to facilitate the creation of valid base triangles)
        // any two consecutive edges shall share a vertex
        jStart = 3;
        jEnd = i;
        jStep = -1;
      }

      for(SimplexId j = jStart; j != jEnd; j += jStep) {

        if(((d[i] > 0) && (d[j] < 0)) || ((d[i] < 0) && (d[j] > 0))) {

          // the edge is crossed by a base triangle
          if(triangleEdgeNumbers[0] == 3) {
            triangleEdges[1][triangleEdgeNumbers[1]] = edgeCounter;
            triangleEdgeNumbers[1]++;
          } else {
            triangleEdges[0][triangleEdgeNumbers[0]] = edgeCounter;
            triangleEdgeNumbers[0]++;
          }
        }

        if((d[i] == 0) && (d[j] == 0)) {
          // special case of a seed tet containing the jacobi edge.
          // the entire edge is on the fiber surface.
          // let's put this edge twice.
          // NOTE: in such a case, we're producing only one triangle.
          // so we just need to add that to the first triangle.
          triangleEdges[0][triangleEdgeNumbers[0]] = edgeCounter;
          triangleEdgeNumbers[0]++;

          triangleEdges[0][triangleEdgeNumbers[0]] = edgeCounter;
          triangleEdgeNumbers[0]++;
        }

        edgeCounter++;
      }
    }

    // post-process in the case of a second base triangle
    if(triangleEdges[1][0] != -1) {
      if(triangleEdges[1][1] == -1) {
        // we need to consider the edges of the first triangle which share
        // a vertex with the second triangle's edge.
        // given the ordering, the following edge pairs are forbidden:
        // (opposite edges of the tetrahedron)
        // 0/5
        // 1/3
        // 2/4
        SimplexId forbiddenEdge = -1;
        switch(triangleEdges[1][0]) {
          case 0:
            forbiddenEdge = 5;
            break;
          case 1:
            forbiddenEdge = 3;
            break;
          case 2:
            forbiddenEdge = 4;
            break;
          case 3:
            forbiddenEdge = 1;
            break;
          case 4:
            forbiddenEdge = 2;
            break;
          case 5:
            forbiddenEdge = 0;
            break;
        }
        for(SimplexId i = 0; i < (SimplexId)triangleEdges[0].size(); i++) {
          if(triangleEdges[0][i] != forbiddenEdge) {
            if(triangleEdges[1][1] != -1) {
              triangleEdges[1][2] = triangleEdges[0][i];
              break;
            } else {
              triangleEdges[1][1] = triangleEdges[0][i];
            }
          }
        }
      }
    }

    // post-process in case of exactly one vertex on the jacobi edge
    for(int i = 0; i < 2; i++) {
      if(triangleEdges[i][0] != -1) {
        // the jacobi vertex has to be the last one
        if(triangleEdges[i][2] == -1) {
          // add whatever edge connected to equalVertexLocalId
          switch(equalVertexLocalId) {
            case 0:
            case 1:
              triangleEdges[i][2] = 0;
              break;
            case 2:
            case 3:
              triangleEdges[i][2] = 5;
              break;
          }
        }
      }
    }

    // 3. crop the resulting triangles to the saddleEdge
    // for each edge recorded previously, get the u,v coordinates of the
    // intersection point.
    // take its barycentric coordinates on the saddle edge.
    // see if it lies in it (in between 0 and 1)
    // figure 7 of the paper
    double d0, d1;
    std::pair<double, double> uv0, uv1;
    std::vector<std::pair<double, double>> uv(3);
    std::vector<double> t(3);

    SimplexId createdVertices = 0;

    for(int i = 0; i < 2; i++) {
      if(triangleEdges[i][0] != -1) {
        // this is a valid triangle, let's go ahead.

        SimplexId lowerVertexNumber = 0;
        SimplexId upperVertexNumber = 0;
        SimplexId greyVertexNumber = 0;
        // iterate over the edges and compute the edge intersection (range)
        for(SimplexId j = 0; j < (SimplexId)triangleEdges[i].size(); j++) {

          SimplexId vertexId0 = 0, vertexId1 = 0;
          if(triangulation) {
            triangulation->getCellVertex(
              tetId, edgeImplicitEncoding_[2 * triangleEdges[i][j]], vertexId0);
            triangulation->getCellVertex(
              tetId, edgeImplicitEncoding_[2 * triangleEdges[i][j] + 1],
              vertexId1);
          } else {
            vertexId0
              = tetList_[5 * tetId + 1
                         + edgeImplicitEncoding_[2 * triangleEdges[i][j]]];
            vertexId1
              = tetList_[5 * tetId + 1
                         + edgeImplicitEncoding_[2 * triangleEdges[i][j] + 1]];
          }

          if((j < (SimplexId)triangleEdges[i].size() - 1)
             && (triangleEdges[i][j] == triangleEdges[i][j + 1])) {

            // special case of a jacobi edge
            uv[j].first = ((dataTypeU *)uField_)[vertexId0];
            uv[j].second = ((dataTypeV *)vField_)[vertexId0];

            uv[j + 1].first = ((dataTypeU *)uField_)[vertexId1];
            uv[j + 1].second = ((dataTypeV *)vField_)[vertexId1];

          } else if((!j)
                    || ((j)
                        && (triangleEdges[i][j] != triangleEdges[i][j - 1]))) {

            // regular intersection case
            d0 = d[edgeImplicitEncoding_[2 * triangleEdges[i][j]]];
            uv0.first = ((dataTypeU *)uField_)[vertexId0];
            uv0.second = ((dataTypeV *)vField_)[vertexId0];

            d1 = d[edgeImplicitEncoding_[2 * triangleEdges[i][j] + 1]];
            uv1.first = ((dataTypeU *)uField_)[vertexId1];
            uv1.second = ((dataTypeV *)vField_)[vertexId1];

            uv[j].first
              = uv0.first + (d0 / (d0 - d1)) * (uv1.first - uv0.first);
            uv[j].second
              = uv0.second + (d0 / (d0 - d1)) * (uv1.second - uv0.second);
          }

          // now determine the line parameterization of this intersection
          // point on the saddle edge
          if(fabs(rangePoint1.first - rangePoint0.first)
             > fabs(rangePoint1.second - rangePoint0.second)) {
            t[j] = (uv[j].first - rangePoint0.first)
                   / (rangePoint1.first - rangePoint0.first);
          } else {
            t[j] = (uv[j].second - rangePoint0.second)
                   / (rangePoint1.second - rangePoint0.second);
          }

          if((t[j] <= 1) && (t[j] >= 0))
            greyVertexNumber++;
          else if(t[j] < 0)
            lowerVertexNumber++;
          else
            upperVertexNumber++;
        }
        // at this point, we know the uv coordinates (and the edge param)
        // of each vertex of the base triangle.
        // we can proceed with the cropping

        // 4. triangulate the result
        if(greyVertexNumber == 3) {
          createdVertices += computeCase0<dataTypeU, dataTypeV>(
            polygonEdgeId, tetId, triangleEdges[i][0], t[0], uv[0].first,
            uv[0].second, triangleEdges[i][1], t[1], uv[1].first, uv[1].second,
            triangleEdges[i][2], t[2], uv[2].first, uv[2].second,
            triangulation);
        } else if(lowerVertexNumber == 3 || upperVertexNumber == 3) {
          // well do nothing (empty triangle)
        } else if((lowerVertexNumber == 1) && (upperVertexNumber == 1)
                  && (greyVertexNumber == 1)) {
          createdVertices += computeCase1<dataTypeU, dataTypeV>(
            polygonEdgeId, tetId, triangleEdges[i][0], t[0], uv[0].first,
            uv[0].second, triangleEdges[i][1], t[1], uv[1].first, uv[1].second,
            triangleEdges[i][2], t[2], uv[2].first, uv[2].second,
            triangulation);
        } else if(((lowerVertexNumber == 2) && (upperVertexNumber == 1))
                  || ((lowerVertexNumber == 1) && (upperVertexNumber == 2))) {
          createdVertices += computeCase2<dataTypeU, dataTypeV>(
            polygonEdgeId, tetId, triangleEdges[i][0], t[0], uv[0].first,
            uv[0].second, triangleEdges[i][1], t[1], uv[1].first, uv[1].second,
            triangleEdges[i][2], t[2], uv[2].first, uv[2].second,
            triangulation);
        } else if((greyVertexNumber == 1)
                  && ((lowerVertexNumber == 2) || (upperVertexNumber == 2))) {
          createdVertices += computeCase3<dataTypeU, dataTypeV>(
            polygonEdgeId, tetId, triangleEdges[i][0], t[0], uv[0].first,
            uv[0].second, triangleEdges[i][1], t[1], uv[1].first, uv[1].second,
            triangleEdges[i][2], t[2], uv[2].first, uv[2].second,
            triangulation);
        } else if(((greyVertexNumber == 2))
                  && ((lowerVertexNumber == 1) || (upperVertexNumber == 1))) {
          createdVertices += computeCase4<dataTypeU, dataTypeV>(
            polygonEdgeId, tetId, triangleEdges[i][0], t[0], uv[0].first,
            uv[0].second, triangleEdges[i][1], t[1], uv[1].first, uv[1].second,
            triangleEdges[i][2], t[2], uv[2].first, uv[2].second,
            triangulation);
        }
      }
    }

    // in case of 2 triangles, remesh locally
    // NOTE: here we just snap vertices together if they are colinear
    // the mergeFiberSurfaces() function will take care of removing any
    // duplicate
    if((triangleEdges[1][0] != -1) && (createdVertices > 3)) {

      std::vector<SimplexId> createdVertexList(createdVertices);
      for(SimplexId i = 0; i < (SimplexId)createdVertices; i++) {
        createdVertexList[i]
          = polygonEdgeVertexLists_[polygonEdgeId]->size() - 1 - i;
      }

      std::vector<bool> snappedVertices(createdVertices, false);

      for(SimplexId i = 0; i < createdVertices; i++) {

        std::vector<SimplexId> colinearVertices;
        if(!snappedVertices[i]) {
          colinearVertices.push_back(i);
          for(SimplexId j = 0; j < createdVertices; j++) {
            if((i != j) && (!snappedVertices[j])) {
              // not the same vertex
              // not snapped already

              if((*polygonEdgeVertexLists_[polygonEdgeId])[createdVertexList[i]]
                   .t_
                 == (*polygonEdgeVertexLists_[polygonEdgeId])
                      [createdVertexList[j]]
                        .t_) {
                colinearVertices.push_back(j);
              }
            }
          }
        }

        if((colinearVertices.size() == 4) || (colinearVertices.size() == 3)) {
          // 3 co-linear vertices with 1 duplicate

          // we just need to find the pair of duplicates and snap both of
          // them to another vertex
          std::pair<SimplexId, SimplexId> minPair;
          double minDistance = -1;
          for(SimplexId j = 0; j < (SimplexId)colinearVertices.size(); j++) {
            for(SimplexId k = 0; k < (SimplexId)colinearVertices.size(); k++) {
              if(j != k) {

                double distance = Geometry::distance(
                  (*polygonEdgeVertexLists_[polygonEdgeId])
                    [createdVertexList[colinearVertices[j]]]
                      .p_.data(),
                  (*polygonEdgeVertexLists_[polygonEdgeId])
                    [createdVertexList[colinearVertices[k]]]
                      .p_.data());

                //                 bool basePointSnap = true;
                //                 for(int l = 0; l < 3; l++){
                //                   if((*polygonEdgeVertexLists_[polygonEdgeId])[
                //                     createdVertexList[colinearVertices[j]]].p_[l]
                //                     !=
                //                   (*polygonEdgeVertexLists_[polygonEdgeId])[
                //                       createdVertexList[colinearVertices[k]]].p_[l]){
                //                     basePointSnap = false;
                //                     break;
                //                   }
                //                 }

                //                 if(!basePointSnap){
                if((minDistance < 0) || (distance < minDistance)) {
                  minDistance = distance;
                  minPair.first = j;
                  minPair.second = k;
                }
                //                 }
              }
            }
          }
          if((minDistance != -1) && (minDistance < prec_dbl)) {
            //           if((minDistance != -1)&&(minDistance <
            //           pointSnappingThreshold_)){
            // snap them to another colinear vertex
            for(SimplexId j = 0; j < (SimplexId)colinearVertices.size(); j++) {
              if((j != minPair.first) && (j != minPair.second)) {
                // snap minPair.first to j

                for(int k = 0; k < 3; k++) {
                  (*polygonEdgeVertexLists_[polygonEdgeId])
                    [createdVertexList[colinearVertices[minPair.first]]]
                      .p_[k]
                    = (*polygonEdgeVertexLists_[polygonEdgeId])
                        [createdVertexList[colinearVertices[j]]]
                          .p_[k];
                }
                (*polygonEdgeVertexLists_[polygonEdgeId])
                  [createdVertexList[colinearVertices[minPair.first]]]
                    .uv_
                  = (*polygonEdgeVertexLists_[polygonEdgeId])
                      [createdVertexList[colinearVertices[j]]]
                        .uv_;

                (*polygonEdgeVertexLists_[polygonEdgeId])
                  [createdVertexList[colinearVertices[minPair.first]]]
                    .t_
                  = (*polygonEdgeVertexLists_[polygonEdgeId])
                      [createdVertexList[colinearVertices[j]]]
                        .t_;

                (*polygonEdgeVertexLists_[polygonEdgeId])
                  [createdVertexList[colinearVertices[minPair.first]]]
                    .isBasePoint_
                  = (*polygonEdgeVertexLists_[polygonEdgeId])
                      [createdVertexList[colinearVertices[j]]]
                        .isBasePoint_;

                (*polygonEdgeVertexLists_[polygonEdgeId])
                  [createdVertexList[colinearVertices[minPair.first]]]
                    .isIntersectionPoint_
                  = (*polygonEdgeVertexLists_[polygonEdgeId])
                      [createdVertexList[colinearVertices[j]]]
                        .isIntersectionPoint_;
                if((*polygonEdgeVertexLists_[polygonEdgeId])
                     [createdVertexList[colinearVertices[j]]]
                       .meshEdge_.first
                   != -1) {
                  (*polygonEdgeVertexLists_[polygonEdgeId])
                    [createdVertexList[colinearVertices[minPair.first]]]
                      .meshEdge_
                    = (*polygonEdgeVertexLists_[polygonEdgeId])
                        [createdVertexList[colinearVertices[j]]]
                          .meshEdge_;
                }

                // snap minPair.second to j
                for(int k = 0; k < 3; k++) {
                  (*polygonEdgeVertexLists_[polygonEdgeId])
                    [createdVertexList[colinearVertices[minPair.second]]]
                      .p_[k]
                    = (*polygonEdgeVertexLists_[polygonEdgeId])
                        [createdVertexList[colinearVertices[j]]]
                          .p_[k];
                }
                (*polygonEdgeVertexLists_[polygonEdgeId])
                  [createdVertexList[colinearVertices[minPair.second]]]
                    .uv_
                  = (*polygonEdgeVertexLists_[polygonEdgeId])
                      [createdVertexList[colinearVertices[j]]]
                        .uv_;

                (*polygonEdgeVertexLists_[polygonEdgeId])
                  [createdVertexList[colinearVertices[minPair.second]]]
                    .t_
                  = (*polygonEdgeVertexLists_[polygonEdgeId])
                      [createdVertexList[colinearVertices[j]]]
                        .t_;

                (*polygonEdgeVertexLists_[polygonEdgeId])
                  [createdVertexList[colinearVertices[minPair.second]]]
                    .isBasePoint_
                  = (*polygonEdgeVertexLists_[polygonEdgeId])
                      [createdVertexList[colinearVertices[j]]]
                        .isBasePoint_;

                (*polygonEdgeVertexLists_[polygonEdgeId])
                  [createdVertexList[colinearVertices[minPair.second]]]
                    .isIntersectionPoint_
                  = (*polygonEdgeVertexLists_[polygonEdgeId])
                      [createdVertexList[colinearVertices[j]]]
                        .isIntersectionPoint_;
                if((*polygonEdgeVertexLists_[polygonEdgeId])
                     [createdVertexList[colinearVertices[j]]]
                       .meshEdge_.first
                   != -1) {
                  (*polygonEdgeVertexLists_[polygonEdgeId])
                    [createdVertexList[colinearVertices[minPair.second]]]
                      .meshEdge_
                    = (*polygonEdgeVertexLists_[polygonEdgeId])
                        [createdVertexList[colinearVertices[j]]]
                          .meshEdge_;
                }

                snappedVertices[colinearVertices[minPair.first]] = true;
                snappedVertices[colinearVertices[minPair.second]] = true;

                // we're done
                break;
              }
            }
          }
        }
      }
    }

    return createdVertices;
  }

  return 0;
}

template <class dataTypeU, class dataTypeV>
inline int ttk::FiberSurface::remeshIntersections() const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(!tetNumber_)
    return -1;
#endif

  // Algorithm
  // 1) loop over each triangle and mark the containing tetrahedron
  // 2) for each tet, in parallel, check for pairwise triangle intersections
  // (based on the triangles' range projection)
  // Given a pair of intersecting triangles:
  // 3) given the point of intersection, take one of its coordinates (u or v)
  // and compute an iso-contour I on both triangles
  // 4) express the barycentric coordinates of I in both triangles, 2 cases:
  // a) the intersection is the entire fiber (old code)
  // b) the intersection is a segment of fiber

  // NOTE:
  // the topological aspect of the code is OK.
  // if any bug, it's very likely to be geometry accuracy related.

  std::vector<std::vector<IntersectionTriangle>> tetIntersections(tetNumber_);
  std::vector<std::vector<Vertex>> tetNewVertices(tetNumber_);

  // fill the information prior to the parallel pass
  for(SimplexId i = 0; i < (SimplexId)polygonEdgeTriangleLists_.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)polygonEdgeTriangleLists_[i]->size();
        j++) {

      SimplexId tetId = (*polygonEdgeTriangleLists_[i])[j].tetId_;

      tetIntersections[tetId].resize(tetIntersections[tetId].size() + 1);

      tetIntersections[tetId].back().caseId_
        = (*polygonEdgeTriangleLists_[i])[j].caseId_;
      tetIntersections[tetId].back().polygonEdgeId_ = i;
      tetIntersections[tetId].back().triangleId_ = j;
      for(int k = 0; k < 3; k++) {
        tetIntersections[tetId].back().vertexIds_[k]
          = (*polygonEdgeTriangleLists_[i])[j].vertexIds_[k];
        tetIntersections[tetId].back().uv_[k]
          = (*globalVertexList_)[tetIntersections[tetId].back().vertexIds_[k]]
              .uv_;
        tetIntersections[tetId].back().t_[k]
          = (*globalVertexList_)[tetIntersections[tetId].back().vertexIds_[k]]
              .t_;
        for(int l = 0; l < 3; l++) {
          tetIntersections[tetId].back().p_[k][l]
            = (*globalVertexList_)[tetIntersections[tetId].back().vertexIds_[k]]
                .p_[l];
        }
      }
      tetIntersections[tetId].back().intersection_.first = -DBL_MAX;
      tetIntersections[tetId].back().intersection_.second = -DBL_MAX;
    }
  }

  std::vector<SimplexId> tetList;
  for(SimplexId i = 0; i < (SimplexId)tetIntersections.size(); i++) {
    if(tetIntersections[i].size() > 1)
      tetList.push_back(i);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)tetList.size(); i++) {
    SimplexId tetId = tetList[i];

    // pre-process by merging nearby vertices...?

    SimplexId newTriangleNumber = 1;
    SimplexId newVertexNumber = 1;

    for(SimplexId j = 0; j < (SimplexId)tetIntersections[tetId].size(); j++) {

      if(j > 1000) {
        this->printWrn("Preventing an infinite loop!");
        this->printWrn("More than 1000 re-meshed triangles in tet #"
                       + std::to_string(tetId));
        this->printWrn("Extra-thin triangles keep on intersecting?!");
        break;
      }

      // if we re-mesh, we add new triangles at the end of the list.
      // there's no need to check intersections with those.
      SimplexId originalTriangleNumber
        = (SimplexId)tetIntersections[tetId].size();

      for(SimplexId k = 0; k < originalTriangleNumber; k++) {

        SimplexId polygonEdgeId0 = tetIntersections[tetId][j].polygonEdgeId_;
        SimplexId polygonEdgeId1 = tetIntersections[tetId][k].polygonEdgeId_;

        if((j != k) && (polygonEdgeId0 != polygonEdgeId1)) {
          // cases 3, 4 and 6 of the fiber surface table (multiple triangle
          // per tet given a single edge). we don't need to re-mesh that.

          // grab the range projection for the triangle j
          std::pair<double, double> edge0point0, edge0point1;

          getTriangleRangeExtremities(
            tetId, j, tetIntersections, edge0point0, edge0point1);

          // now do the same thing for the triangle k
          std::pair<double, double> edge1point0, edge1point1;

          getTriangleRangeExtremities(
            tetId, k, tetIntersections, edge1point0, edge1point1);

          // compute the intersection
          std::pair<double, double> intersection;
          bool hasIntersection = Geometry::computeSegmentIntersection(
            edge0point0.first, edge0point0.second, edge0point1.first,
            edge0point1.second, edge1point0.first, edge1point0.second,
            edge1point1.first, edge1point1.second, intersection.first,
            intersection.second);

          if((hasIntersection)
            // check if that intersection has been registered before
            // in the end, only one intersection per triangle, no matter what
            /*&&(((fabs(tetIntersections[tetId][j].intersection_.first
              - intersection.first) > Geometry::powIntTen(-FLT_DIG))
            ||(fabs(tetIntersections[tetId][j].intersection_.second
              - intersection.second) > Geometry::powIntTen(-FLT_DIG)))
            &&((fabs(tetIntersections[tetId][k].intersection_.first
              - intersection.first) > Geometry::powIntTen(-FLT_DIG))
            ||(fabs(tetIntersections[tetId][k].intersection_.second
              - intersection.second) > Geometry::powIntTen(-FLT_DIG))))*/){

            computeTriangleIntersection(tetId, j, k, polygonEdgeId0,
                                        polygonEdgeId1, intersection,
                                        newVertexNumber, newTriangleNumber,
                                        tetIntersections, tetNewVertices);
          }
        }
      }
    }
  }

  // now copy the new vertices
  for(SimplexId i = 0; i < (SimplexId)tetNewVertices.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)tetNewVertices[i].size(); j++) {

      SimplexId localId = (*globalVertexList_).size();
      tetNewVertices[i][j].localId_ = localId;
      (*globalVertexList_).push_back(tetNewVertices[i][j]);
    }
  }

  for(SimplexId i = 0; i < (SimplexId)tetIntersections.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)tetIntersections[i].size(); j++) {

      if(((tetIntersections[i][j].intersection_.first != -DBL_MAX)
          && (tetIntersections[i][j].intersection_.second != -DBL_MAX))
         || (tetIntersections[i][j].triangleId_ < 0)) {

        SimplexId triangleId = tetIntersections[i][j].triangleId_;

        if(triangleId < 0) {

          // this is a new triangle
          triangleId = (*polygonEdgeTriangleLists_[tetIntersections[i][j]
                                                     .polygonEdgeId_])
                         .size();
          (*polygonEdgeTriangleLists_[tetIntersections[i][j].polygonEdgeId_])
            .resize(triangleId + 1);
          (*polygonEdgeTriangleLists_[tetIntersections[i][j].polygonEdgeId_])
            .back()
            .tetId_
            = i;
          (*polygonEdgeTriangleLists_[tetIntersections[i][j].polygonEdgeId_])
            .back()
            .caseId_
            = tetIntersections[i][j].caseId_;
        }

        for(int k = 0; k < 3; k++) {

          SimplexId vertexId = tetIntersections[i][j].vertexIds_[k];

          if(vertexId < 0) {
            // newly created vertex
            vertexId = tetNewVertices[i][-(vertexId + 1)].localId_;
          }
          (*polygonEdgeTriangleLists_[tetIntersections[i][j]
                                        .polygonEdgeId_])[triangleId]
            .vertexIds_[k]
            = vertexId;
        }
      }
    }
  }

  return 0;
}
