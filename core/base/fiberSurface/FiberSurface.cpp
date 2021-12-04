#include <FiberSurface.h>

using namespace std;
using namespace ttk;

static const float PREC_FLT{powf(10.F, -FLT_DIG)};
static const float PREC_FLT_2{powf(10.F, -FLT_DIG + 2)};
static const double PREC_DBL{Geometry::pow(10.0, -DBL_DIG)};
static const double PREC_DBL_4{Geometry::pow(10.0, -DBL_DIG + 4)};

FiberSurface::FiberSurface() {
  this->setDebugMsgPrefix("FiberSurface");
}

int FiberSurface::getNumberOfCommonVertices(
  const SimplexId &tetId,
  const SimplexId &triangleId0,
  const SimplexId &triangleId1,
  const vector<vector<IntersectionTriangle>> &tetIntersections) const {

  // if the two triangles have at least 2 vertices in common, they are adjacent.
  SimplexId commonVertexNumber = 0;

  for(int i = 0; i < 3; i++) {
    vector<double> p0(3);

    for(int j = 0; j < 3; j++) {
      p0[j] = tetIntersections[tetId][triangleId0].p_[i][j];
    }

    // check if this guy exists in the other triangle
    for(int j = 0; j < 3; j++) {
      vector<double> p1(3);

      bool isTheSame = true;
      for(int k = 0; k < 3; k++) {
        p1[k] = tetIntersections[tetId][triangleId1].p_[j][k];

        if(fabs(p0[k] - p1[k]) > PREC_FLT) {
          isTheSame = false;
          break;
        }
      }
      if(isTheSame) {
        commonVertexNumber++;
        break;
      }
    }
  }

  return commonVertexNumber;
}

int FiberSurface::computeTriangleFiber(
  const SimplexId &tetId,
  const SimplexId &triangleId,
  const pair<double, double> &intersection,
  const vector<vector<IntersectionTriangle>> &tetIntersections,
  vector<double> &pA,
  vector<double> &pB,
  SimplexId &pivotVertexId,
  bool &edgeFiber) const {

  pivotVertexId = -1;

  // let's check first if an edge coincide with the fiber
  for(int i = 0; i < 3; i++) {
    if((fabs(intersection.first
             - tetIntersections[tetId][triangleId].uv_[i].first)
        < PREC_DBL_4)
       && fabs(intersection.second
               - tetIntersections[tetId][triangleId].uv_[i].second)
            < PREC_DBL_4
       && fabs(intersection.first
               - tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].first)
            < PREC_DBL_4
       && fabs(intersection.second
               - tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].second)
            < PREC_DBL_4) {
      // edge 0 - 1 is on the fiber. the pivot is 2
      pivotVertexId = (i + 2) % 3;
      edgeFiber = true;
      break;
    }
  }

  // NOTE: it'd be better here to compute the actual distance
  if(pivotVertexId == -1) {
    for(int i = 0; i < 3; i++) {
      if(tetIntersections[tetId][triangleId].uv_[i].first
         > intersection.first) {
        if((tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].first
            <= intersection.first)
           && (tetIntersections[tetId][triangleId].uv_[(i + 2) % 3].first
               <= intersection.first)) {
          pivotVertexId = i;
          break;
        }
      }
      if(tetIntersections[tetId][triangleId].uv_[i].first
         < intersection.first) {
        if((tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].first
            >= intersection.first)
           && (tetIntersections[tetId][triangleId].uv_[(i + 2) % 3].first
               >= intersection.first)) {
          pivotVertexId = i;
          break;
        }
      }
    }
  }
  if(pivotVertexId == -1) {
    for(int i = 0; i < 3; i++) {
      if(tetIntersections[tetId][triangleId].uv_[i].second
         > intersection.second) {
        if((tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].second
            <= intersection.second)
           && (tetIntersections[tetId][triangleId].uv_[(i + 2) % 3].second
               <= intersection.second)) {
          pivotVertexId = i;
          break;
        }
      }
      if(tetIntersections[tetId][triangleId].uv_[i].second
         < intersection.second) {
        if((tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].second
            >= intersection.second)
           && (tetIntersections[tetId][triangleId].uv_[(i + 2) % 3].second
               >= intersection.second)) {
          pivotVertexId = i;
          break;
        }
      }
    }
  }

  if(pivotVertexId == -1) {
    // initially, several triangles can intersect.
    // after a few iterations, a valid triangle for re-meshing may be
    // subdivided into several sub-triangles. not all of them may
    // intersect with a given query (pivotVertexId == -1)
    return -1;
  }

  // compute the interpolations
  vector<double> baryCentrics0, baryCentrics1;
  vector<double> p(2), p0(2), p1(2), p2(2);

  p[0] = intersection.first;
  p[1] = intersection.second;

  p0[0] = tetIntersections[tetId][triangleId].uv_[pivotVertexId].first;
  p0[1] = tetIntersections[tetId][triangleId].uv_[pivotVertexId].second;

  p1[0]
    = tetIntersections[tetId][triangleId].uv_[(pivotVertexId + 1) % 3].first;
  p1[1]
    = tetIntersections[tetId][triangleId].uv_[(pivotVertexId + 1) % 3].second;

  p2[0]
    = tetIntersections[tetId][triangleId].uv_[(pivotVertexId + 2) % 3].first;
  p2[1]
    = tetIntersections[tetId][triangleId].uv_[(pivotVertexId + 2) % 3].second;

  Geometry::computeBarycentricCoordinates(
    p0.data(), p1.data(), p.data(), baryCentrics0, 2);
  Geometry::computeBarycentricCoordinates(
    p0.data(), p2.data(), p.data(), baryCentrics1, 2);

  pA.resize(3);
  for(int i = 0; i < 3; i++) {
    pA[i] = baryCentrics0[0]
              * tetIntersections[tetId][triangleId].p_[pivotVertexId][i]
            + baryCentrics0[1]
                * tetIntersections[tetId][triangleId]
                    .p_[(pivotVertexId + 1) % 3][i];
  }

  pB.resize(3);
  for(int i = 0; i < 3; i++) {
    pB[i] = baryCentrics1[0]
              * tetIntersections[tetId][triangleId].p_[pivotVertexId][i]
            + baryCentrics1[1]
                * tetIntersections[tetId][triangleId]
                    .p_[(pivotVertexId + 2) % 3][i];
  }

  return 0;
}

int FiberSurface::computeTriangleIntersection(
  const SimplexId &tetId,
  const SimplexId &triangleId0,
  const SimplexId &triangleId1,
  const SimplexId &polygonEdgeId0,
  const SimplexId &polygonEdgeId1,
  const std::pair<double, double> &intersection,
  SimplexId &newVertexNumber,
  SimplexId &newTriangleNumber,
  std::vector<std::vector<IntersectionTriangle>> &tetIntersections,
  std::vector<std::vector<Vertex>> &tetNewVertices) const {

  SimplexId commonVertexNumber = getNumberOfCommonVertices(
    tetId, triangleId0, triangleId1, tetIntersections);

  // make sure the two triangles are not already adjacent
  if(commonVertexNumber == 2) {
    // NOTE: here, we used to quit with only one.
    // however, in high res data-sets you can have triangles that intersect
    // and share one vertex.
    return -1;
  }

  SimplexId pivotVertexIda = -1, pivotVertexIdb = -1;
  vector<double> p0a, p1a, p0b, p1b;

  // extract the fiber in both triangles and see if they match up
  bool edgeFiber0 = false;
  computeTriangleFiber(tetId, triangleId0, intersection, tetIntersections, p0a,
                       p1a, pivotVertexIda, edgeFiber0);

  bool edgeFiber1 = false;
  computeTriangleFiber(tetId, triangleId1, intersection, tetIntersections, p0b,
                       p1b, pivotVertexIdb, edgeFiber1);

  //   if((commonVertexNumber == 1)&&
  //     ((edgeFiber0)||(edgeFiber1))){
  //     // case of adjacent triangles along a vertex.
  //     // one of the two triangles has a fiber along an edge.
  //     return -2;
  //   }

  if(p0a.size() != 3)
    return -1;
  if(p1a.size() != 3)
    return -2;
  if(p0b.size() != 3)
    return -3;
  if(p1b.size() != 3)
    return -4;

  // we need to make sure p0a and p1a are not the same (vertex case)
  bool vertexA = false;
  bool vertexB = false;
  if((fabs(p0a[0] - p1a[0]) < PREC_DBL) && (fabs(p0a[1] - p1a[1]) < PREC_DBL)
     && (fabs(p0a[2] - p1a[2]) < PREC_DBL)) {
    vertexA = true;
  }
  if((fabs(p0b[0] - p1b[0]) < PREC_DBL) && (fabs(p0b[1] - p1b[1]) < PREC_DBL)
     && (fabs(p0b[2] - p1b[2]) < PREC_DBL)) {
    vertexB = true;
  }
  if((vertexA) || (vertexB)) {
    // intersection on triangle vertices
    return -1;
  }

  bool foundA = false, foundB = false;
  vector<double> pA, pB;

  // test if p0a lies in [p0b, p1b]
  if(Geometry::isPointOnSegment(p0a.data(), p0b.data(), p1b.data())) {
    // p0a is in the segment
    pA = p0a;
    foundA = true;
  }

  // test if p1a lies in [p0b, p1b]
  if(Geometry::isPointOnSegment(p1a.data(), p0b.data(), p1b.data())) {
    if(!foundA) {
      pA = p1a;
      foundA = true;
    } else if(!foundB) {
      // check it's far enough from pA
      if((fabs(pA[0] - p1a[0]) > PREC_DBL_4)
         || (fabs(pA[1] - p1a[1]) > PREC_DBL_4)
         || (fabs(pA[2] - p1a[2]) > PREC_DBL_4)) {
        pB = p1a;
        foundB = true;
      }
    }
  }

  // test if p0b lies in [p0a, p1a]
  if(Geometry::isPointOnSegment(p0b.data(), p0a.data(), p1a.data())) {
    if(!foundA) {
      pA = p0b;
      foundA = true;
    } else if(!foundB) {
      // check it's far enough from pA
      if((fabs(pA[0] - p0b[0]) > PREC_DBL_4)
         || (fabs(pA[1] - p0b[1]) > PREC_DBL_4)
         || (fabs(pA[2] - p0b[2]) > PREC_DBL_4)) {
        pB = p0b;
        foundB = true;
      }
    }
  }

  // test if p1b lies in [p0a, p1a]
  if(Geometry::isPointOnSegment(p1b.data(), p0a.data(), p1a.data())) {
    if(!foundA) {
      pA = p1b;
    } else if(!foundB) {
      // check it's far enough from pA
      if((fabs(pA[0] - p1b[0]) > PREC_DBL_4)
         || (fabs(pA[1] - p1b[1]) > PREC_DBL_4)
         || (fabs(pA[2] - p1b[2]) > PREC_DBL_4)) {
        pB = p1b;
      }
    }
  }

  if((!pA.size()) || (!pB.size()))
    return -2;

  if(!edgeFiber0) {
    computeTriangleIntersection(
      tetId, triangleId0, polygonEdgeId0, intersection, pA, pB, pivotVertexIda,
      newVertexNumber, newTriangleNumber, tetIntersections, tetNewVertices);
  }

  if(!edgeFiber1) {
    computeTriangleIntersection(
      tetId, triangleId1, polygonEdgeId1, intersection, pA, pB, pivotVertexIdb,
      newVertexNumber, newTriangleNumber, tetIntersections, tetNewVertices);
  }

  return 0;
}

int FiberSurface::computeTriangleIntersection(
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
  std::vector<std::vector<Vertex>> &tetNewVertices) const {

  // check if the triangle has already been intersected on that fiber
  if((fabs(tetIntersections[tetId][triangleId].intersection_.first
           - intersection.first)
      < PREC_FLT)
     && (fabs(tetIntersections[tetId][triangleId].intersection_.second
              - intersection.second)
         < PREC_FLT)) {

    return -2;
  }

  // check if the fiber is on an edge or not (if so we stop)
  for(int i = 0; i < 3; i++) {
    if((fabs(intersection.first
             - tetIntersections[tetId][triangleId].uv_[i].first)
        < PREC_FLT)
       && fabs(intersection.second
               - tetIntersections[tetId][triangleId].uv_[i].second)
            < PREC_FLT
       && fabs(intersection.first
               - tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].first)
            < PREC_FLT
       && fabs(intersection.second
               - tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].second)
            < PREC_FLT) {
      return -3;
    }
  }

  // 1. compute the barycentric coordinates of pA and pB
  vector<double> barypA, barypB;
  Geometry::computeBarycentricCoordinates(
    tetIntersections[tetId][triangleId].p_[0].data(),
    tetIntersections[tetId][triangleId].p_[1].data(),
    tetIntersections[tetId][triangleId].p_[2].data(), pA.data(), barypA);

  Geometry::computeBarycentricCoordinates(
    tetIntersections[tetId][triangleId].p_[0].data(),
    tetIntersections[tetId][triangleId].p_[1].data(),
    tetIntersections[tetId][triangleId].p_[2].data(), pB.data(), barypB);

  // 2. between the two, find the closest point from the edge
  // [pivotVertexId, (pivotVertexId+2)%3]
  // that's the vertex which minimizes its coordinate [(pivotVertexId+1)%3]
  vector<double> A = pA, B = pB;
  vector<double> baryA = barypA, baryB = barypB;
  if(fabs(barypB[(pivotVertexId + 1) % 3])
     < fabs(barypA[(pivotVertexId + 1) % 3])) {
    // let's swith the two
    A = pB;
    B = pA;
    baryA = barypB;
    baryB = barypA;
  }

  bool isAVertex = false;
  for(int i = 0; i < 3; i++) {
    if(fabs(baryA[i] - 1) < PREC_DBL) {
      isAVertex = true;
      break;
    }
  }
  bool isBVertex = false;
  for(int i = 0; i < 3; i++) {
    if(fabs(baryB[i] - 1) < PREC_DBL) {
      isBVertex = true;
      break;
    }
  }
  if((isAVertex) && (isBVertex))
    return -4;

  // 3. create the two new vertices A and B
  SimplexId vertexIdA = newVertexNumber;
  newVertexNumber++;
  tetNewVertices[tetId].resize(tetNewVertices[tetId].size() + 1);
  for(int i = 0; i < 3; i++)
    tetNewVertices[tetId].back().p_[i] = A[i];
  tetNewVertices[tetId].back().polygonEdgeId_ = polygonEdgeId;
  tetNewVertices[tetId].back().uv_.first
    = baryA[0] * tetIntersections[tetId][triangleId].uv_[0].first
      + baryA[1] * tetIntersections[tetId][triangleId].uv_[1].first
      + baryA[2] * tetIntersections[tetId][triangleId].uv_[2].first;
  tetNewVertices[tetId].back().uv_.second
    = baryA[0] * tetIntersections[tetId][triangleId].uv_[0].second
      + baryA[1] * tetIntersections[tetId][triangleId].uv_[1].second
      + baryA[2] * tetIntersections[tetId][triangleId].uv_[2].second;
  tetNewVertices[tetId].back().t_
    = baryA[0] * tetIntersections[tetId][triangleId].t_[0]
      + baryA[1] * tetIntersections[tetId][triangleId].t_[1]
      + baryA[2] * tetIntersections[tetId][triangleId].t_[2];
  tetNewVertices[tetId].back().isIntersectionPoint_ = true;

  SimplexId vertexIdB = newVertexNumber;
  newVertexNumber++;
  tetNewVertices[tetId].resize(tetNewVertices[tetId].size() + 1);
  for(int i = 0; i < 3; i++)
    tetNewVertices[tetId].back().p_[i] = B[i];
  tetNewVertices[tetId].back().polygonEdgeId_ = polygonEdgeId;
  tetNewVertices[tetId].back().uv_.first
    = baryB[0] * tetIntersections[tetId][triangleId].uv_[0].first
      + baryB[1] * tetIntersections[tetId][triangleId].uv_[1].first
      + baryB[2] * tetIntersections[tetId][triangleId].uv_[2].first;
  tetNewVertices[tetId].back().uv_.second
    = baryB[0] * tetIntersections[tetId][triangleId].uv_[0].second
      + baryB[1] * tetIntersections[tetId][triangleId].uv_[1].second
      + baryB[2] * tetIntersections[tetId][triangleId].uv_[2].second;
  tetNewVertices[tetId].back().t_
    = baryB[0] * tetIntersections[tetId][triangleId].t_[0]
      + baryB[1] * tetIntersections[tetId][triangleId].t_[1]
      + baryB[2] * tetIntersections[tetId][triangleId].t_[2];
  tetNewVertices[tetId].back().isIntersectionPoint_ = true;

  // 4. create the triangle (without saving intersection):
  // pivotVertexId, A, (pivotVertexId+2)%3
  createNewIntersectionTriangle(tetId, triangleId, pivotVertexId, -vertexIdA,
                                (pivotVertexId + 2) % 3, tetNewVertices,
                                newTriangleNumber, tetIntersections);

  // 5. create the triangle (saving intersection)
  // A, B, pivotVertexId+2
  // special case where a, b and p+2 are actually aligned
  // we should detect a colinear triangle and test the opposite diagonal instead
  // a, b, (pivotVertexId+1)%3
  SimplexId ret = createNewIntersectionTriangle(
    tetId, triangleId, -vertexIdA, -vertexIdB, (pivotVertexId + 2) % 3,
    tetNewVertices, newTriangleNumber, tetIntersections, &intersection);
  if(ret == -1) {
    // flip the diagonal
    createNewIntersectionTriangle(
      tetId, triangleId, -vertexIdA, -vertexIdB, (pivotVertexId + 1) % 3,
      tetNewVertices, newTriangleNumber, tetIntersections, &intersection);
    createNewIntersectionTriangle(tetId, triangleId, (pivotVertexId + 1) % 3,
                                  (pivotVertexId + 2) % 3, -vertexIdA,
                                  tetNewVertices, newTriangleNumber,
                                  tetIntersections, &intersection);
  } else {
    // 6. create the triangle (saving intersection)
    // pivotVertexId+1, pivotVertexId+2, B
    createNewIntersectionTriangle(tetId, triangleId, (pivotVertexId + 1) % 3,
                                  (pivotVertexId + 2) % 3, -vertexIdB,
                                  tetNewVertices, newTriangleNumber,
                                  tetIntersections, &intersection);
  }

  // 7. create the triangle (without saving intersection)
  // pivotVertexId, B, pivotVertexId+1
  createNewIntersectionTriangle(tetId, triangleId, pivotVertexId, -vertexIdB,
                                (pivotVertexId + 1) % 3, tetNewVertices,
                                newTriangleNumber, tetIntersections);

  // 8. edit the original triangle to
  // p, A, B (saving intersection)
  tetIntersections[tetId][triangleId].intersection_ = intersection;
  tetIntersections[tetId][triangleId].vertexIds_[0]
    = tetIntersections[tetId][triangleId].vertexIds_[pivotVertexId];
  tetIntersections[tetId][triangleId].uv_[0]
    = tetIntersections[tetId][triangleId].uv_[pivotVertexId];
  tetIntersections[tetId][triangleId].t_[0]
    = tetIntersections[tetId][triangleId].t_[pivotVertexId];
  for(int i = 0; i < 3; i++)
    tetIntersections[tetId][triangleId].p_[0][i]
      = tetIntersections[tetId][triangleId].p_[pivotVertexId][i];
  // vertexA
  tetIntersections[tetId][triangleId].vertexIds_[1] = -vertexIdA;
  tetIntersections[tetId][triangleId].uv_[1]
    = tetNewVertices[tetId][vertexIdA - 1].uv_;
  tetIntersections[tetId][triangleId].t_[1]
    = tetNewVertices[tetId][vertexIdA - 1].t_;
  for(int i = 0; i < 3; i++)
    tetIntersections[tetId][triangleId].p_[1][i]
      = tetNewVertices[tetId][vertexIdA - 1].p_[i];
  // vertexB
  tetIntersections[tetId][triangleId].vertexIds_[2] = -vertexIdB;
  tetIntersections[tetId][triangleId].uv_[2]
    = tetNewVertices[tetId][vertexIdB - 1].uv_;
  tetIntersections[tetId][triangleId].t_[2]
    = tetNewVertices[tetId][vertexIdB - 1].t_;
  for(int i = 0; i < 3; i++)
    tetIntersections[tetId][triangleId].p_[2][i]
      = tetNewVertices[tetId][vertexIdB - 1].p_[i];

  return 0;
}

int FiberSurface::createNewIntersectionTriangle(
  const SimplexId &tetId,
  const SimplexId &triangleId,
  const SimplexId &vertexId0,
  const SimplexId &vertexId1,
  const SimplexId &vertexId2,
  const vector<vector<Vertex>> &tetNewVertices,
  SimplexId &newTriangleNumber,
  vector<vector<IntersectionTriangle>> &tetIntersections,
  const pair<double, double> *intersection) const {

  if(isIntersectionTriangleColinear(tetId, triangleId, tetIntersections,
                                    tetNewVertices, vertexId0, vertexId1,
                                    vertexId2)) {

    return -1;
  }

  tetIntersections[tetId].resize(tetIntersections[tetId].size() + 1);
  tetIntersections[tetId].back().caseId_
    = tetIntersections[tetId][triangleId].caseId_;
  tetIntersections[tetId].back().polygonEdgeId_
    = tetIntersections[tetId][triangleId].polygonEdgeId_;
  tetIntersections[tetId].back().triangleId_ = -newTriangleNumber;

  if(intersection) {
    tetIntersections[tetId].back().intersection_ = *intersection;
  } else {
    tetIntersections[tetId].back().intersection_.first = -DBL_MAX;
    tetIntersections[tetId].back().intersection_.second = -DBL_MAX;
  }

  newTriangleNumber++;

  // process the vertices
  for(int i = 0; i < 3; i++) {

    SimplexId vertexId = vertexId0;
    if(i == 1)
      vertexId = vertexId1;
    if(i == 2)
      vertexId = vertexId2;

    if(vertexId >= 0) {
      // this is not a newly created vertex
      tetIntersections[tetId].back().vertexIds_[i]
        = tetIntersections[tetId][triangleId].vertexIds_[vertexId];
      tetIntersections[tetId].back().uv_[i]
        = tetIntersections[tetId][triangleId].uv_[vertexId];
      tetIntersections[tetId].back().t_[i]
        = tetIntersections[tetId][triangleId].t_[vertexId];
      for(int j = 0; j < 3; j++) {
        tetIntersections[tetId].back().p_[i][j]
          = tetIntersections[tetId][triangleId].p_[vertexId][j];
      }
    } else {
      // this is a newly created vertex
      tetIntersections[tetId].back().vertexIds_[i] = vertexId;
      tetIntersections[tetId].back().uv_[i]
        = tetNewVertices[tetId][(-vertexId) - 1].uv_;
      tetIntersections[tetId].back().t_[i]
        = tetNewVertices[tetId][(-vertexId) - 1].t_;
      for(int j = 0; j < 3; j++) {
        tetIntersections[tetId].back().p_[i][j]
          = tetNewVertices[tetId][(-vertexId) - 1].p_[j];
      }
    }
  }

  return 0;
}

int FiberSurface::flipEdges() const {

  Timer t;

  vector<vector<pair<SimplexId, SimplexId>>> tetTriangles(tetNumber_);
  vector<bool> inQueue(tetNumber_, false);
  vector<SimplexId> tetList;

  for(SimplexId i = 0; i < (SimplexId)polygonEdgeTriangleLists_.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)polygonEdgeTriangleLists_[i]->size();
        j++) {

      SimplexId tetId = (*polygonEdgeTriangleLists_[i])[j].tetId_;

      if(!inQueue[tetId]) {
        inQueue[tetId] = true;
        tetList.push_back(tetId);
      }

      tetTriangles[tetId].push_back(pair<SimplexId, SimplexId>(i, j));
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)tetList.size(); i++) {

    SimplexId tetId = tetList[i];

    if(tetTriangles[tetId].size() >= 2) {

      flipEdges(tetTriangles[tetId]);
    }
  }

  this->printMsg(
    "Performed edge flips", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

int FiberSurface::flipEdges(
  std::vector<std::pair<SimplexId, SimplexId>> &triangles) const {

  for(SimplexId it = 0; it < (SimplexId)triangles.size(); it++) {

    if(it == 2)
      break;

    // avoid infinite loops
    bool hasFlipped = false;

    // go in a greedy manner
    // for a triangle, grab its smallest angle alpha.
    // look at all its neighbors.
    // do the edge flip virtually and evaluate the smallest angle beta.
    // if beta is bigger than alpha do the edge flip

    // heuristic: sort the triangles in order of their minimum angle
    // (make sure we start with bad angles)
    vector<pair<double, pair<SimplexId, SimplexId>>> localTriangles;
    for(SimplexId i = 0; i < (SimplexId)triangles.size(); i++) {
      vector<SimplexId> vertexIds(3);

      vertexIds[0]
        = (*polygonEdgeTriangleLists_[triangles[i].first])[triangles[i].second]
            .vertexIds_[0];
      vertexIds[1]
        = (*polygonEdgeTriangleLists_[triangles[i].first])[triangles[i].second]
            .vertexIds_[1];
      vertexIds[2]
        = (*polygonEdgeTriangleLists_[triangles[i].first])[triangles[i].second]
            .vertexIds_[2];

      vector<double> angles;
      Geometry::computeTriangleAngles(
        (*globalVertexList_)[vertexIds[0]].p_.data(),
        (*globalVertexList_)[vertexIds[1]].p_.data(),
        (*globalVertexList_)[vertexIds[2]].p_.data(), angles);

      double alpha = -1;
      for(int j = 0; j < 3; j++) {
        if((alpha < 0) || (fabs(angles[j]) < alpha))
          alpha = fabs(angles[j]);
      }

      localTriangles.push_back(
        pair<double, pair<SimplexId, SimplexId>>(alpha, triangles[i]));
    }

    const auto FiberSurfaceTriangleCmp
      = [](const pair<double, pair<SimplexId, SimplexId>> &t0,
           const pair<double, pair<SimplexId, SimplexId>> &t1) {
          return t0.first < t1.first;
        };

    sort(localTriangles.begin(), localTriangles.end(), FiberSurfaceTriangleCmp);

    for(SimplexId i = 0; i < (SimplexId)triangles.size(); i++) {
      triangles[i] = localTriangles[i].second;
    }

    for(SimplexId i = 0; i < (SimplexId)triangles.size(); i++) {

      vector<SimplexId> vertexIds(3);

      vertexIds[0]
        = (*polygonEdgeTriangleLists_[triangles[i].first])[triangles[i].second]
            .vertexIds_[0];
      vertexIds[1]
        = (*polygonEdgeTriangleLists_[triangles[i].first])[triangles[i].second]
            .vertexIds_[1];
      vertexIds[2]
        = (*polygonEdgeTriangleLists_[triangles[i].first])[triangles[i].second]
            .vertexIds_[2];

      vector<double> angles;
      Geometry::computeTriangleAngles(
        (*globalVertexList_)[vertexIds[0]].p_.data(),
        (*globalVertexList_)[vertexIds[1]].p_.data(),
        (*globalVertexList_)[vertexIds[2]].p_.data(), angles);

      double alpha = -1;
      for(int j = 0; j < 3; j++) {
        if((alpha < 0) || (fabs(angles[j]) < alpha))
          alpha = fabs(angles[j]);
      }

      // look at neighbors
      for(SimplexId j = 0; j < (SimplexId)triangles.size(); j++) {
        if((i != j)
           // same polygon edge
           && (((*polygonEdgeTriangleLists_[triangles[j].first])[triangles[j]
                                                                   .second]
                  .polygonEdgeId_
                == (*polygonEdgeTriangleLists_[triangles[i].first])[triangles[i]
                                                                      .second]
                     .polygonEdgeId_))) {

          SimplexId commonVertexId0 = -1;
          SimplexId commonVertexId1 = -1;
          SimplexId otherNonCommonVertexId = -1;

          for(int k = 0; k < 3; k++) {
            SimplexId vertexId
              = (*polygonEdgeTriangleLists_[triangles[j].first])[triangles[j]
                                                                   .second]
                  .vertexIds_[k];

            bool isCommon = false;

            for(int l = 0; l < 3; l++) {
              if(vertexId == vertexIds[l]) {
                if(commonVertexId0 == -1) {
                  commonVertexId0 = vertexId;
                } else {
                  commonVertexId1 = vertexId;
                }
                isCommon = true;
              }
            }
            if(!isCommon)
              otherNonCommonVertexId = vertexId;
          }

          if((commonVertexId0 != -1) && (commonVertexId1 != -1)) {

            if((!(*globalVertexList_)[commonVertexId0].isIntersectionPoint_)
               && (!(*globalVertexList_)[commonVertexId1]
                      .isIntersectionPoint_)) {

              SimplexId nonCommonVertexId = -1;
              for(int k = 0; k < 3; k++) {
                if((vertexIds[k] != commonVertexId0)
                   && (vertexIds[k] != commonVertexId1)) {
                  nonCommonVertexId = vertexIds[k];
                  break;
                }
              }

              // we have a neighbor
              // now we want to evaluate the angles of the following triangles
              // nonCommonVertexId, commonVertexId0, otherNonCommonVertexId
              // nonCommonVertexId, commonVertexId1, otherNonCommonVertexId
              // if the min angles of these two guys is bigger than alpha, flip!

              if((nonCommonVertexId != -1) && (otherNonCommonVertexId != -1)) {

                vector<double> beta0angles, beta1angles;

                Geometry::computeTriangleAngles(
                  (*globalVertexList_)[nonCommonVertexId].p_.data(),
                  (*globalVertexList_)[commonVertexId0].p_.data(),
                  (*globalVertexList_)[otherNonCommonVertexId].p_.data(),
                  beta0angles);

                double beta0 = -1;
                for(int k = 0; k < 3; k++) {
                  if((beta0 < 0) || (fabs(beta0angles[j]) < beta0))
                    beta0 = fabs(beta0angles[j]);
                }

                Geometry::computeTriangleAngles(
                  (*globalVertexList_)[nonCommonVertexId].p_.data(),
                  (*globalVertexList_)[commonVertexId1].p_.data(),
                  (*globalVertexList_)[otherNonCommonVertexId].p_.data(),
                  beta1angles);

                double beta1 = -1;
                for(int k = 0; k < 3; k++) {
                  if((beta1 < 0) || (fabs(beta1angles[j]) < beta1))
                    beta1 = fabs(beta1angles[j]);
                }

                if((beta0 > alpha) && (beta1 > alpha)) {
                  // flip!

                  if(isEdgeFlippable(commonVertexId0, commonVertexId1,
                                     nonCommonVertexId,
                                     otherNonCommonVertexId)) {

                    // original triangle:
                    // nonCommonVertexId, otherNonCommonVertexId,
                    // commonVertexId0
                    (*polygonEdgeTriangleLists_[triangles[i]
                                                  .first])[triangles[i].second]
                      .vertexIds_[0]
                      = nonCommonVertexId;
                    (*polygonEdgeTriangleLists_[triangles[i]
                                                  .first])[triangles[i].second]
                      .vertexIds_[1]
                      = otherNonCommonVertexId;
                    (*polygonEdgeTriangleLists_[triangles[i]
                                                  .first])[triangles[i].second]
                      .vertexIds_[2]
                      = commonVertexId0;

                    // other triangle:
                    // nonCommonVertexId, otherNonCommonVertexId,
                    // commonVertexId1
                    (*polygonEdgeTriangleLists_[triangles[j]
                                                  .first])[triangles[j].second]
                      .vertexIds_[0]
                      = nonCommonVertexId;
                    (*polygonEdgeTriangleLists_[triangles[j]
                                                  .first])[triangles[j].second]
                      .vertexIds_[1]
                      = otherNonCommonVertexId;
                    (*polygonEdgeTriangleLists_[triangles[j]
                                                  .first])[triangles[j].second]
                      .vertexIds_[2]
                      = commonVertexId1;

                    hasFlipped = true;
                  }
                }
              }
            }
          }
        }
        if(hasFlipped)
          break;
      }
      if(hasFlipped)
        break;
    }

    if(!hasFlipped)
      break;
  }

  return 0;
}

int FiberSurface::getTriangleRangeExtremities(
  const SimplexId &tetId,
  const SimplexId &triangleId,
  const vector<vector<IntersectionTriangle>> &tetIntersections,
  pair<double, double> &extremity0,
  pair<double, double> &extremity1) const {

  vector<double> p0(2), p1(2), p(2);
  vector<double> baryCentrics;
  bool isInBetween = true;

  // check for edges that project to points first
  for(int i = 0; i < 3; i++) {

    p[0] = tetIntersections[tetId][triangleId].uv_[i].first;
    p[1] = tetIntersections[tetId][triangleId].uv_[i].second;

    p0[0] = tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].first;
    p0[1] = tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].second;

    p1[0] = tetIntersections[tetId][triangleId].uv_[(i + 2) % 3].first;
    p1[1] = tetIntersections[tetId][triangleId].uv_[(i + 2) % 3].second;

    if((fabs(p0[0] - p1[0]) < PREC_FLT) && (fabs(p0[1] - p1[1]) < PREC_FLT)) {
      // one edge of the triangle projects to a point
      extremity0.first = p[0];
      extremity0.second = p[1];

      extremity1.first = p0[0];
      extremity1.second = p0[1];

      return 0;
    }
  }

  for(int i = 0; i < 3; i++) {

    p[0] = tetIntersections[tetId][triangleId].uv_[i].first;
    p[1] = tetIntersections[tetId][triangleId].uv_[i].second;

    p0[0] = tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].first;
    p0[1] = tetIntersections[tetId][triangleId].uv_[(i + 1) % 3].second;

    p1[0] = tetIntersections[tetId][triangleId].uv_[(i + 2) % 3].first;
    p1[1] = tetIntersections[tetId][triangleId].uv_[(i + 2) % 3].second;

    Geometry::computeBarycentricCoordinates(
      p0.data(), p1.data(), p.data(), baryCentrics, 2);

    isInBetween = true;
    for(int j = 0; j < 2; j++) {

      if((baryCentrics[j] < -PREC_FLT) || (baryCentrics[j] > 1 + PREC_FLT)) {
        isInBetween = false;
        break;
      }
    }
    if(isInBetween) {
      extremity0.first = p0[0];
      extremity0.second = p0[1];

      extremity1.first = p1[0];
      extremity1.second = p1[1];

      return 0;
    }
  }

  return 0;
}

bool FiberSurface::hasDuplicatedVertices(const double *p0,
                                         const double *p1,
                                         const double *p2) const {

  if((p0[0] == p1[0]) && (p0[1] == p1[1]) && (p0[2] == p1[2]))
    return true;

  if((p2[0] == p1[0]) && (p2[1] == p1[1]) && (p2[2] == p1[2]))
    return true;

  if((p0[0] == p2[0]) && (p0[1] == p2[1]) && (p0[2] == p2[2]))
    return true;

  return false;
}

int FiberSurface::interpolateBasePoints(const vector<double> &p0,
                                        const pair<double, double> &uv0,
                                        const double &t0,
                                        const vector<double> &p1,
                                        const pair<double, double> &uv1,
                                        const double &t1,
                                        const double &t,
                                        Vertex &v) const {

  for(int j = 0; j < 3; j++) {
    v.p_[j] = p0[j] + ((t - t0) / (t1 - t0)) * (p1[j] - p0[j]);
  }

  // interpolate uv
  v.uv_.first = uv0.first + ((t - t0) / (t1 - t0)) * (uv1.first - uv0.first);
  v.uv_.second
    = uv0.second + ((t - t0) / (t1 - t0)) * (uv1.second - uv0.second);

  v.isBasePoint_ = false;

  return 0;
}

bool FiberSurface::isEdgeAngleCollapsible(
  const SimplexId &source,
  const SimplexId &destination,
  const SimplexId &pivotVertexId,
  const vector<pair<SimplexId, SimplexId>> &starNeighbors) const {

  // NOTE:
  // here I should really look at the triangles' normals more than the angles...

  SimplexId baseId = -1;
  double baseAngle = 0;
  for(SimplexId i = 0; i < (SimplexId)starNeighbors.size(); i++) {
    if(((starNeighbors[i].first == source)
        && (starNeighbors[i].second == destination))
       || ((starNeighbors[i].first == destination)
           && (starNeighbors[i].second == source))) {

      baseAngle = Geometry::angle((*globalVertexList_)[source].p_.data(),
                                  (*globalVertexList_)[pivotVertexId].p_.data(),
                                  (*globalVertexList_)[pivotVertexId].p_.data(),
                                  (*globalVertexList_)[destination].p_.data());
      baseId = i;
      break;
    }
  }

  for(SimplexId i = 0; i < (SimplexId)starNeighbors.size(); i++) {
    if(i != baseId) {
      if((starNeighbors[i].first == source)
         || (starNeighbors[i].first == destination)
         || (starNeighbors[i].second == source)
         || (starNeighbors[i].second == destination)) {

        double localAngle = Geometry::angle(
          (*globalVertexList_)[starNeighbors[i].first].p_.data(),
          (*globalVertexList_)[pivotVertexId].p_.data(),
          (*globalVertexList_)[pivotVertexId].p_.data(),
          (*globalVertexList_)[starNeighbors[i].second].p_.data());
        if(localAngle + baseAngle > 0.9 * M_PI)
          return false;
      }
    }
  }
  return true;
}

bool FiberSurface::isEdgeFlippable(const SimplexId &edgeVertexId0,
                                   const SimplexId &edgeVertexId1,
                                   const SimplexId &otherVertexId0,
                                   const SimplexId &otherVertexId1) const {

  double angle0
    = Geometry::angle((*globalVertexList_)[edgeVertexId0].p_.data(),
                      (*globalVertexList_)[edgeVertexId1].p_.data(),
                      (*globalVertexList_)[edgeVertexId1].p_.data(),
                      (*globalVertexList_)[otherVertexId0].p_.data());

  double angle1
    = Geometry::angle((*globalVertexList_)[edgeVertexId0].p_.data(),
                      (*globalVertexList_)[edgeVertexId1].p_.data(),
                      (*globalVertexList_)[edgeVertexId1].p_.data(),
                      (*globalVertexList_)[otherVertexId1].p_.data());

  if(angle0 + angle1 > 0.9 * M_PI)
    return false;

  // now do the angles at the other extremity of the edge.
  angle0 = Geometry::angle((*globalVertexList_)[edgeVertexId1].p_.data(),
                           (*globalVertexList_)[edgeVertexId0].p_.data(),
                           (*globalVertexList_)[edgeVertexId0].p_.data(),
                           (*globalVertexList_)[otherVertexId0].p_.data());

  angle1 = Geometry::angle((*globalVertexList_)[edgeVertexId1].p_.data(),
                           (*globalVertexList_)[edgeVertexId0].p_.data(),
                           (*globalVertexList_)[edgeVertexId0].p_.data(),
                           (*globalVertexList_)[otherVertexId1].p_.data());

  if(angle0 + angle1 > 0.9 * M_PI)
    return false;

  return true;
}

int FiberSurface::mergeEdges(const double &distanceThreshold) const {

  Timer t;

  // TODO: forbid edge collapse for edges that do not live on the boundary of
  // tetrahedrons....
  // this is the case when the tetId of the two triangles is the same.
  // THAT is not OK

  SimplexId initVertexNumber = (*globalVertexList_).size();

  for(SimplexId it = 0; it < (SimplexId)(*globalVertexList_).size(); it++) {
    // avoid infinite loops

    // make the local list of vertex neighbors
    vector<vector<SimplexId>> vertexNeighbors((*globalVertexList_).size());
    vector<vector<SimplexId>> vertexNeighborsTets((*globalVertexList_).size());
    // for each triangle of the star, list of other two vertices
    vector<vector<pair<SimplexId, SimplexId>>> vertexTriangleNeighbors(
      vertexNeighbors.size());

    for(SimplexId i = 0; i < (SimplexId)polygonEdgeTriangleLists_.size(); i++) {
      for(SimplexId j = 0; j < (SimplexId)polygonEdgeTriangleLists_[i]->size();
          j++) {

        for(int k = 0; k < 3; k++) {
          SimplexId vertexId = (*polygonEdgeTriangleLists_[i])[j].vertexIds_[k];

          vertexTriangleNeighbors[vertexId].resize(
            vertexTriangleNeighbors[vertexId].size() + 1);
          vertexTriangleNeighbors[vertexId].back().first = -1;
          vertexTriangleNeighbors[vertexId].back().second = -1;
          vertexNeighborsTets[vertexId].push_back(
            (*polygonEdgeTriangleLists_[i])[j].tetId_);

          for(int l = 0; l < 3; l++) {
            if(l != k) {
              SimplexId otherVertexId
                = (*polygonEdgeTriangleLists_[i])[j].vertexIds_[l];

              if(vertexTriangleNeighbors[vertexId].back().first == -1) {
                vertexTriangleNeighbors[vertexId].back().first = otherVertexId;
              } else {
                vertexTriangleNeighbors[vertexId].back().second = otherVertexId;
              }

              bool isIn = false;
              for(SimplexId m = 0;
                  m < (SimplexId)vertexNeighbors[vertexId].size(); m++) {
                if(vertexNeighbors[vertexId][m] == otherVertexId) {
                  isIn = true;
                  break;
                }
              }
              if(!isIn) {
                vertexNeighbors[vertexId].push_back(otherVertexId);
              }
            }
          }
        }
      }
    }

    bool hasMerged = false;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)polygonEdgeTriangleLists_.size(); i++) {

      for(SimplexId j = 0; j < (SimplexId)polygonEdgeTriangleLists_[i]->size();
          j++) {

        SimplexId minimizer = -1;
        double minDistance = -1;

        // find the smallest edge on this triangle
        for(int k = 0; k < 3; k++) {

          SimplexId vertexId0
            = (*polygonEdgeTriangleLists_[i])[j].vertexIds_[k];
          SimplexId vertexId1
            = (*polygonEdgeTriangleLists_[i])[j].vertexIds_[(k + 1) % 3];

          bool areAlreadySnapped = true;
          for(int l = 0; l < 3; l++) {
            if((*globalVertexList_)[vertexId0].p_[l]
               != (*globalVertexList_)[vertexId1].p_[l]) {

              areAlreadySnapped = false;
              break;
            }
          }

          if(((*globalVertexList_)[vertexId0].isBasePoint_)
             && ((*globalVertexList_)[vertexId1].isBasePoint_)) {
            areAlreadySnapped = true;
          }

          if(!areAlreadySnapped) {
            double distance
              = Geometry::distance((*globalVertexList_)[vertexId0].p_.data(),
                                   (*globalVertexList_)[vertexId1].p_.data());

            if((minDistance == -1) || (distance < minDistance)) {
              minDistance = distance;
              minimizer = k;
            }
          }
        }

        if((minDistance != -1) && (minDistance < distanceThreshold)) {
          SimplexId vertexId0
            = (*polygonEdgeTriangleLists_[i])[j].vertexIds_[minimizer];
          SimplexId vertexId1 = (*polygonEdgeTriangleLists_[i])[j]
                                  .vertexIds_[(minimizer + 1) % 3];

          // find the number of common neighbors
          vector<SimplexId> commonNeighbors;
          for(SimplexId k = 0; k < (SimplexId)vertexNeighbors[vertexId0].size();
              k++) {
            for(SimplexId l = 0;
                l < (SimplexId)vertexNeighbors[vertexId1].size(); l++) {
              if(vertexNeighbors[vertexId0][k]
                 == vertexNeighbors[vertexId1][l]) {
                commonNeighbors.push_back(vertexNeighbors[vertexId0][k]);
              }
            }
          }

          if(commonNeighbors.size() == 2) {

            // we may create extra non-manifold triangles otherwise
            // plus we don't collapse edges that are strictly inside a tet
            // (only collsapse those on the boundary between two tets)

            SimplexId source = vertexId0, destination = vertexId1;

            if((*globalVertexList_)[destination].isBasePoint_) {
              source = vertexId1;
              destination = vertexId0;
            }

            if(!(*globalVertexList_)[destination].isBasePoint_) {

              // now check the angles
              bool isCollapsible = true;
              for(SimplexId k = 0; k < (SimplexId)commonNeighbors.size(); k++) {
                if(!isEdgeAngleCollapsible(
                     source, destination, commonNeighbors[k],
                     vertexTriangleNeighbors[commonNeighbors[k]])) {
                  isCollapsible = false;
                  break;
                }
              }

              // different tets?
              SimplexId tetId0 = -1, tetId1 = -1;
              for(SimplexId k = 0;
                  k < (SimplexId)vertexNeighborsTets[source].size(); k++) {
                if((vertexTriangleNeighbors[source][k].first == destination)
                   || (vertexTriangleNeighbors[source][k].second
                       == destination)) {
                  if(tetId0 == -1) {
                    tetId0 = vertexNeighborsTets[source][k];
                  } else {
                    tetId1 = vertexNeighborsTets[source][k];
                  }
                }
              }

              if(tetId0 == tetId1) {
                isCollapsible = false;
              }

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
              if(isCollapsible) {
                for(int k = 0; k < 3; k++) {
                  (*globalVertexList_)[destination].p_[k]
                    = (*globalVertexList_)[source].p_[k];
                }

                hasMerged = true;
              }
            }
          }
        }
      }
    }

    if(!hasMerged)
      break;

    // now update the vertices
    mergeVertices(0);
  }

  this->printMsg(
    "Performed edge collapses", 1.0, t.getElapsedTime(), this->threadNumber_);
  this->printMsg(std::vector<std::vector<std::string>>{
    {"#Vertices removed",
     std::to_string(initVertexNumber - (*globalVertexList_).size())}});

  return 0;
}

int FiberSurface::mergeVertices(const double &distanceThreshold) const {

  Timer t;

  vector<Vertex> tmpList = (*globalVertexList_);

  for(SimplexId i = 0; i < (SimplexId)tmpList.size(); i++) {
    tmpList[i].localId_ = i;
  }

  const auto FiberSurfaceVertexComparisonX
    = [](const FiberSurface::Vertex &v0, const FiberSurface::Vertex &v1) {
        if(fabs(v0.p_[0] - v1.p_[0]) < PREC_DBL) {
          // let's consider x coordinates are equal
          if(fabs(v0.p_[1] - v1.p_[1]) < PREC_DBL) {
            // let's consider y coordinates are equal
            if(fabs(v0.p_[2] - v1.p_[2]) < PREC_DBL) {
              // let's consider z coordinates are equal
              // NOTE: the local Id should be sufficient
              return v0.globalId_ < v1.globalId_;
            } else
              return v0.p_[2] < v1.p_[2];
          } else
            return v0.p_[1] < v1.p_[1];
        } else {
          return v0.p_[0] < v1.p_[0];
        }
      };

  const auto FiberSurfaceVertexComparisonY
    = [](const FiberSurface::Vertex &v0, const FiberSurface::Vertex &v1) {
        if(fabs(v0.p_[1] - v1.p_[1]) < PREC_DBL) {
          // let's consider y coordinates are equal
          if(fabs(v0.p_[2] - v1.p_[2]) < PREC_DBL) {
            // let's consider z coordinates are equal
            if(fabs(v0.p_[0] - v1.p_[0]) < PREC_DBL) {
              // let's consider x coordinates are equal
              // NOTE: the local Id should be sufficient
              return v0.globalId_ < v1.globalId_;
            } else
              return v0.p_[0] < v1.p_[0];
          } else
            return v0.p_[2] < v1.p_[2];
        } else {
          return v0.p_[1] < v1.p_[1];
        }
      };

  const auto FiberSurfaceVertexComparisonZ
    = [](const FiberSurface::Vertex &v0, const FiberSurface::Vertex &v1) {
        if(fabs(v0.p_[2] - v1.p_[2]) < PREC_DBL) {
          // let's consider z coordinates are equal
          if(fabs(v0.p_[0] - v1.p_[0]) < PREC_DBL) {
            // let's consider x coordinates are equal
            if(fabs(v0.p_[1] - v1.p_[1]) < PREC_DBL) {
              // let's consider y coordinates are equal
              // NOTE: the local Id should be sufficient
              return v0.globalId_ < v1.globalId_;
            } else
              return v0.p_[1] < v1.p_[1];
          } else
            return v0.p_[0] < v1.p_[0];
        } else {
          return v0.p_[2] < v1.p_[2];
        }
      };

  // now do a parallel sort
  SimplexId uniqueVertexNumber = 0;
  for(int k = 0; k < 3; k++) {

    switch(k) {
      case 0:
        sort(tmpList.begin(), tmpList.end(), FiberSurfaceVertexComparisonX);
        break;

      case 1:
        sort(tmpList.begin(), tmpList.end(), FiberSurfaceVertexComparisonY);
        break;

      case 2:
        sort(tmpList.begin(), tmpList.end(), FiberSurfaceVertexComparisonZ);
        break;
    }

    // NOTE:
    // ethaneDiol dataset
    // 0.019224 in parallel (2 cores HT)
    // 0.041353 in sequential (speedup: x2.15, perfect (HT))

    // now merge the thing
    // 1. Identity duplicates
    // NOTE: order is important, so no parallelism
    // NOTE: points are represented with single precision (floats) so use that
    // as a limit for distances
    uniqueVertexNumber = 0;
    double distance = 0;
    for(SimplexId i = 0; i < (SimplexId)tmpList.size(); i++) {

      bool canMerge = false;

      if(i) {
        distance
          = Geometry::distance(tmpList[i].p_.data(), tmpList[i - 1].p_.data());

        if(distance <= distanceThreshold) {

          // one of the two is -1 (interior vertex)
          //           if((tmpList[i - 1].meshEdge_.first == -1)
          //             ||(tmpList[i].meshEdge_.first == -1)){
          //             canMerge = true;
          //           }
          //
          //           // one vertex in common
          //           if((tmpList[i].meshEdge_.first ==
          //             tmpList[i - 1].meshEdge_.first)
          //             ||
          //             (tmpList[i].meshEdge_.first ==
          //             tmpList[i - 1].meshEdge_.second)){
          //             canMerge = true;
          //           }
          //
          //           // one vertex in common
          //           if((tmpList[i].meshEdge_.second ==
          //             tmpList[i - 1].meshEdge_.first)
          //             ||
          //             (tmpList[i].meshEdge_.second ==
          //             tmpList[i - 1].meshEdge_.second)){
          //             canMerge = true;
          //           }

          // NOTE: still some bugs here in terms of manifoldness.

          if((tmpList[i].meshEdge_.first == tmpList[i - 1].meshEdge_.first)
             && (tmpList[i].meshEdge_.second
                 == tmpList[i - 1].meshEdge_.second)) {
            canMerge = true;
          }
          if((tmpList[i].meshEdge_.first == tmpList[i - 1].meshEdge_.second)
             && (tmpList[i].meshEdge_.second
                 == tmpList[i - 1].meshEdge_.first)) {
            canMerge = true;
          }

          //           canMerge = true;
        }

        if(canMerge) {
          tmpList[i].globalId_ = tmpList[i - 1].globalId_;
          if((tmpList[i].isBasePoint_) || (tmpList[i - 1].isBasePoint_)) {
            tmpList[i].isBasePoint_ = true;
            tmpList[i - 1].isBasePoint_ = true;
          }
          if((tmpList[i].isIntersectionPoint_)
             || (tmpList[i - 1].isIntersectionPoint_)) {
            tmpList[i].isIntersectionPoint_ = true;
            tmpList[i - 1].isIntersectionPoint_ = true;
          }
          for(int j = 0; j < 3; j++) {
            tmpList[i].p_[j] = tmpList[i - 1].p_[j];
          }

          // update meshEdge...
          if((tmpList[i - 1].meshEdge_.first == -1)
             && (tmpList[i].meshEdge_.first != -1)) {
            tmpList[i - 1].meshEdge_ = tmpList[i].meshEdge_;
          }
          if((tmpList[i].meshEdge_.first == -1)
             && (tmpList[i - 1].meshEdge_.first != -1)) {
            tmpList[i].meshEdge_ = tmpList[i - 1].meshEdge_;
          }
        }
      }

      if((!i) || (!canMerge)) {
        tmpList[i].globalId_ = uniqueVertexNumber;
        uniqueVertexNumber++;
      }

      (*globalVertexList_)[tmpList[i].localId_].globalId_
        = tmpList[i].globalId_;
    }
  }

  // update the triangles
  for(SimplexId i = 0; i < (SimplexId)polygonEdgeTriangleLists_.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)polygonEdgeTriangleLists_[i]->size();
        j++) {
      for(int k = 0; k < 3; k++) {
        (*polygonEdgeTriangleLists_[i])[j].vertexIds_[k]
          = (*globalVertexList_)[(*polygonEdgeTriangleLists_[i])[j]
                                   .vertexIds_[k]]
              .globalId_;
      }
    }
  }

  // 2. create the actual global list
  // NOTE: order is no longer important, we can do this in parallel.
  (*globalVertexList_).resize(uniqueVertexNumber);

  // duplicate vertices but only a subset have the right edge information
  // (because they have actually been computed on the edge)
  SimplexId lastId = -1;
  for(SimplexId i = 0; i < (SimplexId)tmpList.size(); i++) {

    if(lastId != -1) {
      if((tmpList[lastId].globalId_ == tmpList[i].globalId_)
         && (tmpList[i].meshEdge_.first != -1)) {

        (*globalVertexList_)[tmpList[lastId].globalId_].meshEdge_
          = tmpList[i].meshEdge_;
        (*globalVertexList_)[tmpList[lastId].globalId_].isBasePoint_ = true;
        lastId = -1;
      }
    }

    if((!i) || (tmpList[i].globalId_ != tmpList[i - 1].globalId_)) {
      (*globalVertexList_)[tmpList[i].globalId_] = tmpList[i];
      if((*globalVertexList_)[tmpList[i].globalId_].meshEdge_.first == -1) {
        lastId = i;
      } else {
        lastId = -1;
      }
    }
  }

  // 3. update the 2-sheets, ignore zero-area triangles
  // and free their vertex memory
  // NOTE: order is no longer important, parallel
  vector<vector<bool>> keepTriangle(polygonEdgeTriangleLists_.size());
  vector<vector<FiberSurface::Triangle>> tmpTriangleLists(
    polygonEdgeTriangleLists_.size());
  for(SimplexId i = 0; i < (SimplexId)polygonEdgeTriangleLists_.size(); i++) {
    tmpTriangleLists[i] = (*polygonEdgeTriangleLists_[i]);
    keepTriangle[i].resize((*polygonEdgeTriangleLists_[i]).size(), true);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)polygonEdgeTriangleLists_.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)tmpTriangleLists[i].size(); j++) {
      for(int k = 0; k < 3; k++) {
        if((k)
           && (tmpTriangleLists[i][j].vertexIds_[k]
               == tmpTriangleLists[i][j].vertexIds_[(k - 1)])) {
          keepTriangle[i][j] = false;
        }
      }
      if(tmpTriangleLists[i][j].vertexIds_[0]
         == tmpTriangleLists[i][j].vertexIds_[2]) {
        keepTriangle[i][j] = false;
      }
    }

    // now copy triangles over with non zero-area triangles
    // NOTE: no need to re-allocate the memory, we know we are not going to use
    // more.
    (*polygonEdgeTriangleLists_[i]).resize(0);
    for(SimplexId j = 0; j < (SimplexId)tmpTriangleLists[i].size(); j++) {

      if(keepTriangle[i][j]) {
        (*polygonEdgeTriangleLists_[i]).push_back(tmpTriangleLists[i][j]);
      }
    }
  }

  this->printMsg(
    "Output made manifold", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

int FiberSurface::snapToBasePoint(const vector<vector<double>> &basePoints,
                                  const vector<pair<double, double>> &uv,
                                  const vector<double> &t,
                                  Vertex &v) const {

  if(!pointSnappingThreshold_)
    return -1;

  SimplexId minimizer = 0;
  double minDistance = -1;

  for(SimplexId i = 0; i < (SimplexId)basePoints.size(); i++) {
    double distance = Geometry::distance(basePoints[i].data(), v.p_.data());
    if((minDistance < 0) || (distance < minDistance)) {
      minDistance = distance;
      minimizer = i;
    }
  }

  if(minDistance < pointSnappingThreshold_) {
    for(int i = 0; i < 3; i++) {
      v.p_[i] = basePoints[minimizer][i];
    }
    v.uv_ = uv[minimizer];
    v.t_ = t[minimizer];
    v.isBasePoint_ = true;
  }

  return 0;
}

int FiberSurface::snapVertexBarycentrics() const {

  vector<bool> inQueue(tetNumber_, false);
  vector<vector<pair<SimplexId, SimplexId>>> tetTriangles(tetNumber_);
  vector<SimplexId> tetList;

  for(SimplexId i = 0; i < (SimplexId)polygonEdgeTriangleLists_.size(); i++) {

    for(SimplexId j = 0; j < (SimplexId)polygonEdgeTriangleLists_[i]->size();
        j++) {

      SimplexId tetId = (*polygonEdgeTriangleLists_[i])[j].tetId_;

      if(!inQueue[tetId]) {
        inQueue[tetId] = true;
        tetList.push_back(tetId);
      }

      tetTriangles[tetId].push_back(pair<SimplexId, SimplexId>(i, j));
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)tetList.size(); i++) {
    snapVertexBarycentrics(tetList[i], tetTriangles[tetList[i]]);
  }

  mergeVertices(0);

  return 0;
}

int FiberSurface::snapVertexBarycentrics(
  const SimplexId &tetId,
  const std::vector<std::pair<SimplexId, SimplexId>> &triangles) const {

  for(SimplexId i = 0; i < (SimplexId)triangles.size(); i++) {

    Triangle *t = &(
      (*polygonEdgeTriangleLists_[triangles[i].first])[triangles[i].second]);

    for(int j = 0; j < 3; j++) {
      SimplexId vertexId = t->vertexIds_[j];

      // check for each triangle of the tet
      double minimum = -DBL_MAX;
      vector<double> minBarycentrics;
      vector<SimplexId> minimizer(3);

      for(int k = 0; k < 2; k++) {
        for(int l = k + 1; l < 3; l++) {
          for(int m = l + 1; m < 4; m++) {

            SimplexId vertexId0 = tetList_[5 * tetId + 1 + k];
            SimplexId vertexId1 = tetList_[5 * tetId + 1 + l];
            SimplexId vertexId2 = tetList_[5 * tetId + 1 + m];

            vector<double> p0(3), p1(3), p2(3);

            for(int n = 0; n < 3; n++) {
              p0[n] = pointSet_[3 * vertexId0 + n];
              p1[n] = pointSet_[3 * vertexId1 + n];
              p2[n] = pointSet_[3 * vertexId2 + n];
            }

            vector<double> barycentrics;
            Geometry::computeBarycentricCoordinates(
              p0.data(), p1.data(), p2.data(),
              (*globalVertexList_)[vertexId].p_.data(), barycentrics);

            if((barycentrics[0] != -1.0) && (barycentrics[1] != -1.0)
               && (barycentrics[2] != -1.0)) {
              // vertexId lies in the current triangle

              double localMin = -1.0;
              for(int n = 0; n < 3; n++) {
                if((localMin == -1.0) || (fabs(barycentrics[n]) < localMin)) {
                  localMin = fabs(barycentrics[n]);
                }
              }

              if((minimum == -DBL_MAX) || (localMin < minimum)) {
                minimum = localMin;
                minBarycentrics = barycentrics;
                minimizer[0] = vertexId0;
                minimizer[1] = vertexId1;
                minimizer[2] = vertexId2;
              }
            }
          }
        }
      }

      if((minimum != -DBL_MAX) && (minimum < PREC_FLT_2)) {
        double sum = 0;
        int numberOfZeros = 0;
        for(int k = 0; k < 3; k++) {
          if(minBarycentrics[k] < PREC_FLT_2) {
            minBarycentrics[k] = 0;
            numberOfZeros++;
          }
          sum += minBarycentrics[k];
        }

        sum = (1 - sum) / numberOfZeros;

        for(int k = 0; k < 3; k++) {
          if(minBarycentrics[k] >= PREC_FLT_2) {
            minBarycentrics[k] += sum;
          }
        }

        // now do the re-location
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
        for(int k = 0; k < 3; k++) {
          (*globalVertexList_)[vertexId].p_[k]
            = minBarycentrics[0] * ((double)pointSet_[3 * minimizer[0] + k])
              + minBarycentrics[1] * ((double)pointSet_[3 * minimizer[1] + k])
              + minBarycentrics[2] * ((double)pointSet_[3 * minimizer[2] + k]);
        }
#endif
      }
    }
  }

  return 0;
}
