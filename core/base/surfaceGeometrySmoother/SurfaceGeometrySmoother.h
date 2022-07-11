/// \ingroup base
/// \class ttk::SurfaceGeometrySmoother
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2022.
///
/// \brief TTK VTK-filter for smoothing meshes on surfaces.
///
/// ttk::GeometrySmoother with a twist!
/// This class smoothes and projects a 1D or a 2D mesh onto a 2D
/// closed triangulated surface.
///
/// \sa ttkSurfaceGeometrySmoother.cpp %for a usage example.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_casting/">Persistent
///   Generators Casting example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_fertility/">Persistent
///   Generators Fertility example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_skull/">Persistent
///   Generators Skull example</a> \n

#pragma once

// base code includes
#include <Triangulation.h>
#include <VisitedMask.h>

#include <stack>
#include <string>

namespace ttk {

  class SurfaceGeometrySmoother : virtual public Debug {

  public:
    SurfaceGeometrySmoother();
    ~SurfaceGeometrySmoother() override = default;

    inline void preconditionTriangulationToSmooth(
      AbstractTriangulation *const triangulation) {
      if(triangulation != nullptr) {
        triangulation->preconditionVertexNeighbors();
      }
    }
    inline void preconditionTriangulationSurface(
      AbstractTriangulation *const triangulation) {
      if(triangulation != nullptr) {
        triangulation->preconditionEdges();
        triangulation->preconditionVertexNeighbors();
        triangulation->preconditionVertexEdges();
        triangulation->preconditionTriangles();
        triangulation->preconditionVertexTriangles();
        triangulation->preconditionEdgeTriangles();
      }
    }

    template <typename triangulationType0, typename triangulationType1>
    int execute(float *const outputCoords,
                const float *const inputCoords,
                const char *const mask,
                const SimplexId *const vertsId,
                const int nIter,
                const triangulationType0 &triangulationToSmooth,
                const triangulationType1 &triangulationSurface) const;

  protected:
    struct Point : public std::array<float, 3> {
      Point operator+(const Point other) const {
        Point res{};
        res[0] = (*this)[0] + other[0];
        res[1] = (*this)[1] + other[1];
        res[2] = (*this)[2] + other[2];
        return res;
      }
      Point operator*(const float scalar) const {
        Point res{};
        res[0] = (*this)[0] * scalar;
        res[1] = (*this)[1] * scalar;
        res[2] = (*this)[2] * scalar;
        return res;
      }
      Point operator-(Point other) const {
        return *this + other * (-1);
      }
      Point operator/(const float scalar) const {
        return (*this * (1.0F / scalar));
      }
      friend std::ostream &operator<<(std::ostream &os, const Point &pt) {
        return os << '(' << pt[0] << " " << pt[1] << " " << pt[2] << ')';
      }
    };

    template <typename triangulationType0, typename triangulationType1>
    int relaxProject(std::vector<Point> &outputPoints,
                     std::vector<Point> &tmpStorage,
                     std::vector<SimplexId> &nearestVertexId,
                     std::vector<bool> &trianglesTested,
                     std::vector<SimplexId> &visitedTriangles,
                     std::vector<float> &dists,
                     const char *const mask,
                     const triangulationType0 &triangulationToSmooth,
                     const triangulationType1 &triangulationSurface) const;

    /**
     * @brief Stores the findProjection() result
     */
    struct ProjectionResult {
      /** Projection coordinates */
      Point pt;
      /** Nearest vertex id in the surface */
      SimplexId nearestVertex;
      /** Number of triangles processed */
      size_t trianglesChecked;
      /** Projection status */
      bool projSuccess;
    };

    /**
     * @brief Stores the findProjection() input
     */
    struct ProjectionInput {
      /** Input point coordinates */
      Point pt;
      /** Nearest vertex in the surface */
      SimplexId nearestVertex;
    };

    template <typename triangulationType>
    ProjectionResult
      findProjection(const ProjectionInput &pi,
                     VisitedMask &trianglesTested,
                     std::vector<float> &dists,
                     std::stack<SimplexId> &trianglesToTest,
                     const triangulationType &triangulation) const;

    /**
     * @brief Computes the barycenter of a given point's neighbors
     *
     * @param[in] a Input point index
     * @param [in] outputPoints Coordinates storage
     * @param[in] triangulationToSmooth To get neighbors
     * @return Neighbors barycenter coordinates
     */
    template <typename triangulationType>
    inline Point
      relax(const SimplexId a,
            std::vector<ttk::SurfaceGeometrySmoother::Point> &outputPoints,
            const triangulationType &triangulationToSmooth) const {
      Point relaxed{outputPoints[a]};
      const auto nneigh{triangulationToSmooth.getVertexNeighborNumber(a)};
      for(SimplexId i = 0; i < nneigh; ++i) {
        SimplexId neigh{};
        triangulationToSmooth.getVertexNeighbor(a, i, neigh);
        relaxed = relaxed + outputPoints[neigh];
      }
      return relaxed * (1.0F / static_cast<float>(nneigh + 1));
    }

    /**
     * @brief Compute euclidian projection in a triangle plane
     *
     * @param[in] p Point to be projected
     * @param[in] a Triangle vertex coordinates
     * @param[in] normTri Triangle normal vector
     * @return Projection coordinates
     */
    inline Point projectOnTrianglePlane(const Point &p,
                                        const Point &a,
                                        const Point &normTri) const {
      const auto ap{p - a};
      return p - normTri * Geometry::dotProduct(normTri.data(), ap.data());
    }

    /**
     * @brief Compute euclidian projection on a 3D segment
     *
     * @param[in] p Point to be projected
     * @param[in] a First segment vertex coordinates
     * @param[in] b Second segment vertex coordinates
     * @return Projection coordinates
     */
    inline Point
      projectOnEdge(const Point &p, const Point &a, const Point &b) const {
      const auto ab{b - a};
      const auto ap{p - a};
      return a
             + ab * Geometry::dotProduct(ap.data(), ab.data())
                 / Geometry::dotProduct(ab.data(), ab.data());
    }

    /**
     * @brief Find nearest vertex on the surface.
     *
     * @param[in] pa Input point
     * @param[in] dists Pre-allocated distances vector
     * @param[in] triangulation Surface triangulation
     * @return Nearest vertex id on the surface
     */
    template <typename triangulationType>
    inline SimplexId
      getNearestSurfaceVertex(const Point &pa,
                              std::vector<float> &dists,
                              const triangulationType &triangulation) const {
      for(SimplexId i = 0; i < triangulation.getNumberOfVertices(); ++i) {
        Point pv{};
        triangulation.getVertexPoint(i, pv[0], pv[1], pv[2]);
        dists[i] = Geometry::distance(pa.data(), pv.data());
      }
      return std::min_element(dists.begin(), dists.end()) - dists.begin();
    }
  };

} // namespace ttk

template <typename triangulationType>
ttk::SurfaceGeometrySmoother::ProjectionResult
  ttk::SurfaceGeometrySmoother::findProjection(
    const ProjectionInput &pi,
    VisitedMask &trianglesTested,
    std::vector<float> &dists,
    std::stack<SimplexId> &trianglesToTest,
    const triangulationType &triangulation) const {

  ProjectionResult res{pi.pt, pi.nearestVertex, 0, false};

  // clean trianglesToTest
  while(!trianglesToTest.empty()) {
    trianglesToTest.pop();
  }

  // init pipeline by checking in the first triangle around selected vertex
  if(triangulation.getVertexTriangleNumber(res.nearestVertex) > 0) {
    SimplexId next{};
    triangulation.getVertexTriangle(res.nearestVertex, 0, next);
    trianglesToTest.push(next);
  }

  while(!trianglesToTest.empty()) {
    const auto curr = trianglesToTest.top();
    trianglesToTest.pop();

    // skip if already tested
    if(trianglesTested.isVisited_[curr]) {
      continue;
    }

    // get triangle vertices
    std::array<SimplexId, 3> tverts{};
    triangulation.getTriangleVertex(curr, 0, tverts[0]);
    triangulation.getTriangleVertex(curr, 1, tverts[1]);
    triangulation.getTriangleVertex(curr, 2, tverts[2]);

    // get coordinates of triangle vertices (lets name is MNO)
    std::array<Point, 3> mno{}; // [0] is M, [1] is N and [2] is O
    triangulation.getVertexPoint(tverts[0], mno[0][0], mno[0][1], mno[0][2]);
    triangulation.getVertexPoint(tverts[1], mno[1][0], mno[1][1], mno[1][2]);
    triangulation.getVertexPoint(tverts[2], mno[2][0], mno[2][1], mno[2][2]);

    // triangle normal: cross product of two edges
    Point crossP{};
    // mn, mo vectors
    const Point mn = mno[1] - mno[0];
    const Point mo = mno[2] - mno[0];
    // compute mn ^ mo
    Geometry::crossProduct(mn.data(), mo.data(), crossP.data());
    // unitary normal vector
    const Point normTri = crossP / Geometry::magnitude(crossP.data());

    res.pt = this->projectOnTrianglePlane(pi.pt, mno[0], normTri);

    // compute barycentric coords of projection
    Point baryCoords{};
    Geometry::computeBarycentricCoordinates(
      mno[0].data(), mno[1].data(), mno[2].data(), res.pt.data(), baryCoords);

    // check if projection in triangle
    bool inTriangle = true;
    static const float PREC_FLT{powf(10, -FLT_DIG)};

    for(auto &coord : baryCoords) {
      if(coord < -PREC_FLT) {
        inTriangle = false;
      }
      if(coord > 1 + PREC_FLT) {
        inTriangle = false;
      }
    }

    // mark triangle as tested
    trianglesTested.insert(curr);
    res.trianglesChecked++;

    if(inTriangle) {
      res.projSuccess = true;
      // should we check if we have the nearest triangle?
      break;
    }

    // extrema values in baryCoords
    const auto extrema
      = std::minmax_element(baryCoords.begin(), baryCoords.end());

    // find the nearest triangle vertices local ids (with the
    // highest/positive values in baryCoords) from proj
    std::array<SimplexId, 2> vid{
      static_cast<SimplexId>(extrema.second - baryCoords.begin()), 0};
    for(size_t j = 0; j < baryCoords.size(); j++) {
      if(j != static_cast<size_t>(extrema.first - baryCoords.begin())
         && j != static_cast<size_t>(extrema.second - baryCoords.begin())) {
        vid[1] = j;
        break;
      }
    }

    // store vertex with highest barycentric coordinate
    res.nearestVertex = tverts[vid[0]];
    const auto secondNearVert{tverts[vid[1]]};

    // get the triangle edge with the two vertices
    SimplexId edge{};
    const auto nEdges{triangulation.getVertexEdgeNumber(res.nearestVertex)};
    for(SimplexId i = 0; i < nEdges; ++i) {
      triangulation.getVertexEdge(res.nearestVertex, i, edge);
      SimplexId v{};
      triangulation.getEdgeVertex(edge, 0, v);
      if(v == res.nearestVertex) {
        triangulation.getEdgeVertex(edge, 1, v);
      }
      if(v == secondNearVert) {
        break;
      }
    }

    // next triangle to visit
    SimplexId next{};
    const auto nTri{triangulation.getEdgeTriangleNumber(edge)};
    for(SimplexId i = 0; i < nTri; ++i) {
      triangulation.getEdgeTriangle(edge, i, next);
      if(next != curr) {
        break;
      }
    }

    const auto nVisited{trianglesTested.visitedIds_.size()};
    if(nVisited > 1 && next == trianglesTested.visitedIds_[nVisited - 2]) {
      // the relaxed point is probably over the edge separating the
      // current and the last visited triangles
      res.pt = this->projectOnEdge(pi.pt, mno[vid[0]], mno[vid[1]]);
      res.projSuccess = true;
      break;
    }

    if(!trianglesTested.isVisited_[next]) {
      trianglesToTest.push(next);
    }
  }

  const size_t maxTrChecked = 100;

  if(res.projSuccess && res.trianglesChecked > maxTrChecked) {
    res.projSuccess = false;
  }

  if(!res.projSuccess) {
    // replace proj by the nearest vertex?
    res.nearestVertex
      = this->getNearestSurfaceVertex(pi.pt, dists, triangulation);
    triangulation.getVertexPoint(
      res.nearestVertex, res.pt[0], res.pt[1], res.pt[2]);
  }

  return res;
}

template <typename triangulationType0, typename triangulationType1>
int ttk::SurfaceGeometrySmoother::relaxProject(
  std::vector<ttk::SurfaceGeometrySmoother::Point> &outputPoints,
  std::vector<ttk::SurfaceGeometrySmoother::Point> &tmpStorage,
  std::vector<SimplexId> &nearestVertexId,
  std::vector<bool> &trianglesTested,
  std::vector<SimplexId> &visitedTriangles,
  std::vector<float> &dists,
  const char *const mask,
  const triangulationType0 &triangulationToSmooth,
  const triangulationType1 &triangulationSurface) const {

  Timer tm;
  std::stack<SimplexId> trianglesToTest{};

  // main loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  firstprivate(trianglesTested, visitedTriangles, dists, trianglesToTest)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints.size(); i++) {

    // skip computation if i in filtered
    if(mask != nullptr && mask[i] == 0) {
      tmpStorage[i] = outputPoints[i];
      continue;
    }
    tmpStorage[i] = this->relax(i, outputPoints, triangulationToSmooth);

    VisitedMask vm{trianglesTested, visitedTriangles};

    // replace curr in outputPoints_ by its projection
    const auto res
      = this->findProjection(ProjectionInput{tmpStorage[i], nearestVertexId[i]},
                             vm, dists, trianglesToTest, triangulationSurface);

    tmpStorage[i] = res.pt;
    nearestVertexId[i] = res.nearestVertex;
  }

  std::swap(outputPoints, tmpStorage);

  this->printMsg("Projected " + std::to_string(outputPoints.size()) + " points",
                 1.0, tm.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  return 0;
}

template <typename triangulationType0, typename triangulationType1>
int ttk::SurfaceGeometrySmoother::execute(
  float *const outputCoords,
  const float *const inputCoords,
  const char *const mask,
  const SimplexId *const vertsId,
  const int nIter,
  const triangulationType0 &triangulationToSmooth,
  const triangulationType1 &triangulationSurface) const {

  const auto nPoints{triangulationToSmooth.getNumberOfVertices()};
  if(triangulationSurface.getDimensionality() != 2) {
    this->printErr("Can only project onto a surface");
    return -1;
  }

  if(triangulationToSmooth.getDimensionality() < 1
     || triangulationToSmooth.getDimensionality() > 2) {
    this->printErr("Can only project a 1D or a 2D triangulated object");
    return -1;
  }

  Timer tm{};
  this->printMsg("Smoothing " + std::to_string(nPoints) + " in "
                 + std::to_string(nIter) + " iterations...");

  // list of triangle IDs already tested
  // (takes more memory to reduce computation time)
  std::vector<bool> trianglesTested(
    triangulationSurface.getNumberOfTriangles(), false);
  std::vector<SimplexId> visitedTriangles{};
  // distance between every mesh point and current point
  std::vector<float> dists(triangulationSurface.getNumberOfVertices());

  // temporary storage
  std::vector<ttk::SurfaceGeometrySmoother::Point> outputPoints(nPoints),
    tmpStorage(nPoints);
  std::vector<SimplexId> nearestVertexId(nPoints);

  // copy input
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nPoints; ++i) {
    outputPoints[i][0] = inputCoords[3 * i + 0];
    outputPoints[i][1] = inputCoords[3 * i + 1];
    outputPoints[i][2] = inputCoords[3 * i + 2];
  }

  // ttkVertexScalarField is optional (helps for instance with
  // MorseSmaleComplex 1-separatrices)
  if(vertsId != nullptr) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < nPoints; ++i) {
      nearestVertexId[i] = vertsId[i];
    }
  } else {
    // generate a ttkVertexScalarField-like point data array using raw
    // euclidian distance between the points to smooth and every
    // vertex of the surface
    Timer tm_nv{};
    this->printMsg("Computing nearest vertices...", debug::Priority::INFO,
                   debug::LineMode::REPLACE);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) firstprivate(dists)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < nPoints; ++i) {
      nearestVertexId[i] = this->getNearestSurfaceVertex(
        outputPoints[i], dists, triangulationSurface);
    }
    this->printMsg("Computed nearest vertices", 1.0, tm_nv.getElapsedTime(),
                   this->threadNumber_);
  }

  for(int i = 0; i < nIter; ++i) {
    this->relaxProject(outputPoints, tmpStorage, nearestVertexId,
                       trianglesTested, visitedTriangles, dists, mask,
                       triangulationToSmooth, triangulationSurface);
  }

  // copy output
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nPoints; ++i) {
    outputCoords[3 * i + 0] = outputPoints[i][0];
    outputCoords[3 * i + 1] = outputPoints[i][1];
    outputCoords[3 * i + 2] = outputPoints[i][2];
  }

  this->printMsg("Smoothed " + std::to_string(nPoints) + " points", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);

  return 0;
}
