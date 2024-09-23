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

#include <numeric>
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
        if(triangulation->getDimensionality() == 2) {
          triangulation->preconditionVertexStars();
        }
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

  protected:
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
  };

} // namespace ttk

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
#ifdef TTK_ENABLE_OPENMP4
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
    const auto res = Geometry::findProjection(
      Geometry::ProjectionInput{i, tmpStorage[i], nearestVertexId[i]}, vm,
      dists, trianglesToTest, false, triangulationToSmooth,
      triangulationSurface);

    tmpStorage[i] = Point{res.pt};
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
  this->printMsg("Smoothing " + std::to_string(nPoints) + " points in "
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
    // euclidean distance between the points to smooth and every
    // vertex of the surface
    Timer tm_nv{};
    this->printMsg("Computing nearest vertices...", debug::Priority::INFO,
                   debug::LineMode::REPLACE);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) firstprivate(dists)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < nPoints; ++i) {
      nearestVertexId[i] = Geometry::getNearestSurfaceVertex(
        outputPoints[i].data(), dists, triangulationSurface);
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
