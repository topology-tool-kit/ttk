/// \ingroup base
/// \class ttk::QuadrangulationSubdivision
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date March 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkQuadrangulationSubdivision.cpp % for a usage example.

#pragma once

// base code includes
#include <Dijkstra.h>
#include <Geometry.h>
#include <Triangulation.h>

#include <array>
#include <cassert>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stack>
#include <tuple>

namespace ttk {

  class QuadrangulationSubdivision : virtual public Debug {

  public:
    QuadrangulationSubdivision() {
      this->setDebugMsgPrefix("QuadrangulationSubdivision");
    }

    inline void setSubdivisionLevel(const unsigned int value) {
      SubdivisionLevel = value;
    }
    inline void setRelaxationIterations(const unsigned int value) {
      RelaxationIterations = value;
    }
    inline void setLockInputExtrema(const bool value) {
      LockInputExtrema = value;
    }
    inline void setLockAllInputVertices(const bool value) {
      LockAllInputVertices = value;
    }
    inline void setReverseProjection(const bool value) {
      ReverseProjection = value;
    }
    inline void setShowResError(const bool value) {
      ShowResError = value;
    }
    inline void setHausdorffLevel(const float value) {
      HausdorffLevel = value;
    }
    inline void setInputQuads(void *const address, unsigned int size) {
      inputQuads_ = static_cast<Quad *>(address);
      inputQuadNumber_ = size;
    }
    inline void setInputVertices(void *const address, unsigned int size) {
      inputVertices_ = static_cast<Point *>(address);
      inputVertexNumber_ = size;
    }
    inline void setInputVertexIdentifiers(void *const address,
                                          unsigned int size) {
      nearestVertexIdentifier_.resize(size);
      auto inputVertexIdentifiers = static_cast<SimplexId *>(address);
      for(size_t i = 0; i < size; i++) {
        nearestVertexIdentifier_[i] = inputVertexIdentifiers[i];
      }
    }
    inline void
      preconditionTriangulation(AbstractTriangulation *const triangl) {
      if(triangl != nullptr) {
        vertexNumber_ = triangl->getNumberOfVertices();
        triangl->preconditionVertexNeighbors();
        triangl->preconditionVertexTriangles();
      }
    }

    template <typename triangulationType = AbstractTriangulation>
    int execute(const triangulationType &triangulation);

    inline float *getPointsBuf() {
      return reinterpret_cast<float *>(outputPoints_.data());
    }
    inline size_t getPointsNumber() const {
      return outputPoints_.size();
    }

  private:
    // vtkPoint instance with interleaved coordinates (AoS)
    struct Point {
      float x;
      float y;
      float z;
      Point operator+(const Point other) const {
        Point res{};
        res.x = x + other.x;
        res.y = y + other.y;
        res.z = z + other.z;
        return res;
      }
      Point operator*(const float scalar) const {
        Point res{};
        res.x = x * scalar;
        res.y = y * scalar;
        res.z = z * scalar;
        return res;
      }
      Point operator-(Point other) const {
        return *this + other * (-1);
      }
      Point operator/(const float scalar) const {
        return (*this * (1.0F / scalar));
      }
      friend std::ostream &operator<<(std::ostream &stream, const Point &pt) {
        stream << pt.x << " " << pt.y << " " << pt.z;
        return stream;
      }
    };

    /**
     * @brief Ad-hoc quad data structure (4 vertex ids)
     */
    using Quad = std::array<LongSimplexId, 4>;

    /**
     * @brief Subdivise a quadrangular mesh
     *
     * Add five new points per existing quadrangle:
     * - four at the middle of the four quadrangle edges
     * - one at the quadrangle barycenter
     *
     * The points are generated once, no duplicate is generated when
     * considering two quadrangles sharing an edge.
     *
     * New quadrangles using these new points are then generated
     * out-of-place in an output vector of quadrangles.
     *
     * @return 0 in case of success
     */
    template <typename triangulationType>
    int subdivise(const triangulationType &triangulation);

    /**
     * @brief Project a generated quadrangle vertex into the
     * triangular input mesh
     *
     * Project a subset of the current quadrangular mesh onto the
     * triangular input mesh.
     *
     * @param[in] filtered Set of indices that should not be projected
     * @param[in] lastIter Indicate last projection iteration for
     * post-processing
     * @return 0 in case of success
     */
    template <typename triangulationType>
    int project(const std::set<size_t> &filtered,
                const triangulationType &triangulation,
                bool lastIter = false);

    /**
     * @brief Relax every generated point of a quadrangular mesh
     *
     * Take every generated point of the current quadrangular mesh,
     * and move its position to the barycenter of its neighbors.
     *
     * @param[in] filtered Set of indices that should not be projected
     * @return 0 in case of success
     */
    int relax(const std::set<size_t> &filtered);

    /**
     * @brief Store for every quad vertex its neighbors
     *
     * Each quad vertex should be linked to four other vertices. This
     * functions stores into the quadNeighbors_ member this relation.
     *
     * @param[in] quads Quadrangular mesh to find neighbors in
     * @param[in] secondNeighbors Also store secondary neighbors (quad third
     * vertex)
     *
     * @return 0 in case of success
     */
    int getQuadNeighbors(const std::vector<Quad> &quads,
                         std::vector<std::set<size_t>> &neighbors,
                         bool secondNeighbors = false) const;

    /**
     * @brief Compute the normal of the quadrangulation at point a
     *
     * @param[in] a input index of quadrangle vertex
     *
     * @return normal to quad surface at point a
     */
    Point getQuadNormal(size_t a) const;

    /**
     * @brief Compute the projection in the nearest triangle
     *
     * @param[in] a input index of quadrangle vertex
     * @param[in] forceReverseProj Try reverse projection
     *
     * @return (coordinates of projection, nearest vertex id, number
     * of triangles checked for result, projection id)
     */
    template <typename triangulationType>
    std::tuple<Point, SimplexId, size_t, SimplexId>
      findProjection(size_t a,
                     bool forceReverseProj,
                     const triangulationType &triangulation) const;

    /**
     * @brief Find the middle of a quad edge using Dijkstra
     *
     * Minimize the sum of the distance to the two edge vertices, and the
     * distance absolute difference
     *
     * @param[in] a First edge vertex index
     * @param[in] b Second edge vertex index
     *
     * @return TTK identifier of potential edge middle
     */
    template <typename triangulationType>
    SimplexId findEdgeMiddle(size_t a,
                             size_t b,
                             const triangulationType &triangulation) const;

    /**
     * @brief Find a quad barycenter using Dijkstra
     *
     * Minimize the sum of the distance to every vertex of the current quad.
     *
     * @param[in] quadVertices Vector of quad vertices point ids in which to
     * find a barycenter
     *
     * @return TTK identifier of potential barycenter
     */
    SimplexId findQuadBary(const std::vector<size_t> &quadVertices) const;

    /**
     * @brief Find input vertices with more than 4 neighbors
     *
     * @param[out] output Output set of input extraordinary point indices
     *
     * @return 0 in case of success
     */
    int findExtraordinaryVertices(std::set<size_t> &output) const;

    /**
     * @brief Compute statistics on generated quadrangles
     *
     * Computes:
     * - quadrangle area
     * - diagonals ratio
     * - ratio between the shortest and the longest edges
     * - ratio between the smallest and the biggest angles
     */
    template <typename triangulationType>
    void quadStatistics(const triangulationType &triangulation);

    /**
     * @brief Clear buffers
     */
    void clearData();

  protected:
    // number of vertices in the mesh
    SimplexId vertexNumber_{};

    // number of subdivisions of the input quadrangles
    unsigned int SubdivisionLevel{1};
    // number of relaxation iterations
    unsigned int RelaxationIterations{10};
    // lock input extrema
    bool LockInputExtrema{false};
    // lock all input vertices
    bool LockAllInputVertices{false};
    // projection method
    bool ReverseProjection{false};
    // display result despite error
    bool ShowResError{false};
    // Hausdorff warning level
    float HausdorffLevel{200.F};

    // number of input quadrangles
    unsigned int inputQuadNumber_{};
    // input quadrangles
    Quad *inputQuads_{};

    // number of input points (quad vertices)
    unsigned int inputVertexNumber_{};
    // input quadrangle vertices (3D coordinates)
    Point *inputVertices_{};

    // array of output quadrangles
    std::vector<Quad> outputQuads_{};
    // array of output quadrangle vertices
    std::vector<Point> outputPoints_{};
    // array mapping quadrangle neighbors
    std::vector<std::set<size_t>> quadNeighbors_{};
    // array of nearest input vertex TTK identifier
    std::vector<SimplexId> nearestVertexIdentifier_{};
    // holds geodesic distance to every other quad vertex sharing a quad
    std::vector<std::vector<float>> vertexDistance_{};

    // array of output quadrangle vertex valences
    std::vector<SimplexId> outputValences_{};
    // array of output quadrangle vertex type
    // 0 - input (critical) point
    // 1 - edge middle
    // 2 - quadrangle barycenter
    std::vector<SimplexId> outputVertType_{};
    // array of output vertex subdivision level
    std::vector<SimplexId> outputSubdivision_{};
    // number of triangles checked per quad vertex for the last projection
    std::vector<SimplexId> trianglesChecked_{};
    // last projection success per quad vertex
    // 0 - not projected (critical point)
    // 1 - projection alongside quadrangle normal
    // 2 - projection alongside triangle normal
    // 3 - failed projection
    std::vector<SimplexId> projSucceeded_{};

    // quadrangles statistics
    std::vector<float> quadArea_{};
    std::vector<float> quadDiagsRatio_{};
    std::vector<float> quadEdgesRatio_{};
    std::vector<float> quadAnglesRatio_{};
    std::vector<float> hausdorff_{};
  };
} // namespace ttk

// if the package is a pure template typename, uncomment the following line
// #include                  <QuadrangulationSubdivision.cpp>

template <typename triangulationType>
ttk::SimplexId ttk::QuadrangulationSubdivision::findEdgeMiddle(
  const size_t a,
  const size_t b,
  const triangulationType &triangulation) const {

  std::vector<SimplexId> midId(this->threadNumber_);
  std::vector<float> minValue(
    this->threadNumber_, std::numeric_limits<float>::infinity());

  // euclidian barycenter of a and b
  Point edgeEuclBary = (outputPoints_[a] + outputPoints_[b]) * 0.5F;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < vertexDistance_[a].size(); ++i) {
#ifdef TTK_ENABLE_OPENMP
    const auto tid = omp_get_thread_num();
#else
    const auto tid = 0;
#endif // TTK_ENABLE_OPENMP
    float m = vertexDistance_[a][i];
    float n = vertexDistance_[b][i];
    // stay on the shortest path between a and b
    float sum = m + n;

    // skip further computation
    if(sum > minValue[tid]) {
      continue;
    }

    if(m != std::numeric_limits<float>::infinity()
       && n != std::numeric_limits<float>::infinity()) {
      // try to get the middle of the shortest path
      sum += std::abs(m - n);
    }

    // get the euclidian distance to AB
    Point curr{};
    triangulation.getVertexPoint(i, curr.x, curr.y, curr.z);
    // try to minimize the euclidian distance to AB too
    sum += Geometry::distance(&curr.x, &edgeEuclBary.x);

    // search for the minimizing index
    if(sum < minValue[tid]) {
      minValue[tid] = sum;
      midId[tid] = i;
    }
  }

#ifdef TTK_ENABLE_OPENMP
  for(int i = 1; i < this->threadNumber_; ++i) {
    if(minValue[i] < minValue[0]) {
      minValue[0] = minValue[i];
      midId[0] = midId[i];
    }
  }
#endif // TTK_ENABLE_OPENMP

  return midId[0];
}

template <typename triangulationType>
std::tuple<ttk::QuadrangulationSubdivision::Point,
           ttk::SimplexId,
           size_t,
           ttk::SimplexId>
  ttk::QuadrangulationSubdivision::findProjection(
    const size_t a,
    const bool forceReverseProj,
    const triangulationType &triangulation) const {

  static const float PREC_FLT{powf(10, -FLT_DIG)};

  // current vertex 3d coordinates
  Point pa = outputPoints_[a];

  Point res{};

  // fallback to euclidian projection code if no normals
  bool doReverseProj = forceReverseProj;

  // quad normal in a
  Point normalsMean{};

  if(doReverseProj) {
    // compute mean of normals
    normalsMean = getQuadNormal(a);

    doReverseProj = !std::isnan(normalsMean.x);
  }

  // found a projection in one triangle
  bool success = false;
  // list of triangle IDs to test to find a potential projection
  std::stack<SimplexId> trianglesToTest;
  // list of triangle IDs already tested
  // (takes more memory to reduce computation time)
  std::vector<bool> trianglesTested(
    triangulation.getNumberOfTriangles(), false);
  // number of triangles tested
  size_t trChecked{0};
  // vertex in triangle with highest barycentric coordinate
  SimplexId nearestVertex = nearestVertexIdentifier_[a];

  // number of triangles around nearest vertex
  SimplexId triangleNumber
    = triangulation.getVertexTriangleNumber(nearestVertex);
  // init pipeline by checking in every triangle around selected vertex
  for(SimplexId j = 0; j < triangleNumber; j++) {
    SimplexId ntid;
    triangulation.getVertexTriangle(nearestVertex, j, ntid);
    trianglesToTest.push(ntid);
  }

  while(!trianglesToTest.empty()) {
    SimplexId i = trianglesToTest.top();
    trianglesToTest.pop();

    // skip if already tested
    if(trianglesTested[i]) {
      continue;
    }

    // get triangle vertices
    std::array<SimplexId, 3> tverts{};
    triangulation.getTriangleVertex(i, 0, tverts[0]);
    triangulation.getTriangleVertex(i, 1, tverts[1]);
    triangulation.getTriangleVertex(i, 2, tverts[2]);

    // get coordinates of triangle vertices
    Point pm{}, pn{}, po{};
    triangulation.getVertexPoint(tverts[0], pm.x, pm.y, pm.z);
    triangulation.getVertexPoint(tverts[1], pn.x, pn.y, pn.z);
    triangulation.getVertexPoint(tverts[2], po.x, po.y, po.z);

    // triangle normal: cross product of two edges
    Point crossP{};
    // mn, mo vectors
    Point mn = pn - pm;
    Point mo = po - pm;
    // compute mn ^ mo
    Geometry::crossProduct(&mn.x, &mo.x, &crossP.x);
    // unitary normal vector
    Point normTri = crossP / Geometry::magnitude(&crossP.x);

    // compute intersection of triangle plane and line (a, normalsMean)
    if(doReverseProj) {

      auto denom = Geometry::dotProduct(&normalsMean.x, &normTri.x);

      // check if triangle plane is parallel to quad normal
      if(std::abs(denom) < PREC_FLT) {
        // skip this iteration after filling pipeline
        trianglesTested[i] = true;
        // fill pipeline with neighboring triangles
        for(auto &vert : tverts) {
          auto ntr = triangulation.getVertexTriangleNumber(vert);
          for(SimplexId j = 0; j < ntr; ++j) {
            SimplexId tid;
            triangulation.getVertexTriangle(vert, j, tid);
            if(tid != i) {
              trianglesToTest.push(tid);
            }
          }
        }
        continue;
      }

      // use formula from Wikipedia: line-plane intersection
      auto tmp = pm - pa;
      auto alpha = Geometry::dotProduct(&tmp.x, &normTri.x) / denom;

      // intersection
      res = pa + normalsMean * alpha;

    }
    // compute euclidian projection of a in triangle plane
    else {

      auto tmp = pa - pm;
      // projection
      res = pa - normTri * Geometry::dotProduct(&normTri.x, &tmp.x);
    }

    // compute barycentric coords of projection
    std::vector<float> baryCoords;
    Geometry::computeBarycentricCoordinates(
      &pm.x, &pn.x, &po.x, &res.x, baryCoords);

    // check if projection in triangle
    bool inTriangle = true;
    for(auto &coord : baryCoords) {
      if(coord < PREC_FLT) {
        inTriangle = false;
      }
      if(coord > 1 + PREC_FLT) {
        inTriangle = false;
      }
    }

    // mark triangle as tested
    trianglesTested[i] = true;
    trChecked++;

    if(inTriangle) {
      success = true;
      // should we check if we have the nearest triangle?
      break;
    }

    // extrema values in baryCoords
    auto extrema = std::minmax_element(baryCoords.begin(), baryCoords.end());

    // find the nearest triangle vertices (with the highest/positive
    // values in baryCoords) from proj
    std::vector<SimplexId> vertices(2);
    vertices[0] = tverts[extrema.second - baryCoords.begin()];
    for(size_t j = 0; j < baryCoords.size(); j++) {
      if(j != static_cast<size_t>(extrema.first - baryCoords.begin())
         && j != static_cast<size_t>(extrema.second - baryCoords.begin())) {
        vertices[1] = tverts[j];
        break;
      }
    }

    // store vertex with highest barycentric coordinate
    nearestVertex = vertices[0];

    // triangles around vertices[0] and vertices[1]
    std::array<std::set<SimplexId>, 2> vertsTriangles{};

    // get triangles around vertices
    for(size_t j = 0; j < vertices.size(); ++j) {
      SimplexId tnum = triangulation.getVertexTriangleNumber(vertices[j]);
      for(SimplexId k = 0; k < tnum; k++) {
        SimplexId tid;
        triangulation.getVertexTriangle(vertices[j], k, tid);
        if(tid == i) {
          continue;
        }
        vertsTriangles[j].insert(tid);
      }
    }

    // triangles to test next
    std::vector<SimplexId> common_triangles;

    // look for triangles sharing the vertices with max values in baryCoords
    std::set_intersection(vertsTriangles[0].begin(), vertsTriangles[0].end(),
                          vertsTriangles[1].begin(), vertsTriangles[1].end(),
                          std::back_inserter(common_triangles));

    for(auto &ntid : common_triangles) {
      if(!trianglesTested[ntid]) {
        trianglesToTest.push(ntid);
      }
    }
  }

  const size_t maxTrChecked = 100;

  if(success && trChecked > maxTrChecked) {
    success = false;
  }

  if(!success) {
    if(!forceReverseProj) {
      return findProjection(a, true, triangulation);
    }
    // replace proj by the nearest vertex?
    std::vector<float> dists(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      Point pv{};
      triangulation.getVertexPoint(i, pv.x, pv.y, pv.z);
      dists[i] = Geometry::distance(&pa.x, &pv.x);
    }
    auto min = std::min_element(dists.begin(), dists.end()) - dists.begin();
    triangulation.getVertexPoint(min, res.x, res.y, res.z);
    nearestVertex = min;
  }

  SimplexId projSucess = success ? (doReverseProj ? 1 : 2) : 3;

  return std::make_tuple(res, nearestVertex, trChecked, projSucess);
}

template <typename triangulationType>
int ttk::QuadrangulationSubdivision::project(
  const std::set<size_t> &filtered,
  const triangulationType &triangulation,
  const bool lastIter) {
  Timer tm;

  if(lastIter) {
    trianglesChecked_.clear();
    projSucceeded_.clear();
    trianglesChecked_.resize(outputPoints_.size());
    projSucceeded_.resize(outputPoints_.size());
  }

  // temp storage for projected points
  std::vector<Point> tmp(outputPoints_.size());

  // main loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_.size(); i++) {

    // skip computation if i in filtered
    if(filtered.find(i) != filtered.end()) {
      tmp[i] = outputPoints_[i];
      continue;
    }

    // replace curr in outputPoints_ by its projection
    auto res = findProjection(i, ReverseProjection, triangulation);

    tmp[i] = std::get<0>(res);
    nearestVertexIdentifier_[i] = std::get<1>(res);

    if(lastIter) {
      // fill in debug info
      trianglesChecked_[i] = std::get<2>(res);
      projSucceeded_[i] = std::get<3>(res);
    }
  }

  outputPoints_ = std::move(tmp);

  this->printMsg("Projected "
                   + std::to_string(outputPoints_.size() - filtered.size())
                   + " points",
                 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::QuadrangulationSubdivision::subdivise(
  const triangulationType &triangulation) {

  using edgeType = std::pair<LongSimplexId, LongSimplexId>;
  using vertexType = std::pair<LongSimplexId, Point>;
  std::map<edgeType, vertexType> processedEdges;

  // temp storage for quad subdivision
  std::vector<Quad> tmp{};

  Timer tm;

  // avoid reallocation in loop, causing invalid pointers
  outputPoints_.reserve(outputPoints_.size() * 5);

  vertexDistance_.resize(outputPoints_.size());

  // get all other vertices sharing a quad
  quadNeighbors_.clear();
  quadNeighbors_.resize(outputPoints_.size());
  getQuadNeighbors(outputQuads_, quadNeighbors_, true);

  // compute shortest distance from every vertex to all other that share a quad
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_.size(); ++i) {

    // skip if already computed on a coarser subdivision
    if(vertexDistance_[i].empty()) {

      // do not propagate on the whole mesh
      std::vector<SimplexId> bounds;
      for(auto &p : quadNeighbors_[i]) {
        bounds.emplace_back(nearestVertexIdentifier_[p]);
      }

      Dijkstra::shortestPath(
        nearestVertexIdentifier_[i], triangulation, vertexDistance_[i], bounds);
    }
  }

  for(auto &q : outputQuads_) {

    auto i = static_cast<size_t>(q[0]);
    auto j = static_cast<size_t>(q[1]);
    auto k = static_cast<size_t>(q[2]);
    auto l = static_cast<size_t>(q[3]);

    // middles of edges
    auto ijid = findEdgeMiddle(i, j, triangulation);
    auto jkid = findEdgeMiddle(j, k, triangulation);
    auto klid = findEdgeMiddle(k, l, triangulation);
    auto liid = findEdgeMiddle(l, i, triangulation);

    Point midij{};
    triangulation.getVertexPoint(ijid, midij.x, midij.y, midij.z);
    Point midjk{};
    triangulation.getVertexPoint(jkid, midjk.x, midjk.y, midjk.z);
    Point midkl{};
    triangulation.getVertexPoint(klid, midkl.x, midkl.y, midkl.z);
    Point midli{};
    triangulation.getVertexPoint(liid, midli.x, midli.y, midli.z);

    std::vector<size_t> quadVertices{i, j, k, l};
    // barycenter TTK identifier
    auto baryid = findQuadBary(quadVertices);
    // barycenter 3D coordinates
    Point bary{};
    triangulation.getVertexPoint(baryid, bary.x, bary.y, bary.z);

    // order edges to avoid duplicates (ij vs. ji)
    auto ij = std::make_pair(std::min(q[0], q[1]), std::max(q[0], q[1]));
    auto jk = std::make_pair(std::min(q[1], q[2]), std::max(q[1], q[2]));
    auto kl = std::make_pair(std::min(q[2], q[3]), std::max(q[2], q[3]));
    auto li = std::make_pair(std::min(q[3], q[0]), std::max(q[3], q[0]));

    auto process_edge_middle
      = [&](const std::pair<LongSimplexId, LongSimplexId> &pair,
            const Point &pt, const SimplexId id) {
          /* check if edge already processed by a neighbor quad */
          if(processedEdges.find(pair) == processedEdges.end()) {
            processedEdges[pair] = std::make_pair(outputPoints_.size(), pt);
            /* add new point 3d coordinates to vector of output points */
            outputPoints_.emplace_back(pt);
            /* new point is an edge middle */
            outputVertType_.emplace_back(1);
            /* store also TTK identifier of triangular mesh vertex */
            nearestVertexIdentifier_.emplace_back(id);
          }
        };

    process_edge_middle(ij, midij, ijid);
    process_edge_middle(jk, midjk, jkid);
    process_edge_middle(kl, midkl, klid);
    process_edge_middle(li, midli, liid);

    // barycenter index in outputPoints_
    auto baryIdx = static_cast<LongSimplexId>(outputPoints_.size());
    outputPoints_.emplace_back(bary);
    outputVertType_.emplace_back(2);
    nearestVertexIdentifier_.emplace_back(baryid);

    // add the four new quads
    tmp.emplace_back(
      Quad{q[0], processedEdges[ij].first, baryIdx, processedEdges[li].first});
    tmp.emplace_back(
      Quad{q[1], processedEdges[jk].first, baryIdx, processedEdges[ij].first});
    tmp.emplace_back(
      Quad{q[2], processedEdges[kl].first, baryIdx, processedEdges[jk].first});
    tmp.emplace_back(
      Quad{q[3], processedEdges[li].first, baryIdx, processedEdges[kl].first});
  }

  // output subdivision level
  auto currSubd = outputSubdivision_.back() + 1;
  auto subdBeg = outputSubdivision_.size();
  outputSubdivision_.resize(outputPoints_.size());
  std::fill(
    outputSubdivision_.begin() + subdBeg, outputSubdivision_.end(), currSubd);

  this->printMsg("Subdivised " + std::to_string(outputQuads_.size())
                   + " quads into " + std::to_string(tmp.size())
                   + " new quads (" + std::to_string(outputPoints_.size())
                   + " points)",
                 1.0, tm.getElapsedTime(), this->threadNumber_);

  outputQuads_ = std::move(tmp);

  return 0;
}

template <typename triangulationType>
void ttk::QuadrangulationSubdivision::quadStatistics(
  const triangulationType &triangulation) {
  Timer tm;

  quadArea_.clear();
  quadArea_.resize(outputQuads_.size());
  quadDiagsRatio_.clear();
  quadDiagsRatio_.resize(outputQuads_.size());
  quadEdgesRatio_.clear();
  quadEdgesRatio_.resize(outputQuads_.size());
  quadAnglesRatio_.clear();
  quadAnglesRatio_.resize(outputQuads_.size());
  hausdorff_.clear();
  hausdorff_.resize(outputPoints_.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputQuads_.size(); ++i) {
    const auto &q = outputQuads_[i];
    Point pi = outputPoints_[q[0]];
    Point pj = outputPoints_[q[1]];
    Point pk = outputPoints_[q[2]];
    Point pl = outputPoints_[q[3]];

    // quadrangle area
    float area0{}, area1{};
    Geometry::computeTriangleArea(&pi.x, &pj.x, &pk.x, area0);
    Geometry::computeTriangleArea(&pi.x, &pk.x, &pl.x, area1);
    quadArea_[i] = area0 + area1;

    // diagonals ratio
    auto diag0 = Geometry::distance(&pi.x, &pk.x);
    auto diag1 = Geometry::distance(&pj.x, &pl.x);
    quadDiagsRatio_[i] = std::min(diag0, diag1) / std::max(diag0, diag1);

    // edges ratio
    std::array<float, 4> edges{
      Geometry::distance(&pi.x, &pj.x), // ij
      Geometry::distance(&pj.x, &pk.x), // jk
      Geometry::distance(&pk.x, &pl.x), // kl
      Geometry::distance(&pl.x, &pi.x), // li
    };
    quadEdgesRatio_[i] = *std::min_element(edges.begin(), edges.end())
                         / *std::max_element(edges.begin(), edges.end());

    // angles ratio
    std::array<float, 4> angles{
      Geometry::angle(&pi.x, &pl.x, &pi.x, &pj.x), // lij
      Geometry::angle(&pj.x, &pi.x, &pj.x, &pk.x), // ijk
      Geometry::angle(&pk.x, &pj.x, &pk.x, &pl.x), // jkl
      Geometry::angle(&pl.x, &pk.x, &pl.x, &pi.x), // kli
    };
    quadAnglesRatio_[i] = *std::min_element(angles.begin(), angles.end())
                          / *std::max_element(angles.begin(), angles.end());
  }

  // compute ratio between quad area and mean quad area

  // global surface area
  float sumArea{};
  for(const auto a : quadArea_) {
    sumArea += a;
  }
  for(auto &a : quadArea_) {
    a *= quadArea_.size() / sumArea;
  }

  // compute the minimal distance from every triangulation point to
  // every quadrangulation point

  std::vector<size_t> triVertsDist(vertexNumber_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < static_cast<size_t>(vertexNumber_); ++i) {
    float minDist{std::numeric_limits<float>::infinity()};
    Point p{};
    triangulation.getVertexPoint(i, p.x, p.y, p.z);

    for(size_t j = 0; j < outputPoints_.size(); ++j) {
      auto dist = Geometry::distance(&p.x, &outputPoints_[j].x);
      if(dist < minDist) {
        minDist = dist;
        triVertsDist[i] = j;
      }
    }
  }

  // compute triangulation bounding box diagonal
  Point pmin{std::numeric_limits<float>::infinity(),
             std::numeric_limits<float>::infinity(),
             std::numeric_limits<float>::infinity()};
  Point pmax{-std::numeric_limits<float>::infinity(),
             -std::numeric_limits<float>::infinity(),
             -std::numeric_limits<float>::infinity()};

  for(size_t i = 0; i < static_cast<size_t>(vertexNumber_); ++i) {
    Point p{};
    triangulation.getVertexPoint(i, p.x, p.y, p.z);
    pmax.x = std::max(pmax.x, p.x);
    pmax.y = std::max(pmax.y, p.y);
    pmax.z = std::max(pmax.z, p.z);
    pmin.x = std::min(pmin.x, p.x);
    pmin.y = std::min(pmin.y, p.y);
    pmin.z = std::min(pmin.z, p.z);
  }

  auto bboxDiag = Geometry::distance(&pmin.x, &pmax.x);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_.size(); ++i) {
    float maxDist{};
    for(size_t j = 0; j < static_cast<size_t>(vertexNumber_); ++j) {
      Point p{};
      triangulation.getVertexPoint(j, p.x, p.y, p.z);

      if(triVertsDist[j] == i) {
        auto dist = Geometry::distance(&p.x, &outputPoints_[i].x);
        if(dist > maxDist) {
          maxDist = dist;
        }
      }
    }
    hausdorff_[i] = maxDist / bboxDiag / vertexNumber_ * 1e8;
  }

  this->printMsg("Computed quad statistics", 1.0, tm.getElapsedTime(),
                 debug::LineMode::NEW, debug::Priority::DETAIL);
}

// main routine
template <typename triangulationType>
int ttk::QuadrangulationSubdivision::execute(
  const triangulationType &triangulation) {

  this->printMsg("Beginning computation...");
  Timer tm;

  // clear output variables
  clearData();

  // ensure consistency of dependent options
  if(LockAllInputVertices) {
    LockInputExtrema = true;
  }

  // store input points (MSC critical points)
  for(size_t i = 0; i < inputVertexNumber_; i++) {
    outputPoints_.emplace_back(inputVertices_[i]);
  }

  // copy input quads into vector
  for(size_t i = 0; i < inputQuadNumber_; i++) {
    outputQuads_.emplace_back(inputQuads_[i]);
  }

  // fill outputInfos_ with input data (critical points)
  outputVertType_.resize(outputPoints_.size());
  std::fill(outputVertType_.begin(), outputVertType_.end(), 0);

  // fill outputSubdivision with input data
  outputSubdivision_.resize(outputPoints_.size());
  std::fill(outputSubdivision_.begin(), outputSubdivision_.end(), 0);

  // vertices to filter from relaxation, projection
  std::set<size_t> filtered{};
  if(!LockAllInputVertices) {
    if(LockInputExtrema) {
      // get extraordinary vertices
      findExtraordinaryVertices(filtered);
    }
  } else {
    // fill vector with all input points indices from 0 to inputVertexNumber_
    for(size_t i = 0; i < inputVertexNumber_; ++i) {
      filtered.insert(i);
    }
  }

  // main loop
  for(size_t i = 0; i < SubdivisionLevel; i++) {
    // subdivise each quadrangle by creating five new points, at the
    // center of each edge (4) and at the barycenter of the four
    // vertices (1).
    subdivise(triangulation);
  }

  // retrieve mapping between every vertex and its neighbors
  quadNeighbors_.clear();
  quadNeighbors_.resize(outputPoints_.size());
  getQuadNeighbors(outputQuads_, quadNeighbors_);

  // "relax" the new points, i.e. replace it by the barycenter of its
  // four neighbors
  for(size_t i = 0; i < RelaxationIterations; i++) {
    relax(filtered);

    // project all points on the nearest triangle (except MSC critical
    // points)
    project(filtered, triangulation, (i == RelaxationIterations - 1));
  }

  // compute valence of every quadrangle vertex
  outputValences_.resize(outputPoints_.size());
  std::transform(
    quadNeighbors_.begin(), quadNeighbors_.end(), outputValences_.begin(),
    [&](const std::set<size_t> &neighbors) { return neighbors.size(); });

  quadStatistics(triangulation);

  bool criterion = false;
  for(size_t i = 0; i < outputPoints_.size(); ++i) {
    if(outputValences_[i] > 4) {
      continue;
    }
    if(hausdorff_[i] > HausdorffLevel) {
      criterion = true;
      break;
    }
  }

  if(criterion) {
    // log, clean & early return
    this->printErr("The output quadrangulation exceeds the provided Haussdorff "
                   "distance tolerance");
    if(!ShowResError) {
      clearData();
      return 1;
    }
  }

  this->printMsg("Produced " + std::to_string(outputQuads_.size())
                   + " quadrangles with " + std::to_string(outputPoints_.size())
                   + " points",
                 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 0;
}
