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
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n

#pragma once

// base code includes
#include <Dijkstra.h>
#include <Geometry.h>
#include <Quadrangulation.h>
#include <Triangulation.h>

#include <limits>
#include <set>

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
        SurfaceGeometrySmoother{}.preconditionTriangulationSurface(triangl);
      }
    }

    template <typename triangulationType = AbstractTriangulation>
    int execute(const triangulationType &triangulation);

  private:
    using Point = Quadrangulation::Point;
    using Quad = Quadrangulation::Quad;

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
    int subdivise(Quadrangulation &qd, const triangulationType &triangulation);

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
    SimplexId findEdgeMiddle(const std::array<SimplexId, 2> &e,
                             const triangulationType &triangulation) const;

    /**
     * @brief Find a quad barycenter using Dijkstra
     *
     * Minimize the sum of the distance to every vertex of the current quad.
     *
     * @param[in] quad Vector of quad vertices point ids in which to
     * find a barycenter
     *
     * @return TTK identifier of potential barycenter
     */
    SimplexId findQuadBary(std::vector<float> &sum, const Quad &quad) const;

    /**
     * @brief Clear buffers
     */
    void clearData();

    template <typename triangulationType>
    float getBoundingBoxDiagonal(const triangulationType &triangulation) const;

    template <typename triangulationType>
    void computeHausdorff(std::vector<float> &hausdorff,
                          const Quadrangulation &qd,
                          const triangulationType &triangulation) const;

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
    FlatJaggedArray quadNeighbors_{};
    // array of nearest input vertex TTK identifier
    std::vector<SimplexId> nearestVertexIdentifier_{};
    // holds geodesic distance to every other quad vertex sharing a quad
    std::vector<std::vector<float>> vertexDistance_{};

    // array of output quadrangle vertex valences
    std::vector<SimplexId> outputValences_{};
    // density around vertices (exp minus euclidean distance between
    // vertex and its closest neighbor)
    std::vector<float> outputDensity_{};
    // quad mesh difformity around vertices (exp minus ratio between
    // smallest and largest euclidean distance to neighbors)
    std::vector<float> outputDifformity_{};
    // array of output quadrangle vertex type
    // 0 - input (critical) point
    // 1 - edge middle
    // 2 - quadrangle barycenter
    std::vector<SimplexId> outputVertType_{};
    // array of output vertex subdivision level
    std::vector<SimplexId> outputSubdivision_{};

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
  const std::array<SimplexId, 2> &e,
  const triangulationType &triangulation) const {

  SimplexId midId{};
  float minValue{std::numeric_limits<float>::infinity()};

  // euclidean barycenter of a and b
  Point edgeEuclBary = (outputPoints_[e[0]] + outputPoints_[e[1]]) * 0.5F;

  for(size_t i = 0; i < vertexDistance_[e[0]].size(); ++i) {
    float const m = vertexDistance_[e[0]][i];
    float const n = vertexDistance_[e[1]][i];
    // stay on the shortest path between a and b
    float sum = m + n;

    // skip further computation
    if(sum > minValue) {
      continue;
    }

    if(m != std::numeric_limits<float>::infinity()
       && n != std::numeric_limits<float>::infinity()) {
      // try to get the middle of the shortest path
      sum += std::abs(m - n);
    }

    // get the euclidean distance to AB
    Point curr{};
    triangulation.getVertexPoint(i, curr[0], curr[1], curr[2]);
    // try to minimize the euclidean distance to AB too
    sum += Geometry::distance(curr.data(), edgeEuclBary.data());

    // search for the minimizing index
    if(sum < minValue) {
      minValue = sum;
      midId = i;
    }
  }

  return midId;
}

template <typename triangulationType>
int ttk::QuadrangulationSubdivision::subdivise(
  Quadrangulation &qd, const triangulationType &triangulation) {

  // temp storage for quad subdivision
  std::vector<Quad> tmp{};

  Timer tm;

  // avoid reallocation in loop, causing invalid pointers
  outputPoints_.reserve(outputPoints_.size() * 5);

  vertexDistance_.resize(outputPoints_.size());

  // set & precondition quadrangulation object
  qd.setInputPoints(this->outputPoints_.size(), this->outputPoints_.data());
  qd.setInputCells(this->outputQuads_.size(), this->outputQuads_.data());
  qd.preconditionEdges();
  qd.preconditionVertexStars();

  // compute shortest distance from every vertex to all other that share a quad
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_.size(); ++i) {

    // skip if already computed on a coarser subdivision
    if(vertexDistance_[i].empty()) {

      // do not propagate on the whole mesh
      std::set<SimplexId> bounds{};
      const auto ns{qd.getVertexStarNumber(i)};
      for(SimplexId j = 0; j < ns; ++j) {
        const auto cid{qd.getVertexStar(i, j)};
        for(const auto v : this->outputQuads_[cid]) {
          if(v == static_cast<LongSimplexId>(i)) {
            continue;
          }
          bounds.emplace(nearestVertexIdentifier_[v]);
        }
      }

      Dijkstra::shortestPath(nearestVertexIdentifier_[i], triangulation,
                             vertexDistance_[i],
                             {bounds.begin(), bounds.end()});
    }
  }

  std::vector<SimplexId> quadBaryId(this->outputQuads_.size());
  std::vector<float> sum{};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) firstprivate(sum)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < this->outputQuads_.size(); ++i) {
    quadBaryId[i] = this->findQuadBary(sum, this->outputQuads_[i]);
  }

  std::vector<SimplexId> edgeMidId(qd.getNumberOfEdges());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < qd.getNumberOfEdges(); ++i) {
    edgeMidId[i] = this->findEdgeMiddle(qd.getEdge(i), triangulation);
  }

  std::vector<SimplexId> processedEdges(qd.getNumberOfEdges(), -1);

  for(size_t a = 0; a < this->outputQuads_.size(); ++a) {
    const auto &q{this->outputQuads_[a]};

    const auto processEdge = [&](const SimplexId e) -> SimplexId {
      if(processedEdges[e] == -1) {
        const auto midab{edgeMidId[e]};
        Point pt{};
        triangulation.getVertexPoint(midab, pt[0], pt[1], pt[2]);
        /* add new point 3d coordinates to vector of output points */
        this->outputPoints_.emplace_back(pt);
        /* new point is an edge middle */
        this->outputVertType_.emplace_back(1);
        /* store also TTK identifier of triangular mesh vertex */
        this->nearestVertexIdentifier_.emplace_back(midab);
        // store in map
        processedEdges[e] = this->outputPoints_.size() - 1;
      }
      return processedEdges[e];
    };

    const auto ij{processEdge(qd.getCellEdge(a, 0))};
    const auto jk{processEdge(qd.getCellEdge(a, 1))};
    const auto kl{processEdge(qd.getCellEdge(a, 2))};
    const auto li{processEdge(qd.getCellEdge(a, 3))};

    // barycenter TTK identifier
    const auto baryid = quadBaryId[a];
    // barycenter 3D coordinates
    Point bary{};
    triangulation.getVertexPoint(baryid, bary[0], bary[1], bary[2]);

    // barycenter index in outputPoints_
    const LongSimplexId baryIdx = outputPoints_.size();
    outputPoints_.emplace_back(bary);
    outputVertType_.emplace_back(2);
    nearestVertexIdentifier_.emplace_back(baryid);

    // add the four new quads
    tmp.emplace_back(Quad{q[0], ij, baryIdx, li});
    tmp.emplace_back(Quad{q[1], jk, baryIdx, ij});
    tmp.emplace_back(Quad{q[2], kl, baryIdx, jk});
    tmp.emplace_back(Quad{q[3], li, baryIdx, kl});
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
float ttk::QuadrangulationSubdivision::getBoundingBoxDiagonal(
  const triangulationType &triangulation) const {

  std::array<float, 3> pmin{std::numeric_limits<float>::max(),
                            std::numeric_limits<float>::max(),
                            std::numeric_limits<float>::max()};
  std::array<float, 3> pmax{std::numeric_limits<float>::min(),
                            std::numeric_limits<float>::min(),
                            std::numeric_limits<float>::min()};

  for(SimplexId i = 0; i < triangulation.getNumberOfVertices(); ++i) {
    std::array<float, 3> p{};
    triangulation.getVertexPoint(i, p[0], p[1], p[2]);
    pmax[0] = std::max(pmax[0], p[0]);
    pmax[1] = std::max(pmax[1], p[1]);
    pmax[2] = std::max(pmax[2], p[2]);
    pmin[0] = std::min(pmin[0], p[0]);
    pmin[1] = std::min(pmin[1], p[1]);
    pmin[2] = std::min(pmin[2], p[2]);
  }

  return Geometry::distance(pmin.data(), pmax.data());
}

template <typename triangulationType>
void ttk::QuadrangulationSubdivision::computeHausdorff(
  std::vector<float> &hausdorff,
  const Quadrangulation &qd,
  const triangulationType &triangulation) const {

  Timer tm{};

  hausdorff.resize(qd.getNumberOfVertices());

  // compute the minimal distance from every triangulation point to
  // every quadrangulation point

  // compute triangulation bounding box diagonal
  const auto bboxDiag = getBoundingBoxDiagonal(triangulation);

  // closest quadrangulation vertex for every triangulation vertex
  std::vector<SimplexId> nearestQuadVert(triangulation.getNumberOfVertices());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->getThreadNumber())
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nearestQuadVert.size(); ++i) {
    float minDist{std::numeric_limits<float>::infinity()};
    std::array<float, 3> p{};
    triangulation.getVertexPoint(i, p[0], p[1], p[2]);

    for(SimplexId j = 0; j < qd.getNumberOfVertices(); ++j) {
      std::array<float, 3> q{};
      qd.getVertexPoint(j, q[0], q[1], q[2]);
      auto dist = Geometry::distance(p.data(), q.data());
      if(dist < minDist) {
        minDist = dist;
        nearestQuadVert[i] = j;
      }
    }
  }

  // nearest triangulation vertices for each quadrangulation vertex
  std::vector<std::vector<SimplexId>> nearestTriVerts(qd.getNumberOfVertices());
  for(SimplexId i = 0; i < triangulation.getNumberOfVertices(); ++i) {
    nearestTriVerts[nearestQuadVert[i]].emplace_back(i);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->getThreadNumber())
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nearestTriVerts.size(); ++i) {
    std::array<float, 3> q{};
    qd.getVertexPoint(i, q[0], q[1], q[2]);
    float maxDist{};
    for(const auto v : nearestTriVerts[i]) {
      std::array<float, 3> p{};
      triangulation.getVertexPoint(v, p[0], p[1], p[2]);
      const auto dist = Geometry::distance(p.data(), q.data());
      if(dist > maxDist) {
        maxDist = dist;
      }
    }
    hausdorff[i] = maxDist / bboxDiag / nearestQuadVert.size() * 1e8;
  }

  this->printMsg("Computed Hausdorff distance", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
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

  Quadrangulation qd{};
  qd.setThreadNumber(this->threadNumber_);
  qd.setDebugLevel(this->debugLevel_);

  // main loop
  for(size_t i = 0; i < SubdivisionLevel; i++) {
    // subdivise each quadrangle by creating five new points, at the
    // center of each edge (4) and at the barycenter of the four
    // vertices (1).
    subdivise(qd, triangulation);
  }

  qd.setInputPoints(this->outputPoints_.size(), this->outputPoints_.data());
  qd.setInputCells(this->outputQuads_.size(), this->outputQuads_.data());

  // also needed by computeStatistics
  qd.preconditionVertexNeighbors();
  qd.preconditionVertexStars();

  if(this->RelaxationIterations > 0) {

    // smoother mask
    std::vector<char> mask(this->outputPoints_.size(), 1);
    if(this->LockAllInputVertices) {
      // all input vertices (before subdivision)
      for(size_t i = 0; i < this->inputVertexNumber_; ++i) {
        mask[i] = 0;
      }
    } else if(this->LockInputExtrema) {
      // extraordinary vertices only (valence != 4)
      for(SimplexId i = 0; i < qd.getNumberOfVertices(); ++i) {
        if(qd.isVertexExtraordinary(i)) {
          mask[i] = 0;
        }
      }
    }

    SurfaceGeometrySmoother worker{};
    worker.setDebugLevel(this->debugLevel_);
    worker.setThreadNumber(this->threadNumber_);
    worker.execute(reinterpret_cast<float *>(this->outputPoints_.data()),
                   reinterpret_cast<float *>(this->outputPoints_.data()),
                   mask.data(), this->nearestVertexIdentifier_.data(),
                   this->RelaxationIterations, qd, triangulation);
  }

  qd.computeStatistics(this->outputValences_, this->outputDensity_,
                       this->outputDifformity_, this->quadArea_,
                       this->quadDiagsRatio_, this->quadEdgesRatio_,
                       this->quadAnglesRatio_);
  this->computeHausdorff(this->hausdorff_, qd, triangulation);

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
