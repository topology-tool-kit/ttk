/// \ingroup base
/// \class ttk::MorseSmaleQuadrangulation
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date March 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkMorseSmaleQuadrangulation.cpp % for a usage example.

#pragma once

// base code includes
#include <BarycentricSubdivision.h>
#include <Dijkstra.h>
#include <Geometry.h>
#include <Triangulation.h>

#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <queue>
#include <set>

namespace ttk {

  class MorseSmaleQuadrangulation : virtual public Debug {
  public:
    MorseSmaleQuadrangulation() {
      this->setDebugMsgPrefix("MorseSmaleQuadrangulation");
    }

    inline void setCriticalPoints(const unsigned int number,
                                  void *const points,
                                  void *const ids,
                                  void *const cellIds,
                                  void *const type) {
      criticalPointsNumber_ = number;
      criticalPoints_ = static_cast<float *>(points);
      criticalPointsIdentifier_ = static_cast<SimplexId *>(ids);
      criticalPointsCellIds_ = static_cast<SimplexId *>(cellIds);
      criticalPointsType_ = static_cast<unsigned char *>(type);
    }

    inline void setSeparatrices(const unsigned int number,
                                void *const cellIds,
                                void *const cellDims,
                                void *const mask,
                                void *const points) {
      separatriceNumber_ = number;
      sepCellIds_ = static_cast<SimplexId *>(cellIds);
      sepCellDims_ = static_cast<unsigned char *>(cellDims);
      sepMask_ = static_cast<unsigned char *>(mask);
      sepPoints_ = static_cast<float *>(points);
    }

    inline void setDualQuadrangulation(const bool input) {
      DualQuadrangulation = input;
    }
    inline void setShowResError(const bool value) {
      ShowResError = value;
    }
    inline void
      preconditionTriangulation(AbstractTriangulation *const triangl) {
      if(triangl != nullptr) {
        triangl->preconditionVertexNeighbors();
        triangl->preconditionVertexTriangles();
        triangl->preconditionBoundaryVertices();
        verticesNumber_ = triangl->getNumberOfVertices();
        bs.preconditionTriangulation(triangl);
      }
    }

    template <typename triangulationType>
    int execute(const triangulationType &triangulation);

  private:
    /**
     * @brief Find the middle of the separatrix specified by its bounds
     *
     * @param[in] a Index in separatrices array of separatrix source
     * @param[in] b Index in separatrices array of separatrix destination
     *
     * @return Index of separatrice source
     */
    template <typename triangulationType>
    size_t findSeparatrixMiddle(const size_t a,
                                const size_t b,
                                const triangulationType &triangulation);

    /**
     * @brief Find the extremities of a set of separatrices
     *
     * @param[in] seps Input vector of separatrices indices
     * @param[out] srcs Output vector of separatrices sources
     * @param[out] dsts Output vector of separatrices destinations
     *
     * @return 0
     */
    int findSepsVertices(const std::vector<size_t> &seps,
                         std::vector<LongSimplexId> &srcs,
                         std::vector<LongSimplexId> &dsts) const;

    /**
     * @brief Perform the quadrangulation
     *
     * The direct quadrangulation links extrema to saddle points to
     * make quadrangles.
     *
     * @param[out] ndegen number of degenerate quadrangles produced
     * @return 0 in case of success
     */
    template <typename triangulationType>
    int quadrangulate(size_t &ndegen, const triangulationType &triangulation);

    /**
     * @brief Perform the dual quadrangulation
     *
     * The dual quadrangulation uses only extrema and no saddle points
     * to output a coarser quadrangulation.
     *
     * @return 0 in case of success
     */
    int dualQuadrangulate();

    /**
     * @brief Subdivise quadrangulation
     *
     * Find duplicate separatrices coming from the same vertices and
     * generate new quads that try to map tubular topologies.
     *
     * @return 0 in case of success
     */
    template <typename triangulationType>
    int subdivise(const triangulationType &triangulation);

    /**
     * @brief Link four separatrices to a cell
     *
     * Perform a breadth-first search from saddle points on a
     * barycentric subdivision of the triangulation to detect the four
     * separatrices around the current cell
     *
     * @return 0
     */
    template <typename triangulationType>
    int detectCellSeps(const triangulationType &triangulation);

    /**
     * @brief Compare the closeness of the input triangulation and the
     * output triangulation
     *
     * @return True if the two surfaces are both close or both open
     */
    template <typename triangulationType>
    bool checkSurfaceCloseness(const triangulationType &triangulation) const;

    /**
     * @brief Clear data
     */
    void clearData();

    /**
     * @brief Ad-hoc quad data structure
     */
    using Quad = std::array<LongSimplexId, 4>;

    /**
     * @brief Subdivise degenerate quads in the quadrangulation
     *
     * @param[out] outputSubd Quad subdivision to be completed
     *
     * @return 0
     */
    template <typename triangulationType>
    int subdiviseDegenerateQuads(std::vector<Quad> &outputSubd,
                                 const triangulationType &triangulation);

    // number of vertices in triangulation
    SimplexId verticesNumber_{};

    // number of critical points from the Morse-Smale complex
    SimplexId criticalPointsNumber_{};
    // critical points 3d coordinates
    float *criticalPoints_{};
    // mapping points id -> cells id
    SimplexId *criticalPointsCellIds_{};
    // mapping point id -> TTK identifier
    SimplexId *criticalPointsIdentifier_{};
    // critical point type: 0 minimum, 1 saddle point, 2 maximum
    unsigned char *criticalPointsType_{};

    // number of separatrices data
    SimplexId separatriceNumber_{};
    // separatrices points cellIds (to be linked to critical points cellIds)
    SimplexId *sepCellIds_{};
    // separatrices mask scalar field (0 for critical points, 1 otherwise)
    unsigned char *sepMask_{};
    // separatrices cell dimension: 0 for vertices, 1 for edges, 2 for triangles
    unsigned char *sepCellDims_{};
    // separatrices points
    float *sepPoints_{};

    // index of separatrices beginnings in separatrices arrays
    std::vector<size_t> sepBegs_{};
    // index of separatrices endings in separatrices arrays
    std::vector<size_t> sepEnds_{};
    // separatrices middles index in output points array
    std::vector<SimplexId> sepMids_{};
    // sub-segmentation of Morse-Smale cells
    std::vector<SimplexId> morseSeg_{};
    // morseSeg_ id -> separatrices that border quads id
    std::vector<std::pair<SimplexId, std::vector<size_t>>> quadSeps_{};

  protected:
    // array of output quads
    std::vector<Quad> outputCells_{};
    // array of output vertices (generated middles of duplicated separatrices)
    std::vector<float> outputPoints_{};
    // array of output vertices identifiers
    std::vector<SimplexId> outputPointsIds_{};
    // 0: critical points, 1: edge middle, 2: quad barycenter
    std::vector<SimplexId> outputPointsTypes_{};
    // for critical points, their id, for sep middles, the sep id and
    // for barycenters the parent quad id
    std::vector<SimplexId> outputPointsCells_{};
    // triangulation subdivision to detect cells
    BarycentricSubdivision bs{};

    // if dual quadrangulation
    bool DualQuadrangulation{false};
    // display result despite error
    bool ShowResError{false};
  };
} // namespace ttk

template <typename triangulationType>
int ttk::MorseSmaleQuadrangulation::detectCellSeps(
  const triangulationType &triangulation) {

  ExplicitTriangulation newT{};
  bs.execute(triangulation, newT);

  newT.setDebugLevel(this->debugLevel_);
  newT.setThreadNumber(this->threadNumber_);
  newT.preconditionEdges();
  newT.preconditionVertexNeighbors();
  newT.preconditionVertexEdges();
  newT.preconditionVertexTriangles();
  newT.preconditionEdgeTriangles();
  newT.preconditionTriangleEdges();

  auto nEdges = triangulation.getNumberOfEdges();

  // store separatrix index on subdivised triangulation edges
  std::vector<SimplexId> edgeOnSep(newT.getNumberOfEdges(), -1);
  // count the number of critical points encountered
  size_t critPoints{0};
  if(sepMask_[0] == 0) {
    critPoints++;
  }

  // id of critical point in new triangulation
  auto critPointId = [&](const SimplexId a) -> SimplexId {
    if(sepCellDims_[a] == 0) {
      return sepCellIds_[a];
    }
    if(sepCellDims_[a] == 1) {
      return verticesNumber_ + sepCellIds_[a];
    }
    if(sepCellDims_[a] == 2) {
      return verticesNumber_ + nEdges + sepCellIds_[a];
    }
    return -1;
  };

  // init first vertex id
  SimplexId prev{critPointId(0)};

  // put separatrices indices onto edges
  // begin loop at one to get the first edge
  for(SimplexId i = 1; i < separatriceNumber_; ++i) {

    // for computing the separatrix id
    if(sepMask_[i] == 0) {
      critPoints++;
      // beginning of a new separatrix
      if(critPoints % 2 == 1) {
        prev = critPointId(i);
        continue;
      }
    }

    // current separatrix id is critPoints // 2
    auto currSepId
      = (critPoints % 2 == 0) ? critPoints / 2 - 1 : critPoints / 2;

    // current point
    SimplexId curr{critPointId(i)};

    // get edge id
    auto ne = newT.getVertexEdgeNumber(curr);
    for(SimplexId j = 0; j < ne; ++j) {
      SimplexId e{};
      newT.getVertexEdge(curr, j, e);
      SimplexId e0{}, e1{};
      newT.getEdgeVertex(e, 0, e0);
      newT.getEdgeVertex(e, 1, e1);

      if((e0 == prev && e1 == curr) || (e1 == prev && e0 == curr)) {
        edgeOnSep[e] = currSepId;
        break;
      }
    }

    prev = curr;
  }

  // if vertex is saddle in the subdivised triangulation
  std::vector<bool> isSaddle(newT.getNumberOfVertices(), false);

  for(SimplexId i = 0; i < criticalPointsNumber_; ++i) {
    // keep only saddle points
    if(criticalPointsType_[i] != 1) {
      continue;
    }
    isSaddle[verticesNumber_ + criticalPointsCellIds_[i]] = true;
  }

  auto getTriangleEdges
    = [&](const SimplexId tr, SimplexId &e0, SimplexId &e1, SimplexId &e2) {
        newT.getTriangleEdge(tr, 0, e0);
        newT.getTriangleEdge(tr, 1, e1);
        newT.getTriangleEdge(tr, 2, e2);
      };

  // get the indices of the two separatrices around a triangle
  // (having a saddle point as vertex)
  auto sepIdAroundTriangle = [&](const SimplexId tr) {
    SimplexId e0{}, e1{}, e2{};
    getTriangleEdges(tr, e0, e1, e2);

    std::set<SimplexId> sepId{};
    if(edgeOnSep[e0] != -1 && edgeOnSep[e1] != -1) {
      sepId.emplace(edgeOnSep[e0]);
      sepId.emplace(edgeOnSep[e1]);
    } else if(edgeOnSep[e1] != -1 && edgeOnSep[e2] != -1) {
      sepId.emplace(edgeOnSep[e1]);
      sepId.emplace(edgeOnSep[e2]);
    } else if(edgeOnSep[e0] != -1 && edgeOnSep[e2] != -1) {
      sepId.emplace(edgeOnSep[e0]);
      sepId.emplace(edgeOnSep[e2]);
    }
    return sepId;
  };

  auto hasTriangleSaddle = [&](const SimplexId tr) {
    SimplexId a, b, c;
    newT.getTriangleVertex(tr, 0, a);
    newT.getTriangleVertex(tr, 1, b);
    newT.getTriangleVertex(tr, 2, c);
    return (isSaddle[a] || isSaddle[b] || isSaddle[c]);
  };

  // propagate from triangles around saddle
  std::vector<SimplexId> processed(newT.getNumberOfTriangles(), -1);

  // store the saddles id
  std::vector<SimplexId> saddles{};

  for(SimplexId i = 0; i < criticalPointsNumber_; ++i) {
    // keep only saddle points
    if(criticalPointsType_[i] == 1) {
      saddles.emplace_back(verticesNumber_ + criticalPointsCellIds_[i]);
    }
  }

  // look around the saddle points
  for(size_t i = 0; i < saddles.size(); ++i) {
    SimplexId saddle = saddles[i];

    auto sadtri = newT.getVertexTriangleNumber(saddle);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId j = 0; j < sadtri; ++j) {
      std::queue<SimplexId> toProcess{};
      SimplexId tr;
      newT.getVertexTriangle(saddle, j, tr);

      std::set<SimplexId> sepIdBeg = sepIdAroundTriangle(tr);
      // current iteration id
      SimplexId iter = i * sadtri + j;

      toProcess.push(tr);

      while(!toProcess.empty()) {
        auto curr = toProcess.front();
        toProcess.pop();

        // mark current triangle, skip already processed
#ifdef TTK_ENABLE_OPENMP
        if(processed[curr] == -1 || processed[curr] > iter) {
#pragma omp atomic write
          processed[curr] = iter;
        }
        // skip if curr marked by thread
        else if(processed[curr] == iter) {
          continue;
        }
        // stop BFS if thread with higher iteration id has reached curr
        else if(processed[curr] < iter) {
          break;
        }
#else
        if(processed[curr] != -1) {
          continue;
        }
        processed[curr] = 0;
#endif // TTK_ENABLE_OPENMP

        // check for saddle at vertices
        bool hasSaddle = hasTriangleSaddle(curr);

        // check for separatrices on edges
        auto sepIdEnd = sepIdAroundTriangle(curr);

        if(hasSaddle && sepIdEnd.size() == 2) {
          // check that seps are different from beginning
          std::vector<size_t> cellSeps{};
          std::set_union(sepIdBeg.begin(), sepIdBeg.end(), sepIdEnd.begin(),
                         sepIdEnd.end(), std::back_inserter(cellSeps));
          if(cellSeps.size() > 2) {
            // found it
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif // TTK_ENABLE_OPENMP
            {
              // keep indices in sync
              quadSeps_.emplace_back(iter, cellSeps);
            }
          }
        }

        // mark vertices from the original mesh with cell id
        for(SimplexId k = 0; k < 3; ++k) {
          SimplexId vert{};
          newT.getTriangleVertex(curr, k, vert);
          if(vert >= verticesNumber_) {
            continue;
          }
          bool vertOnSep = false;
          for(SimplexId l = 0; l < newT.getVertexEdgeNumber(vert); ++l) {
            SimplexId e{};
            newT.getVertexEdge(vert, l, e);
            if(edgeOnSep[e] != -1) {
              vertOnSep = true;
              break;
            }
          }
          if(!vertOnSep) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif // TTK_ENABLE_OPENMP
            morseSeg_[vert] = iter;
          }
        }

        // look for neighboring triangles
        std::array<SimplexId, 3> edges{};
        getTriangleEdges(curr, edges[0], edges[1], edges[2]);
        for(const auto e : edges) {
          // do not cross separatrices
          if(edgeOnSep[e] != -1) {
            continue;
          }
          for(SimplexId k = 0; k < newT.getEdgeTriangleNumber(e); ++k) {
            SimplexId neigh{};
            newT.getEdgeTriangle(e, k, neigh);
            // push only non processed triangles
            if(processed[neigh] == -1) {
              toProcess.push(neigh);
            }
#ifdef TTK_ENABLE_OPENMP
            // push if neigh marked during a previous iteration
            else if(processed[neigh] > iter) {
              toProcess.push(neigh);
            }
            // stop pushing neighbors if curr marked in a newer iteration
            else if(processed[neigh] < iter) {
              break;
            }
#endif // TTK_ENABLE_OPENMP
          }
        }
      }
    }
  }

  if(this->threadNumber_ > 1) {
    // sort quadSeps_ according to cellIds_ to get a deterministic
    // output when filled in parallel
    TTK_PSORT(this->threadNumber_, quadSeps_.begin(), quadSeps_.end());
    // (by default, pairs are sorted by their first element)
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleQuadrangulation::quadrangulate(
  size_t &ndegen, const triangulationType &triangulation) {
  // quadrangle vertices are either extrema or saddle points

  // separatrices bounds indices and cell ids
  std::vector<size_t> sepFlatEdges{};

  for(SimplexId i = 0; i < separatriceNumber_; ++i) {
    if(sepMask_[i] == 1) {
      continue;
    }
    sepFlatEdges.emplace_back(i);
  }

  // number of separatrices
  auto numSeps = sepFlatEdges.size() / 2;

  // clear data members
  sepBegs_.resize(numSeps);
  sepEnds_.resize(numSeps);
  morseSeg_.resize(verticesNumber_);
  std::fill(morseSeg_.begin(), morseSeg_.end(), -1);
  quadSeps_.clear();

  // fill in data arrays
  for(size_t i = 0; i < numSeps; ++i) {
    // separatrices bounds
    sepBegs_[i] = sepFlatEdges[2 * i];
    sepEnds_[i] = sepFlatEdges[2 * i + 1];
  }

  detectCellSeps(triangulation);

  outputCells_.reserve(quadSeps_.size());

  for(const auto &qsp : quadSeps_) {

    const auto &qs{qsp.second};

    std::vector<LongSimplexId> srcs{};
    std::vector<LongSimplexId> dsts{};

    findSepsVertices(qs, srcs, dsts);

    // remove duplicates
    std::set<LongSimplexId> srcs_set(srcs.begin(), srcs.end());
    std::set<LongSimplexId> dsts_set(dsts.begin(), dsts.end());
    srcs.assign(srcs_set.begin(), srcs_set.end());
    dsts.assign(dsts_set.begin(), dsts_set.end());

    // filter out separatrices whose sources are not in contact with
    // the current cell

    bool found = true;

    if(dsts.size() != 2) {
      found = false;
    }

    if(srcs.size() == 2) {
      outputCells_.emplace_back(Quad{dsts[0], srcs[0], dsts[1], srcs[1]});
    } else if(srcs.size() == 1) {
      outputCells_.emplace_back(Quad{dsts[0], srcs[0], dsts[1], srcs[0]});
      ndegen++;
    } else {
      found = false;
    }

    if(!found) {
      this->printMsg("Missing quadrangle", ttk::debug::Priority::DETAIL);
    }
  }

  return 0;
}

template <typename triangulationType>
size_t ttk::MorseSmaleQuadrangulation::findSeparatrixMiddle(
  const size_t a, const size_t b, const triangulationType &triangulation) {

  const int dim = 3;

  std::vector<float> distFromA(b - a + 1);
  std::array<float, dim> prev{}, curr{};

  if(distFromA.empty()) {
    return 0;
  }

  curr[0] = sepPoints_[dim * a];
  curr[1] = sepPoints_[dim * a + 1];
  curr[2] = sepPoints_[dim * a + 2];

  // integrate distances at every point of this separatrix
  for(size_t i = 1; i < b - a + 1; ++i) {
    std::swap(curr, prev);
    curr[0] = sepPoints_[dim * (a + i)];
    curr[1] = sepPoints_[dim * (a + i) + 1];
    curr[2] = sepPoints_[dim * (a + i) + 2];
    distFromA[i]
      = distFromA[i - 1] + ttk::Geometry::distance(&curr[0], &prev[0]);
  }

  auto distAB = distFromA.back();
  for(auto &el : distFromA) {
    el = std::abs(el - distAB / 2.0);
  }

  // index in separatrices point data array of separatrix middle
  auto pos = a + std::min_element(distFromA.begin(), distFromA.end())
             - distFromA.begin();

  // new point!
  outputPoints_.emplace_back(sepPoints_[dim * pos]);
  outputPoints_.emplace_back(sepPoints_[dim * pos + 1]);
  outputPoints_.emplace_back(sepPoints_[dim * pos + 2]);

  SimplexId id = pos;

  // new point identifier (on the triangular mesh)
  switch(sepCellDims_[pos]) {
    case 0:
      outputPointsIds_.emplace_back(sepCellIds_[pos]);
      break;
    case 1: {
      // take the first vertex of the edge
      triangulation.getEdgeVertex(sepCellIds_[pos], 0, id);
      outputPointsIds_.emplace_back(id);
      break;
    }
    case 2: {
      // take the first vertex of the triangle
      triangulation.getTriangleVertex(sepCellIds_[pos], 0, id);
      outputPointsIds_.emplace_back(id);
      break;
    }
    default:
      break;
  }

  outputPointsTypes_.emplace_back(1);

  return id;
}

template <typename triangulationType>
int ttk::MorseSmaleQuadrangulation::subdiviseDegenerateQuads(
  std::vector<Quad> &outputSubd, const triangulationType &triangulation) {

  for(size_t i = 0; i < outputCells_.size(); ++i) {
    auto q = outputCells_[i];
    auto seps = quadSeps_[i].second;

    // don't deal with normal quadrangles
    if(q[1] != q[3]) {
      continue;
    }

    std::vector<LongSimplexId> srcs{};
    std::vector<LongSimplexId> dsts{};

    findSepsVertices(seps, srcs, dsts);

    // identify the extremum that is twice dest
    int count_vi = 0, count_vk = 0;
    for(const auto &s : dsts) {
      if(s == q[0]) {
        count_vi++;
      }
      if(s == q[2]) {
        count_vk++;
      }
    }
    // extremum index
    LongSimplexId vert2Seps = count_vi > count_vk ? q[0] : q[2];
    LongSimplexId vert1Sep = count_vi > count_vk ? q[2] : q[0];
    // the two seps from j to vert2Seps
    std::vector<size_t> borderseps{};
    for(size_t j = 0; j < seps.size(); ++j) {
      if(dsts[j] == vert2Seps && srcs[j] == q[1]) {
        borderseps.emplace_back(seps[j]);
      }
    }

    if(borderseps.size() < 2) {
      continue;
    }

    // find a midpoint between the two extrema on the triangulation

    std::vector<SimplexId> boundi{criticalPointsIdentifier_[q[0]]};
    std::vector<SimplexId> boundk{criticalPointsIdentifier_[q[2]]};
    std::array<std::vector<float>, 6> outputDists{};

    Dijkstra::shortestPath(
      criticalPointsIdentifier_[q[0]], triangulation, outputDists[0], boundk);
    Dijkstra::shortestPath(
      criticalPointsIdentifier_[q[2]], triangulation, outputDists[1], boundi);

    auto inf = std::numeric_limits<float>::infinity();
    std::vector<float> sum(outputDists[0].size(), inf);

    for(size_t j = 0; j < sum.size(); ++j) {
      auto m = outputDists[0][j];
      auto n = outputDists[1][j];
      if(m == inf || n == inf) {
        continue;
      }
      // cost to minimize
      sum[j] = m + n + std::abs(m - n);
    }

    auto insertNewPoint
      = [&](const SimplexId a, const size_t idx, const SimplexId type) {
          float x, y, z;
          triangulation.getVertexPoint(a, x, y, z);
          outputPoints_.emplace_back(x);
          outputPoints_.emplace_back(y);
          outputPoints_.emplace_back(z);
          outputPointsIds_.emplace_back(a);
          outputPointsTypes_.emplace_back(type);
          outputPointsCells_.emplace_back(idx);
          return outputPointsIds_.size() - 1;
        };

    auto v0 = std::min_element(sum.begin(), sum.end()) - sum.begin();
    auto v0Pos = static_cast<LongSimplexId>(insertNewPoint(v0, i, 3));

    // find two other points

    std::vector<SimplexId> bounds{criticalPointsIdentifier_[q[0]],
                                  criticalPointsIdentifier_[q[1]],
                                  criticalPointsIdentifier_[q[2]]};

    auto m0Pos = sepMids_[borderseps[0]];
    auto m1Pos = sepMids_[borderseps[1]];
    auto m0 = outputPointsIds_[m0Pos];
    auto m1 = outputPointsIds_[m1Pos];

    Dijkstra::shortestPath(criticalPointsIdentifier_[vert1Sep], triangulation,
                           outputDists[2], bounds);
    Dijkstra::shortestPath(v0, triangulation, outputDists[3], bounds);
    Dijkstra::shortestPath(m0, triangulation, outputDists[4], bounds);
    Dijkstra::shortestPath(m1, triangulation, outputDists[5], bounds);

    std::fill(sum.begin(), sum.end(), inf);

    for(size_t j = 0; j < sum.size(); ++j) {
      auto m = outputDists[2][j];
      auto n = outputDists[3][j];
      auto o = outputDists[4][j];
      if(m == inf || n == inf || o == inf) {
        continue;
      }
      // cost to minimize
      sum[j] = m + n + o + std::abs(m - n) + std::abs(m - o) + std::abs(n - o);
    }

    auto v1 = std::min_element(sum.begin(), sum.end()) - sum.begin();
    auto v1Pos = static_cast<LongSimplexId>(insertNewPoint(v1, i, 4));
    std::fill(sum.begin(), sum.end(), inf);

    for(size_t j = 0; j < sum.size(); ++j) {
      auto m = outputDists[2][j];
      auto n = outputDists[3][j];
      auto o = outputDists[5][j];
      if(m == inf || n == inf || o == inf) {
        continue;
      }
      // cost to minimize
      sum[j] = m + n + o + std::abs(m - n) + std::abs(m - o) + std::abs(n - o);
    }

    auto v2 = std::min_element(sum.begin(), sum.end()) - sum.begin();
    auto v2Pos = static_cast<LongSimplexId>(insertNewPoint(v2, i, 4));

    outputSubd.emplace_back(Quad{vert2Seps, m0Pos, v1Pos, v0Pos});
    outputSubd.emplace_back(Quad{vert2Seps, m1Pos, v2Pos, v0Pos});
    outputSubd.emplace_back(Quad{q[1], m0Pos, v1Pos, vert1Sep});
    outputSubd.emplace_back(Quad{q[1], m1Pos, v2Pos, vert1Sep});
    outputSubd.emplace_back(Quad{vert1Sep, v1Pos, v0Pos, v2Pos});
  }
  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleQuadrangulation::subdivise(
  const triangulationType &triangulation) {

  // separatrices middles index in output points array
  sepMids_.resize(sepBegs_.size());

  for(size_t i = 0; i < sepMids_.size(); ++i) {
    // separatrices middles
    sepMids_[i] = outputPoints_.size() / 3; // before insertion at next line
    findSeparatrixMiddle(sepBegs_[i], sepEnds_[i], triangulation);
    outputPointsCells_.emplace_back(i);
  }

  // for each output quad, its barycenter position in outputPoints_
  std::vector<size_t> cellBary(outputCells_.size());

  std::array<std::vector<float>, 4> outputDists{};

  // hold quad subdivision
  decltype(outputCells_) outputSubd{};
  outputSubd.reserve(4 * outputCells_.size());

  for(size_t i = 0; i < outputCells_.size(); ++i) {
    auto q = outputCells_[i];
    auto seps = quadSeps_[i].second;

    // skip degenerate case here
    if(q[1] == q[3]) {
      continue;
    }

    std::vector<LongSimplexId> sepMids(seps.size());
    std::vector<SimplexId> midsNearestVertex(seps.size());
    for(size_t j = 0; j < seps.size(); ++j) {
      sepMids[j] = sepMids_[seps[j]];
      midsNearestVertex[j] = outputPointsIds_[sepMids_[seps[j]]];
    }

    // find barycenter of current cell (c.f. QuadrangulationSubdivision.cpp)

    // bound Dijkstra by parent quad vertices
    std::vector<SimplexId> bounds{
      criticalPointsIdentifier_[q[0]], criticalPointsIdentifier_[q[1]],
      criticalPointsIdentifier_[q[2]], criticalPointsIdentifier_[q[3]]};

    // Dijkstra propagation mask
    std::vector<bool> mask(morseSeg_.size(), false);
    // restrict Dijkstra propagation to current cell
    for(size_t j = 0; j < morseSeg_.size(); ++j) {
      if(morseSeg_[j] == quadSeps_[i].first) {
        mask[j] = true;
      }
    }
    // also allow to propagate on separatrices
    for(const auto s : seps) {
      for(size_t j = sepBegs_[s]; j <= sepEnds_[s]; ++j) {
        if(sepCellDims_[j] == 1) {
          auto e = sepCellIds_[j];
          SimplexId e0{}, e1{};
          triangulation.getEdgeVertex(e, 0, e0);
          triangulation.getEdgeVertex(e, 1, e1);
          mask[e0] = true;
          mask[e1] = true;
        }
      }
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t j = 0; j < outputDists.size(); ++j) {
      Dijkstra::shortestPath(midsNearestVertex[j], triangulation,
                             outputDists.at(j), std::vector<SimplexId>(), mask);
    }

    auto inf = std::numeric_limits<float>::infinity();
    std::vector<float> sum(outputDists[0].size(), inf);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t j = 0; j < sum.size(); ++j) {
      // skip if vertex j not in cell i
      if(morseSeg_[j] != quadSeps_[i].first) {
        continue;
      }
      auto m = outputDists[0][j];
      auto n = outputDists[1][j];
      auto o = outputDists[2][j];
      auto p = outputDists[3][j];
      if(m == inf || n == inf || o == inf || p == inf) {
        continue;
      }
      // cost to minimize
      sum[j] = m + n + o + p + std::abs(m - o) + std::abs(n - p);
    }

    size_t verticesInCell
      = verticesNumber_
        - std::count(
          sum.begin(), sum.end(), std::numeric_limits<float>::infinity());

    const size_t thresholdVertsInCell{50};
    if(verticesInCell <= thresholdVertsInCell) {
      this->printMsg("Small cell detected");
    }

    SimplexId baryId{};
    if(verticesInCell == 0) {
      this->printMsg("Barycenter of cell " + std::to_string(i) + " not found");
      // snap bary on sepMids[0]
      baryId = outputPointsIds_[sepMids[0]];
    } else {
      baryId = std::min_element(sum.begin(), sum.end()) - sum.begin();
    }

    LongSimplexId baryPos = outputPointsIds_.size();
    {
      float x, y, z;
      triangulation.getVertexPoint(baryId, x, y, z);
      outputPoints_.emplace_back(x);
      outputPoints_.emplace_back(y);
      outputPoints_.emplace_back(z);
      outputPointsIds_.emplace_back(baryId);
      outputPointsTypes_.emplace_back(2);
      outputPointsCells_.emplace_back(i);
    }

    std::vector<LongSimplexId> srcs{};
    std::vector<LongSimplexId> dsts{};

    findSepsVertices(seps, srcs, dsts);

    auto sepsQuadVertex = [&](const size_t a, const size_t b) {
      if(srcs[a] == srcs[b]) {
        outputSubd.emplace_back(Quad{srcs[a], sepMids[a], baryPos, sepMids[b]});
      }
      if(dsts[a] == dsts[b]) {
        outputSubd.emplace_back(Quad{dsts[a], sepMids[a], baryPos, sepMids[b]});
      }
    };

    // test every pair of current quad seps for a common vertex, if yes, use
    // their middles and the quad barycenter
    sepsQuadVertex(0, 1);
    sepsQuadVertex(0, 2);
    sepsQuadVertex(0, 3);
    sepsQuadVertex(1, 2);
    sepsQuadVertex(1, 3);
    sepsQuadVertex(2, 3);
  }

  subdiviseDegenerateQuads(outputSubd, triangulation);

  // overwrite old quads
  outputCells_ = std::move(outputSubd);

  return 0;
}

template <typename triangulationType>
bool ttk::MorseSmaleQuadrangulation::checkSurfaceCloseness(
  const triangulationType &triangulation) const {

  bool triangulationClosed{true};
  // sweep over all vertices to check if one is on a boundary
  for(SimplexId i = 0; i < verticesNumber_; ++i) {
    if(triangulation.isVertexOnBoundary(i)) {
      triangulationClosed = false;
      break;
    }
  }
  // quadrangles edges -> quadrangles
  std::map<std::pair<LongSimplexId, LongSimplexId>, std::set<size_t>>
    quadEdges{};

  // sweep over quadrangulation edges
  for(size_t i = 0; i < outputCells_.size(); ++i) {
    auto q = outputCells_[i];
    // store edges in order
    if(q[0] < q[1]) {
      quadEdges[std::make_pair(q[0], q[1])].emplace(i);
    } else {
      quadEdges[std::make_pair(q[1], q[0])].emplace(i);
    }
    if(q[1] < q[2]) {
      quadEdges[std::make_pair(q[1], q[2])].emplace(i);
    } else {
      quadEdges[std::make_pair(q[2], q[1])].emplace(i);
    }
    if(q[2] < q[3]) {
      quadEdges[std::make_pair(q[2], q[3])].emplace(i);
    } else {
      quadEdges[std::make_pair(q[3], q[2])].emplace(i);
    }
    if(q[3] < q[0]) {
      quadEdges[std::make_pair(q[3], q[0])].emplace(i);
    } else {
      quadEdges[std::make_pair(q[0], q[3])].emplace(i);
    }
  }

  bool quadrangulationClosed{true};
  for(const auto &e : quadEdges) {
    if(e.second.size() < 2) {
      quadrangulationClosed = false;
      break;
    }
  }

  // each edge should be shared by two different quadrangles
  return triangulationClosed == quadrangulationClosed;
}

// main routine
template <typename triangulationType>
int ttk::MorseSmaleQuadrangulation::execute(
  const triangulationType &triangulation) {

  Timer tm;

  // sanity check
  if(separatriceNumber_ == 0) {
    this->printErr("Unable to perform quadrangulation without separatrices");
    return 1;
  }

  // clear output
  clearData();
  outputPoints_.resize(3 * criticalPointsNumber_);
  outputPointsIds_.resize(criticalPointsNumber_);
  outputPointsTypes_.resize(criticalPointsNumber_);
  outputPointsCells_.resize(criticalPointsNumber_);

  // fill in critical points 3d coordinates and identifiers
  for(SimplexId i = 0; i < criticalPointsNumber_; ++i) {
    outputPoints_[3 * i] = criticalPoints_[3 * i];
    outputPoints_[3 * i + 1] = criticalPoints_[3 * i + 1];
    outputPoints_[3 * i + 2] = criticalPoints_[3 * i + 2];
    outputPointsIds_[i] = criticalPointsIdentifier_[i];
    outputPointsTypes_[i] = 0;
    outputPointsCells_[i] = i;
  }

  // number of degenerate quadrangles
  size_t ndegen = 0;

  // direct quadrangulation with saddle points
  int ret = quadrangulate(ndegen, triangulation);

  if(ret == 0) {
    subdivise(triangulation);
  } else {
    // clean, log & early return
    clearData();
    this->printErr("Unable to generate a quadrangulation from the given "
                   "Morse-Smale complex");
    return 1;
  }

  if(DualQuadrangulation) {
    dualQuadrangulate();
  }

  if(!checkSurfaceCloseness(triangulation)) {
    // log, clean & early return
    this->printErr("Output surface does not match input surface closeness");
    if(!ShowResError) {
      clearData();
      return 1;
    }
  }

  this->printMsg("Produced " + std::to_string(outputCells_.size()) + " ("
                   + std::to_string(ndegen) + " degenerated)",
                 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 0;
}
