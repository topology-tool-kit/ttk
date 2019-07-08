#include <BarycentricSubdivision.h>
#include <Dijkstra.h>
#include <Geometry.h>
#include <SurfaceQuadrangulation.h>

#include <array>
#include <cmath>
#include <numeric>
#include <queue>

#define MODULE_S "[SurfaceQuadrangulation] "

int ttk::SurfaceQuadrangulation::dualQuadrangulate() {

  // quadrangles vertices are only extrema

  // filter sepCellIds_ array according to sepMask_
  std::vector<SimplexId> sepFlatEdges{};

  for(SimplexId i = 0; i < separatriceNumber_; ++i) {
    if(sepMask_[i] == 1) {
      continue;
    }
    sepFlatEdges.emplace_back(sepCellIds_[i]);
  }

  if(sepFlatEdges.size() % 2 != 0) {
    std::stringstream msg;
    msg << MODULE_S "Error: odd number of separatrices edges" << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
    return -1;
  }

  // holds separatrices edges for every separatrix
  std::vector<std::pair<SimplexId, SimplexId>> sepEdges{};

  for(size_t i = 0; i < sepFlatEdges.size() / 2; ++i) {
    sepEdges.emplace_back(
      std::make_pair(sepFlatEdges[2 * i], sepFlatEdges[2 * i + 1]));
  }

  // maps sources (saddle points) to vector of their destinations (extrema)
  std::map<SimplexId, std::vector<SimplexId>> sourceDests{};

  for(auto &p : sepEdges) {
    SimplexId i;
    for(i = 0; i < criticalPointsNumber_; i++) {
      if(p.first == criticalPointsCellIds_[i]) {
        break;
      }
    }
    SimplexId j;
    for(j = 0; j < criticalPointsNumber_; j++) {
      if(p.second == criticalPointsCellIds_[j]) {
        break;
      }
    }

    auto &v = sourceDests[i];
    v.emplace_back(j);
  }

  for(auto &elt : sourceDests) {
    auto extrema = elt.second;
    if(extrema.size() == 4) {
      auto i = extrema[0];
      auto j = i;
      auto k = i;
      auto l = i;
      // filter extrema by nature (minimum: 0 or maximum: 2)
      if(criticalPointsType_[extrema[1]] == criticalPointsType_[i]) {
        j = extrema[2];
        k = extrema[1];
        l = extrema[3];
      } else if(criticalPointsType_[extrema[2]] == criticalPointsType_[i]) {
        j = extrema[1];
        k = extrema[2];
        l = extrema[3];
      } else if(criticalPointsType_[extrema[3]] == criticalPointsType_[i]) {
        j = extrema[2];
        k = extrema[3];
        l = extrema[1];
      }
      outputCells_.emplace_back(4);
      outputCells_.emplace_back(i);
      outputCells_.emplace_back(j);
      outputCells_.emplace_back(k);
      outputCells_.emplace_back(l);
    }
  }

  return 0;
}

// ad-hoc quad data structure (see QuadrangulationSubdivision.h)
struct Quad {
  long long n;
  long long i;
  long long j;
  long long k;
  long long l;
};

int ttk::SurfaceQuadrangulation::detectCellSeps() {
  BarycentricSubdivision bs{};
  Triangulation newT{};

  bs.setupTriangulation(triangulation_);
  bs.setOutputTriangulation(&newT);
  bs.setInputPoints(inputPoints_);
  bs.execute();

  newT.preprocessVertexNeighbors();
  newT.preprocessVertexEdges();
  newT.preprocessVertexTriangles();
  newT.preprocessEdgeTriangles();
  newT.preprocessTriangleEdges();

  auto nVerts = triangulation_->getNumberOfVertices();
  auto nEdges = triangulation_->getNumberOfEdges();

  // store separatrix index on subdivised triangulation vertices
  std::vector<SimplexId> onSep(newT.getNumberOfVertices(), -1);
  // store separatrix index on subdivised triangulation edges
  std::vector<SimplexId> edgeOnSep(newT.getNumberOfEdges(), -1);
  // count the number of critical points encountered
  size_t critPoints{0};
  if(sepMask_[0] == 0) {
    critPoints++;
  }

  // id of critical point in new triangulation
  auto critPointId = [&](const SimplexId a) {
    if(sepCellDims_[a] == 0) {
      return sepCellIds_[a];
    }
    if(sepCellDims_[a] == 1) {
      return nVerts + sepCellIds_[a];
    }
    if(sepCellDims_[a] == 2) {
      return nVerts + nEdges + sepCellIds_[a];
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
    }

    // current separatrix id is critPoints // 2
    auto currSepId = critPoints / 2;

    // current point
    SimplexId curr{critPointId(i)};

    // update value on vertex TODO
    onSep[curr] = currSepId;

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

  // store saddle points id in the subdivised triangulation
  std::set<SimplexId> saddlesId{};

  for(SimplexId i = 0; i < criticalPointsNumber_; ++i) {
    // keep only saddle points
    if(criticalPointsType_[i] != 1) {
      continue;
    }
    saddlesId.emplace(nVerts + criticalPointsCellIds_[i]);
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
      sepId.emplace(edgeOnSep[e0]);
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
    return (saddlesId.find(a) != saddlesId.end()
            || saddlesId.find(b) != saddlesId.end()
            || saddlesId.find(c) != saddlesId.end());
  };

  // look around the saddle points
  for(SimplexId i = 0; i < criticalPointsNumber_; ++i) {
    // keep only saddle points
    if(criticalPointsType_[i] != 1) {
      continue;
    }
    SimplexId saddle = nVerts + criticalPointsCellIds_[i];

    // propagate from triangles around saddle
    std::vector<bool> processed(newT.getNumberOfTriangles(), false);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId j = 0; j < newT.getVertexTriangleNumber(saddle); ++j) {
      std::queue<SimplexId> toProcess{};
      SimplexId tr;
      newT.getVertexTriangle(saddle, j, tr);

      std::set<SimplexId> sepIdBeg = sepIdAroundTriangle(tr);

      toProcess.push(tr);

      while(!toProcess.empty()) {
        auto curr = toProcess.front();
        toProcess.pop();

        // skip already processed
        if(processed[curr]) {
          continue;
        }

        // check for saddle at vertices
        bool hasSaddle = hasTriangleSaddle(curr);

        // check for separatrices on edges
        auto sepIdEnd = sepIdAroundTriangle(curr);

        processed[curr] = true;

        if(hasSaddle && sepIdEnd.size() == 2) {
          // check that seps are different from beginning
          std::vector<size_t> cellSeps{};
          std::set_union(sepIdBeg.begin(), sepIdBeg.end(), sepIdEnd.begin(),
                         sepIdEnd.end(), std::back_inserter(cellSeps));
          if(cellSeps.size() > 2) {
            // found it
            quadSeps_.emplace_back(cellSeps);
            break;
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
            if(!processed[neigh]) {
              toProcess.push(neigh);
            }
          }
        }
      }
    }
  }
  return 0;
}

int ttk::SurfaceQuadrangulation::quadrangulate(size_t &ndegen) {
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
  quadSeps_.clear();

  // fill in data arrays
  for(size_t i = 0; i < numSeps; ++i) {
    // separatrices bounds
    sepBegs_[i] = sepFlatEdges[2 * i];
    sepEnds_[i] = sepFlatEdges[2 * i + 1];
  }

  detectCellSeps();

  outputCells_.reserve(5 * quadSeps_.size());
  auto quads = reinterpret_cast<std::vector<Quad> *>(&outputCells_);

  for(const auto &qs : quadSeps_) {

    std::vector<long long> srcs{};
    std::vector<long long> dsts{};

    findSepsVertices(qs, srcs, dsts);

    // remove duplicates
    std::set<long long> srcs_set(srcs.begin(), srcs.end());
    std::set<long long> dsts_set(dsts.begin(), dsts.end());
    srcs.assign(srcs_set.begin(), srcs_set.end());
    dsts.assign(dsts_set.begin(), dsts_set.end());

    // filter out separatrices whose sources are not in contact with
    // the current cell

    bool found = true;

    if(dsts.size() != 2) {
      found = false;
    }

    if(srcs.size() == 2) {
      quads->emplace_back(Quad{4, dsts[0], srcs[0], dsts[1], srcs[1]});
    } else if(srcs.size() == 1) {
      quads->emplace_back(Quad{4, dsts[0], srcs[0], dsts[1], srcs[0]});
      ndegen++;
    } else {
      found = false;
    }

    if(!found) {
      std::stringstream msg;
      msg << MODULE_S "Missing quadrangle" << std::endl;
      dMsg(std::cout, msg.str(), detailedInfoMsg);
      return 1;
    }
  }

  return 0;
}

size_t ttk::SurfaceQuadrangulation::findSeparatrixMiddle(const size_t a,
                                                         const size_t b) {

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
      triangulation_->getEdgeVertex(sepCellIds_[pos], 0, id);
      outputPointsIds_.emplace_back(id);
      break;
    }
    case 2: {
      // take the first vertex of the triangle
      triangulation_->getTriangleVertex(sepCellIds_[pos], 0, id);
      outputPointsIds_.emplace_back(id);
      break;
    }
    default:
      break;
  }

  outputPointsTypes_.emplace_back(1);

  return id;
}

int ttk::SurfaceQuadrangulation::findSepsVertices(
  const std::vector<size_t> &seps,
  std::vector<long long> &srcs,
  std::vector<long long> &dsts) const {

  srcs.resize(seps.size());
  dsts.resize(seps.size());

  for(size_t i = 0; i < seps.size(); ++i) {
    auto src = sepCellIds_[sepBegs_[seps[i]]];
    auto dst = sepCellIds_[sepEnds_[seps[i]]];
    for(long long j = 0; j < criticalPointsNumber_; ++j) {
      if(criticalPointsCellIds_[j] == src) {
        srcs[i] = j;
      }
      if(criticalPointsCellIds_[j] == dst) {
        dsts[i] = j;
      }
    }
  }

  return 0;
}

std::vector<long long> ttk::SurfaceQuadrangulation::subdiviseDegenerateQuads() {
  // hold quad subdivision
  decltype(outputCells_) outputSubd{};
  outputSubd.reserve(4 * outputCells_.size());
  auto quads = reinterpret_cast<std::vector<Quad> *>(&outputCells_);
  auto qsubd = reinterpret_cast<std::vector<Quad> *>(&outputSubd);

  for(size_t i = 0; i < quads->size(); ++i) {
    auto q = quads->at(i);
    auto seps = quadSeps_[i];

    // don't deal with normal quadrangles
    if(q.j != q.l) {
      continue;
    }

    std::vector<long long> srcs{};
    std::vector<long long> dsts{};

    findSepsVertices(seps, srcs, dsts);

    // identify the extremum that is twice dest
    int count_vi = 0, count_vk = 0;
    for(const auto &s : dsts) {
      if(s == q.i) {
        count_vi++;
      }
      if(s == q.k) {
        count_vk++;
      }
    }
    // extremum index
    long long vert2Seps = count_vi > count_vk ? q.i : q.k;
    long long vert1Sep = count_vi > count_vk ? q.k : q.i;
    // the two seps from j to vert2Seps
    std::vector<size_t> borderseps{};
    for(size_t j = 0; j < seps.size(); ++j) {
      if(dsts[j] == vert2Seps && srcs[j] == q.j) {
        borderseps.emplace_back(seps[j]);
      }
    }

    if(borderseps.size() < 2) {
      continue;
    }

    // find a midpoint between the two extrema on the triangulation

    std::vector<SimplexId> boundi{criticalPointsIdentifier_[q.i]};
    std::vector<SimplexId> boundk{criticalPointsIdentifier_[q.k]};
    std::array<std::vector<float>, 6> outputDists{};

    Dijkstra::shortestPath(
      criticalPointsIdentifier_[q.i], *triangulation_, outputDists[0], boundk);
    Dijkstra::shortestPath(
      criticalPointsIdentifier_[q.k], *triangulation_, outputDists[1], boundi);

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
          triangulation_->getVertexPoint(a, x, y, z);
          outputPoints_.emplace_back(x);
          outputPoints_.emplace_back(y);
          outputPoints_.emplace_back(z);
          outputPointsIds_.emplace_back(a);
          outputPointsTypes_.emplace_back(type);
          outputPointsCells_.emplace_back(idx);
          return outputPointsIds_.size() - 1;
        };

    auto v0 = std::min_element(sum.begin(), sum.end()) - sum.begin();
    auto v0Pos = static_cast<long long>(insertNewPoint(v0, i, 3));

    // find two other points

    std::vector<SimplexId> bounds{criticalPointsIdentifier_[q.i],
                                  criticalPointsIdentifier_[q.j],
                                  criticalPointsIdentifier_[q.k]};

    auto m0Pos = sepMids_[borderseps[0]];
    auto m1Pos = sepMids_[borderseps[1]];
    auto m0 = outputPointsIds_[m0Pos];
    auto m1 = outputPointsIds_[m1Pos];

    Dijkstra::shortestPath(criticalPointsIdentifier_[vert1Sep], *triangulation_,
                           outputDists[2], bounds);
    Dijkstra::shortestPath(v0, *triangulation_, outputDists[3], bounds);
    Dijkstra::shortestPath(m0, *triangulation_, outputDists[4], bounds);
    Dijkstra::shortestPath(m1, *triangulation_, outputDists[5], bounds);

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
    auto v1Pos = static_cast<long long>(insertNewPoint(v1, i, 4));
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
    auto v2Pos = static_cast<long long>(insertNewPoint(v2, i, 4));

    qsubd->emplace_back(Quad{4, vert2Seps, m0Pos, v1Pos, v0Pos});
    qsubd->emplace_back(Quad{4, vert2Seps, m1Pos, v2Pos, v0Pos});
    qsubd->emplace_back(Quad{4, q.j, m0Pos, v1Pos, vert1Sep});
    qsubd->emplace_back(Quad{4, q.j, m1Pos, v2Pos, vert1Sep});
    qsubd->emplace_back(Quad{4, vert1Sep, v1Pos, v0Pos, v2Pos});
  }
  return outputSubd;
}

int ttk::SurfaceQuadrangulation::subdivise() {

  // separatrices middles index in output points array
  sepMids_.resize(sepBegs_.size());

  for(size_t i = 0; i < sepMids_.size(); ++i) {
    // separatrices middles
    sepMids_[i] = outputPoints_.size() / 3; // before insertion at next line
    findSeparatrixMiddle(sepBegs_[i], sepEnds_[i]);
    outputPointsCells_.emplace_back(i);
  }

  // for each output quad, its barycenter position in outputPoints_
  std::vector<size_t> cellBary(outputCells_.size());

  std::array<std::vector<float>, 4> outputDists{};

  // hold quad subdivision
  decltype(outputCells_) outputSubd{};
  outputSubd.reserve(4 * outputCells_.size());
  auto quads = reinterpret_cast<std::vector<Quad> *>(&outputCells_);
  auto qsubd = reinterpret_cast<std::vector<Quad> *>(&outputSubd);

  for(size_t i = 0; i < quads->size(); ++i) {
    auto q = quads->at(i);
    auto seps = quadSeps_[i];

    // skip degenerate case here
    if(q.j == q.l) {
      continue;
    }

    std::vector<long long> sepMids(seps.size());
    std::vector<SimplexId> midsNearestVertex(seps.size());
    for(size_t j = 0; j < seps.size(); ++j) {
      sepMids[j] = sepMids_[seps[j]];
      midsNearestVertex[j] = outputPointsIds_[sepMids_[seps[j]]];
    }

    // find barycenter of current cell (c.f. QuadrangulationSubdivision.cpp)

    // bound Dijkstra by parent quad vertices
    std::vector<SimplexId> bounds{
      criticalPointsIdentifier_[q.i], criticalPointsIdentifier_[q.j],
      criticalPointsIdentifier_[q.k], criticalPointsIdentifier_[q.l]};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t j = 0; j < outputDists.size(); ++j) {
      Dijkstra::shortestPath(
        midsNearestVertex[j], *triangulation_, outputDists.at(j), bounds);
    }

    auto inf = std::numeric_limits<float>::infinity();
    std::vector<float> sum(outputDists[0].size(), inf);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t j = 0; j < sum.size(); ++j) {
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
      = segmentationNumber_
        - std::count(
          sum.begin(), sum.end(), std::numeric_limits<float>::infinity());

    if(verticesInCell == 0) {
      std::stringstream msg;
      msg << MODULE_S "Barycenter in cell " << i << " not found." << std::endl;
      dMsg(std::cout, msg.str(), infoMsg);
    }

    const size_t thresholdVertsInCell{50};
    if(verticesInCell <= thresholdVertsInCell) {
      std::stringstream msg;
      msg << MODULE_S "Small cell detected" << std::endl;
      dMsg(std::cout, msg.str(), infoMsg);
    }

    auto baryId = std::min_element(sum.begin(), sum.end()) - sum.begin();
    long long baryPos = outputPointsIds_.size();
    {
      float x, y, z;
      triangulation_->getVertexPoint(baryId, x, y, z);
      outputPoints_.emplace_back(x);
      outputPoints_.emplace_back(y);
      outputPoints_.emplace_back(z);
      outputPointsIds_.emplace_back(baryId);
      outputPointsTypes_.emplace_back(2);
      outputPointsCells_.emplace_back(i);
    }

    std::vector<long long> srcs{};
    std::vector<long long> dsts{};

    findSepsVertices(seps, srcs, dsts);

    auto sepsQuadVertex = [&](const size_t a, const size_t b) {
      if(srcs[a] == srcs[b]) {
        qsubd->emplace_back(Quad{4, srcs[a], sepMids[a], baryPos, sepMids[b]});
      }
      if(dsts[a] == dsts[b]) {
        qsubd->emplace_back(Quad{4, dsts[a], sepMids[a], baryPos, sepMids[b]});
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

  auto degenSubd = subdiviseDegenerateQuads();
  std::copy(degenSubd.begin(), degenSubd.end(), std::back_inserter(outputSubd));

  // overwrite old quads
  outputCells_ = std::move(outputSubd);

  return 0;
}

bool ttk::SurfaceQuadrangulation::checkSurfaceCloseness() const {
  bool triangulationClosed{true};
  // sweep over all vertices to check if one is on a boundary
  for(size_t i = 0; i < segmentationNumber_; ++i) {
    if(triangulation_->isVertexOnBoundary(i)) {
      triangulationClosed = false;
      break;
    }
  }
  // quadrangles edges -> quadrangles
  std::map<std::pair<long long, long long>, std::set<size_t>> quadEdges{};

  // sweep over quadrangulation edges
  auto quads = reinterpret_cast<const std::vector<Quad> *>(&outputCells_);
  for(size_t i = 0; i < quads->size(); ++i) {
    auto q = quads->at(i);
    // store edges in order
    if(q.i < q.j) {
      quadEdges[std::make_pair(q.i, q.j)].emplace(i);
    } else {
      quadEdges[std::make_pair(q.j, q.i)].emplace(i);
    }
    if(q.j < q.k) {
      quadEdges[std::make_pair(q.j, q.k)].emplace(i);
    } else {
      quadEdges[std::make_pair(q.k, q.j)].emplace(i);
    }
    if(q.k < q.l) {
      quadEdges[std::make_pair(q.k, q.l)].emplace(i);
    } else {
      quadEdges[std::make_pair(q.l, q.k)].emplace(i);
    }
    if(q.l < q.i) {
      quadEdges[std::make_pair(q.l, q.i)].emplace(i);
    } else {
      quadEdges[std::make_pair(q.i, q.l)].emplace(i);
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

void ttk::SurfaceQuadrangulation::clearData() {
  outputCells_.clear();
  outputPoints_.clear();
  outputPointsIds_.clear();
  outputPointsCells_.clear();
}

// main routine
int ttk::SurfaceQuadrangulation::execute() {

  Timer t;

  // sanity check
  if(separatriceNumber_ == 0) {
    std::stringstream msg;
    msg << MODULE_S "Error: cannot perform quadrangulation without separatrices"
        << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
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

  if(dualQuadrangulation_) {
    dualQuadrangulate();
  } else {
    // direct quadrangulation with saddle points
    int ret = quadrangulate(ndegen);

    if(ret == 0) {
      subdivise();
    } else {
      // clean, log & early return
      clearData();
      std::stringstream msg;
      msg << MODULE_S "Error: unable to generate quadrangulation from current "
                      "Morse-Smale complex"
          << std::endl;
      dMsg(std::cout, msg.str(), infoMsg);
      return 1;
    }
  }

  if(!checkSurfaceCloseness()) {
    // clean, log & early return
    clearData();
    std::stringstream msg;
    msg << MODULE_S
      "Error: output surface does not match input surface closeness"
        << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
    return 1;
  }

  // number of produced quads
  size_t quadNumber = outputCells_.size() / 5;

  {
    std::stringstream msg;
    msg << MODULE_S "Produced " << quadNumber << " quadrangles after "
        << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  {
    std::stringstream msg;
    msg << MODULE_S << quadNumber << " quads (" << ndegen << " degenerated)"
        << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}
