#include <Dijkstra.h>
#include <Geometry.h>
#include <SurfaceQuadrangulation.h>
#include <array>
#include <cmath>
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

size_t ttk::SurfaceQuadrangulation::sepFromPoints(const long long src,
                                                  const long long dst) const {
  for(size_t i = 0; i < sepBegs_.size(); ++i) {
    if(sepCellIds_[sepBegs_[i]] == criticalPointsCellIds_[src]
       && sepCellIds_[sepEnds_[i]] == criticalPointsCellIds_[dst]) {
      return i;
    }
  }
  return sepBegs_.size();
}

int ttk::SurfaceQuadrangulation::detectCells(
  const SimplexId src,
  std::vector<SimplexId> &vertexCells,
  std::vector<std::vector<SimplexId>> &cellSeps,
  const std::vector<SimplexId> &vertexSepMask) {

  std::vector<bool> cellMask(vertexCells.size(), false);

  std::queue<SimplexId> toProcess{};
  toProcess.push(src);
  std::vector<SimplexId> borderSeps{};

  auto cellId = segmentation_[src];
  SimplexId newCellId = cellSeps.size();

  auto isCandidate = [&](const SimplexId a) {
    return segmentation_[a] == cellId && vertexSepMask[a] == -1 && !cellMask[a];
  };

  while(!toProcess.empty()) {
    auto curr = toProcess.front();
    toProcess.pop();

    if(!isCandidate(curr)) {
      continue;
    }

    cellMask[curr] = true;

    auto nneigh = triangulation_->getVertexNeighborNumber(curr);
    for(SimplexId j = 0; j < nneigh; ++j) {
      SimplexId next;
      triangulation_->getVertexNeighbor(curr, j, next);

      if(isCandidate(next)) {
        toProcess.push(next);
      }

      // store reached separatrices indices
      if(vertexSepMask[next] != -1) {
        borderSeps.emplace_back(vertexSepMask[next]);
      }
    }
  }

  // histogram of border separatrices indices
  std::map<SimplexId, int> hist{};

  // post-process borderSeps to find the most common separatrices indices
  for(const auto &v : borderSeps) {
    hist[v]++;
  }

  // filter out small cells
  if(hist.size() < 3) {
    return 1;
  }

  // copy sub-segmentation into output
  for(size_t i = 0; i < vertexCells.size(); ++i) {
    if(cellMask[i]) {
      vertexCells[i] = newCellId;
    }
  }

  // map dumped into vector
  std::vector<std::pair<SimplexId, int>> histVec{};
  histVec.reserve(hist.size());

  for(const auto &p : hist) {
    histVec.emplace_back(p);
  }

  // sort by value, descending order
  std::sort(
    histVec.begin(), histVec.end(),
    [&](const std::pair<SimplexId, int> &a,
        const std::pair<SimplexId, int> &b) { return a.second > b.second; });

  std::vector<SimplexId> sepIds(histVec.size());
  for(size_t i = 0; i < histVec.size(); ++i) {
    sepIds[i] = histVec[i].first;
  }

  // return all reached separatrices by importance order
  cellSeps.emplace_back(sepIds);
  // link this cell to MorseSmale Manifold index
  cellId_.emplace_back(segmentation_[src]);

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
  morseSeg_.resize(segmentationNumber_);
  std::fill(morseSeg_.begin(), morseSeg_.end(), -1);
  cellId_.clear();
  quadSeps_.clear();

  // fill in data arrays
  for(size_t i = 0; i < numSeps; ++i) {
    // separatrices bounds
    sepBegs_[i] = sepFlatEdges[2 * i];
    sepEnds_[i] = sepFlatEdges[2 * i + 1];
  }

  // for each vertex on a separatrix, the index of the separatrix
  std::vector<SimplexId> onSep(segmentationNumber_, -1);

  {
    size_t currSep = 0;
    // iterate over separatrices lines to fill in onSep vector
    for(SimplexId i = 0; i < separatriceNumber_; ++i) {
      if(sepMask_[i] == 0) {
        currSep
          = std::find(sepBegs_.begin(), sepBegs_.end(), i) - sepBegs_.begin();
      }
      if(sepCellDims_[i] == 1) {
        SimplexId vertId;
        triangulation_->getEdgeVertex(sepCellIds_[i], 0, vertId);
        onSep[vertId] = static_cast<SimplexId>(currSep);
        triangulation_->getEdgeVertex(sepCellIds_[i], 1, vertId);
        onSep[vertId] = static_cast<SimplexId>(currSep);
      }
    }
  }

  // for each cell, the indices of the bordering separatrices
  std::vector<std::vector<SimplexId>> cellSeps{};

  bool finished = false;
  size_t pos = 0;

  while(true) {
    for(size_t j = pos; j < segmentationNumber_; ++j) {
      if(onSep[j] == -1 && morseSeg_[j] == -1) {
        pos = j;
        break;
      }
      if(j == segmentationNumber_ - 1) {
        finished = true;
      }
    }
    if(finished) {
      break;
    }
    detectCells(pos, morseSeg_, cellSeps, onSep);
    pos++;
  }

  // missing cells?
  // find them and store their bordering separatrices
  auto minCellId
    = *std::min_element(segmentation_, segmentation_ + segmentationNumber_);
  auto maxCellId
    = *std::max_element(segmentation_, segmentation_ + segmentationNumber_);
  for(SimplexId i = minCellId; i <= maxCellId; ++i) {
    if(std::find(cellId_.begin(), cellId_.end(), i) == cellId_.end()) {
      std::set<SimplexId> sepIds{};
      for(size_t j = 0; j < segmentationNumber_; ++j) {
        if(segmentation_[j] == i && onSep[j] != -1) {
          sepIds.emplace(onSep[j]);
        }
      }
      if(sepIds.size() > 2) {
        cellSeps.emplace_back(sepIds.begin(), sepIds.end());
        cellId_.emplace_back(i);
      }
    }
  }

  // hold quad subdivision
  outputCells_.reserve(5 * cellSeps.size());
  auto quads = reinterpret_cast<std::vector<Quad> *>(&outputCells_);

  // for each cell the detected separatrices sources and destinations
  std::vector<std::vector<long long>> cellSrcs{}, cellDsts{};

  for(const auto &c : cellSeps) {

    std::vector<long long> srcs(c.size());
    std::vector<long long> dsts(c.size());

    for(size_t i = 0; i < c.size(); ++i) {
      auto src = sepCellIds_[sepBegs_[c[i]]];
      auto dst = sepCellIds_[sepEnds_[c[i]]];
      for(long long j = 0; j < criticalPointsNumber_; ++j) {
        if(criticalPointsCellIds_[j] == src) {
          srcs[i] = j;
        }
        if(criticalPointsCellIds_[j] == dst) {
          dsts[i] = j;
        }
      }
    }
    cellSrcs.emplace_back(srcs);
    cellDsts.emplace_back(dsts);

    bool found = false;

    // optimisation: try first to only consider the first four separatrices

    if(c.size() >= 4) {
      std::set<SimplexId> srcs_set(srcs.begin(), std::next(srcs.begin(), 4));
      std::set<SimplexId> dsts_set(dsts.begin(), std::next(dsts.begin(), 4));
      if(srcs_set.size() == 2 && dsts_set.size() == 2) {
        auto vi = *dsts_set.begin();
        auto vj = *srcs_set.begin();
        auto vk = *std::next(dsts_set.begin());
        auto vl = *std::next(srcs_set.begin());
        quads->emplace_back(Quad{4, vi, vj, vk, vl});
        found = true;
        continue;
      }
    }

    // normal case: find a pair of extrema and a pair of saddle points
    // with four separatrices linking one another

    // iterate over first dest
    for(size_t i = 0; i < c.size(); ++i) {
      auto vi = dsts[i];

      // iterate over second dest from i + 1
      for(size_t k = i + 1; k < c.size(); ++k) {
        auto vk = dsts[k];
        // skip same extrema type
        if(criticalPointsType_[vi] == criticalPointsType_[vk]) {
          continue;
        }
        // at least four separatrices leading to these two extrema
        if(std::count(dsts.begin(), dsts.end(), vi)
             + std::count(dsts.begin(), dsts.end(), vk)
           < 4) {
          continue;
        }
        // find two (one if degenerate) common sources
        std::set<SimplexId> srcs_i{};
        std::set<SimplexId> srcs_k{};
        std::vector<SimplexId> common_srcs_ik{};
        for(size_t j = 0; j < c.size(); ++j) {
          if(dsts[j] == vi) {
            srcs_i.insert(srcs[j]);
          }
          if(dsts[j] == vk) {
            srcs_k.insert(srcs[j]);
          }
        }
        std::set_intersection(srcs_i.begin(), srcs_i.end(), srcs_k.begin(),
                              srcs_k.end(), std::back_inserter(common_srcs_ik));
        // reorder common_srcs_ik according to srcs vector
        std::vector<SimplexId> common_srcs_ik_ordered{};
        for(const auto &s : srcs) {
          if(std::find(common_srcs_ik.begin(), common_srcs_ik.end(), s)
               != common_srcs_ik.end()
             && std::find(common_srcs_ik_ordered.begin(),
                          common_srcs_ik_ordered.end(), s)
                  == common_srcs_ik_ordered.end()) {
            common_srcs_ik_ordered.emplace_back(s);
          }
        }
        if(common_srcs_ik_ordered.size() > 1) {
          auto vj = common_srcs_ik_ordered[0];
          auto vl = common_srcs_ik_ordered[1];
          quads->emplace_back(Quad{4, vi, vj, vk, vl});
          found = true;
          break;
        }
      }
      if(found) {
        // stop at the first correct quad
        break;
      }
    }

    if(!found) {
      // look at the three first separatrices, try to find a missing fourth
      if(c.size() >= 3) {
        auto vi = dsts[0];
        auto vj = srcs[0];
        decltype(vi) vk{};
        decltype(vj) vl{};
        for(size_t k = 1; k < 3; ++k) {
          if(criticalPointsCellIds_[dsts[k]] != criticalPointsCellIds_[vi]) {
            vk = dsts[k];
            break;
          }
        }
        for(size_t l = 1; l < 3; ++l) {
          if(srcs[l] != vj) {
            vl = srcs[l];
            break;
          }
        }
        // ensure jk, kl and li separatrices exist
        if(sepFromPoints(vj, vk) < numSeps && sepFromPoints(vl, vi) < numSeps
           && sepFromPoints(vl, vk) < numSeps) {
          quads->emplace_back(Quad{4, vi, vj, vk, vl});
          found = true;
          continue;
        }
      }

      // degenerate case:
      // take the first two distinct extrema (separatrices with higher weight)
      for(size_t i = 0; i < c.size(); ++i) {
        auto vi = dsts[0];
        for(size_t k = i + 1; k < dsts.size(); ++k) {
          auto vk = dsts[k];
          // skip same critical point type
          if(criticalPointsCellIds_[vi] == criticalPointsCellIds_[vk]) {
            continue;
          }
          // we need at least 3 separatrices leading to these extrema
          if(std::count(dsts.begin(), dsts.end(), vi)
               + std::count(dsts.begin(), dsts.end(), vk)
             < 3) {
            continue;
          }
          // count saddle point occurences
          int count_vj = 0;
          for(size_t j = i; j < srcs.size(); ++j) {
            if(srcs[j] == srcs[i] && (dsts[j] == vi || dsts[j] == vk)) {
              count_vj++;
            }
          }
          // one saddle point for three separatrices
          if(count_vj < 3) {
            continue;
          }
          auto vj = srcs[i];
          found = true;
          quads->emplace_back(Quad{4, vi, vj, vk, vj});
          ndegen++;
        }
        if(found) {
          break;
        }
      }
    }
    if(!found) {
      std::stringstream msg;
      msg << MODULE_S "Missing quadrangle" << std::endl;
      dMsg(std::cout, msg.str(), detailedInfoMsg);
    }
  }

  // find & order separatrices bordering quadrangles
  for(size_t i = 0; i < quads->size(); ++i) {
    auto c = cellSeps[i];
    auto srcs = cellSrcs[i];
    auto dsts = cellDsts[i];
    auto q = quads->at(i);
    std::vector<size_t> seps{};

    auto findSep = [&](const long long src, const long long dst) {
      for(size_t j = 0; j < c.size(); ++j) {
        if(srcs[j] == src && dsts[j] == dst) {
          seps.emplace_back(c[j]);
          return;
        }
      }
      // if not found in cell seps, try first one from src to dst
      seps.emplace_back(sepFromPoints(src, dst));
    };

    findSep(q.j, q.i);
    findSep(q.j, q.k);
    findSep(q.l, q.k);
    findSep(q.l, q.i);

    quadSeps_.emplace_back(seps);
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

int ttk::SurfaceQuadrangulation::subdivise() {

  // separatrices middles index in output points array
  std::vector<SimplexId> sepMiddle(sepBegs_.size());

  for(size_t i = 0; i < sepMiddle.size(); ++i) {
    // separatrices middles
    sepMiddle[i] = outputPoints_.size() / 3; // before insertion at next line
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

    std::vector<long long> sepMids(seps.size());
    std::vector<SimplexId> midsNearestVertex(seps.size());
    for(size_t j = 0; j < seps.size(); ++j) {
      sepMids[j] = sepMiddle[seps[j]];
      midsNearestVertex[j] = outputPointsIds_[sepMiddle[seps[j]]];
    }

    // find barycenter of current cell (c.f. QuadrangulationSubdivision.cpp)

    // bound Dijkstra by parent quad vertices
    std::vector<SimplexId> bounds{
      criticalPointsIdentifier_[q.i], criticalPointsIdentifier_[q.j],
      criticalPointsIdentifier_[q.k], criticalPointsIdentifier_[q.l]};

    for(size_t j = 0; j < outputDists.size(); ++j) {
      Dijkstra::shortestPath(
        midsNearestVertex[j], *triangulation_, outputDists[j], bounds);
    }

    auto inf = std::numeric_limits<float>::infinity();
    std::vector<float> sum(outputDists[0].size(), inf);

    for(size_t j = 0; j < sum.size(); ++j) {
      // skip if vertex j not in cell i
      if(morseSeg_[j] != static_cast<SimplexId>(i)) {
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
      = segmentationNumber_
        - std::count(
          sum.begin(), sum.end(), std::numeric_limits<float>::infinity());

    if(verticesInCell == 0) {
      std::stringstream msg;
      msg << MODULE_S "Barycenter in cell " << i << " not found." << std::endl;
      dMsg(std::cout, msg.str(), detailedInfoMsg);
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

    qsubd->emplace_back(Quad{4, q.i, sepMids[3], baryPos, sepMids[0]});
    qsubd->emplace_back(Quad{4, q.j, sepMids[0], baryPos, sepMids[1]});
    qsubd->emplace_back(Quad{4, q.k, sepMids[1], baryPos, sepMids[2]});
    qsubd->emplace_back(Quad{4, q.l, sepMids[2], baryPos, sepMids[3]});
  }

  // overwrite old quads
  outputCells_ = std::move(outputSubd);

  return 0;
}

// main routine
int ttk::SurfaceQuadrangulation::execute() {

  Timer t;

  // clear output
  outputCells_.clear();
  outputPoints_.clear();
  outputPointsIds_.clear();
  outputPointsCells_.clear();
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
    outputPointsCells_[i] = criticalPointsCellIds_[i];
  }

  // number of degenerate quadrangles
  size_t ndegen = 0;

  if(dualQuadrangulation_) {
    dualQuadrangulate();
  } else {
    // direct quadrangulation with saddle points
    quadrangulate(ndegen);
    subdivise();
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
