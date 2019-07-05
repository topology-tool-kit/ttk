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

/**
 * @brief Sort map keys according to decreasing values
 */
inline std::vector<size_t> sortHistWeight(const std::map<size_t, int> &hist) {
  // dump map into vector
  std::vector<std::pair<size_t, int>> histVec{};
  histVec.reserve(hist.size());

  for(const auto &p : hist) {
    histVec.emplace_back(p);
  }

  // sort by value, descending order
  std::sort(
    histVec.begin(), histVec.end(),
    [&](const std::pair<size_t, int> &a, const std::pair<size_t, int> &b) {
      return a.second > b.second;
    });

  std::vector<size_t> res(histVec.size());

  // extract identifiers only
  for(size_t i = 0; i < histVec.size(); ++i) {
    res[i] = histVec[i].first;
  }

  return res;
}

int ttk::SurfaceQuadrangulation::detectCells(
  const SimplexId src, const std::vector<SimplexId> &vertexSepMask) {

  std::vector<bool> cellMask(morseSeg_.size(), false);

  std::queue<SimplexId> toProcess{};
  toProcess.push(src);
  std::vector<SimplexId> borderSeps{};

  auto cellId = segmentation_[src];
  SimplexId newCellId = cellSeps_.size();

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
  std::map<size_t, int> hist{};

  // post-process borderSeps to find the most common separatrices indices
  for(const auto &v : borderSeps) {
    hist[v]++;
  }

  // filter out small cells
  if(hist.size() < 3) {
    return 1;
  }

  // copy sub-segmentation into output
  for(size_t i = 0; i < morseSeg_.size(); ++i) {
    if(cellMask[i]) {
      morseSeg_[i] = newCellId;
    }
  }

  auto sepIds = sortHistWeight(hist);

  // return all reached separatrices by importance order
  cellSeps_.emplace_back(sepIds);
  // link this cell to MorseSmale Manifold index
  cellMMId_.emplace_back(segmentation_[src]);
  // add new cell index
  cellId_.emplace_back(newCellId);

  return 0;
}

int ttk::SurfaceQuadrangulation::mergeSmallCells(
  const std::vector<SimplexId> &vertexSepMask) {
  // fuse two (small) cells sharing their cellId and at least three separatrices

  // detect cells sharing the same cell id
  std::vector<std::vector<size_t>> sharedCells{};
  for(size_t i = 0; i < cellMMId_.size(); ++i) {
    // check if i not already included
    bool skip = false;
    for(const auto &s : sharedCells) {
      if(std::find(s.begin(), s.end(), i) != s.end()) {
        skip = true;
      }
    }
    if(skip) {
      continue;
    }

    std::set<size_t> shared_i{i};
    for(size_t j = i + 1; j < cellMMId_.size(); ++j) {
      if(cellMMId_[j] == cellMMId_[i]) {
        shared_i.emplace(j);
      }
    }
    if(shared_i.size() > 1) {
      sharedCells.emplace_back(shared_i.begin(), shared_i.end());
    }
  }

  // store duplicate merged cell indices
  std::vector<size_t> cellsToRemove{};

  for(const auto &c : sharedCells) {
    std::vector<size_t> cellSize(c.size());

    // mesure the size of every cell
    for(size_t i = 0; i < c.size(); ++i) {
      for(size_t j = 0; j < segmentationNumber_; ++j) {
        if(morseSeg_[j] == static_cast<SimplexId>(c[i])) {
          cellSize[i]++;
        }
      }
    }

    const size_t minCellSize{100};
    if(std::accumulate(cellSize.begin(), cellSize.end(), 0UL) > minCellSize) {
      continue;
    }

    // small cells: we should merge everything into c[0]
    for(size_t i = 1; i < c.size(); ++i) {
      for(size_t j = 0; j < segmentationNumber_; ++j) {
        if(morseSeg_[j] == static_cast<SimplexId>(c[i])) {
          morseSeg_[j] = c[0];
        }
      }
    }
    // search for separatrices in the whole area
    std::map<size_t, int> hist{};
    for(size_t i = 0; i < segmentationNumber_; ++i) {
      if(morseSeg_[i] != static_cast<SimplexId>(c[0])) {
        continue;
      }
      auto nneigh = triangulation_->getVertexNeighborNumber(i);
      for(SimplexId j = 0; j < nneigh; ++j) {
        SimplexId next;
        triangulation_->getVertexNeighbor(i, j, next);
        if(vertexSepMask[next] != -1) {
          hist[vertexSepMask[next]]++;
        }
      }
    }
    auto sortedSeps = sortHistWeight(hist);
    cellSeps_[c[0]] = sortedSeps;
    // cleanup merged cells
    for(size_t i = 1; i < c.size(); ++i) {
      cellsToRemove.emplace_back(c[i]);
    }
  }

  // actually remove duplicate merged cells
  if(!cellsToRemove.empty()) {

    // sort cell indices by decreasing order
    std::sort(cellsToRemove.begin(), cellsToRemove.end());
    std::reverse(cellsToRemove.begin(), cellsToRemove.end());

    // remove the end first
    for(const auto i : cellsToRemove) {
      cellId_.erase(std::next(cellId_.begin(), i));
      cellMMId_.erase(std::next(cellMMId_.begin(), i));
      cellSeps_.erase(std::next(cellSeps_.begin(), i));
    }
  }

  return 0;
}

std::set<ttk::SimplexId>
  ttk::SurfaceQuadrangulation::cellsAround(const SimplexId vert) const {

  // set of manifolds indices around vert
  std::set<SimplexId> cellsId{};
  // set of processed neighbors
  std::set<SimplexId> processedNeighbors{};
  // queue of neighbors to process
  std::queue<SimplexId> neighborsToProcess{};
  // maximum number of neighbors to process
  const size_t maxProcessedNeighbors = 50;

  neighborsToProcess.push(vert);

  // search in a neighborhood around vert to store the manifolds
  // indices of the encountered neighbors

  while(!neighborsToProcess.empty()) {
    auto curr = neighborsToProcess.front();
    neighborsToProcess.pop();
    cellsId.insert(morseSeg_[curr]);
    processedNeighbors.insert(curr);
    auto nneigh = triangulation_->getVertexNeighborNumber(curr);
    for(SimplexId j = 0; j < nneigh; ++j) {
      SimplexId next;
      triangulation_->getVertexNeighbor(curr, j, next);
      if(processedNeighbors.find(next) == processedNeighbors.end()) {
        neighborsToProcess.push(next);
      }
    }
    if(processedNeighbors.size() > maxProcessedNeighbors) {
      break;
    }
  }

  return cellsId;
}

int ttk::SurfaceQuadrangulation::sweepOverCells() {
  BarycentricSubdivision bs{};
  Triangulation newT{};

  bs.setupTriangulation(triangulation_);
  bs.setOutputTriangulation(&newT);
  bs.setInputPoints(inputPoints_);
  bs.execute();

  newT.preprocessVertexNeighbors();

  auto nVerts = triangulation_->getNumberOfVertices();
  auto nEdges = triangulation_->getNumberOfEdges();
  auto nTriangles = triangulation_->getNumberOfTriangles();

  // store separatrix index on subdivised triangulation vertices
  std::vector<SimplexId> onSep(newT.getNumberOfVertices(), -1);
  // store separatrix index on subdivised triangulation edges
  std::vector<SimplexId> edgeOnSep(newT.getNumberOfEdges(), -1);
  // count the number of critical points encountered
  size_t critPoints{0};
  if(sepMask_[0] == 0) {
    critPoints++;
  }

  auto critPointId = [&](const SimplexId a) {
    if(sepCellDims_[a] == 0) {
      return sepCellIds_[a];
    } else if(sepCellDims_[a] == 1) {
      return nVerts + sepCellIds_[a];
    } else if(sepCellDims_[a] == 2) {
      return nVerts + nEdges + sepCellIds_[a];
    }
    return -1;
  };

  // init first vertex id
  SimplexId prev{critPointId(0)};

  // begin loop at one to get the first edge
  for(SimplexId i = 1; i < separatriceNumber_; ++i) {

    // for computing the separatrix id
    if(sepMask_[i] == 0) {
      critPoints++;
    }

    // current separatrix id is critPoints // 2
    auto currSepId = critPoints >> 1;

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

  // look around the saddle points
  for(SimplexId i = 0; i < criticalPointsNumber_; ++i) {
    // keep only saddle points
    if(criticalPointsType_[i] != 1) {
      continue;
    }
    SimplexId saddle = nVerts + criticalPointsCellIds_[i];
    auto nn = newT.getVertexNeighborNumber(saddle);
    if(nn != 4) {
      std::cout << "different number of neighbors: " << nn << std::endl;
    }

    std::set<SimplexId> sepsAround{};

    for(SimplexId j = 0; j < nn; ++j) {
      SimplexId neigh;
      newT.getVertexNeighbor(saddle, j, neigh);
      sepsAround.emplace(onSep[neigh]);
    }

    if(sepsAround.size() != 4) {
      std::cout << "different seps ids around " << criticalPointsCellIds_[i]
                << std::endl;
    }

    // propagate from triangles around saddle
    std::vector<bool> processed(nTriangles, false);

    auto nt = newT.getVertexTriangleNumber(saddle);
    for(SimplexId j = 0; j < nt; ++j) {
      std::queue<SimplexId> toProcess{};
      SimplexId tr;
      newT.getVertexTriangle(saddle, j, tr);

      std::set<size_t> sepIds{};
      toProcess.push(tr);
      {
        SimplexId a, b, c;
        newT.getTriangleVertex(tr, 0, a);
        newT.getTriangleVertex(tr, 1, b);
        newT.getTriangleVertex(tr, 2, c);

        if(a == i) {
          sepIds.emplace(onSep[b]);
          sepIds.emplace(onSep[c]);
        } else if(b == i) {
          sepIds.emplace(onSep[a]);
          sepIds.emplace(onSep[c]);
        } else if(c == i) {
          sepIds.emplace(onSep[a]);
          sepIds.emplace(onSep[b]);
        }
      }

      while(!toProcess.empty()) {
        auto curr = toProcess.front();
        toProcess.pop();

        // stop condition
        SimplexId a, b, c;
        newT.getTriangleVertex(tr, 0, a);
        newT.getTriangleVertex(tr, 1, b);
        newT.getTriangleVertex(tr, 2, c);
        if(saddlesId.find(a) != saddlesId.end()) {
        } else if(saddlesId.find(b) != saddlesId.end()) {
        } else if(saddlesId.find(c) != saddlesId.end()) {
        }

        auto ne = newT.getTriangleEdgeNumber(curr);
        for(SimplexId k = 0; k < ne; ++k) {
          SimplexId e;
          newT.getTriangleEdge(curr, k, e);
          SimplexId e0, e1;
          newT.getEdgeVertex(e, 0, e0);
          newT.getEdgeVertex(e, 1, e1);
          // skip edges on separatrices
          if(onSep[e0] != -1 && onSep[e1] != -1) {
            continue;
          }
          // push the triangle
        }

        processed[curr] = true;
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
  morseSeg_.resize(segmentationNumber_);
  std::fill(morseSeg_.begin(), morseSeg_.end(), -1);
  cellId_.clear();
  cellMMId_.clear();
  cellSeps_.clear();
  quadSeps_.clear();

  // fill in data arrays
  for(size_t i = 0; i < numSeps; ++i) {
    // separatrices bounds
    sepBegs_[i] = sepFlatEdges[2 * i];
    sepEnds_[i] = sepFlatEdges[2 * i + 1];
  }

  sweepOverCells();

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
    if(finished || pos >= segmentationNumber_) {
      break;
    }
    detectCells(pos, onSep);
    pos++;
  }

  // missing cells?
  // find them and store their bordering separatrices
  auto minCellId
    = *std::min_element(segmentation_, segmentation_ + segmentationNumber_);
  auto maxCellId
    = *std::max_element(segmentation_, segmentation_ + segmentationNumber_);
  for(SimplexId i = minCellId; i <= maxCellId; ++i) {
    if(std::find(cellMMId_.begin(), cellMMId_.end(), i) == cellMMId_.end()) {
      std::set<SimplexId> sepIds{};
      for(size_t j = 0; j < segmentationNumber_; ++j) {
        if(segmentation_[j] == i && onSep[j] != -1) {
          sepIds.emplace(onSep[j]);
        }
      }
      if(sepIds.size() > 2) {
        cellSeps_.emplace_back(sepIds.begin(), sepIds.end());
        cellMMId_.emplace_back(i);
        cellId_.emplace_back(cellId_.back() + 1);
      }
    }
  }

  mergeSmallCells(onSep);

  // hold quad subdivision
  outputCells_.reserve(5 * cellSeps_.size());
  auto quads = reinterpret_cast<std::vector<Quad> *>(&outputCells_);

  // for each cell the detected separatrices sources and destinations
  std::vector<std::vector<long long>> cellSrcs{}, cellDsts{};

  for(size_t a = 0; a < cellSeps_.size(); ++a) {

    auto &cs = cellSeps_[a];

    std::vector<long long> srcs{};
    std::vector<long long> dsts{};

    findSepsVertices(cs, srcs, dsts);

    // filter out separatrices whose sources are not in contact with
    // the current cell
    if(cs.size() > 5) {
      std::vector<size_t> sepsToRemove{};
      for(size_t i = 0; i < srcs.size(); ++i) {
        auto cellsAroundSrc = cellsAround(criticalPointsIdentifier_[srcs[i]]);
        if(std::find(cellsAroundSrc.begin(), cellsAroundSrc.end(), cellId_[a])
           == cellsAroundSrc.end()) {
          sepsToRemove.emplace_back(i);
        }
      }
      // remove the end first
      std::sort(sepsToRemove.rbegin(), sepsToRemove.rend());
      for(const auto i : sepsToRemove) {
        srcs.erase(std::next(srcs.begin(), i));
        dsts.erase(std::next(dsts.begin(), i));
        cs.erase(std::next(cs.begin(), i));
      }
    }

    cellSrcs.emplace_back(srcs);
    cellDsts.emplace_back(dsts);

    bool found = false;

    // optimisation: try first to only consider the first four separatrices

    if(cs.size() >= 4) {
      auto sbeg = srcs.begin();
      auto send = std::next(srcs.begin(), 4);
      auto dbeg = dsts.begin();
      auto dend = std::next(dsts.begin(), 4);
      std::set<SimplexId> srcs_set(sbeg, send);
      std::set<SimplexId> dsts_set(dbeg, dend);
      if(srcs_set.size() == 2 && dsts_set.size() == 2) {
        auto vi = *dsts_set.begin();
        auto vj = *srcs_set.begin();
        auto vk = *std::next(dsts_set.begin());
        auto vl = *std::next(srcs_set.begin());
        // check that every vertex appears twice
        if(std::count(sbeg, send, vj) == 2 && std::count(sbeg, send, vk) == 2
           && std::count(dbeg, dend, vi) == 2
           && std::count(dbeg, dend, vl) == 2) {
          quads->emplace_back(Quad{4, vi, vj, vk, vl});
          quadSeps_.emplace_back(cs.begin(), std::next(cs.begin(), 4));
          continue;
        }
      }
    }

    // normal case: find a pair of extrema and a pair of saddle points
    // with four separatrices linking one another

    // iterate over first dest
    for(size_t i = 0; i < cs.size(); ++i) {
      auto vi = dsts[i];

      // iterate over second dest from i + 1
      for(size_t k = i + 1; k < cs.size(); ++k) {
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
        for(size_t j = 0; j < cs.size(); ++j) {
          if(dsts[j] == vi) {
            srcs_i.insert(srcs[j]);
          }
          if(dsts[j] == vk) {
            srcs_k.insert(srcs[j]);
          }
        }
        std::set_intersection(srcs_i.begin(), srcs_i.end(), srcs_k.begin(),
                              srcs_k.end(), std::back_inserter(common_srcs_ik));
        if(common_srcs_ik.size() > 2) {
          // reorder common_srcs_ik according to srcs vector
          std::vector<SimplexId> common_srcs_ik_ordered{};
          auto cikb = common_srcs_ik.begin();
          auto cike = common_srcs_ik.end();
          for(const auto &s : srcs) {
            auto cikob = common_srcs_ik_ordered.begin();
            auto cikoe = common_srcs_ik_ordered.end();
            if(std::find(cikb, cike, s) != cike
               && std::find(cikob, cikoe, s) == cikoe) {
              common_srcs_ik_ordered.emplace_back(s);
            }
          }
          common_srcs_ik = std::move(common_srcs_ik_ordered);
        }
        if(common_srcs_ik.size() > 1) {
          auto vj = common_srcs_ik[0];
          auto vl = common_srcs_ik[1];
          quads->emplace_back(Quad{4, vi, vj, vk, vl});
          // find separatrices indices
          std::vector<size_t> seps{};

          auto findSep = [&](const long long src, const long long dst) {
            for(size_t j = 0; j < cs.size(); ++j) {
              if(srcs[j] == src && dsts[j] == dst) {
                seps.emplace_back(cs[j]);
                return;
              }
            }
            // if not found in cell seps, try first one from src to dst
            seps.emplace_back(sepFromPoints(src, dst));
          };

          findSep(vj, vi);
          findSep(vj, vk);
          findSep(vl, vk);
          findSep(vl, vi);
          quadSeps_.emplace_back(seps);
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
      if(cs.size() >= 3) {
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
        auto jk = sepFromPoints(vj, vk);
        auto kl = sepFromPoints(vl, vk);
        auto li = sepFromPoints(vl, vi);
        if(jk < numSeps && kl < numSeps && li < numSeps) {
          quads->emplace_back(Quad{4, vi, vj, vk, vl});
          quadSeps_.emplace_back(std::vector<size_t>{cs[0], jk, kl, li});
          continue;
        }
      }

      // degenerate case:
      // take the first two distinct extrema (separatrices with higher weight)
      for(size_t i = 0; i < cs.size(); ++i) {
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
          // store corresponding seps indices
          std::vector<size_t> seps{};
          for(size_t j = i; j < srcs.size(); ++j) {
            if(srcs[j] == srcs[i] && (dsts[j] == vi || dsts[j] == vk)) {
              seps.emplace_back(cs[j]);
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
          quadSeps_.emplace_back(seps);
          ndegen++;
          break;
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
        midsNearestVertex[j], *triangulation_, outputDists[j], bounds);
    }

    auto inf = std::numeric_limits<float>::infinity();
    std::vector<float> sum(outputDists[0].size(), inf);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
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

    qsubd->emplace_back(Quad{4, q.i, sepMids[3], baryPos, sepMids[0]});
    qsubd->emplace_back(Quad{4, q.j, sepMids[0], baryPos, sepMids[1]});
    qsubd->emplace_back(Quad{4, q.k, sepMids[1], baryPos, sepMids[2]});
    qsubd->emplace_back(Quad{4, q.l, sepMids[2], baryPos, sepMids[3]});
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
