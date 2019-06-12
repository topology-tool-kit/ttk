#include <Dijkstra.h>
#include <Geometry.h>
#include <SurfaceQuadrangulation.h>
#include <cmath>
#include <queue>

std::set<ttk::SimplexId>
  ttk::SurfaceQuadrangulation::manifoldsAround(const SimplexId vert) const {

  // set of manifolds indices around vert
  std::set<SimplexId> manifolds{};
  // set of processed neighbors
  std::set<SimplexId> processedNeighbors{};
  // queue of neighbors to process
  std::queue<SimplexId> neighborsToProcess{};
  // maximum number of neighbors to process
  const size_t maxProcessedNeighbors = 20;

  neighborsToProcess.push(vert);

  // search in a neighborhood around vert to store the manifolds
  // indices of the encountered neighbors

  while(!neighborsToProcess.empty()) {
    auto curr = neighborsToProcess.front();
    neighborsToProcess.pop();
    manifolds.insert(segmentation_[curr]);
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

  return manifolds;
}

std::set<ttk::SimplexId> ttk::SurfaceQuadrangulation::commonManifolds(
  const std::vector<size_t> &verts) const {

  // get the manifold of every neighbors of our input vector points
  std::vector<std::set<SimplexId>> common_manifolds(verts.size());

  // TTK identifiers of input vertices
  std::vector<SimplexId> vertsId(verts.size());

  std::transform(verts.begin(), verts.end(), vertsId.begin(),
                 [&](size_t a) { return criticalPointsIdentifier_[a]; });

  for(size_t i = 0; i < verts.size(); ++i) {
    common_manifolds[i] = manifoldsAround(vertsId[i]);
  }

  // intersect every set to get a unique common manifold
  std::set<SimplexId> curr;
  auto last = common_manifolds[0];

  for(size_t i = 1; i < common_manifolds.size(); ++i) {
    std::set_intersection(last.begin(), last.end(), common_manifolds[i].begin(),
                          common_manifolds[i].end(),
                          std::inserter(curr, curr.begin()));
    std::swap(last, curr);
    curr.clear();
  }

  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] Common manifolds between vertices";
    for(auto &id : vertsId) {
      msg << " " << id;
    }
    msg << ":";
    for(auto &elem : last) {
      msg << " " << elem;
    }
    msg << std::endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);
  }

  return last;
}

int ttk::SurfaceQuadrangulation::dualQuadrangulate(
  const std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &sepEdges) {

  // quadrangles vertices are only extrema

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

int ttk::SurfaceQuadrangulation::quadrangulate(
  const std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &sepEdges,
  size_t &ndegen) {
  // quadrangle vertices are either extrema or saddle points

  // separatrices: destinations (extrema) -> sources (saddle points)
  std::vector<std::set<SimplexId>> sepMappingDests(sepEdges.size());

  // number of separatrices coming out of this edge
  std::vector<size_t> pointSepNumber(criticalPointsNumber_);

  for(auto &p : sepEdges) {
    for(SimplexId i = 0; i < criticalPointsNumber_; i++) {
      if(p.first == criticalPointsCellIds_[i]) {
        for(SimplexId j = 0; j < criticalPointsNumber_; j++) {
          if(p.second == criticalPointsCellIds_[j]) {
            sepMappingDests[j].insert(i);
            pointSepNumber[j]++;
            // should we also use the saddle points valence?
            pointSepNumber[i]++;
          }
        }
      }
    }
  }

  // iterate twice over dests
  for(size_t i = 0; i < sepMappingDests.size(); i++) {
    // skip if no sources
    if(sepMappingDests[i].empty()) {
      continue;
    }
    // begin second loop at i+1 to avoid duplicates and improve
    // performance
    for(size_t k = i + 1; k < sepMappingDests.size(); k++) {
      // skip same dest or if no sources
      if(k == i || sepMappingDests[k].empty()) {
        continue;
      }

      // skip if same type (we are looking for quadrangles containing
      // one minimum, one maximum and two saddle points)
      if(criticalPointsType_[k] == criticalPointsType_[i]) {
        continue;
      }

      // list of common sources to i and k
      std::vector<SimplexId> common_dests;
      std::set_intersection(
        sepMappingDests[i].begin(), sepMappingDests[i].end(),
        sepMappingDests[k].begin(), sepMappingDests[k].end(),
        std::back_inserter(common_dests));

      // find at least two common sources: j and l
      if(common_dests.size() >= 2) {

        // iterate over all possible common sources
        for(size_t m = 0; m < common_dests.size(); m++) {
          // avoid duplicates by beginning at m+1
          for(size_t n = m + 1; n < common_dests.size(); n++) {

            // gotcha!
            size_t j = common_dests[m];
            size_t l = common_dests[n];

            // check for a common shared manifold (looking around
            // saddle points only seems to be sufficient)
            std::vector<size_t> verts{j, l};

            bool validQuad = !commonManifolds(verts).empty();

            // fill output vector
            if(validQuad) {
              outputCells_.emplace_back(4);
              outputCells_.emplace_back(i);
              outputCells_.emplace_back(j);
              outputCells_.emplace_back(k);
              outputCells_.emplace_back(l);
            }
          }
        }
      } else if(common_dests.size() == 1
                && (sepMappingDests[i].size() == 1
                    || sepMappingDests[k].size() == 1)) {
        // we have degenerate quadrangles: i, j, k, j
        size_t j = common_dests[0];
        ndegen++;

        // fill output vector
        outputCells_.emplace_back(4);
        outputCells_.emplace_back(i);
        outputCells_.emplace_back(j);
        outputCells_.emplace_back(k);
        outputCells_.emplace_back(j);
      }
    }
  }

  return 0;
}

int ttk::SurfaceQuadrangulation::postProcess() {
  // post-processing: try to detect missing or extra quadrangles by
  // comparing separatrices number coming out of extrema

  // TODO remove extra quadrangles?

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
  // index of separatrices beginnings in separatrices arrays
  std::vector<size_t> sepBegs(numSeps);
  // index of separatrices endings in separatrices arrays
  std::vector<size_t> sepEnds(numSeps);
  // separatrices middles index in output points array
  std::vector<size_t> sepMiddle(numSeps);
  // separatrices middles nearest vertex id
  std::vector<SimplexId> sepMidNearestVertex(numSeps);
  // store duplicate separatrix id, -1 if no duplicate
  std::vector<SimplexId> sepDup(numSeps, -1);

  // fill in data arrays
  for(size_t i = 0; i < numSeps; ++i) {
    // separatrices bounds
    sepBegs[i] = sepFlatEdges[2 * i];
    sepEnds[i] = sepFlatEdges[2 * i + 1];
    // separatrices middles
    sepMiddle[i] = outputPoints_.size() / 3; // before insertion at next line
    sepMidNearestVertex[i] = findSeparatrixMiddle(sepBegs[i], sepEnds[i]);
  }

  for(size_t i = 0; i < numSeps; ++i) {
    for(size_t j = i + 1; j < numSeps; ++j) {
      if(sepCellIds_[sepBegs[i]] == sepCellIds_[sepBegs[j]]
         && sepCellIds_[sepEnds[i]] == sepCellIds_[sepEnds[j]]) {
        sepDup[i] = j;
        sepDup[j] = i;
      }
    }
  }

  // if output points are on a duplicate separatrix
  std::vector<bool> pointsDupSep(outputPointsIds_.size(), false);

  for(size_t i = 0; i < numSeps; ++i) {
    if(sepDup[i] != -1) {
      for(SimplexId j = 0; j < criticalPointsNumber_; ++j) {
        if(criticalPointsCellIds_[j] == sepCellIds_[sepBegs[i]]
           || criticalPointsCellIds_[j] == sepCellIds_[sepEnds[i]]) {
          pointsDupSep[j] = true;
        }
      }
      pointsDupSep[sepMiddle[i]] = true;
    }
  }

  // ad-hoc quad data structure (see QuadrangulationSubdivision.h)
  struct Quad {
    long long n;
    long long i;
    long long j;
    long long k;
    long long l;
  };

  // for each vertex on a separatrix, the index of the separatrix
  std::vector<SimplexId> onSep(segmentationNumber_, -1);

  {
    size_t currSep = 0;
    // iterate over separatrices lines to fill in onSep vector
    for(SimplexId i = 0; i < separatriceNumber_; ++i) {
      if(sepMask_[i] == 0) {
        currSep
          = std::find(sepBegs.begin(), sepBegs.end(), i) - sepBegs.begin();
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

  // rectified Morse-Smale cells
  std::vector<SimplexId> morseManRect(segmentationNumber_, -1);

  // number of Morse-Smale cells
  auto maxManId
    = *std::max_element(segmentation_, segmentation_ + segmentationNumber_);

  // store current number of quads
  auto nquads = outputCells_.size() / 5;

  // for each cell, the indices of the bordering separatrices
  std::vector<std::vector<SimplexId>> cellSeps{};

  bool notFinished = true;
  size_t pos = 0;

  while(notFinished) {
    for(size_t j = pos; j < segmentationNumber_; ++j) {
      if(onSep[j] == -1 && morseManRect[j] == -1) {
        pos = j;
        break;
      }
      if(j == segmentationNumber_ - 1) {
        notFinished = false;
      }
    }
    detectCells(pos, morseManRect, cellSeps, onSep);
  }

  // hold quad subdivision
  decltype(outputCells_) outQ{};
  outQ.reserve(5 * cellSeps.size());
  auto quads = reinterpret_cast<std::vector<Quad> *>(&outQ);

  for(const auto &c : cellSeps) {

    std::vector<long long> srcs(c.size());
    std::vector<long long> dsts(c.size());

    for(size_t i = 0; i < c.size(); ++i) {
      auto src = sepCellIds_[sepBegs[c[i]]];
      auto dst = sepCellIds_[sepEnds[c[i]]];
      for(long long j = 0; j < criticalPointsNumber_; ++j) {
        if(criticalPointsCellIds_[j] == src) {
          srcs[i] = j;
        }
        if(criticalPointsCellIds_[j] == dst) {
          dsts[i] = j;
        }
      }
    }

    bool found = false;

    // deal here with "small" cells bordered by less than 4 seps
    if(c.size() < 2) {
      continue; // don't bother with those cells
    } else if(c.size() < 4) {
      // pray for finding four different indices!
      std::set<long long> src_set(srcs.begin(), srcs.end());
      std::set<long long> dst_set(dsts.begin(), dsts.end());
      if(src_set.size() >= 2 && dst_set.size() >= 2) {
        // take two first distinct sources and dests
        auto vi = dsts[0];
        auto vk = *std::find_if_not(
          dsts.begin(), dsts.end(), [&](const long long a) { return a == vi; });
        auto vj = srcs[0];
        auto vl = *std::find_if_not(
          srcs.begin(), srcs.end(), [&](const long long a) { return a == vj; });
        quads->emplace_back(Quad{4, vi, vj, vk, vl});
        found = true;
      }
      continue;
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
        if(common_srcs_ik.size() > 1) {
          quads->emplace_back(
            Quad{4, vi, common_srcs_ik[0], vk, common_srcs_ik[1]});
          found = true;
          break;
        }
      }
      if(found) {
        // stop at the first correct quad
        break;
      }
    }

    if(!found) { // not found: degenerate case
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
        }
        if(found) {
          break;
        }
      }
    }
  }

  return 0;

  // rectify Morse-Smale cells indices
  for(size_t i = 0; i < nquads; ++i) {
    auto q = reinterpret_cast<Quad *>(&outputCells_[5 * i]);
    // bounds indices in separatrices array
    std::vector<SimplexId> qVerts{
      criticalPointsCellIds_[q->i], criticalPointsCellIds_[q->j],
      criticalPointsCellIds_[q->k], criticalPointsCellIds_[q->l]};

    auto nDupVerts = std::count_if(
      &q->i, &q->l + 1, [&](const long long a) { return pointsDupSep[a]; });

    if(nDupVerts >= 2) {
      auto saddles = std::vector<size_t>{
        static_cast<size_t>(q->j), static_cast<size_t>(q->l)};
      rectifyManifoldIndex(morseManRect, onSep,
                           *commonManifolds(saddles).begin(), maxManId,
                           cellSeps[i]);
    }
  }

  // fill in the rest of the segmentation
  for(size_t i = 0; i < segmentationNumber_; ++i) {
    if(morseManRect[i] == -1) {
      morseManRect[i] = segmentation_[i];
    }
  }

  // for each output quad, its barycenter position in outputPoints_
  std::vector<size_t> cellBary(outputCells_.size());

  // find separatrix index from vertices
  auto sepFromPoints = [&](const long long src, const long long dst) {
    for(size_t i = 0; i < numSeps; ++i) {
      if(sepCellIds_[sepBegs[i]] == criticalPointsCellIds_[src]
         && sepCellIds_[sepEnds[i]] == criticalPointsCellIds_[dst]) {
        return i;
      }
    }
    return numSeps;
  };

  std::vector<std::vector<float>> outputDists(4);

  // hold quad subdivision
  decltype(outputCells_) quadSubd{};
  quadSubd.reserve(4 * outputCells_.size());
  auto subd = reinterpret_cast<std::vector<Quad> *>(&quadSubd);

  for(size_t i = 0; i < nquads; ++i) {
    auto q = reinterpret_cast<Quad *>(&outputCells_[5 * i]);
    std::vector<size_t> quadSeps
      = {sepFromPoints(q->j, q->i), sepFromPoints(q->j, q->k),
         sepFromPoints(q->l, q->k), sepFromPoints(q->l, q->i)};

    // TODO if dup separatrice

    std::vector<long long> sepMids(quadSeps.size());
    std::vector<SimplexId> midsNearestVertex(quadSeps.size());
    for(size_t j = 0; j < quadSeps.size(); ++j) {
      sepMids[j] = sepMiddle[quadSeps[j]];
      midsNearestVertex[j] = sepMidNearestVertex[quadSeps[j]];
    }

    // find barycenter of current cell (c.f. QuadrangulationSubdivision.cpp)

    std::vector<SimplexId> bounds{
      criticalPointsIdentifier_[q->i], criticalPointsIdentifier_[q->j],
      criticalPointsIdentifier_[q->k], criticalPointsIdentifier_[q->l]};

    for(size_t j = 0; j < quadSeps.size(); ++j) {
      Dijkstra::shortestPath(
        midsNearestVertex[j], *triangulation_, outputDists[j], bounds);
    }

    auto inf = std::numeric_limits<float>::infinity();
    std::vector<float> sum(outputDists[0].size(), inf);

    for(size_t j = 0; j < sum.size(); ++j) {
      auto m = outputDists[0][j];
      auto n = outputDists[1][j];
      auto o = outputDists[2][j];
      auto p = outputDists[3][j];
      if(m == inf || n == inf || o == inf || p == inf) {
        continue;
      }
      // cost to minimize
      sum[j] = m + n + o + p;
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
    }

    subd->emplace_back(Quad{4, q->i, sepMids[3], baryPos, sepMids[0]});
    subd->emplace_back(Quad{4, q->j, sepMids[0], baryPos, sepMids[1]});
    subd->emplace_back(Quad{4, q->k, sepMids[1], baryPos, sepMids[2]});
    subd->emplace_back(Quad{4, q->l, sepMids[2], baryPos, sepMids[3]});
  }

  outputCells_ = std::move(quadSubd);

  return 0;
}

int ttk::SurfaceQuadrangulation::detectCells(
  const SimplexId src,
  std::vector<SimplexId> &vertexCells,
  std::vector<std::vector<SimplexId>> &cellSeps,
  const std::vector<SimplexId> &vertexSepMask) const {

  std::queue<SimplexId> toProcess{};
  toProcess.push(src);
  std::vector<SimplexId> borderSeps{};

  auto cellId = segmentation_[src];
  SimplexId newCellId = cellSeps.size();

  auto isCandidate = [&](const SimplexId a) {
    return segmentation_[a] == cellId && vertexSepMask[a] == -1
           && vertexCells[a] == -1;
  };

  while(!toProcess.empty()) {
    auto curr = toProcess.front();
    toProcess.pop();

    if(!isCandidate(curr)) {
      continue;
    }

    vertexCells[curr] = newCellId;

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

  return 0;
}

int ttk::SurfaceQuadrangulation::rectifyManifoldIndex(
  std::vector<SimplexId> &morseManRect,
  const std::vector<SimplexId> &onSep,
  const SimplexId sharedManifold,
  SimplexId &maxManifoldId,
  std::vector<SimplexId> &sepsIndex) const {

  std::queue<SimplexId> toProcess{};
  std::vector<SimplexId> borderSeps{};

  // look for the first vertex meeting constraints
  for(size_t i = 0; i < segmentationNumber_; ++i) {
    if(segmentation_[i] == sharedManifold && onSep[i] == -1
       && morseManRect[i] == -1) {
      toProcess.push(i);
      break;
    }
  }

  auto isCandidate = [&](const SimplexId a) {
    return segmentation_[a] == sharedManifold && onSep[a] == -1
           && morseManRect[a] == -1;
  };

  // breadth-first search
  while(!toProcess.empty()) {
    auto curr = toProcess.front();
    toProcess.pop();

    if(!isCandidate(curr)) {
      continue;
    }

    morseManRect[curr] = maxManifoldId;

    auto nneigh = triangulation_->getVertexNeighborNumber(curr);
    for(SimplexId j = 0; j < nneigh; ++j) {
      SimplexId next;
      triangulation_->getVertexNeighbor(curr, j, next);

      if(isCandidate(next)) {
        toProcess.push(next);
      }

      // store reached separatrices indices
      if(onSep[next] != -1) {
        borderSeps.emplace_back(onSep[next]);
      }
    }
  }

  // histogram of border separatrices indices
  std::map<SimplexId, int> hist{};

  // post-process borderSeps to find the most common separatrices indices
  for(const auto &v : borderSeps) {
    hist[v]++;
  }

  // reverse histogram to sort map by values
  std::map<int, SimplexId> revHist{};
  for(const auto &p : hist) {
    revHist[p.second] = p.first;
  }

  const size_t sepsPerCell = 4;
  for(auto it = revHist.rbegin(); it != revHist.rend(); ++it) {
    sepsIndex.emplace_back(it->second);
    if(sepsIndex.size() == sepsPerCell) {
      break;
    }
  }

  // update total number of cells
  maxManifoldId++;

  return 0;
}

size_t ttk::SurfaceQuadrangulation::findSeparatrixMiddle(const size_t a,
                                                         const size_t b) {

  const int dim = 3;

  std::vector<float> distFromA(b - a + 1);
  std::array<float, dim> prev{}, curr{};

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

// main routine
int ttk::SurfaceQuadrangulation::execute() {

  Timer t;

  // clear output
  outputCells_.clear();
  outputPoints_.clear();
  outputPointsIds_.clear();
  outputPoints_.resize(3 * criticalPointsNumber_);
  outputPointsIds_.resize(criticalPointsNumber_);
  outputPointsTypes_.resize(criticalPointsNumber_);

  // fill in critical points 3d coordinates and identifiers
  for(SimplexId i = 0; i < criticalPointsNumber_; ++i) {
    outputPoints_[3 * i] = criticalPoints_[3 * i];
    outputPoints_[3 * i + 1] = criticalPoints_[3 * i + 1];
    outputPoints_[3 * i + 2] = criticalPoints_[3 * i + 2];
    outputPointsIds_[i] = criticalPointsIdentifier_[i];
    outputPointsTypes_[i] = 0;
  }

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
    msg << "[SurfaceQuadrangulation] Error: odd number of separatrices edges"
        << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
    return -1;
  }

  // holds separatrices edges for every separatrix
  std::vector<std::pair<SimplexId, SimplexId>> sepEdges{};

  for(size_t i = 0; i < sepFlatEdges.size() / 2; ++i) {
    sepEdges.emplace_back(
      std::make_pair(sepFlatEdges[2 * i], sepFlatEdges[2 * i + 1]));
  }

  // number of degenerate quadrangles
  size_t ndegen = 0;

  if(dualQuadrangulation_) {
    dualQuadrangulate(sepEdges);
  } else {
    // direct quadrangulation with saddle points
    quadrangulate(sepEdges, ndegen);
    postProcess();
  }

  // maximum manifold id
  auto maxSegId
    = *std::max_element(segmentation_, segmentation_ + segmentationNumber_);

  // total number of manifolds
  int nseg = maxSegId + 1;

  // number of produced quads
  size_t quadNumber = outputCells_.size() / 5;

  // print number of quadrangles wrt number of MSC segmentation
  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] " << quadNumber << " quads (" << ndegen
        << " degenerated, " << nseg << " manifolds)" << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] Produced " << quadNumber
        << " quadrangles after " << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  return 0;
}
