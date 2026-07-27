#include <MorseSmaleQuadrangulation.h>

int ttk::MorseSmaleQuadrangulation::findSepsVertices(
  const std::vector<size_t> &seps,
  std::vector<LongSimplexId> &srcs,
  std::vector<LongSimplexId> &dsts) const {

  srcs.resize(seps.size());
  dsts.resize(seps.size());

  for(size_t i = 0; i < seps.size(); ++i) {
    auto src = sepCellIds_[sepBegs_[seps[i]]];
    auto dst = sepCellIds_[sepEnds_[seps[i]]];
    auto src_dim = sepCellDims_[sepBegs_[seps[i]]];
    auto dst_dim = sepCellDims_[sepEnds_[seps[i]]];
    for(LongSimplexId j = 0; j < criticalPointsNumber_; ++j) {
      if(criticalPointsCellIds_[j] == src
         && criticalPointsType_[j] == src_dim) {
        srcs[i] = j;
      }
      if(criticalPointsCellIds_[j] == dst
         && criticalPointsType_[j] == dst_dim) {
        dsts[i] = j;
      }
    }
  }

  return 0;
}

int ttk::MorseSmaleQuadrangulation::dualQuadrangulate() {

  // iterate over separatrices middles to build quadrangles around
  // them: the separatrix vertices and two barycenters

  std::vector<Quad> dualQuads{};

  auto quadNeighbors = [&](const LongSimplexId a) {
    std::set<LongSimplexId> neighs{};
    for(const auto &q : outputCells_) {
      if(a == q[0] || a == q[2]) {
        neighs.emplace(q[1]);
        neighs.emplace(q[3]);
      }
      if(a == q[1] || a == q[3]) {
        neighs.emplace(q[0]);
        neighs.emplace(q[2]);
      }
    }
    return neighs;
  };

  for(size_t i = 0; i < outputPointsIds_.size(); ++i) {
    // only keep sep middles
    if(outputPointsTypes_[i] != 1) {
      continue;
    }

    // get a list of direct neighbors
    auto neighs = quadNeighbors(i);

    // skip sep middle between single extremum and saddle in
    // degenerate cell (point not used)
    if(neighs.empty()) {
      continue;
    }

    // neighs should contain exactly 4 indices, two corresponding to
    // critical points, and two for generated points, mostly cell
    // barycenters
    std::vector<LongSimplexId> crit{}, gen{};
    for(const auto n : neighs) {
      if(n < criticalPointsNumber_) {
        crit.emplace_back(n);
      } else {
        gen.emplace_back(n);
      }
    }

    dualQuads.emplace_back(Quad{crit[0], gen[0], crit[1], gen[1]});
  }

  // for degenerate quadrangles, iterate over the single extremum, and
  // use the three generated points v0, v1 and v2
  for(size_t i = 0; i < outputPointsIds_.size(); ++i) {

    // only keep v0 points in degenerate cells
    if(outputPointsTypes_[i] != 3) {
      continue;
    }

    // v0 has 3 neighbors: v1, v2 and the double extremum

    auto v0neighs = quadNeighbors(i);
    std::vector<LongSimplexId> crit{}, gen{};
    for(const auto n : v0neighs) {
      if(n < criticalPointsNumber_) {
        crit.emplace_back(n);
      } else {
        gen.emplace_back(n);
      }
    }

    // v1 has 3 neighbors: v0, a sep middle and the single extremum

    // follow v1 to get the single extremum
    auto v1neighs = quadNeighbors(gen[0]);
    LongSimplexId single{};
    for(const auto n : v1neighs) {
      if(n < criticalPointsNumber_) {
        single = n;
      }
    }

    // the single extremum has 3 neighbors: v1, v2 and the saddle

    // follow single to find saddle
    auto sneighs = quadNeighbors(single);
    for(const auto n : sneighs) {
      if(n < criticalPointsNumber_) {
        crit.emplace_back(n);
      }
    }

    dualQuads.emplace_back(Quad{crit[0], gen[0], crit[1], gen[1]});
  }

  // clean output points, remove points with no neighbors
  std::vector<bool> hasNeighs(outputPointsIds_.size(), false);
  for(const auto &q : dualQuads) {
    for(const auto &qv : q) {
      hasNeighs[qv] = true;
    }
  }

  std::vector<SimplexId> newPos(outputPointsIds_.size(), -1);
  SimplexId pos{0};
  for(size_t i = 0; i < outputPointsIds_.size(); ++i) {
    if(hasNeighs[i]) {
      newPos[i] = pos++;
    }
  }

  decltype(this->outputPoints_) newPoints(3 * pos);
  decltype(this->outputPointsIds_) newPointsIds(pos);
  decltype(this->outputPointsTypes_) newPointsTypes(pos);
  decltype(this->outputPointsCells_) newPointsCells(pos);

  for(size_t i = 0; i < outputPointsIds_.size(); ++i) {
    if(hasNeighs[i]) {
      newPoints[3 * newPos[i] + 0] = this->outputPoints_[3 * i + 0];
      newPoints[3 * newPos[i] + 1] = this->outputPoints_[3 * i + 1];
      newPoints[3 * newPos[i] + 2] = this->outputPoints_[3 * i + 2];
      newPointsIds[newPos[i]] = this->outputPointsIds_[i];
      newPointsTypes[newPos[i]] = this->outputPointsTypes_[i];
      newPointsCells[newPos[i]] = this->outputPointsCells_[i];
    }
  }

  this->outputPoints_ = std::move(newPoints);
  this->outputPointsIds_ = std::move(newPointsIds);
  this->outputPointsTypes_ = std::move(newPointsTypes);
  this->outputPointsCells_ = std::move(newPointsCells);

  // compact cells
  for(auto &q : dualQuads) {
    for(auto &qv : q) {
      qv = newPos[qv];
    }
  }

  outputCells_ = std::move(dualQuads);

  return 0;
}

void ttk::MorseSmaleQuadrangulation::clearData() {
  outputCells_.clear();
  outputPoints_.clear();
  outputPointsIds_.clear();
  outputPointsCells_.clear();
}
