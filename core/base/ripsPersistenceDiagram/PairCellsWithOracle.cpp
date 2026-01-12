#include <PairCellsWithOracle.h>
#include <ripser.h>

#include <numeric>

ttk::rpd::PairCellsWithOracle::PairCellsWithOracle(
  const PointCloud &points,
  MultidimensionalDiagram const &oracle,
  bool distanceMatrix,
  bool parallelSort)
  : n_(distanceMatrix ? (1 + sqrt(1 + 8 * points[0].size())) / 2
                      : points.size()),
    parallelSort_(parallelSort), oracle_(oracle) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PairCellsWithOracle");

  if(distanceMatrix)
    compressedDM_ = points[0];
  else {
    const unsigned dim = points[0].size();
    for(int i = 1; i < n_; ++i) {
      for(int j = 0; j < i; ++j) {
        double s = 0.;
        for(unsigned d = 0; d < dim; ++d)
          s += (points[i][d] - points[j][d]) * (points[i][d] - points[j][d]);
        compressedDM_.push_back(sqrt(s));
      }
    }
  }
}

ttk::rpd::PairCellsWithOracle::PairCellsWithOracle(
  float *data,
  int n,
  int dim,
  MultidimensionalDiagram const &oracle,
  bool parallelSort)
  : n_(n), parallelSort_(parallelSort), oracle_(oracle) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PairCellsWithOracle");

  for(int i = 1; i < n_; ++i) {
    for(int j = 0; j < i; ++j) {
      double s = 0.;
      for(int d = 0; d < dim; ++d)
        s += (data[dim * i + d] - data[dim * j + d])
             * (data[dim * i + d] - data[dim * j + d]);
      compressedDM_.push_back(sqrt(s));
    }
  }
}

void ttk::rpd::PairCellsWithOracle::callOracle(const PointCloud &points,
                                               MultidimensionalDiagram &oracle,
                                               double threshold,
                                               bool distanceMatrix) {
  ripser::ripser(points, oracle, threshold, 1, distanceMatrix, false, false);

  std::sort(
    oracle[1].begin(), oracle[1].end(),
    [](const PersistencePair &x, const PersistencePair &y) {
      if(x.second.second == y.second.second) { // reverse co-lexicographic order
                                               // to match Ripser's choice
        Simplex const &t1 = x.second.first;
        Simplex const &t2 = y.second.first;
        return !std::lexicographical_compare(
          t1.rbegin(), t1.rend(), t2.rbegin(), t2.rend());
      } else
        return x.second.second < y.second.second;
    });
}

void ttk::rpd::PairCellsWithOracle::run() {
  Timer tm{};
  for(auto const &[b, d] : oracle_[1]) {
    if(d.second < inf)
      bound_ = std::max(bound_, d.second);
  }
  printMsg("Upper bound from the PD oracle: " + std::to_string(bound_));
  initializeWithBound();
  printMsg("Initialized (#p=" + std::to_string(n_)
             + ", #e=" + std::to_string(nEdges_) + ")",
           0, tm.getElapsedTime());
  pairCellsWithOracle();
  printMsg("Complete", 1, tm.getElapsedTime());
}

void ttk::rpd::PairCellsWithOracle::initializeWithBound() {
  // construct edges
  graph_.resize(n_);
  for(id_t i = 0; i < n_; ++i) {
    for(id_t j = i + 1; j < n_; ++j) {
      if(DM(i, j) <= bound_) {
        edgeToIndex_[{i, j}] = edges_.size();
        edges_.push_back({{i, j}, DM(i, j)});
        graph_[i].insert(j);
        graph_[j].insert(i);
      }
    }
  }
  nEdges_ = edges_.size();

  // sort edges
  edgesIndices_.resize(nEdges_);
  edgesOrder_.resize(nEdges_);
  edgesPartner_.resize(nEdges_, -1);
  std::iota(edgesIndices_.begin(), edgesIndices_.end(), 0);
  if(parallelSort_)
    TTK_PSORT(globalThreadNumber_, edgesIndices_.begin(), edgesIndices_.end(),
              [&](id_t i, id_t j) { return edges_[i] < edges_[j]; })
  else
    std::sort(edgesIndices_.begin(), edgesIndices_.end(),
              [&](id_t i, id_t j) { return edges_[i] < edges_[j]; });
  for(unsigned i = 0; i < edgesIndices_.size(); ++i)
    edgesOrder_[edgesIndices_[i]] = i;
}

void ttk::rpd::PairCellsWithOracle::pairCellsWithOracle() {
  for(const auto &[birth, death] : oracle_[1]) {
    if(death.second < inf) {
      const id_t e_id = edgeToIndex_[{birth.first[0], birth.first[1]}];
      const id_t t_id = triangles_.size();
      triangles_.push_back(
        {{death.first[0], death.first[1], death.first[2]}, death.second});
      boundaries_.push_back({edgeToIndex_[{death.first[0], death.first[1]}],
                             edgeToIndex_[{death.first[1], death.first[2]}],
                             edgeToIndex_[{death.first[0], death.first[2]}]});
      cascadeEdges_.emplace_back();
      trianglesPartner_.push_back(e_id);
      edgesPartner_[e_id] = t_id;
      eliminateBoundaryWithOracle(t_id, e_id);
    }
  }
}

void ttk::rpd::PairCellsWithOracle::eliminateBoundaryWithOracle(id_t t_id,
                                                                id_t e_id) {
  BoundaryContainer bc(boundaries_[t_id], nEdges_);
  id_t youngest_id = -1;
  while(youngest_id != e_id) {
    youngest_id = *std::max_element(
      boundaries_[t_id].begin(), boundaries_[t_id].end(),
      [&](id_t i, id_t j) { return edgesOrder_[i] < edgesOrder_[j]; });
    if(youngest_id == e_id)
      return;
    cascadeEdges_[t_id].push_back(youngest_id);
    if(edgesPartner_[youngest_id] != -1)
      bc.exclusiveAddBoundary(boundaries_[edgesPartner_[youngest_id]]);
    else {
      const auto &[i, j] = edges_[youngest_id].e;
      const double d = edges_[youngest_id].d;
      for(id_t const &k : graph_[i]) {
        if(DM(i, k) < d && DM(j, k) < d) {
          bc.exclusiveAddBoundary({youngest_id, edgeToIndex_[std::minmax(i, k)],
                                   edgeToIndex_[std::minmax(j, k)]});
          break;
        }
      }
    }
  }
}

void ttk::rpd::PairCellsWithOracle::getGenerators(
  std::vector<Generator1> &generators) const {
  for(unsigned i = 0; i < triangles_.size(); ++i) {
    EdgeSet boundary;
    for(id_t const &e_id : boundaries_[i])
      boundary.emplace_back(edges_[e_id].e);
    generators.emplace_back(
      boundary,
      std::make_pair(edges_[trianglesPartner_[i]].d, triangles_[i].d));
  }
}

void ttk::rpd::PairCellsWithOracle::getCascades(std::vector<Cascade> &cascades,
                                                EdgeSets3 &critical) const {
  fillRNG(critical);
  for(auto const &c : cascadeEdges_) {
    critical[DEATH1].emplace_back(edges_[c[0]].e);
    Cascade cascade;
    for(id_t const &e_id : c)
      cascade.emplace_back(edges_[e_id].e);
    cascades.push_back(cascade);
  }
}

void ttk::rpd::PairCellsWithOracle::getCascades(EdgeSets4 &critical) const {
  fillRNG(critical);
  std::set<id_t> cascadeSet;
  for(auto const &c : cascadeEdges_) {
    critical[DEATH1].emplace_back(edges_[c[0]].e);
    for(unsigned i = 1; i < c.size(); ++i)
      cascadeSet.insert(c[i]);
  }
  for(const id_t &e_id : cascadeSet)
    critical[CASC1].emplace_back(edges_[e_id].e);
}