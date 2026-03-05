#include <PairCells.h>

#include <numeric>
#include <boost/functional/hash.hpp>

using namespace ttk::rpd;

#ifdef TTK_ENABLE_CGAL
PairCells::PairCells(const std::vector<CGAL::Epick::Point_2> &points,
                     double upperBound,
                     bool parallelSort,
                     bool parallelMatrixConstruction)
  : n_(points.size()), bound_(upperBound), parallelSort_(parallelSort),
    parallelMatrixConstruction_(parallelMatrixConstruction) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PairCells");

  for(int i = 1; i < n_; ++i) {
    for(int j = 0; j < i; ++j)
      compressedDM_.push_back(
        sqrt(CGAL::squared_distance(points[i], points[j])));
  }
}
#endif

PairCells::PairCells(const PointCloud &points,
                     bool distanceMatrix,
                     double upperBound,
                     bool parallelSort,
                     bool parallelMatrixConstruction)
  : n_(distanceMatrix ? (1 + sqrt(1 + 8 * points[0].size())) / 2
                      : points.size()),
    bound_(upperBound), parallelSort_(parallelSort),
    parallelMatrixConstruction_(parallelMatrixConstruction) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PairCells");

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

PairCells::PairCells(float *data,
                     int n,
                     int dim,
                     double upperBound,
                     bool parallelSort,
                     bool parallelMatrixConstruction)
  : n_(n), bound_(upperBound), parallelSort_(parallelSort),
    parallelMatrixConstruction_(parallelMatrixConstruction) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PairCells");

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

void PairCells::run() {
  Timer tm{};
  if(bound_ == inf || bound_ < 0)
    initialize();
  else
    initializeWithBound();
  printMsg("Initialized (#p=" + std::to_string(n_)
             + ", #e=" + std::to_string(nEdges_)
             + ", #t=" + std::to_string(nTriangles_) + ")",
           0, tm.getElapsedTime());
  pairCells();
  printMsg("Complete", 1, tm.getElapsedTime());
}

void PairCells::initialize() {
  // edges
  for(id_t i = 0; i < n_; ++i) {
    for(id_t j = i + 1; j < n_; ++j)
      edges_.push_back({{i, j}, DM(i, j)});
  }
  nEdges_ = n_ * (n_ - 1) / 2;
  executeKruskal();

  // triangles
  nTriangles_ = n_ * (n_ - 1) * (n_ - 2) / 6;
  triangles_.resize(nTriangles_);
  boundaries_.resize(nTriangles_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for if(parallelMatrixConstruction_)
#else
  TTK_FORCE_USE(parallelMatrixConstruction_);
#endif // TTK_ENABLE_OPENMP
  for(id_t i = 0; i < n_ - 2; ++i) {
    const unsigned index_i
      = i * (i * i - 3 * i * (n_ - 1) + 3 * n_ * n_ - 6 * n_ + 2) / 6;
    for(id_t j = i + 1; j < n_ - 1; ++j) {
      const unsigned index_j = (i - j + 1) * (i + j - 2 * n_ + 2) / 2 + index_i;
      for(id_t k = j + 1; k < n_; ++k) {
        const double diameter
          = std::max(DM(i, j), std::max(DM(i, k), DM(j, k)));
        triangles_[index_j + (k - j - 1)] = {{i, j, k}, diameter};
        boundaries_[index_j + (k - j - 1)]
          = {nEdges_ - 1 - (n_ - i) * (n_ - i - 1) / 2 + (j - i),
             nEdges_ - 1 - (n_ - i) * (n_ - i - 1) / 2 + (k - i),
             nEdges_ - 1 - (n_ - j) * (n_ - j - 1) / 2 + (k - j)};
      }
    }
  }

  apparentPairs();
}

void PairCells::initializeWithBound() {
  // edges
  std::unordered_map<Edge, id_t, boost::hash<Edge>> edgeToIndex;
  std::vector<std::set<id_t>> graph(n_); // std::unordered_set?
  for(id_t i = 0; i < n_; ++i) {
    for(id_t j = i + 1; j < n_; ++j) {
      if(DM(i, j) <= bound_) {
        edgeToIndex[{i, j}] = edges_.size();
        edges_.push_back({{i, j}, DM(i, j)});
        graph[i].insert(j); // only directed graph with edges (ij), i<j
      }
    }
  }
  nEdges_ = edges_.size();
  executeKruskal();

  // triangles
  for(const auto &[e, f] : edges_) {
    const auto &[i, j] = e;
    for(auto const &k : graph[j]) {
      if(graph[i].find(k)
         != graph[i].end()) { // c++20 -> std::unordered_set::contains
        triangles_.push_back(
          {{i, j, k}, std::max(DM(i, j), std::max(DM(i, k), DM(j, k)))});
        boundaries_.push_back(
          {edgeToIndex[{i, j}], edgeToIndex[{i, k}], edgeToIndex[{j, k}]});
      }
    }
  }
  nTriangles_ = triangles_.size();

  apparentPairs();
}

void PairCells::executeKruskal() {
  UnionFind UF(n_);
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
  symbolicPerturbation();
  for(unsigned i = 0; i < edgesIndices_.size(); ++i)
    edgesOrder_[edgesIndices_[i]] = i;
  for(const id_t &id : edgesIndices_) {
    const auto &[s, f] = edges_[id];
    if(UF.find(s.first) != UF.find(s.second)) {
      UF.merge(s.first, s.second);
      edgesPartner_[id] = -2;
      if(++nPairedEdges_ == n_ - 1)
        break;
    }
  }
}

void PairCells::symbolicPerturbation(
  double eps) { // could be improved (many calls when lots of collisions)
  bool generic;
  do {
    generic = true;
    for(unsigned i = 0; i + 1 < edgesIndices_.size(); ++i) {
      if(edges_[edgesIndices_[i]].d >= edges_[edgesIndices_[i + 1]].d) {
        generic = false;
        FiltratedEdge &e = edges_[edgesIndices_[i + 1]];
        e.d += eps;
        DM(e.e.first, e.e.second) += eps;
      }
    }
    if(!generic) {
      if(parallelSort_)
        TTK_PSORT(globalThreadNumber_, edgesIndices_.begin(),
                  edgesIndices_.end(),
                  [&](id_t i, id_t j) { return edges_[i] < edges_[j]; })
      else
        std::sort(edgesIndices_.begin(), edgesIndices_.end(),
                  [&](id_t i, id_t j) { return edges_[i] < edges_[j]; });
    }
  } while(!generic);
}

void PairCells::apparentPairs() {
  // arrays initializations
  trianglesIndices_.resize(nTriangles_);
  trianglesPartner_.resize(nTriangles_, -1);
  cascadeEdges_.resize(nTriangles_);

  // full order
  std::iota(trianglesIndices_.begin(), trianglesIndices_.end(), 0);
  if(parallelSort_)
    TTK_PSORT(globalThreadNumber_, trianglesIndices_.begin(),
              trianglesIndices_.end(),
              [&](id_t i, id_t j) { return triangles_[i] < triangles_[j]; })
  else
    std::sort(trianglesIndices_.begin(), trianglesIndices_.end(),
              [&](id_t i, id_t j) { return triangles_[i] < triangles_[j]; });

  // apparent pairs
  id_t t = 0;
  for(const id_t &e : edgesIndices_) {
    const double f = edges_[e].d;
    while(triangles_[trianglesIndices_[t]].d < f)
      ++t;
    if(triangles_[trianglesIndices_[t]].d == f) {
      edgesPartner_[e] = trianglesIndices_[t];
      trianglesPartner_[trianglesIndices_[t]] = e;
      nPairedEdges_++;
    }
  }
}

void PairCells::pairCells() {
  for(const id_t &s : trianglesIndices_) {
    if(trianglesPartner_[s] == -1) {
      const int e = eliminateBoundaries(s);
      if(e != -1) {
        trianglesPartner_[s] = e;
        edgesPartner_[e] = s;
        if(++nPairedEdges_ == nEdges_)
          break;
      }
    }
  }
}

ttk::rpd::id_t PairCells::eliminateBoundaries(id_t s) {
  BoundaryContainer bc(boundaries_[s], nEdges_);
  while(!boundaries_[s].empty()) {
    const id_t e = *std::max_element(
      boundaries_[s].begin(), boundaries_[s].end(),
      [&](id_t i, id_t j) { return edgesOrder_[i] < edgesOrder_[j]; });
    cascadeEdges_[s].push_back(e);
    if(edgesPartner_[e] == -1)
      return e;
    else
      bc.exclusiveAddBoundary(boundaries_[edgesPartner_[e]]);
  }
  return -1;
}

void PairCells::getDiagram(MultidimensionalDiagram &diagrams) const {
  diagrams.resize(2);
  for(const id_t &e : edgesIndices_) {
    if(edgesPartner_[e] == -2)
      diagrams[0].emplace_back(
        FiltratedSimplex{{0}, 0.},
        FiltratedSimplex{{edges_[e].e.first, edges_[e].e.second}, edges_[e].d});
    else if(edgesPartner_[e] > 0
            && edges_[e].d < triangles_[edgesPartner_[e]].d) {
      const FiltratedSimplex birth(
        {edges_[e].e.first, edges_[e].e.second}, edges_[e].d);
      const FiltratedSimplex death(
        {std::get<0>(triangles_[edgesPartner_[e]].t),
         std::get<1>(triangles_[edgesPartner_[e]].t),
         std::get<2>(triangles_[edgesPartner_[e]].t)},
        triangles_[edgesPartner_[e]].d);
      diagrams[1].emplace_back(birth, death);
    }
  }
  diagrams[0].emplace_back(
    FiltratedSimplex{{0}, 0.}, FiltratedSimplex{{-1}, inf});
}

void PairCells::getDiagramAndGenerators(
  MultidimensionalDiagram &diagrams,
  std::vector<Generator1> &generators) const {
  diagrams.resize(2);
  for(const id_t &e : edgesIndices_) {
    if(edgesPartner_[e] == -2)
      diagrams[0].emplace_back(
        FiltratedSimplex{{0}, 0.},
        FiltratedSimplex{{edges_[e].e.first, edges_[e].e.second}, edges_[e].d});
    else if(edgesPartner_[e] > 0
            && edges_[e].d < triangles_[edgesPartner_[e]].d) {
      const FiltratedSimplex birth(
        {edges_[e].e.first, edges_[e].e.second}, edges_[e].d);
      const FiltratedSimplex death(
        {std::get<0>(triangles_[edgesPartner_[e]].t),
         std::get<1>(triangles_[edgesPartner_[e]].t),
         std::get<2>(triangles_[edgesPartner_[e]].t)},
        triangles_[edgesPartner_[e]].d);
      diagrams[1].emplace_back(birth, death);

      EdgeSet boundary;
      for(const id_t &edge : boundaries_[edgesPartner_[e]])
        boundary.emplace_back(edges_[edge].e);
      generators.emplace_back(
        boundary, std::make_pair(edges_[e].d, triangles_[edgesPartner_[e]].d));
    }
  }
  diagrams[0].emplace_back(
    FiltratedSimplex{{0}, 0.}, FiltratedSimplex{{-1}, inf});
}

void PairCells::getCascades(std::vector<Cascade> &cascades,
                            EdgeSets3 &critical) const {
  for(const id_t &e : edgesIndices_) {
    if(edgesPartner_[e] == -2) // MST
      critical[DEATH0].emplace_back(edges_[e].e);
    else if(edgesPartner_[e] > 0
            && edges_[e].d < triangles_[edgesPartner_[e]].d) {
      const Edge &killerEdge = edges_[cascadeEdges_[edgesPartner_[e]][0]].e;
      critical[BIRTH1].emplace_back(edges_[e].e); // RNG edge
      critical[DEATH1].emplace_back(killerEdge); // MML edge
      Cascade cascade = {killerEdge};
      for(unsigned i = 1; i < cascadeEdges_[edgesPartner_[e]].size() - 1; ++i)
        cascade.emplace_back(edges_[cascadeEdges_[edgesPartner_[e]][i]].e);
      cascades.emplace_back(cascade);
    }
  }
}

void PairCells::getCascades(EdgeSets4 &critical) const {
  std::set<id_t> cascadeSet;
  for(const id_t &e : edgesIndices_) {
    if(edgesPartner_[e] == -2) // MST
      critical[DEATH0].emplace_back(edges_[e].e);
    else if(edgesPartner_[e] > 0
            && edges_[e].d < triangles_[edgesPartner_[e]].d) {
      const Edge &killerEdge = edges_[cascadeEdges_[edgesPartner_[e]][0]].e;
      critical[BIRTH1].emplace_back(edges_[e].e); // RNG edge
      critical[DEATH1].emplace_back(killerEdge); // MML edge
      for(unsigned i = 1; i < cascadeEdges_[edgesPartner_[e]].size() - 1; ++i)
        cascadeSet.insert(cascadeEdges_[edgesPartner_[e]][i]);
    }
  }
  for(const id_t &e : cascadeSet)
    critical[CASC1].emplace_back(edges_[e].e);
}

void PairCells::enrichCascades(std::set<Edge> &cascadeSet,
                               EdgeSets4 &critical,
                               std::vector<int> const &globalIndices) const {
  for(const id_t &e : edgesIndices_) {
    if(edgesPartner_[e] > 0 && edges_[e].d < triangles_[edgesPartner_[e]].d) {
      const Edge &killerEdge = edges_[cascadeEdges_[edgesPartner_[e]][0]].e;
      critical[DEATH1].emplace_back(
        globalIndices[killerEdge.first],
        globalIndices[killerEdge.second]); // MML edge
      for(unsigned i = 1; i < cascadeEdges_[edgesPartner_[e]].size() - 1; ++i) {
        const Edge &edge = edges_[cascadeEdges_[edgesPartner_[e]][i]].e;
        cascadeSet.emplace(
          globalIndices[edge.first], globalIndices[edge.second]);
      }
    }
  }
}
