/// \ingroup baseCode
/// \class ttk::PersistentGenerators
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2021.
///
/// \brief TTK %PersistentGenerators processing package.
///
/// ttk::PersistentGenerators uses ttk::discreteMorseSandwich to
/// compute generators for persistence pairs of dimensions 1.
///
/// \b Related \b publication \n
/// "Discrete Morse Sandwich: Fast Computation of Persistence Diagrams for
/// Scalar Data -- An Algorithm and A Benchmark" \n
/// Pierre Guillou, Jules Vidal, Julien Tierny \n
/// IEEE Transactions on Visualization and Computer Graphics, 2023.\n
/// arXiv:2206.13932, 2023.
///
/// \sa ttk::DiscreteMorseSandwich
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_at/">Persistent
///   Generators AT example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_casting/">Persistent
///   Generators Casting example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_darkSky/">Persistent
///   Generators DarkSky example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_fertility/">Persistent
///   Generators Fertility example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_householdAnalysis/">Persistent
///   Generators Household Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_periodicPicture/">Persistent
///   Generators Periodic Picture example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_skull/">Persistent
///   Generators Skull example</a> \n
///

#pragma once

#include <DiscreteMorseSandwich.h>

#include <limits>

namespace ttk {
  class PersistentGenerators : virtual public DiscreteMorseSandwich {
  public:
    PersistentGenerators();

    /**
     * @brief Compute the persistence generators from the discrete gradient
     *
     * @pre @ref buildGradient and @ref preconditionTriangulation
     * should be called prior to this function
     *
     * @param[out] generators Persistent generators
     * @param[out] connComps Generators connected components
     * @param[in] offsets Order field
     * @param[in] triangulation Preconditionned triangulation
     *
     * @return 0 when success
     */
    template <typename triangulationType>
    int computePersistentGenerators(
      std::vector<GeneratorType> &generators,
      std::vector<std::vector<SimplexId>> &connComps,
      const SimplexId *const offsets,
      const triangulationType &triangulation);

  private:
    /**
     * @brief Use Union-Find to get connected components
     *
     * @param[in] edgeNeighs Graph of edge neighbors
     * @param[out] connComp Edges connected component id
     */
    void getConnectedComponents(
      const std::vector<std::array<SimplexId, 2>> &edgeNeighs,
      std::vector<SimplexId> &connComp) const;

    /**
     * @brief Detect connected components of persistent generators
     *
     * @param[out] connComps Connected component id per edge
     * @param[in] generators Pre-computed persistent generators
     * @param[in] triangulation Triangulation
     */
    template <typename triangulationType>
    void findGeneratorsConnectedComponents(
      std::vector<std::vector<SimplexId>> &connComps,
      const std::vector<GeneratorType> &generators,
      const triangulationType &triangulation) const;

    /**
     * @brief Find generators corresponding to infinite pairs of
     * dimension 1 (topological handles)
     *
     * @param[out] generators Persistent generators to compute
     * @param[in] inf1saddles Unpaired 1-saddles (infinite pairs)
     * @param[in] minSadPairs Propagate on simplified topology
     * @param[in] offsets Order field
     * @param[in] triangulation Triangulation
     */
    template <typename triangulationType>
    void findHandlesGenerators(std::vector<GeneratorType> &generators,
                               const std::vector<SimplexId> &inf1Saddles,
                               const std::vector<PersistencePair> &minSadPairs,
                               const SimplexId *const offsets,
                               const triangulationType &triangulation) const;

    /**
     * @brief Clean a topogical handle generator by removing extra branches
     *
     * Performs a Dijkstra shortest path between the two vertices of
     * the critical edge marking the topogical handle
     *
     * @param[in,out] generator Generator to be pruned
     * @param[in,out] genMask Mask generator edges
     * @param[in,out] pathVerts Mask shortest path vertices
     * @param[in,out] preds Vertex predecessor (for Dijkstra)
     * @param[in,out] preds Vertex distance to critical edge (for Dijkstra)
     * @param[in] s1 Critical edge of the topogical handle
     * @param[in] triangulation Triangulation
     */
    template <typename triangulationType>
    void pruneHandleGenerator(std::vector<SimplexId> &generator,
                              std::vector<bool> &genMask,
                              std::vector<bool> &pathVerts,
                              std::vector<SimplexId> &preds,
                              std::vector<size_t> &dists,
                              const SimplexId &s1,
                              const triangulationType &triangulation) const;

  protected:
    bool PruneHandlesGenerators{false};
  };
} // namespace ttk

template <typename triangulationType>
void ttk::PersistentGenerators::findGeneratorsConnectedComponents(
  std::vector<std::vector<SimplexId>> &connComps,
  const std::vector<GeneratorType> &generators,
  const triangulationType &triangulation) const {

  connComps.resize(generators.size());

  // global edge id -> local generator id mapping
  std::vector<SimplexId> isVisited(triangulation.getNumberOfEdges(), -1);

  for(size_t i = 0; i < generators.size(); ++i) {
    const auto &genEdges{generators[i].boundary};
    auto &connComp{connComps[i]};
    if(genEdges.empty()) {
      continue;
    }
    connComp.resize(genEdges.size(), -1);
    // fill isVisited for the local generator
    // (cleanup at iteration end)
    for(size_t j = 0; j < genEdges.size(); ++j) {
      isVisited[genEdges[j]] = j;
    }

    // for each generator edge, the local neighbor ids
    // (assumption: exactly 2 neighboring edges per generator edge)
    std::vector<std::array<SimplexId, 2>> edgeNeighs(genEdges.size());

    for(size_t j = 0; j < genEdges.size(); ++j) {
      const auto edge{genEdges[j]};
      for(size_t k = 0; k < 2; ++k) {
        SimplexId v{};
        triangulation.getEdgeVertex(edge, k, v);
        const auto nneighs{triangulation.getVertexEdgeNumber(v)};
        for(SimplexId l = 0; l < nneighs; ++l) {
          SimplexId neigh{};
          triangulation.getVertexEdge(v, l, neigh);
          const auto neighId{isVisited[neigh]};
          if(neigh == edge || neighId == -1) {
            continue;
          }
          edgeNeighs[j][k] = neighId;
          break; // go to next vertex
        }
      }
    }

    this->getConnectedComponents(edgeNeighs, connComp);

    // clean isVisited
    for(const auto e : genEdges) {
      isVisited[e] = -1;
    }
  }
}

template <typename triangulationType>
void ttk::PersistentGenerators::findHandlesGenerators(
  std::vector<GeneratorType> &generators,
  const std::vector<SimplexId> &inf1Saddles,
  const std::vector<PersistencePair> &minSadPairs,
  const SimplexId *const offsets,
  const triangulationType &triangulation) const {

  Timer tm{};

  // minimum vertex id -> paired 1-saddle edge id mapping
  std::vector<SimplexId> min2sad(triangulation.getNumberOfVertices(), -1);
  for(const auto &p : minSadPairs) {
    min2sad[p.birth] = p.death;
  }

  // global maximum
  const auto globMax{std::distance(
    offsets,
    std::max_element(offsets, offsets + triangulation.getNumberOfVertices()))};

  // storage (allocated once, cleaned after each iteration)
  std::vector<bool> genMask(triangulation.getNumberOfEdges(), false);
  std::vector<size_t> dists(
    triangulation.getNumberOfVertices(), std::numeric_limits<size_t>::max());
  std::vector<SimplexId> preds(triangulation.getNumberOfVertices(), -1);
  std::vector<bool> pathVerts(triangulation.getNumberOfVertices(), false);

  // follow descending 1-separatrices to find representative generator
  for(const auto s1 : inf1Saddles) {
    std::vector<SimplexId> generator{};
    std::set<SimplexId> visited1Sads{};

    std::queue<SimplexId> toPropagate{};
    toPropagate.push(s1);

    while(!toPropagate.empty()) {
      // current 1-saddle
      const auto curr{toPropagate.front()};
      toPropagate.pop();
      visited1Sads.emplace(curr);

      if(curr != s1) {
        generator.emplace_back(curr);
      }

      for(SimplexId i = 0; i < 2; ++i) {
        SimplexId currVert{};
        triangulation.getEdgeVertex(curr, i, currVert);
        std::vector<Cell> vpath{};
        SimplexId min{-1};
        this->dg_.getDescendingPath(Cell{0, currVert}, vpath, triangulation);
        const auto &lastCell = vpath.back();
        if(lastCell.dim_ == 0) {
          min = lastCell.id_;
        }

        if(min == -1) {
          // separatrix leads to nowhere
          continue;
        }

        // store vpath
        for(const auto &c : vpath) {
          if(c.dim_ == 1) {
            generator.emplace_back(c.id_);
          }
        }

        if(min2sad[min] == -1) {
          // global minimum ?
          continue;
        }

        const auto next{min2sad[min]};
        if(visited1Sads.find(next) == visited1Sads.end()) {
          toPropagate.push(next);
        }
      }
    }

    // remove duplicates
    TTK_PSORT(this->threadNumber_, generator.begin(), generator.end());
    const auto last = std::unique(generator.begin(), generator.end());
    generator.erase(last, generator.end());

    if(this->PruneHandlesGenerators) {
      // post-process: Dijkstra from one 1-saddle vertex to the other
      this->pruneHandleGenerator(
        generator, genMask, pathVerts, preds, dists, s1, triangulation);
    }

    // ensure unpaired 1-saddle is first
    generator.emplace_back(s1);
    std::swap(generator[0], generator.back());

    generators.emplace_back(GeneratorType{
      generator, -1,
      std::array<SimplexId, 2>{
        static_cast<SimplexId>(globMax),
        this->dg_.getCellGreaterVertex(Cell{1, s1}, triangulation),
      }});
  }

  this->printMsg("Computed " + std::to_string(inf1Saddles.size())
                   + " topological handle generators",
                 1.0, tm.getElapsedTime(), this->threadNumber_);
}

template <typename triangulationType>
void ttk::PersistentGenerators::pruneHandleGenerator(
  std::vector<SimplexId> &generator,
  std::vector<bool> &genMask,
  std::vector<bool> &pathVerts,
  std::vector<SimplexId> &preds,
  std::vector<size_t> &dists,
  const SimplexId &s1,
  const triangulationType &triangulation) const {

  // store generator vertices for faster cleanup
  std::vector<SimplexId> genVerts{};
  genVerts.reserve(2 * generator.size());

  for(const auto e : generator) {
    genMask[e] = true;
  }

  using PQueueElem = std::pair<size_t, SimplexId>;
  using PQueue = std::priority_queue<PQueueElem, std::vector<PQueueElem>,
                                     std::greater<PQueueElem>>;
  PQueue pq{};
  // try to find the shortest path between seed0 & seed1 (the
  // vertices of the critical edge marking the topological handle)
  // on generator \ {s1}
  SimplexId seed0{}, seed1{};
  triangulation.getEdgeVertex(s1, 0, seed0);
  triangulation.getEdgeVertex(s1, 1, seed1);
  genVerts.emplace_back(seed0);
  genVerts.emplace_back(seed1);
  pq.emplace(0, seed0);
  dists[seed0] = 0;

  while(!pq.empty()) {
    const auto curr{pq.top()};
    genVerts.emplace_back(curr.second);
    pq.pop();

    const auto nneighs{triangulation.getVertexNeighborNumber(curr.second)};
    dists[curr.second] = curr.first;

    for(SimplexId i = 0; i < nneighs; ++i) {
      SimplexId neigh{}, edge{};
      // vertex neighbor & corresponding edge are at the same index
      triangulation.getVertexNeighbor(curr.second, i, neigh);
      triangulation.getVertexEdge(curr.second, i, edge);
      // stay on generator edges
      if(!genMask[edge]) {
        continue;
      }
      // skip critical edge
      if(curr.second == seed0 && neigh == seed1) {
        continue;
      }
      if(dists[neigh] > curr.first + 1) {
        dists[neigh] = curr.first + 1;
        preds[neigh] = curr.second;
        pq.emplace(dists[neigh], neigh);
      }
    }
  }

  // extract shortest path
  SimplexId curr{seed1};
  while(curr != seed0) {
    pathVerts[curr] = true;
    curr = preds[curr];
  }
  pathVerts[seed0] = true;
  pathVerts[seed1] = true;

  std::vector<SimplexId> prunedGenerator{};
  prunedGenerator.reserve(generator.size());
  for(const auto e : generator) {
    SimplexId v0{}, v1{};
    triangulation.getEdgeVertex(e, 0, v0);
    triangulation.getEdgeVertex(e, 1, v1);
    if(pathVerts[v0] && pathVerts[v1]) {
      prunedGenerator.emplace_back(e);
    }
  }
  std::swap(generator, prunedGenerator);

  // cleanup: reset large vectors
  for(const auto e : prunedGenerator) {
    genMask[e] = false;
  }
  for(const auto v : genVerts) {
    preds[v] = -1;
    dists[v] = std::numeric_limits<size_t>::max();
    pathVerts[v] = false;
  }
}

template <typename triangulationType>
int ttk::PersistentGenerators::computePersistentGenerators(
  std::vector<GeneratorType> &generators,
  std::vector<std::vector<SimplexId>> &connComps,
  const SimplexId *const offsets,
  const triangulationType &triangulation) {

  // allocate memory
  this->alloc(triangulation);

  Timer tm{};
  const auto dim = this->dg_.getDimensionality();
  if(dim < 2) {
    this->printWrn("Cannot compute cycles for 1D datasets");
    return 0;
  }
  if(dim == 2) {
    // fix container size for 2D datasets (use
    // eliminateBoundariesSandwich for saddle-max instead of
    // tripletsToPersistencePairs)
    this->critEdges_.resize(triangulation.getNumberOfEdges());
    this->edgeTrianglePartner_.resize(triangulation.getNumberOfEdges(), -1);
    this->onBoundary_.resize(triangulation.getNumberOfEdges(), false);
    this->s2Mapping_.resize(triangulation.getNumberOfTriangles(), -1);
    this->s1Mapping_.resize(triangulation.getNumberOfEdges(), -1);
  }

  // get every critical cell sorted them by dimension
  std::array<std::vector<SimplexId>, 4> criticalCellsByDim{};
  // holds the critical cells order
  auto &critCellsOrder{this->critCellsOrder_};

  this->extractCriticalCells(
    criticalCellsByDim, critCellsOrder, offsets, triangulation, true);

  // if minima are paired
  auto &pairedMinima{this->pairedCritCells_[0]};
  // if 1-saddles are paired
  auto &paired1Saddles{this->pairedCritCells_[1]};
  // if 2-saddles are paired
  auto &paired2Saddles{this->pairedCritCells_[dim - 1]};
  // if maxima are paired
  auto &pairedMaxima{this->pairedCritCells_[dim]};

  // minima - saddle pairs
  // we need the min-saddle pairs to get the infinite pairs of
  // dimension 1 (topological handles)
  std::vector<PersistencePair> minSadPairs{};
  this->getMinSaddlePairs(minSadPairs, pairedMinima, paired1Saddles,
                          criticalCellsByDim[1], critCellsOrder[1], offsets,
                          triangulation);

  if(dim == 3) {
    // saddle - maxima pairs
    // this should fasten the computation of saddle-saddle pairs
    // (sandwich)
    std::vector<PersistencePair> sadMaxPairs{};
    this->getMaxSaddlePairs(
      sadMaxPairs, pairedMaxima, paired2Saddles, criticalCellsByDim[dim - 1],
      critCellsOrder[dim - 1], critCellsOrder[dim], triangulation);
  }

  if(!criticalCellsByDim[1].empty() && !criticalCellsByDim[2].empty()) {
    // compute saddle-saddle pairs, extract generators/cycles
    // (eliminateBoundariesSandwich boundaries right before
    // simplification)
    std::vector<PersistencePair> sadSadPairs{};
    this->getSaddleSaddlePairs(
      sadSadPairs, paired1Saddles, dim == 3 ? paired2Saddles : pairedMaxima,
      true, generators, criticalCellsByDim[1], criticalCellsByDim[2],
      critCellsOrder[1], triangulation);
  }

  // detect topological handles
  std::vector<SimplexId> inf1Saddles{};
  for(const auto s1 : criticalCellsByDim[1]) {
    if(!paired1Saddles[s1]) {
      inf1Saddles.emplace_back(s1);
    }
  }

  if(!inf1Saddles.empty()) {
    this->findHandlesGenerators(
      generators, inf1Saddles, minSadPairs, offsets, triangulation);
  }

  if(!generators.empty()) {
    // detect connected components in cycles
    this->findGeneratorsConnectedComponents(
      connComps, generators, triangulation);
  }

  this->printMsg(
    "Computed " + std::to_string(generators.size()) + " persistent generators",
    1.0, tm.getElapsedTime(), this->threadNumber_);

  // free memory
  this->clear();

  return 0;
}
