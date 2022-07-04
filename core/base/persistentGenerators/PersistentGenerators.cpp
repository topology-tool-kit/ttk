#include <PersistentGenerators.h>
#include <UnionFind.h>

ttk::PersistentGenerators::PersistentGenerators() {
  this->setDebugMsgPrefix("PersistentGenerators");
}

void ttk::PersistentGenerators::getConnectedComponents(
  const std::vector<std::array<SimplexId, 2>> &edgeNeighs,
  std::vector<SimplexId> &connComp) const {

  // use Union-Find to get connected components
  std::vector<ttk::UnionFind> uf(connComp.size());

  // initialize UFs
  for(size_t j = 0; j < edgeNeighs.size(); ++j) {
    uf[j].setRank(j);
  }

  // Union between adjacent edges
  for(size_t j = 0; j < edgeNeighs.size(); ++j) {
    ttk::UnionFind::makeUnion(&uf[j], &uf[edgeNeighs[j][0]]);
    ttk::UnionFind::makeUnion(&uf[j], &uf[edgeNeighs[j][1]]);
  }

  // UF rank -> componend id mapping
  std::vector<SimplexId> connCompIds(edgeNeighs.size(), -1);
  size_t nConnComps{};

  // find connected component ids
  for(size_t j = 0; j < edgeNeighs.size(); ++j) {
    const auto root{uf[j].find()};
    if(root != nullptr) {
      const auto rank{root->getRank()};
      if(connCompIds[rank] == -1) {
        connCompIds[rank] = nConnComps++;
      }
      connComp[j] = connCompIds[rank];
    }
  }
}
