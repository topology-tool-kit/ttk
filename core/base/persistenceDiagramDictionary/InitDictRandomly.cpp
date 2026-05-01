#include <InitDictRandomly.h>
#include <Shuffle.h>

#include <numeric>
#include <random>

using namespace ttk;

void InitRandomDict::execute(std::vector<ttk::DiagramType> &DictDiagrams,
                             const std::vector<ttk::DiagramType> &datas,
                             const int nbAtom,
                             const int seed) {
  int nDiags = datas.size();
  DictDiagrams.resize(nbAtom);
  std::vector<int> indices(nDiags);
  std::iota(indices.begin(), indices.end(), 0);
  std::mt19937 random_engine{};
  random_engine.seed(seed);
  ttk::shuffle(indices, random_engine);
  for(int i = 0; i < nbAtom; ++i) {
    const auto &atom = datas[indices[i]];
    DictDiagrams[i] = atom;
  }
}
