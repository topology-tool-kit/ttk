#include <MorseSmaleComplex.h>

ttk::MorseSmaleComplex::MorseSmaleComplex() {
  this->setDebugMsgPrefix("MorseSmaleComplex");
}

void ttk::MorseSmaleComplex::flattenSeparatricesVectors(
  std::vector<std::vector<Separatrix>> &separatrices,
  std::vector<std::vector<std::vector<ttk::dcg::Cell>>> &separatricesGeometry)
  const {

  std::vector<size_t> partialSizes{0};
  for(const auto &sep : separatrices) {
    partialSizes.emplace_back(partialSizes.back() + sep.size());
  }
  separatrices[0].resize(partialSizes.back());
  separatricesGeometry[0].resize(partialSizes.back());

  for(size_t i = 1; i < separatrices.size(); ++i) {
    const auto offset = partialSizes[i];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t j = 0; j < separatrices[i].size(); ++j) {
      // shift separatrices geometry_
      for(size_t k = 0; k < separatrices[i][j].geometry_.size(); ++k) {
        separatrices[i][j].geometry_[k] += offset;
      }
      // flatten separatrices1 and separatricesGeometry1
      separatrices[0][offset + j] = std::move(separatrices[i][j]);
      separatricesGeometry[0][offset + j]
        = std::move(separatricesGeometry[i][j]);
    }
  }
}
