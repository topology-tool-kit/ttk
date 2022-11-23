#include <MorseSmaleComplex.h>

ttk::MorseSmaleComplex::MorseSmaleComplex() {
  this->setDebugMsgPrefix("MorseSmaleComplex");
}

void ttk::MorseSmaleComplex::flattenSeparatricesVectors(
  std::vector<std::vector<Separatrix>> &separatrices) const {

  std::vector<size_t> partialSizes{0};
  for(const auto &sep : separatrices) {
    partialSizes.emplace_back(partialSizes.back() + sep.size());
  }
  separatrices[0].resize(partialSizes.back());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 1; i < separatrices.size(); ++i) {
    for(size_t j = 0; j < separatrices[i].size(); ++j) {
      const auto o = partialSizes[i] + j;
      // flatten separatrices1 and separatricesGeometry1
      separatrices[0][o].source_ = separatrices[i][j].source_;
      separatrices[0][o].destination_ = separatrices[i][j].destination_;
      separatrices[0][o].geometry_ = std::move(separatrices[i][j].geometry_);
    }
  }
}
