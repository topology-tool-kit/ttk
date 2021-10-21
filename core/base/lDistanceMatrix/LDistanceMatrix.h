/// \ingroup base
/// \class ttk::LDistanceMatrix
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date May 2020

#pragma once

#include <LDistance.h>
#include <Wrapper.h>

#include <string>
#include <vector>

namespace ttk {
  class LDistanceMatrix : virtual public Debug {
  public:
    LDistanceMatrix();

    inline void setDistanceType(const std::string &val) {
      DistanceType = val;
    }

    template <typename T>
    std::vector<std::vector<double>> execute(const std::vector<void *> &inputs,
                                             const size_t nPoints) const;

  protected:
    std::string DistanceType{};
  };
} // namespace ttk

template <typename T>
std::vector<std::vector<double>>
  ttk::LDistanceMatrix::execute(const std::vector<void *> &inputs,
                                const size_t nPoints) const {

  const auto nInputs = inputs.size();
  std::vector<std::vector<double>> distMatrix(nInputs);

  LDistance worker{};
  worker.setThreadNumber(this->threadNumber_);

  // compute matrix upper triangle
  // (some parallelism inside the LDistance computation)
  for(size_t i = 0; i < nInputs; ++i) {
    auto &distCol = distMatrix[i];
    distCol.resize(nInputs);
    // get pointer to scalar field of input i
    const auto inputDataPtr1{static_cast<T *>(inputs[i])};
    for(size_t j = i + 1; j < nInputs; ++j) {
      // get pointer to scalar field of input jc
      const auto inputDataPtr2{static_cast<T *>(inputs[j])};
      // empty output
      T *outputPtr{};
      // call execute
      worker.execute(
        inputDataPtr1, inputDataPtr2, outputPtr, this->DistanceType, nPoints);
      // store result
      distMatrix[i][j] = worker.getResult();
    }
  }

  // distance matrix is symmetric
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; ++i) {
    for(size_t j = i + 1; j < nInputs; ++j) {
      distMatrix[j][i] = distMatrix[i][j];
    }
  }

  return distMatrix;
}
