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
    int execute(std::vector<std::vector<double>> &output,
                const std::vector<const T *> &inputs,
                const size_t nPoints) const;

  protected:
    std::string DistanceType{"2"};
  };
} // namespace ttk

template <typename T>
int ttk::LDistanceMatrix::execute(std::vector<std::vector<double>> &output,
                                  const std::vector<const T *> &inputs,
                                  const size_t nPoints) const {

  const auto nInputs = inputs.size();
  output.resize(nInputs);

  LDistance worker{};
  worker.setThreadNumber(1);
  worker.setPrintRes(false);

  for(size_t i = 0; i < nInputs; ++i) {
    output[i].resize(nInputs);
  }

  // compute matrix upper triangle
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) schedule(dynamic) \
  firstprivate(worker)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; ++i) {
    for(size_t j = i + 1; j < nInputs; ++j) {
      // call execute with nullptr output
      worker.execute(inputs[i], inputs[j], {}, this->DistanceType, nPoints);
      // store result
      output[i][j] = worker.getResult();
      output[j][i] = worker.getResult();
    }
  }

  return 0;
}
