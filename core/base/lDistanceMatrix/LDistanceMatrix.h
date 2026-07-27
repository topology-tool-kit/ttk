/// \ingroup base
/// \class ttk::LDistanceMatrix
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date May 2020
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n

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

    template <typename TIn, typename TOut>
    int execute(std::vector<TOut *> &output,
                const std::vector<const TIn *> &inputs,
                const size_t nPoints) const;

    /// \warning This version of the execute function allocates, in addition to
    /// the output vector (of size inputs.size()^2) output, an extra vector of
    /// pointers of size output.size(). Given the size of the output 2D vector,
    /// this should be negligible.
    template <typename TIn, typename TOut>
    int execute(std::vector<std::vector<TOut>> &output,
                const std::vector<const TIn *> &inputs,
                const size_t nPoints) const;

  protected:
    std::string DistanceType{"2"};
  };
} // namespace ttk

template <typename TIn, typename TOut>
int ttk::LDistanceMatrix::execute(std::vector<TOut *> &output,
                                  const std::vector<const TIn *> &inputs,
                                  const size_t nPoints) const {

  const size_t nInputs = inputs.size();

  if(output.size() != nInputs)
    this->printErr(
      " When using the raw version of execute in LDistanceMatrix module, the "
      "output must be "
      "initialized, with the same number of lines as the number of inputs.");
  for(size_t i = 0; i < nInputs; i++)
    if(output[i] == nullptr)
      this->printErr(" When using the raw version of execute in "
                     "LDistanceMatrix module, the output must be "
                     "fully initialized: each line pointer must not be NULL.");

  LDistance worker{};
  worker.setThreadNumber(1);
  worker.setPrintRes(false);

  // compute matrix upper triangle
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) firstprivate(worker)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; ++i) {
    output[i][i] = 0;
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

template <typename TIn, typename TOut>
int ttk::LDistanceMatrix::execute(std::vector<std::vector<TOut>> &output,
                                  const std::vector<const TIn *> &inputs,
                                  const size_t nPoints) const {

  const auto nInputs = inputs.size();
  output.resize(nInputs);

  std::vector<double *> outputRaw(nInputs);
  for(size_t i = 0; i < nInputs; i++) {
    output[i].resize(nInputs);
    outputRaw[i] = output[i].data();
  }
  return execute(outputRaw, inputs, nPoints);
}
