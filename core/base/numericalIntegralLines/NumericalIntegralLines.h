/// \ingroup base
/// \class ttk::NumericalIntegralLines
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date May 2026
/// \date NumericalIntegralLines extractor wrapping the DiscreteGradient class.
///
/// \brief TTK convenience class wrapping the DiscreteGradient class for
/// the easy extraction of vpaths.
///
/// Given a simplexId and dimension, this class returns a descending (or
/// ascending) vpath started in the given input simplex.
///
/// \sa NumericalIntegralLines.cpp %for an alternative integral line backend.
/// \sa DiscreteGradient.cpp %for the core mechanisms.
/// \sa ttkNumericalIntegralLines.cpp %for a usage example.
///

#pragma once

// base code includes
#include <Triangulation.h>
// std includes

namespace ttk {
  namespace nil {

    class NumericalIntegralLines : virtual public Debug {

    public:
      NumericalIntegralLines();
      ~NumericalIntegralLines() override;

      // template <class triangulationType = ttk::AbstractTriangulation>
      // int execute(triangulationType *triangulation);

      /**
       * @brief Compute a single numerical integral line.
       *
       * @param seed (SimplexId, dimension)
       * @param barycentricWeights Weights for the input seed.
       * @param output Output integral line (vector of 3D points).
       * @param isForawrd Forward or backward line (default: forward).
       */
      template <class dataType, class triangulationType>
        int computeIntegralLine(const triangulationType *triangulation,
          const std::pair<SimplexId, int> &seed,
          const std::array<float, 3> &barycentricWeights,
          std::vector<std::array<float, 3>> &output,
          const bool &isForward = false);

      /**
       * @brief Compute numerical integral lines.
       *
       * @param output Vector storing the output vpaths (1 entry per seed,
       * with possibly multiple v-path per seed).
       * @param isForward Forward or backward vpath (default: forward).
       */
      template <class dataType, class triangulationType>
      int execute(const triangulationType *triangulation,
        const std::vector<std::pair<SimplexId, int>> &seeds,
        std::vector<std::vector<std::array<float, 3>>> &output,
        const bool &isForward = false);

      /**
       * @brief Triangulation preconditioning.
       */
      inline void preconditionTriangulation(AbstractTriangulation *triangulation){
        // precondition simplex2face
        // precondition face2cofacets
      }

      inline void setInputScalarField(const void *const scalars){
        scalars_ = scalars;
      }

    protected:
      int maximumIterationNumber_{1000000000};
      const void *scalars_;
    };
  } // namespace nil
} // namespace ttk

template <class dataType, class triangulationType>
  int ttk::nil::NumericalIntegralLines::computeIntegralLine(
    const triangulationType *triangulation,
    const std::pair<SimplexId, int> &seed,
    const std::array<float, 3> &barycentricWeights,
    std::vector<std::array<float, 3>> &output,
    const bool &isForward){


  for(int i = 0; i < (int) maximumIterationNumber_; i++){

  }

  return 0;
}

template <class dataType, class triangulationType>
int ttk::nil::NumericalIntegralLines::execute(const triangulationType *triangulation,
    const std::vector<std::pair<SimplexId, int>> &seeds,
    std::vector<std::vector<std::array<float, 3>>> &output,
    const bool &isForward){

  Timer t;

  output.resize(seeds.size());

  const std::array<float, 3> barycentricWeights{1/3, 1/3, 1/3};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif
  for(int i = 0; i < (int) seeds.size(); i++){
    computeIntegralLine<dataType, triangulationType>(
      triangulation, seeds[i], barycentricWeights, output[i], isForward);

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
    printMsg("  - Seed-#"
      + std::to_string(seeds[i].first)
      + " (dim: "
      + std::to_string(seeds[i].second)
      + ", f: "
      + std::to_string(isForward)
      + "): "
      + std::to_string(output[i].size()) + " point(s).",
        debug::Priority::DETAIL);
  }

  printMsg("Computed numerical integral line(s) from "
    + std::to_string(output.size())
    + " seed(s)"
    , 1,
    t.getElapsedTime(), threadNumber_);

  return 0;
}
