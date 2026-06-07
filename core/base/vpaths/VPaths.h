/// \ingroup base
/// \class ttk::VPaths
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date May 2026
/// \date VPaths extractor wrapping the DiscreteGradient class.
///
/// \brief TTK convenience class wrapping the DiscreteGradient class for
/// the easy extraction of vpaths.
///
/// Given a simplexId and dimension, this class returns a descending (or
/// ascending) vpath started in the given input simplex.
///
/// \sa VPaths.cpp %for an alternative integral line backend.
/// \sa DiscreteGradient.cpp %for the core mechanisms.
/// \sa ttkVPaths.cpp %for a usage example.
///

#pragma once

// base code includes
#include <DiscreteGradient.h>
#include <Triangulation.h>
// std includes

namespace ttk {
  namespace vp {

    class VPaths : virtual public Debug {

    public:
      VPaths();
      ~VPaths() override;

      // template <class triangulationType = ttk::AbstractTriangulation>
      // int execute(triangulationType *triangulation);

      /**
       * @brief Extract a vpath.
       *
       * @param output Vector storing the output vpaths (1 entry per seed,
       * with possibly multiple v-path per seed).
       * @param isForward Forward or backward vpath (default: forward).
       */
      template <class triangulationType>
      int execute(
        const triangulationType *triangulation,
        const std::vector<ttk::dcg::Cell> &seeds,
        std::vector<std::vector<std::vector<ttk::dcg::Cell>>> &output,
        const bool &isForward = false);

      /**
       * @brief Triangulation preconditioning.
       */
      inline void preconditionTriangulation(AbstractTriangulation *triangulation){

        // see dms precondition
        dcg_.preconditionTriangulation(triangulation);
      }

      inline void setInputOffsets(const SimplexId *const offsets) {
        this->dcg_.setInputOffsets(offsets);
      }

      inline void setInputScalarField(const void *const scalars,
        const size_t &mTime){
        this->dcg_.setInputScalarField(scalars, mTime);
      }

    protected:

      dcg::DiscreteGradient dcg_{};
    };
  } // namespace vp
} // namespace ttk

template <class triangulationType>
int ttk::vp::VPaths::execute(
  const triangulationType *triangulation,
  const std::vector<dcg::Cell> &seeds,
  std::vector<std::vector<std::vector<dcg::Cell>>> &output,
  const bool &isForward){

  // fetching discrete gradient (or pre-computing it)
  dcg_.setDebugLevel(debugLevel_);
  dcg_.setThreadNumber(threadNumber_);
  dcg_.buildGradient(*triangulation, false, nullptr);

  Timer t;

  output.resize(seeds.size());

  /*
   * NOTE:
   * when considering seeds of non-zero dimension, mutliple v-paths may exist
   * for a given seed.
   */

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif
  for(int i = 0; i < (int) seeds.size(); i++){
    if(!isForward){
      dcg_.getAllDescendingPaths(seeds[i], output[i], *triangulation);
    }
    else{
      dcg_.getAllAscendingPaths(seeds[i], output[i], *triangulation);
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
    printMsg("  - Seed-#"
      + std::to_string(seeds[i].id_)
      + " (dim: "
      + std::to_string(seeds[i].dim_)
      + ", f: "
      + std::to_string(isForward)
      + "): "
      + std::to_string(output[i].size()) + " path(s).",
        debug::Priority::DETAIL);
  }

  printMsg("Computed v-path(s) from "
    + std::to_string(output.size()) + " seed(s)", 1,
    t.getElapsedTime(), threadNumber_);

  return 0;
}
