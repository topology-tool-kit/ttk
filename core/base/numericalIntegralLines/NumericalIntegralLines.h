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
  namespace vp {

    class NumericalIntegralLines : virtual public Debug {

    public:
      NumericalIntegralLines();
      ~NumericalIntegralLines() override;

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
        const bool &isForward = false);

      /**
       * @brief Triangulation preconditioning.
       */
      inline void preconditionTriangulation(AbstractTriangulation *triangulation){

        // see dms precondition
      }

      inline void setInputOffsets(const SimplexId *const offsets) {
      }

      inline void setInputScalarField(const void *const scalars,
        const size_t &mTime){
      }

    protected:
    };
  } // namespace vp
} // namespace ttk

template <class triangulationType>
int ttk::vp::NumericalIntegralLines::execute(
  const triangulationType *triangulation,
  const bool &isForward){

  Timer t;

  printMsg("Computed numerical integral line(s) from "
    // + std::to_string(output.size())
    // + " seed(s)"
    , 1,
    t.getElapsedTime(), threadNumber_);

  return 0;
}
