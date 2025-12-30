/// \ingroup base
/// \class ttk::VectorWeightCurve
/// \author Tanner Finken
/// \author Joshua A. Levine
/// \date 2024.
///
/// \brief TTK processing package for the computation of weight curve
///         associated with the discrete vector field topology.
///
/// This package takes a 2D Vector Field as input and computes
/// the number of pairs as a function of weight (i.e. a number
/// similar to persistence).
///
/// These curves provide useful visual clues in order to fine-tune
/// simplification thresholds when displaying fields through
/// TopologicalSkeleton.
///
/// \sa ttk::TopologicalSkeleton
///
/// \b Related \b publication \n
/// "Localized Evaluation for Constructing Discrete Vector Fields" \n
/// Tanner Finken, Julien Tierny, Joshua A. Levine \n
/// IEEE VIS 2024.
///
/// \sa ttk::TopologicalSkeleton
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/discreteVectorFieldTopology/">Discrete
///   Vector Field Topology example</a> \n
///

#pragma once

// base code includes
#include <Debug.h>
#include <VectorSimplification.h>

namespace ttk {
  class VectorWeightCurve : virtual public Debug {
  public:
    VectorWeightCurve();
    VectorSimplification vs_{};

    bool DisplayOrbits{false};
    bool DisplayExtrema{false};
    bool ReverseFullOrbit{true};

    void UpdateTraceFullOrbits() {
      vs_.setFullOrbitSimplification(ReverseFullOrbit);
    }
  };
} // namespace ttk
