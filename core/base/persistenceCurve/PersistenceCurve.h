/// \ingroup base
/// \class ttk::PersistenceCurve
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date September 2016.
///
/// \brief TTK processing package for the computation of persistence curves.
///
/// This package takes a \ref ttk::DiagramType as input and computes
/// the number of pairs as a function of persistence (i.e. the number
/// of pairs whose persistence is higher than a threshold).
///
/// These curves provide useful visual clues in order to fine-tune persistence
/// simplification thresholds.
///
/// \sa ttkPersistenceCurve.cpp %for a usage example.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/BuiltInExample1/">BuiltInExample1
/// </a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/dragon/">Dragon
/// example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morsePersistence/">Morse
///   Persistence example</a> \n

#pragma once

// base code includes
#include <Debug.h>
#include <PersistenceDiagramUtils.h>

namespace ttk {
  class PersistenceCurve : virtual public Debug {
  public:
    PersistenceCurve();

    /**
     * Plot of (Persistence, Number of Pairs)
     */
    using PlotType = std::vector<std::pair<double, SimplexId>>;

    /**
     * @brief Compute the Persistence Curve from the input Diagram
     *
     * @param[out] plots Array of 4 PlotTypes: minimum-saddle plot,
     * saddle-saddle plot, saddle-maximum plot and all pairs plot
     * @param[in] diagram Input Persistence Diagram
     *
     * @return 0 in case of success
     */
    int execute(std::array<PlotType, 4> &plots,
                const DiagramType &diagram) const;
  };
} // namespace ttk
