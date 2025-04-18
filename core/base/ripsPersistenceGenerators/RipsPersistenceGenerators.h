/// \ingroup base
/// \class ttk::RipsPersistenceGenerators
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date June 2024.
///
/// \brief TTK base class that computes 1-dimensional persistence generators in
/// a Rips filtration
///
/// This module defines the %RipsPersistenceGenerators that takes a point
/// cloud and computes 1-dimensional persistence generators, in addition to the
/// associated persistence diagram, of its Rips filtration. It can also compute
/// topologically critical edges and persistent cascades of this Rips filtration
///
/// \sa ttkRipsPersistenceGenerators.cpp %for a usage example.

#pragma once

// ttk common includes
#include <Debug.h>

#include <PairCellsWithOracle.h>

namespace ttk {

  class RipsPersistenceGenerators : virtual public Debug {
  public:
    RipsPersistenceGenerators();
    void execute(const std::vector<std::vector<double>> &points,
                 rpd::MultidimensionalDiagram &diagrams,
                 std::vector<rpd::Generator> &generators) const {
      rpd::PairCellsWithOracle::callOracle(points, diagrams, SimplexMaximumDiameter);
      rpd::PairCellsWithOracle pc(points, diagrams, false, false);
      pc.setDebugLevel(debugLevel_);
      pc.run();
      if (!OutputCascade)
        pc.getGenerators(generators);
      else {
        rpd::EdgeSets4 criticalAndCascade;
        pc.getCascades(criticalAndCascade);
        generators.emplace_back(criticalAndCascade[0], std::make_pair(0., 0.)); // MST
        generators.emplace_back(criticalAndCascade[1], std::make_pair(1., 1.)); // RNG
        generators.emplace_back(criticalAndCascade[2], std::make_pair(2., 2.)); // MML
        generators.emplace_back(criticalAndCascade[3], std::make_pair(3., 3.)); // cascade
      }
    }
  protected:
    double SimplexMaximumDiameter {rpd::inf};
    bool OutputCascade{false};
  };

}