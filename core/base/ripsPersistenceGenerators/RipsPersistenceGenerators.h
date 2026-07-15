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
                 std::vector<rpd::Generator1> &generators) const {
      rpd::PairCellsWithOracle::callOracle(
        points, diagrams, SimplexMaximumDiameter, InputIsDistanceMatrix);
      rpd::PairCellsWithOracle pc(
        points, diagrams, InputIsDistanceMatrix, false);
      pc.setDebugLevel(debugLevel_);
      pc.run();
      if(!OutputCascade)
        pc.getGenerators(generators);
      else {
        rpd::EdgeSets4 criticalAndCascade;
        pc.getCascades(criticalAndCascade);
        generators.emplace_back(
          criticalAndCascade[rpd::DEATH0], std::make_pair(0., 0.));
        generators.emplace_back(
          criticalAndCascade[rpd::BIRTH1], std::make_pair(1., 1.));
        generators.emplace_back(
          criticalAndCascade[rpd::DEATH1], std::make_pair(2., 2.));
        generators.emplace_back(
          criticalAndCascade[rpd::CASC1], std::make_pair(3., 3.));
      }
    }

  protected:
    /** Rips diameter threshold */
    double SimplexMaximumDiameter{rpd::inf};
    /** is input a distance matrix */
    bool InputIsDistanceMatrix{false};
    /** output cascade */
    bool OutputCascade{false};
  };

} // namespace ttk