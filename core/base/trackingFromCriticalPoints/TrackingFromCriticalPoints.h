/// \ingroup base
/// \class ttk::TrackingFromPersistenceDiagrams
/// \author Thomas Daniel <thomas.daniel@lip6.fr>
/// \date January 2025.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/trackingFromCriticalPoints/">Tracking
///   From Critical Points example</a>

#pragma once

// base code includes
#include <DataTypes.h>
#include <PersistenceDiagramUtils.h>
#include <Triangulation.h>

namespace ttk {

  using trackingTuple = std::tuple<int, int, std::vector<SimplexId>>;

  class TrackingFromCriticalPoints : virtual public Debug {

  private:
    double epsilonConstant_{10e-1};
    double epsilonAdapt_{0.5};
    double meshDiameter_{1};
    double tolerance_{10e-3};
    int assignmentMethod_{0};
    double xWeight_{1};
    double yWeight_{1};
    double zWeight_{1};
    double fWeight_{0};
    bool adaptiveDeathBirthCost_{false};

  public:
    TrackingFromCriticalPoints() = default;

    void setMeshDiameter(double r) {
      meshDiameter_ = r;
    }

    void setEpsilon(double e) {
      epsilonConstant_ = e;
    }

    void setEpsilonAdapt(double e) {
      epsilonAdapt_ = e;
    }

    void setTolerance(double t) {
      tolerance_ = t;
    }

    void setAssignmentMethod(int a) {
      if(a == 0 || a == 1) {
        assignmentMethod_ = a;
      }
    }

    void setAdaptDeathBirthCost(bool b) {
      adaptiveDeathBirthCost_ = b;
    }

    void setWeights(double PX, double PY, double PZ, double PF) {
      xWeight_ = PX;
      yWeight_ = PY;
      zWeight_ = PZ;
      fWeight_ = PF;
    }

    double computeBoundingBoxRadius(const DiagramType &d1,
                                    const DiagramType &d2) {
      double maxScalar = d1[0].birth.sfValue;
      double minScalar = d1[0].birth.sfValue;

      for(unsigned int i = 0; i < d1.size(); i++) {
        maxScalar = std::max(maxScalar, d1[i].birth.sfValue);
        maxScalar = std::max(maxScalar, d1[i].death.sfValue);
        minScalar = std::min(minScalar, d1[i].birth.sfValue);
        minScalar = std::min(minScalar, d1[i].death.sfValue);
      }

      for(unsigned int i = 0; i < d2.size(); i++) {
        maxScalar = std::max(maxScalar, d2[i].birth.sfValue);
        maxScalar = std::max(maxScalar, d2[i].death.sfValue);
        minScalar = std::min(minScalar, d2[i].birth.sfValue);
        minScalar = std::min(minScalar, d2[i].death.sfValue);
      }

      return std::sqrt(std::pow(meshDiameter_, 2)
                       + fWeight_ * std::pow(maxScalar - minScalar, 2));
    }

    void
      performMatchings(const std::vector<DiagramType> &persistenceDiagrams,
                       std::vector<std::vector<MatchingType>> &maximaMatchings,
                       std::vector<std::vector<MatchingType>> &sad_1_Matchings,
                       std::vector<std::vector<MatchingType>> &sad_2_Matchings,
                       std::vector<std::vector<MatchingType>> &minimaMatchings,
                       std::vector<std::vector<SimplexId>> &maxMap,
                       std::vector<std::vector<SimplexId>> &sad_1Map,
                       std::vector<std::vector<SimplexId>> &sad_2Map,
                       std::vector<std::vector<SimplexId>> &minMap);
    void performTrackings(
      const std::vector<DiagramType> &persistenceDiagrams,
      const std::vector<std::vector<MatchingType>> &maximaMatchings,
      const std::vector<std::vector<MatchingType>> &sad_1_Matchings,
      const std::vector<std::vector<MatchingType>> &sad_2_Matchings,
      const std::vector<std::vector<MatchingType>> &minimaMatchings,
      const std::vector<std::vector<SimplexId>> &maxMap,
      const std::vector<std::vector<SimplexId>> &sad_1Map,
      const std::vector<std::vector<SimplexId>> &sad_2Map,
      const std::vector<std::vector<SimplexId>> &minMap,
      std::vector<trackingTuple> &allTrackings,
      std::vector<std::vector<double>> &allTrackingsCost,
      std::vector<double> &allTrackingsMeanPersistences,
      unsigned int (&typesArrayLimits)[3]);

  private:
    double computeRelevantPersistence(const DiagramType &d1,
                                      const DiagramType &d2) {
      const auto sp = this->tolerance_;
      const double s = sp > 0.0 && sp < 100.0 ? sp / 100.0 : 0;

      std::vector<double> toSort(d1.size() + d2.size());
      for(size_t i = 0; i < d1.size(); ++i) {
        const auto &t = d1[i];
        toSort[i] = std::abs(t.persistence());
      }
      for(size_t i = 0; i < d2.size(); ++i) {
        const auto &t = d2[i];
        toSort[d1.size() + i] = std::abs(t.persistence());
      }

      const auto minVal = *std::min_element(toSort.begin(), toSort.end());
      const auto maxVal = *std::max_element(toSort.begin(), toSort.end());
      return s * (maxVal - minVal);
    }

    // Compute L_p distance betweem (p,f(p)) and (q,f(q)) where p and q are
    // critical points

    double criticalPointDistance(const std::array<float, 3> &coords_p1,
                                 const double &sfValue_p1,
                                 const std::array<float, 3> &coords_p2,
                                 const double &sfValue_p2,
                                 const int &p);

    // Sort the critical points by types

    void sortCriticalPoint(const DiagramType &d,
                           const double minimumRelevantPersistence,
                           std::vector<std::array<float, 3>> &maxCoords,
                           std::vector<std::array<float, 3>> &sad_1Coords,
                           std::vector<std::array<float, 3>> &sad_2Coords,
                           std::vector<std::array<float, 3>> &minCoords,
                           std::vector<double> &maxScalar,
                           std::vector<double> &sad1Scalar,
                           std::vector<double> &sad_2Scalar,
                           std::vector<double> &minScalar,
                           std::vector<SimplexId> &mapMax,
                           std::vector<SimplexId> &mapSad_1,
                           std::vector<SimplexId> &mapSad_2,
                           std::vector<SimplexId> &mapMin);

    void buildCostMatrix(const std::vector<std::array<float, 3>> &coords_1,
                         const std::vector<double> &sfValues_1,
                         const std::vector<std::array<float, 3>> &coords_2,
                         const std::vector<double> &sfValues_2,
                         const float &costDeathBirth,
                         std::vector<std::vector<double>> &matrix);

    void localToGlobalMatching(const std::vector<int> &startMap,
                               const std::vector<int> &endMap,
                               const std::vector<double> &startPersistence,
                               const std::vector<double> &endPersistence,
                               std::vector<MatchingType> &matchings,
                               std::vector<MatchingType> &matchingsPersistence);

    void assignmentSolver(std::vector<std::vector<double>> &costMatrix,
                          std::vector<ttk::MatchingType> &matching);

    int computeGlobalId(const DiagramType &persistenceDiagram,
                        const CriticalType &type,
                        const SimplexId &id);

    void performTrackingForOneType(
      const std::vector<DiagramType> &persistenceDiagrams,
      const std::vector<std::vector<MatchingType>> &matching,
      const std::vector<std::vector<SimplexId>> &map,
      const CriticalType &currentType,
      std::vector<trackingTuple> &tracking,
      std::vector<std::vector<double>> &trackingCosts,
      std::vector<double> &trackingPersistence);
  };
} // namespace ttk
