/// \ingroup base
/// \class ttk::BottleneckDistance
/// \author Maxime Soler <soler.maxime@total.com>
#pragma once

// base code includes
#include <Debug.h>
#include <PersistenceDiagramUtils.h>

namespace ttk {
  constexpr unsigned long long str2int(const char *str, int h = 0) {
    return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
  }

  class BottleneckDistance : virtual public Debug {
  public:
    BottleneckDistance();

    int execute(const ttk::DiagramType &diag0,
                const ttk::DiagramType &diag1,
                std::vector<MatchingType> &matchings);

    inline void setPersistencePercentThreshold(const double t) {
      Tolerance = t;
    }
    inline void setPX(const double px) {
      PX = px;
    }
    inline void setPY(const double py) {
      PY = py;
    }
    inline void setPZ(const double pz) {
      PZ = pz;
    }
    inline void setPE(const double pe) {
      PE = pe;
    }
    inline void setPS(const double ps) {
      PS = ps;
    }
    inline void setAlgorithm(const std::string &algorithm) {
      DistanceAlgorithm = algorithm;
    }
    inline void setPVAlgorithm(const int algorithm) {
      PVAlgorithm = algorithm;
    }
    inline void setWasserstein(const std::string &wasserstein) {
      WassersteinMetric = wasserstein;
    }

    double getDistance() {
      return distance_;
    }

  protected:
    double distance_{-1.0};
    std::array<double, 3> costs_{};

    std::string WassersteinMetric{"2"};
    std::string DistanceAlgorithm{};
    int PVAlgorithm{-1};
    double Tolerance{1.0};
    double PX{0.0};
    double PY{0.0};
    double PZ{0.0};
    double PE{1.0};
    double PS{1.0};

  private:
    int computeBottleneck(const ttk::DiagramType &d1,
                          const ttk::DiagramType &d2,
                          std::vector<MatchingType> &matchings);

    double computeGeometricalRange(const ttk::DiagramType &CTDiagram1,
                                   const ttk::DiagramType &CTDiagram2) const;

    double computeMinimumRelevantPersistence(
      const ttk::DiagramType &CTDiagram1,
      const ttk::DiagramType &CTDiagram2) const;

    double distanceFunction(const PersistencePair &a,
                            const PersistencePair &b,
                            const int wasserstein) const;
    double diagonalDistanceFunction(const PersistencePair &a,
                                    const int wasserstein) const;

    void computeMinMaxSaddleNumberAndMapping(const ttk::DiagramType &CTDiagram,
                                             int &nbMin,
                                             int &nbMax,
                                             int &nbSaddle,
                                             std::vector<int> &minMap,
                                             std::vector<int> &maxMap,
                                             std::vector<int> &sadMap,
                                             const double zeroThresh) const;

    void buildCostMatrices(const ttk::DiagramType &CTDiagram1,
                           const ttk::DiagramType &CTDiagram2,
                           double zeroThresh,
                           std::vector<std::vector<double>> &minMatrix,
                           std::vector<std::vector<double>> &maxMatrix,
                           std::vector<std::vector<double>> &sadMatrix,
                           bool reverseMin,
                           bool reverseMax,
                           bool reverseSad,
                           int wasserstein) const;

    double buildMappings(const std::vector<MatchingType> &inputMatchings,
                         bool transposeGlobal,
                         bool transposeLocal,
                         std::vector<MatchingType> &outputMatchings,
                         const std::vector<int> &m1,
                         const std::vector<int> &m2,
                         int wasserstein) const;
  };

} // namespace ttk
