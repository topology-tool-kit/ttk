/// \ingroup base
/// \class ttk::BottleneckDistance
/// \author Maxime Soler <soler.maxime@total.com>
#pragma once

// base code includes
#include <Debug.h>
#include <PersistenceDiagramUtils.h>

namespace ttk {

  class BottleneckDistance : virtual public Debug {
  public:
    BottleneckDistance();

    int execute(const ttk::DiagramType &diag0,
                const ttk::DiagramType &diag1,
                std::vector<MatchingType> &matchings,
                const bool usePersistenceMetric);

    inline void setPersistencePercentThreshold(const double t) {
      zeroThreshold_ = t;
    }
    inline void setPX(const double px) {
      px_ = px;
    }
    inline void setPY(const double py) {
      py_ = py;
    }
    inline void setPZ(const double pz) {
      pz_ = pz;
    }
    inline void setPE(const double pe) {
      pe_ = pe;
    }
    inline void setPS(const double ps) {
      ps_ = ps;
    }
    inline void setAlgorithm(const std::string &algorithm) {
      algorithm_ = algorithm;
    }
    inline void setPVAlgorithm(const int algorithm) {
      pvAlgorithm_ = algorithm;
    }
    inline void setWasserstein(const std::string &wasserstein) {
      wasserstein_ = wasserstein;
    }

    double getDistance() {
      return distance_;
    }

  protected:
    double distance_{-1.0};

    std::string wasserstein_{"inf"};
    std::string algorithm_{};
    int pvAlgorithm_{-1};
    double zeroThreshold_{};
    double px_{};
    double py_{};
    double pz_{};
    double pe_{};
    double ps_{};

  private:
    int computeBottleneck(const ttk::DiagramType &d1,
                          const ttk::DiagramType &d2,
                          std::vector<MatchingType> &matchings,
                          bool usePersistenceMetric);

    double computeGeometricalRange(const ttk::DiagramType &CTDiagram1,
                                   const ttk::DiagramType &CTDiagram2,
                                   int d1Size,
                                   int d2Size) const;

    double computeMinimumRelevantPersistence(const ttk::DiagramType &CTDiagram1,
                                             const ttk::DiagramType &CTDiagram2,
                                             int d1Size,
                                             int d2Size) const;

    void computeMinMaxSaddleNumberAndMapping(const ttk::DiagramType &CTDiagram,
                                             int dSize,
                                             int &nbMin,
                                             int &nbMax,
                                             int &nbSaddle,
                                             std::vector<int> &minMap,
                                             std::vector<int> &maxMap,
                                             std::vector<int> &sadMap,
                                             const double zeroThresh) const;

    template <typename distFuncType, typename diagFuncType>
    void buildCostMatrices(const ttk::DiagramType &CTDiagram1,
                           const ttk::DiagramType &CTDiagram2,
                           int d1Size,
                           int d2Size,
                           const distFuncType &distanceFunction,
                           const diagFuncType &diagonalDistanceFunction,
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
