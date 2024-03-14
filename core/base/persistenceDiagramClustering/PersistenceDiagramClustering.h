/// \ingroup base
/// \class ttk::PersistenceDiagramClustering
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \brief TTK processing package for the computation of Wasserstein barycenters
/// and K-Means clusterings of a set of persistence diagrams.
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa ttkPersistenceDiagramClustering
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramClustering/">Persistence
///   Diagram Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramDistance/">Persistence
///   Diagram Distance example</a> \n
///

#pragma once

// base code includes
#include <PDClustering.h>
#include <PersistenceDiagramBarycenter.h>

namespace ttk {

  class PersistenceDiagramClustering : virtual public Debug {

  public:
    PersistenceDiagramClustering() {
      this->setDebugMsgPrefix("PersistenceDiagramClustering");
    }

    ~PersistenceDiagramClustering() override = default;

    std::vector<int> execute(
      std::vector<DiagramType> &intermediateDiagrams,
      std::vector<DiagramType> &centroids,
      std::vector<std::vector<std::vector<MatchingType>>> &all_matchings);

    std::array<double, 3> getDistances() const {
      return this->distances;
    }

    inline void setTimeLimit(const double tl) {
      TimeLimit = tl;
    }
    inline void setForceUseOfAlgorithm(bool data) {
      ForceUseOfAlgorithm = data;
    }
    inline void setUseInterruptible(bool data) {
      UseInterruptible = data;
    }
    inline void setUseProgressive(bool data) {
      UseProgressive = data;
    }
    inline void setUseCustomWeights(bool data) {
      UseCustomWeights = data;
    }
    inline void setCustomWeights(std::vector<double> *pdata) {
      CustomWeights = pdata;
    }

  protected:
    // Critical pairs used for clustering
    // 0:min-saddles ; 1:saddles-saddles ; 2:sad-max ; else : all

    // distance results per pair type
    std::array<double, 3> distances{};

    int DistanceWritingOptions{0};
    int PairTypeClustering{-1};
    bool Deterministic{true};
    int WassersteinMetric{2};

    bool UseProgressive{true};

    bool ForceUseOfAlgorithm{false};
    bool UseInterruptible{true};
    double Alpha{1.0};
    bool UseAdditionalPrecision{false};
    double DeltaLim{0.01};
    double Lambda{1.0};
    double TimeLimit{999999};
    double NonMatchingWeight = 1.0;

    int NumberOfClusters{1};
    bool UseAccelerated{false};
    bool UseKmeansppInit{false};

    int points_added_;
    int points_deleted_;

    std::vector<double> *CustomWeights{};
    bool UseCustomWeights{false};
  };

  void
    computeWeightedBarycenter(std::vector<DiagramType> &intermediateDiagrams,
                              std::vector<double> &weights,
                              DiagramType &barycenter,
                              std::vector<std::vector<MatchingType>> &matchings,
                              const ttk::Debug &dbg,
                              const bool ProgBarycenter);

} // namespace ttk
