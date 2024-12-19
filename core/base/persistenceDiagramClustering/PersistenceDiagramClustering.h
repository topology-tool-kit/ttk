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

    void setTimeLimit(double timeLimit) {
      this->TimeLimit = timeLimit;
    }

    void setForceUseOfAlgorithm(bool forceUseOfAlgorithm) {
      this->ForceUseOfAlgorithm = forceUseOfAlgorithm;
    }

    void setDeltaLim(double DeltaLimNew) {
      this->DeltaLim = DeltaLimNew;
    }

    void setUseAdditionalPrecision(bool Precision) {
      this->UseAdditionalPrecision = Precision;
    }

    void setUseProgressive(bool UseProgressive_) {
      this->UseProgressive = UseProgressive_;
    }

    void setUseInterruptible(bool UseInterruptible_) {
      this->UseInterruptible = UseInterruptible_;
    }

    void setDeterministic(bool Deterministic_) {
      this->Deterministic = Deterministic_;
    }

    void setAlpha(double Alpha_) {
      this->Alpha = Alpha_;
    }

    void setUseAccelerated(bool UseAccelerated_) {
      this->UseAccelerated = UseAccelerated_;
    }

    void setUseKmeansppInit(bool UseKmeansppInit_) {
      this->UseKmeansppInit = UseKmeansppInit_;
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
  };

} // namespace ttk
