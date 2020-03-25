/// \ingroup base
/// \class ttk::PersistenceDiagramDistanceMatrix
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
/// \sa ttkPersistenceDiagramDistanceMatrix

#pragma once

#ifndef BNodeType
#define BNodeType ttk::CriticalType
#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1 ttk::CriticalType::Saddle1
#define BSaddle2 ttk::CriticalType::Saddle2
#define BIdVertex ttk::SimplexId
#endif

// base code includes
//
#include <Wrapper.h>
//
#include <PersistenceDiagram.h>
//
#include <limits>
//
#include <PDDistMat.h>
//

using namespace std;
using namespace ttk;

namespace ttk {
  class PersistenceDiagramDistanceMatrix : public Debug {

  public:
    PersistenceDiagramDistanceMatrix() {
      wasserstein_ = 2;
      use_progressive_ = 1;
      use_kmeanspp_ = 0;
      use_accelerated_ = 0;
      numberOfInputs_ = 0;
      threadNumber_ = 1;
    };

    ~PersistenceDiagramDistanceMatrix(){};

    void execute(std::vector<std::vector<DiagramTuple>> &intermediateDiagrams);

    inline void setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
    }

    inline void setWasserstein(const std::string &wasserstein) {
      wasserstein_ = (wasserstein == "inf") ? -1 : stoi(wasserstein);
    }

    inline void setThreadNumber(const int &ThreadNumber) {
      threadNumber_ = ThreadNumber;
    }

    inline void setUseProgressive(const bool use_progressive) {
      use_progressive_ = use_progressive;
    }

    inline void setAlpha(const double alpha) {
      alpha_ = alpha;
    }
    inline void setLambda(const double lambda) {
      lambda_ = lambda;
    }

    inline void setTimeLimit(const double time_limit) {
      time_limit_ = time_limit;
    }

    inline void setUseKmeansppInit(const bool UseKmeansppInit) {
      use_kmeanspp_ = UseKmeansppInit;
    }

    inline void setUseAccelerated(const bool UseAccelerated) {
      use_accelerated_ = UseAccelerated;
    }

    inline void setNumberOfClusters(const int NumberOfClusters) {
      n_clusters_ = NumberOfClusters;
    }
    inline void setForceUseOfAlgorithm(const bool forceUseOfAlgorithm) {
      forceUseOfAlgorithm_ = forceUseOfAlgorithm;
    }
    inline void setDeterministic(const bool deterministic) {
      deterministic_ = deterministic;
    }
    inline void setPairTypeClustering(const int pairTypeClustering) {
      pairTypeClustering_ = pairTypeClustering;
    }

    inline void setUseDeltaLim(const bool useDeltaLim) {
      useDeltaLim_ = useDeltaLim;
    }

    inline void setDistanceWritingOptions(const int distanceWritingOptions) {
      distanceWritingOptions_ = distanceWritingOptions;
    }

    inline void setDeltaLim(const double deltaLim) {
      deltaLim_ = deltaLim;
    }
    inline void setOutputDistanceMatrix(const bool arg) {
      outputDistanceMatrix_ = arg;
    }
    inline void setUseFullDiagrams(const bool arg) {
      useFullDiagrams_ = arg;
    }
    inline void setPerClusterDistanceMatrix(const bool arg) {
      perClusterDistanceMatrix_ = arg;
    }

    inline const std::vector<std::vector<double>> &&getDiagramsDistMat() {
      return std::move(diagramsDistMat_);
    }

  protected:
    // Critical pairs used for clustering
    // 0:min-saddles ; 1:saddles-saddles ; 2:sad-max ; else : all

    double deltaLim_;
    bool useDeltaLim_;
    int distanceWritingOptions_;
    int pairTypeClustering_;
    bool forceUseOfAlgorithm_;
    bool deterministic_;
    int wasserstein_;
    int n_clusters_;

    int numberOfInputs_;
    int threadNumber_;
    bool use_progressive_;
    bool use_accelerated_;
    bool use_kmeanspp_;
    double alpha_;
    double lambda_;
    double time_limit_;

    int points_added_;
    int points_deleted_;

    std::vector<BidderDiagram<double>> bidder_diagrams_;
    std::vector<GoodDiagram<double>> barycenter_goods_;

    bool outputDistanceMatrix_{false};
    bool useFullDiagrams_{false};
    bool perClusterDistanceMatrix_{false};
    std::vector<std::vector<double>> diagramsDistMat_{};
  };
} // namespace ttk
