/// \ingroup base
/// \class ttk::PDDistMat
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramClustering

#pragma once

#include <Auction.h>
#include <KDTree.h>
//
#include <array>
#include <limits>
//

using namespace std;

namespace ttk {

  using DiagramTuple = std::tuple<ttk::SimplexId,
                                  ttk::CriticalType,
                                  ttk::SimplexId,
                                  ttk::CriticalType,
                                  double,
                                  ttk::SimplexId,
                                  double,
                                  float,
                                  float,
                                  float,
                                  double,
                                  float,
                                  float,
                                  float>;
  using MatchingTuple = std::tuple<ttk::SimplexId, ttk::SimplexId, double>;

  class PDDistMat : public Debug {

  public:
    std::vector<int> execute();

    double getMostPersistent(int type = -1);
    vector<vector<int>> get_centroids_sizes();
    double getLessPersistent(int type = -1);
    std::vector<std::vector<double>> getMinDiagonalPrices();
    std::vector<std::vector<double>> getMinPrices();

    double computeDistance(const BidderDiagram<double> &D1,
                           const BidderDiagram<double> &D2,
                           const double delta_lim);
    double computeDistance(const BidderDiagram<double> D1,
                           const GoodDiagram<double> D2,
                           const double delta_lim);
    double computeDistance(BidderDiagram<double> *const D1,
                           const GoodDiagram<double> *const D2,
                           const double delta_lim);
    double computeDistance(const GoodDiagram<double> &D1,
                           const GoodDiagram<double> &D2,
                           const double delta_lim);

    GoodDiagram<double>
      centroidWithZeroPrices(const GoodDiagram<double> centroid);
    BidderDiagram<double> centroidToDiagram(const GoodDiagram<double> centroid);
    GoodDiagram<double> diagramToCentroid(const BidderDiagram<double> diagram);
    BidderDiagram<double>
      diagramWithZeroPrices(const BidderDiagram<double> diagram);

    void setBidderDiagrams();
    void initializeEmptyClusters();
    void initializeCentroids();
    void initializeCentroidsKMeanspp();
    void initializeAcceleratedKMeans();
    void printDistancesToFile();
    double computeRealCost();

    std::vector<double> enrichCurrentBidderDiagrams(
      std::vector<double> previous_min_persistence,
      std::vector<double> min_persistence,
      std::vector<std::vector<double>> initial_diagonal_prices,
      std::vector<std::vector<double>> initial_off_diagonal_points,
      std::vector<int> min_points_to_add,
      bool add_points_to_barycenter);

    std::vector<std::vector<double>> getDistanceMatrix();
    void getCentroidDistanceMatrix();
    void computeDiagramsDistanceMatrix();

    void updateClusters();
    void invertClusters();
    void invertInverseClusters();

    void acceleratedUpdateClusters();

    inline void resetDosToOriginalValues() {
      do_min_ = original_dos[0];
      do_sad_ = original_dos[1];
      do_max_ = original_dos[2];
    }
    inline void setDiagrams(std::vector<std::vector<DiagramTuple>> *data_min,
                            std::vector<std::vector<DiagramTuple>> *data_saddle,
                            std::vector<std::vector<DiagramTuple>> *data_max) {
      inputDiagramsMin_ = data_min;
      inputDiagramsSaddle_ = data_saddle;
      inputDiagramsMax_ = data_max;
    }

    inline void setDos(bool doMin, bool doSad, bool doMax) {
      do_min_ = doMin;
      do_sad_ = doSad;
      do_max_ = doMax;

      original_dos[0] = do_min_;
      original_dos[1] = do_sad_;
      original_dos[2] = do_max_;
    }

    inline void setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
    }

    inline void setK(const int k) {
      k_ = k;
    }

    inline void setWasserstein(const int &wasserstein) {
      wasserstein_ = wasserstein;
    }

    inline void setThreadNumber(const int &threadNumber) {
      threadNumber_ = threadNumber;
    }

    inline void setUseProgressive(const bool use_progressive) {
      use_progressive_ = use_progressive;
    }

    inline void setKMeanspp(const bool use_kmeanspp) {
      use_kmeanspp_ = use_kmeanspp;
    }

    inline void setUseKDTree(const bool use_kdtree) {
      use_kdtree_ = use_kdtree;
    }

    inline void setAccelerated(const bool use_accelerated) {
      use_accelerated_ = use_accelerated;
    }

    inline void setTimeLimit(const double time_limit) {
      time_limit_ = time_limit;
    }

    inline void setGeometricalFactor(const double geometrical_factor) {
      geometrical_factor_ = geometrical_factor;
    }
    inline void setLambda(const double lambda) {
      lambda_ = lambda;
    }
    inline void setForceUseOfAlgorithm(const bool forceUseOfAlgorithm) {
      forceUseOfAlgorithm_ = forceUseOfAlgorithm;
    }
    inline void setDeterministic(const bool deterministic) {
      deterministic_ = deterministic;
    }

    inline void setUseDeltaLim(const bool UseDeltaLim) {
      UseDeltaLim_ = UseDeltaLim;
      if(UseDeltaLim_) {
        epsilon_min_ = 1e-8;
      } else {
        epsilon_min_ = 5e-5;
      }
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

    inline void printClustering() {
      std::stringstream msg;
      for(int c = 0; c < k_; ++c) {
        msg << "[PersistenceDiagramClustering] Cluster " << c << " = {";
        for(unsigned int idx = 0; idx < clustering_[c].size(); ++idx) {
          if(idx == clustering_[c].size() - 1) {
            msg << clustering_[c][idx] << "}" << std::endl;
          } else {
            msg << clustering_[c][idx] << ", ";
          }
        }
      }
      dMsg(std::cout, msg.str(), infoMsg);
    }

    inline void printOldClustering() {
      std::stringstream msg;
      for(int c = 0; c < k_; ++c) {
        msg << "Cluster " << c << " = {";
        for(unsigned int idx = 0; idx < old_clustering_[c].size(); ++idx) {
          if(idx == old_clustering_[c].size() - 1) {
            msg << old_clustering_[c][idx] << "}" << std::endl;
          } else {
            msg << old_clustering_[c][idx] << ", ";
          }
        }
      }
      dMsg(std::cout, msg.str(), infoMsg);
    }

    inline const std::vector<std::vector<double>> &&
      getDiagramsDistanceMatrix() {
      return std::move(diagramsDistanceMatrix_);
    }

    template <typename type>
    static type abs(const type var) {
      return (var >= 0) ? var : -var;
    }

  protected:
    bool barycenter_inputs_reset_flag;
    bool precision_criterion_{false};
    bool precision_max_{false};
    bool precision_min_{false};
    bool precision_sad_{false};
    bool forceUseOfAlgorithm_{false};
    bool deterministic_{true};
    int wasserstein_{2};
    double geometrical_factor_{1};
    double deltaLim_;
    bool UseDeltaLim_{false};
    int distanceWritingOptions_{0};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double lambda_;

    int k_;
    int numberOfInputs_;
    int threadNumber_{1};
    bool use_progressive_{true};
    bool use_accelerated_;
    bool use_kmeanspp_;
    bool use_kdtree_;
    double time_limit_{std::numeric_limits<double>::max()};

    double epsilon_min_{1e-8};
    std::array<double, 3> epsilon_;
    double cost_;
    double cost_min_{0};
    double cost_sad_{0};
    double cost_max_{0};

    std::vector<std::vector<int>> current_bidder_ids_min_;
    std::vector<std::vector<int>> current_bidder_ids_sad_;
    std::vector<std::vector<int>> current_bidder_ids_max_;
    std::vector<std::vector<DiagramTuple>> *inputDiagramsMin_;
    std::vector<std::vector<DiagramTuple>> *inputDiagramsSaddle_;
    std::vector<std::vector<DiagramTuple>> *inputDiagramsMax_;

    std::array<bool, 3> original_dos;

    bool do_min_;
    std::vector<BidderDiagram<double>> bidder_diagrams_min_;
    std::vector<BidderDiagram<double>> current_bidder_diagrams_min_;
    std::vector<GoodDiagram<double>> centroids_min_;
    std::vector<GoodDiagram<double>> centroids_with_price_min_;

    bool do_sad_;
    std::vector<BidderDiagram<double>> bidder_diagrams_saddle_;
    std::vector<BidderDiagram<double>> current_bidder_diagrams_saddle_;
    std::vector<GoodDiagram<double>> centroids_saddle_;
    std::vector<GoodDiagram<double>> centroids_with_price_saddle_;

    bool do_max_;
    std::vector<BidderDiagram<double>> bidder_diagrams_max_;
    std::vector<BidderDiagram<double>> current_bidder_diagrams_max_;
    std::vector<GoodDiagram<double>> centroids_max_;
    std::vector<GoodDiagram<double>> centroids_with_price_max_;

    std::vector<std::vector<int>> clustering_;
    std::vector<std::vector<int>> old_clustering_;
    std::vector<int> inv_clustering_;

    std::vector<std::vector<int>> centroids_sizes_;

    std::vector<bool> r_;
    std::vector<double> u_;
    std::vector<std::vector<double>> l_;
    std::vector<std::vector<double>> centroidsDistanceMatrix_{};
    std::vector<std::vector<double>> diagramsDistanceMatrix_{};
    bool outputDistanceMatrix_{false};
    bool useFullDiagrams_{false};

    int n_iterations_;
  };
} // namespace ttk
