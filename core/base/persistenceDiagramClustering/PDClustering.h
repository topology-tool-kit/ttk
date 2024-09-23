/// \ingroup base
/// \class ttk::PDClustering
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

#include <KDTree.h>
#include <PDBarycenter.h>
#include <PersistenceDiagramAuction.h>

#include <array>
#include <limits>

namespace ttk {

  class PDClustering : virtual public Debug {

  public:
    PDClustering() {
      threadNumber_ = 1;
      this->setDebugMsgPrefix("PersistenceDiagramClustering");
    }

    ~PDClustering() override = default;

    std::vector<int>
      execute(std::vector<DiagramType> &final_centroids,
              std::vector<std::vector<std::vector<std::vector<MatchingType>>>>
                &all_matchings);

    std::array<double, 3> getDistances() const {
      return {this->cost_min_, this->cost_sad_, this->cost_max_};
    }

    double getMostPersistent(int type = -1);
    std::vector<std::vector<int>> get_centroids_sizes();
    double getLessPersistent(int type = -1);
    std::vector<std::vector<double>> getMinDiagonalPrices();
    std::vector<std::vector<double>> getMinPrices();

    void correctMatchings(
      std::vector<std::vector<std::vector<std::vector<MatchingType>>>>
        &previous_matchings);

    double computeDistance(const BidderDiagram &D1,
                           const BidderDiagram &D2,
                           const double delta_lim);
    double computeDistance(const BidderDiagram &D1,
                           const GoodDiagram &D2,
                           const double delta_lim);
    double computeDistance(BidderDiagram *const D1,
                           const GoodDiagram *const D2,
                           const double delta_lim);
    double computeDistance(const GoodDiagram &D1,
                           const GoodDiagram &D2,
                           const double delta_lim);

    GoodDiagram centroidWithZeroPrices(const GoodDiagram &centroid);
    BidderDiagram centroidToDiagram(const GoodDiagram &centroid);
    GoodDiagram diagramToCentroid(const BidderDiagram &diagram);
    BidderDiagram diagramWithZeroPrices(const BidderDiagram &diagram);

    void setBidderDiagrams();
    void initializeEmptyClusters();
    void initializeCentroids();
    void initializeCentroidsKMeanspp();
    void initializeAcceleratedKMeans();
    void initializeBarycenterComputers(std::vector<double> &min_persistence);
    void printDistancesToFile();
    void printMatchings(std::vector<std::vector<std::vector<MatchingType>>>);
    void printRealDistancesToFile();
    void printPricesToFile(int);
    double computeRealCost();

    std::vector<double> enrichCurrentBidderDiagrams(
      std::vector<double> &previous_min_persistence,
      std::vector<double> &min_persistence,
      std::vector<std::vector<double>> &initial_diagonal_prices,
      std::vector<std::vector<double>> &initial_off_diagonal_prices,
      std::vector<int> &min_points_to_add,
      bool add_points_to_barycenter,
      bool first_enrichment);

    std::vector<std::vector<double>> getDistanceMatrix();
    void getCentroidDistanceMatrix();
    void computeDistanceToCentroid();

    void updateClusters();
    void invertClusters();
    void invertInverseClusters();
    void computeBarycenterForTwoGlobal(
      std::vector<std::vector<std::vector<std::vector<MatchingType>>>> &);

    void computeBarycenterForTwo(
      std::vector<std::vector<MatchingType>> &matchings,
      std::vector<std::vector<int>> &bidders_ids,
      std::vector<BidderDiagram> &current_bidder_diagrams,
      std::vector<BidderDiagram> &bidder_diagrams,
      GoodDiagram &barycenter);

    void acceleratedUpdateClusters();
    std::vector<double> updateCentroidsPosition(
      std::vector<std::vector<double>> *min_price,
      std::vector<std::vector<double>> *min_diag_price,
      std::vector<std::vector<std::vector<std::vector<MatchingType>>>>
        &all_matchings,
      int only_matchings);

    inline void resetDosToOriginalValues() {
      do_min_ = original_dos[0];
      do_sad_ = original_dos[1];
      do_max_ = original_dos[2];
    }
    inline int setDiagrams(std::vector<DiagramType> *data_min,
                           std::vector<DiagramType> *data_saddle,
                           std::vector<DiagramType> *data_max) {
      inputDiagramsMin_ = data_min;
      inputDiagramsSaddle_ = data_saddle;
      inputDiagramsMax_ = data_max;
      return 0;
    }

    inline int setDos(bool doMin, bool doSad, bool doMax) {
      do_min_ = doMin;
      do_sad_ = doSad;
      do_max_ = doMax;

      original_dos[0] = do_min_;
      original_dos[1] = do_sad_;
      original_dos[2] = do_max_;
      return 0;
    }

    inline int setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
      return 0;
    }

    inline int setK(const int k) {
      k_ = k;
      return 0;
    }

    inline void setWasserstein(const int &wasserstein) {
      wasserstein_ = wasserstein;
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

    inline void setNonMatchingWeight(double nonMatchingWeight) {
      nonMatchingWeight_ = nonMatchingWeight;
    }

    inline void printClustering() {
      std::string msg{};
      for(int c = 0; c < k_; ++c) {
        msg.append(" Cluster " + std::to_string(c) + " = {");
        for(size_t idx = 0; idx < clustering_[c].size(); ++idx) {
          if(idx == clustering_[c].size() - 1) {
            msg.append(std::to_string(clustering_[c][idx]) + "}");
            this->printMsg(msg);
            msg = "";
            // msg << clustering_[c][idx] << "}" << std::endl;
          } else {
            msg.append(std::to_string(clustering_[c][idx]) + ", ");
            // msg << clustering_[c][idx] << ", ";
          }
        }
      }
      // cout<<msg.str()<<endl;
    }

    inline void printOldClustering() {
      std::stringstream msg;
      for(int c = 0; c < k_; ++c) {
        msg << "Cluster " << c << " = {";
        for(size_t idx = 0; idx < old_clustering_[c].size(); ++idx) {
          if(idx == old_clustering_[c].size() - 1) {
            msg << old_clustering_[c][idx] << "}" << std::endl;
          } else {
            msg << old_clustering_[c][idx] << ", ";
          }
        }
      }
      this->printMsg(msg.str());
    }

  protected:
    std::vector<PDBarycenter> barycenter_computer_min_{};
    std::vector<PDBarycenter> barycenter_computer_sad_{};
    std::vector<PDBarycenter> barycenter_computer_max_{};

    bool barycenter_inputs_reset_flag;
    bool precision_criterion_{false};
    bool precision_max_{false};
    bool precision_min_{false};
    bool precision_sad_{false};
    bool forceUseOfAlgorithm_{false};
    bool deterministic_{true};
    int wasserstein_{2};
    double geometrical_factor_{1.0};
    double deltaLim_;
    bool UseDeltaLim_{false};
    int distanceWritingOptions_{0};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double lambda_;
    double nonMatchingWeight_ = 1.0;

    int k_;
    int numberOfInputs_;
    bool use_progressive_{true};
    bool use_accelerated_;
    bool use_kmeanspp_;
    bool use_kdtree_;
    double time_limit_{std::numeric_limits<double>::max()};

    double epsilon_min_{1e-8};
    std::array<double, 3> epsilon_;
    double cost_;
    double cost_min_{0.0};
    double cost_sad_{0.0};
    double cost_max_{0.0};

    std::vector<std::vector<int>> current_bidder_ids_min_;
    std::vector<std::vector<int>> current_bidder_ids_sad_;
    std::vector<std::vector<int>> current_bidder_ids_max_;
    std::vector<DiagramType> *inputDiagramsMin_;
    std::vector<DiagramType> *inputDiagramsSaddle_;
    std::vector<DiagramType> *inputDiagramsMax_;

    std::array<bool, 3> original_dos;

    bool do_min_;
    std::vector<BidderDiagram> bidder_diagrams_min_;
    std::vector<BidderDiagram> current_bidder_diagrams_min_;
    std::vector<GoodDiagram> centroids_min_;
    std::vector<GoodDiagram> centroids_with_price_min_;

    bool do_sad_;
    std::vector<BidderDiagram> bidder_diagrams_saddle_;
    std::vector<BidderDiagram> current_bidder_diagrams_saddle_;
    std::vector<GoodDiagram> centroids_saddle_;
    std::vector<GoodDiagram> centroids_with_price_saddle_;

    bool do_max_;
    std::vector<BidderDiagram> bidder_diagrams_max_;
    std::vector<BidderDiagram> current_bidder_diagrams_max_;
    std::vector<GoodDiagram> centroids_max_;
    std::vector<GoodDiagram> centroids_with_price_max_;

    std::vector<std::vector<int>> clustering_;
    std::vector<std::vector<int>> old_clustering_;
    std::vector<int> inv_clustering_;

    std::vector<std::vector<int>> centroids_sizes_;

    std::vector<bool> r_;
    std::vector<double> u_;
    std::vector<std::vector<double>> l_;
    std::vector<std::vector<double>> centroidsDistanceMatrix_{};
    std::vector<double> distanceToCentroid_{};

    int n_iterations_;
  };
} // namespace ttk
