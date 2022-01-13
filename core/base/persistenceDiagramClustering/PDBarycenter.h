/// \ingroup base
/// \class ttk::PDBarycenter
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
#include <PersistenceDiagramAuction.h>

#include <limits>

namespace ttk {
  template <typename dataType>
  class PDBarycenter : public Debug {

  public:
    PDBarycenter() {
      wasserstein_ = 2;
      geometrical_factor_ = 1;
      method_ = "Partial Bidding";
      threadNumber_ = 1;
      use_progressive_ = true;
      time_limit_ = std::numeric_limits<double>::max();
      epsilon_min_ = 1e-5;
      reinit_prices_ = true;
      deterministic_ = false;
      epsilon_decreases_ = true;
      early_stoppage_ = true;
      this->setDebugMsgPrefix("PersistenceDiagramBarycenter");
    }

    ~PDBarycenter() = default;

    std::vector<std::vector<matchingTuple>>
      execute(std::vector<diagramTuple> &barycenter);
    std::vector<std::vector<matchingTuple>>
      executeMunkresBarycenter(std::vector<diagramTuple> &barycenter);
    std::vector<std::vector<matchingTuple>>
      executeAuctionBarycenter(std::vector<diagramTuple> &barycenter);
    std::vector<std::vector<matchingTuple>>
      executePartialBiddingBarycenter(std::vector<diagramTuple> &barycenter);

    void setBidderDiagrams();
    dataType
      enrichCurrentBidderDiagrams(dataType previous_min_persistence,
                                  dataType min_persistence,
                                  std::vector<dataType> initial_diagonal_prices,
                                  std::vector<dataType> initial_prices,
                                  int min_points_to_add,
                                  bool add_points_to_barycenter = true);
    void setInitialBarycenter(dataType min_persistence);
    dataType getMaxPersistence();
    dataType getLowestPersistence();
    dataType getMinimalPrice(int i);
    using KDTreePair = std::pair<typename KDTree<dataType>::KDTreeRoot,
                                 typename KDTree<dataType>::KDTreeMap>;
    KDTreePair getKDTree() const;

    void runMatching(dataType *total_cost,
                     dataType epsilon,
                     std::vector<int> sizes,
                     KDTree<dataType> &kdt,
                     std::vector<KDTree<dataType> *> &correspondance_kdt_map,
                     std::vector<dataType> *min_diag_price,
                     std::vector<dataType> *min_price,
                     std::vector<std::vector<matchingTuple>> *all_matchings,
                     bool use_kdt,
                     int compute_only_distance);

    void runMatchingAuction(
      dataType *total_cost,
      std::vector<int> sizes,
      KDTree<dataType> &kdt,
      std::vector<KDTree<dataType> *> &correspondance_kdt_map,
      std::vector<dataType> *min_diag_price,
      std::vector<std::vector<matchingTuple>> *all_matchings,
      bool use_kdt);

    dataType
      updateBarycenter(std::vector<std::vector<matchingTuple>> &matchings);

    dataType computeRealCost();
    bool isPrecisionObjectiveMet(dataType, int);
    bool hasBarycenterConverged(
      std::vector<std::vector<matchingTuple>> &matchings,
      std::vector<std::vector<matchingTuple>> &previous_matchings);
    std::vector<std::vector<matchingTuple>> correctMatchings(
      std::vector<std::vector<matchingTuple>> previous_matchings);

    bool is_matching_stable();

    dataType getEpsilon(dataType rho);
    dataType getRho(dataType epsilon);

    // 		inline int setDiagram(int idx, void* data){
    // 			if(idx < numberOfInputs_){
    // 			inputData_[idx] = data;
    // 			}
    // 			else{
    // 			return -1;
    // 			}
    // 			return 0;
    // 		}

    inline void setDeterministic(const bool deterministic) {
      deterministic_ = deterministic;
    }

    inline void setMethod(const int &method) {
      if(method == 0) {
        method_ = "Partial Bidding";
      }
      if(method == 1) {
        method_ = "Munkres";
      } else if(method == 2) {
        method_ = "Auction";
      }
    }

    inline int setDiagrams(std::vector<std::vector<diagramTuple>> *data) {
      inputDiagrams_ = data;
      return 0;
    }

    inline int setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
      // precision_objective_.resize(numberOfInputs_);
      precision_.resize(numberOfInputs_);
      return 0;
    }

    inline void setWasserstein(const int &wasserstein) {
      wasserstein_ = wasserstein;
    }

    inline void setUseProgressive(const bool use_progressive) {
      use_progressive_ = use_progressive;
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

    inline void
      setCurrentBidders(std::vector<BidderDiagram<dataType>> &diagrams) {
      current_bidder_diagrams_ = diagrams;
    }

    inline void
      setCurrentBarycenter(std::vector<GoodDiagram<dataType>> &barycenters) {
      barycenter_goods_ = barycenters;
    }

    inline std::vector<BidderDiagram<dataType>> &getCurrentBidders() {
      return current_bidder_diagrams_;
    }

    inline std::vector<GoodDiagram<dataType>> &getCurrentBarycenter() {
      return barycenter_goods_;
    }

    inline void setReinitPrices(const bool reinit_prices) {
      reinit_prices_ = reinit_prices;
    }

    inline void setEpsilonDecreases(const bool epsilon_decreases) {
      epsilon_decreases_ = epsilon_decreases;
    }

    inline void setEarlyStoppage(const bool early_stoppage) {
      early_stoppage_ = early_stoppage;
    }

    inline void setDiagramType(const int &diagramType) {
      diagramType_ = diagramType;
      if(diagramType_ == 0) {
        nt1_ = BLocalMin;
        nt2_ = BSaddle1;
      } else if(diagramType_ == 1) {
        nt1_ = BSaddle1;
        nt2_ = BSaddle2;
      } else {
        nt1_ = BSaddle2;
        nt2_ = BLocalMax;
      }
    }

    dataType getCost() {
      return cost_;
    }

    template <typename type>
    static type abs(const type var) {
      return (var >= 0) ? var : -var;
    }

  protected:
    // std::vector<bool> precision_objective_;
    std::vector<dataType> precision_;

    // to kill any randomness
    bool deterministic_;

    std::string method_;
    int wasserstein_;

    double geometrical_factor_;

    // lambda_ : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda_ = 1 : extremum (min if pair min-sad, max
    // if pair sad-max) lambda_ = 0 : saddle (awful stability) lambda_ = 1/2 :
    // middle of the 2 critical points of the pair (bad stability)
    double lambda_;

    int diagramType_;
    BNodeType nt1_;
    BNodeType nt2_;
    dataType cost_;
    int numberOfInputs_;
    bool use_progressive_;
    double time_limit_;
    double epsilon_min_;
    std::vector<std::vector<diagramTuple>> *inputDiagrams_;

    int points_added_;
    int points_deleted_;

    std::vector<std::vector<dataType>> all_matchings_;
    std::vector<std::vector<dataType>> all_old_matchings_;
    std::vector<BidderDiagram<dataType>> bidder_diagrams_;
    std::vector<BidderDiagram<dataType>> current_bidder_diagrams_;
    std::vector<std::vector<int>> current_bidder_ids_;
    std::vector<GoodDiagram<dataType>> barycenter_goods_;

    bool reinit_prices_;
    bool epsilon_decreases_;
    bool early_stoppage_;
  };
} // namespace ttk

#include <PDBarycenterImpl.h>
