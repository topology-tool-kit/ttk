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
#include <PersistenceDiagramUtils.h>

#include <limits>

namespace ttk {

  class PDBarycenter : public Debug {

  public:
    PDBarycenter() {
      threadNumber_ = 1;
      this->setDebugMsgPrefix("PersistenceDiagramBarycenter");
    }

    ~PDBarycenter() override = default;

    std::vector<std::vector<MatchingType>> execute(DiagramType &barycenter);
    std::vector<std::vector<MatchingType>>
      executeMunkresBarycenter(DiagramType &barycenter);
    std::vector<std::vector<MatchingType>>
      executeAuctionBarycenter(DiagramType &barycenter);
    std::vector<std::vector<MatchingType>>
      executePartialBiddingBarycenter(DiagramType &barycenter);

    void setBidderDiagrams();
    double enrichCurrentBidderDiagrams(
      double previous_min_persistence,
      double min_persistence,
      std::vector<double> &initial_diagonal_prices,
      std::vector<double> &initial_off_diagonal_prices,
      int min_points_to_add,
      bool add_points_to_barycenter = true);
    void setInitialBarycenter(double min_persistence);
    double getMaxPersistence();
    double getLowestPersistence();
    double getMinimalPrice(int i);

    using KDT = PersistenceDiagramAuction::KDT;
    using KDTreePair
      = std::pair<typename KDT::KDTreeRoot, typename KDT::KDTreeMap>;
    KDTreePair getKDTree() const;

    void runMatching(double *total_cost,
                     double epsilon,
                     std::vector<int> &sizes,
                     KDT &kdt,
                     std::vector<KDT *> &correspondence_kdt_map,
                     std::vector<double> *min_diag_price,
                     std::vector<double> *min_price,
                     std::vector<std::vector<MatchingType>> *all_matchings,
                     bool use_kdt,
                     bool actual_distance);

    void
      runMatchingAuction(double *total_cost,
                         std::vector<int> &sizes,
                         KDT &kdt,
                         std::vector<KDT *> &correspondence_kdt_map,
                         std::vector<double> *min_diag_price,
                         std::vector<std::vector<MatchingType>> *all_matchings,
                         bool use_kdt,
                         bool actual_distance);

    double updateBarycenter(std::vector<std::vector<MatchingType>> &matchings);

    double computeRealCost();
    bool isPrecisionObjectiveMet(double, int);
    bool hasBarycenterConverged(
      std::vector<std::vector<MatchingType>> &matchings,
      std::vector<std::vector<MatchingType>> &previous_matchings);
    std::vector<std::vector<MatchingType>> correctMatchings(
      std::vector<std::vector<MatchingType>> &previous_matchings);

    bool is_matching_stable();

    double getEpsilon(double rho);
    double getRho(double epsilon);

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

    inline int setDiagrams(std::vector<DiagramType> *data) {
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

    inline void setGeometricalFactor(const double geometrical_factor) {
      geometrical_factor_ = geometrical_factor;
    }

    inline void setLambda(const double lambda) {
      lambda_ = lambda;
    }

    inline void setCurrentBidders(std::vector<BidderDiagram> &diagrams) {
      current_bidder_diagrams_ = diagrams;
    }

    inline void setCurrentBarycenter(std::vector<GoodDiagram> &barycenters) {
      barycenter_goods_ = barycenters;
    }

    inline std::vector<BidderDiagram> &getCurrentBidders() {
      return current_bidder_diagrams_;
    }

    inline std::vector<GoodDiagram> &getCurrentBarycenter() {
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
        nt1_ = ttk::CriticalType::Local_minimum;
        nt2_ = ttk::CriticalType::Saddle1;
      } else if(diagramType_ == 1) {
        nt1_ = ttk::CriticalType::Saddle1;
        nt2_ = ttk::CriticalType::Saddle2;
      } else {
        nt1_ = ttk::CriticalType::Saddle2;
        nt2_ = ttk::CriticalType::Local_maximum;
      }
    }

    double getCost() {
      return cost_;
    }

  protected:
    // std::vector<bool> precision_objective_;
    std::vector<double> precision_;

    // to kill any randomness
    bool deterministic_{false};

    std::string method_{"Partial Bidding"};
    int wasserstein_{2};

    double geometrical_factor_{1.0};

    // lambda_ : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda_ = 1 : extremum (min if pair min-sad, max
    // if pair sad-max) lambda_ = 0 : saddle (awful stability) lambda_ = 1/2 :
    // middle of the 2 critical points of the pair (bad stability)
    double lambda_;

    int diagramType_;
    ttk::CriticalType nt1_;
    ttk::CriticalType nt2_;
    double cost_;
    int numberOfInputs_;
    bool use_progressive_{true};
    double epsilon_min_{1e-5};
    std::vector<DiagramType> *inputDiagrams_;

    int points_added_;
    int points_deleted_;

    std::vector<std::vector<double>> all_matchings_;
    std::vector<std::vector<double>> all_old_matchings_;
    std::vector<BidderDiagram> bidder_diagrams_;
    std::vector<BidderDiagram> current_bidder_diagrams_;
    std::vector<std::vector<int>> current_bidder_ids_;
    std::vector<GoodDiagram> barycenter_goods_;

    bool reinit_prices_{true};
    bool epsilon_decreases_{true};
    bool early_stoppage_{true};
  };
} // namespace ttk
