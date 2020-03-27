/// \ingroup base
/// \class ttk::PersistenceDiagramDistanceMatrix
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramClustering

#pragma once

#include <array>
#include <limits>

#include <Auction.h>
#include <Wrapper.h>

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

  class PersistenceDiagramDistanceMatrix : public Debug {

  public:
    std::vector<std::vector<double>>
      execute(std::vector<std::vector<DiagramTuple>> &intermediateDiagrams);

    inline void setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
    }
    inline void setWasserstein(const std::string &wasserstein) {
      wasserstein_ = (wasserstein == "inf") ? -1 : stoi(wasserstein);
    }
    inline void setUseKDTree(const bool use_kdtree) {
      use_kdtree_ = use_kdtree;
    }
    inline void setPairTypeClustering(const int pairTypeClustering) {
      pairTypeClustering_ = pairTypeClustering;
    }
    inline void setAlpha(const double alpha) {
      geometrical_factor_ = alpha;
    }
    inline void setLambda(const double lambda) {
      lambda_ = lambda;
    }
    inline void setDeltaLim(const double deltaLim) {
      deltaLim_ = deltaLim;
    }
    inline void setUseFullDiagrams(const bool arg) {
      useFullDiagrams_ = arg;
    }

  protected:
    double getMostPersistent(
      const int type,
      const std::vector<BidderDiagram<double>> &bidder_diags_min,
      const std::vector<BidderDiagram<double>> &bidder_diags_sad,
      const std::vector<BidderDiagram<double>> &bidder_diags_max) const;
    double getLessPersistent(
      const int type,
      const std::vector<BidderDiagram<double>> &bidder_diags_min,
      const std::vector<BidderDiagram<double>> &bidder_diags_sad,
      const std::vector<BidderDiagram<double>> &bidder_diags_max) const;
    double computeDistance(const BidderDiagram<double> &D1,
                           const BidderDiagram<double> &D2,
                           const double delta_lim) const;
    std::vector<std::vector<double>> getDiagramsDistMat(
      const size_t nInputs,
      const bool useFullDiagrams,
      const std::vector<BidderDiagram<double>> &bidder_diags_min,
      const std::vector<BidderDiagram<double>> &bidder_diags_sad,
      const std::vector<BidderDiagram<double>> &bidder_diags_max,
      const std::vector<BidderDiagram<double>> &current_bidder_diags_min,
      const std::vector<BidderDiagram<double>> &current_bidder_diags_sad,
      const std::vector<BidderDiagram<double>> &current_bidder_diags_max) const;
    void setBidderDiagrams(
      std::vector<std::vector<DiagramTuple>> &inputDiagramsMin,
      std::vector<std::vector<DiagramTuple>> &inputDiagramsSad,
      std::vector<std::vector<DiagramTuple>> &inputDiagramsMax,
      std::vector<BidderDiagram<double>> &bidder_diags_min,
      std::vector<BidderDiagram<double>> &bidder_diags_sad,
      std::vector<BidderDiagram<double>> &bidder_diags_max,
      std::vector<BidderDiagram<double>> &current_bidder_diags_min,
      std::vector<BidderDiagram<double>> &current_bidder_diags_sad,
      std::vector<BidderDiagram<double>> &current_bidder_diags_max) const;
    std::array<double, 3> enrichCurrentBidderDiagrams(
      const std::array<double, 3> &previous_min_persistence,
      const std::array<double, 3> &min_persistence,
      const std::array<std::vector<double>, 3> initial_diagonal_prices,
      const std::array<int, 3> min_points_to_add,
      const std::vector<BidderDiagram<double>> &bidder_diags_min,
      const std::vector<BidderDiagram<double>> &bidder_diags_sad,
      const std::vector<BidderDiagram<double>> &bidder_diags_max,
      std::vector<BidderDiagram<double>> &current_bidder_diags_min,
      std::vector<BidderDiagram<double>> &current_bidder_diags_sad,
      std::vector<BidderDiagram<double>> &current_bidder_diags_max) const;

    int wasserstein_{2};
    double geometrical_factor_{1};
    double deltaLim_{0.01};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double lambda_;
    int pairTypeClustering_;
    int numberOfInputs_;
    bool use_kdtree_{true};
    bool useFullDiagrams_{false};
    bool do_min_{true}, do_sad_{true}, do_max_{true};

    std::vector<std::vector<DiagramTuple>> inputDiagramsMin_{};
    std::vector<std::vector<DiagramTuple>> inputDiagramsSaddle_{};
    std::vector<std::vector<DiagramTuple>> inputDiagramsMax_{};

    std::vector<BidderDiagram<double>> bidder_diagrams_min_{};
    std::vector<BidderDiagram<double>> bidder_diagrams_saddle_{};
    std::vector<BidderDiagram<double>> bidder_diagrams_max_{};
    std::vector<BidderDiagram<double>> current_bidder_diagrams_min_{};
    std::vector<BidderDiagram<double>> current_bidder_diagrams_saddle_{};
    std::vector<BidderDiagram<double>> current_bidder_diagrams_max_{};
  };
} // namespace ttk
