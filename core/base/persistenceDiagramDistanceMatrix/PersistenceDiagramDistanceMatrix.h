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

  class PersistenceDiagramDistanceMatrix : virtual public Debug {

  public:
    PersistenceDiagramDistanceMatrix() {
      this->setDebugMsgPrefix("PersistenceDiagramDistanceMatrix");
    }

    std::vector<std::vector<double>> execute(
      std::vector<std::vector<DiagramTuple>> &intermediateDiagrams) const;

    inline void setWasserstein(const int data) {
      Wasserstein = data;
    }
    inline void setDos(const bool min, const bool sad, const bool max) {
      do_min_ = min;
      do_sad_ = sad;
      do_max_ = max;
    }
    inline void setAlpha(const double alpha) {
      Alpha = alpha;
    }
    inline void setLambda(const double lambda) {
      Lambda = lambda;
    }
    inline void setDeltaLim(const double deltaLim) {
      DeltaLim = deltaLim;
    }
    inline void setUseFullDiagrams(const bool arg) {
      UseFullDiagrams = arg;
    }
    inline void setMinPointsToAdd(const size_t data) {
      min_points_to_add_ = data;
    }

  protected:
    double getMostPersistent(
      const std::vector<BidderDiagram<double>> &bidder_diags) const;
    double computeDistance(const BidderDiagram<double> &D1,
                           const BidderDiagram<double> &D2) const;
    void getDiagramsDistMat(
      const size_t nInputs,
      std::vector<std::vector<double>> &distanceMatrix,
      const std::vector<BidderDiagram<double>> &diags_min,
      const std::vector<BidderDiagram<double>> &diags_sad,
      const std::vector<BidderDiagram<double>> &diags_max) const;
    void setBidderDiagrams(
      const size_t nInputs,
      std::vector<std::vector<DiagramTuple>> &inputDiagrams,
      std::vector<BidderDiagram<double>> &bidder_diags,
      std::vector<BidderDiagram<double>> &current_bidder_diags) const;
    double enrichCurrentBidderDiagrams(
      const double prev_min_persistence,
      const size_t min_points_to_add,
      const std::vector<BidderDiagram<double>> &bidder_diags,
      std::vector<BidderDiagram<double>> &current_bidder_diags) const;

    int Wasserstein{2};
    double Alpha{1.0};
    double DeltaLim{0.01};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double Lambda;
    size_t min_points_to_add_{10};
    bool UseFullDiagrams{false};
    bool do_min_{true}, do_sad_{true}, do_max_{true};
  };
} // namespace ttk
