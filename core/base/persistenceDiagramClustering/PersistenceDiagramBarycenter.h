/// \ingroup base
/// \class ttk::PersistenceDiagramBarycenter
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

// base code includes
#include <KDTree.h>
#include <PDBarycenter.h>
#include <PersistenceDiagramAuction.h>
#include <Wrapper.h>

namespace ttk {

  class PersistenceDiagramBarycenter : public Debug {

  public:
    PersistenceDiagramBarycenter() {
      this->threadNumber_ = 1;
      this->setDebugMsgPrefix("PersistenceDiagramBarycenter");
    }

    ~PersistenceDiagramBarycenter() override = default;

    void execute(
      std::vector<DiagramType> &intermediateDiagrams,
      DiagramType &barycenter,
      std::vector<std::vector<std::vector<MatchingType>>> &all_matchings);

    inline void setNumberOfInputs(const int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
    }

    inline void setDeterministic(const bool deterministic) {
      deterministic_ = deterministic;
    }

    inline void setWasserstein(const std::string &wasserstein) {
      wasserstein_ = (wasserstein == "inf") ? -1 : stoi(wasserstein);
    }

    inline void setUseProgressive(const bool use_progressive) {
      if(use_progressive)
        epsilon_decreases_ = true;
      use_progressive_ = use_progressive;
    }

    inline void setAlpha(const double alpha) {
      alpha_ = alpha;
    }

    inline void setLambda(const double lambda) {
      lambda_ = lambda;
    }

    inline void setMethod(const int &method) {
      method_ = method;
    }

    inline void setReinitPrices(const bool reinit_prices) {
      reinit_prices_ = reinit_prices;
    }

    inline void setEpsilonDecreases(const bool epsilon_decreases) {
      if(use_progressive_)
        epsilon_decreases_ = true;
      else
        epsilon_decreases_ = epsilon_decreases;
    }

    inline void setEarlyStoppage(const bool early_stoppage) {
      early_stoppage_ = early_stoppage;
    }

    inline void setNonMatchingWeight(const double nonMatchingWeight) {
      nonMatchingWeight_ = nonMatchingWeight;
    }

    inline void setDeltaLim(double deltaLim) {
      delta_lim_ = deltaLim;
    }

  protected:
    bool deterministic_{true};
    int method_;
    int wasserstein_{2};
    int numberOfInputs_{0};
    bool use_progressive_{true};
    double alpha_{1.0};
    double lambda_{1.0};
    double nonMatchingWeight_ = 1.0;
    double delta_lim_{0.01};

    int points_added_;
    int points_deleted_;

    std::vector<std::vector<double>> all_matchings_;
    std::vector<std::vector<double>> all_old_matchings_;
    std::vector<BidderDiagram> bidder_diagrams_;
    std::vector<GoodDiagram> barycenter_goods_;

    bool reinit_prices_{true};
    bool epsilon_decreases_{true};
    bool early_stoppage_;
  };

} // namespace ttk
