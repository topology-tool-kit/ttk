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
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n

#pragma once

#include <array>

#include <Debug.h>
#include <PersistenceDiagramAuction.h>
#include <PersistenceDiagramUtils.h>

namespace ttk {

  class PersistenceDiagramDistanceMatrix : virtual public Debug {

  public:
    PersistenceDiagramDistanceMatrix() {
      this->setDebugMsgPrefix("PersistenceDiagramDistanceMatrix");
    }

    std::vector<std::vector<double>>
      execute(const std::vector<DiagramType> &intermediateDiagrams,
              const std::array<size_t, 2> &nInputs) const;

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
    inline void setMaxNumberOfPairs(const size_t data) {
      MaxNumberOfPairs = data;
    }
    inline void setMinPersistence(const double data) {
      MinPersistence = data;
    }
    inline void setConstraint(const int data) {
      if(data == 0) {
        this->Constraint = ConstraintType::FULL_DIAGRAMS;
      } else if(data == 1) {
        this->Constraint = ConstraintType::NUMBER_PAIRS;
      } else if(data == 2) {
        this->Constraint = ConstraintType::ABSOLUTE_PERSISTENCE;
      } else if(data == 3) {
        this->Constraint = ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG;
      } else if(data == 4) {
        this->Constraint = ConstraintType::RELATIVE_PERSISTENCE_GLOBAL;
      }
    }

  protected:
    double
      getMostPersistent(const std::vector<BidderDiagram> &bidder_diags) const;
    double computePowerDistance(const BidderDiagram &D1,
                                const BidderDiagram &D2) const;
    void getDiagramsDistMat(const std::array<size_t, 2> &nInputs,
                            std::vector<std::vector<double>> &distanceMatrix,
                            const std::vector<BidderDiagram> &diags_min,
                            const std::vector<BidderDiagram> &diags_sad,
                            const std::vector<BidderDiagram> &diags_max) const;
    void setBidderDiagrams(const size_t nInputs,
                           std::vector<DiagramType> &inputDiagrams,
                           std::vector<BidderDiagram> &bidder_diags) const;

    void enrichCurrentBidderDiagrams(
      const std::vector<BidderDiagram> &bidder_diags,
      std::vector<BidderDiagram> &current_bidder_diags,
      const std::vector<double> &maxDiagPersistence) const;

    int Wasserstein{2};
    double Alpha{1.0};
    double DeltaLim{0.01};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double Lambda;
    size_t MaxNumberOfPairs{20};
    double MinPersistence{0.1};
    bool do_min_{true}, do_sad_{true}, do_max_{true};

    enum class ConstraintType {
      FULL_DIAGRAMS,
      NUMBER_PAIRS,
      ABSOLUTE_PERSISTENCE,
      RELATIVE_PERSISTENCE_PER_DIAG,
      RELATIVE_PERSISTENCE_GLOBAL,
    };
    ConstraintType Constraint{ConstraintType::RELATIVE_PERSISTENCE_GLOBAL};
  };
} // namespace ttk
