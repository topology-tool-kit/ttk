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

#include <Debug.h>

namespace ttk {

  using DiagramTuple = std::tuple<
    /** Vertex Id of low pair element */
    ttk::SimplexId,
    /** Critical Type of low pair element */
    ttk::CriticalType,
    /** Vertex Id of high pair element */
    ttk::SimplexId,
    /** Critical Type of high pair element */
    ttk::CriticalType,
    /** Pair persistence value */
    double,
    /** Pair type */
    ttk::SimplexId,
    /** Pair birth */
    double,
    /** Low pair element 3D coordinates */
    // TODO use std::array<float, 3>
    float,
    float,
    float,
    /** Pair death */
    double,
    /** High pair element 3D coordinates */
    // TODO use std::array<float, 3>
    float,
    float,
    float>;

  using Diagram = std::vector<DiagramTuple>;

  class PersistenceDiagramDistanceMatrix : virtual public Debug {

  public:
    PersistenceDiagramDistanceMatrix() {
      this->setDebugMsgPrefix("PersistenceDiagramDistanceMatrix");
    }

    std::vector<std::vector<double>>
      execute(const std::vector<Diagram> &intermediateDiagrams,
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
    template <typename T>
    double getMostPersistent(const std::vector<T> &bidder_diags) const;
    template <typename T>
    double computePowerDistance(const T &D1, const T &D2) const;
    template <typename T>
    void getDiagramsDistMat(const std::array<size_t, 2> &nInputs,
                            std::vector<std::vector<double>> &distanceMatrix,
                            const std::vector<T> &diags_min,
                            const std::vector<T> &diags_sad,
                            const std::vector<T> &diags_max) const;
    template <typename T>
    void setBidderDiagrams(const size_t nInputs,
                           std::vector<Diagram> &inputDiagrams,
                           std::vector<T> &bidder_diags) const;
    template <typename T>
    void enrichCurrentBidderDiagrams(
      const std::vector<T> &bidder_diags,
      std::vector<T> &current_bidder_diags,
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
