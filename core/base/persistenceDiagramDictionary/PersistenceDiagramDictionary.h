/// \ingroup base
/// \class ttk::PersistenceDiagramDictionary
/// \author Keanu Sisouk <keanu.sisouk@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date Mai 2023
///
/// \b Related \b publication \n
/// "Wasserstein Dictionaries of Persistence Diagrams" \n
/// Keanu Sisouk, Julie Delon and Julien Tierny \n
/// IEEE Transactions on Visualization and Computer Graphics, 2023.
///
/// \sa PersistenceDiagramDictionary

#pragma once

#include "PersistenceDiagramUtils.h"
#include <array>
#include <limits>

#include <ConstrainedGradientDescent.h>
#include <InitDictPersistenceDiagram.h>
#include <InitDictRandomly.h>
#include <PersistenceDiagramAuction.h>
#include <PersistenceDiagramClustering.h>
#include <Wrapper.h>

namespace ttk {
  class PersistenceDiagramDictionary : virtual public Debug {

  public:
    PersistenceDiagramDictionary() {
      this->setDebugMsgPrefix("PersistenceDiagramDictionary");
    }

    void execute(std::vector<ttk::DiagramType> &intermediateDiagrams,
                 const std::vector<ttk::DiagramType> &intermediateAtoms,
                 std::vector<ttk::DiagramType> &dictDiagrams,
                 std::vector<std::vector<double>> &vectorWeights,
                 const int seed,
                 const int numAtom,
                 std::vector<double> &lossTab,
                 std::vector<std::vector<double>> &allLosses,
                 double percent);

    void method(const std::vector<ttk::DiagramType> &intermediateDiagrams,
                std::vector<ttk::DiagramType> &dictDiagrams,
                std::vector<std::vector<double>> &vectorWeights,
                const int numAtom,
                std::vector<double> &lossTab,
                std::vector<std::vector<double>> &allLosses,
                std::vector<std::vector<double>> &histoVectorWeights,
                std::vector<ttk::DiagramType> &histoDictDiagrams,
                bool preWeightOpt,
                double percent,
                bool doCompression);

    enum class BACKEND {
      BORDER_INIT = 0,
      RANDOM_INIT = 1,
      FIRST_DIAGS = 2,
      INPUT_ATOMS = 3,
      GREEDY_INIT = 4
    };

    inline void setUseDimReduct(bool data) {
      DimReductMode_ = data;
    }

    inline void setUseProgApproach(bool data) {
      ProgApproach_ = data;
    }

    inline void setDos(const bool min, const bool sad, const bool max) {
      do_min_ = min;
      do_sad_ = sad;
      do_max_ = max;
    }

    inline void setMinPersistence_(const double data) {
      MinPersistence_ = data;
    }

  protected:
    BACKEND BackEnd{BACKEND::BORDER_INIT};
    double distVect(const std::vector<double> &vec1,
                    const std::vector<double> &vec2) const;

    double
      getMostPersistent(const std::vector<BidderDiagram> &bidder_diags) const;
    double computeDistance(const BidderDiagram &D1,
                           const BidderDiagram &D2,
                           std::vector<ttk::MatchingType> &matching) const;

    void computeGradientWeights(
      std::vector<double> &gradWeights,
      std::vector<Matrix> &hessianList,
      const std::vector<ttk::DiagramType> &dictDiagrams,
      const std::vector<std::vector<ttk::MatchingType>> &matchingsAtoms,
      const ttk::DiagramType &Barycenter,
      const ttk::DiagramType &newData,
      const std::vector<ttk::MatchingType> &matchingsMin,
      const std::vector<ttk::MatchingType> &matchingsMax,
      const std::vector<ttk::MatchingType> &matchingsSad,
      const std::vector<size_t> &indexBaryMin,
      const std::vector<size_t> &indexBaryMax,
      const std::vector<size_t> &indexBarySad,
      const std::vector<size_t> &indexDataMin,
      const std::vector<size_t> &indexDataMax,
      const std::vector<size_t> &indexDataSad,
      const bool doOptimizeAtoms) const;

    void computeGradientAtoms(
      std::vector<Matrix> &gradsAtoms,
      const std::vector<double> &weights,
      const ttk::DiagramType &Barycenter,
      const ttk::DiagramType &newData,
      const std::vector<ttk::MatchingType> &matchingsMin,
      const std::vector<ttk::MatchingType> &matchingsMax,
      const std::vector<ttk::MatchingType> &matchingsSad,
      const std::vector<size_t> &indexBaryMin,
      const std::vector<size_t> &indexBaryMax,
      const std::vector<size_t> &indexBarySad,
      const std::vector<size_t> &indexDataMin,
      const std::vector<size_t> &indexDataMax,
      const std::vector<size_t> &indexDataSad,
      std::vector<int> &checker,
      std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
      ttk::DiagramType &infoToAdd,
      bool doDimReduct) const;

    void computeDirectionsGradWeight(
      const std::vector<std::vector<ttk::MatchingType>> &matchingsAtoms,
      const ttk::DiagramType &Barycenter,
      const ttk::DiagramType &newData,
      const std::vector<ttk::MatchingType> &matchingsCritType,
      const std::vector<size_t> &indexBaryCritType,
      const std::vector<size_t> &indexDataCritType,
      std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
      std::vector<std::array<double, 2>> &directions,
      std::vector<std::array<double, 2>> &data_assigned,
      std::vector<int> &tracker2,
      const bool doOptimizeAtoms) const;

    void computeDirectionsGradAtoms(
      std::vector<Matrix> &gradsAtoms,
      const ttk::DiagramType &Barycenter,
      const std::vector<double> &weights,
      const ttk::DiagramType &newData,
      const std::vector<ttk::MatchingType> &matchingsCritType,
      const std::vector<size_t> &indexBaryCritType,
      const std::vector<size_t> &indexDataCritType,
      std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
      std::vector<std::vector<double>> &directions,
      std::vector<int> &checker,
      std::vector<PersistencePair> &infoToAdd,
      const bool doOptimizeAtoms) const;

    void setBidderDiagrams(const size_t nInputs,
                           std::vector<ttk::DiagramType> &inputDiagrams,
                           std::vector<BidderDiagram> &bidder_diags) const;

    int initDictionary(std::vector<ttk::DiagramType> &dictDiagrams,
                       const std::vector<ttk::DiagramType> &datas,
                       const std::vector<ttk::DiagramType> &inputAtoms,
                       const int nbAtom,
                       bool do_min_,
                       bool do_sad_,
                       bool do_max_,
                       int seed,
                       double percent);

    void gettingBidderDiagrams(
      const std::vector<ttk::DiagramType> &intermediateDiagrams,
      std::vector<ttk::DiagramType> &inputDiagramsMin,
      std::vector<ttk::DiagramType> &inputDiagramsSad,
      std::vector<ttk::DiagramType> &inputDiagramsMax,
      std::vector<BidderDiagram> &bidderDiagramsMin,
      std::vector<BidderDiagram> &bidderDiagramsSad,
      std::vector<BidderDiagram> &bidderDiagramsMax,
      std::vector<std::vector<size_t>> &originIndexMin,
      std::vector<std::vector<size_t>> &originIndexSad,
      std::vector<std::vector<size_t>> &originIndexMax,
      bool insertOriginIndexMode) const;

    void computeAllDistances(
      std::vector<ttk::DiagramType> &barycentersList,
      const size_t nDiag,
      std::vector<ttk::DiagramType> &barycentersListMin,
      std::vector<ttk::DiagramType> &barycentersListSad,
      std::vector<ttk::DiagramType> &barycentersListMax,
      std::vector<BidderDiagram> &bidderBarycentersListMin,
      std::vector<BidderDiagram> &bidderBarycentersListSad,
      std::vector<BidderDiagram> &bidderBarycentersListMax,
      std::vector<std::vector<size_t>> &originIndexBarysMin,
      std::vector<std::vector<size_t>> &originIndexBarysSad,
      std::vector<std::vector<size_t>> &originIndexBarysMax,
      std::vector<BidderDiagram> &bidderDiagramsMin,
      std::vector<BidderDiagram> &bidderDiagramsMax,
      std::vector<BidderDiagram> &bidderDiagramsSad,
      std::vector<std::vector<ttk::MatchingType>> &matchingsDatasMin,
      std::vector<std::vector<ttk::MatchingType>> &matchingsDatasMax,
      std::vector<std::vector<ttk::MatchingType>> &matchingsDatasSad,
      std::vector<double> &allLossesAtEpoch,
      bool firstDistComputation) const;

    void controlAtomsSize(
      const std::vector<ttk::DiagramType> &intermediateDiagrams,
      std::vector<ttk::DiagramType> &dictDiagrams) const;

    double getMaxPers(const ttk::DiagramType &data);

    int Wasserstein{2};
    double Alpha{1.0};
    double DeltaLim{0.01};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double Lambda{1.0};
    size_t MaxNumberOfPairs{20};
    double MinPersistence_{0.1};

    double CompressionFactor{1.5};
    bool do_min_{true}, do_sad_{true}, do_max_{true};

    int maxLag2_;

    int MaxEpoch_;
    bool MaxEigenValue_{true};
    bool OptimizeWeights_{true};
    bool OptimizeAtoms_{true};

    bool CreationFeatures_{true};
    bool Fusion_{false};
    bool ProgBarycenter_{false};

    bool sortedForTest_{false};
    bool ProgApproach_{false};
    bool StopCondition_{true};

    bool CompressionMode_{false};
    bool DimReductMode_{false};
  };
} // namespace ttk
