/// \ingroup base
/// \class ttk::TopologicalOptimization
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Mohamed Amine Kissi <mohamed.kissi@lip6.fr>
/// \date March 2024

#pragma once

#ifdef TTK_ENABLE_TORCH
#include <torch/optim.h>
#include <torch/torch.h>
#endif

// base code includes
#include <Debug.h>
#include <PersistenceDiagram.h>
#include <PersistenceDiagramClustering.h>
#include <Triangulation.h>

namespace ttk {

  class TopologicalOptimization : virtual public Debug {
  public:
    TopologicalOptimization();

    template <typename dataType, typename triangulationType>
    int execute(const dataType *const inputScalars,
                dataType *const outputScalars,
                SimplexId *const inputOffsets,
                triangulationType *triangulation,
                const ttk::DiagramType &constraintDiagram) const;

    inline int preconditionTriangulation(AbstractTriangulation *triangulation) {
      if(triangulation) {
        vertexNumber_ = triangulation->getNumberOfVertices();
        triangulation->preconditionVertexNeighbors();
      }
      return 0;
    }

    /*
      This function allows us to retrieve the indices of the critical points
      that we must modify in order to match our current diagram to our target
      diagram.
    */
    template <typename dataType, typename triangulationType>
    void getIndices(
      triangulationType *triangulation,
      SimplexId *&inputOffsets,
      dataType *const inputScalars,
      const ttk::DiagramType &constraintDiagram,
      int epoch,
      std::vector<SimplexId> &listAllIndicesToChange,
      std::vector<std::vector<SimplexId>> &pair2MatchedPair,
      std::vector<std::vector<SimplexId>> &pair2Delete,
      std::vector<SimplexId> &pairChangeMatchingPair,
      std::vector<SimplexId> &birthPairToDeleteCurrentDiagram,
      std::vector<double> &birthPairToDeleteTargetDiagram,
      std::vector<SimplexId> &deathPairToDeleteCurrentDiagram,
      std::vector<double> &deathPairToDeleteTargetDiagram,
      std::vector<SimplexId> &birthPairToChangeCurrentDiagram,
      std::vector<double> &birthPairToChangeTargetDiagram,
      std::vector<SimplexId> &deathPairToChangeCurrentDiagram,
      std::vector<double> &deathPairToChangeTargetDiagram,
      std::vector<std::vector<SimplexId>> &currentVertex2PairsCurrentDiagram,
      std::vector<int> &vertexInHowManyPairs) const;

/*
  This function allows you to copy the values of a pytorch tensor
  to a vector in an optimized way.
*/
#ifdef TTK_ENABLE_TORCH
    int tensorToVectorFast(const torch::Tensor &tensor,
                           std::vector<double> &result) const;
#endif

    inline void setUseFastPersistenceUpdate(bool UseFastPersistenceUpdate) {
      useFastPersistenceUpdate_ = UseFastPersistenceUpdate;
    }

    inline void setFastAssignmentUpdate(bool FastAssignmentUpdate) {
      fastAssignmentUpdate_ = FastAssignmentUpdate;
    }

    inline void setEpochNumber(int EpochNumber) {
      epochNumber_ = EpochNumber;
    }

    inline void setPDCMethod(int PDCMethod) {
      pdcMethod_ = PDCMethod;
    }

    inline void setMethodOptimization(int methodOptimization) {
      methodOptimization_ = methodOptimization;
    }

    inline void setFinePairManagement(int finePairManagement) {
      finePairManagement_ = finePairManagement;
    }

    inline void setChooseLearningRate(int chooseLearningRate) {
      chooseLearningRate_ = chooseLearningRate;
    }

    inline void setLearningRate(double learningRate) {
      learningRate_ = learningRate;
    }

    inline void setAlpha(double alpha) {
      alpha_ = alpha;
    }

    inline void setCoefStopCondition(double coefStopCondition) {
      coefStopCondition_ = coefStopCondition;
    }

    inline void
      setOptimizationWithoutMatching(bool optimizationWithoutMatching) {
      optimizationWithoutMatching_ = optimizationWithoutMatching;
    }

    inline void setThresholdMethod(int thresholdMethod) {
      thresholdMethod_ = thresholdMethod;
    }

    inline void setThresholdPersistence(double thresholdPersistence) {
      thresholdPersistence_ = thresholdPersistence;
    }

    inline void setLowerThreshold(int lowerThreshold) {
      lowerThreshold_ = lowerThreshold;
    }

    inline void setUpperThreshold(int upperThreshold) {
      upperThreshold_ = upperThreshold;
    }

    inline void setPairTypeToDelete(int pairTypeToDelete) {
      pairTypeToDelete_ = pairTypeToDelete;
    }

    inline void setConstraintAveraging(bool ConstraintAveraging) {
      constraintAveraging_ = ConstraintAveraging;
    }

    inline void setPrintFrequency(int printFrequency) {
      printFrequency_ = printFrequency;
    }

  protected:
    SimplexId vertexNumber_{};
    int epochNumber_;

    // enable the fast update of the persistence diagram
    bool useFastPersistenceUpdate_;

    // enable the fast update of the pair assignments between the target diagram
    bool fastAssignmentUpdate_;

    // if pdcMethod_ == 0 then we use Progressive approach
    // if pdcMethod_ == 1 then we use Classical Auction approach
    int pdcMethod_;

    // if methodOptimization_ == 0 then we use Direct gradient descent
    // if methodOptimization_ == 1 then we use Adam
    int methodOptimization_;

    // if finePairManagement_ == 0 then we let the algorithm choose
    // if finePairManagement_ == 1 then we fill the domain
    // if finePairManagement_ == 2 then we cut the domain
    int finePairManagement_;

    // Adam
    bool chooseLearningRate_;
    double learningRate_;

    // Direct gradient descent
    // alpha_ : the gradient step size
    double alpha_;

    // Stopping criterion: when the loss becomes less than a percentage
    // coefStopCondition_ (e.g. coefStopCondition_ = 0.01 => 1%) of the original
    // loss (between input diagram and simplified diagram)
    double coefStopCondition_;

    // Optimization without matching (OWM)
    bool optimizationWithoutMatching_;

    // [OWM] if thresholdMethod_ == 0 : threshold on persistence
    // [OWM] if thresholdMethod_ == 1 : threshold on pair type
    int thresholdMethod_;

    // [OWM] thresholdPersistence_ : The threshold value on persistence.
    double thresholdPersistence_;

    // [OWM] lowerThreshold_ : The lower threshold on pair type
    int lowerThreshold_;

    // [OWM] upperThreshold_ : The upper threshold on pair type
    int upperThreshold_;

    // [OWM] pairTypeToDelete_ : Remove only pairs of type pairTypeToDelete_
    int pairTypeToDelete_;

    bool constraintAveraging_;

    int printFrequency_{10};
  };

} // namespace ttk

#ifdef TTK_ENABLE_TORCH
class PersistenceGradientDescent : public torch::nn::Module,
                                   public ttk::TopologicalOptimization {
public:
  PersistenceGradientDescent(torch::Tensor X_tensor) : torch::nn::Module() {
    X = register_parameter("X", X_tensor, true);
  }
  torch::Tensor X;
};

#endif

/*
  This function allows us to retrieve the indices of the critical points
  that we must modify in order to match our current diagram to our target
  diagram.
*/
template <typename dataType, typename triangulationType>
void ttk::TopologicalOptimization::getIndices(
  triangulationType *triangulation,
  SimplexId *&inputOffsets,
  dataType *const inputScalars,
  const ttk::DiagramType &constraintDiagram,
  int epoch,
  std::vector<SimplexId> &listAllIndicesToChange,
  std::vector<std::vector<SimplexId>> &pair2MatchedPair,
  std::vector<std::vector<SimplexId>> &pair2Delete,
  std::vector<SimplexId> &pairChangeMatchingPair,
  std::vector<SimplexId> &birthPairToDeleteCurrentDiagram,
  std::vector<double> &birthPairToDeleteTargetDiagram,
  std::vector<SimplexId> &deathPairToDeleteCurrentDiagram,
  std::vector<double> &deathPairToDeleteTargetDiagram,
  std::vector<SimplexId> &birthPairToChangeCurrentDiagram,
  std::vector<double> &birthPairToChangeTargetDiagram,
  std::vector<SimplexId> &deathPairToChangeCurrentDiagram,
  std::vector<double> &deathPairToChangeTargetDiagram,
  std::vector<std::vector<SimplexId>> &currentVertex2PairsCurrentDiagram,
  std::vector<int> &vertexInHowManyPairs) const {

  //=========================================
  //            Lazy Gradient
  //=========================================

  bool needUpdateDefaultValue
    = (useFastPersistenceUpdate_ ? (epoch == 0 || epoch < 0 ? true : false)
                                 : true);
  std::vector<bool> needUpdate(vertexNumber_, needUpdateDefaultValue);
  if(useFastPersistenceUpdate_) {
    /*
      There is a 10% loss of performance
    */
    this->printMsg(
      "Get Indices | UseFastPersistenceUpdate_", debug::Priority::DETAIL);

    if(not(epoch == 0 || epoch < 0)) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(size_t index = 0; index < listAllIndicesToChange.size(); index++) {
        if(listAllIndicesToChange[index] == 1) {
          needUpdate[index] = true;

          // Find all the neighbors of the vertex
          int vertexNumber = triangulation->getVertexNeighborNumber(index);
          for(int i = 0; i < vertexNumber; i++) {
            SimplexId vertexNeighborId = -1;
            triangulation->getVertexNeighbor(index, i, vertexNeighborId);
            needUpdate[vertexNeighborId] = true;
          }
        }
      }
    }
  }

  SimplexId count = std::count(needUpdate.begin(), needUpdate.end(), true);

  this->printMsg(
    "Get Indices | The number of vertices that need to be updated is: "
      + std::to_string(count),
    debug::Priority::DETAIL);

  //=========================================
  //     Compute the persistence diagram
  //=========================================
  ttk::Timer timePersistenceDiagram;

  ttk::PersistenceDiagram diagram;
  std::vector<ttk::PersistencePair> diagramOutput;
  ttk::preconditionOrderArray<dataType>(
    vertexNumber_, inputScalars, inputOffsets, threadNumber_);
  diagram.setDebugLevel(debugLevel_);
  diagram.setThreadNumber(threadNumber_);
  diagram.preconditionTriangulation(triangulation);

  if(useFastPersistenceUpdate_) {
    diagram.execute(
      diagramOutput, inputScalars, 0, inputOffsets, triangulation, &needUpdate);
  } else {
    diagram.execute(
      diagramOutput, inputScalars, epoch, inputOffsets, triangulation);
  }

  //=====================================
  //          Matching Pairs
  //=====================================

  if(optimizationWithoutMatching_) {
    for(SimplexId i = 0; i < (SimplexId)diagramOutput.size(); i++) {
      auto pair = diagramOutput[i];
      if((thresholdMethod_ == 0)
         && (pair.persistence() < thresholdPersistence_)) {
        birthPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.birth.id));
        birthPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
        deathPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.death.id));
        deathPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
      } else if((thresholdMethod_ == 1)
                && ((pair.dim < lowerThreshold_)
                    || (pair.dim > upperThreshold_))) {
        birthPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.birth.id));
        birthPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
        deathPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.death.id));
        deathPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
      } else if((thresholdMethod_ == 2) && (pair.dim == pairTypeToDelete_)) {
        birthPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.birth.id));
        birthPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
        deathPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.death.id));
        deathPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
      }
    }
  } else if(fastAssignmentUpdate_) {

    std::vector<std::vector<SimplexId>> vertex2PairsCurrentDiagram(
      vertexNumber_, std::vector<SimplexId>());
    for(SimplexId i = 0; i < (SimplexId)diagramOutput.size(); i++) {
      auto &pair = diagramOutput[i];
      vertex2PairsCurrentDiagram[pair.birth.id].push_back(i);
      vertex2PairsCurrentDiagram[pair.death.id].push_back(i);
      vertexInHowManyPairs[pair.birth.id]++;
      vertexInHowManyPairs[pair.death.id]++;
    }

    std::vector<std::vector<SimplexId>> vertex2PairsTargetDiagram(
      vertexNumber_, std::vector<SimplexId>());
    for(SimplexId i = 0; i < (SimplexId)constraintDiagram.size(); i++) {
      auto &pair = constraintDiagram[i];
      vertex2PairsTargetDiagram[pair.birth.id].push_back(i);
      vertex2PairsTargetDiagram[pair.death.id].push_back(i);
    }

    std::vector<std::vector<SimplexId>> matchedPairs;
    for(SimplexId i = 0; i < (SimplexId)constraintDiagram.size(); i++) {
      auto &pair = constraintDiagram[i];

      SimplexId birthId = -1;
      SimplexId deathId = -1;

      if(pairChangeMatchingPair[i] == 1) {
        birthId = pair2MatchedPair[i][0];
        deathId = pair2MatchedPair[i][1];
      } else {
        birthId = pair.birth.id;
        deathId = pair.death.id;
      }

      if(epoch == 0) {
        for(auto &idPairBirth : vertex2PairsCurrentDiagram[birthId]) {
          for(auto &idPairDeath : vertex2PairsCurrentDiagram[deathId]) {
            if(idPairBirth == idPairDeath) {
              matchedPairs.push_back({i, idPairBirth});
            }
          }
        }
      } else if((vertex2PairsCurrentDiagram[birthId].size() == 1)
                && (vertex2PairsCurrentDiagram[deathId].size() == 1)) {
        if(vertex2PairsCurrentDiagram[birthId][0]
           == vertex2PairsCurrentDiagram[deathId][0]) {
          matchedPairs.push_back({i, vertex2PairsCurrentDiagram[deathId][0]});
        }
      }
    }

    std::vector<SimplexId> matchingPairCurrentDiagram(
      (SimplexId)diagramOutput.size(), -1);
    std::vector<SimplexId> matchingPairTargetDiagram(
      (SimplexId)constraintDiagram.size(), -1);

    for(auto &match : matchedPairs) {
      auto &indicePairTargetDiagram = match[0];
      auto &indicePairCurrentDiagram = match[1];

      auto &pairCurrentDiagram = diagramOutput[indicePairCurrentDiagram];
      auto &pairTargetDiagram = constraintDiagram[indicePairTargetDiagram];

      pair2MatchedPair[indicePairTargetDiagram][0]
        = pairCurrentDiagram.birth.id;
      pair2MatchedPair[indicePairTargetDiagram][1]
        = pairCurrentDiagram.death.id;

      matchingPairCurrentDiagram[indicePairCurrentDiagram] = 1;
      matchingPairTargetDiagram[indicePairTargetDiagram] = 1;

      SimplexId valueBirthPairToChangeCurrentDiagram
        = (SimplexId)(pairCurrentDiagram.birth.id);
      SimplexId valueDeathPairToChangeCurrentDiagram
        = (SimplexId)(pairCurrentDiagram.death.id);

      double valueBirthPairToChangeTargetDiagram
        = pairTargetDiagram.birth.sfValue;
      double valueDeathPairToChangeTargetDiagram
        = pairTargetDiagram.death.sfValue;

      birthPairToChangeCurrentDiagram.push_back(
        valueBirthPairToChangeCurrentDiagram);
      birthPairToChangeTargetDiagram.push_back(
        valueBirthPairToChangeTargetDiagram);
      deathPairToChangeCurrentDiagram.push_back(
        valueDeathPairToChangeCurrentDiagram);
      deathPairToChangeTargetDiagram.push_back(
        valueDeathPairToChangeTargetDiagram);
    }

    ttk::DiagramType thresholdCurrentDiagram{};
    for(SimplexId i = 0; i < (SimplexId)diagramOutput.size(); i++) {
      auto &pair = diagramOutput[i];

      if((pair2Delete[pair.birth.id].size() == 1)
         && (pair2Delete[pair.death.id].size() == 1)
         && (pair2Delete[pair.birth.id] == pair2Delete[pair.death.id])) {

        birthPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.birth.id));
        birthPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
        deathPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.death.id));
        deathPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
        continue;
      }
      if(matchingPairCurrentDiagram[i] == -1) {
        thresholdCurrentDiagram.push_back(pair);
      }
    }

    ttk::DiagramType thresholdConstraintDiagram{};
    std::vector<SimplexId> pairIndiceLocal2Global{};
    for(SimplexId i = 0; i < (SimplexId)constraintDiagram.size(); i++) {
      auto &pair = constraintDiagram[i];

      if(matchingPairTargetDiagram[i] == -1) {
        thresholdConstraintDiagram.push_back(pair);
        pairIndiceLocal2Global.push_back(i);
      }
    }

    this->printMsg("Get Indices | thresholdCurrentDiagram.size(): "
                     + std::to_string(thresholdCurrentDiagram.size()),
                   debug::Priority::DETAIL);

    this->printMsg("Get Indices | thresholdConstraintDiagram.size(): "
                     + std::to_string(thresholdConstraintDiagram.size()),
                   debug::Priority::DETAIL);

    if(thresholdConstraintDiagram.size() == 0) {
      for(SimplexId i = 0; i < (SimplexId)thresholdCurrentDiagram.size(); i++) {
        auto &pair = thresholdCurrentDiagram[i];

        if(!constraintAveraging_) {

          // If the point pair.birth.id is in a signal pair
          // AND If the point pair.death.id is not in a signal pair
          // Then we only modify the pair.death.id
          if((vertex2PairsTargetDiagram[pair.birth.id].size() >= 1)
             && (vertex2PairsTargetDiagram[pair.death.id].size() == 0)) {
            deathPairToDeleteCurrentDiagram.push_back(
              static_cast<SimplexId>(pair.death.id));
            deathPairToDeleteTargetDiagram.push_back(
              (pair.birth.sfValue + pair.death.sfValue) / 2);
            continue;
          }

          // If the point pair.death.id is in a signal pair
          // AND If the point pair.birth.id is not in a signal pair
          // Then we only modify the pair.birth.id
          if((vertex2PairsTargetDiagram[pair.birth.id].size() == 0)
             && (vertex2PairsTargetDiagram[pair.death.id].size() >= 1)) {
            birthPairToDeleteCurrentDiagram.push_back(
              static_cast<SimplexId>(pair.birth.id));
            birthPairToDeleteTargetDiagram.push_back(
              (pair.birth.sfValue + pair.death.sfValue) / 2);
            continue;
          }

          // If the point pair.birth.id is in a signal pair
          // AND If the point pair.death.id is in a signal pair
          // Then we do not modify either point
          if((vertex2PairsTargetDiagram[pair.birth.id].size() >= 1)
             || (vertex2PairsTargetDiagram[pair.death.id].size() >= 1)) {
            continue;
          }
        }

        birthPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.birth.id));
        birthPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
        deathPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.death.id));
        deathPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);

        pair2Delete[pair.birth.id].push_back(i);
        pair2Delete[pair.death.id].push_back(i);
      }
    } else {

      ttk::Timer timePersistenceDiagramClustering;

      ttk::PersistenceDiagramClustering persistenceDiagramClustering;
      PersistenceDiagramBarycenter pdBarycenter{};
      std::vector<ttk::DiagramType> intermediateDiagrams{
        thresholdConstraintDiagram, thresholdCurrentDiagram};
      std::vector<std::vector<std::vector<ttk::MatchingType>>> allMatchings;
      std::vector<ttk::DiagramType> centroids{};

      if(pdcMethod_ == 0) {
        persistenceDiagramClustering.setDebugLevel(debugLevel_);
        persistenceDiagramClustering.setThreadNumber(threadNumber_);
        // setDeterministic ==> Deterministic algorithm
        persistenceDiagramClustering.setDeterministic(true);
        // setUseProgressive ==> Compute Progressive Barycenter
        persistenceDiagramClustering.setUseProgressive(true);
        // setUseInterruptible ==> Interruptible algorithm
        persistenceDiagramClustering.setUseInterruptible(false);
        // // setTimeLimit ==> Maximal computation time (s)
        persistenceDiagramClustering.setTimeLimit(0.01);
        // setUseAdditionalPrecision ==> Force minimum precision on matchings
        persistenceDiagramClustering.setUseAdditionalPrecision(true);
        // setDeltaLim ==> Minimal relative precision
        persistenceDiagramClustering.setDeltaLim(1e-5);
        // setUseAccelerated ==> Use Accelerated KMeans
        persistenceDiagramClustering.setUseAccelerated(false);
        // setUseKmeansppInit ==> KMeanspp Initialization
        persistenceDiagramClustering.setUseKmeansppInit(false);

        std::vector<int> clusterIds = persistenceDiagramClustering.execute(
          intermediateDiagrams, centroids, allMatchings);
      } else {

        centroids.resize(1);
        const auto wassersteinMetric = std::to_string(2);
        pdBarycenter.setWasserstein(wassersteinMetric);
        pdBarycenter.setMethod(2);
        pdBarycenter.setNumberOfInputs(2);
        pdBarycenter.setDeterministic(true);
        pdBarycenter.setUseProgressive(true);
        pdBarycenter.setDebugLevel(debugLevel_);
        pdBarycenter.setThreadNumber(threadNumber_);
        pdBarycenter.setAlpha(1);
        pdBarycenter.setLambda(1);
        pdBarycenter.execute(intermediateDiagrams, centroids[0], allMatchings);
      }

      std::vector<std::vector<SimplexId>> allPairsSelected{};
      std::vector<std::vector<SimplexId>> matchingsBlockPairs(
        centroids[0].size());

      for(auto i = 1; i >= 0; --i) {
        std::vector<ttk::MatchingType> &matching = allMatchings[0][i];

        const auto &diag{intermediateDiagrams[i]};

        for(SimplexId j = 0; j < (SimplexId)matching.size(); j++) {

          const auto &m{matching[j]};
          const auto &bidderId{std::get<0>(m)};
          const auto &goodId{std::get<1>(m)};

          if((goodId == -1) | (bidderId == -1)) {
            continue;
          }

          if(diag[bidderId].persistence() != 0) {
            if(i == 1) {
              matchingsBlockPairs[goodId].push_back(bidderId);
            } else if(matchingsBlockPairs[goodId].size() > 0) {
              matchingsBlockPairs[goodId].push_back(bidderId);
            }
            allPairsSelected.push_back(
              {diag[bidderId].birth.id, diag[bidderId].death.id});
          }
        }
      }

      std::vector<ttk::PersistencePair> pairsToErase{};

      std::map<std::vector<SimplexId>, SimplexId> currentToTarget;
      for(auto &pair : allPairsSelected) {
        currentToTarget[{pair[0], pair[1]}] = 1;
      }

      for(auto &pair : intermediateDiagrams[1]) {
        if(pair.isFinite != 0) {
          if(!(currentToTarget.count({pair.birth.id, pair.death.id}) > 0)) {
            pairsToErase.push_back(pair);
          }
        }
      }

      for(auto &pair : pairsToErase) {

        if(!constraintAveraging_) {

          // If the point pair.birth.id is in a signal pair
          // AND If the point pair.death.id is not in a signal pair
          // Then we only modify the pair.death.id
          if((vertex2PairsTargetDiagram[pair.birth.id].size() >= 1)
             && (vertex2PairsTargetDiagram[pair.death.id].size() == 0)) {
            deathPairToDeleteCurrentDiagram.push_back(
              static_cast<SimplexId>(pair.death.id));
            deathPairToDeleteTargetDiagram.push_back(
              (pair.birth.sfValue + pair.death.sfValue) / 2);
            continue;
          }

          // If the point pair.death.id is in a signal pair
          // AND If the point pair.birth.id is not in a signal pair
          // Then we only modify the pair.birth.id
          if((vertex2PairsTargetDiagram[pair.birth.id].size() == 0)
             && (vertex2PairsTargetDiagram[pair.death.id].size() >= 1)) {
            birthPairToDeleteCurrentDiagram.push_back(
              static_cast<SimplexId>(pair.birth.id));
            birthPairToDeleteTargetDiagram.push_back(
              (pair.birth.sfValue + pair.death.sfValue) / 2);
            continue;
          }

          // If the point pair.birth.id is in a signal pair
          // AND If the point pair.death.id is in a signal pair
          // Then we do not modify either point
          if((vertex2PairsTargetDiagram[pair.birth.id].size() >= 1)
             || (vertex2PairsTargetDiagram[pair.death.id].size() >= 1)) {
            continue;
          }
        }

        birthPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.birth.id));
        birthPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
        deathPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(pair.death.id));
        deathPairToDeleteTargetDiagram.push_back(
          (pair.birth.sfValue + pair.death.sfValue) / 2);
      }

      for(const auto &entry : matchingsBlockPairs) {
        // Delete pairs that have no equivalence
        if(entry.size() == 1) {

          if(!constraintAveraging_) {
            // If the point thresholdCurrentDiagram[entry[0]].birth.id is in a
            // signal pair AND If the point
            // thresholdCurrentDiagram[entry[0]].death.id is not in a signal
            // pair Then we only modify the
            // thresholdCurrentDiagram[entry[0]].death.id
            if((vertex2PairsTargetDiagram[thresholdCurrentDiagram[entry[0]]
                                            .birth.id]
                  .size()
                >= 1)
               && (vertex2PairsTargetDiagram[thresholdCurrentDiagram[entry[0]]
                                               .death.id]
                     .size()
                   == 0)) {
              deathPairToDeleteCurrentDiagram.push_back(static_cast<SimplexId>(
                thresholdCurrentDiagram[entry[0]].death.id));
              deathPairToDeleteTargetDiagram.push_back(
                (thresholdCurrentDiagram[entry[0]].birth.sfValue
                 + thresholdCurrentDiagram[entry[0]].death.sfValue)
                / 2);
              continue;
            }

            // If the point thresholdCurrentDiagram[entry[0]].death.id is in a
            // signal pair AND If the point
            // thresholdCurrentDiagram[entry[0]].birth.id is not in a signal
            // pair Then we only modify the
            // thresholdCurrentDiagram[entry[0]].birth.id
            if((vertex2PairsTargetDiagram[thresholdCurrentDiagram[entry[0]]
                                            .birth.id]
                  .size()
                == 0)
               && (vertex2PairsTargetDiagram[thresholdCurrentDiagram[entry[0]]
                                               .death.id]
                     .size()
                   >= 1)) {
              birthPairToDeleteCurrentDiagram.push_back(static_cast<SimplexId>(
                thresholdCurrentDiagram[entry[0]].birth.id));
              birthPairToDeleteTargetDiagram.push_back(
                (thresholdCurrentDiagram[entry[0]].birth.sfValue
                 + thresholdCurrentDiagram[entry[0]].death.sfValue)
                / 2);
              continue;
            }

            // If the point thresholdCurrentDiagram[entry[0]].birth.id is in a
            // signal pair AND If the point
            // thresholdCurrentDiagram[entry[0]].death.id is in a signal pair
            // Then we do not modify either point
            if((vertex2PairsTargetDiagram[thresholdCurrentDiagram[entry[0]]
                                            .birth.id]
                  .size()
                >= 1)
               || (vertex2PairsTargetDiagram[thresholdCurrentDiagram[entry[0]]
                                               .death.id]
                     .size()
                   >= 1)) {
              continue;
            }
          }

          birthPairToDeleteCurrentDiagram.push_back(
            static_cast<SimplexId>(thresholdCurrentDiagram[entry[0]].birth.id));
          birthPairToDeleteTargetDiagram.push_back(
            (thresholdCurrentDiagram[entry[0]].birth.sfValue
             + thresholdCurrentDiagram[entry[0]].death.sfValue)
            / 2);
          deathPairToDeleteCurrentDiagram.push_back(
            static_cast<SimplexId>(thresholdCurrentDiagram[entry[0]].death.id));
          deathPairToDeleteTargetDiagram.push_back(
            (thresholdCurrentDiagram[entry[0]].birth.sfValue
             + thresholdCurrentDiagram[entry[0]].death.sfValue)
            / 2);
          continue;
        } else if(entry.empty())
          continue;

        SimplexId valueBirthPairToChangeCurrentDiagram
          = static_cast<SimplexId>(thresholdCurrentDiagram[entry[0]].birth.id);
        SimplexId valueDeathPairToChangeCurrentDiagram
          = static_cast<SimplexId>(thresholdCurrentDiagram[entry[0]].death.id);

        double valueBirthPairToChangeTargetDiagram
          = thresholdConstraintDiagram[entry[1]].birth.sfValue;
        double valueDeathPairToChangeTargetDiagram
          = thresholdConstraintDiagram[entry[1]].death.sfValue;

        pair2MatchedPair[pairIndiceLocal2Global[entry[1]]][0]
          = thresholdCurrentDiagram[entry[0]].birth.id;
        pair2MatchedPair[pairIndiceLocal2Global[entry[1]]][1]
          = thresholdCurrentDiagram[entry[0]].death.id;

        pairChangeMatchingPair[pairIndiceLocal2Global[entry[1]]] = 1;

        birthPairToChangeCurrentDiagram.push_back(
          valueBirthPairToChangeCurrentDiagram);
        birthPairToChangeTargetDiagram.push_back(
          valueBirthPairToChangeTargetDiagram);
        deathPairToChangeCurrentDiagram.push_back(
          valueDeathPairToChangeCurrentDiagram);
        deathPairToChangeTargetDiagram.push_back(
          valueDeathPairToChangeTargetDiagram);
      }
    }
  }
  //=====================================//
  //            Basic Matching          //
  //=====================================//
  else {
    this->printMsg(
      "Get Indices | Compute Wasserstein distance: ", debug::Priority::DETAIL);

    if(epoch == 0) {
      for(SimplexId i = 0; i < (SimplexId)diagramOutput.size(); i++) {
        auto &pair = diagramOutput[i];
        currentVertex2PairsCurrentDiagram[pair.birth.id].push_back(i);
        currentVertex2PairsCurrentDiagram[pair.death.id].push_back(i);
      }
    } else {
      std::vector<std::vector<SimplexId>> newVertex2PairsCurrentDiagram(
        vertexNumber_, std::vector<SimplexId>());

      for(SimplexId i = 0; i < (SimplexId)diagramOutput.size(); i++) {
        auto &pair = diagramOutput[i];
        newVertex2PairsCurrentDiagram[pair.birth.id].push_back(i);
        newVertex2PairsCurrentDiagram[pair.death.id].push_back(i);
      }

      currentVertex2PairsCurrentDiagram = newVertex2PairsCurrentDiagram;
    }

    std::vector<std::vector<SimplexId>> vertex2PairsCurrentDiagram(
      vertexNumber_, std::vector<SimplexId>());
    for(SimplexId i = 0; i < (SimplexId)diagramOutput.size(); i++) {
      auto &pair = diagramOutput[i];
      vertex2PairsCurrentDiagram[pair.birth.id].push_back(i);
      vertex2PairsCurrentDiagram[pair.death.id].push_back(i);
      vertexInHowManyPairs[pair.birth.id]++;
      vertexInHowManyPairs[pair.death.id]++;
    }

    std::vector<std::vector<SimplexId>> vertex2PairsTargetDiagram(
      vertexNumber_, std::vector<SimplexId>());
    for(SimplexId i = 0; i < (SimplexId)constraintDiagram.size(); i++) {
      auto &pair = constraintDiagram[i];
      vertex2PairsTargetDiagram[pair.birth.id].push_back(i);
      vertex2PairsTargetDiagram[pair.death.id].push_back(i);
    }

    //=========================================
    //     Compute wasserstein distance
    //=========================================
    ttk::Timer timePersistenceDiagramClustering;

    ttk::PersistenceDiagramClustering persistenceDiagramClustering;
    PersistenceDiagramBarycenter pdBarycenter{};
    std::vector<ttk::DiagramType> intermediateDiagrams{
      constraintDiagram, diagramOutput};
    std::vector<ttk::DiagramType> centroids;
    std::vector<std::vector<std::vector<ttk::MatchingType>>> allMatchings;

    if(pdcMethod_ == 0) {
      persistenceDiagramClustering.setDebugLevel(debugLevel_);
      persistenceDiagramClustering.setThreadNumber(threadNumber_);
      // SetForceUseOfAlgorithm ==> Force the progressive approch if 2 inputs
      persistenceDiagramClustering.setForceUseOfAlgorithm(false);
      // setDeterministic ==> Deterministic algorithm
      persistenceDiagramClustering.setDeterministic(true);
      // setUseProgressive ==> Compute Progressive Barycenter
      persistenceDiagramClustering.setUseProgressive(true);
      // setUseInterruptible ==> Interruptible algorithm
      // persistenceDiagramClustering.setUseInterruptible(true);
      persistenceDiagramClustering.setUseInterruptible(false);
      // // setTimeLimit ==> Maximal computation time (s)
      persistenceDiagramClustering.setTimeLimit(0.01);
      // setUseAdditionalPrecision ==> Force minimum precision on matchings
      persistenceDiagramClustering.setUseAdditionalPrecision(true);
      // setDeltaLim ==> Minimal relative precision
      persistenceDiagramClustering.setDeltaLim(0.00000001);
      // setUseAccelerated ==> Use Accelerated KMeans
      persistenceDiagramClustering.setUseAccelerated(false);
      // setUseKmeansppInit ==> KMeanspp Initialization
      persistenceDiagramClustering.setUseKmeansppInit(false);

      std::vector<int> clusterIds = persistenceDiagramClustering.execute(
        intermediateDiagrams, centroids, allMatchings);
    } else {
      centroids.resize(1);
      const auto wassersteinMetric = std::to_string(2);
      pdBarycenter.setWasserstein(wassersteinMetric);
      pdBarycenter.setMethod(2);
      pdBarycenter.setNumberOfInputs(2);
      pdBarycenter.setDeterministic(true);
      pdBarycenter.setUseProgressive(true);
      pdBarycenter.setDebugLevel(debugLevel_);
      pdBarycenter.setThreadNumber(threadNumber_);
      pdBarycenter.setAlpha(1);
      pdBarycenter.setLambda(1);
      pdBarycenter.execute(intermediateDiagrams, centroids[0], allMatchings);
    }

    this->printMsg(
      "Get Indices | Persistence Diagram Clustering Time: "
        + std::to_string(timePersistenceDiagramClustering.getElapsedTime()),
      debug::Priority::DETAIL);

    //=========================================
    //             Find matched pairs
    //=========================================

    std::vector<std::vector<SimplexId>> allPairsSelected{};
    std::vector<std::vector<std::vector<double>>> matchingsBlock(
      centroids[0].size());
    std::vector<std::vector<ttk::PersistencePair>> matchingsBlockPairs(
      centroids[0].size());

    for(auto i = 1; i >= 0; --i) {
      std::vector<ttk::MatchingType> &matching = allMatchings[0][i];

      const auto &diag{intermediateDiagrams[i]};

      for(SimplexId j = 0; j < (SimplexId)matching.size(); j++) {

        const auto &m{matching[j]};
        const auto &bidderId{std::get<0>(m)};
        const auto &goodId{std::get<1>(m)};

        if((goodId == -1) | (bidderId == -1))
          continue;

        if(diag[bidderId].persistence() != 0) {
          matchingsBlock[goodId].push_back(
            {static_cast<double>(diag[bidderId].birth.id),
             static_cast<double>(diag[bidderId].death.id),
             diag[bidderId].persistence()});
          if(i == 1) {
            matchingsBlockPairs[goodId].push_back(diag[bidderId]);
          } else if(matchingsBlockPairs[goodId].size() > 0) {
            matchingsBlockPairs[goodId].push_back(diag[bidderId]);
          }
          allPairsSelected.push_back(
            {diag[bidderId].birth.id, diag[bidderId].death.id});
        }
      }
    }

    std::vector<ttk::PersistencePair> pairsToErase{};

    std::map<std::vector<SimplexId>, SimplexId> currentToTarget;
    for(auto &pair : allPairsSelected) {
      currentToTarget[{pair[0], pair[1]}] = 1;
    }

    for(auto &pair : intermediateDiagrams[1]) {
      if(pair.isFinite != 0) {
        if(!(currentToTarget.count({pair.birth.id, pair.death.id}) > 0)) {
          pairsToErase.push_back(pair);
        }
      }
    }

    for(auto &pair : pairsToErase) {

      if(!constraintAveraging_) {

        // If the point pair.birth.id is in a signal pair
        // AND If the point pair.death.id is not in a signal pair
        // Then we only modify the pair.death.id
        if((vertex2PairsTargetDiagram[pair.birth.id].size() >= 1)
           && (vertex2PairsTargetDiagram[pair.death.id].size() == 0)) {
          deathPairToDeleteCurrentDiagram.push_back(
            static_cast<SimplexId>(pair.death.id));
          deathPairToDeleteTargetDiagram.push_back(
            (pair.birth.sfValue + pair.death.sfValue) / 2);
          continue;
        }

        // If the point pair.death.id is in a signal pair
        // AND If the point pair.birth.id is not in a signal pair
        // Then we only modify the pair.birth.id
        if((vertex2PairsTargetDiagram[pair.birth.id].size() == 0)
           && (vertex2PairsTargetDiagram[pair.death.id].size() >= 1)) {
          birthPairToDeleteCurrentDiagram.push_back(
            static_cast<SimplexId>(pair.birth.id));
          birthPairToDeleteTargetDiagram.push_back(
            (pair.birth.sfValue + pair.death.sfValue) / 2);
          continue;
        }

        // If the point pair.birth.id is in a signal pair
        // AND If the point pair.death.id is in a signal pair
        // Then we do not modify either point
        if((vertex2PairsTargetDiagram[pair.birth.id].size() >= 1)
           || (vertex2PairsTargetDiagram[pair.death.id].size() >= 1)) {
          continue;
        }
      }

      birthPairToDeleteCurrentDiagram.push_back(
        static_cast<SimplexId>(pair.birth.id));
      birthPairToDeleteTargetDiagram.push_back(
        (pair.birth.sfValue + pair.death.sfValue) / 2);
      deathPairToDeleteCurrentDiagram.push_back(
        static_cast<SimplexId>(pair.death.id));
      deathPairToDeleteTargetDiagram.push_back(
        (pair.birth.sfValue + pair.death.sfValue) / 2);
    }

    for(const auto &entry : matchingsBlockPairs) {
      // Delete pairs that have no equivalence
      if(entry.size() == 1) {
        birthPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(entry[0].birth.id));
        birthPairToDeleteTargetDiagram.push_back(
          (entry[0].birth.sfValue + entry[0].death.sfValue) / 2);
        deathPairToDeleteCurrentDiagram.push_back(
          static_cast<SimplexId>(entry[0].death.id));
        deathPairToDeleteTargetDiagram.push_back(
          (entry[0].birth.sfValue + entry[0].death.sfValue) / 2);
        continue;
      } else if(entry.empty())
        continue;

      SimplexId valueBirthPairToChangeCurrentDiagram
        = static_cast<SimplexId>(entry[0].birth.id);
      SimplexId valueDeathPairToChangeCurrentDiagram
        = static_cast<SimplexId>(entry[0].death.id);

      double valueBirthPairToChangeTargetDiagram = entry[1].birth.sfValue;
      double valueDeathPairToChangeTargetDiagram = entry[1].death.sfValue;

      birthPairToChangeCurrentDiagram.push_back(
        valueBirthPairToChangeCurrentDiagram);
      birthPairToChangeTargetDiagram.push_back(
        valueBirthPairToChangeTargetDiagram);
      deathPairToChangeCurrentDiagram.push_back(
        valueDeathPairToChangeCurrentDiagram);
      deathPairToChangeTargetDiagram.push_back(
        valueDeathPairToChangeTargetDiagram);
    }
  }
}

/*
  This function allows you to copy the values of a pytorch tensor
  to a vector in an optimized way.
*/
#ifdef TTK_ENABLE_TORCH
int ttk::TopologicalOptimization::tensorToVectorFast(
  const torch::Tensor &tensor, std::vector<double> &result) const {
  TORCH_CHECK(
    tensor.dtype() == torch::kDouble, "The tensor must be of double type");
  const double *dataPtr = tensor.data_ptr<double>();
  result.assign(dataPtr, dataPtr + tensor.numel());

  return 0;
}
#endif

template <typename dataType, typename triangulationType>
int ttk::TopologicalOptimization::execute(
  const dataType *const inputScalars,
  dataType *const outputScalars,
  SimplexId *const inputOffsets,
  triangulationType *triangulation,
  const ttk::DiagramType &constraintDiagram) const {

  Timer t;
  double stoppingCondition = 0;
  bool enableTorch = true;

  if(methodOptimization_ == 1) {
#ifndef TTK_ENABLE_TORCH
    this->printWrn("Adam unavailable (Torch not found).");
    this->printWrn("Using direct gradient descent.");
    enableTorch = false;
#endif
  }

  //=======================
  //    Copy input data
  //=======================
  std::vector<double> dataVector(vertexNumber_);
  SimplexId *inputOffsetsCopie = inputOffsets;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId k = 0; k < vertexNumber_; ++k) {
    outputScalars[k] = inputScalars[k];
    dataVector[k] = inputScalars[k];
    if(std::isnan((double)outputScalars[k]))
      outputScalars[k] = 0;
  }

  //===============================
  //        Normalize the data
  //===============================

  dataType minVal = *std::min_element(dataVector.begin(), dataVector.end());
  dataType maxVal = *std::max_element(dataVector.begin(), dataVector.end());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < dataVector.size(); ++i) {
    dataVector[i] = (dataVector[i] - minVal) / (maxVal - minVal);
  }

  ttk::DiagramType normalizedConstraintDiagram(constraintDiagram.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)constraintDiagram.size(); i++) {
    auto pair = constraintDiagram[i];
    pair.birth.sfValue = (pair.birth.sfValue - minVal) / (maxVal - minVal);
    pair.death.sfValue = (pair.death.sfValue - minVal) / (maxVal - minVal);
    normalizedConstraintDiagram[i] = pair;
  }

  std::vector<double> losses;
  std::vector<double> inputScalarsX(vertexNumber_);

  //========================================
  //          Direct gradient descent
  //========================================
  if((methodOptimization_ == 0) || !(enableTorch)) {
    std::vector<SimplexId> listAllIndicesToChangeSmoothing(vertexNumber_, 0);
    std::vector<std::vector<SimplexId>> pair2MatchedPair(
      constraintDiagram.size(), std::vector<SimplexId>(2));
    std::vector<SimplexId> pairChangeMatchingPair(constraintDiagram.size(), -1);
    std::vector<std::vector<SimplexId>> pair2Delete(
      vertexNumber_, std::vector<SimplexId>());
    std::vector<std::vector<SimplexId>> currentVertex2PairsCurrentDiagram(
      vertexNumber_, std::vector<SimplexId>());

    for(int it = 0; it < epochNumber_; it++) {

      if(it % printFrequency_ == 0) {
        debugLevel_ = 3;
      } else {
        debugLevel_ = 0;
      }

      this->printMsg("DirectGradientDescent - iteration #" + std::to_string(it),
                     debug::Priority::PERFORMANCE);

      // pairs to change
      std::vector<SimplexId> birthPairToChangeCurrentDiagram{};
      std::vector<double> birthPairToChangeTargetDiagram{};
      std::vector<SimplexId> deathPairToChangeCurrentDiagram{};
      std::vector<double> deathPairToChangeTargetDiagram{};

      // pairs to delete
      std::vector<SimplexId> birthPairToDeleteCurrentDiagram{};
      std::vector<double> birthPairToDeleteTargetDiagram{};
      std::vector<SimplexId> deathPairToDeleteCurrentDiagram{};
      std::vector<double> deathPairToDeleteTargetDiagram{};

      std::vector<int> vertexInHowManyPairs(vertexNumber_, 0);

      getIndices(
        triangulation, inputOffsetsCopie, dataVector.data(),
        normalizedConstraintDiagram, it, listAllIndicesToChangeSmoothing,
        pair2MatchedPair, pair2Delete, pairChangeMatchingPair,
        birthPairToDeleteCurrentDiagram, birthPairToDeleteTargetDiagram,
        deathPairToDeleteCurrentDiagram, deathPairToDeleteTargetDiagram,
        birthPairToChangeCurrentDiagram, birthPairToChangeTargetDiagram,
        deathPairToChangeCurrentDiagram, deathPairToChangeTargetDiagram,
        currentVertex2PairsCurrentDiagram, vertexInHowManyPairs);
      std::fill(listAllIndicesToChangeSmoothing.begin(),
                listAllIndicesToChangeSmoothing.end(), 0);

      //==========================================================================
      //    Retrieve the indices for the pairs that we want to send diagonally
      //==========================================================================
      double lossDeletePairs = 0;

      std::vector<SimplexId> &indexBirthPairToDelete
        = birthPairToDeleteCurrentDiagram;
      std::vector<double> &targetValueBirthPairToDelete
        = birthPairToDeleteTargetDiagram;
      std::vector<SimplexId> &indexDeathPairToDelete
        = deathPairToDeleteCurrentDiagram;
      std::vector<double> &targetValueDeathPairToDelete
        = deathPairToDeleteTargetDiagram;

      this->printMsg("DirectGradientDescent - Number of pairs to delete: "
                       + std::to_string(indexBirthPairToDelete.size()),
                     debug::Priority::DETAIL);

      std::vector<int> vertexInCellMultiple(vertexNumber_, -1);
      std::vector<std::vector<double>> vertexToTargetValue(
        vertexNumber_, std::vector<double>());

      if(indexBirthPairToDelete.size() == indexDeathPairToDelete.size()) {
        for(size_t i = 0; i < indexBirthPairToDelete.size(); i++) {
          lossDeletePairs += std::pow(dataVector[indexBirthPairToDelete[i]]
                                        - targetValueBirthPairToDelete[i],
                                      2)
                             + std::pow(dataVector[indexDeathPairToDelete[i]]
                                          - targetValueDeathPairToDelete[i],
                                        2);
          SimplexId indexMax = indexBirthPairToDelete[i];
          SimplexId indexSelle = indexDeathPairToDelete[i];

          if(!(finePairManagement_ == 2) && !(finePairManagement_ == 1)) {
            if(constraintAveraging_) {
              if(vertexInHowManyPairs[indexMax] == 1) {
                dataVector[indexMax]
                  = dataVector[indexMax]
                    - alpha_ * 2
                        * (dataVector[indexMax]
                           - targetValueBirthPairToDelete[i]);
                listAllIndicesToChangeSmoothing[indexMax] = 1;
              } else {
                vertexInCellMultiple[indexMax] = 1;
                vertexToTargetValue[indexMax].push_back(
                  targetValueBirthPairToDelete[i]);
              }

              if(vertexInHowManyPairs[indexSelle] == 1) {
                dataVector[indexSelle]
                  = dataVector[indexSelle]
                    - alpha_ * 2
                        * (dataVector[indexSelle]
                           - targetValueDeathPairToDelete[i]);
                listAllIndicesToChangeSmoothing[indexSelle] = 1;
              } else {
                vertexInCellMultiple[indexSelle] = 1;
                vertexToTargetValue[indexSelle].push_back(
                  targetValueDeathPairToDelete[i]);
              }
            } else {
              dataVector[indexMax] = dataVector[indexMax]
                                     - alpha_ * 2
                                         * (dataVector[indexMax]
                                            - targetValueBirthPairToDelete[i]);
              dataVector[indexSelle]
                = dataVector[indexSelle]
                  - alpha_ * 2
                      * (dataVector[indexSelle]
                         - targetValueDeathPairToDelete[i]);
              listAllIndicesToChangeSmoothing[indexMax] = 1;
              listAllIndicesToChangeSmoothing[indexSelle] = 1;
            }
          } else if(finePairManagement_ == 1) {
            if(constraintAveraging_) {
              if(vertexInHowManyPairs[indexSelle] == 1) {
                dataVector[indexSelle]
                  = dataVector[indexSelle]
                    - alpha_ * 2
                        * (dataVector[indexSelle]
                           - targetValueDeathPairToDelete[i]);
                listAllIndicesToChangeSmoothing[indexSelle] = 1;
              } else {
                vertexInCellMultiple[indexSelle] = 1;
                vertexToTargetValue[indexSelle].push_back(
                  targetValueDeathPairToDelete[i]);
              }
            } else {
              dataVector[indexSelle]
                = dataVector[indexSelle]
                  - alpha_ * 2
                      * (dataVector[indexSelle]
                         - targetValueDeathPairToDelete[i]);
              listAllIndicesToChangeSmoothing[indexSelle] = 1;
            }
          } else if(finePairManagement_ == 2) {
            if(constraintAveraging_) {
              if(vertexInHowManyPairs[indexMax] == 1) {
                dataVector[indexMax]
                  = dataVector[indexMax]
                    - alpha_ * 2
                        * (dataVector[indexMax]
                           - targetValueBirthPairToDelete[i]);
                listAllIndicesToChangeSmoothing[indexMax] = 1;
              } else {
                vertexInCellMultiple[indexMax] = 1;
                vertexToTargetValue[indexMax].push_back(
                  targetValueBirthPairToDelete[i]);
              }
            } else {
              dataVector[indexMax] = dataVector[indexMax]
                                     - alpha_ * 2
                                         * (dataVector[indexMax]
                                            - targetValueBirthPairToDelete[i]);
              listAllIndicesToChangeSmoothing[indexMax] = 1;
            }
          }
        }
      } else {
        for(size_t i = 0; i < indexBirthPairToDelete.size(); i++) {
          lossDeletePairs += std::pow(dataVector[indexBirthPairToDelete[i]]
                                        - targetValueBirthPairToDelete[i],
                                      2);
          SimplexId indexMax = indexBirthPairToDelete[i];

          if(!(finePairManagement_ == 1)) {
            if(constraintAveraging_) {
              if(vertexInHowManyPairs[indexMax] == 1) {
                dataVector[indexMax]
                  = dataVector[indexMax]
                    - alpha_ * 2
                        * (dataVector[indexMax]
                           - targetValueBirthPairToDelete[i]);
                listAllIndicesToChangeSmoothing[indexMax] = 1;
              } else {
                vertexInCellMultiple[indexMax] = 1;
                vertexToTargetValue[indexMax].push_back(
                  targetValueBirthPairToDelete[i]);
              }
            } else {
              dataVector[indexMax] = dataVector[indexMax]
                                     - alpha_ * 2
                                         * (dataVector[indexMax]
                                            - targetValueBirthPairToDelete[i]);
              listAllIndicesToChangeSmoothing[indexMax] = 1;
            }
          } else { // finePairManagement_ == 1
            continue;
          }
        }

        for(size_t i = 0; i < indexDeathPairToDelete.size(); i++) {
          lossDeletePairs += std::pow(dataVector[indexDeathPairToDelete[i]]
                                        - targetValueDeathPairToDelete[i],
                                      2);
          SimplexId indexSelle = indexDeathPairToDelete[i];

          if(!(finePairManagement_ == 2)) {
            if(constraintAveraging_) {
              if(vertexInHowManyPairs[indexSelle] == 1) {
                dataVector[indexSelle]
                  = dataVector[indexSelle]
                    - alpha_ * 2
                        * (dataVector[indexSelle]
                           - targetValueDeathPairToDelete[i]);
                listAllIndicesToChangeSmoothing[indexSelle] = 1;
              } else {
                vertexInCellMultiple[indexSelle] = 1;
                vertexToTargetValue[indexSelle].push_back(
                  targetValueDeathPairToDelete[i]);
              }
            } else {
              dataVector[indexSelle]
                = dataVector[indexSelle]
                  - alpha_ * 2
                      * (dataVector[indexSelle]
                         - targetValueDeathPairToDelete[i]);
              listAllIndicesToChangeSmoothing[indexSelle] = 1;
            }
          } else { // finePairManagement_ == 2
            continue;
          }
        }
      }

      this->printMsg("DirectGradientDescent - Loss Delete Pairs: "
                       + std::to_string(lossDeletePairs),
                     debug::Priority::PERFORMANCE);
      //==========================================================================
      //      Retrieve the indices for the pairs that we want to change
      //==========================================================================
      double lossChangePairs = 0;

      std::vector<SimplexId> &indexBirthPairToChange
        = birthPairToChangeCurrentDiagram;
      std::vector<double> &targetValueBirthPairToChange
        = birthPairToChangeTargetDiagram;
      std::vector<SimplexId> &indexDeathPairToChange
        = deathPairToChangeCurrentDiagram;
      std::vector<double> &targetValueDeathPairToChange
        = deathPairToChangeTargetDiagram;

      for(size_t i = 0; i < indexBirthPairToChange.size(); i++) {
        lossChangePairs += std::pow(dataVector[indexBirthPairToChange[i]]
                                      - targetValueBirthPairToChange[i],
                                    2)
                           + std::pow(dataVector[indexDeathPairToChange[i]]
                                        - targetValueDeathPairToChange[i],
                                      2);

        SimplexId indexMax = indexBirthPairToChange[i];
        SimplexId indexSelle = indexDeathPairToChange[i];

        if(constraintAveraging_) {
          if(vertexInHowManyPairs[indexMax] == 1) {
            dataVector[indexMax]
              = dataVector[indexMax]
                - alpha_ * 2
                    * (dataVector[indexMax] - targetValueBirthPairToChange[i]);
            listAllIndicesToChangeSmoothing[indexMax] = 1;
          } else {
            vertexInCellMultiple[indexMax] = 1;
            vertexToTargetValue[indexMax].push_back(
              targetValueBirthPairToChange[i]);
          }

          if(vertexInHowManyPairs[indexSelle] == 1) {
            dataVector[indexSelle] = dataVector[indexSelle]
                                     - alpha_ * 2
                                         * (dataVector[indexSelle]
                                            - targetValueDeathPairToChange[i]);
            listAllIndicesToChangeSmoothing[indexSelle] = 1;
          } else {
            vertexInCellMultiple[indexSelle] = 1;
            vertexToTargetValue[indexSelle].push_back(
              targetValueDeathPairToChange[i]);
          }
        } else {
          dataVector[indexMax]
            = dataVector[indexMax]
              - alpha_ * 2
                  * (dataVector[indexMax] - targetValueBirthPairToChange[i]);
          dataVector[indexSelle]
            = dataVector[indexSelle]
              - alpha_ * 2
                  * (dataVector[indexSelle] - targetValueDeathPairToChange[i]);
          listAllIndicesToChangeSmoothing[indexMax] = 1;
          listAllIndicesToChangeSmoothing[indexSelle] = 1;
        }
      }

      this->printMsg("DirectGradientDescent - Loss Change Pairs: "
                       + std::to_string(lossChangePairs),
                     debug::Priority::PERFORMANCE);

      if(constraintAveraging_) {
        for(SimplexId i = 0; i < (SimplexId)vertexInCellMultiple.size(); i++) {
          double averageTargetValue = 0;

          if(vertexInCellMultiple[i] == 1) {
            for(auto targetValue : vertexToTargetValue[i]) {
              averageTargetValue += targetValue;
            }
            averageTargetValue
              = averageTargetValue / (int)vertexToTargetValue[i].size();

            dataVector[i] = dataVector[i]
                            - alpha_ * 2 * (dataVector[i] - averageTargetValue);
            listAllIndicesToChangeSmoothing[i] = 1;
          }
        }
      }

      //==================================
      //          Stop Condition
      //==================================

      if(it == 0) {
        stoppingCondition
          = coefStopCondition_ * (lossDeletePairs + lossChangePairs);
      }

      if(((lossDeletePairs + lossChangePairs) <= stoppingCondition))
        break;
    }

//========================================================
//              De-normalize data &  Update output data
//========================================================
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId k = 0; k < vertexNumber_; ++k) {
      outputScalars[k] = dataVector[k] * (maxVal - minVal) + minVal;
    }
  }

//=======================================
//           Adam Optimization
//=======================================
#ifdef TTK_ENABLE_TORCH
  else if(methodOptimization_ == 1) {
    //=====================================================
    //          Initialization of model parameters
    //=====================================================
    torch::Tensor F
      = torch::from_blob(dataVector.data(), {SimplexId(dataVector.size())},
                         torch::dtype(torch::kDouble))
          .to(torch::kDouble);
    PersistenceGradientDescent model(F);

    torch::optim::Adam optimizer(model.parameters(), learningRate_);

    //=======================================
    //            Optimization
    //=======================================

    std::vector<std::vector<SimplexId>> pair2MatchedPair(
      constraintDiagram.size(), std::vector<SimplexId>(2));
    std::vector<SimplexId> pairChangeMatchingPair(constraintDiagram.size(), -1);
    std::vector<SimplexId> listAllIndicesToChange(vertexNumber_, 0);
    std::vector<std::vector<SimplexId>> pair2Delete(
      vertexNumber_, std::vector<SimplexId>());
    std::vector<std::vector<SimplexId>> currentVertex2PairsCurrentDiagram(
      vertexNumber_, std::vector<SimplexId>());

    for(int i = 0; i < epochNumber_; i++) {

      if(i % printFrequency_ == 0) {
        debugLevel_ = 3;
      } else {
        debugLevel_ = 0;
      }

      this->printMsg(
        "Adam - epoch: " + std::to_string(i), debug::Priority::PERFORMANCE);

      ttk::Timer timeOneIteration;

      // Update the tensor with the new optimized values
      tensorToVectorFast(model.X.to(torch::kDouble), inputScalarsX);

      // pairs to change
      std::vector<SimplexId> birthPairToChangeCurrentDiagram{};
      std::vector<double> birthPairToChangeTargetDiagram{};
      std::vector<SimplexId> deathPairToChangeCurrentDiagram{};
      std::vector<double> deathPairToChangeTargetDiagram{};

      // pairs to delete
      std::vector<SimplexId> birthPairToDeleteCurrentDiagram{};
      std::vector<double> birthPairToDeleteTargetDiagram{};
      std::vector<SimplexId> deathPairToDeleteCurrentDiagram{};
      std::vector<double> deathPairToDeleteTargetDiagram{};

      std::vector<int> vertexInHowManyPairs(vertexNumber_, 0);

      // Retrieve the indices of the critical points that we must modify in
      // order to match our current diagram to our target diagram.
      getIndices(
        triangulation, inputOffsetsCopie, inputScalarsX.data(),
        normalizedConstraintDiagram, i, listAllIndicesToChange,
        pair2MatchedPair, pair2Delete, pairChangeMatchingPair,
        birthPairToDeleteCurrentDiagram, birthPairToDeleteTargetDiagram,
        deathPairToDeleteCurrentDiagram, deathPairToDeleteTargetDiagram,
        birthPairToChangeCurrentDiagram, birthPairToChangeTargetDiagram,
        deathPairToChangeCurrentDiagram, deathPairToChangeTargetDiagram,
        currentVertex2PairsCurrentDiagram, vertexInHowManyPairs);

      std::fill(
        listAllIndicesToChange.begin(), listAllIndicesToChange.end(), 0);
      //==========================================================================
      //    Retrieve the indices for the pairs that we want to send diagonally
      //==========================================================================

      torch::Tensor valueOfXDeleteBirth = torch::index_select(
        model.X, 0, torch::tensor(birthPairToDeleteCurrentDiagram));
      auto valueDeleteBirth = torch::from_blob(
        birthPairToDeleteTargetDiagram.data(),
        {static_cast<SimplexId>(birthPairToDeleteTargetDiagram.size())},
        torch::kDouble);
      torch::Tensor valueOfXDeleteDeath = torch::index_select(
        model.X, 0, torch::tensor(deathPairToDeleteCurrentDiagram));
      auto valueDeleteDeath = torch::from_blob(
        deathPairToDeleteTargetDiagram.data(),
        {static_cast<SimplexId>(deathPairToDeleteTargetDiagram.size())},
        torch::kDouble);

      torch::Tensor lossDeletePairs = torch::zeros({1}, torch::kDouble);
      if(!(finePairManagement_ == 2) && !(finePairManagement_ == 1)) {
        lossDeletePairs
          = torch::sum(torch::pow(valueOfXDeleteBirth - valueDeleteBirth, 2));
        lossDeletePairs
          = lossDeletePairs
            + torch::sum(torch::pow(valueOfXDeleteDeath - valueDeleteDeath, 2));
      } else if(finePairManagement_ == 1) {
        lossDeletePairs
          = torch::sum(torch::pow(valueOfXDeleteDeath - valueDeleteDeath, 2));
      } else if(finePairManagement_ == 2) {
        lossDeletePairs
          = torch::sum(torch::pow(valueOfXDeleteBirth - valueDeleteBirth, 2));
      }

      this->printMsg("Adam - Loss Delete Pairs: "
                       + std::to_string(lossDeletePairs.item<double>()),
                     debug::Priority::PERFORMANCE);

      //==========================================================================
      //      Retrieve the indices for the pairs that we want to change
      //==========================================================================

      torch::Tensor valueOfXChangeBirth = torch::index_select(
        model.X, 0, torch::tensor(birthPairToChangeCurrentDiagram));
      auto valueChangeBirth = torch::from_blob(
        birthPairToChangeTargetDiagram.data(),
        {static_cast<SimplexId>(birthPairToChangeTargetDiagram.size())},
        torch::kDouble);
      torch::Tensor valueOfXChangeDeath = torch::index_select(
        model.X, 0, torch::tensor(deathPairToChangeCurrentDiagram));
      auto valueChangeDeath = torch::from_blob(
        deathPairToChangeTargetDiagram.data(),
        {static_cast<SimplexId>(deathPairToChangeTargetDiagram.size())},
        torch::kDouble);

      auto lossChangePairs
        = torch::sum((torch::pow(valueOfXChangeBirth - valueChangeBirth, 2)
                      + torch::pow(valueOfXChangeDeath - valueChangeDeath, 2)));

      this->printMsg("Adam - Loss Change Pairs: "
                       + std::to_string(lossChangePairs.item<double>()),
                     debug::Priority::PERFORMANCE);

      //====================================
      //      Definition of final loss
      //====================================

      auto loss = lossDeletePairs + lossChangePairs;

      this->printMsg("Adam - Loss: " + std::to_string(loss.item<double>()),
                     debug::Priority::PERFORMANCE);

      //==========================================
      //            Back Propagation
      //==========================================

      losses.push_back(loss.item<double>());

      ttk::Timer timeBackPropagation;
      optimizer.zero_grad();
      loss.backward();
      optimizer.step();

      //==========================================
      //         Modified index checking
      //==========================================

      // On trouve les indices qui ont changé
      std::vector<double> NewinputScalarsX(vertexNumber_);
      tensorToVectorFast(model.X.to(torch::kDouble), NewinputScalarsX);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(SimplexId k = 0; k < vertexNumber_; ++k) {
        double diff = NewinputScalarsX[k] - inputScalarsX[k];
        if(diff != 0) {
          listAllIndicesToChange[k] = 1;
        }
      }

      //=======================================
      //              Stop condition
      //=======================================
      if(i == 0) {
        stoppingCondition = coefStopCondition_ * loss.item<double>();
      }

      if(loss.item<double>() < stoppingCondition)
        break;
    }

//============================================
//   De-normalize data &  Update output data
//============================================
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId k = 0; k < vertexNumber_; ++k) {
      outputScalars[k]
        = model.X[k].item().to<double>() * (maxVal - minVal) + minVal;
      if(std::isnan((double)outputScalars[k]))
        outputScalars[k] = 0;
    }
  }
#endif
  //========================================
  //            Information display
  //========================================
  debugLevel_ = 3;

  // Total execution time
  double time = t.getElapsedTime();

  // Number Pairs Constraint Diagram
  SimplexId numberPairsConstraintDiagram = (SimplexId)constraintDiagram.size();
  this->printMsg("Number of constrained pairs: "
                   + std::to_string(numberPairsConstraintDiagram),
                 debug::Priority::PERFORMANCE);

  this->printMsg("Stopping condition: " + std::to_string(stoppingCondition),
                 debug::Priority::PERFORMANCE);

  this->printMsg("Scalar field optimized", 1.0, time, this->threadNumber_);

  return 0;
}
