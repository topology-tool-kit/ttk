#include "PersistenceDiagramUtils.h"
#include <algorithm>
#include <cmath>
#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif

#include <PersistenceDiagramDictionary.h>

using namespace ttk;

void PersistenceDiagramDictionary::execute(
  std::vector<ttk::DiagramType> &intermediateDiagrams,
  const std::vector<ttk::DiagramType> &intermediateAtoms,
  std::vector<ttk::DiagramType> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  const int seed,
  const int numAtom,
  std::vector<double> &lossTab,
  std::vector<std::vector<double>> &allLosses,
  double percent) {

  if(!ProgApproach_) {
    // Regular approach
    for(size_t i = 0; i < intermediateDiagrams.size(); ++i) {
      if(sortedForTest_) {
        auto &diag = intermediateDiagrams[i];
        std::sort(diag.begin(), diag.end(),
                  [](const PersistencePair &t1, const PersistencePair &t2) {
                    return (t1.death.sfValue - t1.birth.sfValue)
                           > (t2.death.sfValue - t2.birth.sfValue);
                  });
      }
    }

    bool doCompression = false;
    if(CompressionMode_) {
      doCompression = true;
    }

    std::vector<BidderDiagram> trueBidderDiagramMin(
      intermediateDiagrams.size());
    std::vector<BidderDiagram> trueBidderDiagramSad(
      intermediateDiagrams.size());
    std::vector<BidderDiagram> trueBidderDiagramMax(
      intermediateDiagrams.size());

    std::vector<std::vector<double>> histoVectorWeights(
      intermediateDiagrams.size());
    std::vector<ttk::DiagramType> histoDictDiagrams(numAtom);
    this->maxLag2_ = 5;

    Timer tm_init{};
    bool preWeightOpt = false;
    this->printMsg("Regular approach:");
    initDictionary(dictDiagrams, intermediateDiagrams, intermediateAtoms,
                   numAtom, this->do_min_, this->do_sad_, this->do_max_, seed,
                   percent);
    this->printMsg("Initialization computed ", 1, tm_init.getElapsedTime(),
                   threadNumber_, debug::LineMode::NEW);

    controlAtomsSize(intermediateDiagrams, dictDiagrams);

    method(intermediateDiagrams, dictDiagrams, vectorWeights, numAtom, lossTab,
           allLosses, histoVectorWeights, histoDictDiagrams, preWeightOpt,
           percent, doCompression);
  } else {
    // Multi scale approach
    bool doCompression = false;
    for(size_t i = 0; i < intermediateDiagrams.size(); ++i) {
      auto &diag = intermediateDiagrams[i];
      std::sort(diag.begin(), diag.end(),
                [](const PersistencePair &t1, const PersistencePair &t2) {
                  return (t1.death.sfValue - t1.birth.sfValue)
                         > (t2.death.sfValue - t2.birth.sfValue);
                });
    }

    int start = 20;
    double stop = percent;
    std::vector<double> percentages;
    for(int value = start; value > stop; value -= 5) {
      percentages.push_back(static_cast<double>(value) / 100.);
    }
    percentages.push_back(static_cast<double>(percent) / 100.);
    std::vector<std::vector<double>> histoVectorWeights(
      intermediateDiagrams.size());
    std::vector<ttk::DiagramType> histoDictDiagrams(numAtom);
    std::vector<ttk::DiagramType> dataTemp(intermediateDiagrams.size());
    bool preWeightOpt = true;
    Timer tm_method{};
    for(size_t j = 0; j < 1; ++j) {
      double percentage = percentages[j];
      this->maxLag2_ = 0;
      this->printMsg("First step multi-scale approach:");
      for(size_t i = 0; i < intermediateDiagrams.size(); ++i) {
        auto &diag = intermediateDiagrams[i];
        auto &t = diag[0];
        double maxPers = getMaxPers(diag);
        dataTemp[i].push_back(t);
        for(size_t p = 1; p < diag.size(); ++p) {
          auto &t2 = diag[p];
          if(percentage * maxPers <= (t2.death.sfValue - t2.birth.sfValue)) {
            dataTemp[i].push_back(t2);
          } else {
            continue;
          }
        }
      }

      Timer tm_init{};
      initDictionary(dictDiagrams, dataTemp, intermediateAtoms, numAtom,
                     this->do_min_, this->do_sad_, this->do_max_, seed,
                     percent);
      this->printMsg("Initialization computed ", 1, tm_init.getElapsedTime(),
                     threadNumber_, debug::LineMode::NEW);

      method(dataTemp, dictDiagrams, vectorWeights, numAtom, lossTab, allLosses,
             histoVectorWeights, histoDictDiagrams, preWeightOpt, percent,
             doCompression);
    }

    int counter = 0;
    for(size_t j = 1; j < percentages.size(); ++j) {
      if(j == percentages.size() - static_cast<size_t>(1) && CompressionMode_) {
        doCompression = true;
      }
      double percentage = percentages[j];
      double previousPerc = percentages[j - 1];
      if(j < percentages.size() - 1) {
        this->maxLag2_ = 0;
      } else {
        this->maxLag2_ = 5;
      }
      if(j < percentages.size() - 1) {
        this->printMsg("New step multi-scale approach:");
      } else {
        this->printMsg("Final step multi-scale approach:");
      }

      for(size_t i = 0; i < intermediateDiagrams.size(); ++i) {
        auto &diag = intermediateDiagrams[i];
        double maxPers = getMaxPers(diag);
        for(size_t p = 0; p < diag.size(); ++p) {
          auto &t2 = diag[p];
          if(percentage * maxPers <= (t2.death.sfValue - t2.birth.sfValue)
             && (t2.death.sfValue - t2.birth.sfValue)
                  < previousPerc * maxPers) {
            dataTemp[i].push_back(t2);
            counter += 1;
          }
        }
      }
      if(counter == 0) {
        continue;
      }
      method(dataTemp, dictDiagrams, vectorWeights, numAtom, lossTab, allLosses,
             histoVectorWeights, histoDictDiagrams, preWeightOpt, percent,
             doCompression);
    }

    this->printMsg(
      "Total time", 1.0, tm_method.getElapsedTime(), threadNumber_);
  }
}

void PersistenceDiagramDictionary::method(
  const std::vector<ttk::DiagramType> &intermediateDiagrams,
  std::vector<ttk::DiagramType> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  const int numAtom,
  std::vector<double> &lossTab,
  std::vector<std::vector<double>> &allLosses,
  std::vector<std::vector<double>> &histoVectorWeights,
  std::vector<ttk::DiagramType> &histoDictDiagrams,
  bool preWeightOpt,
  double percent,
  bool doCompression) {

  Timer tm{};

  bool doOptimizeAtoms = false;
  bool doOptimizeWeights = false;

  if(OptimizeWeights_) {
    printMsg("Weight Optimization activated");
    doOptimizeWeights = true;
  } else {
    printWrn("Weight Optimization desactivated");
  }
  if(OptimizeAtoms_) {
    doOptimizeAtoms = true;
    printMsg("Atom Optimization activated");
  } else {
    printWrn("Atom Optimization desactivated");
  }

  const auto nDiags = intermediateDiagrams.size();

  // tracking the original indices
  std::vector<ttk::DiagramType> inputDiagramsMin(nDiags);
  std::vector<ttk::DiagramType> inputDiagramsSad(nDiags);
  std::vector<ttk::DiagramType> inputDiagramsMax(nDiags);

  std::vector<BidderDiagram> bidderDiagramsMin{};
  std::vector<BidderDiagram> bidderDiagramsSad{};
  std::vector<BidderDiagram> bidderDiagramsMax{};

  std::vector<std::vector<size_t>> originIndexDatasMin(nDiags);
  std::vector<std::vector<size_t>> originIndexDatasSad(nDiags);
  std::vector<std::vector<size_t>> originIndexDatasMax(nDiags);

  gettingBidderDiagrams(
    intermediateDiagrams, inputDiagramsMin, inputDiagramsSad, inputDiagramsMax,
    bidderDiagramsMin, bidderDiagramsSad, bidderDiagramsMax,
    originIndexDatasMin, originIndexDatasSad, originIndexDatasMax, true);

  std::vector<ttk::DiagramType> barycentersList(nDiags);
  std::vector<std::vector<std::vector<ttk::MatchingType>>> allMatchingsAtoms(
    nDiags);
  std::vector<std::vector<ttk::MatchingType>> matchingsDatasMin(nDiags);
  std::vector<std::vector<ttk::MatchingType>> matchingsDatasSad(nDiags);
  std::vector<std::vector<ttk::MatchingType>> matchingsDatasMax(nDiags);
  ConstrainedGradientDescent gradActor;
  double loss;
  int lag = 0;
  int lag2 = 0;
  int lag3 = 0;
  int lagLimit = 10;
  int minEpoch = 20;
  int maxEpoch = MaxEpoch_;
  bool cond = true;
  int epoch = 0;
  int nbEpochPrevious = lossTab.size();
  std::vector<size_t> initSizes(dictDiagrams.size());
  for(size_t j = 0; j < dictDiagrams.size(); ++j) {
    initSizes[j] = dictDiagrams[j].size();
  }

  double factEquiv = static_cast<double>(numAtom);
  double step = 1. / (2. * 2. * factEquiv);
  gradActor.setStep(factEquiv);

  // BUFFERS
  std::vector<std::vector<int>> bufferHistoAllEpochLife(dictDiagrams.size());
  std::vector<std::vector<bool>> bufferHistoAllBoolLife(dictDiagrams.size());
  std::vector<std::vector<bool>> bufferCheckUnderDiag(dictDiagrams.size());
  std::vector<std::vector<bool>> bufferCheckDiag(dictDiagrams.size());
  std::vector<std::vector<bool>> bufferCheckAboveGlobal(dictDiagrams.size());

  std::vector<std::vector<int>> histoAllEpochLife(dictDiagrams.size());
  std::vector<std::vector<bool>> histoAllBoolLife(dictDiagrams.size());
  std::vector<std::vector<bool>> checkUnderDiag(dictDiagrams.size());
  std::vector<std::vector<bool>> checkDiag(dictDiagrams.size());
  std::vector<std::vector<bool>> checkAboveGlobal(dictDiagrams.size());

  std::vector<double> allLossesAtEpoch(nDiags, 0.);
  std::vector<double> trueAllLossesAtEpoch(nDiags, 0.);

  while(epoch < maxEpoch && cond) {
    loss = 0.;
    Timer tm_it{};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    //////////////////////////////WEIGHTS///////////////////////////////////
    for(size_t i = 0; i < nDiags; ++i) {
      auto &barycenter = barycentersList[i];
      std::vector<double> &weight = vectorWeights[i];
      std::vector<std::vector<ttk::MatchingType>> &matchings
        = allMatchingsAtoms[i];
      computeWeightedBarycenter(
        dictDiagrams, weight, barycenter, matchings, *this, ProgBarycenter_);
    }
    this->printMsg(
      "Computed 1st Barycenters for epoch " + std::to_string(epoch),
      epoch / static_cast<double>(maxEpoch), tm_it.getElapsedTime(),
      threadNumber_, debug::LineMode::NEW, debug::Priority::DETAIL);

    //   "====================BARYCENTER FINISHED======================");

    // tracking the original indices
    std::vector<ttk::DiagramType> barycentersListMin(nDiags);
    std::vector<ttk::DiagramType> barycentersListSad(nDiags);
    std::vector<ttk::DiagramType> barycentersListMax(nDiags);

    std::vector<BidderDiagram> bidderBarycentersListMin{};
    std::vector<BidderDiagram> bidderBarycentersListSad{};
    std::vector<BidderDiagram> bidderBarycentersListMax{};

    std::vector<std::vector<size_t>> originIndexBarysMin(nDiags);
    std::vector<std::vector<size_t>> originIndexBarysSad(nDiags);
    std::vector<std::vector<size_t>> originIndexBarysMax(nDiags);

    computeAllDistances(
      barycentersList, nDiags, barycentersListMin, barycentersListSad,
      barycentersListMax, bidderBarycentersListMin, bidderBarycentersListSad,
      bidderBarycentersListMax, originIndexBarysMin, originIndexBarysSad,
      originIndexBarysMax, bidderDiagramsMin, bidderDiagramsMax,
      bidderDiagramsSad, matchingsDatasMin, matchingsDatasMax,
      matchingsDatasSad, allLossesAtEpoch, true);

    for(size_t p = 0; p < nDiags; ++p) {
      loss += allLossesAtEpoch[p];
    }

    for(size_t p = 0; p < nDiags; ++p) {
      if(!ProgApproach_) {
        allLosses[p].push_back(allLossesAtEpoch[p]);
      } else {
        allLosses[p].push_back(trueAllLossesAtEpoch[p]);
      }
    }

    lossTab.push_back(loss);

    printMsg(
      " Epoch " + std::to_string(epoch) + ", loss = " + std::to_string(loss), 1,
      threadNumber_, ttk::debug::LineMode::REPLACE);

    if(preWeightOpt && OptimizeAtoms_) {
      if(epoch < 10) {
        // Pre optimization of barycentric weights
        doOptimizeAtoms = false;
      } else {
        doOptimizeAtoms = true;
      }
    }

    if(loss < 1e-7) {
      cond = false;
      doOptimizeWeights = false;
      doOptimizeAtoms = false;
    }

    double mini
      = *std::min_element(lossTab.begin() + nbEpochPrevious, lossTab.end() - 1);
    if(loss <= mini) {
      for(size_t p = 0; p < dictDiagrams.size(); ++p) {
        const auto atom = dictDiagrams[p];
        histoDictDiagrams[p] = atom;

        const auto &histoEpochAtom = histoAllEpochLife[p];
        const auto &histoBoolAtom = histoAllBoolLife[p];
        const auto &boolUnderDiag = checkUnderDiag[p];
        const auto &boolDiag = checkDiag[p];
        const auto &boolAboveGlobal = checkAboveGlobal[p];

        bufferHistoAllEpochLife[p] = histoEpochAtom;
        bufferHistoAllBoolLife[p] = histoBoolAtom;
        bufferCheckUnderDiag[p] = boolUnderDiag;
        bufferCheckDiag[p] = boolDiag;
        bufferCheckAboveGlobal[p] = boolAboveGlobal;
      }
      for(size_t p = 0; p < nDiags; ++p) {
        const auto weights = vectorWeights[p];
        histoVectorWeights[p] = weights;
      }
      lag = 0;
      lag3 = 0;
    } else {
      if(epoch > minEpoch) {
        lag += 1;
      } else {
        lag = 0;
      }
    }

    if((epoch > minEpoch)
       && (lossTab[epoch + nbEpochPrevious]
             / lossTab[epoch + nbEpochPrevious - 1]
           > 0.99)) {
      if(lossTab[epoch + nbEpochPrevious]
         < lossTab[epoch + nbEpochPrevious - 1]) {
        if(lag2 == this->maxLag2_) {
          this->printMsg("Loss not decreasing enough");
          if(StopCondition_) {
            for(size_t p = 0; p < dictDiagrams.size(); ++p) {
              const auto atom = histoDictDiagrams[p];
              dictDiagrams[p] = atom;
            }
            for(size_t p = 0; p < nDiags; ++p) {
              const auto weights = histoVectorWeights[p];
              vectorWeights[p] = weights;
            }
            doOptimizeWeights = false;
            doOptimizeAtoms = false;
            cond = false;
          }
        } else {
          lag2 += 1;
        }
      }
    } else {
      lag2 = 0;
    }

    if(epoch > minEpoch && lag > lagLimit) {

      if(StopCondition_) {
        for(size_t p = 0; p < dictDiagrams.size(); ++p) {
          const auto atom = histoDictDiagrams[p];
          dictDiagrams[p] = atom;
        }
        for(size_t p = 0; p < nDiags; ++p) {
          const auto weights = histoVectorWeights[p];
          vectorWeights[p] = weights;
        }
        this->printMsg("Minimum not passed");
        doOptimizeWeights = false;
        doOptimizeAtoms = false;

        cond = false;
      }
    }

    if(cond && (epoch > 0) && (lossTab[epoch + nbEpochPrevious] > 2. * mini)) {
      lag3 += 1;
      lag2 = 0;

      if(lag3 > 2) {
        this->printMsg("Loss increasing too much, reducing step and recompute "
                       "Barycenters and matchings");
        for(size_t p = 0; p < dictDiagrams.size(); ++p) {
          const auto atom = histoDictDiagrams[p];
          dictDiagrams[p] = atom;

          const auto &bufferHistoEpochAtom = bufferHistoAllEpochLife[p];
          const auto &bufferHistoBoolAtom = bufferHistoAllBoolLife[p];
          const auto &bufferBoolUnderDiag = bufferCheckUnderDiag[p];
          const auto &bufferBoolDiag = bufferCheckDiag[p];
          const auto &bufferBoolAboveGlobal = bufferCheckAboveGlobal[p];

          histoAllEpochLife[p] = bufferHistoEpochAtom;
          histoAllBoolLife[p] = bufferHistoBoolAtom;
          checkUnderDiag[p] = bufferBoolUnderDiag;
          checkDiag[p] = bufferBoolDiag;
          checkAboveGlobal[p] = bufferBoolAboveGlobal;
        }

        for(size_t p = 0; p < nDiags; ++p) {
          const auto weights = histoVectorWeights[p];
          vectorWeights[p] = weights;
        }
        gradActor.reduceStep();
        step = step / 2.;
        lag3 = 0;

        barycentersList.clear();
        barycentersList.resize(nDiags);
        allMatchingsAtoms.clear();
        allMatchingsAtoms.resize(nDiags);

        barycentersListMin.clear();
        barycentersListSad.clear();
        barycentersListMax.clear();

        barycentersListMin.resize(nDiags);
        barycentersListSad.resize(nDiags);
        barycentersListMax.resize(nDiags);

        bidderBarycentersListMin.clear();
        bidderBarycentersListSad.clear();
        bidderBarycentersListMax.clear();

        originIndexBarysMin.clear();
        originIndexBarysSad.clear();
        originIndexBarysMax.clear();

        originIndexBarysMin.resize(nDiags);
        originIndexBarysSad.resize(nDiags);
        originIndexBarysMax.resize(nDiags);

        matchingsDatasMin.clear();
        matchingsDatasSad.clear();
        matchingsDatasMax.clear();

        matchingsDatasMin.resize(nDiags);
        matchingsDatasSad.resize(nDiags);
        matchingsDatasMax.resize(nDiags);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
       //////////////////////////////WEIGHTS///////////////////////////////////
        for(size_t i = 0; i < nDiags; ++i) {
          auto &barycenter = barycentersList[i];
          std::vector<double> &weight = vectorWeights[i];
          std::vector<std::vector<ttk::MatchingType>> &matchings
            = allMatchingsAtoms[i];
          computeWeightedBarycenter(dictDiagrams, weight, barycenter, matchings,
                                    *this, ProgBarycenter_);
        }

        computeAllDistances(
          barycentersList, nDiags, barycentersListMin, barycentersListSad,
          barycentersListMax, bidderBarycentersListMin,
          bidderBarycentersListSad, bidderBarycentersListMax,
          originIndexBarysMin, originIndexBarysSad, originIndexBarysMax,
          bidderDiagramsMin, bidderDiagramsMax, bidderDiagramsSad,
          matchingsDatasMin, matchingsDatasMax, matchingsDatasSad,
          allLossesAtEpoch, false);
      }
    }

    std::vector<std::vector<Matrix>> allHessianLists(nDiags);
    std::vector<std::vector<double>> gradWeightsList(nDiags);
    Timer tm_opt1{};

    // WEIGHT OPTIMIZATION
    if(doOptimizeWeights) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nDiags; ++i) {
        auto &gradWeights = gradWeightsList[i];
        const auto &matchingsAtoms = allMatchingsAtoms[i];
        const auto &Barycenter = barycentersList[i];
        const auto &Data = intermediateDiagrams[i];
        std::vector<Matrix> &hessianList = allHessianLists[i];
        const std::vector<ttk::MatchingType> &matchingsMin
          = matchingsDatasMin[i];
        const std::vector<ttk::MatchingType> &matchingsMax
          = matchingsDatasMax[i];
        const std::vector<ttk::MatchingType> &matchingsSad
          = matchingsDatasSad[i];
        const std::vector<size_t> &indexBaryMin = originIndexBarysMin[i];
        const std::vector<size_t> &indexBarySad = originIndexBarysSad[i];
        const std::vector<size_t> &indexBaryMax = originIndexBarysMax[i];
        const std::vector<size_t> &indexDataMin = originIndexDatasMin[i];
        const std::vector<size_t> &indexDataSad = originIndexDatasSad[i];
        const std::vector<size_t> &indexDataMax = originIndexDatasMax[i];
        std::vector<double> &weights = vectorWeights[i];
        computeGradientWeights(gradWeights, hessianList, dictDiagrams,
                               matchingsAtoms, Barycenter, Data, matchingsMin,
                               matchingsMax, matchingsSad, indexBaryMin,
                               indexBaryMax, indexBarySad, indexDataMin,
                               indexDataMax, indexDataSad, doOptimizeAtoms);
        gradActor.executeWeightsProjected(
          hessianList, weights, gradWeights, MaxEigenValue_);
      }

      this->printMsg("Computed 1st opt for epoch " + std::to_string(epoch),
                     epoch / static_cast<double>(maxEpoch),
                     tm_opt1.getElapsedTime(), threadNumber_,
                     debug::LineMode::NEW, debug::Priority::DETAIL);
    }

    for(size_t p = 0; p < nDiags; ++p) {
      allLossesAtEpoch[p] = 0.;
      trueAllLossesAtEpoch[p] = 0.;
    }
    barycentersList.clear();
    barycentersList.resize(nDiags);
    allMatchingsAtoms.clear();
    allMatchingsAtoms.resize(nDiags);
    ////////////////////////////////ATOM////////////////////////////////////////
    if(doOptimizeAtoms) {
      Timer tm_it2{};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nDiags; ++i) {
        auto &barycenter = barycentersList[i];
        std::vector<double> &weight = vectorWeights[i];
        std::vector<std::vector<ttk::MatchingType>> &matchings
          = allMatchingsAtoms[i];
        computeWeightedBarycenter(
          dictDiagrams, weight, barycenter, matchings, *this, ProgBarycenter_);
      }
      this->printMsg(
        "Computed 2nd Barycenters for epoch " + std::to_string(epoch),
        epoch / static_cast<double>(maxEpoch), tm_it2.getElapsedTime(),
        threadNumber_, debug::LineMode::NEW, debug::Priority::DETAIL);

      barycentersListMin.clear();
      barycentersListSad.clear();
      barycentersListMax.clear();

      barycentersListMin.resize(nDiags);
      barycentersListSad.resize(nDiags);
      barycentersListMax.resize(nDiags);

      bidderBarycentersListMin.clear();
      bidderBarycentersListSad.clear();
      bidderBarycentersListMax.clear();

      originIndexBarysMin.clear();
      originIndexBarysSad.clear();
      originIndexBarysMax.clear();

      originIndexBarysMin.resize(nDiags);
      originIndexBarysSad.resize(nDiags);
      originIndexBarysMax.resize(nDiags);

      matchingsDatasMin.clear();
      matchingsDatasSad.clear();
      matchingsDatasMax.clear();

      matchingsDatasMin.resize(nDiags);
      matchingsDatasSad.resize(nDiags);
      matchingsDatasMax.resize(nDiags);

      computeAllDistances(
        barycentersList, nDiags, barycentersListMin, barycentersListSad,
        barycentersListMax, bidderBarycentersListMin, bidderBarycentersListSad,
        bidderBarycentersListMax, originIndexBarysMin, originIndexBarysSad,
        originIndexBarysMax, bidderDiagramsMin, bidderDiagramsMax,
        bidderDiagramsSad, matchingsDatasMin, matchingsDatasMax,
        matchingsDatasSad, allLossesAtEpoch, false);

      std::vector<std::vector<std::vector<std::array<double, 2>>>>
        allPairToAddToGradList(nDiags);
      std::vector<ttk::DiagramType> allInfoToAdd(nDiags);
      std::vector<std::vector<Matrix>> gradsAtomsList(nDiags);
      std::vector<std::vector<int>> checkerAtomsList(nDiags);

      bool doDimReduct = false;
      if(DimReductMode_ && numAtom <= 3) {
        doDimReduct = true;
      } else {
        doDimReduct = false;
      }

      Timer tm_opt2{};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nDiags; ++i) {
        auto &pairToAddGradList = allPairToAddToGradList[i];
        auto &infoToAdd = allInfoToAdd[i];
        auto &gradsAtoms = gradsAtomsList[i];
        auto &checkerAtoms = checkerAtomsList[i];
        const auto &Barycenter = barycentersList[i];
        const auto &Data = intermediateDiagrams[i];
        const std::vector<ttk::MatchingType> &matchingsMin
          = matchingsDatasMin[i];
        const std::vector<ttk::MatchingType> &matchingsMax
          = matchingsDatasMax[i];
        const std::vector<ttk::MatchingType> &matchingsSad
          = matchingsDatasSad[i];
        const std::vector<size_t> &indexBaryMin = originIndexBarysMin[i];
        const std::vector<size_t> &indexBarySad = originIndexBarysSad[i];
        const std::vector<size_t> &indexBaryMax = originIndexBarysMax[i];
        const std::vector<size_t> &indexDataMin = originIndexDatasMin[i];
        const std::vector<size_t> &indexDataSad = originIndexDatasSad[i];
        const std::vector<size_t> &indexDataMax = originIndexDatasMax[i];
        const std::vector<double> &weights = vectorWeights[i];
        computeGradientAtoms(
          gradsAtoms, weights, Barycenter, Data, matchingsMin, matchingsMax,
          matchingsSad, indexBaryMin, indexBaryMax, indexBarySad, indexDataMin,
          indexDataMax, indexDataSad, checkerAtoms, pairToAddGradList,
          infoToAdd, doDimReduct);
      }

      std::vector<double> maxiDeath(numAtom);
      std::vector<double> minBirth(numAtom);
      for(int j = 0; j < numAtom; ++j) {
        auto &atom = dictDiagrams[j];
        auto &temp = atom[0];
        maxiDeath[j] = temp.death.sfValue;
        minBirth[j] = temp.birth.sfValue;
      }

      std::vector<std::vector<std::vector<int>>> allProjectionsList(nDiags);
      std::vector<ttk::DiagramType> allFeaturesToAdd(nDiags);
      std::vector<std::vector<std::array<double, 2>>> allProjLocations(nDiags);
      std::vector<std::vector<std::vector<double>>>
        allVectorForProjContributions(nDiags);
      std::vector<ttk::DiagramType> allTrueFeaturesToAdd(numAtom);
      std::vector<std::vector<std::vector<int>>> allTrueProj(numAtom);
      std::vector<std::vector<std::array<double, 2>>> allTrueProjLoc(numAtom);

      for(size_t i = 0; i < nDiags; ++i) {
        auto &pairToAddGradList = allPairToAddToGradList[i];
        auto &infoToAdd = allInfoToAdd[i];
        auto &projForDiag = allProjectionsList[i];
        auto &vectorForProjContrib = allVectorForProjContributions[i];
        auto &featuresToAdd = allFeaturesToAdd[i];
        auto &projLocations = allProjLocations[i];
        auto &gradsAtoms = gradsAtomsList[i];
        const auto &matchingsAtoms = allMatchingsAtoms[i];
        const auto &Barycenter = barycentersList[i];
        const auto &checkerAtoms = checkerAtomsList[i];
        gradActor.executeAtoms(
          dictDiagrams, matchingsAtoms, Barycenter, gradsAtoms, checkerAtoms,
          projForDiag, featuresToAdd, projLocations, vectorForProjContrib,
          pairToAddGradList, infoToAdd);
      }

      if(CreationFeatures_) {
        for(size_t i = 0; i < nDiags; ++i) {
          auto &projForDiag = allProjectionsList[i];
          auto &featuresToAdd = allFeaturesToAdd[i];
          auto &projLocations = allProjLocations[i];
          auto &vectorForProjContrib = allVectorForProjContributions[i];
          for(size_t j = 0; j < projForDiag.size(); ++j) {
            auto &t = featuresToAdd[j];
            std::array<double, 2> &pair = projLocations[j];
            std::vector<double> &vectorContrib = vectorForProjContrib[j];
            std::vector<int> &projAndIndex = projForDiag[j];
            std::vector<int> proj(numAtom);
            for(int m = 0; m < numAtom; ++m) {
              proj[m] = projAndIndex[m];
            }
            int atomIndex = static_cast<int>(projAndIndex[numAtom]);
            bool lenNull = allTrueProj[atomIndex].size() == 0;
            if(lenNull) {
              const CriticalType c1 = t.birth.type;
              const CriticalType c2 = t.death.type;
              const SimplexId idTemp = t.dim;
              pair[0] = pair[0] - step * vectorContrib[0];
              pair[1] = pair[1] - step * vectorContrib[1];
              if(pair[0] > pair[1]) {
                pair[1] = pair[0];
              }
              if(pair[0] < minBirth[atomIndex]) {
                continue;
              }
              if(pair[1] - pair[0] < 1e-7) {
                continue;
              }

              PersistencePair newPair{CriticalVertex{0, c1, pair[0], {}},
                                      CriticalVertex{0, c2, pair[1], {}},
                                      idTemp, true};
              allTrueProj[atomIndex].push_back(proj);
              allTrueFeaturesToAdd[atomIndex].push_back(newPair);
              allTrueProjLoc[atomIndex].push_back(pair);
            } else {
              bool ralph = true;
              size_t index = 0;
              if(Fusion_) {
                for(size_t n = 0; n < allTrueProj[atomIndex].size(); ++n) {
                  auto &projStocked = allTrueProj[atomIndex][n];
                  auto &projLocStocked = allTrueProj[atomIndex][n];
                  double distance
                    = sqrt(pow((pair[0] - projLocStocked[0]), 2)
                           + pow((pair[1] - projLocStocked[1]), 2));
                  if(proj == projStocked && distance < 1e-3) {
                    ralph = false;
                    index = n;
                    break;
                  }
                }
              }
              if(ralph) {

                const CriticalType c1 = t.birth.type;
                const CriticalType c2 = t.death.type;
                const SimplexId idTemp = t.dim;
                pair[0] = pair[0] - step * vectorContrib[0];
                pair[1] = pair[1] - step * vectorContrib[1];
                if(pair[0] > pair[1]) {
                  pair[1] = pair[0];
                }
                if(pair[0] < minBirth[atomIndex]) {
                  continue;
                }
                if(pair[1] - pair[0] < 1e-7) {
                  continue;
                }

                PersistencePair newPair{CriticalVertex{0, c1, pair[0], {}},
                                        CriticalVertex{0, c2, pair[1], {}},
                                        idTemp, true};
                allTrueProj[atomIndex].push_back(proj);
                allTrueFeaturesToAdd[atomIndex].push_back(newPair);
                allTrueProjLoc[atomIndex].push_back(pair);

              } else {
                auto &tReal = allTrueFeaturesToAdd[atomIndex][index];
                tReal.birth.sfValue
                  = tReal.birth.sfValue - step * vectorContrib[0];
                tReal.death.sfValue
                  = tReal.death.sfValue - step * vectorContrib[1];
                if(tReal.birth.sfValue > tReal.death.sfValue) {
                  tReal.death.sfValue = tReal.birth.sfValue;
                }
              }
            }
          }
        }

        if(!doCompression) {
          for(int i = 0; i < numAtom; ++i) {
            auto &atom = dictDiagrams[i];
            auto &histoEpochAtom = histoAllEpochLife[i];
            auto &histoBoolAtom = histoAllBoolLife[i];
            auto &boolUnderDiag = checkUnderDiag[i];
            auto &boolDiag = checkDiag[i];
            auto initSize = initSizes[i];
            auto &boolAboveGlobal = checkAboveGlobal[i];
            if(histoEpochAtom.size() > 0) {
              for(size_t j = 0; j < histoEpochAtom.size(); ++j) {
                auto &t = atom[initSize + j];
                histoEpochAtom[j] += 1;
                histoBoolAtom[j] = t.death.sfValue - t.birth.sfValue
                                   < 0.1 * (percent / 100.) * maxiDeath[i];
                boolDiag[j] = t.death.sfValue - t.birth.sfValue < 1e-6;
                boolUnderDiag[j] = t.death.sfValue < t.birth.sfValue;
                boolAboveGlobal[j] = t.birth.sfValue > maxiDeath[i];
              }
            }
          }

          for(int i = 0; i < numAtom; ++i) {
            auto &atom = dictDiagrams[i];
            auto &histoEpochAtom = histoAllEpochLife[i];
            auto &histoBoolAtom = histoAllBoolLife[i];
            auto &boolUnderDiag = checkUnderDiag[i];
            auto &boolDiag = checkDiag[i];
            auto &trueFeaturesToAdd = allTrueFeaturesToAdd[i];
            auto &boolAboveGlobal = checkAboveGlobal[i];
            for(size_t j = 0; j < trueFeaturesToAdd.size(); ++j) {
              auto &t = trueFeaturesToAdd[j];
              atom.push_back(t);
              histoEpochAtom.push_back(0);
              histoBoolAtom.push_back(t.death.sfValue - t.birth.sfValue
                                      < 0.1 * (percent / 100.) * maxiDeath[i]);
              boolDiag.push_back(t.death.sfValue - t.birth.sfValue < 1e-6);
              boolUnderDiag.push_back(t.death.sfValue < t.birth.sfValue);
              boolAboveGlobal.push_back(t.birth.sfValue > maxiDeath[i]);
            }
          }

          std::vector<std::vector<size_t>> allIndicesToDelete(numAtom);
          for(int i = 0; i < numAtom; ++i) {
            auto &indicesAtomToDelete = allIndicesToDelete[i];
            auto &histoEpochAtom = histoAllEpochLife[i];
            auto &histoBoolAtom = histoAllBoolLife[i];
            auto &boolUnderDiag = checkUnderDiag[i];
            auto &boolDiag = checkDiag[i];
            auto &boolAboveGlobal = checkAboveGlobal[i];
            for(size_t j = 0; j < histoEpochAtom.size(); ++j) {
              if(boolUnderDiag[j] || boolDiag[j] || boolAboveGlobal[j]
                 || (histoEpochAtom[j] > 5 && histoBoolAtom[j])) {
                indicesAtomToDelete.push_back(j);
              }
            }
          }

          for(int i = 0; i < numAtom; ++i) {
            auto &atom = dictDiagrams[i];
            auto &histoEpochAtom = histoAllEpochLife[i];
            auto &histoBoolAtom = histoAllBoolLife[i];
            auto &indicesAtomToDelete = allIndicesToDelete[i];
            auto &boolUnderDiag = checkUnderDiag[i];
            auto &boolDiag = checkDiag[i];
            auto &boolAboveGlobal = checkAboveGlobal[i];
            auto initSize = initSizes[i];
            // size_t initAtomSize = atom.size() -histoEpochAtom.size();
            if(static_cast<int>(indicesAtomToDelete.size()) > 0) {
              for(int j = static_cast<int>(indicesAtomToDelete.size()) - 1;
                  j >= 0; j--) {
                atom.erase(atom.begin() + initSize + indicesAtomToDelete[j]);
                histoEpochAtom.erase(histoEpochAtom.begin()
                                     + indicesAtomToDelete[j]);
                histoBoolAtom.erase(histoBoolAtom.begin()
                                    + indicesAtomToDelete[j]);
                boolUnderDiag.erase(boolUnderDiag.begin()
                                    + indicesAtomToDelete[j]);
                boolDiag.erase(boolDiag.begin() + indicesAtomToDelete[j]);
                boolAboveGlobal.erase(boolAboveGlobal.begin()
                                      + indicesAtomToDelete[j]);
              }
            } else {
              continue;
            }
          }
        } else {
          for(int i = 0; i < numAtom; ++i) {
            auto &atom = dictDiagrams[i];
            auto &trueFeaturesToAdd = allTrueFeaturesToAdd[i];
            for(size_t j = 0; j < trueFeaturesToAdd.size(); ++j) {
              auto &t = trueFeaturesToAdd[j];
              atom.push_back(t);
            }
          }
          controlAtomsSize(intermediateDiagrams, dictDiagrams);
        }
      } else {
        if(doCompression) {
          if(epoch > -1) {
            controlAtomsSize(intermediateDiagrams, dictDiagrams);
          }
        }
      }
      // Deleting unallowed pairs:
      for(int i = 0; i < numAtom; ++i) {
        auto &atom = dictDiagrams[i];
        auto &globalPair = atom[0];
        for(size_t j = 0; j < atom.size(); ++j) {
          auto &t = atom[j];

          if(t.death.sfValue > globalPair.death.sfValue) {
            t.death.sfValue = globalPair.death.sfValue;
          }

          if(t.birth.sfValue > t.death.sfValue) {
            t.death.sfValue = t.birth.sfValue;
          }

          if(t.birth.sfValue < globalPair.birth.sfValue) {
            t.birth.sfValue = globalPair.birth.sfValue;
          }
        }
      }
      this->printMsg("Computed 2nd opt for epoch " + std::to_string(epoch),
                     epoch / static_cast<double>(maxEpoch),
                     tm_opt2.getElapsedTime(), threadNumber_,
                     debug::LineMode::NEW, debug::Priority::DETAIL);

      // ATOM OPTIMIZATION
    }
    epoch += 1;
    barycentersList.clear();
    barycentersList.resize(nDiags);
    allMatchingsAtoms.clear();
    allMatchingsAtoms.resize(nDiags);
  }
  printMsg(
    " Epoch " + std::to_string(epoch) + ", loss = " + std::to_string(loss), 1,
    threadNumber_);

  printMsg("Loss returned "
           + std::to_string(*std::min_element(
             lossTab.begin() + nbEpochPrevious, lossTab.end()))
           + " at Epoch "
           + std::to_string(
             std::min_element(lossTab.begin() + nbEpochPrevious, lossTab.end())
             - lossTab.begin()));

  for(size_t p = 0; p < dictDiagrams.size(); ++p) {
    auto atom = histoDictDiagrams[p];
    dictDiagrams[p] = atom;
  }
  for(size_t p = 0; p < nDiags; ++p) {
    auto weights = histoVectorWeights[p];
    vectorWeights[p] = weights;
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);
}

double PersistenceDiagramDictionary::distVect(
  const std::vector<double> &vec1, const std::vector<double> &vec2) const {

  double dist = 0.;
  for(size_t i = 0; i < vec1.size(); ++i) {
    dist = dist + (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
  }
  return std::sqrt(dist);
}

double PersistenceDiagramDictionary::getMostPersistent(
  const std::vector<BidderDiagram> &bidder_diags) const {

  double max_persistence = 0;

  for(unsigned int i = 0; i < bidder_diags.size(); ++i) {
    for(size_t j = 0; j < bidder_diags[i].size(); ++j) {
      const double persistence = bidder_diags[i].at(j).getPersistence();
      if(persistence > max_persistence) {
        max_persistence = persistence;
      }
    }
  }

  return max_persistence;
}

double PersistenceDiagramDictionary::computeDistance(
  const BidderDiagram &D1,
  const BidderDiagram &D2,
  std::vector<ttk::MatchingType> &matching) const {

  GoodDiagram D2_bis{};
  for(size_t i = 0; i < D2.size(); i++) {
    const Bidder &b = D2.at(i);
    Good g(b.x_, b.y_, b.isDiagonal(), D2_bis.size());
    g.SetCriticalCoordinates(b.coords_);
    g.setPrice(0);
    D2_bis.emplace_back(g);
  }

  PersistenceDiagramAuction auction(
    this->Wasserstein, this->Alpha, this->Lambda, this->DeltaLim, true);
  auction.BuildAuctionDiagrams(D1, D2_bis);
  double loss = auction.run(matching);
  return loss;
}

void PersistenceDiagramDictionary::computeGradientWeights(
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
  const bool doOptimizeAtoms) const {

  // initialization
  std::vector<std::vector<std::array<double, 2>>> gradBuffersList(
    Barycenter.size());
  std::vector<std::vector<std::array<double, 2>>> pairToAddGradList;
  for(size_t i = 0; i < gradBuffersList.size(); ++i) {
    gradBuffersList[i].resize(matchingsAtoms.size());
  }
  std::vector<std::array<double, 2>> directions(Barycenter.size());
  std::vector<std::array<double, 2>> dataAssigned(Barycenter.size());
  gradWeights.resize(dictDiagrams.size());
  for(size_t i = 0; i < dictDiagrams.size(); ++i) {
    gradWeights[i] = 0.;
  }

  std::vector<std::vector<int>> checker(Barycenter.size());
  for(size_t j = 0; j < Barycenter.size(); ++j) {
    checker[j].resize(matchingsAtoms.size());
  }
  std::vector<int> tracker(Barycenter.size(), 0);
  std::vector<int> tracker2(Barycenter.size(), 0);

  // computing gradients
  for(size_t i = 0; i < matchingsAtoms.size(); ++i) {
    for(size_t j = 0; j < matchingsAtoms[i].size(); ++j) {
      const ttk::MatchingType &t = matchingsAtoms[i][j];
      // Id in atom
      const SimplexId Id1 = std::get<0>(t);
      // Id in barycenter
      const SimplexId Id2 = std::get<1>(t);
      // if(Id2 < 0) {
      if(Id2 < 0 || static_cast<int>(gradBuffersList.size()) <= Id2
         || static_cast<int>(dictDiagrams[i].size()) <= Id1) {
        continue;
      } else if(Id1 < 0) {
        const PersistencePair &t3 = Barycenter[Id2];
        auto &point = gradBuffersList[Id2][i];
        const double birthBarycenter = t3.birth.sfValue;
        const double deathBarycenter = t3.death.sfValue;
        const double birth_death_atom
          = birthBarycenter + (deathBarycenter - birthBarycenter) / 2.;
        point[0] = birth_death_atom;
        point[1] = birth_death_atom;
        checker[Id2][i] = i;
        tracker[Id2] = 1;
      } else {
        const PersistencePair &t2 = dictDiagrams[i][Id1];
        auto &point = gradBuffersList[Id2][i];
        const double birth_atom = t2.birth.sfValue;
        const double death_atom = t2.death.sfValue;
        point[0] = birth_atom;
        point[1] = death_atom;
        checker[Id2][i] = i;
        tracker[Id2] = 1;
      }
    }
  }

  // Compute directions for min diagram
  computeDirectionsGradWeight(matchingsAtoms, Barycenter, newData, matchingsMin,
                              indexBaryMin, indexDataMin, pairToAddGradList,
                              directions, dataAssigned, tracker2,
                              doOptimizeAtoms);

  // Compute directions for max diagram
  computeDirectionsGradWeight(matchingsAtoms, Barycenter, newData, matchingsMax,
                              indexBaryMax, indexDataMax, pairToAddGradList,
                              directions, dataAssigned, tracker2,
                              doOptimizeAtoms);

  // Compute directions for sad diagram
  computeDirectionsGradWeight(matchingsAtoms, Barycenter, newData, matchingsSad,
                              indexBarySad, indexDataSad, pairToAddGradList,
                              directions, dataAssigned, tracker2,
                              doOptimizeAtoms);

  std::vector<int> temp(pairToAddGradList.size(), 1);
  gradBuffersList.insert(
    gradBuffersList.end(), pairToAddGradList.begin(), pairToAddGradList.end());
  tracker2.insert(tracker2.end(), temp.begin(), temp.end());
  tracker.insert(tracker.end(), temp.begin(), temp.end());
  std::vector<int> temp2(matchingsAtoms.size());
  for(size_t j = 0; j < matchingsAtoms.size(); ++j) {
    temp2.push_back(static_cast<int>(j));
  }
  for(size_t j = 0; j < pairToAddGradList.size(); ++j) {
    checker.push_back(temp2);
  }

  for(size_t i = 0; i < gradBuffersList.size(); ++i) {
    const auto &data_point = dataAssigned[i];
    for(size_t j = 0; j < checker[i].size(); ++j) {
      auto &point = gradBuffersList[i][checker[i][j]];
      point[0] -= data_point[0];
      point[1] -= data_point[1];
    }
  }

  for(size_t i = 0; i < gradBuffersList.size(); ++i) {
    if(tracker[i] == 0 || tracker2[i] == 0) {
      continue;
    } else {
      for(size_t j = 0; j < checker[i].size(); ++j) {
        const auto &point = gradBuffersList[i][checker[i][j]];
        const auto &direction = directions[i];
        gradWeights[checker[i][j]]
          += -2 * (point[0] * direction[0] + point[1] * direction[1]);
      }
    }
  }
  hessianList.resize(gradBuffersList.size());
  for(size_t i = 0; i < gradBuffersList.size(); ++i) {
    Matrix &hessian = hessianList[i];
    hessian.resize(checker[i].size());
    for(size_t j = 0; j < checker[i].size(); ++j) {
      auto &line = hessian[j];
      line.resize(checker[i].size());
      const auto &point = gradBuffersList[i][checker[i][j]];
      for(size_t q = 0; q < checker[i].size(); ++q) {
        const auto &point_temp = gradBuffersList[i][checker[i][q]];
        line[q] = point[0] * point_temp[0] + point[1] * point_temp[1];
      }
    }
  }
}

void PersistenceDiagramDictionary::computeGradientAtoms(
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
  std::vector<PersistencePair> &infoToAdd,
  bool doDimReduct) const {

  gradsAtoms.resize(Barycenter.size());
  for(size_t i = 0; i < Barycenter.size(); ++i) {
    gradsAtoms[i].resize(weights.size());
  }

  checker.resize(Barycenter.size());
  for(size_t i = 0; i < Barycenter.size(); ++i) {
    checker[i] = 0;
  }
  std::vector<std::vector<double>> directions(Barycenter.size());

  // Compute directions for min diagram
  computeDirectionsGradAtoms(gradsAtoms, Barycenter, weights, newData,
                             matchingsMin, indexBaryMin, indexDataMin,
                             pairToAddGradList, directions, checker, infoToAdd,
                             doDimReduct);

  // Compute directions for max diagram
  computeDirectionsGradAtoms(gradsAtoms, Barycenter, weights, newData,
                             matchingsMax, indexBaryMax, indexDataMax,
                             pairToAddGradList, directions, checker, infoToAdd,
                             doDimReduct);

  // Compute directions for sad diagram
  computeDirectionsGradAtoms(gradsAtoms, Barycenter, weights, newData,
                             matchingsSad, indexBarySad, indexDataSad,
                             pairToAddGradList, directions, checker, infoToAdd,
                             doDimReduct);

  for(size_t i = 0; i < Barycenter.size(); ++i) {
    if(checker[i] == 0) {
      continue;
    } else {
      for(size_t j = 0; j < weights.size(); ++j) {
        std::vector<double> temp(2);
        const std::vector<double> &direction = directions[i];
        temp[0] = -2 * weights[j] * direction[0];
        temp[1] = -2 * weights[j] * direction[1];
        gradsAtoms[i][j] = temp;
      }
    }
  }
}

void PersistenceDiagramDictionary::setBidderDiagrams(
  const size_t nInputs,
  std::vector<ttk::DiagramType> &inputDiagrams,
  std::vector<BidderDiagram> &bidder_diags) const {

  bidder_diags.resize(nInputs);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; i++) {
    auto &diag = inputDiagrams[i];
    auto &bidders = bidder_diags[i];

    for(size_t j = 0; j < diag.size(); j++) {
      // Add bidder to bidders
      Bidder b(diag[j], j, this->Lambda);
      b.setPositionInAuction(bidders.size());
      bidders.emplace_back(b);
      if(b.isDiagonal() || b.x_ == b.y_) {
        this->printMsg("Diagonal point in diagram " + std::to_string(i) + "!",
                       ttk::debug::Priority::DETAIL);
      }
    }
  }
}

int PersistenceDiagramDictionary::initDictionary(
  std::vector<ttk::DiagramType> &dictDiagrams,
  const std::vector<ttk::DiagramType> &datas,
  const std::vector<ttk::DiagramType> &inputAtoms,
  const int nbAtom,
  bool do_min,
  bool do_sad,
  bool do_max,
  int seed,
  double percent) {
  switch(this->BackEnd) {
    case BACKEND::INPUT_ATOMS: {
      if(static_cast<int>(inputAtoms.size()) != nbAtom) {
        printErr("number of atoms don't match, going with "
                 + std::to_string(inputAtoms.size()) + " atoms");
      }
      for(size_t i = 0; i < inputAtoms.size(); i++) {
        const auto &t = inputAtoms[i];
        double maxPers = getMaxPers(t);

        dictDiagrams.push_back(t);
        dictDiagrams[i].erase(
          std::remove_if(dictDiagrams[i].begin(), dictDiagrams[i].end(),
                         [maxPers](ttk::PersistencePair &p) {
                           return (p.death.sfValue - p.birth.sfValue)
                                  < 0. * maxPers;
                         }),
          dictDiagrams[i].end());
      }
      break;
    }

    case BACKEND::BORDER_INIT: {
      ttk::InitDictPersistenceDiagram initializer;
      initializer.setThreadNumber(this->threadNumber_);
      initializer.execute(dictDiagrams, datas, nbAtom, do_min, do_sad, do_max);
      break;
    }

    case BACKEND::RANDOM_INIT:
      ttk::InitRandomDict{}.execute(dictDiagrams, datas, nbAtom, seed); // TODO
      break;

    case BACKEND::FIRST_DIAGS: {
      for(int i = 0; i < nbAtom; ++i) {
        const auto &t = datas[i];
        dictDiagrams.push_back(t);
      }
      break;
    }
    case BACKEND::GREEDY_INIT: {
      dictDiagrams.resize(datas.size());
      for(size_t i = 0; i < datas.size(); ++i) {
        dictDiagrams[i] = datas[i];
      }
      while(static_cast<int>(dictDiagrams.size()) != nbAtom) {
        std::vector<double> allEnergy(dictDiagrams.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(size_t j = 0; j < dictDiagrams.size(); ++j) {
          std::vector<ttk::DiagramType> dictTemp;
          std::vector<ttk::DiagramType> dataAlone;
          std::vector<double> lossTabTemp;
          std::vector<std::vector<double>> allLossesTemp(1);

          for(size_t p = 0; p < dictDiagrams.size(); ++p) {
            if(p != j) {
              dictTemp.push_back(dictDiagrams[p]);
            } else {
              dataAlone.push_back(dictDiagrams[p]);
            }
          }
          std::vector<std::vector<double>> weightsTemp(1);
          std::vector<double> weights(
            dictTemp.size(), 1. / (static_cast<int>(dictTemp.size()) * 1.));
          weightsTemp[0] = weights;
          std::vector<std::vector<double>> histoVectorWeights(1);
          std::vector<ttk::DiagramType> histoDictDiagrams(dictTemp.size());
          bool doCompression = false;
          this->method(dataAlone, dictTemp, weightsTemp,
                       static_cast<int>(dictTemp.size()), lossTabTemp,
                       allLossesTemp, histoVectorWeights, histoDictDiagrams,
                       false, percent, doCompression);
          double min_loss
            = *std::min_element(lossTabTemp.begin(), lossTabTemp.end());
          allEnergy[j] = min_loss;
        }
        int index = std::min_element(allEnergy.begin(), allEnergy.end())
                    - allEnergy.begin();
        dictDiagrams.erase(dictDiagrams.begin() + index);
      }
      break;
    }
    default:
      break;
  }
  return 0;
}

void PersistenceDiagramDictionary::gettingBidderDiagrams(
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
  bool insertOriginIndexMode) const {

  size_t nDiags = intermediateDiagrams.size();

  // Create diagrams for min, saddle and max persistence pairs
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDiags; i++) {
    const auto &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      const PersistencePair &t = CTDiagram[j];
      const ttk::CriticalType nt1 = t.birth.type;
      const ttk::CriticalType nt2 = t.death.type;
      const double pers = t.persistence();
      if(pers > 0) {
        if(nt1 == CriticalType::Local_minimum
           && nt2 == CriticalType::Local_maximum) {
          inputDiagramsMin[i].emplace_back(t);
          if(insertOriginIndexMode) {
            originIndexMin[i].push_back(j);
          }
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax[i].emplace_back(t);
            if(insertOriginIndexMode) {
              originIndexMax[i].push_back(j);
            }
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            inputDiagramsMin[i].emplace_back(t);
            if(insertOriginIndexMode) {
              originIndexMin[i].push_back(j);
            }
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            inputDiagramsSad[i].emplace_back(t);
            if(insertOriginIndexMode) {
              originIndexSad[i].push_back(j);
            }
          }
        }
      }
    }
  }

  if(this->do_min_) {
    setBidderDiagrams(nDiags, inputDiagramsMin, bidderDiagramsMin);
  }
  if(this->do_sad_) {
    setBidderDiagrams(nDiags, inputDiagramsSad, bidderDiagramsSad);
  }
  if(this->do_max_) {
    setBidderDiagrams(nDiags, inputDiagramsMax, bidderDiagramsMax);
  }
}

double PersistenceDiagramDictionary::getMaxPers(const ttk::DiagramType &data) {
  double maxPers = 0.;
  for(size_t j = 0; j < data.size(); ++j) {
    auto &t = data[j];
    double pers = t.death.sfValue - t.birth.sfValue;
    if(pers > maxPers) {
      maxPers = pers;
    }
  }

  return maxPers;
}

void PersistenceDiagramDictionary::controlAtomsSize(
  const std::vector<ttk::DiagramType> &intermediateDiagrams,
  std::vector<ttk::DiagramType> &dictDiagrams) const {

  size_t m = dictDiagrams.size();
  int globalSize = 0;

  for(size_t j = 0; j < intermediateDiagrams.size(); ++j) {
    auto &data = intermediateDiagrams[j];
    globalSize += static_cast<int>(data.size());
  }

  int dictSize = 0;

  for(size_t j = 0; j < m; ++j) {
    auto &atom = dictDiagrams[j];
    dictSize += static_cast<int>(atom.size());
  }

  if(static_cast<double>(dictSize)
     > (1. / this->CompressionFactor) * static_cast<double>(globalSize)) {
    double factor
      = (1. / this->CompressionFactor)
        * (static_cast<double>(globalSize) / static_cast<double>(dictSize));
    std::vector<std::vector<double>> tempDictPersistencePairs(m);

    for(size_t j = 0; j < m; ++j) {
      auto &temp = tempDictPersistencePairs[j];
      auto &atom = dictDiagrams[j];
      temp.resize(atom.size());
      for(size_t p = 0; p < atom.size(); ++p) {
        auto &t = atom[p];
        temp[p] = t.persistence();
      }
    }

    for(size_t j = 0; j < m; ++j) {
      auto &temp = tempDictPersistencePairs[j];
      std::sort(temp.begin(), temp.end(), std::greater<double>());
    }

    std::vector<double> persThresholds(m);
    for(size_t j = 0; j < m; ++j) {
      auto &temp = tempDictPersistencePairs[j];
      int index
        = static_cast<int>(floor(factor * static_cast<double>(temp.size())));
      persThresholds[j] = temp[index];
    }

    for(size_t j = 0; j < m; ++j) {
      double persThreshold = persThresholds[j];
      auto &atom = dictDiagrams[j];
      atom.erase(std::remove_if(atom.begin(), atom.end(),
                                [persThreshold](ttk::PersistencePair &t) {
                                  return (t.death.sfValue - t.birth.sfValue)
                                         < persThreshold;
                                }),
                 atom.end());
    }
  }
}

void PersistenceDiagramDictionary::computeDirectionsGradWeight(
  const std::vector<std::vector<ttk::MatchingType>> &matchingsAtoms,
  const ttk::DiagramType &Barycenter,
  const ttk::DiagramType &newData,
  const std::vector<ttk::MatchingType> &matchingsCritType,
  const std::vector<size_t> &indexBaryCritType,
  const std::vector<size_t> &indexDataCritType,
  std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
  std::vector<std::array<double, 2>> &directions,
  std::vector<std::array<double, 2>> &dataAssigned,
  std::vector<int> &tracker2,
  const bool doOptimizeAtoms) const {

  size_t m = matchingsCritType.size();
  // int k = 0;
  for(size_t i = 0; i < m; ++i) {
    const ttk::MatchingType &t = matchingsCritType[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);
    if(Id2 < 0) {

      if(Id1 < 0) {
        continue;
      } else {
        if(doOptimizeAtoms && CreationFeatures_ && ProgApproach_) {
          const PersistencePair &t2 = newData[indexDataCritType[Id1]];
          const double birthData = t2.birth.sfValue;
          const double deathData = t2.death.sfValue;
          const double birthDeathBarycenter
            = birthData + (deathData - birthData) / 2.;
          std::array<double, 2> direction;
          direction[0] = birthData - birthDeathBarycenter;
          direction[1] = deathData - birthDeathBarycenter;

          std::vector<std::array<double, 2>> newPairs(matchingsAtoms.size());
          for(size_t j = 0; j < matchingsAtoms.size(); ++j) {
            std::array<double, 2> pair{
              birthDeathBarycenter, birthDeathBarycenter};
            newPairs[j] = pair;
          }
          // pair to add later from the diagonal
          pairToAddGradList.push_back(newPairs);
          dataAssigned.push_back({birthData, deathData});
          directions.push_back(direction);
        } else {
          continue;
        }
      }
    } else {
      const PersistencePair &t3 = Barycenter[indexBaryCritType[Id2]];
      const double birthBarycenter = t3.birth.sfValue;
      const double deathBarycenter = t3.death.sfValue;
      auto &direction = directions[indexBaryCritType[Id2]];
      if(Id1 < 0) {
        // If matching on the diagonal
        const double birthDeathData
          = birthBarycenter + (deathBarycenter - birthBarycenter) / 2.;
        direction[0] = birthDeathData - birthBarycenter;
        direction[1] = birthDeathData - deathBarycenter;
        dataAssigned[indexBaryCritType[Id2]] = {birthDeathData, birthDeathData};

      } else {
        const PersistencePair &t2 = newData[indexDataCritType[Id1]];
        const double birthData = t2.birth.sfValue;
        const double deathData = t2.death.sfValue;
        direction[0] = birthData - birthBarycenter;
        direction[1] = deathData - deathBarycenter;
        dataAssigned[indexBaryCritType[Id2]] = {birthData, deathData};
      }
      tracker2[indexBaryCritType[Id2]] = 1;
    }
  }
}

void PersistenceDiagramDictionary::computeDirectionsGradAtoms(
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
  const bool doDimReduct) const {

  for(size_t i = 0; i < matchingsCritType.size(); ++i) {
    const ttk::MatchingType &t = matchingsCritType[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);
    int maxWeights
      = std::max_element(weights.begin(), weights.end()) - weights.begin();
    if(Id2 < 0) {
      if(Id1 < 0) {
        continue;
      } else {
        if(CreationFeatures_) {
          const PersistencePair &t2 = newData[indexDataCritType[Id1]];
          const double birthData = t2.birth.sfValue;
          const double deathData = t2.death.sfValue;
          const double birthDeathBarycenter
            = birthData + (deathData - birthData) / 2.;
          std::vector<double> direction(2);
          direction[0] = birthData - birthDeathBarycenter;
          direction[1] = deathData - birthDeathBarycenter;
          std::vector<std::vector<double>> temp3(weights.size());
          std::vector<double> temp2(2);
          if(CompressionMode_ && !doDimReduct) {
            for(size_t j = 0; j < weights.size(); ++j) {
              if(j == static_cast<size_t>(maxWeights)) {
                temp2[0] = -2 * weights[j] * direction[0];
                temp2[1] = -2 * weights[j] * direction[1];
              } else {
                temp2[0] = 0.;
                temp2[1] = 0.;
              }
              temp3[j] = temp2;
            }
          } else {
            for(size_t j = 0; j < weights.size(); ++j) {
              temp2[0] = -2 * weights[j] * direction[0];
              temp2[1] = -2 * weights[j] * direction[1];
              temp3[j] = temp2;
            }
          }
          gradsAtoms.push_back(temp3);
          checker.push_back(1);
          std::vector<std::array<double, 2>> newPairs(weights.size());
          for(size_t j = 0; j < weights.size(); ++j) {
            std::array<double, 2> pair{
              birthDeathBarycenter, birthDeathBarycenter};
            newPairs[j] = pair;
          }
          // pair to add later from the diagonal
          pairToAddGradList.push_back(newPairs);
          infoToAdd.push_back(t2);

        } else {
          continue;
        }
      }
    } else {
      const PersistencePair &t3 = Barycenter[indexBaryCritType[Id2]];
      const double birthBarycenter = t3.birth.sfValue;
      const double deathBarycenter = t3.death.sfValue;
      std::vector<double> direction(2);
      if(Id1 < 0) {
        const double birthDeathData
          = birthBarycenter + (deathBarycenter - birthBarycenter) / 2.;
        direction[0] = birthDeathData - birthBarycenter;
        direction[1] = birthDeathData - deathBarycenter;

      } else {
        const PersistencePair &t2 = newData[indexDataCritType[Id1]];
        const double birthData = t2.birth.sfValue;
        const double deathData = t2.death.sfValue;
        direction[0] = birthData - birthBarycenter;
        direction[1] = deathData - deathBarycenter;
      }
      directions[indexBaryCritType[Id2]] = std::move(direction);
      checker[indexBaryCritType[Id2]] = 1;
    }
  }
}

void PersistenceDiagramDictionary::computeAllDistances(
  std::vector<ttk::DiagramType> &barycentersList,
  const size_t nDiags,
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
  bool firstDistComputation) const {

  gettingBidderDiagrams(barycentersList, barycentersListMin, barycentersListSad,
                        barycentersListMax, bidderBarycentersListMin,
                        bidderBarycentersListSad, bidderBarycentersListMax,
                        originIndexBarysMin, originIndexBarysSad,
                        originIndexBarysMax, true);

  // Compute distance and matchings
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDiags; ++i) {
    std::vector<ttk::MatchingType> matchingMin;
    std::vector<ttk::MatchingType> matchingSad;
    std::vector<ttk::MatchingType> matchingMax;

    if(this->do_min_) {
      auto &barycentermin = bidderBarycentersListMin[i];
      auto &datamin = bidderDiagramsMin[i];

      if(firstDistComputation) {
        allLossesAtEpoch[i]
          += computeDistance(datamin, barycentermin, matchingMin);

      } else {
        computeDistance(datamin, barycentermin, matchingMin);
      }
    }
    if(this->do_max_) {
      auto &barycentermax = bidderBarycentersListMax[i];
      auto &datamax = bidderDiagramsMax[i];

      if(firstDistComputation) {
        allLossesAtEpoch[i]
          += computeDistance(datamax, barycentermax, matchingMax);

      } else {
        computeDistance(datamax, barycentermax, matchingMax);
      }
    }
    if(this->do_sad_) {
      auto &barycentersListad = bidderBarycentersListSad[i];
      auto &datasad = bidderDiagramsSad[i];

      if(firstDistComputation) {
        allLossesAtEpoch[i]
          += computeDistance(datasad, barycentersListad, matchingSad);

      } else {
        computeDistance(datasad, barycentersListad, matchingSad);
      }
    }
    matchingsDatasMin[i] = std::move(matchingMin);
    matchingsDatasSad[i] = std::move(matchingSad);
    matchingsDatasMax[i] = std::move(matchingMax);
  }
}