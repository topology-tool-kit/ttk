#include <MergeTreeNeuralNetwork.h>
#include <cmath>

#ifdef TTK_ENABLE_TORCH
using namespace torch::indexing;
#endif

ttk::MergeTreeNeuralNetwork::MergeTreeNeuralNetwork() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("MergeTreeNeuralNetwork");
}

#ifdef TTK_ENABLE_TORCH
//  -----------------------------------------------------------------------
//  --- Init
//  -----------------------------------------------------------------------
void ttk::MergeTreeNeuralNetwork::initInputBasis(
  unsigned int l,
  unsigned int layerNoAxes,
  std::vector<mtu::TorchMergeTree<float>> &tmTrees,
  std::vector<mtu::TorchMergeTree<float>> &tmTrees2,
  std::vector<bool> &ttkNotUsed(isTrain),
  std::vector<std::vector<torch::Tensor>> &allAlphasInit) {
  // TODO is there a way to avoid copy of merge trees?
  std::vector<ftm::MergeTree<float>> trees, trees2;
  for(unsigned int i = 0; i < tmTrees.size(); ++i) {
    trees.emplace_back(tmTrees[i].mTree);
    if(useDoubleInput_)
      trees2.emplace_back(tmTrees2[i].mTree);
  }

  // - Compute origin
  printMsg("Compute origin...", debug::Priority::DETAIL);
  Timer t_origin;
  std::vector<double> inputToBaryDistances;
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    baryMatchings, baryMatchings2;
  if(l != 0 or not layers_[0].getOrigin().tensor.defined()) {
    double sizeLimit = (l == 0 ? barycenterSizeLimitPercent_ : 0);
    unsigned int maxNoPairs
      = (l == 0 ? 0 : layers_[l - 1].getOriginPrime().tensor.sizes()[0] / 2);
    unsigned int maxNoPairs2
      = (l == 0 or not useDoubleInput_
           ? 0
           : layers_[l - 1].getOrigin2Prime().tensor.sizes()[0] / 2);
    layers_[l].initInputBasisOrigin(trees, trees2, sizeLimit, maxNoPairs,
                                    maxNoPairs2, inputToBaryDistances,
                                    baryMatchings, baryMatchings2);
    if(l == 0) {
      baryMatchings_L0_ = baryMatchings;
      baryMatchings2_L0_ = baryMatchings2;
      inputToBaryDistances_L0_ = inputToBaryDistances;
    }
  } else {
    baryMatchings = baryMatchings_L0_;
    baryMatchings2 = baryMatchings2_L0_;
    inputToBaryDistances = inputToBaryDistances_L0_;
  }
  printMsg("Compute origin time", 1, t_origin.getElapsedTime(), threadNumber_,
           debug::LineMode::NEW, debug::Priority::DETAIL);

  // - Compute vectors
  printMsg("Compute vectors...", debug::Priority::DETAIL);
  Timer t_vectors;
  std::vector<torch::Tensor> allAlphasInitT(tmTrees.size());
  layers_[l].initInputBasisVectors(
    tmTrees, tmTrees2, trees, trees2, layerNoAxes, allAlphasInitT,
    inputToBaryDistances, baryMatchings, baryMatchings2);
  for(unsigned int i = 0; i < allAlphasInitT.size(); ++i)
    allAlphasInit[i][l] = allAlphasInitT[i];
  printMsg("Compute vectors time", 1, t_vectors.getElapsedTime(), threadNumber_,
           debug::LineMode::NEW, debug::Priority::DETAIL);
}

void ttk::MergeTreeNeuralNetwork::initOutputBasis(
  unsigned int l,
  double layerOriginPrimeSizePercent,
  std::vector<mtu::TorchMergeTree<float>> &tmTrees,
  std::vector<mtu::TorchMergeTree<float>> &tmTrees2,
  std::vector<bool> &ttkNotUsed(isTrain)) {
  std::vector<ftm::FTMTree_MT *> ftmTrees(tmTrees.size()),
    ftmTrees2(tmTrees2.size());
  for(unsigned int i = 0; i < tmTrees.size(); ++i)
    ftmTrees[i] = &(tmTrees[i].mTree.tree);
  for(unsigned int i = 0; i < tmTrees2.size(); ++i)
    ftmTrees2[i] = &(tmTrees2[i].mTree.tree);
  auto sizeMetric = getSizeLimitMetric(ftmTrees);
  auto sizeMetric2 = getSizeLimitMetric(ftmTrees2);
  auto getDim = [](double _sizeMetric, double _percent) {
    unsigned int dim = std::max((int)(_sizeMetric * _percent / 100.0), 2) * 2;
    return dim;
  };

  unsigned int dim = getDim(sizeMetric, layerOriginPrimeSizePercent);
  dim = std::min(dim, (unsigned int)layers_[l].getOrigin().tensor.sizes()[0]);
  unsigned int dim2 = getDim(sizeMetric2, layerOriginPrimeSizePercent);
  if(useDoubleInput_)
    dim2
      = std::min(dim2, (unsigned int)layers_[l].getOrigin2().tensor.sizes()[0]);
  auto baseTensor = (l == 0 ? layers_[0].getOrigin().tensor
                            : layers_[l - 1].getOriginPrime().tensor);
  layers_[l].initOutputBasis(dim, dim2, baseTensor);
}

bool ttk::MergeTreeNeuralNetwork::initGetReconstructed(
  unsigned int l,
  unsigned int layerNoAxes,
  double layerOriginPrimeSizePercent,
  std::vector<mtu::TorchMergeTree<float>> &trees,
  std::vector<mtu::TorchMergeTree<float>> &trees2,
  std::vector<bool> &isTrain,
  std::vector<mtu::TorchMergeTree<float>> &recs,
  std::vector<mtu::TorchMergeTree<float>> &recs2,
  std::vector<std::vector<torch::Tensor>> &allAlphasInit) {
  printMsg("Get reconstructed", debug::Priority::DETAIL);
  recs.resize(trees.size());
  recs2.resize(trees.size());
  unsigned int i = 0;
  unsigned int noReset = 0;
  while(i < trees.size()) {
    layers_[l].outputBasisReconstruction(
      allAlphasInit[i][l], recs[i], recs2[i], activateOutputInit_);
    if(recs[i].mTree.tree.getRealNumberOfNodes() == 0) {
      bool fullReset = initResetOutputBasis(
        l, layerNoAxes, layerOriginPrimeSizePercent, trees, trees2, isTrain);
      if(fullReset)
        return true;
      i = 0;
      ++noReset;
      if(noReset >= 100) {
        printWrn("[initParameters] noReset >= 100");
        return true;
      }
    }
    ++i;
  }
  return false;
}

void ttk::MergeTreeNeuralNetwork::initStep(
  std::vector<mtu::TorchMergeTree<float>> &trees,
  std::vector<mtu::TorchMergeTree<float>> &trees2,
  std::vector<bool> &isTrain) {
  layers_.clear();

  float bestError = std::numeric_limits<float>::max();
  std::vector<torch::Tensor> bestVSTensor, bestVSPrimeTensor, bestVS2Tensor,
    bestVS2PrimeTensor, bestLatentCentroids;
  std::vector<mtu::TorchMergeTree<float>> bestOrigins, bestOriginsPrime,
    bestOrigins2, bestOrigins2Prime;
  std::vector<std::vector<torch::Tensor>> bestAlphasInit;
  for(unsigned int n = 0; n < noInit_; ++n) {
    // Init parameters
    float error = initParameters(trees, trees2, isTrain, (noInit_ != 1));
    // Save best parameters
    if(noInit_ != 1) {
      std::stringstream ss;
      ss << "Init error = " << error;
      printMsg(ss.str());
      if(error < bestError) {
        bestError = error;
        copyParams(bestOrigins, bestOriginsPrime, bestVSTensor,
                   bestVSPrimeTensor, bestOrigins2, bestOrigins2Prime,
                   bestVS2Tensor, bestVS2PrimeTensor, allAlphas_,
                   bestAlphasInit, true);
        copyCustomParams(true);
      }
    }
  }
  // TODO this copy can be avoided if initParameters takes dummy tensors to fill
  // as parameters and then copy to the member tensors when a better init is
  // found.
  if(noInit_ != 1) {
    // Put back best parameters
    std::stringstream ss;
    ss << "Best init error = " << bestError;
    printMsg(ss.str());
    copyParams(bestOrigins, bestOriginsPrime, bestVSTensor, bestVSPrimeTensor,
               bestOrigins2, bestOrigins2Prime, bestVS2Tensor,
               bestVS2PrimeTensor, bestAlphasInit, allAlphas_, false);
    copyCustomParams(false);
  }

  for(unsigned int l = 0; l < noLayers_; ++l) {
    layers_[l].requires_grad(true);

    // Print
    printMsg(debug::Separator::L2);
    std::stringstream ss;
    ss << "Layer " << l;
    printMsg(ss.str());
    if(isTreeHasBigValues(layers_[l].getOrigin().mTree, bigValuesThreshold_)) {
      ss.str("");
      ss << "origins_[" << l << "] has big values!" << std::endl;
      printMsg(ss.str());
      printPairs(layers_[l].getOrigin().mTree);
    }
    if(isTreeHasBigValues(
         layers_[l].getOriginPrime().mTree, bigValuesThreshold_)) {
      ss.str("");
      ss << "originsPrime_[" << l << "] has big values!" << std::endl;
      printMsg(ss.str());
      printPairs(layers_[l].getOriginPrime().mTree);
    }
    ss.str("");
    ss << "vS size   = " << layers_[l].getVSTensor().sizes();
    printMsg(ss.str());
    ss.str("");
    ss << "vS' size  = " << layers_[l].getVSPrimeTensor().sizes();
    printMsg(ss.str());
    if(trees2.size() != 0) {
      ss.str("");
      ss << "vS2 size  = " << layers_[l].getVS2Tensor().sizes();
      printMsg(ss.str());
      ss.str("");
      ss << "vS2' size = " << layers_[l].getVS2PrimeTensor().sizes();
      printMsg(ss.str());
    }
  }
}

void ttk::MergeTreeNeuralNetwork::passLayerParameters(
  MergeTreeNeuralLayer &layer) {
  layer.setDropout(dropout_);
  layer.setEuclideanVectorsInit(euclideanVectorsInit_);
  layer.setRandomAxesInit(randomAxesInit_);
  layer.setInitBarycenterRandom(initBarycenterRandom_);
  layer.setInitBarycenterOneIter(initBarycenterOneIter_);
  layer.setInitOriginPrimeStructByCopy(initOriginPrimeStructByCopy_);
  layer.setInitOriginPrimeValuesByCopy(initOriginPrimeValuesByCopy_);
  layer.setInitOriginPrimeValuesByCopyRandomness(
    initOriginPrimeValuesByCopyRandomness_);
  layer.setActivate(activate_);
  layer.setActivationFunction(activationFunction_);
  layer.setUseGpu(useGpu_);
  layer.setBigValuesThreshold(bigValuesThreshold_);

  layer.setDeterministic(deterministic_);
  layer.setNumberOfProjectionSteps(k_);
  layer.setBarycenterSizeLimitPercent(barycenterSizeLimitPercent_);
  layer.setProbabilisticVectorsInit(probabilisticVectorsInit_);

  layer.setNormalizedWasserstein(normalizedWasserstein_);
  layer.setAssignmentSolver(assignmentSolverID_);
  layer.setNodePerTask(nodePerTask_);
  layer.setUseDoubleInput(useDoubleInput_);
  layer.setJoinSplitMixtureCoefficient(mixtureCoefficient_);
  layer.setIsPersistenceDiagram(isPersistenceDiagram_);

  layer.setDebugLevel(debugLevel_);
  layer.setThreadNumber(threadNumber_);
}

//  ---------------------------------------------------------------------------
//  --- Forward
//  ---------------------------------------------------------------------------
bool ttk::MergeTreeNeuralNetwork::forwardOneData(
  mtu::TorchMergeTree<float> &tree,
  mtu::TorchMergeTree<float> &tree2,
  unsigned int treeIndex,
  unsigned int k,
  std::vector<torch::Tensor> &alphasInit,
  mtu::TorchMergeTree<float> &out,
  mtu::TorchMergeTree<float> &out2,
  std::vector<torch::Tensor> &dataAlphas,
  std::vector<mtu::TorchMergeTree<float>> &outs,
  std::vector<mtu::TorchMergeTree<float>> &outs2,
  bool train) {
  outs.resize(noLayers_ - 1);
  outs2.resize(noLayers_ - 1);
  dataAlphas.resize(noLayers_);
  for(unsigned int l = 0; l < noLayers_; ++l) {
    auto &treeToUse = (l == 0 ? tree : outs[l - 1]);
    auto &tree2ToUse = (l == 0 ? tree2 : outs2[l - 1]);
    auto &outToUse = (l != noLayers_ - 1 ? outs[l] : out);
    auto &out2ToUse = (l != noLayers_ - 1 ? outs2[l] : out2);
    bool reset = layers_[l].forward(treeToUse, tree2ToUse, k, alphasInit[l],
                                    outToUse, out2ToUse, dataAlphas[l], train);
    if(reset)
      return true;
    // Update recs
    auto updateRecs
      = [this, &treeIndex, &l](
          std::vector<std::vector<mtu::TorchMergeTree<float>>> &recs,
          mtu::TorchMergeTree<float> &outT) {
          if(recs[treeIndex].size() > noLayers_)
            mtu::copyTorchMergeTree<float>(outT, recs[treeIndex][l + 1]);
          else {
            mtu::TorchMergeTree<float> tmt;
            mtu::copyTorchMergeTree<float>(outT, tmt);
            recs[treeIndex].emplace_back(tmt);
          }
        };
    updateRecs(recs_, outToUse);
    if(useDoubleInput_)
      updateRecs(recs2_, out2ToUse);
  }
  return false;
}

bool ttk::MergeTreeNeuralNetwork::forwardStep(
  std::vector<mtu::TorchMergeTree<float>> &trees,
  std::vector<mtu::TorchMergeTree<float>> &trees2,
  std::vector<unsigned int> &indexes,
  std::vector<bool> &isTrain,
  unsigned int k,
  std::vector<std::vector<torch::Tensor>> &allAlphasInit,
  bool computeError,
  std::vector<mtu::TorchMergeTree<float>> &outs,
  std::vector<mtu::TorchMergeTree<float>> &outs2,
  std::vector<std::vector<torch::Tensor>> &bestAlphas,
  std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts,
  std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts2,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &matchings,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &matchings2,
  float &loss,
  float &testLoss) {
  loss = 0;
  testLoss = 0;
  outs.resize(trees.size());
  outs2.resize(trees.size());
  bestAlphas.resize(trees.size());
  layersOuts.resize(trees.size());
  layersOuts2.resize(trees.size());
  matchings.resize(trees.size());
  if(useDoubleInput_)
    matchings2.resize(trees2.size());
  mtu::TorchMergeTree<float> dummyTMT;
  bool reset = false;
  unsigned int noTrainLoss = 0, noTestLoss = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic)                                \
  num_threads(this->threadNumber_) if(parallelize_) reduction(|| : reset) \
  reduction(+ : loss)
#endif
  for(unsigned int ind = 0; ind < indexes.size(); ++ind) {
    unsigned int i = indexes[ind];
    auto &tree2ToUse = (trees2.size() == 0 ? dummyTMT : trees2[i]);
    bool dReset = forwardOneData(trees[i], tree2ToUse, i, k, allAlphasInit[i],
                                 outs[i], outs2[i], bestAlphas[i],
                                 layersOuts[i], layersOuts2[i], isTrain[i]);
    if(computeError) {
      float iLoss
        = computeOneLoss(trees[i], outs[i], trees2[i], outs2[i], matchings[i],
                         matchings2[i], bestAlphas[i], i);
      if(isTrain[i]) {
        loss += iLoss;
        ++noTrainLoss;
      } else {
        testLoss += iLoss;
        ++noTestLoss;
      }
    }
    if(dReset)
      reset = reset || dReset;
  }
  if(noTrainLoss != 0)
    loss /= noTrainLoss;
  if(noTestLoss != 0)
    testLoss /= noTestLoss;
  return reset;
}

bool ttk::MergeTreeNeuralNetwork::forwardStep(
  std::vector<mtu::TorchMergeTree<float>> &trees,
  std::vector<mtu::TorchMergeTree<float>> &trees2,
  std::vector<unsigned int> &indexes,
  unsigned int k,
  std::vector<std::vector<torch::Tensor>> &allAlphasInit,
  bool computeError,
  std::vector<mtu::TorchMergeTree<float>> &outs,
  std::vector<mtu::TorchMergeTree<float>> &outs2,
  std::vector<std::vector<torch::Tensor>> &bestAlphas,
  std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts,
  std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts2,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &matchings,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &matchings2,
  float &loss) {
  std::vector<bool> isTrain(trees.size(), false);
  float tempLoss;
  return forwardStep(trees, trees2, indexes, isTrain, k, allAlphasInit,
                     computeError, outs, outs2, bestAlphas, layersOuts,
                     layersOuts2, matchings, matchings2, tempLoss, loss);
}

//  ---------------------------------------------------------------------------
//  --- Projection
//  ---------------------------------------------------------------------------
void ttk::MergeTreeNeuralNetwork::projectionStep() {
  for(unsigned int l = 0; l < noLayers_; ++l)
    layers_[l].projectionStep();
}

//  -----------------------------------------------------------------------
//  --- Convergence
//  -----------------------------------------------------------------------
bool ttk::MergeTreeNeuralNetwork::isBestLoss(float loss,
                                             float &minLoss,
                                             unsigned int &cptBlocked) {
  bool isBestEnergy = false;
  if(loss + ENERGY_COMPARISON_TOLERANCE < minLoss) {
    minLoss = loss;
    cptBlocked = 0;
    isBestEnergy = true;
  }
  return isBestEnergy;
}

bool ttk::MergeTreeNeuralNetwork::convergenceStep(float loss,
                                                  float &oldLoss,
                                                  float &minLoss,
                                                  unsigned int &cptBlocked) {
  double tol = oldLoss / 125.0;
  bool converged = std::abs(loss - oldLoss) < std::abs(tol);
  oldLoss = loss;
  if(not converged) {
    cptBlocked += (minLoss < loss) ? 1 : 0;
    converged = (cptBlocked >= 10 * 10);
    if(converged)
      printMsg("Blocked!", debug::Priority::DETAIL);
  }
  return converged;
}

//  -----------------------------------------------------------------------
//  --- Main Functions
//  -----------------------------------------------------------------------
void ttk::MergeTreeNeuralNetwork::fit(
  std::vector<ftm::MergeTree<float>> &trees,
  std::vector<ftm::MergeTree<float>> &trees2) {
  torch::set_num_threads(1);
  if(useGpu_) {
    if(torch::cuda::device_count() > 0 and torch::cuda::is_available())
      printMsg("Computation with GPU support.");
    else {
      printMsg("Disabling GPU support because no device were found.");
      useGpu_ = false;
      // TODO cache useGpu parameter to be in accordance with ParaView GUI
    }
  } else {
    printMsg("Computation without GPU support.");
  }
  //  ----- Determinism
  if(deterministic_) {
    int m_seed = 0;
    bool m_torch_deterministic = true;
    srand(m_seed);
    torch::manual_seed(m_seed);
    at::globalContext().setDeterministicCuDNN(m_torch_deterministic ? true
                                                                    : false);
    if(not useGpu_)
      at::globalContext().setDeterministicAlgorithms(
        m_torch_deterministic ? true : false, true);
  }

  //  ----- Testing
  for(unsigned int i = 0; i < trees.size(); ++i) {
    for(unsigned int n = 0; n < trees[i].tree.getNumberOfNodes(); ++n) {
      if(trees[i].tree.isNodeAlone(n))
        continue;
      auto birthDeath = trees[i].tree.template getBirthDeath<float>(n);
      bigValuesThreshold_
        = std::max(std::abs(std::get<0>(birthDeath)), bigValuesThreshold_);
      bigValuesThreshold_
        = std::max(std::abs(std::get<1>(birthDeath)), bigValuesThreshold_);
    }
  }
  bigValuesThreshold_ *= 100;

  // ----- Convert MergeTree to TorchMergeTree
  std::vector<mtu::TorchMergeTree<float>> torchTrees, torchTrees2;
  mergeTreesToTorchTrees(trees, torchTrees, normalizedWasserstein_);
  mergeTreesToTorchTrees(trees2, torchTrees2, normalizedWasserstein_);
  if(useGpu_) {
    for(unsigned i = 0; i < torchTrees.size(); ++i)
      torchTrees[i].tensor = torchTrees[i].tensor.cuda();
    for(unsigned i = 0; i < torchTrees2.size(); ++i)
      torchTrees2[i].tensor = torchTrees2[i].tensor.cuda();
  }

  auto initRecs = [](std::vector<std::vector<mtu::TorchMergeTree<float>>> &recs,
                     std::vector<mtu::TorchMergeTree<float>> &torchTreesT) {
    recs.clear();
    recs.resize(torchTreesT.size());
    for(unsigned int i = 0; i < torchTreesT.size(); ++i) {
      mtu::TorchMergeTree<float> tmt;
      mtu::copyTorchMergeTree<float>(torchTreesT[i], tmt);
      recs[i].emplace_back(tmt);
    }
  };
  initRecs(recs_, torchTrees);
  if(useDoubleInput_)
    initRecs(recs2_, torchTrees2);

  // --- Train/Test Split
  unsigned int trainSize = std::min(
    std::max((int)(trees.size() * trainTestSplit_), 1), (int)trees.size());
  std::vector<unsigned int> trainIndexes(trees.size()), testIndexes;
  std::iota(trainIndexes.begin(), trainIndexes.end(), 0);
  std::random_device rd;
  std::default_random_engine rng(deterministic_ ? 0 : rd());
  bool trainTestSplitted = trainSize != trees.size();
  if(trainTestSplitted) {
    if(shuffleBeforeSplit_)
      std::shuffle(trainIndexes.begin(), trainIndexes.end(), rng);
    testIndexes.insert(
      testIndexes.end(), trainIndexes.begin() + trainSize, trainIndexes.end());
    trainIndexes.resize(trainSize);
  }
  std::vector<bool> isTrain(trees.size(), true);
  for(auto &ind : testIndexes)
    isTrain[ind] = false;

  // ----- Custom Init
  customInit(torchTrees, torchTrees2);

  // ----- Init Model Parameters
  Timer t_init;
  initStep(torchTrees, torchTrees2, isTrain);
  printMsg("Init", 1, t_init.getElapsedTime(), threadNumber_);

  // --- Init optimizer
  std::vector<torch::Tensor> parameters;
  for(unsigned int l = 0; l < noLayers_; ++l) {
    parameters.emplace_back(layers_[l].getOrigin().tensor);
    parameters.emplace_back(layers_[l].getOriginPrime().tensor);
    parameters.emplace_back(layers_[l].getVSTensor());
    parameters.emplace_back(layers_[l].getVSPrimeTensor());
    if(trees2.size() != 0) {
      parameters.emplace_back(layers_[l].getOrigin2().tensor);
      parameters.emplace_back(layers_[l].getOrigin2Prime().tensor);
      parameters.emplace_back(layers_[l].getVS2Tensor());
      parameters.emplace_back(layers_[l].getVS2PrimeTensor());
    }
  }
  addCustomParameters(parameters);

  torch::optim::Optimizer *optimizer;
  // - Init Adam
  auto adamOptions = torch::optim::AdamOptions(gradientStepSize_);
  adamOptions.betas(std::make_tuple(beta1_, beta2_));
  auto adamOptimizer = torch::optim::Adam(parameters, adamOptions);
  // - Init SGD optimizer
  auto sgdOptions = torch::optim::SGDOptions(gradientStepSize_);
  auto sgdOptimizer = torch::optim::SGD(parameters, sgdOptions);
  // -Init RMSprop optimizer
  auto rmspropOptions = torch::optim::RMSpropOptions(gradientStepSize_);
  auto rmspropOptimizer = torch::optim::RMSprop(parameters, rmspropOptions);
  // - Set optimizer pointer
  switch(optimizer_) {
    case 1:
      optimizer = &sgdOptimizer;
      break;
    case 2:
      optimizer = &rmspropOptimizer;
      break;
    case 0:
    default:
      optimizer = &adamOptimizer;
  }

  // --- Print train/test split
  if(trainTestSplitted) {
    std::stringstream ss;
    ss << "trainSize = " << trainIndexes.size() << " / " << trees.size();
    printMsg(ss.str());
    ss.str("");
    ss << "testSize = " << testIndexes.size() << " / " << trees.size();
    printMsg(ss.str());
  }

  // --- Init batches indexes
  unsigned int batchSize
    = std::min(std::max((int)(trainIndexes.size() * batchSize_), 1),
               (int)trainIndexes.size());
  std::stringstream ssBatch;
  ssBatch << "batchSize = " << batchSize;
  printMsg(ssBatch.str());
  unsigned int noBatch = trainIndexes.size() / batchSize
                         + ((trainIndexes.size() % batchSize) != 0 ? 1 : 0);
  std::vector<std::vector<unsigned int>> allIndexes(noBatch);
  if(noBatch == 1) {
    // Yes, trees.size() below is correct and it is not trainIndexes.size(), the
    // goal is to forward everyone (even test data) if noBatch == 1 to benefit
    // from full parallelism, but only train data will be used for backward.
    allIndexes[0].resize(trees.size());
    std::iota(allIndexes[0].begin(), allIndexes[0].end(), 0);
  }

  // ----- Testing
  originsNoZeroGrad_.resize(noLayers_);
  originsPrimeNoZeroGrad_.resize(noLayers_);
  vSNoZeroGrad_.resize(noLayers_);
  vSPrimeNoZeroGrad_.resize(noLayers_);
  for(unsigned int l = 0; l < noLayers_; ++l) {
    originsNoZeroGrad_[l] = 0;
    originsPrimeNoZeroGrad_[l] = 0;
    vSNoZeroGrad_[l] = 0;
    vSPrimeNoZeroGrad_[l] = 0;
  }
  if(useDoubleInput_) {
    origins2NoZeroGrad_.resize(noLayers_);
    origins2PrimeNoZeroGrad_.resize(noLayers_);
    vS2NoZeroGrad_.resize(noLayers_);
    vS2PrimeNoZeroGrad_.resize(noLayers_);
    for(unsigned int l = 0; l < noLayers_; ++l) {
      origins2NoZeroGrad_[l] = 0;
      origins2PrimeNoZeroGrad_[l] = 0;
      vS2NoZeroGrad_[l] = 0;
      vS2PrimeNoZeroGrad_[l] = 0;
    }
  }

  // ----- Init Variables
  unsigned int k = k_;
  float oldLoss, minLoss, minTestLoss;
  std::vector<float> minCustomLoss;
  unsigned int cptBlocked, iteration = 0;
  auto initLoop = [&]() {
    oldLoss = -1;
    minLoss = std::numeric_limits<float>::max();
    minTestLoss = std::numeric_limits<float>::max();
    cptBlocked = 0;
    iteration = 0;
  };
  initLoop();
  int convWinSize = 5;
  int noConverged = 0, noConvergedToGet = 10;
  std::vector<float> gapLosses, gapTestLosses;
  std::vector<std::vector<float>> gapCustomLosses;
  float windowLoss = 0;

  double assignmentTime = 0.0, updateTime = 0.0, projectionTime = 0.0,
         lossTime = 0.0;

  int bestIteration = 0;
  std::vector<torch::Tensor> bestVSTensor, bestVSPrimeTensor, bestVS2Tensor,
    bestVS2PrimeTensor;
  std::vector<mtu::TorchMergeTree<float>> bestOrigins, bestOriginsPrime,
    bestOrigins2, bestOrigins2Prime;
  std::vector<std::vector<torch::Tensor>> bestAlphasInit;
  std::vector<std::vector<mtu::TorchMergeTree<float>>> bestRecs, bestRecs2;
  double bestTime = 0;

  auto printLoss = [this, trainTestSplitted](
                     float loss, float testLoss, std::vector<float> &customLoss,
                     int iterationT, int iterationTT, double time,
                     const debug::Priority &priority = debug::Priority::INFO) {
    std::stringstream prefix;
    prefix << (priority == debug::Priority::VERBOSE ? "Iter " : "Best ");
    std::stringstream ssBestLoss;
    ssBestLoss << prefix.str() << "loss is " << loss << " (iteration "
               << iterationT << " / " << iterationTT << ") at time " << time;
    printMsg(ssBestLoss.str(), priority);
    if(trainTestSplitted) {
      ssBestLoss.str("");
      ssBestLoss << prefix.str() << "test loss is " << testLoss;
      printMsg(ssBestLoss.str(), priority);
    }
    printCustomLosses(customLoss, prefix, priority);
  };

  auto copyAlphas = [this](std::vector<std::vector<torch::Tensor>> &alphas,
                           std::vector<unsigned int> &indexes) {
    for(unsigned int ind = 0; ind < indexes.size(); ++ind) {
      unsigned int i = indexes[ind];
      for(unsigned int j = 0; j < alphas[i].size(); ++j)
        mtu::copyTensor(alphas[i][j], allAlphas_[i][j]);
    }
  };

  // ----- Algorithm
  Timer t_alg;
  bool converged = false;
  while(not converged) {
    if(iteration % iterationGap_ == 0) {
      std::stringstream ss;
      ss << "Iteration " << iteration;
      printMsg(debug::Separator::L2);
      printMsg(ss.str());
    }

    bool forwardReset = false;
    std::vector<float> iterationLosses, iterationTestLosses;
    std::vector<std::vector<float>> iterationCustomLosses;
    if(noBatch != 1) {
      std::vector<unsigned int> indexes = trainIndexes;
      std::shuffle(std::begin(indexes), std::end(indexes), rng);
      for(unsigned int i = 0; i < allIndexes.size(); ++i) {
        unsigned int noProcessed = batchSize * i;
        unsigned int remaining = trainIndexes.size() - noProcessed;
        unsigned int size = std::min(batchSize, remaining);
        allIndexes[i].resize(size);
        for(unsigned int j = 0; j < size; ++j)
          allIndexes[i][j] = indexes[noProcessed + j];
      }
    }
    for(unsigned batchNum = 0; batchNum < allIndexes.size(); ++batchNum) {
      auto &indexes = allIndexes[batchNum];

      // --- Forward
      Timer t_assignment;
      std::vector<mtu::TorchMergeTree<float>> outs, outs2;
      std::vector<std::vector<torch::Tensor>> bestAlphas;
      std::vector<std::vector<mtu::TorchMergeTree<float>>> layersOuts,
        layersOuts2;
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        matchings, matchings2;
      float loss, testLoss;
      bool computeError = true;
      forwardReset
        = forwardStep(torchTrees, torchTrees2, indexes, isTrain, k, allAlphas_,
                      computeError, outs, outs2, bestAlphas, layersOuts,
                      layersOuts2, matchings, matchings2, loss, testLoss);
      if(forwardReset)
        break;
      copyAlphas(bestAlphas, indexes);
      assignmentTime += t_assignment.getElapsedTime();

      // --- Loss
      Timer t_loss;
      gapLosses.emplace_back(loss);
      iterationLosses.emplace_back(loss);
      if(noBatch == 1 and trainTestSplitted) {
        gapTestLosses.emplace_back(testLoss);
        iterationTestLosses.emplace_back(testLoss);
      }
      std::vector<torch::Tensor> torchCustomLoss;
      computeCustomLosses(layersOuts, layersOuts2, bestAlphas, indexes, isTrain,
                          iteration, gapCustomLosses, iterationCustomLosses,
                          torchCustomLoss);
      lossTime += t_loss.getElapsedTime();

      // --- Backward
      Timer t_update;
      backwardStep(torchTrees, outs, matchings, torchTrees2, outs2, matchings2,
                   bestAlphas, *optimizer, indexes, isTrain, torchCustomLoss);
      updateTime += t_update.getElapsedTime();

      // --- Projection
      Timer t_projection;
      projectionStep();
      projectionTime += t_projection.getElapsedTime();
    } // end batch

    if(noBatch != 1 and trainTestSplitted) {
      std::vector<mtu::TorchMergeTree<float>> outs, outs2;
      std::vector<std::vector<torch::Tensor>> bestAlphas;
      std::vector<std::vector<mtu::TorchMergeTree<float>>> layersOuts,
        layersOuts2;
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        matchings, matchings2;
      float loss, testLoss;
      bool computeError = true;
      forwardStep(torchTrees, torchTrees2, testIndexes, isTrain, k, allAlphas_,
                  computeError, outs, outs2, bestAlphas, layersOuts,
                  layersOuts2, matchings, matchings2, loss, testLoss);
      copyAlphas(bestAlphas, testIndexes);
      gapTestLosses.emplace_back(testLoss);
      iterationTestLosses.emplace_back(testLoss);
      std::vector<torch::Tensor> torchCustomLoss;
      computeCustomLosses(layersOuts, layersOuts2, bestAlphas, testIndexes,
                          isTrain, iteration, gapCustomLosses,
                          iterationCustomLosses, torchCustomLoss);
    }

    if(forwardReset) {
      // TODO better manage reset by init new parameters and start again for
      // example (should not happen anymore)
      printWrn("Forward reset!");
      break;
    }

    // --- Get iteration loss
    // TODO an approximation is made here if batch size != 1 because the
    // iteration loss that will be printed will not be exact, we need to do a
    // forward step and compute loss with the whole dataset
    float iterationLoss = torch::tensor(iterationLosses).mean().item<float>();
    float iterationTestLoss
      = torch::tensor(iterationTestLosses).mean().item<float>();
    std::vector<float> iterationCustomLoss;
    float iterationTotalLoss = computeIterationTotalLoss(
      iterationLoss, iterationCustomLosses, iterationCustomLoss);
    printLoss(iterationTotalLoss, iterationTestLoss, iterationCustomLoss,
              iteration, iteration,
              t_alg.getElapsedTime() - t_allVectorCopy_time_,
              debug::Priority::VERBOSE);

    // --- Update best parameters
    bool isBest = false;
    if(not trainTestSplitted)
      isBest = isBestLoss(iterationTotalLoss, minLoss, cptBlocked);
    else {
      // TODO generalize these lines when accuracy is not the metric computed or
      // evaluated
      if(minCustomLoss.empty())
        isBest = true;
      else {
        float minusAcc = -iterationCustomLoss[1];
        float minMinusAcc = -minCustomLoss[1];
        isBest = isBestLoss(minusAcc, minMinusAcc, cptBlocked);
      }
    }
    if(isBest) {
      Timer t_copy;
      bestIteration = iteration;
      copyParams(bestOrigins, bestOriginsPrime, bestVSTensor, bestVSPrimeTensor,
                 bestOrigins2, bestOrigins2Prime, bestVS2Tensor,
                 bestVS2PrimeTensor, allAlphas_, bestAlphasInit, true);
      copyCustomParams(true);
      copyParams(recs_, bestRecs);
      copyParams(recs2_, bestRecs2);
      t_allVectorCopy_time_ += t_copy.getElapsedTime();
      bestTime = t_alg.getElapsedTime() - t_allVectorCopy_time_;
      minCustomLoss = iterationCustomLoss;
      if(trainTestSplitted) {
        minLoss = iterationTotalLoss;
        minTestLoss = iterationTestLoss;
      }
      printLoss(minLoss, minTestLoss, minCustomLoss, bestIteration, iteration,
                bestTime, debug::Priority::DETAIL);
    }

    // --- Convergence
    windowLoss += iterationTotalLoss;
    if((iteration + 1) % convWinSize == 0) {
      windowLoss /= convWinSize;
      converged = convergenceStep(windowLoss, oldLoss, minLoss, cptBlocked);
      windowLoss = 0;
      if(converged) {
        ++noConverged;
      } else
        noConverged = 0;
      converged = noConverged >= noConvergedToGet;
      if(converged and iteration < minIteration_)
        printMsg("convergence is detected but iteration < minIteration_",
                 debug::Priority::DETAIL);
      if(iteration < minIteration_)
        converged = false;
      if(converged)
        break;
    }

    // --- Print
    if(iteration % iterationGap_ == 0) {
      printMsg("Assignment", 1, assignmentTime, threadNumber_);
      printMsg("Loss", 1, lossTime, threadNumber_);
      printMsg("Update", 1, updateTime, threadNumber_);
      printMsg("Projection", 1, projectionTime, threadNumber_);
      assignmentTime = 0.0;
      lossTime = 0.0;
      updateTime = 0.0;
      projectionTime = 0.0;
      float loss = torch::tensor(gapLosses).mean().item<float>();
      gapLosses.clear();
      float testLoss = torch::tensor(gapTestLosses).mean().item<float>();
      gapTestLosses.clear();
      if(trainTestSplitted) {
        std::stringstream ss;
        ss << "Test Loss = " << testLoss;
        printMsg(ss.str());
      }
      printGapLoss(loss, gapCustomLosses);

      // Verify grad and big values (testing)
      for(unsigned int l = 0; l < noLayers_; ++l) {
        std::stringstream ss;
        if(originsNoZeroGrad_[l] != 0)
          ss << originsNoZeroGrad_[l] << " originsNoZeroGrad_[" << l << "]"
             << std::endl;
        if(originsPrimeNoZeroGrad_[l] != 0)
          ss << originsPrimeNoZeroGrad_[l] << " originsPrimeNoZeroGrad_[" << l
             << "]" << std::endl;
        if(vSNoZeroGrad_[l] != 0)
          ss << vSNoZeroGrad_[l] << " vSNoZeroGrad_[" << l << "]" << std::endl;
        if(vSPrimeNoZeroGrad_[l] != 0)
          ss << vSPrimeNoZeroGrad_[l] << " vSPrimeNoZeroGrad_[" << l << "]"
             << std::endl;
        originsNoZeroGrad_[l] = 0;
        originsPrimeNoZeroGrad_[l] = 0;
        vSNoZeroGrad_[l] = 0;
        vSPrimeNoZeroGrad_[l] = 0;
        if(useDoubleInput_) {
          if(origins2NoZeroGrad_[l] != 0)
            ss << origins2NoZeroGrad_[l] << " origins2NoZeroGrad_[" << l << "]"
               << std::endl;
          if(origins2PrimeNoZeroGrad_[l] != 0)
            ss << origins2PrimeNoZeroGrad_[l] << " origins2PrimeNoZeroGrad_["
               << l << "]" << std::endl;
          if(vS2NoZeroGrad_[l] != 0)
            ss << vS2NoZeroGrad_[l] << " vS2NoZeroGrad_[" << l << "]"
               << std::endl;
          if(vS2PrimeNoZeroGrad_[l] != 0)
            ss << vS2PrimeNoZeroGrad_[l] << " vS2PrimeNoZeroGrad_[" << l << "]"
               << std::endl;
          origins2NoZeroGrad_[l] = 0;
          origins2PrimeNoZeroGrad_[l] = 0;
          vS2NoZeroGrad_[l] = 0;
          vS2PrimeNoZeroGrad_[l] = 0;
        }
        if(isTreeHasBigValues(
             layers_[l].getOrigin().mTree, bigValuesThreshold_))
          ss << "origins_[" << l << "] has big values!" << std::endl;
        if(isTreeHasBigValues(
             layers_[l].getOriginPrime().mTree, bigValuesThreshold_))
          ss << "originsPrime_[" << l << "] has big values!" << std::endl;
        if(ss.rdbuf()->in_avail() != 0)
          printMsg(ss.str(), debug::Priority::DETAIL);
      }
    }

    ++iteration;
    if(maxIteration_ != 0 and iteration >= maxIteration_) {
      printMsg("iteration >= maxIteration_", debug::Priority::DETAIL);
      break;
    }
  }
  printMsg(debug::Separator::L2);
  printLoss(
    minLoss, minTestLoss, minCustomLoss, bestIteration, iteration, bestTime);
  printMsg(debug::Separator::L2);
  bestLoss_ = (trainTestSplitted ? minTestLoss : minLoss);

  Timer t_copy;
  copyParams(bestOrigins, bestOriginsPrime, bestVSTensor, bestVSPrimeTensor,
             bestOrigins2, bestOrigins2Prime, bestVS2Tensor, bestVS2PrimeTensor,
             bestAlphasInit, allAlphas_, false);
  copyCustomParams(false);
  copyParams(bestRecs, recs_);
  copyParams(bestRecs2, recs2_);
  t_allVectorCopy_time_ += t_copy.getElapsedTime();
  printMsg("Copy time", 1, t_allVectorCopy_time_, threadNumber_);
}

//  ---------------------------------------------------------------------------
//  --- End Functions
//  ---------------------------------------------------------------------------
void ttk::MergeTreeNeuralNetwork::computeTrackingInformation(
  unsigned int endLayer) {
  originsMatchings_.resize(endLayer);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
  for(unsigned int l = 0; l < endLayer; ++l) {
    auto &tree1
      = (l == 0 ? layers_[0].getOrigin() : layers_[l - 1].getOriginPrime());
    auto &tree2
      = (l == 0 ? layers_[0].getOriginPrime() : layers_[l].getOriginPrime());
    bool isCalled = true;
    float distance;
    computeOneDistance<float>(tree1.mTree, tree2.mTree, originsMatchings_[l],
                              distance, isCalled, useDoubleInput_);
  }

  // Data matchings
  ++endLayer;
  dataMatchings_.resize(endLayer);
  for(unsigned int l = 0; l < endLayer; ++l) {
    dataMatchings_[l].resize(recs_.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < recs_.size(); ++i) {
      bool isCalled = true;
      float distance;
      auto &origin
        = (l == 0 ? layers_[0].getOrigin() : layers_[l - 1].getOriginPrime());
      computeOneDistance<float>(origin.mTree, recs_[i][l].mTree,
                                dataMatchings_[l][i], distance, isCalled,
                                useDoubleInput_);
    }
  }

  // Reconst matchings
  reconstMatchings_.resize(recs_.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
  for(unsigned int i = 0; i < recs_.size(); ++i) {
    bool isCalled = true;
    float distance;
    auto l = recs_[i].size() - 1;
    computeOneDistance<float>(recs_[i][0].mTree, recs_[i][l].mTree,
                              reconstMatchings_[i], distance, isCalled,
                              useDoubleInput_);
  }
}

void ttk::MergeTreeNeuralNetwork::computeCorrelationMatrix(
  std::vector<ftm::MergeTree<float>> &trees, unsigned int layer) {
  std::vector<std::vector<double>> allTs;
  auto noGeod = allAlphas_[0][layer].sizes()[0];
  allTs.resize(noGeod);
  for(unsigned int i = 0; i < noGeod; ++i) {
    allTs[i].resize(allAlphas_.size());
    for(unsigned int j = 0; j < allAlphas_.size(); ++j)
      allTs[i][j] = allAlphas_[j][layer][i].item<double>();
  }
  computeBranchesCorrelationMatrix(
    layers_[0].getOrigin().mTree, trees, dataMatchings_[0], allTs,
    branchesCorrelationMatrix_, persCorrelationMatrix_);
}

void ttk::MergeTreeNeuralNetwork::createScaledAlphas(
  std::vector<std::vector<torch::Tensor>> &alphas,
  std::vector<std::vector<torch::Tensor>> &scaledAlphas) {
  scaledAlphas.clear();
  scaledAlphas.resize(
    alphas.size(), std::vector<torch::Tensor>(alphas[0].size()));
  for(unsigned int l = 0; l < alphas[0].size(); ++l) {
    torch::Tensor scale = layers_[l].getVSTensor().pow(2).sum(0).sqrt();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < alphas.size(); ++i) {
      scaledAlphas[i][l] = alphas[i][l] * scale.reshape({-1, 1});
    }
  }
}

void ttk::MergeTreeNeuralNetwork::createScaledAlphas() {
  createScaledAlphas(allAlphas_, allScaledAlphas_);
}

void ttk::MergeTreeNeuralNetwork::createActivatedAlphas() {
  allActAlphas_ = allAlphas_;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
  for(unsigned int i = 0; i < allActAlphas_.size(); ++i)
    for(unsigned int j = 0; j < allActAlphas_[i].size(); ++j)
      allActAlphas_[i][j] = activation(allActAlphas_[i][j]);
  createScaledAlphas(allActAlphas_, allActScaledAlphas_);
}

//  ---------------------------------------------------------------------------
//  --- Utils
//  ---------------------------------------------------------------------------
void ttk::MergeTreeNeuralNetwork::copyParams(
  std::vector<mtu::TorchMergeTree<float>> &origins,
  std::vector<mtu::TorchMergeTree<float>> &originsPrime,
  std::vector<torch::Tensor> &vS,
  std::vector<torch::Tensor> &vSPrime,
  std::vector<mtu::TorchMergeTree<float>> &origins2,
  std::vector<mtu::TorchMergeTree<float>> &origins2Prime,
  std::vector<torch::Tensor> &vS2,
  std::vector<torch::Tensor> &vS2Prime,
  std::vector<std::vector<torch::Tensor>> &srcAlphas,
  std::vector<std::vector<torch::Tensor>> &dstAlphas,
  bool get) {
  dstAlphas.resize(srcAlphas.size(), std::vector<torch::Tensor>(noLayers_));
  if(get) {
    origins.resize(noLayers_);
    originsPrime.resize(noLayers_);
    vS.resize(noLayers_);
    vSPrime.resize(noLayers_);
    if(useDoubleInput_) {
      origins2.resize(noLayers_);
      origins2Prime.resize(noLayers_);
      vS2.resize(noLayers_);
      vS2Prime.resize(noLayers_);
    }
  }
  for(unsigned int l = 0; l < noLayers_; ++l) {
    layers_[l].copyParams(origins[l], originsPrime[l], vS[l], vSPrime[l],
                          origins2[l], origins2Prime[l], vS2[l], vS2Prime[l],
                          get);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < srcAlphas.size(); ++i)
      mtu::copyTensor(srcAlphas[i][l], dstAlphas[i][l]);
  }
}

void ttk::MergeTreeNeuralNetwork::copyParams(
  std::vector<std::vector<mtu::TorchMergeTree<float>>> &src,
  std::vector<std::vector<mtu::TorchMergeTree<float>>> &dst) {
  dst.resize(src.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
  for(unsigned int i = 0; i < src.size(); ++i) {
    dst[i].resize(src[i].size());
    for(unsigned int j = 0; j < src[i].size(); ++j)
      mtu::copyTorchMergeTree(src[i][j], dst[i][j]);
  }
}

void ttk::MergeTreeNeuralNetwork::getAlphasTensor(
  std::vector<std::vector<torch::Tensor>> &alphas,
  std::vector<unsigned int> &indexes,
  std::vector<bool> &toGet,
  unsigned int layerIndex,
  torch::Tensor &alphasOut) {
  unsigned int beg = 0;
  while(not toGet[indexes[beg]])
    ++beg;
  alphasOut = alphas[indexes[beg]][layerIndex].transpose(0, 1);
  for(unsigned int ind = beg + 1; ind < indexes.size(); ++ind) {
    if(not toGet[indexes[ind]])
      continue;
    alphasOut = torch::cat(
      {alphasOut, alphas[indexes[ind]][layerIndex].transpose(0, 1)});
  }
}

void ttk::MergeTreeNeuralNetwork::getAlphasTensor(
  std::vector<std::vector<torch::Tensor>> &alphas,
  std::vector<unsigned int> &indexes,
  unsigned int layerIndex,
  torch::Tensor &alphasOut) {
  std::vector<bool> toGet(indexes.size(), true);
  getAlphasTensor(alphas, indexes, toGet, layerIndex, alphasOut);
}

void ttk::MergeTreeNeuralNetwork::getAlphasTensor(
  std::vector<std::vector<torch::Tensor>> &alphas,
  unsigned int layerIndex,
  torch::Tensor &alphasOut) {
  std::vector<unsigned int> indexes(alphas.size());
  std::iota(indexes.begin(), indexes.end(), 0);
  getAlphasTensor(alphas, indexes, layerIndex, alphasOut);
}

//  ---------------------------------------------------------------------------
//  --- Testing
//  ---------------------------------------------------------------------------
void ttk::MergeTreeNeuralNetwork::checkZeroGrad(unsigned int l,
                                                bool checkOutputBasis) {
  if(not layers_[l].getOrigin().tensor.grad().defined()
     or not layers_[l].getOrigin().tensor.grad().count_nonzero().is_nonzero())
    ++originsNoZeroGrad_[l];
  if(not layers_[l].getVSTensor().grad().defined()
     or not layers_[l].getVSTensor().grad().count_nonzero().is_nonzero())
    ++vSNoZeroGrad_[l];
  if(checkOutputBasis) {
    if(not layers_[l].getOriginPrime().tensor.grad().defined()
       or not layers_[l]
                .getOriginPrime()
                .tensor.grad()
                .count_nonzero()
                .is_nonzero())
      ++originsPrimeNoZeroGrad_[l];
    if(not layers_[l].getVSPrimeTensor().grad().defined()
       or not layers_[l].getVSPrimeTensor().grad().count_nonzero().is_nonzero())
      ++vSPrimeNoZeroGrad_[l];
  }
  if(useDoubleInput_) {
    if(not layers_[l].getOrigin2().tensor.grad().defined()
       or not layers_[l]
                .getOrigin2()
                .tensor.grad()
                .count_nonzero()
                .is_nonzero())
      ++origins2NoZeroGrad_[l];
    if(not layers_[l].getVS2Tensor().grad().defined()
       or not layers_[l].getVS2Tensor().grad().count_nonzero().is_nonzero())
      ++vS2NoZeroGrad_[l];
    if(checkOutputBasis) {
      if(not layers_[l].getOrigin2Prime().tensor.grad().defined()
         or not layers_[l]
                  .getOrigin2Prime()
                  .tensor.grad()
                  .count_nonzero()
                  .is_nonzero())
        ++origins2PrimeNoZeroGrad_[l];
      if(not layers_[l].getVS2PrimeTensor().grad().defined()
         or not layers_[l]
                  .getVS2PrimeTensor()
                  .grad()
                  .count_nonzero()
                  .is_nonzero())
        ++vS2PrimeNoZeroGrad_[l];
    }
  }
}

bool ttk::MergeTreeNeuralNetwork::isTreeHasBigValues(
  const ftm::MergeTree<float> &mTree, float threshold) {
  bool found = false;
  for(unsigned int n = 0; n < mTree.tree.getNumberOfNodes(); ++n) {
    if(mTree.tree.isNodeAlone(n))
      continue;
    auto birthDeath = mTree.tree.template getBirthDeath<float>(n);
    if(std::abs(std::get<0>(birthDeath)) > threshold
       or std::abs(std::get<1>(birthDeath)) > threshold) {
      found = true;
      break;
    }
  }
  return found;
}
#endif

//  ---------------------------------------------------------------------------
//  --- Main Functions
//  ---------------------------------------------------------------------------
void ttk::MergeTreeNeuralNetwork::execute(
  std::vector<ftm::MergeTree<float>> &trees,
  std::vector<ftm::MergeTree<float>> &trees2) {
#ifndef TTK_ENABLE_TORCH
  TTK_FORCE_USE(trees);
  TTK_FORCE_USE(trees2);
  printErr("This module requires Torch.");
#else
#ifdef TTK_ENABLE_OPENMP
  int ompNested = omp_get_nested();
  omp_set_nested(1);
#endif
  // makeExponentialExample(trees, trees2);

  // --- Preprocessing
  Timer t_preprocess;
  preprocessingTrees<float>(trees, treesNodeCorr_);
  if(trees2.size() != 0)
    preprocessingTrees<float>(trees2, trees2NodeCorr_);
  printMsg("Preprocessing", 1, t_preprocess.getElapsedTime(), threadNumber_);
  useDoubleInput_ = (trees2.size() != 0);

  // --- Fit neural network
  Timer t_total;
  fit(trees, trees2);
  auto totalTime = t_total.getElapsedTime() - t_allVectorCopy_time_;
  printMsg(debug::Separator::L1);
  printMsg("Total time", 1, totalTime, threadNumber_);

  // --- End functions
  Timer t_end;
  createScaledAlphas();
  createActivatedAlphas();
  executeEndFunction(trees, trees2);
  printMsg("End functions", 1, t_end.getElapsedTime(), threadNumber_);

  // --- Postprocessing
  if(createOutput_) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < trees.size(); ++i)
      postprocessingPipeline<float>(&(trees[i].tree));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < trees2.size(); ++i)
      postprocessingPipeline<float>(&(trees2[i].tree));

    originsCopy_.resize(layers_.size());
    originsPrimeCopy_.resize(layers_.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int l = 0; l < layers_.size(); ++l) {
      mtu::copyTorchMergeTree<float>(layers_[l].getOrigin(), originsCopy_[l]);
      mtu::copyTorchMergeTree<float>(
        layers_[l].getOriginPrime(), originsPrimeCopy_[l]);
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int l = 0; l < originsCopy_.size(); ++l) {
      fillMergeTreeStructure(originsCopy_[l]);
      postprocessingPipeline<float>(&(originsCopy_[l].mTree.tree));
      fillMergeTreeStructure(originsPrimeCopy_[l]);
      postprocessingPipeline<float>(&(originsPrimeCopy_[l].mTree.tree));
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < recs_.size(); ++i) {
      for(unsigned int j = 0; j < recs_[i].size(); ++j) {
        fixTreePrecisionScalars(recs_[i][j].mTree);
        postprocessingPipeline<float>(&(recs_[i][j].mTree.tree));
      }
    }
  }

  if(not isPersistenceDiagram_) {
    for(unsigned int l = 0; l < originsMatchings_.size(); ++l) {
      auto &tree1 = (l == 0 ? originsCopy_[0] : originsPrimeCopy_[l - 1]);
      auto &tree2 = (l == 0 ? originsPrimeCopy_[0] : originsPrimeCopy_[l]);
      convertBranchDecompositionMatching<float>(
        &(tree1.mTree.tree), &(tree2.mTree.tree), originsMatchings_[l]);
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < recs_.size(); ++i) {
      for(unsigned int l = 0; l < dataMatchings_.size(); ++l) {
        auto &origin = (l == 0 ? originsCopy_[0] : originsPrimeCopy_[l - 1]);
        convertBranchDecompositionMatching<float>(&(origin.mTree.tree),
                                                  &(recs_[i][l].mTree.tree),
                                                  dataMatchings_[l][i]);
      }
    }
    for(unsigned int i = 0; i < reconstMatchings_.size(); ++i) {
      auto l = recs_[i].size() - 1;
      convertBranchDecompositionMatching<float>(&(recs_[i][0].mTree.tree),
                                                &(recs_[i][l].mTree.tree),
                                                reconstMatchings_[i]);
    }
  }
#ifdef TTK_ENABLE_OPENMP
  omp_set_nested(ompNested);
#endif
#endif
}
