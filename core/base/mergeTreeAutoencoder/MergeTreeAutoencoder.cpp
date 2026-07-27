#include <MergeTreeAutoencoder.h>
#include <cmath>

#ifdef TTK_ENABLE_TORCH
using namespace torch::indexing;
#endif

ttk::MergeTreeAutoencoder::MergeTreeAutoencoder() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("MergeTreeAutoencoder");
}

#ifdef TTK_ENABLE_TORCH
void ttk::MergeTreeAutoencoder::initClusteringLossParameters() {
  unsigned int l = getLatentLayerIndex();
  unsigned int noCentroids
    = std::set<unsigned int>(clusterAsgn_.begin(), clusterAsgn_.end()).size();
  latentCentroids_.resize(noCentroids);
  for(unsigned int c = 0; c < noCentroids; ++c) {
    unsigned int firstIndex = std::numeric_limits<unsigned int>::max();
    for(unsigned int i = 0; i < clusterAsgn_.size(); ++i) {
      if(clusterAsgn_[i] == c) {
        firstIndex = i;
        break;
      }
    }
    if(firstIndex >= allAlphas_.size()) {
      printWrn("no data found for cluster " + std::to_string(c));
      // TODO init random centroid
    }
    latentCentroids_[c] = allAlphas_[firstIndex][l].detach().clone();
    float noData = 1;
    for(unsigned int i = 0; i < allAlphas_.size(); ++i) {
      if(i == firstIndex)
        continue;
      if(clusterAsgn_[i] == c) {
        latentCentroids_[c] += allAlphas_[i][l];
        ++noData;
      }
    }
    latentCentroids_[c] /= torch::tensor(noData);
    latentCentroids_[c] = latentCentroids_[c].detach();
    latentCentroids_[c].requires_grad_(true);
  }
}

bool ttk::MergeTreeAutoencoder::initResetOutputBasis(
  unsigned int l,
  unsigned int layerNoAxes,
  double layerOriginPrimeSizePercent,
  std::vector<mtu::TorchMergeTree<float>> &trees,
  std::vector<mtu::TorchMergeTree<float>> &trees2,
  std::vector<bool> &isTrain) {
  printMsg("Reset output basis", debug::Priority::DETAIL);
  if((noLayers_ == 2 and l == 1) or noLayers_ == 1) {
    initOutputBasisSpecialCase(l, layerNoAxes, trees, trees2);
  } else if(l < (unsigned int)(noLayers_ / 2)) {
    initOutputBasis(l, layerOriginPrimeSizePercent, trees, trees2, isTrain);
  } else {
    printErr("recs[i].mTree.tree.getRealNumberOfNodes() == 0");
    std::stringstream ssT;
    ssT << "layer " << l;
    printWrn(ssT.str());
    return true;
  }
  return false;
}

void ttk::MergeTreeAutoencoder::initOutputBasisSpecialCase(
  unsigned int l,
  unsigned int layerNoAxes,
  std::vector<mtu::TorchMergeTree<float>> &trees,
  std::vector<mtu::TorchMergeTree<float>> &trees2) {
  // - Compute Origin
  printMsg("Compute output basis origin", debug::Priority::DETAIL);
  layers_[l].setOriginPrime(layers_[0].getOrigin());
  if(useDoubleInput_)
    layers_[l].setOrigin2Prime(layers_[0].getOrigin2());
  // - Compute vectors
  printMsg("Compute output basis vectors", debug::Priority::DETAIL);
  if(layerNoAxes != layers_[0].getVSTensor().sizes()[1]) {
    // TODO is there a way to avoid copy of merge trees?
    std::vector<ftm::MergeTree<float>> treesToUse, trees2ToUse;
    for(unsigned int i = 0; i < trees.size(); ++i) {
      treesToUse.emplace_back(trees[i].mTree);
      if(useDoubleInput_)
        trees2ToUse.emplace_back(trees2[i].mTree);
    }
    std::vector<torch::Tensor> allAlphasInitT(trees.size());
    layers_[l].initInputBasisVectors(
      trees, trees2, treesToUse, trees2ToUse, layerNoAxes, allAlphasInitT,
      inputToBaryDistances_L0_, baryMatchings_L0_, baryMatchings2_L0_, false);
  } else {
    layers_[l].setVSPrimeTensor(layers_[0].getVSTensor());
    if(useDoubleInput_)
      layers_[l].setVS2PrimeTensor(layers_[0].getVS2Tensor());
  }
}

float ttk::MergeTreeAutoencoder::initParameters(
  std::vector<mtu::TorchMergeTree<float>> &trees,
  std::vector<mtu::TorchMergeTree<float>> &trees2,
  std::vector<bool> &isTrain,
  bool computeError) {
  // ----- Init variables
  // noLayers_ = number of encoder layers + number of decoder layers + the
  // latent layer + the output layer
  noLayers_ = encoderNoLayers_ * 2 + 1 + 1;
  if(encoderNoLayers_ <= -1)
    noLayers_ = 1;
  std::vector<double> layersOriginPrimeSizePercent(noLayers_);
  std::vector<unsigned int> layersNoAxes(noLayers_);
  if(noLayers_ <= 2) {
    layersNoAxes[0] = numberOfAxes_;
    layersOriginPrimeSizePercent[0] = latentSpaceOriginPrimeSizePercent_;
    if(noLayers_ == 2) {
      layersNoAxes[1] = inputNumberOfAxes_;
      layersOriginPrimeSizePercent[1] = barycenterSizeLimitPercent_;
    }
  } else {
    for(unsigned int l = 0; l < noLayers_ / 2; ++l) {
      double alpha = (double)(l) / (noLayers_ / 2 - 1);
      unsigned int noAxes
        = (1 - alpha) * inputNumberOfAxes_ + alpha * numberOfAxes_;
      layersNoAxes[l] = noAxes;
      layersNoAxes[noLayers_ - 1 - l] = noAxes;
      double originPrimeSizePercent
        = (1 - alpha) * inputOriginPrimeSizePercent_
          + alpha * latentSpaceOriginPrimeSizePercent_;
      layersOriginPrimeSizePercent[l] = originPrimeSizePercent;
      layersOriginPrimeSizePercent[noLayers_ - 1 - l] = originPrimeSizePercent;
    }
    if(scaleLayerAfterLatent_)
      layersNoAxes[noLayers_ / 2]
        = (layersNoAxes[noLayers_ / 2 - 1] + layersNoAxes[noLayers_ / 2 + 1])
          / 2.0;
  }

  // ----- Resize parameters
  layers_.resize(noLayers_);
  for(unsigned int l = 0; l < layers_.size(); ++l) {
    initOriginPrimeValuesByCopy_
      = trackingLossWeight_ != 0
        and l < (trackingLossDecoding_ ? noLayers_ : getLatentLayerIndex() + 1);
    initOriginPrimeValuesByCopyRandomness_ = trackingLossInitRandomness_;
    passLayerParameters(layers_[l]);
  }

  // ----- Compute parameters of each layer
  bool fullSymmetricAE = fullSymmetricAE_;

  std::vector<mtu::TorchMergeTree<float>> recs, recs2;
  std::vector<std::vector<torch::Tensor>> allAlphasInit(
    trees.size(), std::vector<torch::Tensor>(noLayers_));
  for(unsigned int l = 0; l < noLayers_; ++l) {
    printMsg(debug::Separator::L2, debug::Priority::DETAIL);
    std::stringstream ss;
    ss << "Init Layer " << l;
    printMsg(ss.str(), debug::Priority::DETAIL);

    // --- Init Input Basis
    if(l < (unsigned int)(noLayers_ / 2) or not fullSymmetricAE
       or (noLayers_ <= 2 and not fullSymmetricAE)) {
      auto &treesToUse = (l == 0 ? trees : recs);
      auto &trees2ToUse = (l == 0 ? trees2 : recs2);
      initInputBasis(
        l, layersNoAxes[l], treesToUse, trees2ToUse, isTrain, allAlphasInit);
    } else {
      // - Copy output tensors of the opposite layer (full symmetric init)
      printMsg(
        "Copy output tensors of the opposite layer", debug::Priority::DETAIL);
      unsigned int middle = noLayers_ / 2;
      unsigned int l_opp = middle - (l - middle + 1);
      layers_[l].setOrigin(layers_[l_opp].getOriginPrime());
      layers_[l].setVSTensor(layers_[l_opp].getVSPrimeTensor());
      if(trees2.size() != 0) {
        if(fullSymmetricAE) {
          layers_[l].setOrigin2(layers_[l_opp].getOrigin2Prime());
          layers_[l].setVS2Tensor(layers_[l_opp].getVS2PrimeTensor());
        }
      }
      for(unsigned int i = 0; i < trees.size(); ++i)
        allAlphasInit[i][l] = allAlphasInit[i][l_opp];
    }

    // --- Init Output Basis
    if((noLayers_ == 2 and l == 1) or noLayers_ == 1) {
      // -- Special case
      initOutputBasisSpecialCase(l, layersNoAxes[l], trees, trees2);
    } else if(l < (unsigned int)(noLayers_ / 2)) {
      initOutputBasis(
        l, layersOriginPrimeSizePercent[l], trees, trees2, isTrain);
    } else {
      // - Copy input tensors of the opposite layer (symmetric init)
      printMsg(
        "Copy input tensors of the opposite layer", debug::Priority::DETAIL);
      unsigned int middle = noLayers_ / 2;
      unsigned int l_opp = middle - (l - middle + 1);
      layers_[l].setOriginPrime(layers_[l_opp].getOrigin());
      if(trees2.size() != 0)
        layers_[l].setOrigin2Prime(layers_[l_opp].getOrigin2());
      if(l == (unsigned int)(noLayers_) / 2 and scaleLayerAfterLatent_) {
        unsigned int dim2
          = (trees2.size() != 0 ? layers_[l].getOrigin2Prime().tensor.sizes()[0]
                                : 0);
        layers_[l].initOutputBasisVectors(
          layers_[l].getOriginPrime().tensor.sizes()[0], dim2);
      } else {
        layers_[l].setVSPrimeTensor(layers_[l_opp].getVSTensor());
        if(trees2.size() != 0)
          layers_[l].setVS2PrimeTensor(layers_[l_opp].getVS2Tensor());
      }
    }

    // --- Get reconstructed
    bool fullReset = initGetReconstructed(
      l, layersNoAxes[l], layersOriginPrimeSizePercent[l], trees, trees2,
      isTrain, recs, recs2, allAlphasInit);
    if(fullReset)
      return std::numeric_limits<float>::max();
  }
  allAlphas_ = allAlphasInit;

  // Init clustering parameters if needed
  if(clusteringLossWeight_ != 0)
    initClusteringLossParameters();

  // Compute error
  float error = 0.0, recLoss = 0.0;
  if(computeError) {
    printMsg("Compute error", debug::Priority::DETAIL);
    std::vector<unsigned int> indexes(trees.size());
    std::iota(indexes.begin(), indexes.end(), 0);
    // TODO forward only if necessary
    unsigned int k = k_;
    std::vector<std::vector<torch::Tensor>> bestAlphas;
    std::vector<std::vector<mtu::TorchMergeTree<float>>> layersOuts,
      layersOuts2;
    std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
      matchings, matchings2;
    bool reset = forwardStep(trees, trees2, indexes, k, allAlphasInit,
                             computeError, recs, recs2, bestAlphas, layersOuts,
                             layersOuts2, matchings, matchings2, recLoss);
    if(reset) {
      printWrn("[initParameters] forwardStep reset");
      return std::numeric_limits<float>::max();
    }
    error = recLoss * reconstructionLossWeight_;
    if(metricLossWeight_ != 0) {
      torch::Tensor metricLoss;
      computeMetricLoss(layersOuts, layersOuts2, allAlphas_, distanceMatrix_,
                        indexes, metricLoss);
      baseRecLoss_ = std::numeric_limits<double>::max();
      metricLoss *= metricLossWeight_
                    * getCustomLossDynamicWeight(recLoss, baseRecLoss_);
      error += metricLoss.item<float>();
    }
    if(clusteringLossWeight_ != 0) {
      torch::Tensor clusteringLoss, asgn;
      computeClusteringLoss(allAlphas_, indexes, clusteringLoss, asgn);
      baseRecLoss_ = std::numeric_limits<double>::max();
      clusteringLoss *= clusteringLossWeight_
                        * getCustomLossDynamicWeight(recLoss, baseRecLoss_);
      error += clusteringLoss.item<float>();
    }
    if(trackingLossWeight_ != 0) {
      torch::Tensor trackingLoss;
      computeTrackingLoss(trackingLoss);
      trackingLoss *= trackingLossWeight_;
      error += trackingLoss.item<float>();
    }
  }
  return error;
}

//  ---------------------------------------------------------------------------
//  --- Backward
//  ---------------------------------------------------------------------------
bool ttk::MergeTreeAutoencoder::backwardStep(
  std::vector<mtu::TorchMergeTree<float>> &trees,
  std::vector<mtu::TorchMergeTree<float>> &outs,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &matchings,
  std::vector<mtu::TorchMergeTree<float>> &trees2,
  std::vector<mtu::TorchMergeTree<float>> &outs2,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &matchings2,
  std::vector<std::vector<torch::Tensor>> &ttkNotUsed(alphas),
  torch::optim::Optimizer &optimizer,
  std::vector<unsigned int> &indexes,
  std::vector<bool> &ttkNotUsed(isTrain),
  std::vector<torch::Tensor> &torchCustomLoss) {
  double totalLoss = 0;
  bool retainGraph = (metricLossWeight_ != 0 or clusteringLossWeight_ != 0
                      or trackingLossWeight_ != 0);
  if(reconstructionLossWeight_ != 0
     or (customLossDynamicWeight_ and retainGraph)) {
    std::vector<torch::Tensor> outTensors(indexes.size()),
      reorderedTensors(indexes.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int ind = 0; ind < indexes.size(); ++ind) {
      unsigned int i = indexes[ind];
      torch::Tensor reorderedTensor;
      dataReorderingGivenMatching(
        outs[i], trees[i], matchings[i], reorderedTensor);
      auto outTensor = outs[i].tensor;
      if(useDoubleInput_) {
        torch::Tensor reorderedTensor2;
        dataReorderingGivenMatching(
          outs2[i], trees2[i], matchings2[i], reorderedTensor2);
        outTensor = torch::cat({outTensor, outs2[i].tensor});
        reorderedTensor = torch::cat({reorderedTensor, reorderedTensor2});
      }
      outTensors[ind] = outTensor;
      reorderedTensors[ind] = reorderedTensor;
    }
    for(unsigned int ind = 0; ind < indexes.size(); ++ind) {
      auto loss = torch::nn::functional::mse_loss(
        outTensors[ind], reorderedTensors[ind]);
      // Same as next loss with a factor of 1 / n where n is the number of nodes
      // in the output
      // auto loss = (outTensors[ind] - reorderedTensors[ind]).pow(2).sum();
      totalLoss += loss.item<float>();
      loss *= reconstructionLossWeight_;
      loss.backward({}, retainGraph);
    }
  }
  if(metricLossWeight_ != 0) {
    bool retainGraphMetricLoss
      = (clusteringLossWeight_ != 0 or trackingLossWeight_ != 0);
    torchCustomLoss[0] *= metricLossWeight_
                          * getCustomLossDynamicWeight(
                            totalLoss / indexes.size(), baseRecLoss2_);
    torchCustomLoss[0].backward({}, retainGraphMetricLoss);
  }
  if(clusteringLossWeight_ != 0) {
    bool retainGraphClusteringLoss = (trackingLossWeight_ != 0);
    torchCustomLoss[1] *= clusteringLossWeight_
                          * getCustomLossDynamicWeight(
                            totalLoss / indexes.size(), baseRecLoss2_);
    torchCustomLoss[1].backward({}, retainGraphClusteringLoss);
  }
  if(trackingLossWeight_ != 0) {
    torchCustomLoss[2] *= trackingLossWeight_;
    torchCustomLoss[2].backward();
  }

  for(unsigned int l = 0; l < noLayers_; ++l)
    checkZeroGrad(l);

  optimizer.step();
  optimizer.zero_grad();
  return false;
}

//  ---------------------------------------------------------------------------
//  --- Convergence
//  ---------------------------------------------------------------------------
float ttk::MergeTreeAutoencoder::computeOneLoss(
  mtu::TorchMergeTree<float> &tree,
  mtu::TorchMergeTree<float> &out,
  mtu::TorchMergeTree<float> &tree2,
  mtu::TorchMergeTree<float> &out2,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching2,
  std::vector<torch::Tensor> &ttkNotUsed(alphas),
  unsigned int ttkNotUsed(treeIndex)) {
  float loss = 0;
  bool isCalled = true;
  float distance;
  computeOneDistance<float>(
    out.mTree, tree.mTree, matching, distance, isCalled, useDoubleInput_);
  if(useDoubleInput_) {
    float distance2;
    computeOneDistance<float>(out2.mTree, tree2.mTree, matching2, distance2,
                              isCalled, useDoubleInput_, false);
    distance = mixDistances<float>(distance, distance2);
  }
  loss += distance * distance;
  return loss;
}

//  ---------------------------------------------------------------------------
//  --- Main Functions
//  ---------------------------------------------------------------------------
void ttk::MergeTreeAutoencoder::customInit(
  std::vector<mtu::TorchMergeTree<float>> &torchTrees,
  std::vector<mtu::TorchMergeTree<float>> &torchTrees2) {
  baseRecLoss_ = std::numeric_limits<double>::max();
  baseRecLoss2_ = std::numeric_limits<double>::max();
  // ----- Init Metric Loss
  if(metricLossWeight_ != 0)
    getDistanceMatrix(torchTrees, torchTrees2, distanceMatrix_);
}

void ttk::MergeTreeAutoencoder::addCustomParameters(
  std::vector<torch::Tensor> &parameters) {
  if(clusteringLossWeight_ != 0)
    for(unsigned int i = 0; i < latentCentroids_.size(); ++i)
      parameters.emplace_back(latentCentroids_[i]);
}

void ttk::MergeTreeAutoencoder::computeCustomLosses(
  std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts,
  std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts2,
  std::vector<std::vector<torch::Tensor>> &bestAlphas,
  std::vector<unsigned int> &indexes,
  std::vector<bool> &ttkNotUsed(isTrain),
  unsigned int ttkNotUsed(iteration),
  std::vector<std::vector<float>> &gapCustomLosses,
  std::vector<std::vector<float>> &iterationCustomLosses,
  std::vector<torch::Tensor> &torchCustomLoss) {
  if(gapCustomLosses.empty())
    gapCustomLosses.resize(3);
  if(iterationCustomLosses.empty())
    iterationCustomLosses.resize(3);
  torchCustomLoss.resize(3);
  // - Metric Loss
  if(metricLossWeight_ != 0) {
    computeMetricLoss(layersOuts, layersOuts2, bestAlphas, distanceMatrix_,
                      indexes, torchCustomLoss[0]);
    float metricLossF = torchCustomLoss[0].item<float>();
    gapCustomLosses[0].emplace_back(metricLossF);
    iterationCustomLosses[0].emplace_back(metricLossF);
  }
  // - Clustering Loss
  if(clusteringLossWeight_ != 0) {
    torch::Tensor asgn;
    computeClusteringLoss(bestAlphas, indexes, torchCustomLoss[1], asgn);
    float clusteringLossF = torchCustomLoss[1].item<float>();
    gapCustomLosses[1].emplace_back(clusteringLossF);
    iterationCustomLosses[1].emplace_back(clusteringLossF);
  }
  // - Tracking Loss
  if(trackingLossWeight_ != 0) {
    computeTrackingLoss(torchCustomLoss[2]);
    float trackingLossF = torchCustomLoss[2].item<float>();
    gapCustomLosses[2].emplace_back(trackingLossF);
    iterationCustomLosses[2].emplace_back(trackingLossF);
  }
}

float ttk::MergeTreeAutoencoder::computeIterationTotalLoss(
  float iterationLoss,
  std::vector<std::vector<float>> &iterationCustomLosses,
  std::vector<float> &iterationCustomLoss) {
  iterationCustomLoss.emplace_back(iterationLoss);
  float iterationTotalLoss = reconstructionLossWeight_ * iterationLoss;
  // Metric
  float iterationMetricLoss = 0;
  if(metricLossWeight_ != 0) {
    iterationMetricLoss
      = torch::tensor(iterationCustomLosses[0]).mean().item<float>();
    iterationTotalLoss
      += metricLossWeight_
         * getCustomLossDynamicWeight(iterationLoss, baseRecLoss_)
         * iterationMetricLoss;
  }
  iterationCustomLoss.emplace_back(iterationMetricLoss);
  // Clustering
  float iterationClusteringLoss = 0;
  if(clusteringLossWeight_ != 0) {
    iterationClusteringLoss
      = torch::tensor(iterationCustomLosses[1]).mean().item<float>();
    iterationTotalLoss
      += clusteringLossWeight_
         * getCustomLossDynamicWeight(iterationLoss, baseRecLoss_)
         * iterationClusteringLoss;
  }
  iterationCustomLoss.emplace_back(iterationClusteringLoss);
  // Tracking
  float iterationTrackingLoss = 0;
  if(trackingLossWeight_ != 0) {
    iterationTrackingLoss
      = torch::tensor(iterationCustomLosses[2]).mean().item<float>();
    iterationTotalLoss += trackingLossWeight_ * iterationTrackingLoss;
  }
  iterationCustomLoss.emplace_back(iterationTrackingLoss);
  return iterationTotalLoss;
}

void ttk::MergeTreeAutoencoder::printCustomLosses(
  std::vector<float> &customLoss,
  std::stringstream &prefix,
  const debug::Priority &priority) {
  if(priority != debug::Priority::VERBOSE)
    prefix.str("");
  std::stringstream ssBestLoss;
  if(metricLossWeight_ != 0 or clusteringLossWeight_ != 0
     or trackingLossWeight_ != 0) {
    ssBestLoss.str("");
    ssBestLoss << "- Rec. " << prefix.str() << "loss   = " << customLoss[0];
    printMsg(ssBestLoss.str(), priority);
  }
  if(metricLossWeight_ != 0) {
    ssBestLoss.str("");
    ssBestLoss << "- Metric " << prefix.str() << "loss = " << customLoss[1];
    printMsg(ssBestLoss.str(), priority);
  }
  if(clusteringLossWeight_ != 0) {
    ssBestLoss.str("");
    ssBestLoss << "- Clust. " << prefix.str() << "loss = " << customLoss[2];
    printMsg(ssBestLoss.str(), priority);
  }
  if(trackingLossWeight_ != 0) {
    ssBestLoss.str("");
    ssBestLoss << "- Track. " << prefix.str() << "loss = " << customLoss[3];
    printMsg(ssBestLoss.str(), priority);
  }
}

void ttk::MergeTreeAutoencoder::printGapLoss(
  float loss, std::vector<std::vector<float>> &gapCustomLosses) {
  std::stringstream ss;
  ss << "Rec. loss   = " << loss;
  printMsg(ss.str());
  if(metricLossWeight_ != 0) {
    float metricLoss = torch::tensor(gapCustomLosses[0]).mean().item<float>();
    gapCustomLosses[0].clear();
    ss.str("");
    ss << "Metric loss = " << metricLoss;
    printMsg(ss.str());
  }
  if(clusteringLossWeight_ != 0) {
    float clusteringLoss
      = torch::tensor(gapCustomLosses[1]).mean().item<float>();
    gapCustomLosses[1].clear();
    ss.str("");
    ss << "Clust. loss = " << clusteringLoss;
    printMsg(ss.str());
  }
  if(trackingLossWeight_ != 0) {
    float trackingLoss = torch::tensor(gapCustomLosses[2]).mean().item<float>();
    gapCustomLosses[2].clear();
    ss.str("");
    ss << "Track. loss = " << trackingLoss;
    printMsg(ss.str());
  }
}

//  ---------------------------------------------------------------------------
//  --- Custom Losses
//  ---------------------------------------------------------------------------
double ttk::MergeTreeAutoencoder::getCustomLossDynamicWeight(double recLoss,
                                                             double &baseLoss) {
  baseLoss = std::min(recLoss, baseLoss);
  if(customLossDynamicWeight_)
    return baseLoss;
  else
    return 1.0;
}

void ttk::MergeTreeAutoencoder::computeMetricLoss(
  std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts,
  std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts2,
  std::vector<std::vector<torch::Tensor>> alphas,
  std::vector<std::vector<float>> &baseDistanceMatrix,
  std::vector<unsigned int> &indexes,
  torch::Tensor &metricLoss) {
  auto layerIndex = getLatentLayerIndex();
  std::vector<std::vector<torch::Tensor>> losses(
    layersOuts.size(), std::vector<torch::Tensor>(layersOuts.size()));

  std::vector<mtu::TorchMergeTree<float> *> trees, trees2;
  for(unsigned int ind = 0; ind < indexes.size(); ++ind) {
    unsigned int i = indexes[ind];
    trees.emplace_back(&(layersOuts[i][layerIndex]));
    if(useDoubleInput_)
      trees2.emplace_back(&(layersOuts2[i][layerIndex]));
  }

  std::vector<std::vector<torch::Tensor>> outDistMat;
  torch::Tensor coefDistMat;
  if(customLossSpace_) {
    getDifferentiableDistanceMatrix(trees, trees2, outDistMat);
  } else {
    std::vector<std::vector<torch::Tensor>> scaledAlphas;
    createScaledAlphas(alphas, scaledAlphas);
    torch::Tensor latentAlphas;
    getAlphasTensor(scaledAlphas, indexes, layerIndex, latentAlphas);
    if(customLossActivate_)
      latentAlphas = activation(latentAlphas);
    coefDistMat = torch::cdist(latentAlphas, latentAlphas).pow(2);
  }

  torch::Tensor maxLoss = torch::tensor(0);
  metricLoss = torch::tensor(0);
  float div = 0;
  for(unsigned int ind = 0; ind < indexes.size(); ++ind) {
    unsigned int i = indexes[ind];
    for(unsigned int ind2 = ind + 1; ind2 < indexes.size(); ++ind2) {
      unsigned int j = indexes[ind2];
      torch::Tensor loss;
      torch::Tensor toCompare
        = (customLossSpace_ ? outDistMat[i][j] : coefDistMat[ind][ind2]);
      loss = torch::nn::MSELoss()(
        torch::tensor(baseDistanceMatrix[i][j]), toCompare);
      metricLoss = metricLoss + loss;
      maxLoss = torch::max(loss, maxLoss);
      ++div;
    }
  }
  metricLoss = metricLoss / torch::tensor(div);
  if(normalizeMetricLoss_)
    metricLoss /= maxLoss;
}

void ttk::MergeTreeAutoencoder::computeClusteringLoss(
  std::vector<std::vector<torch::Tensor>> &alphas,
  std::vector<unsigned int> &indexes,
  torch::Tensor &clusteringLoss,
  torch::Tensor &asgn) {
  // Compute distance matrix
  unsigned int layerIndex = getLatentLayerIndex();
  torch::Tensor latentAlphas;
  getAlphasTensor(alphas, indexes, layerIndex, latentAlphas);
  if(customLossActivate_)
    latentAlphas = activation(latentAlphas);
  torch::Tensor centroids = latentCentroids_[0].transpose(0, 1);
  for(unsigned int i = 1; i < latentCentroids_.size(); ++i)
    centroids = torch::cat({centroids, latentCentroids_[i].transpose(0, 1)});
  torch::Tensor dist = torch::cdist(latentAlphas, centroids);

  // Compute softmax and one hot real asgn
  dist = dist * -clusteringLossTemp_;
  asgn = torch::nn::Softmax(1)(dist);
  std::vector<float> clusterAsgn;
  for(unsigned int ind = 0; ind < indexes.size(); ++ind) {
    clusterAsgn.emplace_back(clusterAsgn_[indexes[ind]]);
  }
  torch::Tensor realAsgn = torch::tensor(clusterAsgn).to(torch::kInt64);
  realAsgn
    = torch::nn::functional::one_hot(realAsgn, asgn.sizes()[1]).to(torch::kF32);

  // Compute KL div.
  clusteringLoss = torch::nn::KLDivLoss(
    torch::nn::KLDivLossOptions().reduction(torch::kBatchMean))(asgn, realAsgn);
}

void ttk::MergeTreeAutoencoder::computeTrackingLoss(
  torch::Tensor &trackingLoss) {
  unsigned int latentLayerIndex = getLatentLayerIndex() + 1;
  auto endLayer = (trackingLossDecoding_ ? noLayers_ : latentLayerIndex);
  std::vector<torch::Tensor> losses(endLayer);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
  for(unsigned int l = 0; l < endLayer; ++l) {
    auto &tree1
      = (l == 0 ? layers_[0].getOrigin() : layers_[l - 1].getOriginPrime());
    auto &tree2
      = (l == 0 ? layers_[0].getOriginPrime() : layers_[l].getOriginPrime());
    torch::Tensor tensorDist;
    bool isCalled = true, doSqrt = false;
    getDifferentiableDistance(tree1, tree2, tensorDist, isCalled, doSqrt);
    losses[l] = tensorDist;
  }
  trackingLoss = torch::tensor(0, torch::kFloat32);
  for(unsigned int i = 0; i < losses.size(); ++i)
    trackingLoss += losses[i];
}

//  ---------------------------------------------------------------------------
//  --- End Functions
//  ---------------------------------------------------------------------------
void ttk::MergeTreeAutoencoder::createCustomRecs() {
  if(customAlphas_.empty())
    return;

  bool initByTreesAlphas = not allAlphas_.empty();
  std::vector<torch::Tensor> allTreesAlphas;
  if(initByTreesAlphas) {
    allTreesAlphas.resize(allAlphas_[0].size());
    for(unsigned int l = 0; l < allTreesAlphas.size(); ++l) {
      allTreesAlphas[l] = allAlphas_[0][l].reshape({-1, 1});
      for(unsigned int i = 1; i < allAlphas_.size(); ++i)
        allTreesAlphas[l]
          = torch::cat({allTreesAlphas[l], allAlphas_[i][l]}, 1);
      allTreesAlphas[l] = allTreesAlphas[l].transpose(0, 1);
    }
  }

  unsigned int latLayer = getLatentLayerIndex();
  customRecs_.resize(customAlphas_.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
  for(unsigned int i = 0; i < customAlphas_.size(); ++i) {
    torch::Tensor alphas = torch::tensor(customAlphas_[i]).reshape({-1, 1});

    torch::Tensor alphasWeight;
    if(initByTreesAlphas) {
      auto driver = "gelsd";
      alphasWeight = std::get<0>(torch::linalg_lstsq(
                                   allTreesAlphas[latLayer].transpose(0, 1),
                                   alphas, c10::nullopt, driver))
                       .transpose(0, 1);
    }

    // Reconst latent
    std::vector<mtu::TorchMergeTree<float>> outs, outs2;
    auto noOuts = noLayers_ - latLayer;
    outs.resize(noOuts);
    outs2.resize(noOuts);
    mtu::TorchMergeTree<float> out, out2;
    layers_[latLayer].outputBasisReconstruction(alphas, outs[0], outs2[0]);
    // Decoding
    unsigned int k = 32;
    for(unsigned int l = latLayer + 1; l < noLayers_; ++l) {
      unsigned int noIter = (initByTreesAlphas ? 1 : 32);
      std::vector<torch::Tensor> allAlphasInit(noIter);
      torch::Tensor maxNorm;
      for(unsigned int j = 0; j < allAlphasInit.size(); ++j) {
        allAlphasInit[j]
          = torch::randn({layers_[l].getVSTensor().sizes()[1], 1});
        auto norm = torch::linalg_vector_norm(
          allAlphasInit[j], 2, 0, false, c10::nullopt);
        if(j == 0 or maxNorm.item<float>() < norm.item<float>())
          maxNorm = norm;
      }
      for(unsigned int j = 0; j < allAlphasInit.size(); ++j)
        allAlphasInit[j] /= maxNorm;
      float bestDistance = std::numeric_limits<float>::max();
      auto outIndex = l - latLayer;
      mtu::TorchMergeTree<float> outToUse;
      for(unsigned int j = 0; j < noIter; ++j) {
        torch::Tensor alphasInit, dataAlphas;
        if(initByTreesAlphas) {
          alphasInit
            = torch::matmul(alphasWeight, allTreesAlphas[l]).transpose(0, 1);
        } else {
          alphasInit = allAlphasInit[j];
        }
        float distance;
        layers_[l].forward(outs[outIndex - 1], outs2[outIndex - 1], k,
                           alphasInit, outToUse, outs2[outIndex], dataAlphas,
                           distance);
        if(distance < bestDistance) {
          bestDistance = distance;
          mtu::copyTorchMergeTree<float>(
            outToUse, (l != noLayers_ - 1 ? outs[outIndex] : customRecs_[i]));
        }
      }
    }
  }

  customMatchings_.resize(customRecs_.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
  for(unsigned int i = 0; i < customRecs_.size(); ++i) {
    bool isCalled = true;
    float distance;
    computeOneDistance<float>(layers_[0].getOrigin().mTree,
                              customRecs_[i].mTree, customMatchings_[i],
                              distance, isCalled, useDoubleInput_);
  }

  mtu::TorchMergeTree<float> originCopy;
  mtu::copyTorchMergeTree<float>(layers_[0].getOrigin(), originCopy);
  postprocessingPipeline<float>(&(originCopy.mTree.tree));
  for(unsigned int i = 0; i < customRecs_.size(); ++i) {
    fixTreePrecisionScalars(customRecs_[i].mTree);
    postprocessingPipeline<float>(&(customRecs_[i].mTree.tree));
    if(not isPersistenceDiagram_) {
      convertBranchDecompositionMatching<float>(&(originCopy.mTree.tree),
                                                &(customRecs_[i].mTree.tree),
                                                customMatchings_[i]);
    }
  }
}

//  ---------------------------------------------------------------------------
//  --- Utils
//  ---------------------------------------------------------------------------
unsigned int ttk::MergeTreeAutoencoder::getLatentLayerIndex() {
  unsigned int idx = noLayers_ / 2 - 1;
  if(idx > noLayers_) // unsigned negativeness
    idx = 0;
  return idx;
}

void ttk::MergeTreeAutoencoder::copyCustomParams(bool get) {
  auto &srcLatentCentroids = (get ? latentCentroids_ : bestLatentCentroids_);
  auto &dstLatentCentroids = (!get ? latentCentroids_ : bestLatentCentroids_);
  dstLatentCentroids.resize(srcLatentCentroids.size());
  for(unsigned int i = 0; i < dstLatentCentroids.size(); ++i)
    mtu::copyTensor(srcLatentCentroids[i], dstLatentCentroids[i]);
}

//  ---------------------------------------------------------------------------
//  --- Main Functions
//  ---------------------------------------------------------------------------
void ttk::MergeTreeAutoencoder::executeEndFunction(
  std::vector<ftm::MergeTree<float>> &trees,
  std::vector<ftm::MergeTree<float>> &ttkNotUsed(trees2)) {
  // Tracking
  computeTrackingInformation(getLatentLayerIndex() + 1);
  // Correlation
  computeCorrelationMatrix(trees, getLatentLayerIndex());
  // Custom recs
  createCustomRecs();
}
#endif
