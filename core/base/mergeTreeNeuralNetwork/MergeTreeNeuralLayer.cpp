#include <MergeTreeNeuralLayer.h>
#include <cmath>

#ifdef TTK_ENABLE_TORCH
using namespace torch::indexing;
#endif

ttk::MergeTreeNeuralLayer::MergeTreeNeuralLayer() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("MergeTreeNeuralLayer");
}

#ifdef TTK_ENABLE_TORCH
//  -----------------------------------------------------------------------
//  --- Getter/Setter
//  -----------------------------------------------------------------------
const ttk::mtu::TorchMergeTree<float> &
  ttk::MergeTreeNeuralLayer::getOrigin() const {
  return origin_;
}

const ttk::mtu::TorchMergeTree<float> &
  ttk::MergeTreeNeuralLayer::getOriginPrime() const {
  return originPrime_;
}

const ttk::mtu::TorchMergeTree<float> &
  ttk::MergeTreeNeuralLayer::getOrigin2() const {
  return origin2_;
}

const ttk::mtu::TorchMergeTree<float> &
  ttk::MergeTreeNeuralLayer::getOrigin2Prime() const {
  return origin2Prime_;
}

const torch::Tensor &ttk::MergeTreeNeuralLayer::getVSTensor() const {
  return vSTensor_;
}

const torch::Tensor &ttk::MergeTreeNeuralLayer::getVSPrimeTensor() const {
  return vSPrimeTensor_;
}

const torch::Tensor &ttk::MergeTreeNeuralLayer::getVS2Tensor() const {
  return vS2Tensor_;
}

const torch::Tensor &ttk::MergeTreeNeuralLayer::getVS2PrimeTensor() const {
  return vS2PrimeTensor_;
}

void ttk::MergeTreeNeuralLayer::setOrigin(
  const mtu::TorchMergeTree<float> &tmt) {
  mtu::copyTorchMergeTree(tmt, origin_);
}

void ttk::MergeTreeNeuralLayer::setOriginPrime(
  const mtu::TorchMergeTree<float> &tmt) {
  mtu::copyTorchMergeTree(tmt, originPrime_);
}

void ttk::MergeTreeNeuralLayer::setOrigin2(
  const mtu::TorchMergeTree<float> &tmt) {
  mtu::copyTorchMergeTree(tmt, origin2_);
}

void ttk::MergeTreeNeuralLayer::setOrigin2Prime(
  const mtu::TorchMergeTree<float> &tmt) {
  mtu::copyTorchMergeTree(tmt, origin2Prime_);
}

void ttk::MergeTreeNeuralLayer::setVSTensor(const torch::Tensor &vS) {
  mtu::copyTensor(vS, vSTensor_);
}

void ttk::MergeTreeNeuralLayer::setVSPrimeTensor(const torch::Tensor &vS) {
  mtu::copyTensor(vS, vSPrimeTensor_);
}

void ttk::MergeTreeNeuralLayer::setVS2Tensor(const torch::Tensor &vS) {
  mtu::copyTensor(vS, vS2Tensor_);
}

void ttk::MergeTreeNeuralLayer::setVS2PrimeTensor(const torch::Tensor &vS) {
  mtu::copyTensor(vS, vS2PrimeTensor_);
}

//  ---------------------------------------------------------------------------
//  --- Init
//  ---------------------------------------------------------------------------
void ttk::MergeTreeNeuralLayer::initOutputBasisTreeStructure(
  mtu::TorchMergeTree<float> &originPrime,
  bool isJT,
  mtu::TorchMergeTree<float> &baseOrigin) {
  // ----- Create scalars vector
  torch::Tensor originTensor = originPrime.tensor;
  if(!originTensor.device().is_cpu())
    originTensor = originTensor.cpu();
  std::vector<float> scalarsVector(
    originTensor.data_ptr<float>(),
    originTensor.data_ptr<float>() + originTensor.numel());
  unsigned int noNodes = scalarsVector.size() / 2;
  std::vector<std::vector<ftm::idNode>> childrenFinal(noNodes);

  // ----- Init tree structure and modify scalars if necessary
  if(isPersistenceDiagram_) {
    for(unsigned int i = 2; i < scalarsVector.size(); i += 2)
      childrenFinal[0].emplace_back(i / 2);
  } else {
    // --- Fix or swap min-max pair
    float maxPers = std::numeric_limits<float>::lowest();
    unsigned int indMax = 0;
    for(unsigned int i = 0; i < scalarsVector.size(); i += 2) {
      if(maxPers < (scalarsVector[i + 1] - scalarsVector[i])) {
        maxPers = (scalarsVector[i + 1] - scalarsVector[i]);
        indMax = i;
      }
    }
    if(indMax != 0) {
      float temp = scalarsVector[0];
      scalarsVector[0] = scalarsVector[indMax];
      scalarsVector[indMax] = temp;
      temp = scalarsVector[1];
      scalarsVector[1] = scalarsVector[indMax + 1];
      scalarsVector[indMax + 1] = temp;
    }
    ftm::idNode refNode = 0;
    for(unsigned int i = 2; i < scalarsVector.size(); i += 2) {
      ftm::idNode node = i / 2;
      adjustNestingScalars(scalarsVector, node, refNode);
    }

    if(not initOriginPrimeStructByCopy_
       or (int) noNodes > baseOrigin.mTree.tree.getRealNumberOfNodes()) {
      // --- Get possible children and parent relations
      std::vector<std::vector<ftm::idNode>> parents(noNodes), children(noNodes);
      for(unsigned int i = 0; i < scalarsVector.size(); i += 2) {
        for(unsigned int j = i; j < scalarsVector.size(); j += 2) {
          if(i == j)
            continue;
          unsigned int iN = i / 2, jN = j / 2;
          if(scalarsVector[i] <= scalarsVector[j]
             and scalarsVector[i + 1] >= scalarsVector[j + 1]) {
            // - i is parent of j
            parents[jN].emplace_back(iN);
            children[iN].emplace_back(jN);
          } else if(scalarsVector[i] >= scalarsVector[j]
                    and scalarsVector[i + 1] <= scalarsVector[j + 1]) {
            // - j is parent of i
            parents[iN].emplace_back(jN);
            children[jN].emplace_back(iN);
          }
        }
      }
      createBalancedBDT(parents, children, scalarsVector, childrenFinal);
    } else {
      ftm::MergeTree<float> mTreeTemp
        = ftm::copyMergeTree<float>(baseOrigin.mTree);
      bool useBD = true;
      keepMostImportantPairs<float>(&(mTreeTemp.tree), noNodes, useBD);
      torch::Tensor reshaped = torch::tensor(scalarsVector).reshape({-1, 2});
      torch::Tensor order = torch::argsort(
        (reshaped.index({Slice(), 1}) - reshaped.index({Slice(), 0})), -1,
        true);
      std::vector<unsigned int> nodeCorr(mTreeTemp.tree.getNumberOfNodes(), 0);
      unsigned int nodeNum = 1;
      std::queue<ftm::idNode> queue;
      queue.emplace(mTreeTemp.tree.getRoot());
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();
        std::vector<ftm::idNode> children;
        mTreeTemp.tree.getChildren(node, children);
        for(auto &child : children) {
          queue.emplace(child);
          unsigned int tNode = nodeCorr[node];
          nodeCorr[child] = order[nodeNum].item<int>();
          ++nodeNum;
          unsigned int tChild = nodeCorr[child];
          childrenFinal[tNode].emplace_back(tChild);
          adjustNestingScalars(scalarsVector, tChild, tNode);
        }
      }
    }
  }

  // ----- Create new tree
  originPrime.mTree = ftm::createEmptyMergeTree<float>(scalarsVector.size());
  ftm::FTMTree_MT *tree = &(originPrime.mTree.tree);
  if(isJT) {
    for(unsigned int i = 0; i < scalarsVector.size(); i += 2) {
      float temp = scalarsVector[i];
      scalarsVector[i] = scalarsVector[i + 1];
      scalarsVector[i + 1] = temp;
    }
  }
  ftm::setTreeScalars<float>(originPrime.mTree, scalarsVector);

  // ----- Create tree structure
  originPrime.nodeCorr.clear();
  originPrime.nodeCorr.assign(
    scalarsVector.size(), std::numeric_limits<unsigned int>::max());
  for(unsigned int i = 0; i < scalarsVector.size(); i += 2) {
    tree->makeNode(i);
    tree->makeNode(i + 1);
    tree->getNode(i)->setOrigin(i + 1);
    tree->getNode(i + 1)->setOrigin(i);
    originPrime.nodeCorr[i] = (unsigned int)(i / 2);
  }
  for(unsigned int i = 0; i < scalarsVector.size(); i += 2) {
    unsigned int node = i / 2;
    for(auto &child : childrenFinal[node])
      tree->makeSuperArc(child * 2, i);
  }
  mtu::getParentsVector(originPrime.mTree, originPrime.parentsOri);

  if(isTreeHasBigValues(originPrime.mTree, bigValuesThreshold_)) {
    std::stringstream ss;
    ss << originPrime.mTree.tree.printPairsFromTree<float>(true).str()
       << std::endl;
    ss << "isTreeHasBigValues(originPrime.mTree)" << std::endl;
    ss << "pause" << std::endl;
    printMsg(ss.str());
    std::cin.get();
  }
}

void ttk::MergeTreeNeuralLayer::initOutputBasis(
  const unsigned int dim,
  const unsigned int dim2,
  const torch::Tensor &baseTensor) {
  unsigned int originSize = origin_.tensor.sizes()[0];
  unsigned int origin2Size = 0;
  if(useDoubleInput_)
    origin2Size = origin2_.tensor.sizes()[0];

  // --- Compute output basis origin
  printMsg("Compute output basis origin", debug::Priority::DETAIL);
  auto initOutputBasisOrigin = [this, &baseTensor](
                                 torch::Tensor &w,
                                 mtu::TorchMergeTree<float> &tmt,
                                 mtu::TorchMergeTree<float> &baseTmt) {
    // - Create scalars
    torch::nn::init::xavier_normal_(w);
    torch::Tensor baseTmtTensor = baseTmt.tensor;
    if(normalizedWasserstein_)
      // Work on unnormalized tensor
      mtu::mergeTreeToTorchTensor(baseTmt.mTree, baseTmtTensor, false);
    torch::Tensor b
      = torch::full({w.sizes()[0], 1}, 0.01,
                    torch::TensorOptions().device(baseTmtTensor.device()));
    tmt.tensor = (torch::matmul(w, baseTmtTensor) + b);
    // - Shift to keep mean birth and max pers
    mtu::meanBirthMaxPersShift(tmt.tensor, baseTmtTensor);
    // - Shift to avoid diagonal points
    mtu::belowDiagonalPointsShift(tmt.tensor, baseTmtTensor);
    //
    if(initOriginPrimeValuesByCopy_) {
      auto baseTensorDiag = baseTensor.reshape({-1, 2});
      auto basePersDiag = (baseTensorDiag.index({Slice(), 1})
                           - baseTensorDiag.index({Slice(), 0}));
      auto tmtTensorDiag = tmt.tensor.reshape({-1, 2});
      auto persDiag = (tmtTensorDiag.index({Slice(1, None), 1})
                       - tmtTensorDiag.index({Slice(1, None), 0}));
      int noK = std::min(baseTensorDiag.sizes()[0], tmtTensorDiag.sizes()[0]);
      auto topVal = baseTensorDiag.index({std::get<1>(basePersDiag.topk(noK))});
      auto indexes = std::get<1>(persDiag.topk(noK - 1)) + 1;
      auto zeros
        = torch::zeros(1, torch::TensorOptions().device(indexes.device()));
      indexes = torch::cat({zeros, indexes}).to(torch::kLong);
      if(initOriginPrimeValuesByCopyRandomness_ != 0) {
        topVal = (1 - initOriginPrimeValuesByCopyRandomness_) * topVal
                 + initOriginPrimeValuesByCopyRandomness_
                     * tmtTensorDiag.index({indexes});
      }
      tmtTensorDiag.index_put_({indexes}, topVal);
    }
    // - Create tree structure
    initOutputBasisTreeStructure(
      tmt, baseTmt.mTree.tree.isJoinTree<float>(), baseTmt);
    if(normalizedWasserstein_)
      // Normalize tensor
      mtu::mergeTreeToTorchTensor(tmt.mTree, tmt.tensor, true);
    // - Projection
    interpolationProjection(tmt);
  };
  torch::Tensor w = torch::zeros(
    {dim, originSize}, torch::TensorOptions().device(origin_.tensor.device()));
  initOutputBasisOrigin(w, originPrime_, origin_);
  torch::Tensor w2;
  if(useDoubleInput_) {
    w2 = torch::zeros({dim2, origin2Size},
                      torch::TensorOptions().device(origin2_.tensor.device()));
    initOutputBasisOrigin(w2, origin2Prime_, origin2_);
  }

  // --- Compute output basis vectors
  printMsg("Compute output basis vectors", debug::Priority::DETAIL);
  initOutputBasisVectors(w, w2);
}

void ttk::MergeTreeNeuralLayer::initOutputBasisVectors(torch::Tensor &w,
                                                       torch::Tensor &w2) {
  vSPrimeTensor_ = torch::matmul(w, vSTensor_);
  if(useDoubleInput_)
    vS2PrimeTensor_ = torch::matmul(w2, vS2Tensor_);
  if(normalizedWasserstein_) {
    mtu::normalizeVectors(originPrime_.tensor, vSPrimeTensor_);
    if(useDoubleInput_)
      mtu::normalizeVectors(origin2Prime_.tensor, vS2PrimeTensor_);
  }
}

void ttk::MergeTreeNeuralLayer::initOutputBasisVectors(unsigned int dim,
                                                       unsigned int dim2) {
  unsigned int originSize = origin_.tensor.sizes()[0];
  unsigned int origin2Size = 0;
  if(useDoubleInput_)
    origin2Size = origin2_.tensor.sizes()[0];
  torch::Tensor w = torch::zeros({dim, originSize});
  torch::nn::init::xavier_normal_(w);
  torch::Tensor w2 = torch::zeros({dim2, origin2Size});
  torch::nn::init::xavier_normal_(w2);
  initOutputBasisVectors(w, w2);
}

void ttk::MergeTreeNeuralLayer::initInputBasisOrigin(
  std::vector<ftm::MergeTree<float>> &treesToUse,
  std::vector<ftm::MergeTree<float>> &trees2ToUse,
  double barycenterSizeLimitPercent,
  unsigned int barycenterMaxNoPairs,
  unsigned int barycenterMaxNoPairs2,
  std::vector<double> &inputToBaryDistances,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &baryMatchings,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &baryMatchings2) {
  int barycenterInitIndex = -1;
  if(initBarycenterRandom_) {
    std::random_device rd;
    std::default_random_engine rng(deterministic_ ? 0 : rd());
    barycenterInitIndex
      = std::uniform_int_distribution<>(0, treesToUse.size() - 1)(rng);
  }
  int maxNoPairs = (initBarycenterRandom_ ? barycenterMaxNoPairs : 0);
  computeOneBarycenter<float>(treesToUse, origin_.mTree, baryMatchings,
                              inputToBaryDistances, barycenterSizeLimitPercent,
                              maxNoPairs, barycenterInitIndex,
                              initBarycenterOneIter_, useDoubleInput_, true);
  if(not initBarycenterRandom_ and barycenterMaxNoPairs > 0)
    keepMostImportantPairs<float>(
      &(origin_.mTree.tree), barycenterMaxNoPairs, true);
  if(useDoubleInput_) {
    std::vector<double> baryDistances2;
    int maxNoPairs2 = (initBarycenterRandom_ ? barycenterMaxNoPairs2 : 0);
    computeOneBarycenter<float>(trees2ToUse, origin2_.mTree, baryMatchings2,
                                baryDistances2, barycenterSizeLimitPercent,
                                maxNoPairs2, barycenterInitIndex,
                                initBarycenterOneIter_, useDoubleInput_, false);
    if(not initBarycenterRandom_ and barycenterMaxNoPairs2 > 0)
      keepMostImportantPairs<float>(
        &(origin2_.mTree.tree), barycenterMaxNoPairs2, true);
    for(unsigned int i = 0; i < inputToBaryDistances.size(); ++i)
      inputToBaryDistances[i]
        = mixDistances(inputToBaryDistances[i], baryDistances2[i]);
  }

  mtu::getParentsVector(origin_.mTree, origin_.parentsOri);
  mtu::mergeTreeToTorchTensor<float>(
    origin_.mTree, origin_.tensor, origin_.nodeCorr, normalizedWasserstein_);
  if(useGpu_)
    origin_.tensor = origin_.tensor.cuda();
  if(useDoubleInput_) {
    mtu::getParentsVector(origin2_.mTree, origin2_.parentsOri);
    mtu::mergeTreeToTorchTensor<float>(origin2_.mTree, origin2_.tensor,
                                       origin2_.nodeCorr,
                                       normalizedWasserstein_);
    if(useGpu_)
      origin2_.tensor = origin2_.tensor.cuda();
  }
}

void ttk::MergeTreeNeuralLayer::initInputBasisVectors(
  std::vector<mtu::TorchMergeTree<float>> &tmTrees,
  std::vector<mtu::TorchMergeTree<float>> &tmTrees2,
  std::vector<ftm::MergeTree<float>> &trees,
  std::vector<ftm::MergeTree<float>> &trees2,
  unsigned int noVectors,
  std::vector<torch::Tensor> &allAlphasInit,
  std::vector<double> &inputToBaryDistances,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &baryMatchings,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &baryMatchings2,
  mtu::TorchMergeTree<float> &origin,
  mtu::TorchMergeTree<float> &origin2,
  torch::Tensor &vSTensor,
  torch::Tensor &vS2Tensor,
  bool useInputBasis) {
  if(randomAxesInit_) {
    auto initRandomAxes = [&noVectors](mtu::TorchMergeTree<float> &originT,
                                       torch::Tensor &axes) {
      torch::Tensor w = torch::zeros({noVectors, originT.tensor.sizes()[0]});
      torch::nn::init::xavier_normal_(w);
      axes = torch::linalg_pinv(w);
    };
    initRandomAxes(origin, vSTensor);
    if(useGpu_)
      vSTensor = vSTensor.cuda();
    if(useDoubleInput_) {
      initRandomAxes(origin2, vS2Tensor);
      if(useGpu_)
        vS2Tensor = vS2Tensor.cuda();
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < trees.size(); ++i)
      allAlphasInit[i] = torch::randn({noVectors, 1});
    return;
  }

  // --- Initialized vectors projection function to avoid collinearity
  auto initializedVectorsProjection
    = [=](int ttkNotUsed(_axeNumber),
          ftm::MergeTree<float> &ttkNotUsed(_barycenter),
          std::vector<std::vector<double>> &_v,
          std::vector<std::vector<double>> &ttkNotUsed(_v2),
          std::vector<std::vector<std::vector<double>>> &_vS,
          std::vector<std::vector<std::vector<double>>> &ttkNotUsed(_v2s),
          ftm::MergeTree<float> &ttkNotUsed(_barycenter2),
          std::vector<std::vector<double>> &ttkNotUsed(_trees2V),
          std::vector<std::vector<double>> &ttkNotUsed(_trees2V2),
          std::vector<std::vector<std::vector<double>>> &ttkNotUsed(_trees2Vs),
          std::vector<std::vector<std::vector<double>>> &ttkNotUsed(_trees2V2s),
          bool ttkNotUsed(_useSecondInput),
          unsigned int ttkNotUsed(_noProjectionStep)) {
        std::vector<double> scaledV, scaledVSi;
        Geometry::flattenMultiDimensionalVector(_v, scaledV);
        Geometry::scaleVector(
          scaledV, 1.0 / Geometry::magnitude(scaledV), scaledV);
        for(unsigned int i = 0; i < _vS.size(); ++i) {
          Geometry::flattenMultiDimensionalVector(_vS[i], scaledVSi);
          Geometry::scaleVector(
            scaledVSi, 1.0 / Geometry::magnitude(scaledVSi), scaledVSi);
          auto prod = Geometry::dotProduct(scaledV, scaledVSi);
          double tol = 0.01;
          if(prod <= -1.0 + tol or prod >= 1.0 - tol) {
            // Reset vector to initialize it again
            for(unsigned int j = 0; j < _v.size(); ++j)
              for(unsigned int k = 0; k < _v[j].size(); ++k)
                _v[j][k] = 0;
            break;
          }
        }
        return 0;
      };

  // --- Init vectors
  std::vector<std::vector<double>> inputToAxesDistances;
  std::vector<std::vector<std::vector<double>>> vS, v2s, trees2Vs, trees2V2s;
  std::stringstream ss;
  for(unsigned int vecNum = 0; vecNum < noVectors; ++vecNum) {
    ss.str("");
    ss << "Compute vectors " << vecNum;
    printMsg(ss.str(), debug::Priority::VERBOSE);
    std::vector<std::vector<double>> v1, v2, trees2V1, trees2V2;
    int newVectorOffset = 0;
    bool projectInitializedVectors = true;
    int bestIndex = MergeTreeAxesAlgorithmBase::initVectors<float>(
      vecNum, origin.mTree, trees, origin2.mTree, trees2, v1, v2, trees2V1,
      trees2V2, newVectorOffset, inputToBaryDistances, baryMatchings,
      baryMatchings2, inputToAxesDistances, vS, v2s, trees2Vs, trees2V2s,
      projectInitializedVectors, initializedVectorsProjection);
    vS.emplace_back(v1);
    v2s.emplace_back(v2);
    trees2Vs.emplace_back(trees2V1);
    trees2V2s.emplace_back(trees2V2);

    ss.str("");
    ss << "bestIndex = " << bestIndex;
    printMsg(ss.str(), debug::Priority::VERBOSE);

    // Update inputToAxesDistances
    printMsg("Update inputToAxesDistances", debug::Priority::VERBOSE);
    inputToAxesDistances.resize(1, std::vector<double>(trees.size()));
    if(bestIndex == -1 and normalizedWasserstein_) {
      mtu::normalizeVectors(origin, vS[vS.size() - 1]);
      if(useDoubleInput_)
        mtu::normalizeVectors(origin2, trees2Vs[vS.size() - 1]);
    }
    mtu::axisVectorsToTorchTensor(origin.mTree, vS, vSTensor);
    if(useGpu_)
      vSTensor = vSTensor.cuda();
    if(useDoubleInput_) {
      mtu::axisVectorsToTorchTensor(origin2.mTree, trees2Vs, vS2Tensor);
      if(useGpu_)
        vS2Tensor = vS2Tensor.cuda();
    }
    mtu::TorchMergeTree<float> dummyTmt;
    std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
      dummyBaryMatching2;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < trees.size(); ++i) {
      auto &tmt2ToUse = (not useDoubleInput_ ? dummyTmt : tmTrees2[i]);
      if(not euclideanVectorsInit_) {
        unsigned int k = k_;
        auto newAlpha = torch::ones({1, 1});
        if(bestIndex == -1) {
          newAlpha = torch::zeros({1, 1});
        }
        allAlphasInit[i] = (allAlphasInit[i].defined()
                              ? torch::cat({allAlphasInit[i], newAlpha})
                              : newAlpha);
        torch::Tensor bestAlphas;
        bool isCalled = true;
        inputToAxesDistances[0][i]
          = assignmentOneData(tmTrees[i], tmt2ToUse, k, allAlphasInit[i],
                              bestAlphas, isCalled, useInputBasis);
        allAlphasInit[i] = bestAlphas.detach();
      } else {
        auto &baryMatching2ToUse
          = (not useDoubleInput_ ? dummyBaryMatching2 : baryMatchings2[i]);
        torch::Tensor alphas;
        computeAlphas(tmTrees[i], origin, vSTensor, origin, baryMatchings[i],
                      tmt2ToUse, origin2, vS2Tensor, origin2,
                      baryMatching2ToUse, alphas);
        mtu::TorchMergeTree<float> interpolated, interpolated2;
        getMultiInterpolation(origin, vSTensor, alphas, interpolated);
        if(useDoubleInput_)
          getMultiInterpolation(origin2, vS2Tensor, alphas, interpolated2);
        torch::Tensor tensorDist;
        bool doSqrt = true;
        getDifferentiableDistanceFromMatchings(
          interpolated, tmTrees[i], interpolated2, tmt2ToUse, baryMatchings[i],
          baryMatching2ToUse, tensorDist, doSqrt);
        inputToAxesDistances[0][i] = tensorDist.item<double>();
        allAlphasInit[i] = alphas.detach();
      }
    }
  }
}

void ttk::MergeTreeNeuralLayer::initInputBasisVectors(
  std::vector<mtu::TorchMergeTree<float>> &tmTrees,
  std::vector<mtu::TorchMergeTree<float>> &tmTrees2,
  std::vector<ftm::MergeTree<float>> &trees,
  std::vector<ftm::MergeTree<float>> &trees2,
  unsigned int noVectors,
  std::vector<torch::Tensor> &allAlphasInit,
  std::vector<double> &inputToBaryDistances,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &baryMatchings,
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
    &baryMatchings2,
  bool useInputBasis) {
  mtu::TorchMergeTree<float> &origin = (useInputBasis ? origin_ : originPrime_);
  mtu::TorchMergeTree<float> &origin2
    = (useInputBasis ? origin2_ : origin2Prime_);
  torch::Tensor &vSTensor = (useInputBasis ? vSTensor_ : vSPrimeTensor_);
  torch::Tensor &vS2Tensor = (useInputBasis ? vS2Tensor_ : vS2PrimeTensor_);

  initInputBasisVectors(tmTrees, tmTrees2, trees, trees2, noVectors,
                        allAlphasInit, inputToBaryDistances, baryMatchings,
                        baryMatchings2, origin, origin2, vSTensor, vS2Tensor,
                        useInputBasis);
}

void ttk::MergeTreeNeuralLayer::requires_grad(const bool requireGrad) {
  origin_.tensor.requires_grad_(requireGrad);
  originPrime_.tensor.requires_grad_(requireGrad);
  vSTensor_.requires_grad_(requireGrad);
  vSPrimeTensor_.requires_grad_(requireGrad);
  if(useDoubleInput_) {
    origin2_.tensor.requires_grad_(requireGrad);
    origin2Prime_.tensor.requires_grad_(requireGrad);
    vS2Tensor_.requires_grad_(requireGrad);
    vS2PrimeTensor_.requires_grad_(requireGrad);
  }
}

void ttk::MergeTreeNeuralLayer::cuda() {
  origin_.tensor = origin_.tensor.cuda();
  originPrime_.tensor = originPrime_.tensor.cuda();
  vSTensor_ = vSTensor_.cuda();
  vSPrimeTensor_ = vSPrimeTensor_.cuda();
  if(useDoubleInput_) {
    origin2_.tensor = origin2_.tensor.cuda();
    origin2Prime_.tensor = origin2Prime_.tensor.cuda();
    vS2Tensor_ = vS2Tensor_.cuda();
    vS2PrimeTensor_ = vS2PrimeTensor_.cuda();
  }
}

//  ---------------------------------------------------------------------------
//  --- Interpolation
//  ---------------------------------------------------------------------------
void ttk::MergeTreeNeuralLayer::interpolationDiagonalProjection(
  mtu::TorchMergeTree<float> &interpolation) {
  torch::Tensor diagTensor = interpolation.tensor.reshape({-1, 2});
  if(interpolation.tensor.requires_grad())
    diagTensor = diagTensor.detach();

  torch::Tensor birthTensor = diagTensor.index({Slice(), 0});
  torch::Tensor deathTensor = diagTensor.index({Slice(), 1});

  torch::Tensor indexer = (birthTensor > deathTensor);

  torch::Tensor allProj = (birthTensor + deathTensor) / 2.0;
  allProj = allProj.index({indexer});
  allProj = allProj.reshape({-1, 1});

  diagTensor.index_put_({indexer}, allProj);
}

void ttk::MergeTreeNeuralLayer::interpolationNestingProjection(
  mtu::TorchMergeTree<float> &interpolation) {
  torch::Tensor diagTensor = interpolation.tensor.reshape({-1, 2});
  if(interpolation.tensor.requires_grad())
    diagTensor = diagTensor.detach();

  torch::Tensor birthTensor = diagTensor.index({Slice(1, None), 0});
  torch::Tensor deathTensor = diagTensor.index({Slice(1, None), 1});

  torch::Tensor birthIndexer = (birthTensor < 0);
  torch::Tensor deathIndexer = (deathTensor < 0);
  birthTensor.index_put_(
    {birthIndexer}, torch::zeros_like(birthTensor.index({birthIndexer})));
  deathTensor.index_put_(
    {deathIndexer}, torch::zeros_like(deathTensor.index({deathIndexer})));

  birthIndexer = (birthTensor > 1);
  deathIndexer = (deathTensor > 1);
  birthTensor.index_put_(
    {birthIndexer}, torch::ones_like(birthTensor.index({birthIndexer})));
  deathTensor.index_put_(
    {deathIndexer}, torch::ones_like(deathTensor.index({deathIndexer})));
}

void ttk::MergeTreeNeuralLayer::interpolationProjection(
  mtu::TorchMergeTree<float> &interpolation) {
  interpolationDiagonalProjection(interpolation);
  if(normalizedWasserstein_)
    interpolationNestingProjection(interpolation);

  ftm::MergeTree<float> interpolationNew;
  bool noRoot = mtu::torchTensorToMergeTree<float>(
    interpolation, normalizedWasserstein_, interpolationNew);
  if(noRoot)
    printWrn("[interpolationProjection] no root found");
  interpolation.mTree = copyMergeTree(interpolationNew);

  persistenceThresholding<float>(&(interpolation.mTree.tree), 0.001);

  if(isPersistenceDiagram_ and isThereMissingPairs(interpolation))
    printWrn("[getMultiInterpolation] missing pairs");
}

void ttk::MergeTreeNeuralLayer::getMultiInterpolation(
  const mtu::TorchMergeTree<float> &origin,
  const torch::Tensor &vS,
  torch::Tensor &alphas,
  mtu::TorchMergeTree<float> &interpolation) {
  mtu::copyTorchMergeTree<float>(origin, interpolation);
  interpolation.tensor = origin.tensor + torch::matmul(vS, alphas);
  interpolationProjection(interpolation);
}

//  ---------------------------------------------------------------------------
//  --- Forward
//  ---------------------------------------------------------------------------
void ttk::MergeTreeNeuralLayer::getAlphasOptimizationTensors(
  mtu::TorchMergeTree<float> &tree,
  mtu::TorchMergeTree<float> &origin,
  torch::Tensor &vSTensor,
  mtu::TorchMergeTree<float> &interpolated,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
  torch::Tensor &reorderedTreeTensor,
  torch::Tensor &deltaOrigin,
  torch::Tensor &deltaA,
  torch::Tensor &originTensor_f,
  torch::Tensor &vSTensor_f) {
  // Create matching indexing
  std::vector<int> tensorMatching;
  mtu::getTensorMatching(interpolated, tree, matching, tensorMatching);

  torch::Tensor indexes = torch::tensor(tensorMatching);
  torch::Tensor projIndexer = (indexes == -1).reshape({-1, 1});

  dataReorderingGivenMatching(
    origin, tree, projIndexer, indexes, reorderedTreeTensor, deltaOrigin);

  // Create axes projection given matching
  deltaA = vSTensor.transpose(0, 1).reshape({vSTensor.sizes()[1], -1, 2});
  deltaA = (deltaA.index({Slice(), Slice(), 0})
            + deltaA.index({Slice(), Slice(), 1}))
           / 2.0;
  deltaA = torch::stack({deltaA, deltaA}, 2);
  if(!deltaA.device().is_cpu())
    projIndexer = projIndexer.to(deltaA.device());
  deltaA = deltaA * projIndexer;
  deltaA = deltaA.reshape({vSTensor.sizes()[1], -1}).transpose(0, 1);

  //
  originTensor_f = origin.tensor;
  vSTensor_f = vSTensor;
}

void ttk::MergeTreeNeuralLayer::computeAlphas(
  mtu::TorchMergeTree<float> &tree,
  mtu::TorchMergeTree<float> &origin,
  torch::Tensor &vSTensor,
  mtu::TorchMergeTree<float> &interpolated,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
  mtu::TorchMergeTree<float> &tree2,
  mtu::TorchMergeTree<float> &origin2,
  torch::Tensor &vS2Tensor,
  mtu::TorchMergeTree<float> &interpolated2,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching2,
  torch::Tensor &alphasOut) {
  torch::Tensor reorderedTreeTensor, deltaOrigin, deltaA, originTensor_f,
    vSTensor_f;
  getAlphasOptimizationTensors(tree, origin, vSTensor, interpolated, matching,
                               reorderedTreeTensor, deltaOrigin, deltaA,
                               originTensor_f, vSTensor_f);

  if(useDoubleInput_) {
    torch::Tensor reorderedTree2Tensor, deltaOrigin2, deltaA2, origin2Tensor_f,
      vS2Tensor_f;
    getAlphasOptimizationTensors(tree2, origin2, vS2Tensor, interpolated2,
                                 matching2, reorderedTree2Tensor, deltaOrigin2,
                                 deltaA2, origin2Tensor_f, vS2Tensor_f);
    vSTensor_f = torch::cat({vSTensor_f, vS2Tensor_f});
    deltaA = torch::cat({deltaA, deltaA2});
    reorderedTreeTensor
      = torch::cat({reorderedTreeTensor, reorderedTree2Tensor});
    originTensor_f = torch::cat({originTensor_f, origin2Tensor_f});
    deltaOrigin = torch::cat({deltaOrigin, deltaOrigin2});
  }

  torch::Tensor r_axes = vSTensor_f - deltaA;
  torch::Tensor r_data = reorderedTreeTensor - originTensor_f + deltaOrigin;

  // Pseudo inverse
  auto driver = "gelsd";
  bool is_cpu = r_axes.device().is_cpu();
  auto device = r_axes.device();
  if(!is_cpu) {
    r_axes = r_axes.cpu();
    r_data = r_data.cpu();
  }
  alphasOut
    = std::get<0>(torch::linalg_lstsq(r_axes, r_data, c10::nullopt, driver));
  if(!is_cpu)
    alphasOut = alphasOut.to(device);

  alphasOut.reshape({-1, 1});
}

float ttk::MergeTreeNeuralLayer::assignmentOneData(
  mtu::TorchMergeTree<float> &tree,
  mtu::TorchMergeTree<float> &tree2,
  unsigned int k,
  torch::Tensor &alphasInit,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &bestMatching,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &bestMatching2,
  torch::Tensor &bestAlphas,
  bool isCalled,
  bool useInputBasis) {
  mtu::TorchMergeTree<float> &origin = (useInputBasis ? origin_ : originPrime_);
  mtu::TorchMergeTree<float> &origin2
    = (useInputBasis ? origin2_ : origin2Prime_);
  torch::Tensor &vSTensor = (useInputBasis ? vSTensor_ : vSPrimeTensor_);
  torch::Tensor &vS2Tensor = (useInputBasis ? vS2Tensor_ : vS2PrimeTensor_);

  torch::Tensor alphas, oldAlphas;
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching, matching2;
  float bestDistance = std::numeric_limits<float>::max();
  mtu::TorchMergeTree<float> interpolated, interpolated2;
  unsigned int i = 0;
  auto reset = [&]() {
    alphasInit = torch::randn_like(alphas);
    i = 0;
  };
  unsigned int noUpdate = 0;
  unsigned int noReset = 0;
  while(i < k) {
    if(i == 0) {
      if(alphasInit.defined())
        alphas = alphasInit;
      else
        alphas = torch::zeros({vSTensor.sizes()[1], 1});
    } else {
      computeAlphas(tree, origin, vSTensor, interpolated, matching, tree2,
                    origin2, vS2Tensor, interpolated2, matching2, alphas);
      if(oldAlphas.defined() and alphas.defined() and alphas.equal(oldAlphas)
         and i != 1) {
        break;
      }
    }
    mtu::copyTensor(alphas, oldAlphas);
    getMultiInterpolation(origin, vSTensor, alphas, interpolated);
    if(useDoubleInput_)
      getMultiInterpolation(origin2, vS2Tensor, alphas, interpolated2);
    if(interpolated.mTree.tree.getRealNumberOfNodes() == 0
       or (useDoubleInput_
           and interpolated2.mTree.tree.getRealNumberOfNodes() == 0)) {
      ++noReset;
      if(noReset >= 100)
        printWrn("[assignmentOneData] noReset >= 100");
      reset();
      continue;
    }
    float distance;
    computeOneDistance<float>(interpolated.mTree, tree.mTree, matching,
                              distance, isCalled, useDoubleInput_);
    if(useDoubleInput_) {
      float distance2;
      computeOneDistance<float>(interpolated2.mTree, tree2.mTree, matching2,
                                distance2, isCalled, useDoubleInput_, false);
      distance = mixDistances<float>(distance, distance2);
    }
    if(distance < bestDistance and i != 0) {
      bestDistance = distance;
      bestMatching = matching;
      bestMatching2 = matching2;
      bestAlphas = alphas;
      noUpdate += 1;
    }
    i += 1;
  }
  if(noUpdate == 0)
    printErr("[assignmentOneData] noUpdate ==  0");
  return bestDistance;
}

float ttk::MergeTreeNeuralLayer::assignmentOneData(
  mtu::TorchMergeTree<float> &tree,
  mtu::TorchMergeTree<float> &tree2,
  unsigned int k,
  torch::Tensor &alphasInit,
  torch::Tensor &bestAlphas,
  bool isCalled,
  bool useInputBasis) {
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> bestMatching,
    bestMatching2;
  return assignmentOneData(tree, tree2, k, alphasInit, bestMatching,
                           bestMatching2, bestAlphas, isCalled, useInputBasis);
}

void ttk::MergeTreeNeuralLayer::outputBasisReconstruction(
  torch::Tensor &alphas,
  mtu::TorchMergeTree<float> &out,
  mtu::TorchMergeTree<float> &out2,
  bool activate,
  bool train) {
  if(not activate_)
    activate = false;
  torch::Tensor act = (activate ? activation(alphas) : alphas);
  if(dropout_ != 0.0 and train) {
    torch::nn::Dropout model(torch::nn::DropoutOptions().p(dropout_));
    act = model(act);
  }
  getMultiInterpolation(originPrime_, vSPrimeTensor_, act, out);
  if(useDoubleInput_)
    getMultiInterpolation(origin2Prime_, vS2PrimeTensor_, act, out2);
}

bool ttk::MergeTreeNeuralLayer::forward(mtu::TorchMergeTree<float> &tree,
                                        mtu::TorchMergeTree<float> &tree2,
                                        unsigned int k,
                                        torch::Tensor &alphasInit,
                                        mtu::TorchMergeTree<float> &out,
                                        mtu::TorchMergeTree<float> &out2,
                                        torch::Tensor &bestAlphas,
                                        float &bestDistance,
                                        bool train) {
  bool goodOutput = false;
  int noReset = 0;
  while(not goodOutput) {
    bool isCalled = true;
    bestDistance
      = assignmentOneData(tree, tree2, k, alphasInit, bestAlphas, isCalled);
    outputBasisReconstruction(bestAlphas, out, out2, true, train);
    goodOutput = (out.mTree.tree.getRealNumberOfNodes() != 0
                  and (not useDoubleInput_
                       or out2.mTree.tree.getRealNumberOfNodes() != 0));
    if(not goodOutput) {
      ++noReset;
      if(noReset >= 100) {
        printWrn("[forwardOneLayer] noReset >= 100");
        return true;
      }
      alphasInit = torch::randn_like(alphasInit);
    }
  }
  return false;
}

bool ttk::MergeTreeNeuralLayer::forward(mtu::TorchMergeTree<float> &tree,
                                        mtu::TorchMergeTree<float> &tree2,
                                        unsigned int k,
                                        torch::Tensor &alphasInit,
                                        mtu::TorchMergeTree<float> &out,
                                        mtu::TorchMergeTree<float> &out2,
                                        torch::Tensor &bestAlphas,
                                        bool train) {
  float bestDistance;
  return forward(
    tree, tree2, k, alphasInit, out, out2, bestAlphas, bestDistance, train);
}

//  ---------------------------------------------------------------------------
//  --- Projection
//  ---------------------------------------------------------------------------
void ttk::MergeTreeNeuralLayer::projectionStep() {
  auto projectTree = [this](mtu::TorchMergeTree<float> &tmt) {
    interpolationProjection(tmt);
    tmt.tensor = tmt.tensor.detach();
    tmt.tensor.requires_grad_(true);
  };
  projectTree(origin_);
  projectTree(originPrime_);
  if(useDoubleInput_) {
    projectTree(origin2_);
    projectTree(origin2Prime_);
  }
}

//  ---------------------------------------------------------------------------
//  --- Utils
//  ---------------------------------------------------------------------------
void ttk::MergeTreeNeuralLayer::copyParams(
  mtu::TorchMergeTree<float> &origin,
  mtu::TorchMergeTree<float> &originPrime,
  torch::Tensor &vS,
  torch::Tensor &vSPrime,
  mtu::TorchMergeTree<float> &origin2,
  mtu::TorchMergeTree<float> &origin2Prime,
  torch::Tensor &vS2,
  torch::Tensor &vS2Prime,
  bool get) {

  // Source
  mtu::TorchMergeTree<float> &srcOrigin = (get ? origin_ : origin);
  mtu::TorchMergeTree<float> &srcOriginPrime
    = (get ? originPrime_ : originPrime);
  torch::Tensor &srcVS = (get ? vSTensor_ : vS);
  torch::Tensor &srcVSPrime = (get ? vSPrimeTensor_ : vSPrime);
  mtu::TorchMergeTree<float> &srcOrigin2 = (get ? origin2_ : origin2);
  mtu::TorchMergeTree<float> &srcOrigin2Prime
    = (get ? origin2Prime_ : origin2Prime);
  torch::Tensor &srcVS2 = (get ? vS2Tensor_ : vS2);
  torch::Tensor &srcVS2Prime = (get ? vS2PrimeTensor_ : vS2Prime);

  // Destination
  mtu::TorchMergeTree<float> &dstOrigin = (!get ? origin_ : origin);
  mtu::TorchMergeTree<float> &dstOriginPrime
    = (!get ? originPrime_ : originPrime);
  torch::Tensor &dstVS = (!get ? vSTensor_ : vS);
  torch::Tensor &dstVSPrime = (!get ? vSPrimeTensor_ : vSPrime);
  mtu::TorchMergeTree<float> &dstOrigin2 = (!get ? origin2_ : origin2);
  mtu::TorchMergeTree<float> &dstOrigin2Prime
    = (!get ? origin2Prime_ : origin2Prime);
  torch::Tensor &dstVS2 = (!get ? vS2Tensor_ : vS2);
  torch::Tensor &dstVS2Prime = (!get ? vS2PrimeTensor_ : vS2Prime);

  // Copy
  mtu::copyTorchMergeTree(srcOrigin, dstOrigin);
  mtu::copyTorchMergeTree(srcOriginPrime, dstOriginPrime);
  mtu::copyTensor(srcVS, dstVS);
  mtu::copyTensor(srcVSPrime, dstVSPrime);
  if(useDoubleInput_) {
    mtu::copyTorchMergeTree(srcOrigin2, dstOrigin2);
    mtu::copyTorchMergeTree(srcOrigin2Prime, dstOrigin2Prime);
    mtu::copyTensor(srcVS2, dstVS2);
    mtu::copyTensor(srcVS2Prime, dstVS2Prime);
  }
}

void ttk::MergeTreeNeuralLayer::adjustNestingScalars(
  std::vector<float> &scalarsVector, ftm::idNode node, ftm::idNode refNode) {
  float birth = scalarsVector[refNode * 2];
  float death = scalarsVector[refNode * 2 + 1];
  auto getSign = [](float v) { return (v > 0 ? 1 : -1); };
  auto getPrecValue = [&getSign](float v, bool opp = false) {
    return v * (1 + (opp ? -1 : 1) * getSign(v) * 1e-6);
  };
  // Shift scalars
  if(scalarsVector[node * 2 + 1] > getPrecValue(death, true)) {
    float diff = scalarsVector[node * 2 + 1] - getPrecValue(death, true);
    scalarsVector[node * 2] -= diff;
    scalarsVector[node * 2 + 1] -= diff;
  } else if(scalarsVector[node * 2] < getPrecValue(birth)) {
    float diff = getPrecValue(birth) - scalarsVector[node * 2];
    scalarsVector[node * 2] += getPrecValue(diff);
    scalarsVector[node * 2 + 1] += getPrecValue(diff);
  }
  // Cut scalars
  if(scalarsVector[node * 2] < getPrecValue(birth))
    scalarsVector[node * 2] = getPrecValue(birth);
  if(scalarsVector[node * 2 + 1] > getPrecValue(death, true))
    scalarsVector[node * 2 + 1] = getPrecValue(death, true);
}

void ttk::MergeTreeNeuralLayer::createBalancedBDT(
  std::vector<std::vector<ftm::idNode>> &parents,
  std::vector<std::vector<ftm::idNode>> &children,
  std::vector<float> &scalarsVector,
  std::vector<std::vector<ftm::idNode>> &childrenFinal) {
  // ----- Some variables
  unsigned int noNodes = scalarsVector.size() / 2;
  childrenFinal.resize(noNodes);
  int mtLevel = ceil(log(noNodes * 2) / log(2)) + 1;
  int bdtLevel = mtLevel - 1;
  int noDim = bdtLevel;

  // ----- Get node levels
  std::vector<int> nodeLevels(noNodes, -1);
  std::queue<ftm::idNode> queueLevels;
  std::vector<int> noChildDone(noNodes, 0);
  for(unsigned int i = 0; i < children.size(); ++i) {
    if(children[i].size() == 0) {
      queueLevels.emplace(i);
      nodeLevels[i] = 1;
    }
  }
  while(!queueLevels.empty()) {
    ftm::idNode node = queueLevels.front();
    queueLevels.pop();
    for(auto &parent : parents[node]) {
      ++noChildDone[parent];
      nodeLevels[parent] = std::max(nodeLevels[parent], nodeLevels[node] + 1);
      if(noChildDone[parent] >= (int)children[parent].size())
        queueLevels.emplace(parent);
    }
  }

  // ----- Sort heuristic lambda
  auto sortChildren = [this, &parents, &scalarsVector, &noNodes](
                        ftm::idNode nodeOrigin, std::vector<bool> &nodeDone,
                        std::vector<std::vector<ftm::idNode>> &childrenT) {
    double refPers = scalarsVector[1] - scalarsVector[0];
    auto getRemaining = [&nodeDone](std::vector<ftm::idNode> &vec) {
      unsigned int remaining = 0;
      for(auto &e : vec)
        remaining += (not nodeDone[e]);
      return remaining;
    };
    std::vector<unsigned int> parentsRemaining(noNodes, 0),
      childrenRemaining(noNodes, 0);
    for(auto &child : childrenT[nodeOrigin]) {
      parentsRemaining[child] = getRemaining(parents[child]);
      childrenRemaining[child] = getRemaining(childrenT[child]);
    }
    TTK_PSORT(
      threadNumber_, childrenT[nodeOrigin].begin(), childrenT[nodeOrigin].end(),
      [&](ftm::idNode nodeI, ftm::idNode nodeJ) {
        double persI = scalarsVector[nodeI * 2 + 1] - scalarsVector[nodeI * 2];
        double persJ = scalarsVector[nodeJ * 2 + 1] - scalarsVector[nodeJ * 2];
        return parentsRemaining[nodeI] + childrenRemaining[nodeI]
                 - persI / refPers * noNodes
               < parentsRemaining[nodeJ] + childrenRemaining[nodeJ]
                   - persJ / refPers * noNodes;
      });
  };

  // ----- Greedy approach to find balanced BDT structures
  const auto findStructGivenDim =
    [&children, &noNodes, &nodeLevels](
      ftm::idNode _nodeOrigin, int _dimToFound, bool _searchMaxDim,
      std::vector<bool> &_nodeDone, std::vector<bool> &_dimFound,
      std::vector<std::vector<ftm::idNode>> &_childrenFinalOut) {
      // --- Recursive lambda
      auto findStructGivenDimImpl =
        [&children, &noNodes, &nodeLevels](
          ftm::idNode nodeOrigin, int dimToFound, bool searchMaxDim,
          std::vector<bool> &nodeDone, std::vector<bool> &dimFound,
          std::vector<std::vector<ftm::idNode>> &childrenFinalOut,
          auto &findStructGivenDimRef) mutable {
          childrenFinalOut.resize(noNodes);
          // - Find structures
          int dim = (searchMaxDim ? dimToFound - 1 : 0);
          unsigned int i = 0;
          //
          auto searchMaxDimReset = [&i, &dim, &nodeDone]() {
            --dim;
            i = 0;
            unsigned int noDone = 0;
            for(auto done : nodeDone)
              if(done)
                ++noDone;
            return noDone == nodeDone.size() - 1; // -1 for root
          };
          while(i < children[nodeOrigin].size()) {
            auto child = children[nodeOrigin][i];
            // Skip if child was already processed
            if(nodeDone[child]) {
              // If we have processed all children while searching for max
              // dim then restart at the beginning to find a lower dim
              if(searchMaxDim and i == children[nodeOrigin].size() - 1) {
                if(searchMaxDimReset())
                  break;
              } else
                ++i;
              continue;
            }
            if(dim == 0) {
              // Base case
              childrenFinalOut[nodeOrigin].emplace_back(child);
              nodeDone[child] = true;
              dimFound[0] = true;
              if(dimToFound <= 1 or searchMaxDim)
                return true;
              ++dim;
            } else {
              // General case
              std::vector<std::vector<ftm::idNode>> childrenFinalDim;
              std::vector<bool> nodeDoneDim;
              std::vector<bool> dimFoundDim(dim);
              bool found = false;
              if(nodeLevels[child] > dim) {
                nodeDoneDim = nodeDone;
                found = findStructGivenDimRef(child, dim, false, nodeDoneDim,
                                              dimFoundDim, childrenFinalDim,
                                              findStructGivenDimRef);
              }
              if(found) {
                dimFound[dim] = true;
                childrenFinalOut[nodeOrigin].emplace_back(child);
                for(unsigned int j = 0; j < childrenFinalDim.size(); ++j)
                  for(auto &e : childrenFinalDim[j])
                    childrenFinalOut[j].emplace_back(e);
                nodeDone[child] = true;
                for(unsigned int j = 0; j < nodeDoneDim.size(); ++j)
                  nodeDone[j] = nodeDone[j] || nodeDoneDim[j];
                // Return if it is the last dim to found
                if(dim == dimToFound - 1 and not searchMaxDim)
                  return true;
                // Reset index if we search for the maximum dim
                if(searchMaxDim) {
                  if(searchMaxDimReset())
                    break;
                } else {
                  ++dim;
                }
                continue;
              } else if(searchMaxDim and i == children[nodeOrigin].size() - 1) {
                // If we have processed all children while searching for max dim
                // then restart at the beginning to find a lower dim
                if(searchMaxDimReset())
                  break;
                continue;
              }
            }
            ++i;
          }
          return false;
        };
      return findStructGivenDimImpl(_nodeOrigin, _dimToFound, _searchMaxDim,
                                    _nodeDone, _dimFound, _childrenFinalOut,
                                    findStructGivenDimImpl);
    };
  std::vector<bool> dimFound(noDim - 1, false);
  std::vector<bool> nodeDone(noNodes, false);
  for(unsigned int i = 0; i < children.size(); ++i)
    sortChildren(i, nodeDone, children);
  Timer t_find;
  ftm::idNode startNode = 0;
  findStructGivenDim(startNode, noDim, true, nodeDone, dimFound, childrenFinal);

  // ----- Greedy approach to create non found structures
  const auto createStructGivenDim =
    [this, &children, &noNodes, &findStructGivenDim, &nodeLevels](
      int _nodeOrigin, int _dimToCreate, std::vector<bool> &_nodeDone,
      ftm::idNode &_structOrigin, std::vector<float> &_scalarsVectorOut,
      std::vector<std::vector<ftm::idNode>> &_childrenFinalOut) {
      // --- Recursive lambda
      auto createStructGivenDimImpl =
        [this, &children, &noNodes, &findStructGivenDim, &nodeLevels](
          int nodeOrigin, int dimToCreate, std::vector<bool> &nodeDoneImpl,
          ftm::idNode &structOrigin, std::vector<float> &scalarsVectorOut,
          std::vector<std::vector<ftm::idNode>> &childrenFinalOut,
          auto &createStructGivenDimRef) mutable {
          // Deduction of auto lambda type
          if(false)
            return;
          // - Find structures of lower dimension
          int dimToFound = dimToCreate - 1;
          std::vector<std::vector<std::vector<ftm::idNode>>> childrenFinalT(2);
          std::array<ftm::idNode, 2> structOrigins;
          for(unsigned int n = 0; n < 2; ++n) {
            bool found = false;
            for(unsigned int i = 0; i < children[nodeOrigin].size(); ++i) {
              auto child = children[nodeOrigin][i];
              if(nodeDoneImpl[child])
                continue;
              if(dimToFound != 0) {
                if(nodeLevels[child] > dimToFound) {
                  std::vector<bool> dimFoundT(dimToFound, false);
                  childrenFinalT[n].clear();
                  childrenFinalT[n].resize(noNodes);
                  std::vector<bool> nodeDoneImplFind = nodeDoneImpl;
                  found = findStructGivenDim(child, dimToFound, false,
                                             nodeDoneImplFind, dimFoundT,
                                             childrenFinalT[n]);
                }
              } else
                found = true;
              if(found) {
                structOrigins[n] = child;
                nodeDoneImpl[child] = true;
                for(unsigned int j = 0; j < childrenFinalT[n].size(); ++j) {
                  for(auto &e : childrenFinalT[n][j]) {
                    childrenFinalOut[j].emplace_back(e);
                    nodeDoneImpl[e] = true;
                  }
                }
                break;
              }
            } // end for children[nodeOrigin]
            if(not found) {
              if(dimToFound <= 0) {
                structOrigins[n] = std::numeric_limits<ftm::idNode>::max();
                continue;
              }
              childrenFinalT[n].clear();
              childrenFinalT[n].resize(noNodes);
              createStructGivenDimRef(
                nodeOrigin, dimToFound, nodeDoneImpl, structOrigins[n],
                scalarsVectorOut, childrenFinalT[n], createStructGivenDimRef);
              for(unsigned int j = 0; j < childrenFinalT[n].size(); ++j) {
                for(auto &e : childrenFinalT[n][j]) {
                  if(e == structOrigins[n])
                    continue;
                  childrenFinalOut[j].emplace_back(e);
                }
              }
            }
          } // end for n
          // - Combine both structures
          if(structOrigins[0] == std::numeric_limits<ftm::idNode>::max()
             and structOrigins[1] == std::numeric_limits<ftm::idNode>::max()) {
            structOrigin = std::numeric_limits<ftm::idNode>::max();
            return;
          }
          bool firstIsParent = true;
          if(structOrigins[0] == std::numeric_limits<ftm::idNode>::max())
            firstIsParent = false;
          else if(structOrigins[1] == std::numeric_limits<ftm::idNode>::max())
            firstIsParent = true;
          else if(scalarsVectorOut[structOrigins[1] * 2 + 1]
                    - scalarsVectorOut[structOrigins[1] * 2]
                  > scalarsVectorOut[structOrigins[0] * 2 + 1]
                      - scalarsVectorOut[structOrigins[0] * 2])
            firstIsParent = false;
          structOrigin = (firstIsParent ? structOrigins[0] : structOrigins[1]);
          ftm::idNode modOrigin
            = (firstIsParent ? structOrigins[1] : structOrigins[0]);
          childrenFinalOut[nodeOrigin].emplace_back(structOrigin);
          if(modOrigin != std::numeric_limits<ftm::idNode>::max()) {
            childrenFinalOut[structOrigin].emplace_back(modOrigin);
            std::queue<std::array<ftm::idNode, 2>> queue;
            queue.emplace(std::array<ftm::idNode, 2>{modOrigin, structOrigin});
            while(!queue.empty()) {
              auto &nodeAndParent = queue.front();
              ftm::idNode node = nodeAndParent[0];
              ftm::idNode parent = nodeAndParent[1];
              queue.pop();
              adjustNestingScalars(scalarsVectorOut, node, parent);
              // Push children
              for(auto &child : childrenFinalOut[node])
                queue.emplace(std::array<ftm::idNode, 2>{child, node});
            }
          }
          return;
        };
      return createStructGivenDimImpl(
        _nodeOrigin, _dimToCreate, _nodeDone, _structOrigin, _scalarsVectorOut,
        _childrenFinalOut, createStructGivenDimImpl);
    };
  for(unsigned int i = 0; i < children.size(); ++i)
    sortChildren(i, nodeDone, children);
  Timer t_create;
  for(unsigned int i = 0; i < dimFound.size(); ++i) {
    if(dimFound[i])
      continue;
    ftm::idNode structOrigin;
    createStructGivenDim(
      startNode, i, nodeDone, structOrigin, scalarsVector, childrenFinal);
  }
}

//  ---------------------------------------------------------------------------
//  --- Testing
//  ---------------------------------------------------------------------------
bool ttk::MergeTreeNeuralLayer::isTreeHasBigValues(ftm::MergeTree<float> &mTree,
                                                   float threshold) {
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
