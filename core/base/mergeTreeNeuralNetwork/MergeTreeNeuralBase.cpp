#include <MergeTreeNeuralBase.h>
#include <cmath>

#ifdef TTK_ENABLE_TORCH
using namespace torch::indexing;
#endif

ttk::MergeTreeNeuralBase::MergeTreeNeuralBase() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("MergeTreeNeuralBase");
}

#ifdef TTK_ENABLE_TORCH
//  -----------------------------------------------------------------------
//  --- Setter
//  -----------------------------------------------------------------------
void ttk::MergeTreeNeuralBase::setDropout(const double dropout) {
  dropout_ = dropout;
}

void ttk::MergeTreeNeuralBase::setEuclideanVectorsInit(
  const bool euclideanVectorsInit) {
  euclideanVectorsInit_ = euclideanVectorsInit;
}

void ttk::MergeTreeNeuralBase::setRandomAxesInit(const bool randomAxesInit) {
  randomAxesInit_ = randomAxesInit;
}

void ttk::MergeTreeNeuralBase::setInitBarycenterRandom(
  const bool initBarycenterRandom) {
  initBarycenterRandom_ = initBarycenterRandom;
}

void ttk::MergeTreeNeuralBase::setInitBarycenterOneIter(
  const bool initBarycenterOneIter) {
  initBarycenterOneIter_ = initBarycenterOneIter;
}

void ttk::MergeTreeNeuralBase::setInitOriginPrimeStructByCopy(
  const bool initOriginPrimeStructByCopy) {
  initOriginPrimeStructByCopy_ = initOriginPrimeStructByCopy;
}

void ttk::MergeTreeNeuralBase::setInitOriginPrimeValuesByCopy(
  const bool initOriginPrimeValuesByCopy) {
  initOriginPrimeValuesByCopy_ = initOriginPrimeValuesByCopy;
}

void ttk::MergeTreeNeuralBase::setInitOriginPrimeValuesByCopyRandomness(
  const double initOriginPrimeValuesByCopyRandomness) {
  initOriginPrimeValuesByCopyRandomness_
    = initOriginPrimeValuesByCopyRandomness;
}

void ttk::MergeTreeNeuralBase::setActivate(const bool activate) {
  activate_ = activate;
}

void ttk::MergeTreeNeuralBase::setActivationFunction(
  const unsigned int activationFunction) {
  activationFunction_ = activationFunction;
}

void ttk::MergeTreeNeuralBase::setUseGpu(const bool useGpu) {
  useGpu_ = useGpu;
}

void ttk::MergeTreeNeuralBase::setBigValuesThreshold(
  const float bigValuesThreshold) {
  bigValuesThreshold_ = bigValuesThreshold;
}

//  -----------------------------------------------------------------------
//  --- Utils
//  -----------------------------------------------------------------------
torch::Tensor ttk::MergeTreeNeuralBase::activation(torch::Tensor &in) {
  torch::Tensor act;
  switch(activationFunction_) {
    case 1:
      act = torch::nn::LeakyReLU()(in);
      break;
    case 0:
    default:
      act = torch::nn::ReLU()(in);
  }
  return act;
}

void ttk::MergeTreeNeuralBase::fixTreePrecisionScalars(
  ftm::MergeTree<float> &mTree) {
  double eps = 1e-6;
  auto shiftSubtree
    = [&mTree, &eps](ftm::idNode node, ftm::idNode birthNodeParent,
                     ftm::idNode deathNodeParent, std::vector<float> &scalars,
                     bool invalidBirth, bool invalidDeath) {
        std::queue<ftm::idNode> queue;
        queue.emplace(node);
        while(!queue.empty()) {
          ftm::idNode nodeT = queue.front();
          queue.pop();
          auto birthDeathNode = mTree.tree.getBirthDeathNode<float>(node);
          auto birthNode = std::get<0>(birthDeathNode);
          auto deathNode = std::get<1>(birthDeathNode);
          if(invalidBirth)
            scalars[birthNode] = scalars[birthNodeParent] + 2 * eps;
          if(invalidDeath)
            scalars[deathNode] = scalars[deathNodeParent] - 2 * eps;
          std::vector<ftm::idNode> children;
          mTree.tree.getChildren(nodeT, children);
          for(auto &child : children)
            queue.emplace(child);
        }
      };
  std::vector<float> scalars;
  getTreeScalars(mTree, scalars);
  std::queue<ftm::idNode> queue;
  auto root = mTree.tree.getRoot();
  queue.emplace(root);
  while(!queue.empty()) {
    ftm::idNode node = queue.front();
    queue.pop();
    auto birthDeathNode = mTree.tree.getBirthDeathNode<float>(node);
    auto birthNode = std::get<0>(birthDeathNode);
    auto deathNode = std::get<1>(birthDeathNode);
    auto birthDeathNodeParent
      = mTree.tree.getBirthDeathNode<float>(mTree.tree.getParentSafe(node));
    auto birthNodeParent = std::get<0>(birthDeathNodeParent);
    auto deathNodeParent = std::get<1>(birthDeathNodeParent);
    bool invalidBirth = (scalars[birthNode] <= scalars[birthNodeParent] + eps);
    bool invalidDeath = (scalars[deathNode] >= scalars[deathNodeParent] - eps);
    if(!mTree.tree.isRoot(node) and (invalidBirth or invalidDeath))
      shiftSubtree(node, birthNodeParent, deathNodeParent, scalars,
                   invalidBirth, invalidDeath);
    std::vector<ftm::idNode> children;
    mTree.tree.getChildren(node, children);
    for(auto &child : children)
      queue.emplace(child);
  }
  ftm::setTreeScalars<float>(mTree, scalars);
}

void ttk::MergeTreeNeuralBase::printPairs(const ftm::MergeTree<float> &mTree,
                                          bool useBD) {
  std::stringstream ss;
  if(mTree.tree.getRealNumberOfNodes() != 0)
    ss = mTree.tree.template printPairsFromTree<float>(useBD);
  else {
    std::vector<bool> nodeDone(mTree.tree.getNumberOfNodes(), false);
    for(unsigned int i = 0; i < mTree.tree.getNumberOfNodes(); ++i) {
      if(nodeDone[i])
        continue;
      std::tuple<ftm::idNode, ftm::idNode, float> pair
        = std::make_tuple(i, mTree.tree.getNode(i)->getOrigin(),
                          mTree.tree.getNodePersistence<float>(i));
      ss << std::get<0>(pair) << " ("
         << mTree.tree.getValue<float>(std::get<0>(pair)) << ") _ ";
      ss << std::get<1>(pair) << " ("
         << mTree.tree.getValue<float>(std::get<1>(pair)) << ") _ ";
      ss << std::get<2>(pair) << std::endl;
      nodeDone[i] = true;
      nodeDone[mTree.tree.getNode(i)->getOrigin()] = true;
    }
  }
  ss << std::endl;
  std::cout << ss.str();
}

//  -----------------------------------------------------------------------
//  --- Distance
//  -----------------------------------------------------------------------
void ttk::MergeTreeNeuralBase::getDistanceMatrix(
  const std::vector<mtu::TorchMergeTree<float>> &tmts,
  std::vector<std::vector<float>> &distanceMatrix,
  bool useDoubleInput,
  bool isFirstInput) {
  distanceMatrix.clear();
  distanceMatrix.resize(tmts.size(), std::vector<float>(tmts.size(), 0));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) if(parallelize_) \
  shared(distanceMatrix, tmts)
  {
#pragma omp single nowait
    {
#endif
      for(unsigned int i = 0; i < tmts.size(); ++i) {
        for(unsigned int j = i + 1; j < tmts.size(); ++j) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task UNTIED() shared(distanceMatrix, tmts) firstprivate(i, j)
          {
#endif
            std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching;
            float distance;
            bool isCalled = true;
            computeOneDistance(tmts[i].mTree, tmts[j].mTree, matching, distance,
                               isCalled, useDoubleInput, isFirstInput);
            distance = distance * distance;
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
#ifdef TTK_ENABLE_OPENMP
          } // pragma omp task
#endif
        }
      }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
    } // pragma omp single nowait
  } // pragma omp parallel
#endif
}

void ttk::MergeTreeNeuralBase::getDistanceMatrix(
  const std::vector<mtu::TorchMergeTree<float>> &tmts,
  const std::vector<mtu::TorchMergeTree<float>> &tmts2,
  std::vector<std::vector<float>> &distanceMatrix) {
  getDistanceMatrix(tmts, distanceMatrix, useDoubleInput_);
  if(useDoubleInput_) {
    std::vector<std::vector<float>> distanceMatrix2;
    getDistanceMatrix(tmts2, distanceMatrix2, useDoubleInput_, false);
    mixDistancesMatrix<float>(distanceMatrix, distanceMatrix2);
  }
}

void ttk::MergeTreeNeuralBase::getDifferentiableDistanceFromMatchings(
  const mtu::TorchMergeTree<float> &tree1,
  const mtu::TorchMergeTree<float> &tree2,
  const mtu::TorchMergeTree<float> &tree1_2,
  const mtu::TorchMergeTree<float> &tree2_2,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings,
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings2,
  torch::Tensor &tensorDist,
  bool doSqrt) {
  torch::Tensor reorderedITensor, reorderedJTensor;
  dataReorderingGivenMatching(
    tree1, tree2, matchings, reorderedITensor, reorderedJTensor);
  if(useDoubleInput_) {
    torch::Tensor reorderedI2Tensor, reorderedJ2Tensor;
    dataReorderingGivenMatching(
      tree1_2, tree2_2, matchings2, reorderedI2Tensor, reorderedJ2Tensor);
    reorderedITensor = torch::cat({reorderedITensor, reorderedI2Tensor});
    reorderedJTensor = torch::cat({reorderedJTensor, reorderedJ2Tensor});
  }
  tensorDist = (reorderedITensor - reorderedJTensor).pow(2).sum();
  if(doSqrt)
    tensorDist = tensorDist.sqrt();
}

void ttk::MergeTreeNeuralBase::getDifferentiableDistance(
  const mtu::TorchMergeTree<float> &tree1,
  const mtu::TorchMergeTree<float> &tree2,
  const mtu::TorchMergeTree<float> &tree1_2,
  const mtu::TorchMergeTree<float> &tree2_2,
  torch::Tensor &tensorDist,
  bool isCalled,
  bool doSqrt) {
  std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matchings,
    matchings2;
  float distance;
  computeOneDistance<float>(
    tree1.mTree, tree2.mTree, matchings, distance, isCalled, useDoubleInput_);
  if(useDoubleInput_) {
    float distance2;
    computeOneDistance<float>(tree1_2.mTree, tree2_2.mTree, matchings2,
                              distance2, isCalled, useDoubleInput_, false);
  }
  getDifferentiableDistanceFromMatchings(
    tree1, tree2, tree1_2, tree2_2, matchings, matchings2, tensorDist, doSqrt);
}

void ttk::MergeTreeNeuralBase::getDifferentiableDistance(
  const mtu::TorchMergeTree<float> &tree1,
  const mtu::TorchMergeTree<float> &tree2,
  torch::Tensor &tensorDist,
  bool isCalled,
  bool doSqrt) {
  mtu::TorchMergeTree<float> tree1_2, tree2_2;
  getDifferentiableDistance(
    tree1, tree2, tree1_2, tree2_2, tensorDist, isCalled, doSqrt);
}

void ttk::MergeTreeNeuralBase::getDifferentiableDistanceMatrix(
  const std::vector<mtu::TorchMergeTree<float> *> &trees,
  const std::vector<mtu::TorchMergeTree<float> *> &trees2,
  std::vector<std::vector<torch::Tensor>> &outDistMat) {
  outDistMat.resize(trees.size(), std::vector<torch::Tensor>(trees.size()));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) if(parallelize_) \
  shared(trees, trees2, outDistMat)
  {
#pragma omp single nowait
    {
#endif
      for(unsigned int i = 0; i < trees.size(); ++i) {
        outDistMat[i][i] = torch::tensor(0);
        for(unsigned int j = i + 1; j < trees.size(); ++j) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task UNTIED() shared(trees, trees2, outDistMat) firstprivate(i, j)
          {
#endif
            bool isCalled = true;
            bool doSqrt = false;
            torch::Tensor tensorDist;
            getDifferentiableDistance(*(trees[i]), *(trees[j]), *(trees2[i]),
                                      *(trees2[j]), tensorDist, isCalled,
                                      doSqrt);
            outDistMat[i][j] = tensorDist;
            outDistMat[j][i] = tensorDist;
#ifdef TTK_ENABLE_OPENMP
          } // pragma omp task
#endif
        }
      }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
    } // pragma omp single nowait
  } // pragma omp parallel
#endif
}
#endif
