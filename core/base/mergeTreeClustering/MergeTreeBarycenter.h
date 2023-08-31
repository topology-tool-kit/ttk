/// \ingroup base
/// \class MergeTreeBarycenter
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// This module defines the %MergeTreeBarycenter class that computes
/// the barycenter of an ensemble of merge trees.
///
/// \b Related \b publication \n
/// "Wasserstein Distances, Geodesics and Barycenters of Merge Trees" \n
/// Mathieu Pont, Jules Vidal, Julie Delon, Julien Tierny.\n
/// Proc. of IEEE VIS 2021.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021

#pragma once

#include <random>

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include "MergeTreeBase.h"
#include "MergeTreeDistance.h"

namespace ttk {

  /**
   * The MergeTreeBarycenter class that computes
   * the barycenter of an ensemble of merge trees.
   */
  class MergeTreeBarycenter : virtual public Debug, public MergeTreeBase {

  protected:
    double tol_ = 0.0;
    bool addNodes_ = true;
    bool deterministic_ = true;
    bool isCalled_ = false;
    bool progressiveBarycenter_ = false;
    double progressiveSpeedDivisor_ = 4.0;
    double alpha_ = 0.5;
    unsigned int barycenterMaximumNumberOfPairs_ = 0;
    double barycenterSizeLimitPercent_ = 0.0;

    double allDistanceTime_ = 0;

    double addDeletedNodesTime_ = 0;

    bool preprocess_ = true;
    bool postprocess_ = true;

    // Output
    std::vector<double> finalDistances_;

  public:
    MergeTreeBarycenter() {
      this->setDebugMsgPrefix(
        "MergeTreeBarycenter"); // inherited from Debug: prefix will be printed
                                // at the beginning of every msg
#ifdef TTK_ENABLE_OPENMP4
      omp_set_nested(1);
#endif
    }
    ~MergeTreeBarycenter() override = default;

    void setTol(double tolT) {
      tol_ = tolT;
    }

    void setAddNodes(bool addNodesT) {
      addNodes_ = addNodesT;
    }

    void setDeterministic(bool deterministicT) {
      deterministic_ = deterministicT;
    }

    void setProgressiveBarycenter(bool progressive) {
      progressiveBarycenter_ = progressive;
    }

    void setProgressiveSpeedDivisor(double progSpeed) {
      progressiveSpeedDivisor_ = progSpeed;
    }

    void setIsCalled(bool ic) {
      isCalled_ = ic;
    }

    double getAllDistanceTime() {
      return allDistanceTime_;
    }

    double getAddDeletedNodesTime() {
      return addDeletedNodesTime_;
    }

    void setAlpha(double alpha) {
      alpha_ = alpha;
    }

    void setBarycenterMaximumNumberOfPairs(unsigned int maxi) {
      barycenterMaximumNumberOfPairs_ = maxi;
    }

    void setBarycenterSizeLimitPercent(double percent) {
      barycenterSizeLimitPercent_ = percent;
    }

    void setPreprocess(bool preproc) {
      preprocess_ = preproc;
    }

    void setPostprocess(bool postproc) {
      postprocess_ = postproc;
    }

    std::vector<double> getFinalDistances() {
      return finalDistances_;
    }

    /**
     * Implementation of the algorithm.
     */
    // ------------------------------------------------------------------------
    // Initialization
    // ------------------------------------------------------------------------
    template <class dataType>
    void getDistanceMatrix(std::vector<ftm::FTMTree_MT *> &trees,
                           std::vector<ftm::FTMTree_MT *> &trees2,
                           std::vector<std::vector<double>> &distanceMatrix,
                           bool useDoubleInput = false,
                           bool isFirstInput = true) {
      distanceMatrix.clear();
      distanceMatrix.resize(trees.size(), std::vector<double>(trees.size(), 0));
#ifdef TTK_ENABLE_OPENMP4
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
      for(unsigned int i = 0; i < trees.size(); ++i)
        for(unsigned int j = i + 1; j < trees.size(); ++j) {
          std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching;
          dataType distance;
          computeOneDistance<dataType>(trees[i], trees2[j], matching, distance,
                                       useDoubleInput, isFirstInput);
          distanceMatrix[i][j] = distance;
          distanceMatrix[j][i] = distance;
        }
    }

    template <class dataType>
    void getDistanceMatrix(std::vector<ftm::FTMTree_MT *> &trees,
                           std::vector<std::vector<double>> &distanceMatrix,
                           bool useDoubleInput = false,
                           bool isFirstInput = true) {
      getDistanceMatrix<dataType>(
        trees, trees, distanceMatrix, useDoubleInput, isFirstInput);
    }

    template <class dataType>
    void getSizeLimitedTrees(
      std::vector<ftm::FTMTree_MT *> &trees,
      unsigned int barycenterMaximumNumberOfPairs,
      double sizeLimitPercent,
      std::vector<ftm::MergeTree<dataType>> &mTreesLimited) {
      mTreesLimited.resize(trees.size());
      for(unsigned int i = 0; i < trees.size(); ++i) {
        mTreesLimited[i] = ftm::copyMergeTree<dataType>(trees[i]);
        limitSizeBarycenter(mTreesLimited[i], trees,
                            barycenterMaximumNumberOfPairs, sizeLimitPercent);
        ftm::cleanMergeTree<dataType>(mTreesLimited[i]);
      }
    }

    template <class dataType>
    void getSizeLimitedDistanceMatrix(
      std::vector<ftm::FTMTree_MT *> &trees,
      std::vector<std::vector<double>> &distanceMatrix,
      unsigned int barycenterMaximumNumberOfPairs,
      double sizeLimitPercent,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      std::vector<ftm::MergeTree<dataType>> mTreesLimited;
      getSizeLimitedTrees<dataType>(
        trees, barycenterMaximumNumberOfPairs, sizeLimitPercent, mTreesLimited);
      std::vector<ftm::FTMTree_MT *> treesLimited;
      ftm::mergeTreeToFTMTree<dataType>(mTreesLimited, treesLimited);
      getDistanceMatrix<dataType>(
        trees, treesLimited, distanceMatrix, useDoubleInput, isFirstInput);
    }

    template <class dataType>
    void getParametrizedDistanceMatrix(
      std::vector<ftm::FTMTree_MT *> &trees,
      std::vector<std::vector<double>> &distanceMatrix,
      unsigned int barycenterMaximumNumberOfPairs,
      double sizeLimitPercent,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      if(barycenterMaximumNumberOfPairs <= 0 and sizeLimitPercent <= 0.0)
        getDistanceMatrix<dataType>(
          trees, distanceMatrix, useDoubleInput, isFirstInput);
      else
        getSizeLimitedDistanceMatrix<dataType>(
          trees, distanceMatrix, barycenterMaximumNumberOfPairs,
          sizeLimitPercent, useDoubleInput, isFirstInput);
    }

    template <class dataType>
    int getBestInitTreeIndex(std::vector<ftm::FTMTree_MT *> &trees,
                             std::vector<ftm::FTMTree_MT *> &trees2,
                             unsigned int barycenterMaximumNumberOfPairs,
                             double sizeLimitPercent,
                             bool distMinimizer = true) {
      std::vector<std::vector<double>> distanceMatrix, distanceMatrix2;
      bool const useDoubleInput = (trees2.size() != 0);
      getParametrizedDistanceMatrix<dataType>(trees, distanceMatrix,
                                              barycenterMaximumNumberOfPairs,
                                              sizeLimitPercent, useDoubleInput);
      if(trees2.size() != 0)
        getParametrizedDistanceMatrix<dataType>(
          trees2, distanceMatrix2, barycenterMaximumNumberOfPairs,
          sizeLimitPercent, useDoubleInput, false);

      int bestIndex = -1;
      dataType bestValue
        = distMinimizer ? std::numeric_limits<dataType>::max() : 0;
      std::vector<int> sizes(trees.size());
      for(unsigned int i = 0; i < trees.size(); ++i) {
        dataType value = 0;
        for(unsigned int j = 0; j < distanceMatrix[i].size(); ++j)
          value += (not useDoubleInput ? distanceMatrix[i][j]
                                       : mixDistances(distanceMatrix[i][j],
                                                      distanceMatrix2[i][j]));
        if((distMinimizer and value < bestValue)
           or (not distMinimizer and value > bestValue)) {
          bestIndex = i;
          bestValue = value;
        }
        sizes[i] = -value;
        sizes[i] *= (distMinimizer) ? 1 : -1;
      }
      if(not deterministic_) {
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::discrete_distribution<int> distribution(
          sizes.begin(), sizes.end());
        bestIndex = distribution(generator);
      }
      return bestIndex;
    }

    template <class dataType>
    int getBestInitTreeIndex(std::vector<ftm::FTMTree_MT *> &trees,
                             std::vector<ftm::FTMTree_MT *> &trees2,
                             double sizeLimitPercent,
                             bool distMinimizer = true) {
      return getBestInitTreeIndex<dataType>(trees, trees2,
                                            barycenterMaximumNumberOfPairs_,
                                            sizeLimitPercent, distMinimizer);
    }

    template <class dataType>
    int getBestInitTreeIndex(std::vector<ftm::FTMTree_MT *> &trees,
                             bool distMinimizer = true) {
      std::vector<ftm::FTMTree_MT *> trees2;
      return getBestInitTreeIndex<dataType>(
        trees, trees2, barycenterMaximumNumberOfPairs_,
        barycenterSizeLimitPercent_, distMinimizer);
    }

    template <class dataType>
    void initBarycenterTree(std::vector<ftm::FTMTree_MT *> &trees,
                            ftm::MergeTree<dataType> &baryTree,
                            bool distMinimizer = true) {
      int const bestIndex
        = getBestInitTreeIndex<dataType>(trees, distMinimizer);
      baryTree = ftm::copyMergeTree<dataType>(trees[bestIndex], true);
      limitSizeBarycenter(baryTree, trees);
    }

    // ------------------------------------------------------------------------
    // Update
    // ------------------------------------------------------------------------
    template <class dataType>
    ftm::idNode getNodesAndScalarsToAdd(
      ftm::MergeTree<dataType> &ttkNotUsed(mTree1),
      ftm::idNode nodeId1,
      ftm::FTMTree_MT *tree2,
      ftm::idNode nodeId2,
      std::vector<dataType> &newScalarsVector,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> &nodesToProcess,
      ftm::idNode nodeCpt,
      int i) {
      // Get nodes and scalars to add
      std::queue<std::tuple<ftm::idNode, ftm::idNode>> queue;
      queue.emplace(nodeId2, nodeId1);
      nodesToProcess.emplace_back(nodeId2, nodeId1, i);
      while(!queue.empty()) {
        auto queueTuple = queue.front();
        queue.pop();
        ftm::idNode const node = std::get<0>(queueTuple);
        // Get scalars
        newScalarsVector.push_back(
          tree2->getValue<dataType>(tree2->getNode(node)->getOrigin()));
        newScalarsVector.push_back(tree2->getValue<dataType>(node));
        // Process children
        std::vector<ftm::idNode> children;
        tree2->getChildren(node, children);
        for(auto child : children) {
          queue.emplace(child, nodeCpt + 1);
          nodesToProcess.emplace_back(child, nodeCpt + 1, i);
        }
        nodeCpt += 2; // we will add two nodes (birth and death)
      }

      return nodeCpt;
    }

    template <class dataType>
    void addNodes(
      ftm::MergeTree<dataType> &mTree1,
      int noTrees,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> &nodesToProcess,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>
        &nodesProcessed) {
      ftm::FTMTree_MT *tree1 = &(mTree1.tree);

      // Add nodes
      nodesProcessed.clear();
      nodesProcessed.resize(noTrees);
      for(auto processTuple : nodesToProcess) {
        ftm::idNode const parent = std::get<1>(processTuple);
        ftm::idNode const nodeTree1 = tree1->getNumberOfNodes();
        int const index = std::get<2>(processTuple);
        nodesProcessed[index].emplace_back(
          nodeTree1 + 1, std::get<0>(processTuple));
        // Make node and its origin
        tree1->makeNode(nodeTree1);
        tree1->makeNode(nodeTree1 + 1);
        tree1->setParent(nodeTree1 + 1, parent);
        tree1->getNode(nodeTree1)->setOrigin(nodeTree1 + 1);
        tree1->getNode(nodeTree1 + 1)->setOrigin(nodeTree1);
      }
    }

    template <class dataType>
    void updateNodesAndScalars(
      ftm::MergeTree<dataType> &mTree1,
      int noTrees,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> &nodesToProcess,
      std::vector<dataType> &newScalarsVector,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>
        &nodesProcessed) {
      ftm::FTMTree_MT *tree1 = &(mTree1.tree);

      // Create new tree
      ftm::MergeTree<dataType> mTreeNew
        = ftm::createEmptyMergeTree<dataType>(newScalarsVector.size());
      ftm::setTreeScalars<dataType>(mTreeNew, newScalarsVector);
      ftm::FTMTree_MT *treeNew = &(mTreeNew.tree);

      // Copy the old tree structure
      treeNew->copyMergeTreeStructure(tree1);

      // Add nodes in the other trees
      addNodes<dataType>(mTreeNew, noTrees, nodesToProcess, nodesProcessed);

      // Copy new tree
      mTree1 = mTreeNew;
    }

    template <class dataType>
    void updateBarycenterTreeStructure(
      std::vector<ftm::FTMTree_MT *> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings) {
      ftm::FTMTree_MT *baryTree = &(baryMergeTree.tree);
      ftm::idNode const baryTreeRoot = baryTree->getRoot();

      // Init matching matrix
      // m[i][j] contains the node in the barycenter matched to the jth node of
      // the ith tree
      std::vector<std::vector<ftm::idNode>> matrixMatchings(trees.size());
      std::vector<bool> baryMatched(baryTree->getNumberOfNodes(), false);
      for(unsigned int i = 0; i < matchings.size(); ++i) {
        auto matching = matchings[i];
        matrixMatchings[i].resize(trees[i]->getNumberOfNodes(),
                                  std::numeric_limits<ftm::idNode>::max());
        for(auto match : matching) {
          matrixMatchings[i][std::get<1>(match)] = std::get<0>(match);
          baryMatched[std::get<0>(match)] = true;
        }
      }

      // Iterate through trees to get the nodes to add in the barycenter
      std::vector<std::vector<ftm::idNode>> nodesToAdd(trees.size());
      for(unsigned int i = 0; i < trees.size(); ++i) {
        ftm::idNode const root = trees[i]->getRoot();
        std::queue<ftm::idNode> queue;
        queue.emplace(root);
        while(!queue.empty()) {
          ftm::idNode const node = queue.front();
          queue.pop();
          bool processChildren = true;
          // if node in trees[i] is not matched
          if(matrixMatchings[i][node]
             == std::numeric_limits<ftm::idNode>::max()) {
            if(not keepSubtree_) {
              processChildren = false;
              nodesToAdd[i].push_back(node);
            } else {
              // not todo manage if keepSubtree=true (not important since it is
              // not a valid merge tree)
              printErr(
                "barycenter with keepSubtree_=true is not implemented yet");
            }
          }
          if(processChildren) {
            std::vector<ftm::idNode> children;
            trees[i]->getChildren(node, children);
            for(auto child : children)
              if(not(trees[i]->isThereOnlyOnePersistencePair()
                     and trees[i]->isLeaf(child)))
                queue.emplace(child);
          }
        }
      }

      bool foundRootNotMatched = false;
      for(unsigned int i = 0; i < trees.size(); ++i)
        foundRootNotMatched |= baryTree->isNodeIdInconsistent(
          matrixMatchings[i][trees[i]->getRoot()]);
      if(foundRootNotMatched)
        printWrn("[updateBarycenterTreeStructure] an input tree has its root "
                 "not matched.");

      // Delete nodes that are not matched in the barycenter
      for(unsigned int i = 0; i < baryTree->getNumberOfNodes(); ++i)
        if(not baryMatched[i])
          baryTree->deleteNode(i);

      if(not keepSubtree_) {
        // Add scalars and nodes not present in the barycenter
        ftm::idNode nodeCpt = baryTree->getNumberOfNodes();
        std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> nodesToProcess;
        std::vector<dataType> newScalarsVector;
        ftm::getTreeScalars<dataType>(baryMergeTree, newScalarsVector);
        for(unsigned int i = 0; i < nodesToAdd.size(); ++i) {
          for(auto node : nodesToAdd[i]) {
            ftm::idNode parent
              = matrixMatchings[i][trees[i]->getParentSafe(node)];
            if(matchings[i].size() == 0)
              parent = baryTreeRoot;

            if((baryTree->isNodeIdInconsistent(parent)
                or baryTree->isNodeAlone(parent))
               and matchings[i].size() != 0) {
              std::stringstream ss;
              ss << trees[i]->getParentSafe(node) << " _ " << node;
              printMsg(ss.str());
              printMsg(trees[i]->printTree().str());
              printMsg(trees[i]->printPairsFromTree<dataType>(true).str());
              printMatching(matchings[i]);
              std::stringstream ss2;
              ss2 << "parent " << parent;
              printMsg(ss2.str());
            }
            /*if(isRoot(trees[i], node))
              parent = baryTree->getRoot();*/
            std::vector<dataType> addedScalars;
            nodeCpt = getNodesAndScalarsToAdd<dataType>(
              baryMergeTree, parent, trees[i], node, addedScalars,
              nodesToProcess, nodeCpt, i);
            newScalarsVector.insert(
              newScalarsVector.end(), addedScalars.begin(), addedScalars.end());
          }
        }
        if(addNodes_) {
          std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>
            nodesProcessed;
          updateNodesAndScalars<dataType>(baryMergeTree, trees.size(),
                                          nodesToProcess, newScalarsVector,
                                          nodesProcessed);
          for(unsigned int i = 0; i < matchings.size(); ++i) {
            std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
              nodesProcessedT;
            for(auto tup : nodesProcessed[i])
              nodesProcessedT.emplace_back(
                std::get<0>(tup), std::get<1>(tup), -1);
            matchings[i].insert(matchings[i].end(), nodesProcessedT.begin(),
                                nodesProcessedT.end());
          }
        }
      } else {
        // not todo manage if keepSubtree=true (not important since it is not a
        // valid merge tree)
        printErr("barycenter with keepSubtree_=true is not implemented yet");
      }
    }

    template <class dataType>
    std::tuple<dataType, dataType>
      getParametrizedBirthDeath(ftm::FTMTree_MT *tree1, ftm::idNode nodeId1) {
      std::tuple<dataType, dataType> birthDeath;
      // Normalized Wasserstein
      if(normalizedWasserstein_)
        birthDeath = getNormalizedBirthDeathDouble<dataType>(tree1, nodeId1);
      // Classical Wasserstein
      else
        birthDeath = tree1->getBirthDeath<dataType>(nodeId1);
      return birthDeath;
    }

    template <class dataType>
    std::tuple<dataType, dataType>
      interpolation(ftm::MergeTree<dataType> &baryMergeTree,
                    ftm::idNode nodeId,
                    std::vector<dataType> &newScalarsVector,
                    std::vector<ftm::FTMTree_MT *> &trees,
                    std::vector<ftm::idNode> &nodes,
                    std::vector<double> &alphas) {
      ftm::FTMTree_MT *baryTree = &(baryMergeTree.tree);
      dataType mu_max = getMinMaxLocalFromVector<dataType>(
        baryTree, nodeId, newScalarsVector, false);
      dataType mu_min = getMinMaxLocalFromVector<dataType>(
        baryTree, nodeId, newScalarsVector);
      double newBirth = 0, newDeath = 0;

      // Compute projection
      double tempBirth = 0, tempDeath = 0;
      double alphaSum = 0;
      for(unsigned int i = 0; i < trees.size(); ++i)
        if(nodes[i] != std::numeric_limits<ftm::idNode>::max())
          alphaSum += alphas[i];
      for(unsigned int i = 0; i < trees.size(); ++i) {
        // if node is matched in trees[i]
        if(nodes[i] != std::numeric_limits<ftm::idNode>::max()) {
          auto iBirthDeath
            = getParametrizedBirthDeath<dataType>(trees[i], nodes[i]);
          double tTempBirth = 0, tTempDeath = 0;
          tTempBirth += std::get<0>(iBirthDeath);
          tTempDeath += std::get<1>(iBirthDeath);
          tempBirth += tTempBirth * alphas[i] / alphaSum;
          tempDeath += tTempDeath * alphas[i] / alphaSum;
        }
      }
      double const projec = (tempBirth + tempDeath) / 2;

      // Compute newBirth and newDeath
      for(unsigned int i = 0; i < trees.size(); ++i) {
        double iBirth = projec, iDeath = projec;
        // if node is matched in trees[i]
        if(nodes[i] != std::numeric_limits<ftm::idNode>::max()) {
          auto iBirthDeath
            = getParametrizedBirthDeath<dataType>(trees[i], nodes[i]);
          iBirth = std::get<0>(iBirthDeath);
          iDeath = std::get<1>(iBirthDeath);
        }
        newBirth += alphas[i] * iBirth;
        newDeath += alphas[i] * iDeath;
      }
      if(normalizedWasserstein_) {
        newBirth = newBirth * (mu_max - mu_min) + mu_min;
        newDeath = newDeath * (mu_max - mu_min) + mu_min;
      }

      return std::make_tuple(newBirth, newDeath);
    }

    template <class dataType>
    std::tuple<dataType, dataType>
      interpolationAdded(ftm::FTMTree_MT *tree,
                         ftm::idNode nodeId,
                         double alpha,
                         ftm::MergeTree<dataType> &baryMergeTree,
                         ftm::idNode nodeB,
                         std::vector<dataType> &newScalarsVector) {
      ftm::FTMTree_MT *baryTree = &(baryMergeTree.tree);
      dataType mu_max = getMinMaxLocalFromVector<dataType>(
        baryTree, nodeB, newScalarsVector, false);
      dataType mu_min
        = getMinMaxLocalFromVector<dataType>(baryTree, nodeB, newScalarsVector);

      auto birthDeath = getParametrizedBirthDeath<dataType>(tree, nodeId);
      double newBirth = std::get<0>(birthDeath);
      double newDeath = std::get<1>(birthDeath);
      double const projec = (newBirth + newDeath) / 2;

      newBirth = alpha * newBirth + (1 - alpha) * projec;
      newDeath = alpha * newDeath + (1 - alpha) * projec;

      if(normalizedWasserstein_) {
        newBirth = newBirth * (mu_max - mu_min) + mu_min;
        newDeath = newDeath * (mu_max - mu_min) + mu_min;
      }

      dataType newBirthT = newBirth;
      dataType newDeathT = newDeath;
      return std::make_tuple(newBirthT, newDeathT);
    }

    template <class dataType>
    void updateBarycenterTreeScalars(
      std::vector<ftm::FTMTree_MT *> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<double> &alphas,
      unsigned int indexAddedNodes,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings) {
      ftm::FTMTree_MT *baryTree = &(baryMergeTree.tree);
      bool const isJT = baryTree->isJoinTree<dataType>();

      // Init matching matrix
      // m[i][j] contains the node in trees[j] matched to the node i in the
      // barycenter
      std::vector<std::vector<ftm::idNode>> baryMatching(
        baryTree->getNumberOfNodes(),
        std::vector<ftm::idNode>(
          trees.size(), std::numeric_limits<ftm::idNode>::max()));
      std::vector<int> nodesAddedTree(baryTree->getNumberOfNodes(), -1);
      for(unsigned int i = 0; i < matchings.size(); ++i) {
        auto matching = matchings[i];
        for(auto match : matching) {
          baryMatching[std::get<0>(match)][i] = std::get<1>(match);
          if(std::get<0>(match)
             >= indexAddedNodes) // get the tree of this added node
            nodesAddedTree[std::get<0>(match)] = i;
        }
      }

      // Interpolate scalars
      std::vector<dataType> newScalarsVector(baryTree->getNumberOfNodes());
      ftm::idNode const root = baryTree->getRoot();
      std::queue<ftm::idNode> queue;
      queue.emplace(root);
      while(!queue.empty()) {
        ftm::idNode const node = queue.front();
        queue.pop();
        std::tuple<dataType, dataType> newBirthDeath;
        if(node < indexAddedNodes) {
          newBirthDeath
            = interpolation<dataType>(baryMergeTree, node, newScalarsVector,
                                      trees, baryMatching[node], alphas);
        } else {
          int const i = nodesAddedTree[node];
          ftm::idNode const nodeT = baryMatching[node][i];
          newBirthDeath = interpolationAdded<dataType>(
            trees[i], nodeT, alphas[i], baryMergeTree, node, newScalarsVector);
        }
        dataType nodeScalar
          = (isJT ? std::get<1>(newBirthDeath) : std::get<0>(newBirthDeath));
        dataType nodeOriginScalar
          = (isJT ? std::get<0>(newBirthDeath) : std::get<1>(newBirthDeath));
        newScalarsVector[node] = nodeScalar;
        newScalarsVector[baryTree->getNode(node)->getOrigin()]
          = nodeOriginScalar;
        std::vector<ftm::idNode> children;
        baryTree->getChildren(node, children);
        for(auto child : children)
          queue.emplace(child);
      }

      if(baryMergeTree.tree.isFullMerge()) {
        auto mergedRootOrigin = baryTree->getMergedRootOrigin<dataType>();
        dataType mergedRootOriginScalar = 0.0;
        for(unsigned int i = 0; i < trees.size(); ++i)
          mergedRootOriginScalar += trees[i]->getValue<dataType>(
            trees[i]->getMergedRootOrigin<dataType>());
        mergedRootOriginScalar /= trees.size();
        newScalarsVector[mergedRootOrigin] = mergedRootOriginScalar;
      }

      setTreeScalars(baryMergeTree, newScalarsVector);

      std::vector<ftm::idNode> deletedNodesT;
      persistenceThresholding<dataType>(
        &(baryMergeTree.tree), 0, deletedNodesT);
      limitSizeBarycenter(baryMergeTree, trees);
      ftm::cleanMergeTree<dataType>(baryMergeTree);
    }

    template <class dataType>
    void updateBarycenterTree(
      std::vector<ftm::FTMTree_MT *> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<double> &alphas,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings) {
      int const indexAddedNodes = baryMergeTree.tree.getNumberOfNodes();
      updateBarycenterTreeStructure<dataType>(trees, baryMergeTree, matchings);
      updateBarycenterTreeScalars<dataType>(
        trees, baryMergeTree, alphas, indexAddedNodes, matchings);
    }

    // ------------------------------------------------------------------------
    // Assignment
    // ------------------------------------------------------------------------
    template <class dataType>
    void computeOneDistance(
      ftm::FTMTree_MT *tree,
      ftm::FTMTree_MT *baryTree,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      dataType &distance,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      // Timer t_distance;
      MergeTreeDistance mergeTreeDistance;
      mergeTreeDistance.setDebugLevel(std::min(debugLevel_, 2));
      mergeTreeDistance.setPreprocess(false);
      mergeTreeDistance.setPostprocess(false);
      mergeTreeDistance.setBranchDecomposition(true);
      mergeTreeDistance.setNormalizedWasserstein(normalizedWasserstein_);
      mergeTreeDistance.setKeepSubtree(keepSubtree_);
      mergeTreeDistance.setAssignmentSolver(assignmentSolverID_);
      mergeTreeDistance.setIsCalled(true);
      mergeTreeDistance.setThreadNumber(this->threadNumber_);
      mergeTreeDistance.setDistanceSquaredRoot(true); // squared root
      mergeTreeDistance.setNodePerTask(nodePerTask_);
      if(useDoubleInput) {
        double const weight = mixDistancesMinMaxPairWeight(isFirstInput);
        mergeTreeDistance.setMinMaxPairWeight(weight);
      }
      /*if(progressiveBarycenter_){
        mergeTreeDistance.setAuctionNoRounds(1);
        mergeTreeDistance.setAuctionEpsilonDiviser(NoIteration-1);
      }*/
      distance
        = mergeTreeDistance.computeDistance<dataType>(baryTree, tree, matching);
      std::stringstream ss, ss2;
      ss << "distance tree : " << distance;
      printMsg(ss.str(), debug::Priority::VERBOSE);
      ss2 << "distance²tree : " << distance * distance;
      printMsg(ss2.str(), debug::Priority::VERBOSE);

      // auto t_distance_time = t_distance.getElapsedTime();
      // allDistanceTime_ += t_distance_time;
    }

    template <class dataType>
    void computeOneDistance(
      ftm::FTMTree_MT *tree,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      dataType &distance,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      computeOneDistance<dataType>(tree, &(baryMergeTree.tree), matching,
                                   distance, useDoubleInput, isFirstInput);
    }

    template <class dataType>
    void computeOneDistance(
      ftm::MergeTree<dataType> &baryMergeTree,
      ftm::MergeTree<dataType> &baryMergeTree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      dataType &distance,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      computeOneDistance<dataType>(&(baryMergeTree.tree), baryMergeTree2,
                                   matching, distance, useDoubleInput,
                                   isFirstInput);
    }

    template <class dataType>
    void assignment(
      std::vector<ftm::FTMTree_MT *> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<dataType> &distances,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      if(not isCalled_)
        assignmentPara(trees, baryMergeTree, matchings, distances,
                       useDoubleInput, isFirstInput);
      else
        assignmentTask(trees, baryMergeTree, matchings, distances,
                       useDoubleInput, isFirstInput);
    }

    template <class dataType>
    void assignmentPara(
      std::vector<ftm::FTMTree_MT *> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<dataType> &distances,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
#ifdef TTK_ENABLE_OPENMP4
#pragma omp parallel num_threads(this->threadNumber_) \
  shared(baryMergeTree) if(parallelize_)
      {
#pragma omp single nowait
#endif
        assignmentTask(trees, baryMergeTree, matchings, distances,
                       useDoubleInput, isFirstInput);
#ifdef TTK_ENABLE_OPENMP4
      } // pragma omp parallel
#endif
    }

    template <class dataType>
    void assignmentTask(
      std::vector<ftm::FTMTree_MT *> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<dataType> &distances,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      for(unsigned int i = 0; i < trees.size(); ++i)
#ifdef TTK_ENABLE_OPENMP4
#pragma omp task firstprivate(i) UNTIED() \
  shared(baryMergeTree, matchings, distances)
#endif
        computeOneDistance<dataType>(trees[i], baryMergeTree, matchings[i],
                                     distances[i], useDoubleInput,
                                     isFirstInput);
#ifdef TTK_ENABLE_OPENMP4
#pragma omp taskwait
#endif
    }

    // ------------------------------------------------------------------------
    // Progressivity
    // ------------------------------------------------------------------------
    template <class dataType>
    unsigned int
      persistenceScaling(std::vector<ftm::FTMTree_MT *> &trees,
                         std::vector<ftm::MergeTree<dataType>> &mergeTrees,
                         std::vector<ftm::FTMTree_MT *> &oriTrees,
                         int iterationNumber,
                         std::vector<std::vector<ftm::idNode>> &deletedNodes) {
      deletedNodes.clear();
      deletedNodes.resize(oriTrees.size());
      unsigned int noTreesUnscaled = 0;

      // Scale trees
      for(unsigned int i = 0; i < oriTrees.size(); ++i) {
        double persistenceThreshold = 50.0;
        if(iterationNumber != -1) {
          // Get number of pairs in scaled merge tree
          int const noPairs = mergeTrees[i].tree.getRealNumberOfNodes();

          // Get pairs in original merge tree
          std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> pairs;
          oriTrees[i]->getPersistencePairsFromTree<dataType>(
            pairs, branchDecomposition_);

          // Compute new persistence threshold
          double const multiplier
            = (progressiveSpeedDivisor_ < 1e-6
                 ? 1.
                 : iterationNumber / progressiveSpeedDivisor_);
          int const decrement = multiplier * pairs.size() / 10;
          int thresholdIndex = pairs.size() - noPairs - std::max(decrement, 2);
          thresholdIndex = std::max(thresholdIndex, 0);
          const double persistence = std::get<2>(pairs[thresholdIndex]);
          persistenceThreshold
            = persistence / std::get<2>(pairs.back()) * 100.0;
          if(thresholdIndex == 0) {
            persistenceThreshold = 0.;
            ++noTreesUnscaled;
          }
        }
        if(persistenceThreshold != 0.) {
          ftm::MergeTree<dataType> mt
            = ftm::copyMergeTree<dataType>(oriTrees[i]);
          persistenceThresholding<dataType>(
            &(mt.tree), persistenceThreshold, deletedNodes[i]);
          if(mergeTrees.size() == 0)
            mergeTrees.resize(oriTrees.size());
          mergeTrees[i] = mt;
          trees[i] = &(mt.tree);
        } else {
          trees[i] = oriTrees[i];
        }
      }

      printTreesStats(trees);

      return noTreesUnscaled;
    }

    template <class dataType>
    void addScaledDeletedNodesCost(
      std::vector<ftm::FTMTree_MT *> &oriTrees,
      std::vector<std::vector<ftm::idNode>> &deletedNodes,
      std::vector<dataType> &distances) {
      for(unsigned int i = 0; i < oriTrees.size(); ++i)
        for(auto node : deletedNodes[i])
          distances[i] += deleteCost<dataType>(oriTrees[i], node);
    }

    // ------------------------------------------------------------------------
    // Main Functions
    // ------------------------------------------------------------------------
    template <class dataType>
    void computeBarycenter(
      std::vector<ftm::FTMTree_MT *> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<double> &alphas,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &finalMatchings,
      bool finalAsgnDoubleInput = false,
      bool finalAsgnFirstInput = true) {
      Timer t_bary;

      ftm::FTMTree_MT *baryTree = &(baryMergeTree.tree);

      // Persistence scaling
      std::vector<ftm::FTMTree_MT *> oriTrees;
      std::vector<ftm::MergeTree<dataType>> scaledMergeTrees;
      std::vector<std::vector<ftm::idNode>> deletedNodes;
      if(progressiveBarycenter_) {
        oriTrees.insert(oriTrees.end(), trees.begin(), trees.end());
        persistenceScaling<dataType>(
          trees, scaledMergeTrees, oriTrees, -1, deletedNodes);
        std::vector<ftm::idNode> deletedNodesT;
        persistenceThresholding<dataType>(baryTree, 50, deletedNodesT);
      }
      bool treesUnscaled = false;

      // Print bary stats
      printBaryStats(baryTree);

      // Run
      bool converged = false;
      dataType frechetEnergy = -1;
      dataType minFrechet = std::numeric_limits<dataType>::max();
      int cptBlocked = 0;
      int NoIteration = 0;
      while(not converged) {
        ++NoIteration;

        printMsg(debug::Separator::L2);
        std::stringstream ss;
        ss << "Iteration " << NoIteration;
        printMsg(ss.str());

        // --- Assignment
        std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
          matchings(trees.size());
        std::vector<dataType> distances(trees.size(), -1);
        Timer t_assignment;
        assignment<dataType>(trees, baryMergeTree, matchings, distances);
        Timer t_addDeletedNodes;
        if(progressiveBarycenter_)
          addScaledDeletedNodesCost<dataType>(
            oriTrees, deletedNodes, distances);
        addDeletedNodesTime_ += t_addDeletedNodes.getElapsedTime();
        auto t_assignment_time
          = t_assignment.getElapsedTime() - t_addDeletedNodes.getElapsedTime();
        printMsg("Assignment", 1, t_assignment_time, this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::INFO);

        // --- Update
        Timer t_update;
        updateBarycenterTree<dataType>(trees, baryMergeTree, alphas, matchings);
        auto t_update_time = t_update.getElapsedTime();
        baryTree = &(baryMergeTree.tree);
        printMsg("Update", 1, t_update_time, this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::INFO);

        // --- Check convergence
        dataType currentFrechetEnergy = 0;
        for(unsigned int i = 0; i < trees.size(); ++i)
          currentFrechetEnergy += alphas[i] * distances[i] * distances[i];
        auto frechetDiff
          = std::abs((double)(frechetEnergy - currentFrechetEnergy));
        converged = (frechetDiff <= tol_);
        converged = converged and (not progressiveBarycenter_ or treesUnscaled);
        frechetEnergy = currentFrechetEnergy;
        tol_ = frechetEnergy / 125.0;

        std::stringstream ss4;
        auto barycenterTime = t_bary.getElapsedTime() - addDeletedNodesTime_;
        printMsg("Total", 1, barycenterTime, this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::INFO);
        printBaryStats(baryTree);
        ss4 << "Frechet energy : " << frechetEnergy;
        printMsg(ss4.str());

        minFrechet = std::min(minFrechet, frechetEnergy);
        if(not converged and (not progressiveBarycenter_ or treesUnscaled)) {
          cptBlocked = (minFrechet < frechetEnergy) ? cptBlocked + 1 : 0;
          converged = (cptBlocked >= 10);
        }

        // --- Persistence scaling
        if(progressiveBarycenter_) {
          unsigned int const noTreesUnscaled = persistenceScaling<dataType>(
            trees, scaledMergeTrees, oriTrees, NoIteration, deletedNodes);
          treesUnscaled = (noTreesUnscaled == oriTrees.size());
        }
      }

      // Final processing
      printMsg(debug::Separator::L2);
      printMsg("Final assignment");

      std::vector<dataType> distances(trees.size(), -1);
      assignment<dataType>(trees, baryMergeTree, finalMatchings, distances,
                           finalAsgnDoubleInput, finalAsgnFirstInput);
      for(auto dist : distances)
        finalDistances_.push_back(dist);
      dataType currentFrechetEnergy = 0;
      for(unsigned int i = 0; i < trees.size(); ++i)
        currentFrechetEnergy += alphas[i] * distances[i] * distances[i];

      std::stringstream ss, ss2;
      ss << "Frechet energy : " << currentFrechetEnergy;
      printMsg(ss.str());
      auto barycenterTime = t_bary.getElapsedTime() - addDeletedNodesTime_;
      printMsg("Total", 1, barycenterTime, this->threadNumber_,
               debug::LineMode::NEW, debug::Priority::INFO);
      // std::cout << "Bary Distance Time = " << allDistanceTime_ << std::endl;

      if(trees.size() == 2 and not isCalled_)
        verifyBarycenterTwoTrees<dataType>(
          trees, baryMergeTree, finalMatchings, distances);

      // Persistence (un)scaling
      if(progressiveBarycenter_) {
        scaledMergeTrees.clear();
        trees.clear();
        trees.insert(trees.end(), oriTrees.begin(), oriTrees.end());
      }
    }

    template <class dataType>
    void execute(
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<double> &alphas,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &finalMatchings,
      ftm::MergeTree<dataType> &baryMergeTree,
      bool finalAsgnDoubleInput = false,
      bool finalAsgnFirstInput = true) {
      // --- Preprocessing
      if(preprocess_) {
        treesNodeCorr_.resize(trees.size());
        for(unsigned int i = 0; i < trees.size(); ++i)
          preprocessingPipeline<dataType>(trees[i], epsilonTree2_,
                                          epsilon2Tree2_, epsilon3Tree2_,
                                          branchDecomposition_, useMinMaxPair_,
                                          cleanTree_, treesNodeCorr_[i]);
        printTreesStats(trees);
      }

      // --- Init barycenter
      std::vector<ftm::FTMTree_MT *> treesT;
      ftm::mergeTreeToFTMTree<dataType>(trees, treesT);
      initBarycenterTree<dataType>(treesT, baryMergeTree);

      // --- Execute
      computeBarycenter<dataType>(treesT, baryMergeTree, alphas, finalMatchings,
                                  finalAsgnDoubleInput, finalAsgnFirstInput);

      // --- Postprocessing
      if(postprocess_) {
        std::vector<int> const allRealNodes(trees.size());
        for(unsigned int i = 0; i < trees.size(); ++i) {
          postprocessingPipeline<dataType>(treesT[i]);
        }

        // fixMergedRootOriginBarycenter<dataType>(baryMergeTree);
        postprocessingPipeline<dataType>(&(baryMergeTree.tree));
        for(unsigned int i = 0; i < trees.size(); ++i) {
          convertBranchDecompositionMatching<dataType>(
            &(baryMergeTree.tree), treesT[i], finalMatchings[i]);
        }
      }
    }

    template <class dataType>
    void execute(
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &finalMatchings,
      ftm::MergeTree<dataType> &baryMergeTree,
      bool finalAsgnDoubleInput = false,
      bool finalAsgnFirstInput = true) {
      std::vector<double> alphas;
      if(trees.size() != 2) {
        for(unsigned int i = 0; i < trees.size(); ++i)
          alphas.push_back(1.0 / trees.size());
      } else {
        alphas.push_back(alpha_);
        alphas.push_back(1 - alpha_);
      }

      execute<dataType>(trees, alphas, finalMatchings, baryMergeTree,
                        finalAsgnDoubleInput, finalAsgnFirstInput);
    }

    // ------------------------------------------------------------------------
    // Preprocessing
    // ------------------------------------------------------------------------
    template <class dataType>
    void limitSizePercent(ftm::MergeTree<dataType> &bary,
                          std::vector<ftm::FTMTree_MT *> &trees,
                          double percent,
                          bool useBD) {
      auto metric = getSizeLimitMetric(trees);
      unsigned int const newNoNodes = metric * percent / 100.0;
      keepMostImportantPairs<dataType>(&(bary.tree), newNoNodes, useBD);

      unsigned int const noNodesAfter = bary.tree.getRealNumberOfNodes();
      if(bary.tree.isFullMerge() and noNodesAfter > newNoNodes * 1.1 + 1
         and noNodesAfter > 3) {
        std::cout << "metric = " << metric << std::endl;
        std::cout << "newNoNodes = " << newNoNodes << std::endl;
        std::cout << "noNodesAfter = " << noNodesAfter << std::endl;
      }
    }

    template <class dataType>
    void limitSizeBarycenter(ftm::MergeTree<dataType> &bary,
                             std::vector<ftm::FTMTree_MT *> &trees,
                             unsigned int barycenterMaximumNumberOfPairs,
                             double percent,
                             bool useBD = true) {
      if(barycenterMaximumNumberOfPairs > 0)
        keepMostImportantPairs<dataType>(
          &(bary.tree), barycenterMaximumNumberOfPairs, useBD);
      if(percent > 0)
        limitSizePercent(bary, trees, percent, useBD);
    }
    template <class dataType>
    void limitSizeBarycenter(ftm::MergeTree<dataType> &bary,
                             std::vector<ftm::FTMTree_MT *> &trees,
                             double percent,
                             bool useBD = true) {
      limitSizeBarycenter(
        bary, trees, barycenterMaximumNumberOfPairs_, percent, useBD);
    }
    template <class dataType>
    void limitSizeBarycenter(ftm::MergeTree<dataType> &bary,
                             std::vector<ftm::FTMTree_MT *> &trees,
                             bool useBD = true) {
      limitSizeBarycenter(bary, trees, barycenterMaximumNumberOfPairs_,
                          barycenterSizeLimitPercent_, useBD);
    }

    // ------------------------------------------------------------------------
    // Postprocessing
    // ------------------------------------------------------------------------
    template <class dataType>
    void fixMergedRootOriginBarycenter(ftm::MergeTree<dataType> &barycenter) {
      if(not barycenter.tree.isFullMerge())
        return;

      ftm::FTMTree_MT *tree = &(barycenter.tree);
      auto tup = fixMergedRootOrigin<dataType>(tree);
      int maxIndex = std::get<0>(tup);
      dataType oldOriginValue = std::get<1>(tup);

      // Verify that scalars are consistent
      ftm::idNode const treeRoot = tree->getRoot();
      std::vector<dataType> newScalarsVector;
      ftm::getTreeScalars<dataType>(tree, newScalarsVector);
      bool isJT = tree->isJoinTree<dataType>();
      if((isJT and tree->getValue<dataType>(maxIndex) > oldOriginValue)
         or (not isJT
             and tree->getValue<dataType>(maxIndex) < oldOriginValue)) {
        newScalarsVector[treeRoot] = newScalarsVector[maxIndex];
        newScalarsVector[maxIndex] = oldOriginValue;
      } else
        newScalarsVector[treeRoot] = oldOriginValue;
      setTreeScalars(barycenter, newScalarsVector);
    }

    // ------------------------------------------------------------------------
    // Utils
    // ------------------------------------------------------------------------
    void printBaryStats(ftm::FTMTree_MT *baryTree,
                        const debug::Priority &priority
                        = debug::Priority::INFO) {
      auto noNodesT = baryTree->getNumberOfNodes();
      auto noNodes = baryTree->getRealNumberOfNodes();
      std::stringstream ss;
      ss << "Barycenter number of nodes : " << noNodes << " / " << noNodesT;
      printMsg(ss.str(), priority);
    }

    // ------------------------------------------------------------------------
    // Testing
    // ------------------------------------------------------------------------
    template <class dataType>
    void verifyBarycenterTwoTrees(
      std::vector<ftm::FTMTree_MT *> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &finalMatchings,
      std::vector<dataType> distances) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching;
      dataType distance;
      computeOneDistance(trees[0], trees[1], matching, distance);
      if(distance != (distances[0] + distances[1])) {
        std::stringstream ss, ss2, ss3, ss4;
        ss << "distance T1 T2    : " << distance;
        printMsg(ss.str());
        ss2 << "distance T1 T' T2 : " << distances[0] + distances[1];
        printMsg(ss2.str());
        ss3 << "distance T1 T'    : " << distances[0];
        printMsg(ss3.str());
        ss4 << "distance T' T2    : " << distances[1];
        printMsg(ss4.str());
      }
      return;

      auto baryTree = &(baryMergeTree.tree);
      std::vector<std::vector<ftm::idNode>> baryMatched(
        baryTree->getNumberOfNodes(),
        std::vector<ftm::idNode>(
          trees.size(), std::numeric_limits<ftm::idNode>::max()));
      for(unsigned int i = 0; i < finalMatchings.size(); ++i)
        for(auto match : finalMatchings[i])
          baryMatched[std::get<0>(match)][i] = std::get<1>(match);

      std::queue<ftm::idNode> queue;
      queue.emplace(baryTree->getRoot());
      while(!queue.empty()) {
        auto node = queue.front();
        queue.pop();
        std::vector<dataType> costs(trees.size());
        for(unsigned int i = 0; i < trees.size(); ++i)
          if(baryMatched[node][i] != std::numeric_limits<ftm::idNode>::max())
            costs[i] = relabelCost<dataType>(
              baryTree, node, trees[i], baryMatched[node][i]);
          else
            costs[i] = deleteCost<dataType>(baryTree, node);
        dataType cost = 0;
        if(baryMatched[node][0] != std::numeric_limits<ftm::idNode>::max()
           and baryMatched[node][1] != std::numeric_limits<ftm::idNode>::max())
          cost = relabelCost<dataType>(
            trees[0], baryMatched[node][0], trees[1], baryMatched[node][1]);
        else if(baryMatched[node][0] == std::numeric_limits<ftm::idNode>::max())
          cost = deleteCost<dataType>(trees[1], baryMatched[node][1]);
        else if(baryMatched[node][1] == std::numeric_limits<ftm::idNode>::max())
          cost = deleteCost<dataType>(trees[0], baryMatched[node][0]);
        else
          printErr("problem");
        costs[0] = std::sqrt(costs[0]);
        costs[1] = std::sqrt(costs[1]);
        cost = std::sqrt(cost);
        if(std::abs((double)(costs[0] - costs[1])) > 1e-7) {
          printMsg(debug::Separator::L1);
          std::stringstream ss, ss2, ss3, ss4;
          ss << "cost T' T0    : " << costs[0];
          printMsg(ss.str());
          ss2 << "cost T' T1    : " << costs[1];
          printMsg(ss2.str());
          ss3 << "cost T0 T1    : " << cost;
          printMsg(ss2.str());
          ss4 << "cost T0 T' T1 : " << costs[0] + costs[1];
          printMsg(ss4.str());
          if(std::abs((double)((costs[0] + costs[1]) - cost)) > 1e-7) {
            std::stringstream ss5;
            ss5 << "diff          : "
                << std::abs((double)((costs[0] + costs[1]) - cost));
            printMsg(ss5.str());
          }
          std::stringstream ss6;
          ss6 << "diff2         : " << std::abs((double)(costs[0] - costs[1]));
          printMsg(ss.str());
          // baryTree->printNode2<dataType>(node);
          // baryTree->printNode2<dataType>(baryTree->getParentSafe(node));
          for(unsigned int i = 0; i < 2; ++i)
            if(baryMatched[node][i]
               != std::numeric_limits<ftm::idNode>::max()) {
              printMsg(
                trees[i]->printNode2<dataType>(baryMatched[node][i]).str());
              printMsg(trees[i]
                         ->printNode2<dataType>(
                           trees[i]->getParentSafe(baryMatched[node][i]))
                         .str());
            }
        }
        std::vector<ftm::idNode> children;
        baryTree->getChildren(node, children);
        for(auto child : children)
          queue.emplace(child);
      }
    }

  }; // MergeTreeBarycenter class

} // namespace ttk
