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

#ifndef _MERGETREEBARYCENTER_H
#define _MERGETREEBARYCENTER_H

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
#ifdef TTK_ENABLE_OPENMP
      omp_set_nested(1);
#endif
    };
    ~MergeTreeBarycenter(){};

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
    // ----------------------------------------
    // Initialization
    // ----------------------------------------
    template <class dataType>
    void getDistanceMatrix(std::vector<ftm::FTMTree_MT *> &trees,
                           std::vector<std::vector<double>> &distanceMatrix) {
      distanceMatrix = std::vector<std::vector<double>>(
        trees.size(), std::vector<double>(trees.size(), 0));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
      for(unsigned int i = 0; i < trees.size(); ++i)
        for(unsigned int j = i + 1; j < trees.size(); ++j) {
          std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching;
          dataType distance;
          computeOneDistance<dataType>(trees[i], trees[j], matching, distance);
          distanceMatrix[i][j] = distance;
          distanceMatrix[j][i] = distance;
        }
    }

    template <class dataType>
    int getBestInitTreeIndex(std::vector<ftm::FTMTree_MT *> &trees,
                             std::vector<ftm::FTMTree_MT *> &trees2,
                             bool distMinimizer = true) {
      std::vector<std::vector<double>> distanceMatrix, distanceMatrix2;
      if(distMinimizer) {
        getDistanceMatrix<dataType>(trees, distanceMatrix);
        if(trees2.size() != 0)
          getDistanceMatrix<dataType>(trees2, distanceMatrix2);
      }

      int bestIndex = -1;
      dataType bestValue
        = distMinimizer ? std::numeric_limits<dataType>::max() : 0;
      std::vector<int> sizes(trees.size());
      for(unsigned int i = 0; i < trees.size(); ++i) {
        dataType value = 0;
        if(distMinimizer) {
          for(auto v : distanceMatrix[i])
            value += v;
          if(trees2.size() != 0)
            for(auto v : distanceMatrix2[i])
              value += v;
        } else {
          value = trees[i]->getRealNumberOfNodes();
          if(trees2.size() != 0)
            value += trees2[i]->getRealNumberOfNodes();
        }
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
                             bool distMinimizer = true) {
      std::vector<ftm::FTMTree_MT *> trees2;
      return getBestInitTreeIndex<dataType>(trees, trees2, distMinimizer);
    }

    template <class dataType>
    void initBarycenterTree(std::vector<ftm::FTMTree_MT *> &trees,
                            MergeTree<dataType> &baryTree,
                            bool distMinimizer = true) {
      int bestIndex = getBestInitTreeIndex<dataType>(trees, distMinimizer);
      baryTree = copyMergeTree<dataType>(trees[bestIndex], true);
    }

    // ----------------------------------------
    // Update
    // ----------------------------------------
    template <class dataType>
    ftm::idNode getNodesAndScalarsToAdd(
      MergeTree<dataType> &mTree1,
      ftm::idNode nodeId1,
      ftm::FTMTree_MT *tree2,
      ftm::idNode nodeId2,
      std::vector<dataType> &newScalarsVector,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> &nodesToProcess,
      ftm::idNode nodeCpt,
      int i) {
      // Get nodes and scalars to add
      std::queue<std::tuple<ftm::idNode, ftm::idNode>> queue;
      queue.emplace(std::make_tuple(nodeId2, nodeId1));
      nodesToProcess.push_back(std::make_tuple(nodeId2, nodeId1, i));
      while(!queue.empty()) {
        auto queueTuple = queue.front();
        queue.pop();
        ftm::idNode node = std::get<0>(queueTuple);
        // Get scalars
        newScalarsVector.push_back(
          tree2->getValue<dataType>(tree2->getNode(node)->getOrigin()));
        newScalarsVector.push_back(tree2->getValue<dataType>(node));
        // Process children
        std::vector<ftm::idNode> children;
        tree2->getChildren(node, children);
        for(auto child : children) {
          queue.emplace(std::make_tuple(child, nodeCpt + 1));
          nodesToProcess.push_back(std::make_tuple(child, nodeCpt + 1, i));
        }
        nodeCpt += 2; // we will add two nodes (birth and death)
      }

      return nodeCpt;
    }

    template <class dataType>
    void addNodes(
      MergeTree<dataType> &mTree1,
      int noTrees,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> &nodesToProcess,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>
        &nodesProcessed) {
      ftm::FTMTree_MT *tree1 = &(mTree1.tree);

      // Add nodes
      nodesProcessed
        = std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>(
          noTrees);
      for(auto processTuple : nodesToProcess) {
        ftm::idNode parent = std::get<1>(processTuple);
        ftm::idNode nodeTree1 = tree1->getNumberOfNodes();
        int index = std::get<2>(processTuple);
        nodesProcessed[index].push_back(
          std::make_tuple(nodeTree1 + 1, std::get<0>(processTuple)));
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
      MergeTree<dataType> &mTree1,
      int noTrees,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> &nodesToProcess,
      std::vector<dataType> &newScalarsVector,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>
        &nodesProcessed) {
      ftm::FTMTree_MT *tree1 = &(mTree1.tree);

      // Create new tree
      MergeTree<dataType> mTreeNew
        = createEmptyMergeTree<dataType>(newScalarsVector.size());
      setTreeScalars<dataType>(mTreeNew, newScalarsVector);
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
      MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings) {
      ftm::FTMTree_MT *baryTree = &(baryMergeTree.tree);
      ftm::idNode baryTreeRoot = baryTree->getRoot();

      // Init matching matrix
      // m[i][j] contains the node in the barycenter matched to the jth node of
      // the ith tree
      std::vector<std::vector<ftm::idNode>> matrixMatchings(trees.size());
      std::vector<bool> baryMatched(baryTree->getNumberOfNodes(), false);
      for(unsigned int i = 0; i < matchings.size(); ++i) {
        auto matching = matchings[i];
        std::vector<ftm::idNode> matchingT(trees[i]->getNumberOfNodes(), -1);
        for(auto match : matching) {
          matchingT[std::get<1>(match)] = std::get<0>(match);
          baryMatched[std::get<0>(match)] = true;
        }
        matrixMatchings[i].insert(
          matrixMatchings[i].end(), matchingT.begin(), matchingT.end());
      }

      // Iterate through trees to get the nodes to add in the barycenter
      std::vector<std::vector<ftm::idNode>> nodesToAdd(trees.size());
      for(unsigned int i = 0; i < trees.size(); ++i) {
        ftm::idNode root = trees[i]->getRoot();
        std::queue<ftm::idNode> queue;
        queue.emplace(root);
        while(!queue.empty()) {
          ftm::idNode node = queue.front();
          queue.pop();
          bool processChildren = true;
          if((int)matrixMatchings[i][node]
             == -1) { // if node in trees[i] is not matched
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

      // Delete nodes that are not matched in the barycenter
      for(unsigned int i = 0; i < baryTree->getNumberOfNodes(); ++i)
        if(not baryMatched[i])
          baryTree->deleteNode(i);

      if(not keepSubtree_) {
        // Add scalars and nodes not present in the barycenter
        ftm::idNode nodeCpt = baryTree->getNumberOfNodes();
        std::vector<std::tuple<ftm::idNode, ftm::idNode, int>> nodesToProcess;
        std::vector<dataType> newScalarsVector;
        getTreeScalars<dataType>(baryMergeTree, newScalarsVector);
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
              trees[i]->printTree();
              trees[i]->printPairsFromTree<dataType>(true);
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
              nodesProcessedT.push_back(
                std::make_tuple(std::get<0>(tup), std::get<1>(tup), -1));
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
    std::tuple<double, double>
      getParametrizedBirthDeath(ftm::FTMTree_MT *tree1,
                                ftm::idNode nodeId1,
                                ftm::FTMTree_MT *tree2 = nullptr,
                                ftm::idNode nodeId2 = ftm::nullNodes) {
      std::tuple<double, double> birthDeath;
      // Normalized Wasserstein
      if(normalizedWasserstein_ and not rescaledWasserstein_)
        birthDeath = getNormalizedBirthDeathDouble<dataType>(tree1, nodeId1);
      // Rescaled Wasserstein
      else if(normalizedWasserstein_ and rescaledWasserstein_)
        birthDeath
          = getRescaledBirthDeath<dataType>(tree1, nodeId1, tree2, nodeId2);
      // Classical Wasserstein
      else
        birthDeath = tree1->getBirthDeath<dataType>(nodeId1);
      return birthDeath;
    }

    template <class dataType>
    std::tuple<dataType, dataType> getParametrizedBirthDeathFromVector(
      ftm::FTMTree_MT *tree1,
      ftm::idNode nodeId1,
      ftm::FTMTree_MT *tree2,
      ftm::idNode nodeId2,
      std::vector<dataType> &newScalarsVector) {
      if(normalizedWasserstein_ and rescaledWasserstein_)
        return getRescaledBirthDeathFromVector<dataType>(
          tree1, nodeId1, tree2, nodeId2, newScalarsVector);
      return getParametrizedBirthDeath<dataType>(
        tree1, nodeId1, tree2, nodeId2);
    }

    template <class dataType>
    std::tuple<dataType, dataType>
      interpolation(MergeTree<dataType> &baryMergeTree,
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

      // Compute projection h
      double tempBirth = 0, tempDeath = 0;
      int offDiagonal = 0;
      double alphaSum = 0;
      for(unsigned int i = 0; i < trees.size(); ++i)
        if((int)nodes[i] != -1)
          alphaSum += alphas[i];
      for(unsigned int i = 0; i < trees.size(); ++i) {
        if((int)nodes[i] != -1) { // if node is matched in trees[i]
          auto iBirthDeath = getParametrizedBirthDeathFromVector<dataType>(
            trees[i], nodes[i], baryTree, nodeId, newScalarsVector);
          double tTempBirth = 0, tTempDeath = 0;
          tTempBirth += std::get<0>(iBirthDeath);
          tTempDeath += std::get<1>(iBirthDeath);
          if(normalizedWasserstein_ and rescaledWasserstein_) {
            auto newMinMax = getNewMinMaxFromVector<dataType>(
              trees[i], nodes[i], baryTree, nodeId, newScalarsVector);
            tTempBirth /= (std::get<1>(newMinMax) - std::get<0>(newMinMax));
            tTempDeath /= (std::get<1>(newMinMax) - std::get<0>(newMinMax));
          }
          tempBirth += tTempBirth * alphas[i] / alphaSum;
          tempDeath += tTempDeath * alphas[i] / alphaSum;
          ++offDiagonal;
        }
      }
      double projec = (tempBirth + tempDeath) / 2;

      // Compute newBirth and newDeath
      dataType divisor = 0;
      for(unsigned int i = 0; i < trees.size(); ++i) {
        double iBirth = projec, iDeath = projec;
        if((int)nodes[i] != -1) { // if node is matched in trees[i]
          auto iBirthDeath = getParametrizedBirthDeathFromVector<dataType>(
            trees[i], nodes[i], baryTree, nodeId, newScalarsVector);
          iBirth = std::get<0>(iBirthDeath);
          iDeath = std::get<1>(iBirthDeath);
        }
        if(normalizedWasserstein_ and rescaledWasserstein_) {
          dataType beta_max = 1, beta_min = 0;
          auto newMinMax
            = ((int)nodes[i] == -1)
                ? getNewMinMax<dataType>(baryTree, nodeId, baryTree, nodeId)
                : getNewMinMaxFromVector<dataType>(
                  trees[i], nodes[i], baryTree, nodeId, newScalarsVector);
          if((int)nodes[i] == -1) {
            beta_max = mu_max;
            beta_min = mu_min;
            iBirth *= (beta_max - beta_min);
            iDeath *= (beta_max - beta_min);
          } else {
            beta_min = std::get<0>(newMinMax);
            beta_max = std::get<1>(newMinMax);
          }
          iBirth *= (beta_max - beta_min);
          iDeath *= (beta_max - beta_min);
          divisor += alphas[i] * (beta_max - beta_min) * (beta_max - beta_min);
        }
        newBirth += alphas[i] * iBirth;
        newDeath += alphas[i] * iDeath;
      }
      if(normalizedWasserstein_ and rescaledWasserstein_) {
        newBirth /= divisor;
        newDeath /= divisor;
      }
      if(normalizedWasserstein_ or rescaledWasserstein_) {
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
                         MergeTree<dataType> &baryMergeTree,
                         ftm::idNode nodeB,
                         std::vector<dataType> &newScalarsVector) {
      ftm::FTMTree_MT *baryTree = &(baryMergeTree.tree);
      dataType mu_max = getMinMaxLocalFromVector<dataType>(
        baryTree, nodeB, newScalarsVector, false);
      dataType mu_min
        = getMinMaxLocalFromVector<dataType>(baryTree, nodeB, newScalarsVector);

      auto birthDeath = getParametrizedBirthDeathFromVector<dataType>(
        tree, nodeId, baryTree, nodeB, newScalarsVector);
      double newBirth = std::get<0>(birthDeath);
      double newDeath = std::get<1>(birthDeath);
      double projec = (newBirth + newDeath) / 2;

      dataType beta_min = 0, beta_max = 0, divisor = 1;
      if(normalizedWasserstein_ and rescaledWasserstein_) {
        auto newMinMax = getNewMinMaxFromVector<dataType>(
          tree, nodeId, baryTree, nodeB, newScalarsVector);
        beta_min = std::get<0>(newMinMax);
        beta_max = std::get<1>(newMinMax);
        newBirth *= (beta_max - beta_min);
        newDeath *= (beta_max - beta_min);
        projec = projec * (mu_max - mu_min) * (mu_max - mu_min)
                 / (beta_max - beta_min);
        divisor = alpha * (beta_max - beta_min) * (beta_max - beta_min)
                  + (1 - alpha) * (mu_max - mu_min) * (mu_max - mu_min);
      }

      newBirth = alpha * newBirth + (1 - alpha) * projec;
      newDeath = alpha * newDeath + (1 - alpha) * projec;

      if(normalizedWasserstein_ and rescaledWasserstein_) {
        newBirth /= divisor;
        newDeath /= divisor;
      }

      if(normalizedWasserstein_ or rescaledWasserstein_) {
        newBirth = newBirth * (mu_max - mu_min) + mu_min;
        newDeath = newDeath * (mu_max - mu_min) + mu_min;
      }

      dataType newBirthT = newBirth;
      dataType newDeathT = newDeath;
      return std::make_tuple(newBirthT, newDeathT);
    }

    template <class dataType>
    void purgeBarycenter(MergeTree<dataType> &baryMergeTree,
                         std::vector<std::vector<ftm::idNode>> &baryMatching,
                         std::vector<ftm::FTMTree_MT *> &trees,
                         std::vector<double> &alphas) {
      ftm::FTMTree_MT *baryTree = &(baryMergeTree.tree);
      std::vector<bool> nodesProcessed(baryTree->getNumberOfNodes(), false);
      std::vector<dataType> nodesMatchingCost(baryTree->getNumberOfNodes(), 0);
      std::vector<dataType> nodesDestructCost(baryTree->getNumberOfNodes(), 0);
      std::vector<ftm::idNode> leaves;
      baryTree->getLeavesFromTree(leaves);
      std::queue<ftm::idNode> queue;
      for(auto leaf : leaves)
        queue.emplace(leaf);
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();
        for(unsigned int i = 0; i < trees.size(); ++i) {
          dataType newMatchingCost = alphas[i];
          dataType newDestructCost = alphas[i];
          if((int)baryMatching[node][i] != -1) {
            newMatchingCost *= relabelCost<dataType>(
              baryTree, node, trees[i], baryMatching[node][i]);
            newDestructCost
              *= deleteCost<dataType>(trees[i], baryMatching[node][i]);
          } else {
            newMatchingCost *= deleteCost<dataType>(baryTree, node);
            newDestructCost *= 0;
          }
          nodesMatchingCost[node] += newMatchingCost;
          nodesDestructCost[node] += newDestructCost;
        }
        std::vector<ftm::idNode> children;
        baryTree->getChildren(node, children);
        for(auto child : children) {
          nodesMatchingCost[node] += nodesMatchingCost[child];
          nodesDestructCost[node] += nodesDestructCost[child];
        }
        nodesProcessed[node] = true;
        if(not nodesProcessed[baryTree->getParentSafe(node)])
          queue.emplace(baryTree->getParentSafe(node));

        // Destruct subtree if better
        if(nodesDestructCost[node] < nodesMatchingCost[node]) {
          baryTree->deleteSubtree(node);
          nodesDestructCost[node] = 0;
          nodesMatchingCost[node] = 0;
        }
      }
    }

    template <class dataType>
    void updateBarycenterTreeScalars(
      std::vector<ftm::FTMTree_MT *> &trees,
      MergeTree<dataType> &baryMergeTree,
      std::vector<double> &alphas,
      int indexAddedNodes,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings) {
      ftm::FTMTree_MT *baryTree = &(baryMergeTree.tree);
      bool isJT = baryTree->isJoinTree<dataType>();

      // Init matching matrix
      // m[i][j] contains the node in trees[j] matched to the node i in the
      // barycenter
      std::vector<std::vector<ftm::idNode>> baryMatching(
        baryTree->getNumberOfNodes(),
        std::vector<ftm::idNode>(trees.size(), -1));
      std::vector<int> nodesAddedTree(baryTree->getNumberOfNodes(), -1);
      for(unsigned int i = 0; i < matchings.size(); ++i) {
        auto matching = matchings[i];
        for(auto match : matching) {
          baryMatching[std::get<0>(match)][i] = std::get<1>(match);
          if((int)std::get<0>(match)
             >= indexAddedNodes) // get the tree of this added node
            nodesAddedTree[std::get<0>(match)] = i;
        }
      }

      // Interpolate scalars
      std::vector<dataType> newScalarsVector(baryTree->getNumberOfNodes());
      ftm::idNode root = baryTree->getRoot();
      std::queue<ftm::idNode> queue;
      queue.emplace(root);
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();
        std::tuple<dataType, dataType> newBirthDeath;
        if((int)node < indexAddedNodes) {
          newBirthDeath
            = interpolation<dataType>(baryMergeTree, node, newScalarsVector,
                                      trees, baryMatching[node], alphas);
        } else {
          int i = nodesAddedTree[node];
          ftm::idNode nodeT = baryMatching[node][i];
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

      setTreeScalars(baryMergeTree, newScalarsVector);
      if(normalizedWasserstein_ and rescaledWasserstein_)
        purgeBarycenter<dataType>(baryMergeTree, baryMatching, trees, alphas);
      std::vector<ftm::idNode> deletedNodesT;
      persistenceThresholding<dataType>(
        &(baryMergeTree.tree), 0, deletedNodesT);
      if(not isCalled_)
        cleanMergeTree<dataType>(baryMergeTree);
    }

    int getNumberOfRoots(ftm::FTMTree_MT *tree) {
      int noRoots = 0;
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
        noRoots += (tree->isRoot(i) and not tree->isLeaf(i)) ? 1 : 0;
      return noRoots;
    }

    template <class dataType>
    void updateBarycenterTree(
      std::vector<ftm::FTMTree_MT *> &trees,
      MergeTree<dataType> &baryMergeTree,
      std::vector<double> &alphas,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings) {
      int indexAddedNodes = baryMergeTree.tree.getNumberOfNodes();
      updateBarycenterTreeStructure<dataType>(trees, baryMergeTree, matchings);
      updateBarycenterTreeScalars<dataType>(
        trees, baryMergeTree, alphas, indexAddedNodes, matchings);
    }

    template <class dataType>
    void computeOneDistance(
      ftm::FTMTree_MT *tree,
      ftm::FTMTree_MT *baryTree,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      dataType &distance) {
      // Timer t_distance;
      MergeTreeDistance mergeTreeDistance;
      mergeTreeDistance.setDebugLevel(2);
      mergeTreeDistance.setProgressiveComputation(false);
      mergeTreeDistance.setPreprocess(false);
      mergeTreeDistance.setPostprocess(false);
      mergeTreeDistance.setBranchDecomposition(true);
      mergeTreeDistance.setNormalizedWasserstein(normalizedWasserstein_);
      mergeTreeDistance.setNormalizedWassersteinReg(normalizedWassersteinReg_);
      mergeTreeDistance.setRescaledWasserstein(rescaledWasserstein_);
      mergeTreeDistance.setKeepSubtree(keepSubtree_);
      mergeTreeDistance.setAssignmentSolver(assignmentSolverID_);
      mergeTreeDistance.setIsCalled(true);
      mergeTreeDistance.setThreadNumber(this->threadNumber_);
      mergeTreeDistance.setDistanceSquared(true); // squared root
      mergeTreeDistance.setNodePerTask(nodePerTask_);
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

    // ----------------------------------------
    // Assignment
    // ----------------------------------------
    template <class dataType>
    void computeOneDistance(
      ftm::FTMTree_MT *tree,
      MergeTree<dataType> &baryMergeTree,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      dataType &distance) {
      computeOneDistance<dataType>(
        tree, &(baryMergeTree.tree), matching, distance);
    }

    template <class dataType>
    void computeOneDistance(
      MergeTree<dataType> &baryMergeTree,
      MergeTree<dataType> &baryMergeTree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      dataType &distance) {
      computeOneDistance<dataType>(
        &(baryMergeTree.tree), baryMergeTree2, matching, distance);
    }

    template <class dataType>
    void assignment(
      std::vector<ftm::FTMTree_MT *> &trees,
      MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<dataType> &distances) {
      if(not isCalled_)
        assignmentPara(trees, baryMergeTree, matchings, distances);
      else
        assignmentTask(trees, baryMergeTree, matchings, distances);
    }

    template <class dataType>
    void assignmentPara(
      std::vector<ftm::FTMTree_MT *> &trees,
      MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<dataType> &distances) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) \
  shared(baryMergeTree) if(parallelize_)
      {
#pragma omp single nowait
#endif
        assignmentTask(trees, baryMergeTree, matchings, distances);
#ifdef TTK_ENABLE_OPENMP
      } // pragma omp parallel
#endif
    }

    template <class dataType>
    void assignmentTask(
      std::vector<ftm::FTMTree_MT *> &trees,
      MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<dataType> &distances) {
      for(unsigned int i = 0; i < trees.size(); ++i)
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(i) \
  untied shared(baryMergeTree, matchings, distances)
#endif
        computeOneDistance<dataType>(
          trees[i], baryMergeTree, matchings[i], distances[i]);
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
    }

    // ----------------------------------------
    //
    // ----------------------------------------
    template <class dataType>
    unsigned int
      persistenceScaling(std::vector<ftm::FTMTree_MT *> &trees,
                         std::vector<MergeTree<dataType>> &mergeTrees,
                         std::vector<ftm::FTMTree_MT *> &oriTrees,
                         int iterationNumber,
                         std::vector<std::vector<ftm::idNode>> &deletedNodes) {
      deletedNodes = std::vector<std::vector<ftm::idNode>>(oriTrees.size());
      unsigned int noTreesUnscaled = 0;

      // Scale trees
      for(unsigned int i = 0; i < oriTrees.size(); ++i) {
        double persistenceThreshold = 50.0;
        if(iterationNumber != -1) {
          // Get number of pairs in scaled merge tree
          int noPairs = mergeTrees[i].tree.getRealNumberOfNodes();

          // Get pairs in original merge tree
          std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> pairs;
          oriTrees[i]->getPersistencePairsFromTree<dataType>(
            pairs, branchDecomposition_);

          // Compute new persistence threshold
          double multiplier = (progressiveSpeedDivisor_ < 1e-6
                                 ? 1.
                                 : iterationNumber / progressiveSpeedDivisor_);
          int decrement = multiplier * pairs.size() / 10;
          int thresholdIndex = pairs.size() - noPairs - std::max(decrement, 2);
          thresholdIndex = std::max(thresholdIndex, 0);
          dataType persistence = std::get<2>(pairs[thresholdIndex]);
          persistenceThreshold
            = persistence / std::get<2>(pairs[pairs.size() - 1]) * 100;
          if(thresholdIndex == 0) {
            persistenceThreshold = 0.;
            ++noTreesUnscaled;
          }
        }
        if(persistenceThreshold != 0.) {
          MergeTree<dataType> mt = copyMergeTree<dataType>(oriTrees[i]);
          persistenceThresholding<dataType>(
            &(mt.tree), persistenceThreshold, deletedNodes[i]);
          if(mergeTrees.size() == 0)
            mergeTrees = std::vector<MergeTree<dataType>>(oriTrees.size());
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

    void printBaryStats(ftm::FTMTree_MT *baryTree) {
      auto noNodesT = baryTree->getNumberOfNodes();
      auto noNodes = baryTree->getRealNumberOfNodes();
      std::stringstream ss;
      ss << "Barycenter number of nodes : " << noNodes << " / " << noNodesT;
      printMsg(ss.str());
    }

    // ----------------------------------------
    // Main Functions
    // ----------------------------------------
    template <class dataType>
    void computeBarycenter(
      std::vector<ftm::FTMTree_MT *> &trees,
      MergeTree<dataType> &baryMergeTree,
      std::vector<double> &alphas,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &finalMatchings) {
      Timer t_bary;

      ftm::FTMTree_MT *baryTree = &(baryMergeTree.tree);

      // Persistence scaling
      std::vector<ftm::FTMTree_MT *> oriTrees;
      std::vector<MergeTree<dataType>> scaledMergeTrees;
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
          unsigned int noTreesUnscaled = persistenceScaling<dataType>(
            trees, scaledMergeTrees, oriTrees, NoIteration, deletedNodes);
          treesUnscaled = (noTreesUnscaled == oriTrees.size());
        }
      }

      // Final processing
      printMsg(debug::Separator::L2);
      printMsg("Final assignment");

      std::vector<dataType> distances(trees.size(), -1);
      assignment<dataType>(trees, baryMergeTree, finalMatchings, distances);
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
      std::vector<MergeTree<dataType>> &trees,
      std::vector<double> &alphas,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &finalMatchings,
      MergeTree<dataType> &baryMergeTree) {
      // --- Preprocessing
      if(preprocess_) {
        treesNodeCorr_ = std::vector<std::vector<int>>(trees.size());
        for(unsigned int i = 0; i < trees.size(); ++i)
          preprocessingPipeline<dataType>(trees[i], epsilonTree2_,
                                          epsilon2Tree2_, epsilon3Tree2_,
                                          branchDecomposition_, useMinMaxPair_,
                                          cleanTree_, treesNodeCorr_[i]);
        printTreesStats(trees);
      }

      // --- Init barycenter
      std::vector<ftm::FTMTree_MT *> treesT;
      mergeTreeToFTMTree<dataType>(trees, treesT);
      initBarycenterTree<dataType>(treesT, baryMergeTree);

      // --- Execute
      computeBarycenter<dataType>(
        treesT, baryMergeTree, alphas, finalMatchings);

      // --- Postprocessing
      if(postprocess_) {
        std::vector<int> allRealNodes(trees.size());
        for(unsigned int i = 0; i < trees.size(); ++i) {
          postprocessingPipeline<dataType>(treesT[i]);
        }

        fixMergedRootOriginBarycenter<dataType>(baryMergeTree);
        postprocessingPipeline<dataType>(&(baryMergeTree.tree));
        for(unsigned int i = 0; i < trees.size(); ++i) {
          convertBranchDecompositionMatching<dataType>(
            &(baryMergeTree.tree), treesT[i], finalMatchings[i]);
        }
      }
    }

    template <class dataType>
    void execute(
      std::vector<MergeTree<dataType>> &trees,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &finalMatchings,
      MergeTree<dataType> &baryMergeTree) {
      std::vector<double> alphas;
      if(trees.size() != 2) {
        for(unsigned int i = 0; i < trees.size(); ++i)
          alphas.push_back(1.0 / trees.size());
      } else {
        alphas.push_back(alpha_);
        alphas.push_back(1 - alpha_);
      }

      execute<dataType>(trees, alphas, finalMatchings, baryMergeTree);
    }

    // ----------------------------------------
    // Postprocessing
    // ----------------------------------------
    template <class dataType>
    void fixMergedRootOriginBarycenter(MergeTree<dataType> &barycenter) {
      if(not barycenter.tree.isFullMerge())
        return;

      ftm::FTMTree_MT *tree = &(barycenter.tree);
      auto tup = fixMergedRootOrigin<dataType>(tree);
      int maxIndex = std::get<0>(tup);
      dataType oldOriginValue = std::get<1>(tup);

      // Verify that scalars are consistent
      ftm::idNode treeRoot = tree->getRoot();
      std::vector<dataType> newScalarsVector;
      getTreeScalars<dataType>(tree, newScalarsVector);
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

    // ----------------------------------------
    // Testing
    // ----------------------------------------
    template <class dataType>
    void verifyBarycenterTwoTrees(
      std::vector<ftm::FTMTree_MT *> &trees,
      MergeTree<dataType> &baryMergeTree,
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
        std::vector<ftm::idNode>(trees.size(), -1));
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
          if((int)baryMatched[node][i] != -1)
            costs[i] = relabelCost<dataType>(
              baryTree, node, trees[i], baryMatched[node][i]);
          else
            costs[i] = deleteCost<dataType>(baryTree, node);
        dataType cost = 0;
        if((int)baryMatched[node][0] != -1 and (int) baryMatched[node][1] != -1)
          cost = relabelCost<dataType>(
            trees[0], baryMatched[node][0], trees[1], baryMatched[node][1]);
        else if((int)baryMatched[node][0] == -1)
          cost = deleteCost<dataType>(trees[1], baryMatched[node][1]);
        else if((int)baryMatched[node][1] == -1)
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
            if((int)baryMatched[node][i] != -1) {
              trees[i]->printNode2<dataType>(baryMatched[node][i]);
              trees[i]->printNode2<dataType>(
                trees[i]->getParentSafe(baryMatched[node][i]));
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

#endif
