/// \ingroup base
/// \class MergeTreeBase
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// \b Related \b publication \n
/// "Wasserstein Distances, Geodesics and Barycenters of Merge Trees" \n
/// Mathieu Pont, Jules Vidal, Julie Delon, Julien Tierny.\n
/// Proc. of IEEE VIS 2021.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021

#pragma once

#include <AssignmentSolver.h>
#include <FTMNode.h>
#include <FTMTree.h>
#include <FTMTreePPUtils.h>
#include <FTMTreeUtils.h>

#include "MergeTreeUtils.h"

namespace ttk {

  class MergeTreeBase : virtual public Debug {
  protected:
    int assignmentSolverID_ = 0;
    bool epsilon1UseFarthestSaddle_ = false;
    double epsilonTree1_ = 0;
    double epsilonTree2_ = 0;
    double epsilon2Tree1_ = 100;
    double epsilon2Tree2_ = 100;
    double epsilon3Tree1_ = 100;
    double epsilon3Tree2_ = 100;
    double persistenceThreshold_ = 0;
    bool barycenterMergeTree_ = false;
    bool useMinMaxPair_ = true;
    bool deleteMultiPersPairs_ = false;

    bool branchDecomposition_ = true;
    int wassersteinPower_ = 2;
    bool normalizedWasserstein_ = true;
    bool keepSubtree_ = false;
    double nonMatchingWeight_ = 1.0;

    bool distanceSquaredRoot_ = true; // squared root
    bool useFullMerge_ = false;

    bool isPersistenceDiagram_ = false;
    bool convertToDiagram_ = false;

    // Double input
    double mixtureCoefficient_ = 0.5;
    bool useDoubleInput_ = false;

    // Old
    bool parallelize_ = true;
    int nodePerTask_ = 32;
    bool cleanTree_ = true;

    // Clean correspondence
    std::vector<std::vector<int>> treesNodeCorr_;

  public:
    MergeTreeBase() {
      this->setDebugMsgPrefix(
        "MergeTreeBase"); // inherited from Debug: prefix will be printed
                          // at the beginning of every msg
    }

    void setAssignmentSolver(int assignmentSolver) {
      assignmentSolverID_ = assignmentSolver;
    }

    void setEpsilon1UseFarthestSaddle(bool b) {
      epsilon1UseFarthestSaddle_ = b;
    }

    void setEpsilonTree1(double epsilon) {
      epsilonTree1_ = epsilon;
    }

    void setEpsilonTree2(double epsilon) {
      epsilonTree2_ = epsilon;
    }

    void setEpsilon2Tree1(double epsilon) {
      epsilon2Tree1_ = epsilon;
    }

    void setEpsilon2Tree2(double epsilon) {
      epsilon2Tree2_ = epsilon;
    }

    void setEpsilon3Tree1(double epsilon) {
      epsilon3Tree1_ = epsilon;
    }

    void setEpsilon3Tree2(double epsilon) {
      epsilon3Tree2_ = epsilon;
    }

    void setPersistenceThreshold(double pt) {
      persistenceThreshold_ = pt;
    }

    void setParallelize(bool para) {
      parallelize_ = para;
    }

    void setNodePerTask(int npt) {
      nodePerTask_ = npt;
    }

    void setBranchDecomposition(bool useBD) {
      branchDecomposition_ = useBD;
    }

    void setNormalizedWasserstein(bool normalizedWasserstein) {
      normalizedWasserstein_ = normalizedWasserstein;
    }

    void setKeepSubtree(bool keepSubtree) {
      keepSubtree_ = keepSubtree;
    }

    void setNonMatchingWeight(double weight) {
      nonMatchingWeight_ = weight;
    }

    void setBarycenterMergeTree(bool imt) {
      barycenterMergeTree_ = imt;
    }

    void setDistanceSquaredRoot(bool distanceSquaredRoot) {
      distanceSquaredRoot_ = distanceSquaredRoot;
    }

    void setUseMinMaxPair(bool useMinMaxPair) {
      useMinMaxPair_ = useMinMaxPair;
    }

    void setDeleteMultiPersPairs(bool deleteMultiPersPairsT) {
      deleteMultiPersPairs_ = deleteMultiPersPairsT;
    }

    void setCleanTree(bool clean) {
      cleanTree_ = clean;
    }

    void setIsPersistenceDiagram(bool isPD) {
      isPersistenceDiagram_ = isPD;
    }

    std::vector<std::vector<int>> getTreesNodeCorr() {
      return treesNodeCorr_;
    }

    // ------------------------------------------------------------------------
    // Double Input
    // ------------------------------------------------------------------------
    double mixDistancesMinMaxPairWeight(bool isFirstInput) {
      return (
        mixtureCoefficient_ == 0.0 or mixtureCoefficient_ == 1.0
          ? (isFirstInput ? mixtureCoefficient_ : (1.0 - mixtureCoefficient_))
          : (isFirstInput ? 1.0 / std::pow(mixDistancesWeight(isFirstInput), 2)
                          : 0.0));
    }

    double mixDistancesWeight(bool isFirstInput) {
      return (isFirstInput ? std::min(mixtureCoefficient_ * 2, 1.0)
                           : std::min(-mixtureCoefficient_ * 2 + 2, 1.0));
    }

    template <class dataType>
    double mixDistances(dataType distance1, dataType distance2) {
      return mixDistancesWeight(true) * distance1
             + mixDistancesWeight(false) * distance2;
    }

    template <class dataType>
    void
      mixDistancesMatrix(std::vector<std::vector<dataType>> &distanceMatrix,
                         std::vector<std::vector<dataType>> &distanceMatrix2) {
      for(unsigned int i = 0; i < distanceMatrix.size(); ++i)
        for(unsigned int j = 0; j < distanceMatrix[i].size(); ++j)
          distanceMatrix[i][j] = mixDistances<dataType>(
            distanceMatrix[i][j], distanceMatrix2[i][j]);
    }

    // ------------------------------------------------------------------------
    // Tree Preprocessing
    // ------------------------------------------------------------------------
    // Epsilon 1 processing
    template <class dataType>
    void mergeSaddle(ftm::FTMTree_MT *tree,
                     double epsilon,
                     std::vector<std::vector<ftm::idNode>> &treeNodeMerged,
                     bool mergeByPersistence = false) {
      bool fullMerge = (epsilon == 100);
      fullMerge &= useFullMerge_;

      treeNodeMerged.clear();
      treeNodeMerged.resize(tree->getNumberOfNodes());

      if(mergeByPersistence)
        ftm::computePersistencePairs<dataType>(
          tree); // need to have the pairing (if merge by persistence)

      // Compute epsilon value
      dataType maxValue = tree->getValue<dataType>(0);
      dataType minValue = tree->getValue<dataType>(0);
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i) {
        if(!tree->isRoot(i) and !tree->isLeaf(i)) {
          dataType iValue = tree->getValue<dataType>(i);
          if(mergeByPersistence) {
            maxValue = (maxValue < iValue) ? iValue : maxValue;
            minValue = (minValue > iValue) ? iValue : minValue;
          } else {
            ftm::idNode const parent = tree->getParentSafe(i);
            dataType parentValue = tree->getValue<dataType>(parent);
            dataType tempMax = std::max(iValue, parentValue);
            dataType tempMin = std::min(iValue, parentValue);
            if((tempMax - tempMin) > (maxValue - minValue)) {
              maxValue = tempMax;
              minValue = tempMin;
            }
          }
        }
      }
      double const epsilonOri = epsilon;
      epsilon = (maxValue - minValue) * epsilon / 100;

      // For Farthest Saddle option
      if(epsilon1UseFarthestSaddle_)
        epsilon = tree->getMaximumPersistence<dataType>() * epsilonOri / 100;
      bool isJT = tree->isJoinTree<dataType>();
      auto isFarthest = [&](ftm::idNode a, ftm::idNode b) {
        return (isJT
                and tree->getValue<dataType>(a) > tree->getValue<dataType>(b))
               or (not isJT
                   and tree->getValue<dataType>(a)
                         < tree->getValue<dataType>(b));
      };
      std::vector<ftm::idNode> farthestSaddle(tree->getNumberOfNodes());
      for(unsigned int i = 0; i < farthestSaddle.size(); ++i)
        farthestSaddle[i] = i;

      // --- Merge saddle
      // Create stack
      std::stack<int> nodeStack;
      std::queue<ftm::idNode> queue;
      queue.emplace(tree->getRoot());
      while(!queue.empty()) {
        ftm::idNode const node = queue.front();
        queue.pop();
        nodeStack.emplace(node);
        std::vector<ftm::idNode> children;
        tree->getChildren(node, children);
        for(auto child : children)
          queue.emplace(child);
      }
      // Iterate through nodes
      while(!nodeStack.empty()) {
        ftm::idNode const nodeId = nodeStack.top();
        nodeStack.pop();
        if(!tree->isRoot(nodeId) and !tree->isLeaf(nodeId)) {
          ftm::idNode const parentNodeId = tree->getParentSafe(nodeId);
          dataType nodeValue = tree->getValue<dataType>(nodeId);
          if(epsilon1UseFarthestSaddle_)
            nodeValue = tree->getValue<dataType>(farthestSaddle[nodeId]);
          dataType parentNodeValue = tree->getValue<dataType>(parentNodeId);
          dataType diffValue = std::max(nodeValue, parentNodeValue)
                               - std::min(nodeValue, parentNodeValue);
          if(diffValue <= epsilon) {
            ftm::idNode nodeIdToDelete, nodeIdToKeep;
            if(mergeByPersistence) {
              auto birthDeath1 = tree->getBirthDeath<dataType>(nodeId);
              auto birthDeath2 = tree->getBirthDeath<dataType>(parentNodeId);
              dataType pers1
                = std::get<1>(birthDeath1) - std::get<0>(birthDeath1);
              dataType pers2
                = std::get<1>(birthDeath2) - std::get<0>(birthDeath2);
              nodeIdToDelete = (pers1 > pers2) ? parentNodeId : nodeId;
              nodeIdToKeep = (pers1 > pers2) ? nodeId : parentNodeId;
              if(nodeIdToDelete == parentNodeId)
                nodeStack.emplace(nodeId);
            } else {
              nodeIdToDelete = nodeId;
              nodeIdToKeep = parentNodeId;
            }
            // Manage nodeMerged vector of vector
            for(auto node : treeNodeMerged[nodeIdToDelete]) {
              treeNodeMerged[nodeIdToKeep].push_back(node);
              if(isFarthest(farthestSaddle[nodeIdToKeep],
                            tree->getNode(node)->getOrigin()))
                farthestSaddle[nodeIdToKeep] = tree->getNode(node)->getOrigin();
            }
            treeNodeMerged[nodeIdToKeep].push_back(
              tree->getNode(nodeIdToDelete)->getOrigin());
            if(isFarthest(farthestSaddle[nodeIdToKeep], nodeIdToDelete))
              farthestSaddle[nodeIdToKeep] = nodeIdToDelete;
            treeNodeMerged[nodeIdToDelete].clear();
            // Delete node
            tree->deleteNode(nodeIdToDelete);
          }
        }
      }

      if(fullMerge) {
        auto root = tree->getRoot();
        tree->getNode(root)->setOrigin(root);
      }
    }

    // Epsilon 2 and 3 processing
    template <class dataType>
    void persistenceMerging(ftm::FTMTree_MT *tree,
                            double epsilon2,
                            double epsilon3 = 100) {
      bool fullMerge = (epsilon2 == 0);
      fullMerge &= useFullMerge_;
      epsilon2 /= 100;
      epsilon3 /= 100;
      dataType maxPers = tree->getMaximumPersistence<dataType>();

      std::queue<ftm::idNode> queue;
      queue.emplace(tree->getRoot());
      while(!queue.empty()) {
        ftm::idNode const node = queue.front();
        queue.pop();
        ftm::idNode const nodeParent = tree->getParentSafe(node);
        if(!tree->isRoot(node)) {
          const double nodePers = tree->getNodePersistence<dataType>(node);
          const double nodeParentPers
            = tree->getNodePersistence<dataType>(nodeParent);
          if(nodePers / nodeParentPers > epsilon2
             and nodePers / maxPers < epsilon3)
            tree->setParent(node, tree->getParentSafe(nodeParent));
        }
        std::vector<ftm::idNode> children;
        tree->getChildren(node, children);
        for(auto child : children)
          queue.emplace(child);
      }

      if(fullMerge) {
        auto root = tree->getRoot();
        if(tree->getNode(root)->getOrigin() != (int)root) {
          tree->setParent(tree->getNode(root)->getOrigin(), root);
          tree->getNode(root)->setOrigin(root);
        }
      }
    }

    template <class dataType>
    void keepMostImportantPairs(ftm::FTMTree_MT *tree, int n, bool useBD) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> pairs;
      tree->getPersistencePairsFromTree(pairs, useBD);
      n = std::max(n, 2); // keep at least 2 pairs
      int const index = std::max((int)(pairs.size() - n), 0);
      dataType threshold = std::get<2>(pairs[index]) * (1.0 - 1e-6)
                           / tree->getMaximumPersistence<dataType>() * 100.0;
      persistenceThresholding<dataType>(tree, threshold);
    }

    template <class dataType>
    void persistenceThresholding(ftm::FTMTree_MT *tree,
                                 double persistenceThresholdT,
                                 std::vector<ftm::idNode> &deletedNodes) {
      ftm::idNode const treeRoot = tree->getRoot();
      dataType maxPers = tree->getMaximumPersistence<dataType>();
      dataType threshold = persistenceThresholdT / 100 * maxPers;

      dataType secondMax = tree->getSecondMaximumPersistence<dataType>();
      bool keepOneZeroPersistencePair = (secondMax == 0 or maxPers == 0);
      if(threshold >= secondMax)
        threshold = (1.0 - 1e-6) * secondMax;

      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i) {
        if(tree->isRoot(i))
          continue;
        dataType nodePers = tree->getNodePersistence<dataType>(i);
        if(nodePers == 0 and keepOneZeroPersistencePair
           and tree->getParentSafe(i) == treeRoot) {
          keepOneZeroPersistencePair = false;
          continue;
        }
        if((nodePers == 0 or nodePers <= threshold
            or not tree->isNodeOriginDefined(i))) {
          tree->deleteNode(i);
          deletedNodes.push_back(i);
          ftm::idNode const nodeOrigin = tree->getNode(i)->getOrigin();
          if(tree->isNodeOriginDefined(i)
             and tree->getNode(nodeOrigin)->getOrigin() == (int)i) {
            tree->deleteNode(nodeOrigin);
            deletedNodes.push_back(nodeOrigin);
          }
        }
      }
    }

    template <class dataType>
    void persistenceThresholding(ftm::FTMTree_MT *tree,
                                 std::vector<ftm::idNode> &deletedNodes) {
      persistenceThresholding<dataType>(
        tree, persistenceThreshold_, deletedNodes);
    }

    template <class dataType>
    void persistenceThresholding(ftm::FTMTree_MT *tree,
                                 double persistenceThresholdT) {
      std::vector<ftm::idNode> deletedNodes;
      persistenceThresholding<dataType>(
        tree, persistenceThresholdT, deletedNodes);
    }

    template <class dataType>
    void persistenceThresholding(ftm::FTMTree_MT *tree) {
      std::vector<ftm::idNode> deletedNodes;
      persistenceThresholding<dataType>(
        tree, persistenceThreshold_, deletedNodes);
    }

    template <class dataType>
    void verifyOrigins(ftm::FTMTree_MT *tree) {
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
        if(not tree->isNodeAlone(i) and not tree->isNodeOriginDefined(i)) {
          std::stringstream ss;
          std::vector<ftm::idNode> children;
          tree->getChildren(i, children);
          ss << i << " has no origin (scalar=" << tree->getValue<dataType>(i)
             << ", noChildren=" << children.size()
             << ", parent=" << tree->getParentSafe(i) << ")";
          printMsg(ss.str());
          if(!tree->isRoot(i))
            tree->deleteNode(i);
          else {
            std::stringstream ss2;
            ss2 << "the root has no origin!";
            printErr(ss2.str());
          }
        }
    }

    template <class dataType>
    void preprocessTree(ftm::FTMTree_MT *tree,
                        bool deleteInconsistentNodes = true) {
      if(deleteInconsistentNodes) {
        // Manage inconsistent critical points
        // Critical points with same scalar value than parent
        for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
          if(!tree->isNodeAlone(i) and !tree->isRoot(i)
             and tree->getValue<dataType>(tree->getParentSafe(i))
                   == tree->getValue<dataType>(i)) {
            /*printMsg("[preprocessTree] " + std::to_string(i)
                     + " has same scalar value than parent (will be
               deleted).");*/
            tree->deleteNode(i);
          }
        // Valence 2 nodes
        for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
          if(tree->getNode(i)->getNumberOfUpSuperArcs() == 1
             and tree->getNode(i)->getNumberOfDownSuperArcs() == 1) {
            /*printMsg("[preprocessTree] " + std::to_string(i)
                     + " has 1 up arc and 1 down arc (will be deleted).");*/
            tree->deleteNode(i);
          }
      }

      // Compute persistence pairs
      if(not isPersistenceDiagram_ or convertToDiagram_) {
        auto pairs = ftm::computePersistencePairs<dataType>(tree);
        // Verify pairs
        verifyOrigins<dataType>(tree);
      }
    }

    template <class dataType>
    ftm::FTMTree_MT *computeBranchDecomposition(
      ftm::FTMTree_MT *tree,
      std::vector<std::vector<ftm::idNode>> &treeNodeMerged) {
      ftm::FTMTree_MT *treeNew = tree;

      ftm::idNode const root = treeNew->getRoot();

      // Manage when there is only one pair
      if(tree->isThereOnlyOnePersistencePair()) {
        ftm::idNode const rootOrigin = treeNew->getNode(root)->getOrigin();
        treeNew->getNode(rootOrigin)->setOrigin(rootOrigin);
        return treeNew;
      }

      // Manage multi persistence pairing
      std::vector<std::vector<ftm::idNode>> treeMultiPers;
      tree->getMultiPersOriginsVectorFromTree(treeMultiPers);

      // General case
      std::vector<bool> nodeDone(tree->getNumberOfNodes(), false);
      std::queue<ftm::idNode> queueNodes;
      queueNodes.emplace(root);
      while(!queueNodes.empty()) {
        ftm::idNode const node = queueNodes.front();
        queueNodes.pop();
        ftm::idNode const nodeOrigin = treeNew->getNode(node)->getOrigin();
        if(node == nodeOrigin
           or treeNew->getNodeLevel(node) > treeNew->getNodeLevel(nodeOrigin))
          continue;

        // Init vector with all origins
        std::vector<std::tuple<ftm::idNode, int>> vecOrigins;
        for(auto nodeMergedOrigin : treeNodeMerged[node]) {
          vecOrigins.emplace_back(nodeMergedOrigin, 0);
          for(auto multiPersOrigin :
              treeMultiPers[tree->getNode(nodeMergedOrigin)->getOrigin()])
            vecOrigins.emplace_back(multiPersOrigin, 1);
        }
        if(not tree->isNodeMerged(node))
          for(auto multiPersOrigin : treeMultiPers[node])
            vecOrigins.emplace_back(multiPersOrigin, 1);
        vecOrigins.emplace_back(nodeOrigin, 2);

        bool splitRoot = (vecOrigins.size() != 1 and treeNew->isRoot(node));
        splitRoot = false; // disabled

        // Process each origin
        for(auto stackTuple : vecOrigins) {
          ftm::idNode const nodeOriginT = std::get<0>(stackTuple);
          int const nodeOriginTID = std::get<1>(stackTuple);
          if(nodeDone[nodeOriginT]
             and nodeDone[tree->getNode(nodeOriginT)->getOrigin()])
            continue;
          nodeDone[nodeOriginT] = true;
          nodeDone[tree->getNode(nodeOriginT)->getOrigin()] = true;

          // Manage new parent
          ftm::idNode newParent = node;
          // - if merged node
          if(nodeOriginTID == 0) {
            newParent = treeNew->getNode(nodeOriginT)->getOrigin();
            treeNew->setParent(newParent, treeNew->getParentSafe(node));
            // - if multi pers node or nodeOrigin and splitRoot
          } else if(nodeOriginTID == 1 or (nodeOriginTID == 2 and splitRoot)) {
            newParent = nodeOriginT;
          }

          // Set nodes in the branch as childrens of the node
          ftm::idNode parentNodeOrigin = treeNew->getParentSafe(nodeOriginT);
          while(parentNodeOrigin != node) {
            ftm::idNode const oldParentNodeOrigin
              = treeNew->getParentSafe(parentNodeOrigin);
            treeNew->setParent(parentNodeOrigin, newParent);
            parentNodeOrigin = oldParentNodeOrigin;
          }

          if(nodeOriginTID == 1 or (nodeOriginTID == 2 and splitRoot))
            treeNew->setParent(newParent, treeNew->getParentSafe(node));
          else // if(nodeOriginTID != 1) // if not a multi pers node
            // Delete the other node of the pair
            treeNew->deleteNode(nodeOriginT);
          if(nodeOriginTID == 2 and splitRoot)
            tree->getNode(node)->setOrigin(node);

          // Push childrens of the node to the stack to process them
          std::vector<ftm::idNode> childrenNode;
          treeNew->getChildren(newParent, childrenNode);
          for(ftm::idNode const children : childrenNode)
            if(!treeNew->isLeaf(children))
              queueNodes.emplace(children);
        }
      }

      // Verify inconsistency
      // verifyBranchDecompositionInconsistency<dataType>(treeNew);

      return treeNew;
    }

    template <class dataType>
    void dontUseMinMaxPair(ftm::FTMTree_MT *tree) {
      ftm::idNode const treeRoot = tree->getRoot();
      // Full merge case, search for the origin
      if(tree->getNode(treeRoot)->getOrigin() == (int)treeRoot) {
        ftm::idNode const nodeIdToDelete
          = tree->getMergedRootOrigin<dataType>();
        if(nodeIdToDelete != treeRoot
           and not tree->isNodeIdInconsistent(nodeIdToDelete)) {
          if(tree->isThereOnlyOnePersistencePair())
            tree->getNode(nodeIdToDelete)->setOrigin(nodeIdToDelete);
          else
            tree->deleteNode(nodeIdToDelete);
        }
        // Classic case
      } else {
        ftm::idNode const rootOrigin = tree->getNode(treeRoot)->getOrigin();
        if(tree->isThereOnlyOnePersistencePair())
          tree->getNode(rootOrigin)->setOrigin(rootOrigin);
        else
          tree->deleteNode(rootOrigin);
      }

      tree->getNode(treeRoot)->setOrigin(treeRoot);
    }

    void verifyPairsTree(ftm::FTMTree_MT *tree) {
      int cptBug = 0;
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i) {
        if(not tree->isNodeOriginDefined(i)) {
          std::stringstream ss;
          ss << i << " _ " << tree->getNode(i)->getOrigin() << " / "
             << tree->getNumberOfNodes();
          printMsg(ss.str());
          if(tree->isNodeAlone(i))
            printMsg("alone");
          cptBug++;
        }
      }
      std::stringstream ss;
      ss << cptBug;
      printMsg(ss.str());
    }

    template <class dataType>
    void deleteMultiPersPairs(ftm::FTMTree_MT *tree, bool useBD) {
      auto multiPersOrigins = tree->getMultiPersOrigins<dataType>(useBD);
      for(auto origin : multiPersOrigins)
        tree->deleteNode(origin);
    }

    template <class dataType>
    void preprocessingPipeline(ftm::MergeTree<dataType> &mTree,
                               double epsilonTree,
                               double epsilon2Tree,
                               double epsilon3Tree,
                               bool branchDecompositionT,
                               bool useMinMaxPairT,
                               bool cleanTreeT,
                               double persistenceThreshold,
                               std::vector<int> &nodeCorr,
                               bool deleteInconsistentNodes = true) {
      Timer t_proc;

      ftm::FTMTree_MT *tree = &(mTree.tree);

      preprocessTree<dataType>(tree, deleteInconsistentNodes);

      // - Delete null persistence pairs and persistence thresholding
      persistenceThresholding<dataType>(tree, persistenceThreshold);

      // - Merge saddle points according epsilon
      std::vector<std::vector<ftm::idNode>> treeNodeMerged(
        tree->getNumberOfNodes());
      if(not isPersistenceDiagram_ or convertToDiagram_) {
        if(epsilonTree != 0)
          mergeSaddle<dataType>(tree, epsilonTree, treeNodeMerged);
      }

      // - Compute branch decomposition
      // verifyPairsTree(tree);
      if(branchDecompositionT
         and (not isPersistenceDiagram_ or convertToDiagram_))
        tree = computeBranchDecomposition<dataType>(tree, treeNodeMerged);

      // - Delete multi pers pairs
      if(deleteMultiPersPairs_)
        deleteMultiPersPairs<dataType>(tree, branchDecompositionT);

      // - Remove min max pair
      // verifyPairsTree(tree);
      if(not useMinMaxPairT)
        dontUseMinMaxPair<dataType>(tree);

      // - Epsilon 2 and 3 processing
      if(branchDecompositionT and not isPersistenceDiagram_)
        persistenceMerging<dataType>(tree, epsilon2Tree, epsilon3Tree);

      // - Tree cleaning (remove unused nodes)
      if(cleanTreeT) {
        ftm::cleanMergeTree<dataType>(mTree, nodeCorr, branchDecompositionT);
        tree = &(mTree.tree);
        reverseNodeCorr(tree, nodeCorr);
      }

      // - Root number verification
      if(tree->getNumberOfRoot() != 1)
        printErr("preprocessingPipeline tree->getNumberOfRoot() != 1");

      // - Time printing
      // verifyPairsTree(tree);
      auto t_preproc_time = t_proc.getElapsedTime();
      std::stringstream ss;
      ss << "TIME PREPROC.   = " << t_preproc_time;
      printMsg(ss.str(), debug::Priority::VERBOSE);
    }

    template <class dataType>
    void preprocessingPipeline(ftm::MergeTree<dataType> &mTree,
                               double epsilonTree,
                               double epsilon2Tree,
                               double epsilon3Tree,
                               bool branchDecompositionT,
                               bool useMinMaxPairT,
                               bool cleanTreeT,
                               std::vector<int> &nodeCorr,
                               bool deleteInconsistentNodes = true) {
      preprocessingPipeline<dataType>(
        mTree, epsilonTree, epsilon2Tree, epsilon3Tree, branchDecompositionT,
        useMinMaxPairT, cleanTreeT, persistenceThreshold_, nodeCorr,
        deleteInconsistentNodes);
    }

    void reverseNodeCorr(ftm::FTMTree_MT *tree, std::vector<int> &nodeCorr) {
      std::vector<int> newNodeCorr(tree->getNumberOfNodes());
      for(unsigned int i = 0; i < nodeCorr.size(); ++i)
        if(nodeCorr[i] >= 0 && nodeCorr[i] < (int)newNodeCorr.size())
          newNodeCorr[nodeCorr[i]] = i;
      nodeCorr = newNodeCorr;
    }

    template <class dataType>
    void mtFlattening(ftm::MergeTree<dataType> &mt) {
      ftm::FTMTree_MT *tree = &(mt.tree);
      ttk::ftm::computePersistencePairs<dataType>(tree);
      persistenceThresholding<dataType>(tree);
      std::vector<std::vector<ftm::idNode>> treeNodeMerged;
      mergeSaddle<dataType>(tree, 100.0, treeNodeMerged);
      computeBranchDecomposition<dataType>(tree, treeNodeMerged);
    }

    template <class dataType>
    void mtsFlattening(std::vector<ftm::MergeTree<dataType>> &mts) {
      for(auto &mt : mts)
        mtFlattening(mt);
    }

    double getSizeLimitMetric(std::vector<ftm::FTMTree_MT *> &trees) {
      std::array<double, 3> stats;
      getTreesStats(trees, stats);
      auto meanNodes = stats[0];
      unsigned int const n = trees.size();
      return meanNodes * n;
    }

    // ------------------------------------------------------------------------
    // Tree Postprocessing
    // ------------------------------------------------------------------------
    template <class dataType>
    void copyMinMaxPair(ftm::MergeTree<dataType> &mTree1,
                        ftm::MergeTree<dataType> &mTree2,
                        bool setOrigins = false) {
      // Get min max pair
      ftm::FTMTree_MT *tree1 = &(mTree1.tree);
      ftm::idNode root = tree1->getRoot();
      dataType newMax = tree1->getValue<dataType>(root);
      ftm::idNode rootOrigin = tree1->getNode(root)->getOrigin();
      dataType newMin = tree1->getValue<dataType>(rootOrigin);

      // Update tree
      ftm::idNode root2 = mTree2.tree.getRoot();
      std::vector<dataType> newScalarsVector;
      ftm::getTreeScalars<dataType>(mTree2, newScalarsVector);
      newScalarsVector[root2] = newMax;

      auto root2Origin = mTree2.tree.getNode(root2)->getOrigin();
      if(root2Origin == (int)root2)
        root2Origin = mTree2.tree.template getMergedRootOrigin<dataType>();
      if(mTree2.tree.isNodeIdInconsistent(root2Origin))
        newScalarsVector.push_back(newMin);
      else
        newScalarsVector[root2Origin] = newMin;

      // Set new scalars
      ftm::setTreeScalars<dataType>(mTree2, newScalarsVector);

      // Create root origin if not already there
      ftm::FTMTree_MT *treeNew = &(mTree2.tree);
      if(mTree2.tree.isNodeIdInconsistent(root2Origin)) {
        root2Origin = treeNew->getNumberOfNodes();
        treeNew->makeNode(root2Origin);
      }

      // Manage new origins
      if(setOrigins) {
        treeNew->getNode(root2Origin)->setOrigin(root2);
        treeNew->getNode(root2)->setOrigin(root2Origin);
      }
    }

    template <class dataType>
    std::tuple<int, dataType> fixMergedRootOrigin(ftm::FTMTree_MT *tree) {
      if(not tree->isFullMerge())
        return std::make_tuple(-1, -1);

      // Get node of the min max pair
      int maxIndex = tree->getMergedRootOrigin<dataType>();

      // Link node of the min max pair with the root
      ftm::idNode const treeRoot = tree->getRoot();
      dataType oldOriginValue
        = tree->getValue<dataType>(tree->getNode(maxIndex)->getOrigin());
      tree->getNode(maxIndex)->setOrigin(treeRoot);

      return std::make_tuple(maxIndex, oldOriginValue);
    }

    // TODO fix bug when one multi pers. pairs is moved up with epsilon 2 and 3
    // but not its brothers
    template <class dataType>
    void branchDecompositionToTree(ftm::FTMTree_MT *tree) {
      ftm::idNode const treeRoot = tree->getRoot();

      // Get original tree message
      std::stringstream const oriPrintTree = tree->printTree();
      std::stringstream const oriPrintPairs
        = tree->printPairsFromTree<dataType>(true);
      std::stringstream const oriPrintMultiPers
        = tree->printMultiPersPairsFromTree<dataType>(true);

      // One pair case
      if(tree->isThereOnlyOnePersistencePair()) {
        ftm::idNode const treeRootOrigin = tree->getNode(treeRoot)->getOrigin();
        tree->getNode(treeRootOrigin)->setOrigin(treeRoot);
        return;
      }

      // Manage full merge and dontuseMinMaxPair_
      bool const isFM = tree->isFullMerge();
      if(isFM) {
        ftm::idNode const mergedRootOrigin
          = tree->getMergedRootOrigin<dataType>();
        if(not tree->isNodeIdInconsistent(mergedRootOrigin)
           and mergedRootOrigin != treeRoot)
          tree->getNode(treeRoot)->setOrigin(mergedRootOrigin);
        else {
          printErr("branchDecompositionToTree mergedRootOrigin inconsistent");
        }
      }

      // Some functions
      bool isJT = tree->isJoinTree<dataType>();
      auto comp = [&](const std::tuple<ftm::idNode, dataType> &a,
                      const std::tuple<ftm::idNode, dataType> &b) {
        return isJT ? std::get<1>(a) > std::get<1>(b)
                    : std::get<1>(a) < std::get<1>(b);
      };
      auto getIndexNotMultiPers = [&](int index, ftm::FTMTree_MT *treeT,
                                      std::vector<ftm::idNode> &children) {
        while(index >= 0 and treeT->isMultiPersPair(children[index]))
          --index;
        return index;
      };

      // Branch Decomposition To Tree
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> nodeParent;
      std::queue<ftm::idNode> queue;
      queue.emplace(treeRoot);
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();
        auto nodeOrigin = tree->getNode(node)->getOrigin();
        if(tree->isLeaf(node)) {
          if(tree->isNodeAlone(nodeOrigin)) {
            if(not isFM)
              nodeParent.emplace_back(nodeOrigin, node);
            else
              nodeParent.emplace_back(node, nodeOrigin);
          } else if(tree->isMultiPersPair(node)) {
            nodeParent.emplace_back(node, nodeOrigin);
          }
          continue;
        }

        // Get children and sort them by scalar values
        std::vector<ftm::idNode> childrenOri;
        tree->getChildren(node, childrenOri);
        std::vector<ftm::idNode> children = childrenOri;
        std::vector<std::tuple<ftm::idNode, dataType>> childrenScalars;
        for(unsigned int i = 0; i < children.size(); ++i) {
          if(isFM and (int) children[i] != nodeOrigin)
            children[i] = tree->getNode(children[i])->getOrigin();
          childrenScalars.push_back(std::make_tuple(
            children[i], tree->getValue<dataType>(children[i])));
        }
        std::sort(std::begin(childrenScalars), std::end(childrenScalars), comp);
        children.clear();
        for(unsigned int i = 0; i < childrenScalars.size(); ++i)
          children.push_back(std::get<0>(childrenScalars[i]));

        // Get new parent of children
        for(unsigned int i = 1; i < children.size(); ++i) {
          if(tree->isMultiPersPair(children[i]))
            continue;
          int const index = getIndexNotMultiPers(i - 1, tree, children);
          if(index >= 0)
            nodeParent.emplace_back(children[i], children[index]);
        }

        bool const multiPersPair
          = tree->getNode(nodeOrigin)->getOrigin() != (int)node;
        if(not multiPersPair) {
          if(not isFM) {
            int const index
              = getIndexNotMultiPers(children.size() - 1, tree, children);
            nodeParent.emplace_back(nodeOrigin, children[index]);
          } else
            nodeParent.emplace_back(children[0], node);
        } else {
          // std::cout << "branchDecompositionToTree multiPersPair" <<
          // std::endl;
          nodeParent.emplace_back(children[0], nodeOrigin);
          int index = getIndexNotMultiPers(children.size() - 1, tree, children);
          if(index < 0) { // should not be possible
            printErr("[branchDecompositionToTree] index < 0");
            index = 0;
          }
          nodeParent.emplace_back(node, children[index]);
        }

        // Push children to the queue
        for(auto child : childrenOri)
          queue.emplace(child);
      }

      // Set new parents for each node
      for(auto nodeParentT : nodeParent)
        tree->setParent(std::get<0>(nodeParentT), std::get<1>(nodeParentT));

      // Verify that the tree is correct
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
        if(tree->getNode(i)->getNumberOfDownSuperArcs() == 1
           and tree->getNode(i)->getNumberOfUpSuperArcs() == 1) {
          printMsg(oriPrintPairs.str());
          printMsg(oriPrintMultiPers.str());
          printMsg(oriPrintTree.str());
          printMsg(tree->printTree().str());
          std::stringstream ss;
          auto iOrigin = tree->getNode(i)->getOrigin();
          ss << i << " _ " << iOrigin;
          if(tree->getNode(iOrigin)->getOrigin() != int(i))
            ss << " _ " << tree->getNode(iOrigin)->getOrigin() << " _ "
               << tree->getNode(tree->getNode(iOrigin)->getOrigin())
                    ->getOrigin();
          printMsg(ss.str());
          printErr("[branchDecompositionToTree] 1 up arc and 1 down arc");
        }
    }

    // For not branch decomposition tree
    template <class dataType>
    void putBackMergedNodes(ftm::FTMTree_MT *tree) {
      bool isJT = tree->isJoinTree<dataType>();
      std::queue<ftm::idNode> queue;
      queue.emplace(tree->getRoot());
      while(!queue.empty()) {
        ftm::idNode const node = queue.front();
        queue.pop();
        ftm::idNode const nodeOrigin = tree->getNode(node)->getOrigin();
        if(!tree->isLeaf(node)) {
          std::vector<ftm::idNode> children;
          tree->getChildren(node, children);
          std::vector<ftm::idNode> lowestNodes;
          ftm::idNode branchOrigin = nodeOrigin;
          for(auto child : children) {
            ftm::idNode lowestNode = tree->getLowestNode<dataType>(child);
            lowestNodes.push_back(lowestNode);
            ftm::idNode const lowestNodeOrigin
              = tree->getNode(lowestNode)->getOrigin();
            if(not tree->isNodeAlone(lowestNodeOrigin)
               and lowestNodeOrigin != node)
              branchOrigin = lowestNode;
          }
          for(size_t i = 0; i < children.size(); ++i) {
            ftm::idNode lowestNodeOrigin
              = tree->getNode(lowestNodes[i])->getOrigin();
            if(branchOrigin == lowestNodes[i] or lowestNodeOrigin == node)
              continue;
            dataType lowestNodeOriginVal
              = tree->getValue<dataType>(lowestNodeOrigin);
            ftm::idNode branchOriginT = branchOrigin;
            ftm::idNode const branchRoot
              = tree->getNode(branchOrigin)->getOrigin();
            while(branchRoot != branchOriginT) {
              dataType val
                = tree->getValue<dataType>(tree->getParentSafe(branchOriginT));
              if((val > lowestNodeOriginVal and isJT)
                 or (val < lowestNodeOriginVal and not isJT))
                break;
              branchOriginT = tree->getParentSafe(branchOriginT);
            }
            tree->setParent(
              lowestNodeOrigin, tree->getParentSafe(branchOriginT));
            tree->setParent(branchOriginT, lowestNodeOrigin);
            tree->setParent(children[i], lowestNodeOrigin);
          }
        }
        std::vector<ftm::idNode> children;
        tree->getChildren(node, children);
        for(auto child : children)
          queue.emplace(child);
      }
    }

    template <class dataType>
    void postprocessingPipeline(ftm::FTMTree_MT *tree) {
      // if(not branchDecomposition_ or not useMinMaxPair)
      // fixMergedRootOrigin<dataType>(tree);
      if(tree->isFullMerge()) {
        auto mergedRootOrigin = tree->getMergedRootOrigin<dataType>();
        if(not tree->isNodeIdInconsistent(mergedRootOrigin))
          tree->getNode(tree->getRoot())->setOrigin(mergedRootOrigin);
        else
          printErr(
            "[postprocessingPipeline] mergedRootOrigin inconsistent id.");
      }
      if(branchDecomposition_) {
        if(not isPersistenceDiagram_ and tree->getRealNumberOfNodes() != 0)
          branchDecompositionToTree<dataType>(tree);
      } else
        putBackMergedNodes<dataType>(tree);
    }

    // ------------------------------------------------------------------------
    // Output Matching
    // ------------------------------------------------------------------------
    template <class dataType>
    void convertBranchDecompositionMatching(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
        &outputMatching) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> toAdd;
      for(auto mTuple : outputMatching) {
        ftm::idNode const node1 = std::get<0>(mTuple);
        ftm::idNode const node2 = std::get<1>(mTuple);
        double const cost = std::get<2>(mTuple);
        ftm::idNode const node1Origin = tree1->getNode(node1)->getOrigin();
        ftm::idNode const node2Origin = tree2->getNode(node2)->getOrigin();

        int const node1Level = tree1->getNodeLevel(node1);
        int const node1OriginLevel = tree1->getNodeLevel(node1Origin);
        int const node2Level = tree2->getNodeLevel(node2);
        int const node2OriginLevel = tree2->getNodeLevel(node2Origin);

        ftm::idNode const node1Higher
          = (node1Level > node1OriginLevel) ? node1 : node1Origin;
        ftm::idNode const node1Lower
          = (node1Level > node1OriginLevel) ? node1Origin : node1;
        ftm::idNode const node2Higher
          = (node2Level > node2OriginLevel) ? node2 : node2Origin;
        ftm::idNode const node2Lower
          = (node2Level > node2OriginLevel) ? node2Origin : node2;

        if(((tree1->isRoot(node1Higher) and tree1->isFullMerge())
            or (tree2->isRoot(node2Higher) and tree2->isFullMerge())))
          continue;

        if(!tree1->isNodeAlone(node1Higher)
           and !tree2->isNodeAlone(node2Higher))
          toAdd.emplace_back(node1Higher, node2Higher, cost);
        if(!tree1->isNodeAlone(node1Lower) and !tree2->isNodeAlone(node2Lower))
          toAdd.emplace_back(node1Lower, node2Lower, cost);
      }
      outputMatching.clear();
      outputMatching.insert(outputMatching.end(), toAdd.begin(), toAdd.end());
    }

    template <class dataType>
    void convertBranchDecompositionMatching(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
        realOutputMatching(outputMatching.size());
      for(size_t i = 0; i < outputMatching.size(); ++i) {
        const auto &tup{outputMatching[i]};
        realOutputMatching[i] = {std::get<0>(tup), std::get<1>(tup), 0.0};
      }

      convertBranchDecompositionMatching<dataType>(
        tree1, tree2, realOutputMatching);

      outputMatching.clear();
      for(auto tup : realOutputMatching)
        outputMatching.emplace_back(std::get<0>(tup), std::get<1>(tup));
    }

    template <class dataType>
    void identifyRealMatching(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, bool>> &realMatching) {
      for(std::tuple<ftm::idNode, ftm::idNode> mTuple : outputMatching) {
        ftm::idNode tree1Node = std::get<0>(mTuple);
        ftm::idNode tree2Node = std::get<1>(mTuple);
        dataType relabelCostVal
          = relabelCostOnly<dataType>(tree1, tree1Node, tree2, tree2Node);
        dataType deleteInsertCostVal = deleteCost<dataType>(tree1, tree1Node)
                                       + insertCost<dataType>(tree2, tree2Node);
        bool isRealMatching = (relabelCostVal <= deleteInsertCostVal);
        realMatching.emplace_back(tree1Node, tree2Node, isRealMatching);
      }
    }

    // ------------------------------------------------------------------------
    // Edit Costs
    // ------------------------------------------------------------------------
    template <class dataType>
    dataType computeDistance(
      dataType x1, dataType x2, dataType y1, dataType y2, double power = 2) {
      if(power <= 0)
        return std::max(
          std::abs((double)(x1 - y1)), std::abs((double)(x2 - y2)));
      else
        return std::pow(std::abs((double)(x1 - y1)), power)
               + std::pow(std::abs((double)(x2 - y2)), power);
    }

    template <class dataType>
    dataType deleteCost(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
      dataType cost = 0;
      dataType newMin = 0.0, newMax = 1.0;
      // Get birth/death
      auto birthDeath
        = normalizedWasserstein_
            ? getNormalizedBirthDeath<dataType>(tree, nodeId, newMin, newMax)
            : tree->getBirthDeath<dataType>(nodeId);
      dataType birth = std::get<0>(birthDeath);
      dataType death = std::get<1>(birthDeath);
      dataType projec = (birth + death) / 2;
      // Compute delete cost
      cost = computeDistance<dataType>(
        birth, death, projec, projec, wassersteinPower_);
      // Divide cost by two if not branch decomposition and not merged
      /*if(! branchDecomposition_ and ! tree->isNodeMerged(nodeId))
        cost /= 2;*/
      cost *= nonMatchingWeight_;

      return cost;
    }

    template <class dataType>
    dataType insertCost(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
      return deleteCost<dataType>(tree, nodeId);
    }

    template <class dataType>
    dataType relabelCostOnly(ftm::FTMTree_MT *tree1,
                             ftm::idNode nodeId1,
                             ftm::FTMTree_MT *tree2,
                             ftm::idNode nodeId2) {
      dataType cost = 0;
      dataType newMin = 0.0, newMax = 1.0;
      // Get birth/death of the first tree
      auto birthDeath1
        = normalizedWasserstein_
            ? getNormalizedBirthDeath<dataType>(tree1, nodeId1, newMin, newMax)
            : tree1->getBirthDeath<dataType>(nodeId1);
      dataType birth1 = std::get<0>(birthDeath1);
      dataType death1 = std::get<1>(birthDeath1);
      // Get birth/death of the second tree
      auto birthDeath2
        = normalizedWasserstein_
            ? getNormalizedBirthDeath<dataType>(tree2, nodeId2, newMin, newMax)
            : tree2->getBirthDeath<dataType>(nodeId2);
      dataType birth2 = std::get<0>(birthDeath2);
      dataType death2 = std::get<1>(birthDeath2);
      // Compute relabel cost
      cost = computeDistance<dataType>(
        birth1, death1, birth2, death2, wassersteinPower_);
      // Divide cost by two if not branch decomposition and not merged
      /*bool merged = isNodeMerged(tree1, nodeId1) or isNodeMerged(tree2,
      nodeId2); if(! branchDecomposition_ and ! merged) cost /= 2;*/

      return cost;
    }

    template <class dataType>
    dataType relabelCost(ftm::FTMTree_MT *tree1,
                         ftm::idNode nodeId1,
                         ftm::FTMTree_MT *tree2,
                         ftm::idNode nodeId2) {
      // Full merge case and only one persistence pair case
      if(tree1->getNode(nodeId1)->getOrigin() == (int)nodeId1
         or tree2->getNode(nodeId2)->getOrigin() == (int)nodeId2)
        return 0;

      // Compute relabel cost
      dataType cost = relabelCostOnly<dataType>(tree1, nodeId1, tree2, nodeId2);

      if(keepSubtree_) {
        // Compute deleteInsert cost
        dataType deleteInsertCost = deleteCost<dataType>(tree1, nodeId1)
                                    + insertCost<dataType>(tree2, nodeId2);
        if(deleteInsertCost < cost)
          cost = deleteInsertCost;
      }

      return cost;
    }

    // ------------------------------------------------------------------------
    // Utils
    // ------------------------------------------------------------------------
    void getParamNames(std::vector<std::string> &paramNames) {
      paramNames = std::vector<std::string>{"epsilon1",
                                            "epsilon2",
                                            "epsilon3",
                                            "persistenceThreshold",
                                            "branchDecomposition",
                                            "normalizedWasserstein",
                                            "keepSubtree",
                                            "isPersistenceDiagram",
                                            "deleteMultiPersPairs",
                                            "epsilon1UseFarthestSaddle",
                                            "mixtureCoefficient"};
    }

    double getParamValueFromName(std::string &paramName) {
      double value = 0.0;
      if(paramName == "epsilon1")
        value = epsilonTree1_;
      else if(paramName == "epsilon2")
        value = epsilon2Tree1_;
      else if(paramName == "epsilon3")
        value = epsilon3Tree1_;
      else if(paramName == "persistenceThreshold")
        value = persistenceThreshold_;
      else if(paramName == "branchDecomposition")
        value = branchDecomposition_;
      else if(paramName == "normalizedWasserstein")
        value = normalizedWasserstein_;
      else if(paramName == "keepSubtree")
        value = keepSubtree_;
      else if(paramName == "isPersistenceDiagram")
        value = isPersistenceDiagram_;
      else if(paramName == "deleteMultiPersPairs")
        value = deleteMultiPersPairs_;
      else if(paramName == "epsilon1UseFarthestSaddle")
        value = epsilon1UseFarthestSaddle_;
      else if(paramName == "mixtureCoefficient")
        value = mixtureCoefficient_;
      return value;
    }

    void setParamValueFromName(std::string &paramName, double value) {
      if(paramName == "epsilon1")
        epsilonTree1_ = value;
      else if(paramName == "epsilon2")
        epsilon2Tree1_ = value;
      else if(paramName == "epsilon3")
        epsilon3Tree1_ = value;
      else if(paramName == "persistenceThreshold")
        persistenceThreshold_ = value;
      else if(paramName == "branchDecomposition")
        branchDecomposition_ = value;
      else if(paramName == "normalizedWasserstein")
        normalizedWasserstein_ = value;
      else if(paramName == "keepSubtree")
        keepSubtree_ = value;
      else if(paramName == "isPersistenceDiagram")
        isPersistenceDiagram_ = value;
      else if(paramName == "deleteMultiPersPairs")
        deleteMultiPersPairs_ = value;
      else if(paramName == "epsilon1UseFarthestSaddle")
        epsilon1UseFarthestSaddle_ = value;
      else if(paramName == "mixtureCoefficient")
        mixtureCoefficient_ = value;
    }

    void getTreesStats(std::vector<ftm::FTMTree_MT *> &trees,
                       std::array<double, 3> &stats) {
      double avgNodes = 0, avgNodesT = 0;
      double avgDepth = 0;
      for(unsigned int i = 0; i < trees.size(); ++i) {
        auto noNodesT = trees[i]->getNumberOfNodes();
        auto noNodes = trees[i]->getRealNumberOfNodes();
        avgNodes += noNodes;
        avgNodesT += noNodesT;
        avgDepth += trees[i]->getTreeDepth();
      }
      avgNodes /= trees.size();
      avgNodesT /= trees.size();
      avgDepth /= trees.size();

      stats = {avgNodes, avgNodesT, avgDepth};
    }

    void printTreesStats(std::vector<ftm::FTMTree_MT *> &trees) {
      std::array<double, 3> stats;
      getTreesStats(trees, stats);
      int avgNodes = stats[0], avgNodesT = stats[1];
      double const avgDepth = stats[2];
      std::stringstream ss;
      ss << trees.size() << " trees average [node: " << avgNodes << " / "
         << avgNodesT << ", depth: " << avgDepth << "]";
      printMsg(ss.str());
    }

    template <class dataType>
    void printTreesStats(std::vector<ftm::MergeTree<dataType>> &trees) {
      std::vector<ftm::FTMTree_MT *> treesT;
      ftm::mergeTreeToFTMTree<dataType>(trees, treesT);
      printTreesStats(treesT);
    }

    template <class dataType>
    void printTableVector(std::vector<std::vector<dataType>> &table) {
      std::streamsize const ssize = std::cout.precision();
      std::stringstream ss;
      ss << "      ";
      for(unsigned int j = 0; j < table[0].size(); ++j)
        ss << j - 1 << "    ";
      printMsg(ss.str());
      ss.str("");
      ss.clear();
      for(unsigned int i = 0; i < table.size(); ++i) {
        ss << std::setw(3) << std::setfill('0') << std::internal << i - 1
           << " | ";
        for(unsigned int j = 0; j < table[0].size(); ++j) {
          ss << std::fixed << std::setprecision(2) << table[i][j] << " ";
        }
        printMsg(ss.str());
        printMsg("");
      }
      std::cout.precision(ssize);
      printMsg(debug::Separator::L2);
    }

    template <class dataType>
    void printTable(dataType *table, int nRows, int nCols) {
      std::vector<std::vector<dataType>> vec(nRows, std::vector<dataType>());
      for(int i = 0; i < nRows; ++i)
        for(int j = 0; j < nCols; ++j)
          vec[i].push_back(table[i * nCols + j]);
      printTableVector<dataType>(vec);
    }

    void printMatching(std::vector<MatchingType> &matchings) {
      printMsg(debug::Separator::L2);
      for(const auto &mTuple : matchings) {
        std::stringstream ss;
        ss << std::get<0>(mTuple) << " - " << std::get<1>(mTuple) << " - "
           << std::get<2>(mTuple);
        printMsg(ss.str());
      }
      printMsg(debug::Separator::L2);
    }

    void printMatching(
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings) {
      printMsg(debug::Separator::L2);
      for(auto mTuple : matchings) {
        std::stringstream ss;
        ss << std::get<0>(mTuple) << " - " << std::get<1>(mTuple);
        printMsg(ss.str());
      }
      printMsg(debug::Separator::L2);
    }

    void printMatching(
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> &matchings) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matchingsT(
        matchings.size());
      for(size_t i = 0; i < matchings.size(); ++i) {
        const auto &tup{matchings[i]};
        matchingsT[i] = {std::get<0>(tup), std::get<1>(tup), 0.0};
      }
      printMatching(matchingsT);
    }

    template <class dataType>
    void printPairs(
      std::vector<std::tuple<SimplexId, SimplexId, dataType>> &treePairs) {
      for(auto pair : treePairs) {
        std::stringstream const ss;
        ss << std::get<0>(pair) << " _ " << std::get<1>(pair) << " _ "
           << std::get<2>(pair);
        printMsg(ss.str());
      }
      printMsg(debug::Separator::L2);
    }

    template <class dataType>
    void printOutputMatching(
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching,
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      bool computeCosts = true) {
      dataType cost = 0;
      std::vector<bool> tree1Done(tree1->getNumberOfNodes(), false);
      std::vector<bool> tree2Done(tree2->getNumberOfNodes(), false);
      std::stringstream ss;
      for(std::tuple<ftm::idNode, ftm::idNode> matching : outputMatching) {
        ftm::idNode node0 = std::get<0>(matching);
        ftm::idNode node1 = std::get<1>(matching);
        ftm::idNode node0Origin = tree1->getNode(node0)->getOrigin();
        ftm::idNode node1Origin = tree2->getNode(node1)->getOrigin();
        ss << node0 << " - " << node1 << " _ [ ";
        ss << "f(" << node0 << ")=" << tree1->getValue<dataType>(node0)
           << " _ ";
        ss << "g(" << node1 << ")=" << tree2->getValue<dataType>(node1) << " ]"
           << " _ [ ";
        ss << "f(" << node0Origin
           << ")=" << tree1->getValue<dataType>(node0Origin) << " _ ";
        ss << "g(" << node1Origin
           << ")=" << tree2->getValue<dataType>(node1Origin) << " ] ";

        if(computeCosts) {
          dataType tempCost = relabelCost<dataType>(tree1, node0, tree2, node1);
          dataType tempCost2
            = relabelCostOnly<dataType>(tree1, node0, tree2, node1);
          ss << "cost = " << tempCost << " (" << tempCost2 << ")" << std::endl;
          cost += tempCost;
        } else
          ss << std::endl;
        tree1Done[node0] = true;
        tree2Done[node1] = true;
      }

      for(unsigned int i = 0; i < tree1->getNumberOfNodes(); ++i)
        if(not tree1Done[i] and not tree1->isNodeAlone(i)) {
          ftm::idNode nodeOrigin = tree1->getNode(i)->getOrigin();
          ss << "T1 " << i << " _ [ f(" << i
             << ") = " << tree1->getValue<dataType>(i);
          ss << "_ f(" << nodeOrigin
             << ") = " << tree1->getValue<dataType>(nodeOrigin);
          ss << "]";
          if(computeCosts) {
            dataType tempCost = deleteCost<dataType>(tree1, i);
            ss << " _ cost = " << tempCost << std::endl;
            cost += tempCost;
          } else
            ss << std::endl;
        }
      for(unsigned int i = 0; i < tree2->getNumberOfNodes(); ++i)
        if(not tree2Done[i] and not tree2->isNodeAlone(i)) {
          ftm::idNode nodeOrigin = tree2->getNode(i)->getOrigin();
          ss << "T2 " << i << " _ [ g(" << i
             << ") = " << tree2->getValue<dataType>(i);
          ss << "_ g(" << nodeOrigin
             << ") = " << tree2->getValue<dataType>(nodeOrigin);
          ss << "]";
          if(computeCosts) {
            dataType tempCost = deleteCost<dataType>(tree2, i);
            ss << " _ cost = " << tempCost << std::endl;
            cost += tempCost;
          } else
            ss << std::endl;
        }
      if(computeCosts)
        ss << "total cost = " << cost << " (" << std::sqrt(cost) << ")"
           << std::endl;

      printMsg(ss.str());
      printMsg(debug::Separator::L2);
    }
  }; // MergeTreeBase class

} // namespace ttk
