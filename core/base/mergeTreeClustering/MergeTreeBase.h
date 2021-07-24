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

#ifndef _MERGETREEBASE_H
#define _MERGETREEBASE_H

#include <AssignmentSolver.h>
#include <FTMNode.h>
#include <FTMTree.h>
#include <FTMTreePPUtils.h>
#include <FTMTreeUtils.h>

#include "MergeTreeUtils.h"

using namespace ttk;

class MergeTreeBase : virtual public Debug {
protected:
  int assignmentSolverID_ = 0;
  bool epsilon1UseFarthestSaddle_ = false;
  double epsilonTree1_ = 5;
  double epsilonTree2_ = 5;
  double epsilon2Tree1_ = 95;
  double epsilon2Tree2_ = 95;
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

  bool distanceSquared_ = true;
  bool useFullMerge_ = false;

  // Old
  bool progressiveComputation_ = false;
  bool rescaledWasserstein_ = false;
  double normalizedWassersteinReg_ = 0.;
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
  };

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

  void setRescaledWasserstein(bool rescaledWasserstein) {
    rescaledWasserstein_ = rescaledWasserstein;
  }

  void setNormalizedWassersteinReg(double normalizedWassersteinReg) {
    normalizedWassersteinReg_ = normalizedWassersteinReg;
  }

  void setKeepSubtree(bool keepSubtree) {
    keepSubtree_ = keepSubtree;
  }

  void setProgressiveComputation(bool progressive) {
    progressiveComputation_ = progressive;
    /*if(progressiveComputation_)
      Preprocess = false;*/
  }

  void setBarycenterMergeTree(bool imt) {
    barycenterMergeTree_ = imt;
  }

  void setDistanceSquared(bool distanceSquared) {
    distanceSquared_ = distanceSquared;
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

  std::vector<std::vector<int>> getTreesNodeCorr() {
    return treesNodeCorr_;
  }

  // --------------------------------------------------------------------------------
  // Tree Preprocessing
  // --------------------------------------------------------------------------------
  // Epsilon 1 processing
  template <class dataType>
  void mergeSaddle(ftm::FTMTree_MT *tree,
                   double epsilon,
                   std::vector<std::vector<ftm::idNode>> &treeNodeMerged,
                   bool mergeByPersistence = false) {
    bool fullMerge = (epsilon == 100);
    fullMerge &= useFullMerge_;

    if(mergeByPersistence)
      computePersistencePairs<dataType>(
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
          ftm::idNode parent = tree->getParentSafe(i);
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
    double epsilonOri = epsilon;
    epsilon = (maxValue - minValue) * epsilon / 100;

    // For Farthest Saddle option
    if(epsilon1UseFarthestSaddle_)
      epsilon = tree->getMaximumPersistence<dataType>() * epsilonOri / 100;
    bool isJT = tree->isJoinTree<dataType>();
    auto isFarthest = [&](ftm::idNode a, ftm::idNode b) {
      return (isJT
              and tree->getValue<dataType>(a) > tree->getValue<dataType>(b))
             or (not isJT
                 and tree->getValue<dataType>(a) < tree->getValue<dataType>(b));
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
      ftm::idNode node = queue.front();
      queue.pop();
      nodeStack.emplace(node);
      std::vector<idNode> children;
      tree->getChildren(node, children);
      for(auto child : children)
        queue.emplace(child);
    }
    // Iterate through nodes
    while(!nodeStack.empty()) {
      ftm::idNode nodeId = nodeStack.top();
      nodeStack.pop();
      if(!tree->isRoot(nodeId) and !tree->isLeaf(nodeId)) {
        ftm::idNode parentNodeId = tree->getParentSafe(nodeId);
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
      ftm::idNode node = queue.front();
      queue.pop();
      ftm::idNode nodeParent = tree->getParentSafe(node);
      if(!tree->isRoot(node)) {
        dataType nodePers = tree->getNodePersistence<dataType>(node);
        dataType nodeParentPers
          = tree->getNodePersistence<dataType>(nodeParent);
        if(nodePers / nodeParentPers > epsilon2
           and nodePers / maxPers < epsilon3)
          tree->setParent(node, tree->getParentSafe(nodeParent));
      }
      std::vector<idNode> children;
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
  void persistenceThresholding(ftm::FTMTree_MT *tree,
                               double persistenceThresholdT,
                               std::vector<ftm::idNode> &deletedNodes) {
    dataType threshold
      = persistenceThresholdT / 100 * tree->getMaximumPersistence<dataType>();

    dataType secondMax = tree->getSecondMaximumPersistence<dataType>();
    if(threshold >= secondMax)
      threshold = 0.99 * secondMax;

    for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i) {
      dataType nodePers = tree->getNodePersistence<dataType>(i);
      if((nodePers == 0 or nodePers <= threshold
          or not tree->isNodeOriginDefined(i))
         and !tree->isRoot(i)) {
        tree->deleteNode(i);
        deletedNodes.push_back(i);
        ftm::idNode nodeOrigin = tree->getNode(i)->getOrigin();
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
  void verifyOrigins(ftm::FTMTree_MT *tree) {
    for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
      if(not tree->isNodeAlone(i) and not tree->isNodeOriginDefined(i)) {
        std::stringstream ss;
        std::vector<idNode> children;
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
                      double epsilon,
                      std::vector<std::vector<ftm::idNode>> &treeNodeMerged) {
    // Manage inconsistent critical points
    // Critical points with same scalar value than parent
    for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
      if(!tree->isNodeAlone(i) and !tree->isRoot(i)
         and tree->getValue<dataType>(tree->getParentSafe(i))
               == tree->getValue<dataType>(i))
        tree->deleteNode(i);
    // Valence 2 nodes
    for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
      if(tree->getNode(i)->getNumberOfUpSuperArcs() == 1
         and tree->getNode(i)->getNumberOfDownSuperArcs() == 1)
        tree->deleteNode(i);

    // Compute persistence pairs
    auto pairs = computePersistencePairs<dataType>(tree);
    // Verify pairs
    verifyOrigins<dataType>(tree);

    // Delete null persistence pairs and persistence thresholding
    std::vector<ftm::idNode> deletedNodes;
    persistenceThresholding<dataType>(tree, deletedNodes);

    // Merge saddle points according epsilon
    if(epsilon != 0)
      mergeSaddle<dataType>(tree, epsilon, treeNodeMerged);
  }

  template <class dataType>
  ftm::FTMTree_MT *computeBranchDecomposition(
    ftm::FTMTree_MT *tree,
    std::vector<std::vector<ftm::idNode>> &treeNodeMerged) {
    ftm::FTMTree_MT *treeNew = tree;

    ftm::idNode root = treeNew->getRoot();

    // Manage when there is only one pair
    if(tree->isThereOnlyOnePersistencePair()) {
      ftm::idNode rootOrigin = treeNew->getNode(root)->getOrigin();
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
      ftm::idNode node = queueNodes.front();
      queueNodes.pop();
      ftm::idNode nodeOrigin = treeNew->getNode(node)->getOrigin();
      if(node == nodeOrigin
         or treeNew->getNodeLevel(node) > treeNew->getNodeLevel(nodeOrigin))
        continue;

      // Init vector with all origins
      std::vector<std::tuple<ftm::idNode, int>> vecOrigins;
      for(auto nodeMergedOrigin : treeNodeMerged[node]) {
        vecOrigins.push_back(std::make_tuple(nodeMergedOrigin, 0));
        for(auto multiPersOrigin :
            treeMultiPers[tree->getNode(nodeMergedOrigin)->getOrigin()])
          vecOrigins.push_back(std::make_tuple(multiPersOrigin, 1));
      }
      if(not tree->isNodeMerged(node))
        for(auto multiPersOrigin : treeMultiPers[node])
          vecOrigins.push_back(std::make_tuple(multiPersOrigin, 1));
      vecOrigins.push_back(std::make_tuple(nodeOrigin, 2));

      bool splitRoot = (vecOrigins.size() != 1 and treeNew->isRoot(node));
      splitRoot = false; // disabled

      // Process each origin
      for(auto stackTuple : vecOrigins) {
        ftm::idNode nodeOriginT = std::get<0>(stackTuple);
        int nodeOriginTID = std::get<1>(stackTuple);
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
          ftm::idNode oldParentNodeOrigin
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
        for(ftm::idNode children : childrenNode)
          if(!treeNew->isLeaf(children))
            queueNodes.emplace(children);
      }
    }

    // Verify inconsistency
    // verifyBranchDecompositionInconsistency<dataType>(treeNew);

    return treeNew;
  }

  template <class dataType>
  ftm::idNode getMergedRootOrigin(ftm::FTMTree_MT *tree) {
    ftm::idNode treeRoot = tree->getRoot();
    bool isJt = tree->isJoinTree<dataType>();
    ftm::idNode mergedRootOrigin = treeRoot;
    for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
      if(tree->getNode(i)->getOrigin() == (int)treeRoot and i != treeRoot)
        if((isJt
            and tree->getValue<dataType>(i)
                  < tree->getValue<dataType>(mergedRootOrigin))
           or (not isJt
               and tree->getValue<dataType>(i)
                     > tree->getValue<dataType>(mergedRootOrigin)))
          mergedRootOrigin = i;
    return mergedRootOrigin;
  }

  template <class dataType>
  void dontUseMinMaxPair(ftm::FTMTree_MT *tree) {
    ftm::idNode treeRoot = tree->getRoot();
    // Full merge case, search for the origin
    if(tree->getNode(treeRoot)->getOrigin() == (int)treeRoot) {
      ftm::idNode nodeIdToDelete = getMergedRootOrigin<dataType>(tree);
      if(nodeIdToDelete != treeRoot
         and not tree->isNodeIdInconsistent(nodeIdToDelete)) {
        if(tree->isThereOnlyOnePersistencePair())
          tree->getNode(nodeIdToDelete)->setOrigin(nodeIdToDelete);
        else
          tree->deleteNode(nodeIdToDelete);
      }
      // Classic case
    } else {
      ftm::idNode rootOrigin = tree->getNode(treeRoot)->getOrigin();
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
  ftm::FTMTree_MT *preprocessingPipeline(MergeTree<dataType> &mTree,
                                         double epsilonTree,
                                         double epsilon2Tree,
                                         double epsilon3Tree,
                                         bool branchDecompositionT,
                                         bool useMinMaxPairT,
                                         bool cleanTreeT,
                                         std::vector<int> &nodeCorr) {
    Timer t_proc;

    FTMTree_MT *tree = &(mTree.tree);

    std::vector<std::vector<ftm::idNode>> treeNodeMerged(
      tree->getNumberOfNodes());

    preprocessTree<dataType>(tree, epsilonTree, treeNodeMerged);
    ftm::FTMTree_MT *treeOld = tree;

    // verifyPairsTree(tree);
    if(branchDecompositionT)
      tree = computeBranchDecomposition<dataType>(tree, treeNodeMerged);

    if(deleteMultiPersPairs_)
      deleteMultiPersPairs<dataType>(tree, branchDecompositionT);

    // verifyPairsTree(tree);
    if(not useMinMaxPairT)
      dontUseMinMaxPair<dataType>(tree);

    if(branchDecompositionT)
      persistenceMerging<dataType>(tree, epsilon2Tree, epsilon3Tree);

    if(cleanTreeT) {
      cleanMergeTree<dataType>(mTree, nodeCorr, branchDecompositionT);
      tree = &(mTree.tree);
      reverseNodeCorr(tree, nodeCorr);
    }

    if(tree->getNumberOfRoot() != 1) {
      printErr("preprocessingPipeline tree->getNumberOfRoot() != 1");
    }

    // verifyPairsTree(tree);
    auto t_preproc_time = t_proc.getElapsedTime();
    std::stringstream ss;
    ss << "TIME PREPROC.   = " << t_preproc_time;
    printMsg(ss.str(), debug::Priority::DETAIL);

    return treeOld;
  }

  void reverseNodeCorr(ftm::FTMTree_MT *tree, std::vector<int> &nodeCorr) {
    std::vector<int> newNodeCorr(tree->getNumberOfNodes());
    for(unsigned int i = 0; i < nodeCorr.size(); ++i)
      if(nodeCorr[i] >= 0 && nodeCorr[i] < (int)newNodeCorr.size())
        newNodeCorr[nodeCorr[i]] = i;
    nodeCorr = newNodeCorr;
  }

  // --------------------------------------------------------------------------------
  // Tree Postprocessing
  // --------------------------------------------------------------------------------
  template <class dataType>
  std::tuple<int, dataType> fixMergedRootOrigin(ftm::FTMTree_MT *tree) {
    if(not tree->isFullMerge())
      return std::make_tuple(-1, -1);

    // Get node of the min max pair
    dataType maxPers = 0;
    int maxIndex = -1;
    for(unsigned int j = 0; j < tree->getNumberOfNodes(); ++j) {
      if(not tree->isNodeAlone(j)) {
        dataType nodePers = tree->getNodePersistence<dataType>(j);
        if(nodePers > maxPers) {
          maxPers = nodePers;
          maxIndex = j;
        }
      }
    }

    // Link node of the min max pair with the root
    ftm::idNode treeRoot = tree->getRoot();
    dataType oldOriginValue
      = tree->getValue<dataType>(tree->getNode(maxIndex)->getOrigin());
    tree->getNode(maxIndex)->setOrigin(treeRoot);

    return std::make_tuple(maxIndex, oldOriginValue);
  }

  template <class dataType>
  void branchDecompositionToTree(ftm::FTMTree_MT *tree) {
    ftm::idNode treeRoot = tree->getRoot();

    // One pair case
    if(tree->isThereOnlyOnePersistencePair()) {
      idNode treeRootOrigin = tree->getNode(treeRoot)->getOrigin();
      tree->getNode(treeRootOrigin)->setOrigin(treeRoot);
      return;
    }

    // Manage full merge and dontuseMinMaxPair_
    bool isFM = tree->isFullMerge();
    if(isFM) {
      ftm::idNode mergedRootOrigin = getMergedRootOrigin<dataType>(tree);
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
      while(treeT->isMultiPersPair(children[index]))
        --index;
      index = std::max(0, index);
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
            nodeParent.push_back(std::make_tuple(nodeOrigin, node));
          else
            nodeParent.push_back(std::make_tuple(node, nodeOrigin));
        } else if(tree->isMultiPersPair(node)) {
          nodeParent.push_back(std::make_tuple(node, nodeOrigin));
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
        childrenScalars.push_back(
          std::make_tuple(children[i], tree->getValue<dataType>(children[i])));
      }
      std::sort(std::begin(childrenScalars), std::end(childrenScalars), comp);
      children.clear();
      for(unsigned int i = 0; i < childrenScalars.size(); ++i)
        children.push_back(std::get<0>(childrenScalars[i]));

      // Get new parent of children
      for(unsigned int i = 1; i < children.size(); ++i) {
        if(tree->isMultiPersPair(children[i]))
          continue;
        int index = getIndexNotMultiPers(i - 1, tree, children);
        nodeParent.push_back(std::make_tuple(children[i], children[index]));
      }

      bool multiPersPair = tree->getNode(nodeOrigin)->getOrigin() != (int)node;
      if(not multiPersPair) {
        if(not isFM) {
          int index = getIndexNotMultiPers(children.size() - 1, tree, children);
          nodeParent.push_back(std::make_tuple(nodeOrigin, children[index]));
        } else
          nodeParent.push_back(std::make_tuple(children[0], node));
      } else {
        // std::cout << "branchDecompositionToTree multiPersPair" << std::endl;
        nodeParent.push_back(std::make_tuple(children[0], nodeOrigin));
        int index = getIndexNotMultiPers(children.size() - 1, tree, children);
        nodeParent.push_back(std::make_tuple(node, children[index]));
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
        tree->printTree();
        std::stringstream ss;
        ss << i << " _ " << tree->getNode(i)->getOrigin();
        printMsg(ss.str());
        printErr("1 up arc and 1 down arc");
      }
  }

  // For not branch decomposition tree
  template <class dataType>
  void putBackMergedNodes(ftm::FTMTree_MT *tree) {
    bool isJT = tree->isJoinTree<dataType>();
    std::queue<ftm::idNode> queue;
    queue.emplace(tree->getRoot());
    while(!queue.empty()) {
      ftm::idNode node = queue.front();
      queue.pop();
      ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
      if(!tree->isLeaf(node)) {
        std::vector<ftm::idNode> children;
        tree->getChildren(node, children);
        std::vector<ftm::idNode> lowestNodes;
        ftm::idNode branchOrigin = nodeOrigin;
        for(auto child : children) {
          ftm::idNode lowestNode = tree->getLowestNode<dataType>(child);
          lowestNodes.push_back(lowestNode);
          ftm::idNode lowestNodeOrigin = tree->getNode(lowestNode)->getOrigin();
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
          ftm::idNode branchRoot = tree->getNode(branchOrigin)->getOrigin();
          while(branchRoot != branchOriginT) {
            dataType val
              = tree->getValue<dataType>(tree->getParentSafe(branchOriginT));
            if((val > lowestNodeOriginVal and isJT)
               or (val < lowestNodeOriginVal and not isJT))
              break;
            branchOriginT = tree->getParentSafe(branchOriginT);
          }
          tree->setParent(lowestNodeOrigin, tree->getParentSafe(branchOriginT));
          tree->setParent(branchOriginT, lowestNodeOrigin);
          tree->setParent(children[i], lowestNodeOrigin);
        }
      }
      std::vector<idNode> children;
      tree->getChildren(node, children);
      for(auto child : children)
        queue.emplace(child);
    }
  }

  template <class dataType>
  void postprocessingPipeline(ftm::FTMTree_MT *tree) {
    // TODO fix dont use min max pair
    if(not branchDecomposition_ or not useMinMaxPair_)
      fixMergedRootOrigin<dataType>(tree);
    if(branchDecomposition_)
      branchDecompositionToTree<dataType>(tree);
    else
      putBackMergedNodes<dataType>(tree);
  }

  // --------------------------------------------------------------------------------
  // Output Matching
  // --------------------------------------------------------------------------------
  template <class dataType>
  void computeMatching(
    ftm::FTMTree_MT *tree1,
    ftm::FTMTree_MT *tree2,
    std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
    std::vector<std::vector<std::vector<std::tuple<int, int>>>>
      &forestBackTable,
    std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &outputMatching,
    int startR,
    int startC) {
    std::queue<std::tuple<int, int, bool>> backQueue;
    backQueue.emplace(std::make_tuple(startR, startC, true));
    while(!backQueue.empty()) {
      std::tuple<int, int, bool> elem = backQueue.front();
      backQueue.pop();
      bool useTreeTable = std::get<2>(elem);
      int i = std::get<0>(elem);
      int j = std::get<1>(elem);

      if(useTreeTable) {
        int tupleI = std::get<0>(treeBackTable[i][j]);
        int tupleJ = std::get<1>(treeBackTable[i][j]);
        if(tupleI != 0 && tupleJ != 0) {
          useTreeTable = (tupleI != i || tupleJ != j);
          backQueue.emplace(std::make_tuple(tupleI, tupleJ, useTreeTable));
          if(not useTreeTable) { // We have matched i and j
            ftm::idNode tree1Node = tupleI - 1;
            ftm::idNode tree2Node = tupleJ - 1;
            double cost = 0;
            dataType costT
              = relabelCost<dataType>(tree1, tree1Node, tree2, tree2Node);
            cost = static_cast<double>(costT);
            outputMatching.push_back(
              std::make_tuple(tree1Node, tree2Node, cost));
          }
        }
      } else {
        for(std::tuple<int, int> forestBackElem : forestBackTable[i][j]) {
          int tupleI = std::get<0>(forestBackElem);
          int tupleJ = std::get<1>(forestBackElem);
          if(tupleI != 0 && tupleJ != 0) {
            useTreeTable = (tupleI != i && tupleJ != j);
            backQueue.emplace(std::make_tuple(tupleI, tupleJ, useTreeTable));
          }
        }
      }
    }
  }

  template <class dataType>
  void convertBranchDecompositionMatching(
    ftm::FTMTree_MT *tree1,
    ftm::FTMTree_MT *tree2,
    std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &outputMatching) {
    std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> toAdd;
    for(auto mTuple : outputMatching) {
      ftm::idNode node1 = std::get<0>(mTuple);
      ftm::idNode node2 = std::get<1>(mTuple);
      double cost = std::get<2>(mTuple);
      ftm::idNode node1Origin = tree1->getNode(node1)->getOrigin();
      ftm::idNode node2Origin = tree2->getNode(node2)->getOrigin();

      int node1Level = tree1->getNodeLevel(node1);
      int node1OriginLevel = tree1->getNodeLevel(node1Origin);
      int node2Level = tree2->getNodeLevel(node2);
      int node2OriginLevel = tree2->getNodeLevel(node2Origin);

      ftm::idNode node1Higher
        = (node1Level > node1OriginLevel) ? node1 : node1Origin;
      ftm::idNode node1Lower
        = (node1Level > node1OriginLevel) ? node1Origin : node1;
      ftm::idNode node2Higher
        = (node2Level > node2OriginLevel) ? node2 : node2Origin;
      ftm::idNode node2Lower
        = (node2Level > node2OriginLevel) ? node2Origin : node2;

      if(((tree1->isRoot(node1Higher) and tree1->isFullMerge())
          or (tree2->isRoot(node2Higher) and tree2->isFullMerge())))
        continue;

      if(!tree1->isNodeAlone(node1Higher) and !tree2->isNodeAlone(node2Higher))
        toAdd.push_back(std::make_tuple(node1Higher, node2Higher, cost));
      if(!tree1->isNodeAlone(node1Lower) and !tree2->isNodeAlone(node2Lower))
        toAdd.push_back(std::make_tuple(node1Lower, node2Lower, cost));
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
      realOutputMatching;
    for(auto tup : outputMatching)
      realOutputMatching.push_back(
        std::make_tuple(std::get<0>(tup), std::get<1>(tup), 0));

    convertBranchDecompositionMatching<dataType>(
      tree1, tree2, realOutputMatching);

    outputMatching.clear();
    for(auto tup : realOutputMatching)
      outputMatching.push_back(
        std::make_tuple(std::get<0>(tup), std::get<1>(tup)));
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
      realMatching.push_back(
        std::make_tuple(tree1Node, tree2Node, isRealMatching));
    }
  }

  // --------------------------------------------------------------------------------
  // Edit Costs
  // --------------------------------------------------------------------------------
  template <class dataType>
  dataType computeDistance(
    dataType x1, dataType x2, dataType y1, dataType y2, double power = 2) {
    if(power <= 0)
      return std::max(std::abs((double)(x1 - y1)), std::abs((double)(x2 - y2)));
    else
      return std::pow(std::abs((double)(x1 - y1)), power)
             + std::pow(std::abs((double)(x2 - y2)), power);
  }

  template <class dataType>
  dataType regularizationCost(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
    dataType shiftMin1 = getMinMaxLocal<dataType>(tree, nodeId);
    dataType shiftMax1 = getMinMaxLocal<dataType>(tree, nodeId, false);
    return computeDistance<dataType>(
      shiftMin1, 0, shiftMax1, 0, wassersteinPower_);
  }

  template <class dataType>
  dataType deleteCost(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
    dataType cost = 0;
    dataType newMin = 0.0, newMax = 1.0;
    if(normalizedWasserstein_ and rescaledWasserstein_) {
      auto newMinMax = getNewMinMax<dataType>(tree, nodeId, tree, nodeId);
      newMin = std::get<0>(newMinMax);
      newMax = std::get<1>(newMinMax);
    }
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
    // Regularize
    if(normalizedWasserstein_ and normalizedWassersteinReg_ != 0)
      cost += normalizedWassersteinReg_
              * regularizationCost<dataType>(tree, nodeId);

    return cost;
  }

  template <class dataType>
  dataType insertCost(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
    return deleteCost<dataType>(tree, nodeId);
  }

  template <class dataType>
  dataType regularizationCost2(ftm::FTMTree_MT *tree1,
                               ftm::idNode nodeId1,
                               ftm::FTMTree_MT *tree2,
                               ftm::idNode nodeId2) {
    dataType shiftMin1 = getMinMaxLocal<dataType>(tree1, nodeId1);
    dataType shiftMax1 = getMinMaxLocal<dataType>(tree1, nodeId1, false);
    dataType shiftMin2 = getMinMaxLocal<dataType>(tree2, nodeId2);
    dataType shiftMax2 = getMinMaxLocal<dataType>(tree2, nodeId2, false);
    return computeDistance<dataType>(
      shiftMin1, shiftMax1, shiftMin2, shiftMax2, wassersteinPower_);
  }

  template <class dataType>
  dataType relabelCostOnly(ftm::FTMTree_MT *tree1,
                           ftm::idNode nodeId1,
                           ftm::FTMTree_MT *tree2,
                           ftm::idNode nodeId2) {
    dataType cost = 0;
    dataType newMin = 0.0, newMax = 1.0;
    if(normalizedWasserstein_ and rescaledWasserstein_) {
      auto newMinMax = getNewMinMax<dataType>(tree1, nodeId1, tree2, nodeId2);
      newMin = std::get<0>(newMinMax);
      newMax = std::get<1>(newMinMax);
    }
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
    // Regularize
    if(normalizedWasserstein_ and normalizedWassersteinReg_ != 0)
      cost += normalizedWassersteinReg_
              * regularizationCost2<dataType>(tree1, nodeId1, tree2, nodeId2);

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
    // Compute deleteInsert cost
    dataType deleteInsertCost = deleteCost<dataType>(tree1, nodeId1)
                                + insertCost<dataType>(tree2, nodeId2);

    if(keepSubtree_ and deleteInsertCost < cost)
      cost = deleteInsertCost;

    return cost;
  }

  // --------------------------------------------------------------------------------
  // Edit Distance Dynamic Programming Equations
  // --------------------------------------------------------------------------------
  template <class dataType>
  void computeEquation8(ftm::FTMTree_MT *tree1,
                        ftm::idNode nodeI,
                        int i,
                        std::vector<std::vector<dataType>> &treeTable,
                        std::vector<std::vector<dataType>> &forestTable) {
    std::vector<ftm::idNode> children;
    tree1->getChildren(nodeI, children);
    forestTable[i][0] = 0;
    for(ftm::idNode child : children)
      forestTable[i][0] += treeTable[child + 1][0];
  }

  template <class dataType>
  void computeEquation9(ftm::FTMTree_MT *tree1,
                        ftm::idNode nodeI,
                        int i,
                        std::vector<std::vector<dataType>> &treeTable,
                        std::vector<std::vector<dataType>> &forestTable) {
    treeTable[i][0] = forestTable[i][0] + deleteCost<dataType>(tree1, nodeI);
  }

  template <class dataType>
  void computeEquation10(ftm::FTMTree_MT *tree2,
                         ftm::idNode nodeJ,
                         int j,
                         std::vector<std::vector<dataType>> &treeTable,
                         std::vector<std::vector<dataType>> &forestTable) {
    std::vector<ftm::idNode> children;
    tree2->getChildren(nodeJ, children);
    forestTable[0][j] = 0;
    for(ftm::idNode child : children)
      forestTable[0][j] += treeTable[0][child + 1];
  }

  template <class dataType>
  void computeEquation11(ftm::FTMTree_MT *tree2,
                         ftm::idNode nodeJ,
                         int j,
                         std::vector<std::vector<dataType>> &treeTable,
                         std::vector<std::vector<dataType>> &forestTable) {
    treeTable[0][j] = forestTable[0][j] + insertCost<dataType>(tree2, nodeJ);
  }

  // Compute first or second term of equation 12 or 13
  template <class dataType>
  std::tuple<dataType, ftm::idNode>
    computeTerm1_2(std::vector<ftm::idNode> &childrens,
                   int ind,
                   std::vector<std::vector<dataType>> &table,
                   bool computeTerm1) {
    dataType tempMin = (childrens.size() == 0)
                         ? ((computeTerm1) ? table[ind][0] : table[0][ind])
                         : std::numeric_limits<dataType>::max();
    ftm::idNode bestIdNode = 0;
    for(ftm::idNode children : childrens) {
      children += 1;
      dataType temp;
      if(computeTerm1) {
        temp = table[ind][children] - table[0][children];
      } else {
        temp = table[children][ind] - table[children][0];
      }
      if(temp < tempMin) {
        tempMin = temp;
        bestIdNode = children;
      }
    }
    return std::make_tuple(tempMin, bestIdNode);
  }

  template <class dataType>
  void computeEquation12(
    ftm::FTMTree_MT *tree1,
    ftm::FTMTree_MT *tree2,
    int i,
    int j,
    ftm::idNode nodeI,
    ftm::idNode nodeJ,
    std::vector<std::vector<dataType>> &treeTable,
    std::vector<std::vector<dataType>> &forestTable,
    std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
    std::vector<ftm::idNode> &children1,
    std::vector<ftm::idNode> &children2) {
    dataType treeTerm1, treeTerm2, treeTerm3;
    std::tuple<dataType, ftm::idNode> treeCoTerm1, treeCoTerm2;
    // Term 1
    treeCoTerm1 = computeTerm1_2<dataType>(children2, i, treeTable, true);
    treeTerm1 = treeTable[0][j] + std::get<0>(treeCoTerm1);

    // Term 2
    treeCoTerm2 = computeTerm1_2<dataType>(children1, j, treeTable, false);
    treeTerm2 = treeTable[i][0] + std::get<0>(treeCoTerm2);

    // Term 3
    treeTerm3
      = forestTable[i][j] + relabelCost<dataType>(tree1, nodeI, tree2, nodeJ);

    // Compute table value
    treeTable[i][j] = keepSubtree_
                        ? std::min(std::min(treeTerm1, treeTerm2), treeTerm3)
                        : treeTerm3;

    // Add backtracking information
    if(treeTable[i][j] == treeTerm3) {
      treeBackTable[i][j] = std::make_tuple(i, j);
    } else if(treeTable[i][j] == treeTerm2) {
      treeBackTable[i][j] = std::make_tuple(std::get<1>(treeCoTerm2), j);
    } else {
      treeBackTable[i][j] = std::make_tuple(i, std::get<1>(treeCoTerm1));
    }
  }

  // --------------------------------------------------------------------------------
  // Assignment Problem Functions
  // --------------------------------------------------------------------------------
  template <class dataType>
  dataType
    postprocessAssignment(std::vector<asgnMatchingTuple> &matchings,
                          std::vector<ftm::idNode> &children1,
                          std::vector<ftm::idNode> &children2,
                          std::vector<std::tuple<int, int>> &forestAssignment) {
    dataType cost = 0;
    for(asgnMatchingTuple mTuple : matchings) {
      cost += std::get<2>(mTuple);
      if(std::get<0>(mTuple) >= (int)children1.size()
         || std::get<1>(mTuple) >= (int)children2.size())
        continue;
      int tableId1 = children1[std::get<0>(mTuple)] + 1;
      int tableId2 = children2[std::get<1>(mTuple)] + 1;
      forestAssignment.push_back(std::make_tuple(tableId1, tableId2));
    }
    return cost;
  }

  // --------------------------------------------------------------------------------
  // Utils
  // --------------------------------------------------------------------------------
  void printTreesStats(std::vector<ftm::FTMTree_MT *> &trees) {
    int avgNodes = 0, avgNodesT = 0;
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
    std::stringstream ss;
    ss << trees.size() << " trees average [node: " << avgNodes << " / "
       << avgNodesT << ", depth: " << avgDepth << "]";
    printMsg(ss.str());
  }

  template <class dataType>
  void printTreesStats(std::vector<MergeTree<dataType>> &trees) {
    std::vector<ftm::FTMTree_MT *> treesT;
    mergeTreeToFTMTree<dataType>(trees, treesT);
    printTreesStats(treesT);
  }

  template <class dataType>
  void printTableVector(std::vector<std::vector<dataType>> &table) {
    std::streamsize ssize = std::cout.precision();
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

  void printMatching(std::vector<asgnMatchingTuple> &matchings) {
    printMsg(debug::Separator::L2);
    for(asgnMatchingTuple mTuple : matchings) {
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
    std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matchingsT;
    for(auto tup : matchings)
      matchingsT.push_back(
        std::make_tuple(std::get<0>(tup), std::get<1>(tup), 0));
    printMatching(matchingsT);
  }

  template <class dataType>
  void printPairs(
    std::vector<std::tuple<SimplexId, SimplexId, dataType>> &treePairs) {
    for(auto pair : treePairs) {
      std::stringstream ss;
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
      ss << "f(" << node0 << ")=" << tree1->getValue<dataType>(node0) << " _ ";
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

#endif
