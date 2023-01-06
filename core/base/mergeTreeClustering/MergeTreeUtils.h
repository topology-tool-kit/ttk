/// \ingroup base
/// \class MergeTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///

#pragma once

#include <FTMTreeUtils.h>
#include <stack>

namespace ttk {

  template <class dataType>
  bool isDeathHigher(ftm::FTMTree_MT *tree1,
                     ftm::idNode tree1Node,
                     ftm::FTMTree_MT *tree2,
                     ftm::idNode tree2Node) {
    dataType tree1NodeDeath
      = std::get<1>(tree1->getBirthDeath<dataType>(tree1Node));
    dataType tree2NodeDeath
      = std::get<1>(tree2->getBirthDeath<dataType>(tree2Node));
    return tree1NodeDeath > tree2NodeDeath;
  }

  // --------------------
  // Normalized Wasserstein
  // --------------------

  template <class dataType>
  dataType getMinMaxLocal(ftm::FTMTree_MT *tree,
                          ftm::idNode nodeId,
                          bool getMin = true) {
    auto nodeIdParent = tree->getParentSafe(nodeId);

    if(tree->notNeedToNormalize(nodeId))
      return getMin ? 0.0 : 1.0;

    // General case
    std::tuple<dataType, dataType> birthDeath
      = tree->getBirthDeath<dataType>(nodeIdParent);

    // Verify inconsistency
    if(tree->isParentInconsistent<dataType>(nodeId)) {
      // printErr("inconsistency");
      // if(tree->getNumberOfNodes() < 1000)
      tree->printTree();
      // printPairsFromTree<dataType>(tree, true);
      tree->printNode2<dataType>(nodeId);
      tree->printNode2<dataType>(nodeIdParent);
    }

    return getMin ? std::get<0>(birthDeath) : std::get<1>(birthDeath);
  }

  template <class dataType>
  std::tuple<dataType, dataType> getMinMaxLocalT(ftm::FTMTree_MT *tree,
                                                 ftm::idNode nodeId) {
    dataType min = getMinMaxLocal<dataType>(tree, nodeId);
    dataType max = getMinMaxLocal<dataType>(tree, nodeId, false);
    return std::make_tuple(min, max);
  }

  template <class dataType>
  std::tuple<double, double>
    getNormalizedBirthDeathDouble(ftm::FTMTree_MT *tree,
                                  ftm::idNode nodeId,
                                  dataType newMin = 0.0,
                                  dataType newMax = 1.0) {
    auto birthDeath = tree->getBirthDeath<dataType>(nodeId);
    double birth = std::get<0>(birthDeath);
    double death = std::get<1>(birthDeath);
    dataType shiftMin = getMinMaxLocal<dataType>(tree, nodeId);
    dataType shiftMax = getMinMaxLocal<dataType>(tree, nodeId, false);
    birth = (newMax - newMin) * (birth - shiftMin)
            / (shiftMax - shiftMin); // + newMin;
    death = (newMax - newMin) * (death - shiftMin)
            / (shiftMax - shiftMin); // + newMin;
    return std::make_tuple(birth, death);
  }

  template <class dataType>
  std::tuple<dataType, dataType> getNormalizedBirthDeath(ftm::FTMTree_MT *tree,
                                                         ftm::idNode nodeId,
                                                         dataType newMin = 0.0,
                                                         dataType newMax
                                                         = 1.0) {
    auto birthDeath = tree->getBirthDeath<dataType>(nodeId);
    dataType birth = std::get<0>(birthDeath);
    dataType death = std::get<1>(birthDeath);
    dataType shiftMin = getMinMaxLocal<dataType>(tree, nodeId);
    dataType shiftMax = getMinMaxLocal<dataType>(tree, nodeId, false);
    birth = (newMax - newMin) * (birth - shiftMin)
            / (shiftMax - shiftMin); // + newMin;
    death = (newMax - newMin) * (death - shiftMin)
            / (shiftMax - shiftMin); // + newMin;
    return std::make_tuple(birth, death);
  }

  // --------------------
  // Rescaled Wasserstein (old)
  // --------------------

  template <class dataType>
  std::tuple<dataType, dataType> getNewMinMax(ftm::FTMTree_MT *tree1,
                                              ftm::idNode nodeId1,
                                              ftm::FTMTree_MT *tree2,
                                              ftm::idNode nodeId2) {
    dataType shiftMin1 = getMinMaxLocal<dataType>(tree1, nodeId1);
    dataType shiftMax1 = getMinMaxLocal<dataType>(tree1, nodeId1, false);
    dataType shiftMin2 = getMinMaxLocal<dataType>(tree2, nodeId2);
    dataType shiftMax2 = getMinMaxLocal<dataType>(tree2, nodeId2, false);
    return std::make_tuple(
      (shiftMin1 + shiftMin2) / 2.0, (shiftMax1 + shiftMax2) / 2.0);
  }

  template <class dataType>
  std::tuple<dataType, dataType> getRescaledBirthDeath(ftm::FTMTree_MT *tree1,
                                                       ftm::idNode nodeId1,
                                                       ftm::FTMTree_MT *tree2,
                                                       ftm::idNode nodeId2) {
    auto newMinMax = getNewMinMax<dataType>(tree1, nodeId1, tree2, nodeId2);
    return getNormalizedBirthDeath<dataType>(
      tree1, nodeId1, std::get<0>(newMinMax), std::get<1>(newMinMax));
  }

  template <class dataType>
  dataType getMinMaxLocalFromVector(ftm::FTMTree_MT *tree,
                                    ftm::idNode nodeId,
                                    std::vector<dataType> &scalarsVector,
                                    bool getMin = true) {
    auto nodeIdParent = tree->getParentSafe(nodeId);

    if(tree->notNeedToNormalize(nodeId))
      return getMin ? 0.0 : 1.0;

    // General case
    dataType death = scalarsVector[nodeIdParent];
    dataType birth = scalarsVector[tree->getNode(nodeIdParent)->getOrigin()];
    if(death < birth) {
      dataType temp = death;
      death = birth;
      birth = temp;
    }

    return getMin ? birth : death;
  }

  template <class dataType>
  std::tuple<dataType, dataType>
    getMinMaxLocalFromVectorT(ftm::FTMTree_MT *tree,
                              ftm::idNode nodeId,
                              std::vector<dataType> &scalarsVector) {
    dataType min
      = getMinMaxLocalFromVector<dataType>(tree, nodeId, scalarsVector);
    dataType max
      = getMinMaxLocalFromVector<dataType>(tree, nodeId, scalarsVector, false);
    return std::make_tuple(min, max);
  }

  template <class dataType>
  std::tuple<dataType, dataType>
    getNewMinMaxFromVector(ftm::FTMTree_MT *tree1,
                           ftm::idNode nodeId1,
                           ftm::FTMTree_MT *tree2,
                           ftm::idNode nodeId2,
                           std::vector<dataType> &scalarsVector) {
    dataType shiftMin1 = getMinMaxLocal<dataType>(tree1, nodeId1);
    dataType shiftMax1 = getMinMaxLocal<dataType>(tree1, nodeId1, false);
    dataType shiftMin2
      = getMinMaxLocalFromVector<dataType>(tree2, nodeId2, scalarsVector);
    dataType shiftMax2 = getMinMaxLocalFromVector<dataType>(
      tree2, nodeId2, scalarsVector, false);
    return std::make_tuple(
      (shiftMin1 + shiftMin2) / 2.0, (shiftMax1 + shiftMax2) / 2.0);
  }

  template <class dataType>
  std::tuple<dataType, dataType>
    getRescaledBirthDeathFromVector(ftm::FTMTree_MT *tree1,
                                    ftm::idNode nodeId1,
                                    ftm::FTMTree_MT *tree2,
                                    ftm::idNode nodeId2,
                                    std::vector<dataType> &scalarsVector) {
    auto newMinMax = getNewMinMaxFromVector<dataType>(
      tree1, nodeId1, tree2, nodeId2, scalarsVector);
    return getNormalizedBirthDeath<dataType>(
      tree1, nodeId1, std::get<0>(newMinMax), std::get<1>(newMinMax));
  }

  // --------------------
  // Testing
  // --------------------
  template <class dataType>
  std::shared_ptr<ftm::FTMTree_MT>
    makeFakeTree(dataType *nodesScalar,
                 std::vector<SimplexId> &nodes,
                 std::vector<std::tuple<ftm::idNode, ftm::idNode>> &arcs) {
    // Init Scalars
    auto scalars = std::make_shared<ftm::Scalars>();
    scalars->size = nodes.size();
    scalars->values = (void *)nodesScalar;

    // Init Tree
    auto params = std::make_shared<ftm::Params>();
    auto tree
      = std::make_shared<ftm::FTMTree_MT>(params, scalars, ftm::Join_Split);
    tree->makeAlloc();

    // Add Nodes
    for(auto i : nodes)
      tree->makeNode(i);

    // Add Arcs
    for(std::tuple<ftm::idNode, ftm::idNode> arc : arcs)
      tree->makeSuperArc(std::get<0>(arc), std::get<1>(arc)); // (down, Up)

    return tree;
  }

  template <class dataType>
  ftm::MergeTree<dataType>
    makeFakeMergeTree(std::vector<dataType> &scalarsVector,
                      std::vector<std::tuple<ftm::idNode, ftm::idNode>> &arcs) {
    ftm::MergeTree<dataType> mergeTree
      = ttk::ftm::createEmptyMergeTree<dataType>(scalarsVector.size());
    ttk::ftm::setTreeScalars<dataType>(mergeTree, scalarsVector);
    ftm::FTMTree_MT *tree = &(mergeTree.tree);

    // Add Nodes
    for(unsigned int i = 0; i < scalarsVector.size(); ++i)
      tree->makeNode(i);

    // Add Arcs
    for(std::tuple<ftm::idNode, ftm::idNode> arc : arcs)
      tree->makeSuperArc(std::get<0>(arc), std::get<1>(arc)); // (down, Up)

    return mergeTree;
  }

} // namespace ttk
