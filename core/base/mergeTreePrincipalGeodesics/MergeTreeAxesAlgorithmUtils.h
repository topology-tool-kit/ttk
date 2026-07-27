#pragma once

#include <FTMTree_MT.h>
#include <string>

namespace ttk {
  namespace axa {
    //----------------------------------------------------------------------------
    // Utils
    //----------------------------------------------------------------------------
    // v[i] contains the node in tree matched to the node i in barycenter
    template <class dataType>
    void getMatchingVector(
      const ftm::MergeTree<dataType> &barycenter,
      const ftm::MergeTree<dataType> &tree,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings,
      std::vector<ftm::idNode> &matchingVector) {
      matchingVector.clear();
      matchingVector.resize(barycenter.tree.getNumberOfNodes(),
                            std::numeric_limits<ftm::idNode>::max());
      for(unsigned int j = 0; j < matchings.size(); ++j) {
        auto &match0 = std::get<0>(matchings[j]);
        auto &match1 = std::get<1>(matchings[j]);
        if(match0 < barycenter.tree.getNumberOfNodes()
           and match1 < tree.tree.getNumberOfNodes())
          matchingVector[match0] = match1;
      }
    }

    // v[i] contains the node in barycenter matched to the node i in tree
    template <class dataType>
    void getInverseMatchingVector(
      const ftm::MergeTree<dataType> &barycenter,
      const ftm::MergeTree<dataType> &tree,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings,
      std::vector<ftm::idNode> &matchingVector) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> invMatchings(
        matchings.size());
      for(unsigned int i = 0; i < matchings.size(); ++i)
        invMatchings[i] = std::make_tuple(std::get<1>(matchings[i]),
                                          std::get<0>(matchings[i]),
                                          std::get<2>(matchings[i]));
      getMatchingVector(tree, barycenter, invMatchings, matchingVector);
    }

    void reverseMatchingVector(unsigned int noNodes,
                               std::vector<ftm::idNode> &matchingVector,
                               std::vector<ftm::idNode> &invMatchingVector);

    template <class dataType>
    void reverseMatchingVector(ftm::MergeTree<dataType> &tree,
                               std::vector<ftm::idNode> &matchingVector,
                               std::vector<ftm::idNode> &invMatchingVector) {
      reverseMatchingVector(
        tree.tree.getNumberOfNodes(), matchingVector, invMatchingVector);
    }

    // m[i][j] contains the node in trees[j] matched to the node i in the
    // barycenter
    template <class dataType>
    void getMatchingMatrix(
      const ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<std::vector<ftm::idNode>> &matchingMatrix) {
      matchingMatrix.clear();
      matchingMatrix.resize(
        barycenter.tree.getNumberOfNodes(),
        std::vector<ftm::idNode>(
          trees.size(), std::numeric_limits<ftm::idNode>::max()));
      for(unsigned int i = 0; i < trees.size(); ++i) {
        std::vector<ftm::idNode> matchingVector;
        getMatchingVector<dataType>(
          barycenter, trees[i], matchings[i], matchingVector);
        for(unsigned int j = 0; j < matchingVector.size(); ++j)
          matchingMatrix[j][i] = matchingVector[j];
      }
    }

    //----------------------------------------------------------------------------
    // Output Utils
    //----------------------------------------------------------------------------
    void zeroPadding(std::string &colName,
                     const size_t numberCols,
                     const size_t colIdx);

    std::string getTableCoefficientName(int noAxes, int axeNum);

    std::string getTableCoefficientNormName(int noAxes, int axeNum);

    std::string getTableVectorName(
      int noAxes, int axeNum, int vId, int vComp, bool isSecondInput = false);

    std::string getTableCorrelationName(int noAxes, int axeNum);

    std::string getTableCorrelationPersName(int noAxes, int axeNum);

    std::string getTableTreeName(int noTrees, int treeNum);
  } // namespace axa
} // namespace ttk
