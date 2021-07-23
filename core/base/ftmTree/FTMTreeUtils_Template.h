/// \ingroup base
/// \class FTMTreeUtils_Template
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _FTMTREEUTILS_TEMPLATE_H
#define _FTMTREEUTILS_TEMPLATE_H

#pragma once

namespace ttk {
  namespace ftm {

    // --------------------
    // Is
    // --------------------
    template <class dataType>
    bool FTMTree_MT::isJoinTree() {
      auto root = this->getRoot();
      std::vector<idNode> rootChildren;
      this->getChildren(root, rootChildren);
      idNode child = rootChildren[0];
      if(this->isFullMerge()) {
        dataType min = std::numeric_limits<dataType>::max();
        for(unsigned int i = 0; i < this->getNumberOfNodes(); ++i) {
          dataType value = this->getValue<dataType>(i);
          if(not this->isNodeAlone(i) and value < min) {
            min = value;
            child = i;
          }
        }
      }
      return this->getValue<dataType>(root) > this->getValue<dataType>(child);
    }

    template <class dataType>
    bool FTMTree_MT::isImportantPair(idNode nodeId, double threshold) {
      dataType rootPers = this->getNodePersistence<dataType>(this->getRoot());
      if(threshold > 1)
        threshold /= 100.0;
      threshold = rootPers * threshold;
      return this->getNodePersistence<dataType>(nodeId) > threshold;
    }

    template <class dataType>
    bool FTMTree_MT::isParentInconsistent(idNode nodeId) {
      auto parentBirthDeath
        = this->getBirthDeath<dataType>(this->getParentSafe(nodeId));
      dataType parentBirth = std::get<0>(parentBirthDeath);
      dataType parentDeath = std::get<1>(parentBirthDeath);
      auto birthDeath = this->getBirthDeath<dataType>(nodeId);
      dataType birth = std::get<0>(birthDeath);
      dataType death = std::get<1>(birthDeath);
      bool parentInconsistent = parentDeath < death or parentBirth > birth;
      return parentInconsistent;
    }

    template <class dataType>
    bool FTMTree_MT::verifyBranchDecompositionInconsistency() {
      bool inconsistency = false;
      std::queue<idNode> queue;
      queue.emplace(this->getRoot());
      while(!queue.empty()) {
        idNode node = queue.front();
        queue.pop();
        if(!this->isRoot(node) and this->isParentInconsistent<dataType>(node)) {
          printErr("inconsistency");
          this->printNode2<dataType>(node);
          this->printNode2<dataType>(this->getParentSafe(node));
          inconsistency = true;
        }
        std::vector<idNode> children;
        this->getChildren(node, children);
        for(idNode child : children)
          queue.emplace(child);
      }
      return inconsistency;
    }

    // --------------------
    // Get
    // --------------------
    template <class dataType>
    idNode FTMTree_MT::getLowestNode(idNode nodeStart) {
      idNode lowestNode = nodeStart;
      bool isJT = this->isJoinTree<dataType>();
      dataType bestVal = isJT ? std::numeric_limits<dataType>::max()
                              : std::numeric_limits<dataType>::lowest();
      std::queue<idNode> queue;
      queue.emplace(nodeStart);
      while(!queue.empty()) {
        idNode node = queue.front();
        queue.pop();
        dataType val = this->getValue<dataType>(node);
        if((val < bestVal and isJT) or (val > bestVal and not isJT)) {
          lowestNode = node;
          bestVal = val;
        }
        std::vector<idNode> children;
        this->getChildren(node, children);
        for(idNode child : children)
          queue.emplace(child);
      }
      return lowestNode;
    }

    // --------------------
    // Persistence
    // --------------------
    template <class dataType>
    std::tuple<dataType, dataType> FTMTree_MT::getBirthDeath(idNode nodeId) {
      idNode originId = this->getNode(nodeId)->getOrigin();
      if(this->isNodeOriginDefined(
           nodeId)) { // Avoid error if origin is not defined
        dataType pers1 = this->getValue<dataType>(nodeId);
        dataType pers2 = this->getValue<dataType>(originId);
        dataType birth = std::min(pers1, pers2);
        dataType death = std::max(pers1, pers2);
        return std::make_tuple(birth, death);
      }
      return std::make_tuple(0.0, 0.0);
    }

    template <class dataType>
    dataType FTMTree_MT::getBirth(idNode nodeId) {
      return std::get<0>(this->getBirthDeath<dataType>(nodeId));
    }

    template <class dataType>
    dataType FTMTree_MT::getNodePersistence(idNode nodeId) {
      std::tuple<dataType, dataType> birthDeath
        = this->getBirthDeath<dataType>(nodeId);
      return std::get<1>(birthDeath) - std::get<0>(birthDeath);
    }

    template <class dataType>
    dataType FTMTree_MT::getMaximumPersistence() {
      idNode root = this->getRoot();
      bool fullMerge = this->isFullMerge();

      // Classic case
      if(not fullMerge)
        return this->getNodePersistence<dataType>(this->getRoot());

      // Full merge case
      dataType maxPers = std::numeric_limits<dataType>::lowest();
      for(unsigned int i = 0; i < this->getNumberOfNodes(); ++i)
        if(not this->isNodeAlone(i) and this->isNodeOriginDefined(i)
           and this->getNode(i)->getOrigin() == (int)root)
          maxPers = std::max(maxPers, this->getNodePersistence<dataType>(i));

      return maxPers;
    }

    template <class dataType>
    dataType FTMTree_MT::getSecondMaximumPersistence() {
      idNode root = this->getRoot();
      dataType pers = std::numeric_limits<dataType>::lowest();
      for(unsigned int i = 0; i < this->getNumberOfNodes(); ++i) {
        if(this->isLeaf(i) and not this->isNodeAlone(i)
           and this->isNodeOriginDefined(i)) {
          idNode nodeOrigin = this->getNode(i)->getOrigin();
          if(not(nodeOrigin == root
                 and this->getNode(nodeOrigin)->getOrigin() == (int)i))
            pers = std::max(pers, this->getNodePersistence<dataType>(i));
        }
      }

      return pers;
    }

    template <class dataType>
    void FTMTree_MT::getPersistencePairsFromTree(
      std::vector<std::tuple<idNode, idNode, dataType>> &pairs, bool useBD) {
      std::vector<idNode> nodes;
      if(useBD) {
        for(unsigned int i = 0; i < this->getNumberOfNodes(); ++i)
          if(!this->isNodeAlone(i) and this->getNode(i)->getOrigin() != (int)i)
            nodes.push_back(i);
      } else
        this->getLeavesFromTree(nodes);
      for(auto node : nodes) {
        auto pers = this->getNodePersistence<dataType>(node);
        pairs.push_back(
          std::make_tuple(node, this->getNode(node)->getOrigin(), pers));
      }
      auto comp = [&](const std::tuple<idNode, idNode, dataType> a,
                      const std::tuple<idNode, idNode, dataType> b) {
        return std::get<2>(a) < std::get<2>(b);
      };
      sort(pairs.begin(), pairs.end(), comp);
    }

    template <class dataType>
    std::vector<idNode> FTMTree_MT::getMultiPersOrigins(bool useBD) {
      std::vector<idNode> multiPersOrigins;

      std::vector<std::tuple<idNode, idNode, dataType>> pairs;
      this->getPersistencePairsFromTree(pairs, useBD);
      // std::vector<idNode> origins(this->getNumberOfNodes(), -1);
      std::vector<std::vector<idNode>> origins(this->getNumberOfNodes());
      std::vector<bool> birthFound(this->getNumberOfNodes(), false);
      for(auto pair : pairs) {
        idNode nodeBirth = std::get<0>(pair);
        idNode nodeDeath = std::get<1>(pair);

        origins[nodeDeath].push_back(nodeBirth);
        birthFound[nodeBirth] = true;
      }

      for(unsigned int i = 0; i < origins.size(); ++i)
        if(birthFound[i])
          for(auto node : origins[i])
            multiPersOrigins.push_back(node);

      return multiPersOrigins;
    }

    // --------------------
    // Utils
    // --------------------
    template <class dataType>
    void FTMTree_MT::printNode2(idNode nodeId) {
      auto birthDeath = this->getBirthDeath<dataType>(nodeId);
      std::stringstream ss;
      ss << "nodeId = " << nodeId << " (" << std::get<0>(birthDeath)
         << ") _ originId = " << this->getNode(nodeId)->getOrigin() << " ("
         << std::get<1>(birthDeath) << ")";
      printMsg(ss.str());
    }

    template <class dataType>
    void FTMTree_MT::printTreeScalars(bool printNodeAlone) {
      std::streamsize sSize = std::cout.precision();
      for(unsigned int i = 0; i < this->getNumberOfNodes(); ++i) {
        idNode iOrigin
          = this->isNodeOriginDefined(i) ? this->getNode(i)->getOrigin() : i;
        if(printNodeAlone
           or (not printNodeAlone
               and (not this->isNodeAlone(i)
                    or not this->isNodeAlone(iOrigin)))) {
          std::stringstream ss;
          ss << i << " _ " << std::setprecision(12)
             << this->getValue<dataType>(i);
          printMsg(ss.str());
        }
      }
      printMsg(debug::Separator::L2);
      std::cout.precision(sSize);
    }

    template <class dataType>
    std::stringstream FTMTree_MT::printPairsFromTree(bool useBD,
                                                     bool printPairs,
                                                     bool doPrint) {
      std::stringstream ss;
      std::vector<std::tuple<idNode, idNode, dataType>> pairs;
      this->getPersistencePairsFromTree(pairs, useBD);
      ss << "size=" << pairs.size() << std::endl;
      if(printPairs)
        for(auto pair : pairs) {
          ss << std::get<0>(pair) << " ("
             << this->getValue<dataType>(std::get<0>(pair)) << ") _ ";
          ss << std::get<1>(pair) << " ("
             << this->getValue<dataType>(std::get<1>(pair)) << ") _ ";
          ss << std::get<2>(pair) << std::endl;
        }

      if(doPrint) {
        printMsg(ss.str());
        printMsg(debug::Separator::L2);
      }
      return ss;
    }

    template <class dataType>
    void FTMTree_MT::printMultiPersPairsFromTree(bool useBD, bool printPairs) {
      std::vector<std::tuple<idNode, idNode, dataType>> pairs;
      this->getPersistencePairsFromTree(pairs, useBD);
      std::vector<int> noOrigin(this->getNumberOfNodes(), 0);
      int noMultiPers = 0;
      for(auto pair : pairs) {
        noOrigin[std::get<0>(pair)]++;
        noMultiPers += (noOrigin[std::get<0>(pair)] > 1) ? 1 : 0;
        noOrigin[std::get<1>(pair)]++;
        noMultiPers += (noOrigin[std::get<1>(pair)] > 1) ? 1 : 0;
      }
      std::stringstream ss;
      ss << "Number of multi pers pairs : " << noMultiPers << std::endl;
      if(printPairs) {
        auto multiPers = this->getMultiPersOrigins<dataType>(useBD);
        for(auto node : multiPers)
          ss << node << std::endl;
      }
      printMsg(ss.str());
      printMsg(debug::Separator::L2);
    }

  } // namespace ftm
} // namespace ttk

#endif
