/// \ingroup base
/// \class ttk::FTMTreePP
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2017-08-03
///
///\brief TTK processing package that add persistance pairs features
// to the contour tree package.
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa ttkPersistenceDiagram.cpp %for a usage example.

#ifndef FTMTREE_PP_H
#define FTMTREE_PP_H

#include "FTMTree.h"

namespace ttk {
  namespace ftm {
    /**
     * Compute the persistence pairs of a function on a triangulation.
     * TTK assumes that the input dataset is made of only one connected
     * component.
     */
    class FTMTreePP : public FTMTree {
    private:
      std::vector<AtomicUF *> nodesUF_;

    public:
      FTMTreePP();
      virtual ~FTMTreePP();

      template <typename scalarType>
      void computePersistencePairs(
        std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs,
        const bool jt);

    protected:
      template <typename scalarType>
      void computePairs(
        ftm::FTMTree_MT *tree,
        std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs);

      template <typename scalarType>
      void sortPairs(
        ftm::FTMTree_MT *tree,
        std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs);

      void addPendingNode(const idNode parentNode, const idNode toAdd) {
        // Trick, we use the arc list to maintaint nodes
        // coming to this UF.
        nodesUF_[parentNode]->find()->addArcToClose(toAdd);
      }

      idNode countPendingNode(const idNode current) {
        return nodesUF_[current]->find()->getOpenedArcs().size();
      }

      template <typename scalarType>
      SimplexId getMostPersistVert(const idNode current,
                                   ftm::FTMTree_MT *tree) {
        SimplexId minVert = tree->getNode(current)->getVertexId();
        AtomicUF *uf = nodesUF_[current]->find();

        for(const auto nodeid : uf->getOpenedArcs()) {
          const SimplexId vtmp = nodesUF_[nodeid]->find()->getExtrema();
          if(tree->compLower(vtmp, minVert)) {
            minVert = vtmp;
          }
        }

        return minVert;
      }

      void clearPendingNodes(const idNode current) {
        nodesUF_[current]->find()->clearOpenedArcs();
      }

      template <typename scalarType>
      void createPairs(
        const idNode current,
        std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs,
        ftm::FTMTree_MT *tree,
        const SimplexId mp) {
        AtomicUF *uf = nodesUF_[current]->find();
        const SimplexId curVert = tree->getNode(current)->getVertexId();
        const scalarType curVal = getValue<scalarType>(curVert);

        for(const auto nodeid : uf->getOpenedArcs()) {
          const SimplexId tmpVert = nodesUF_[nodeid]->find()->getExtrema();
          AtomicUF::makeUnion(uf, nodesUF_[nodeid]);
          if(tmpVert != mp) {
            const scalarType tmpVal = getValue<scalarType>(tmpVert);
            if(scalars_->isLower(tmpVert, curVert)) {
              pairs.emplace_back(tmpVert, curVert, curVal - tmpVal);
            } else {
              pairs.emplace_back(tmpVert, curVert, tmpVal - curVal);
            }
          }
        }
      }
    };
  } // namespace ftm
} // namespace ttk

template <typename scalarType>
void ttk::ftm::FTMTreePP::computePersistencePairs(
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs,
  const bool jt) {
  ftm::FTMTree_MT *tree = jt ? getJoinTree() : getSplitTree();

  nodesUF_.clear();
  pairs.clear();
  nodesUF_.resize(tree->getNumberOfNodes(), nullptr);
  pairs.reserve(tree->getNumberOfLeaves());

  const auto nbNodes = tree->getNumberOfNodes();
  for(idNode nid = 0; nid < nbNodes; ++nid) {
    nodesUF_[nid] = new AtomicUF(tree->getNode(nid)->getVertexId());
  }

  computePairs<scalarType>(tree, pairs);

  sortPairs<scalarType>(tree, pairs);

  // destruct
  for(AtomicUF *uf : nodesUF_) {
    delete uf;
    uf = nullptr;
  }
}

template <typename scalarType>
void ttk::ftm::FTMTreePP::computePairs(
  ftm::FTMTree_MT *tree,
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs) {
  auto getParentNode = [&](const idNode current) {
    const idSuperArc parentArc = tree->getNode(current)->getUpSuperArcId(0);
    return tree->getSuperArc(parentArc)->getUpNodeId();
  };

  ::std::queue<idNode> toSee;

  // start at the leaves
  auto &vectLeaves = tree->getLeaves();
  for(auto nid : vectLeaves) {
    toSee.emplace(nid);
  }

  while(!toSee.empty()) {
    const idNode current = toSee.front();
    toSee.pop();

    if(!tree->getNode(current)->getNumberOfUpSuperArcs()) {
      createPairs<scalarType>(current, pairs, tree, ftm::nullVertex);
      clearPendingNodes(current);
      continue;
    } else {
      clearPendingNodes(current);
    }

    const idNode parentNode = getParentNode(current);
    addPendingNode(parentNode, current);

    if(countPendingNode(parentNode)
       == tree->getNode(parentNode)->getNumberOfDownSuperArcs()) {
      const SimplexId mostPersist
        = getMostPersistVert<scalarType>(parentNode, tree);
      createPairs<scalarType>(parentNode, pairs, tree, mostPersist);
      nodesUF_[parentNode]->find()->setExtrema(mostPersist);
      toSee.push(parentNode);
    }
  }
}

template <typename scalarType>
void ttk::ftm::FTMTreePP::sortPairs(
  ftm::FTMTree_MT *tree,
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs) {
  auto comp = [&](const std::tuple<SimplexId, SimplexId, scalarType> a,
                  const std::tuple<SimplexId, SimplexId, scalarType> b) {
    return std::get<2>(a) < std::get<2>(b);
  };

  sort(pairs.begin(), pairs.end(), comp);
}

#endif // FTMTREE_PP_H
