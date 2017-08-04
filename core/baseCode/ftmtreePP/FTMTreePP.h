/// \ingroup baseCode
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

namespace ttk
{
   namespace ftm
   {
      class FTMTreePP : public FTMTree
      {
        private:
         vector<AtomicUF*> nodesUF_;

        public:
         FTMTreePP();
         virtual ~FTMTreePP();

         template <typename scalarType>
         void computePersistencePairs(vector<tuple<idVertex, idVertex, scalarType>>& pairs,
                                      const bool jt);

        protected:
         template <typename scalarType>
         void computePairs(ftm::FTMTree_MT* tree,
                           vector<tuple<idVertex, idVertex, scalarType>>& pairs);

         template <typename scalarType>
         void sortPairs(ftm::FTMTree_MT* tree,
                        vector<tuple<idVertex, idVertex, scalarType>>& pairs);

         void addPendingNode(const idNode parent, const idNode toAdd)
         {
            // Trick, we use the arc list to maintaint nodes
            // coming to this UF.
            nodesUF_[parent]->addArcToClose(toAdd);
         }

         idNode countPendingNode(const idNode current)
         {
            return nodesUF_[current]->getOpenedArcs().size();
         }

         template <typename scalarType>
         idVertex getMostPersistVert(const idNode current, ftm::FTMTree_MT* tree)
         {
            idVertex  minVert = tree->getNode(current)->getVertexId();
            AtomicUF* uf      = nodesUF_[current];

            for (const auto nodeid : uf->getOpenedArcs()) {
               const idVertex vtmp = nodesUF_[nodeid]->getExtrema();
               if (tree->compLower(vtmp, minVert)) {
                  minVert = vtmp;
               }
            }

            return minVert;
         }

         void clearPendingNodes(const idNode current)
         {
            nodesUF_[current]->clearOpenedArcs();
         }

         template <typename scalarType>
         void createPairs(const idNode current,
                          vector<tuple<idVertex, idVertex, scalarType>>& pairs,
                          ftm::FTMTree_MT* tree,
                          const idVertex   mp)
         {
            AtomicUF*        uf      = nodesUF_[current];
            const idVertex   curVert = tree->getNode(current)->getVertexId();
            const scalarType curVal  = getValue<scalarType>(curVert);

            for (const auto nodeid : uf->getOpenedArcs()) {
               AtomicUF::makeUnion(uf, nodesUF_[nodeid]);
               const idVertex tmpVert = nodesUF_[nodeid]->getExtrema();
               if (tmpVert != mp) {
                  const scalarType tmpVal = getValue<scalarType>(tmpVert);
                  if (scalars_->isLower(tmpVert, curVert)) {
                     pairs.emplace_back(tmpVert, curVert, curVal - tmpVal);
                  } else {
                     pairs.emplace_back(tmpVert, curVert, tmpVal - curVal);
                  }
               }
            }
         }

      };
   }
}

template <typename scalarType>
void ftm::FTMTreePP::computePersistencePairs(vector<tuple<idVertex, idVertex, scalarType>>& pairs,
                                             const bool jt)
{
   ftm::FTMTree_MT* tree = jt ? getJoinTree() : getSplitTree();

   nodesUF_.clear();
   pairs.clear();
   nodesUF_.resize(tree->getNumberOfNodes(), nullptr);
   pairs.reserve(tree->getNumberOfLeaves());

   const auto nbNodes = tree->getNumberOfNodes();
   for (idNode nid = 0; nid < nbNodes; ++nid) {
      nodesUF_[nid] = new AtomicUF(tree->getNode(nid)->getVertexId());
   }

   computePairs<scalarType>(tree, pairs);

   sortPairs<scalarType>(tree, pairs);

   // destruct
   for (AtomicUF* uf : nodesUF_) {
      delete uf;
      uf = nullptr;
   }
}

template <typename scalarType>
void ftm::FTMTreePP::computePairs(ftm::FTMTree_MT* tree,
                                  vector<tuple<idVertex, idVertex, scalarType>>& pairs)
{
   auto getParentNode = [&](const idNode current) {
      const idSuperArc parentArc = tree->getNode(current)->getUpSuperArcId(0);
      return tree->getSuperArc(parentArc)->getUpNodeId();
   };

   std::queue<idNode> toSee;

   // start at the leaves
   const auto& vectLeaves = tree->getLeaves();
   for (const idNode nid : vectLeaves) {
      toSee.emplace(nid);
   }

   while (!toSee.empty()) {
      const idNode current = toSee.front();
      toSee.pop();

      if (!tree->getNode(current)->getNumberOfUpSuperArcs()) {
         createPairs<scalarType>(current, pairs, tree, ftm::nullVertex);
         clearPendingNodes(current);
         continue;
      } else {
         clearPendingNodes(current);
      }

      const idNode parent = getParentNode(current);
      addPendingNode(parent, current);

      if (countPendingNode(parent) == tree->getNode(parent)->getNumberOfDownSuperArcs()) {
         const idVertex mostPersist = getMostPersistVert<scalarType>(parent, tree);
         createPairs<scalarType>(parent, pairs, tree, mostPersist);
         nodesUF_[parent]->setExtrema(mostPersist);
         toSee.push(parent);
      }
   }
}

template <typename scalarType>
void ftm::FTMTreePP::sortPairs(ftm::FTMTree_MT* tree,
                               vector<tuple<idVertex, idVertex, scalarType>>& pairs)
{
   auto comp = [&](const tuple<idVertex, idVertex, scalarType> a,
                   const tuple<idVertex, idVertex, scalarType> b) {
       return get<2>(a) < get<2>(b);
   };

   sort(pairs.begin(), pairs.end(), comp);
}

#endif  // FTMTREE_PP_H
