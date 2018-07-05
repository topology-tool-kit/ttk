/// \ingroup base
/// \class ttk::ftr::Graph
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-25
///
/// \brief TTK %FTRGraph graph skeleton
///
/// This class manage nodes and arcs of the graph structure
/// along with their Segmentation
///
/// \sa ttk::FTRGraph

#pragma once

#include "AtomicVector.h"
#include "DataTypesFTR.h"
#include "FTRCommon.h"
#include "Node.h"
#include "SuperArc.h"
#include "Scalars.h"
#include "Propagation.h"

#ifndef TTK_ENABLE_KAMIKAZE
#include<iostream>
#endif

#include <vector>

namespace ttk
{
   namespace ftr
   {

      class Graph : public Allocable
      {
        private:
         AtomicVector<std::tuple<idVertex, bool>> leaves_;
         AtomicVector<Node>                       nodes_;
         AtomicVector<SuperArc>                   arcs_;

         std::vector<AtomicVector<idSegmentation>> segmentation_;

        public:

         // For openmp capture, we need direct access...
         std::vector<valence> valDown_, valUp_;

         Graph();
         Graph(Graph&& other)      = default;
         Graph(const Graph& other) = delete;
         virtual ~Graph();

         Graph& operator=(Graph&& other) {
            if (this != &other) {
               leaves_       = std::move(other.leaves_);
               nodes_        = std::move(other.nodes_);
               arcs_         = std::move(other.arcs_);
               segmentation_ = std::move(other.segmentation_);
               valDown_      = std::move(other.valDown_);
               valUp_        = std::move(other.valUp_);
            }
            return *this;
         }

         Graph& operator=(Graph& other) = delete;

         // Accessor on structure
         // ---------------------

         idNode getNumberOfNodes(void) const
         {
            return nodes_.size();
         }

         idSuperArc getNumberOfArcs(void) const
         {
            return arcs_.size();
         }

         idSuperArc getNumberOfVisibleArcs(void) const
         {
            idSuperArc res = 0;
            for (const auto& arc : arcs_) {
               if (arc.isVisible())
                  ++res;
            }
            return res;
         }

         idNode getNumberOfLeaves(void) const
         {
            return leaves_.size();
         }

         idVertex getLeaf(const idNode id) const
         {
            return std::get<0>(leaves_[id]);
         }

         bool isLeafFromMin(const idNode id) const
         {
            return std::get<1>(leaves_[id]);
         }

         const Node& getNode(const idNode id) const
         {
           return nodes_[id];
         }

         Node& getNode(const idNode id)
         {
           return nodes_[id];
         }

         const SuperArc& getArc(const idSuperArc id) const
         {
            return arcs_[id];
         }

         SuperArc& getArc(const idSuperArc id)
         {
            return arcs_[id];
         }

         bool isVisited(const idVertex v) const
         {
            return !segmentation_[v].empty();
         }

         void visit(const idVertex v, const idSegmentation id, bool regular = true)
         {
            if (regular) {
               segmentation_[v].emplace_back(id);
            } else {
               segmentation_[v].emplace_back(-id-1);
            }
         }

         const AtomicVector<idSegmentation> visit(const idVertex v) const
         {
            return segmentation_[v];
         }

         bool hasVisited(const idVertex v, const idSuperArc id) const
         {
            for(const idSegmentation tmp :  segmentation_[v]){
               if (tmp == (idSegmentation)id) {
                  return true;
               }
            }
            return false;
         }

         bool haveSeen(const idVertex v, const idSuperArc rpz) const
         {
            for(const idSegmentation tmp :  segmentation_[v]){
               if (tmp >= 0 && getArc(tmp).getPropagation() &&
                   getArc(tmp).getPropagation()->getRpz() == rpz) {
                  return true;
               }
            }
            return false;
         }

         bool hasArcEndingHere(const idVertex saddle, const idVertex regular, VertCompFN comp) const
         {
            const idNode lookupNode = getNodeId(saddle);
            for(const idSegmentation tmp :  segmentation_[regular]){
               if (tmp >= 0 && getArc(tmp).getUpNodeId() == lookupNode) {
                  const idNode downNodeId = getArc(tmp).getDownNodeId();
                  const idVertex downVert = getNode(downNodeId).getVertexIdentifier();
                  if (comp(saddle, downVert)) {
                     return true;
                  }
               }
            }
            return false;
         }

         bool isNode(const idVertex v) const
         {
            for (const idSegmentation tmp : segmentation_[v]) {
               if (tmp < 0) {
                  return true;
               }
            }
            return false;
         }

         bool isArc(const idVertex v) const
         {
            for (const idSegmentation tmp : segmentation_[v]) {
               if (tmp >= 0) {
                  return true;
               }
            }
            return false;
         }

         bool isNode(const idVertex v, const idNode id) const
         {
            for(const idSegmentation tmp :  segmentation_[v]){
               if ((-tmp - 1) == (idSegmentation)id) {
                  return true;
               }
            }
            return false;
         }

         idNode getNodeId(const idVertex v) const
         {
            for (const idSegmentation tmp : segmentation_[v]) {
               if (tmp < 0) {
                  return -tmp - 1;
               }
            }
            return nullNode;
         }

         idSuperArc getFirstArcId(const idVertex v) const
         {
            for (const idSegmentation tmp : segmentation_[v]) {
               if (tmp >= 0) {
                  return tmp;
               }
            }
            return nullSuperArc;
         }

#ifndef TTK_ENABLE_KAMIKAZE
         valence getNbVisit(const idVertex v) const
         {
            // Debug only, costly
            valence s = 0;
            for(const idSuperArc tmp :  segmentation_[v]){
               std::ignore = tmp;
               ++s;
            }
            return s;
         }
#endif

         // direct access for openmp capture
         const valence& valUp(const idVertex v) const
         {
            return valUp_[v];
         }

         valence& valUp(const idVertex v)
         {
            return valUp_[v];
         }

         const valence& valDown(const idVertex v) const
         {
            return valDown_[v];
         }

         valence& valDown(const idVertex v)
         {
            return valDown_[v];
         }

         // Build structure
         // ---------------

         void addLeaf(const idVertex v, bool isMax)
         {
            leaves_.emplace_back(std::make_tuple(v, isMax));
         }

         idNode makeNode(const idVertex v)
         {
            if (isNode(v)) {
               return getNodeId(v);
            }

            const idNode newNode = nodes_.getNext();
            nodes_[newNode].setVerterIdentifier(v);
            visit(v, newNode, false);
            return newNode;
         }

         idSuperArc openArc(const idNode downId, Propagation * p = nullptr)
         {
             idSuperArc newArc = arcs_.getNext();
             arcs_[newArc].setDownNodeId(downId);
             if (p) {
                arcs_[newArc].setPropagation(p);
             }
             return newArc;
         }

         void closeArc(const idSuperArc arc, const idNode upId)
         {
            arcs_[arc].setUpNodeId(upId);
         }

         void makeArc(const idNode downId, const idNode upId)
         {
            arcs_.emplace_back(SuperArc{downId, upId});
         }

         idSuperArc makeHiddenArc(Propagation* const lp)
         {
            idSuperArc newArc = arcs_.getNext();
            arcs_[newArc].hide();
            arcs_[newArc].setPropagation(lp);
            return newArc;
         }

         // Process
         // -------

         // sort leaves vector by scalar value,
         // can be done in parallel
         template<typename ScalarType>
         void sortLeaves(const Scalars<ScalarType>* s, const bool parallel = false){
            auto compare_fun = [&](std::tuple<idVertex, bool> a, std::tuple<idVertex, bool> b) {
               return s->isLower(std::get<0>(a), std::get<0>(b));
            };
            if(parallel){
               ::ttk::ftr::parallel_sort<decltype(leaves_.begin()), std::tuple<idVertex, bool>>(
                   leaves_.begin(), leaves_.end(), compare_fun);
            } else {
               ::ttk::ftr::sort<decltype(leaves_.begin()), std::tuple<idVertex, bool>>(
                   leaves_.begin(), leaves_.end(), compare_fun);
            }
         }

         template<typename ScalarType>
         void shuffleLeaves(const Scalars<ScalarType>* s)
         {
            std::random_shuffle(leaves_.begin(), leaves_.end());
         }

         // Merge porpagations under a join saddle point inside the given prop
         // Only concern propagation growing in the same direction
         void mergeAtSaddle(const idNode saddleId, Propagation* const localProp);

         // Link nodes to arcs when arc are completely created
         void arcs2nodes(VertCompFN comp);

         // Tools
         // -----

         std::string print(const int verbosity) const;

         std::string printArc(const idSuperArc arcId) const;

         std::string printNode(const idNode nodeId) const;

         std::string printVisit(const idVertex v) const;

         // DEBUG function
         std::string printVisit() const;

         // Initialize functions
         // --------------------

         void alloc() override;

         void init() override;

      };

   }
}

