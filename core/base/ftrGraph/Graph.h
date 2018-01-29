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
#include "DataTypes.h"
#include "FTRCommon.h"
#include "Node.h"
#include "SuperArc.h"

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
         AtomicVector<idVertex> leaves_;
         AtomicVector<Node>     nodes_;
         AtomicVector<SuperArc> arcs_;
         std::vector<graphElmt> vert2Graph_;

        public:
         Graph();
         virtual ~Graph();

         // Accessor on structure
         // ---------------------

         idNode getNumberOfNodes(void) const
         {
            return nodes_.size();
         }

         idNode getNumberOfLeaves(void) const
         {
            return leaves_.size();
         }

         const Node& getNode(const idNode id) const
         {
           return nodes_[id];
         }

         const SuperArc& getArc(const idSuperArc id) const
         {
            return arcs_[id];
         }

         bool isVisied(const idVertex v) const
         {
            return vert2Graph_[v] == nullGraphElmt;
         }

         bool isNode(const idVertex v) const
         {
#ifndef TTK_ENABLE_KAMIKAZE
            if (!isVisied(v)) {
               std::cerr << "[FTR Graph]: Node is not visited: " << v << std::endl;
            }
#endif
            return vert2Graph_[v] < 0;
         }

         bool isArc(const idVertex v) const
         {
#ifndef TTK_ENABLE_KAMIKAZE
            if (!isVisied(v)) {
               std::cerr << "[FTR Graph]: Node is not visited: " << v << std::endl;
            }
#endif
            return vert2Graph_[v] >= 0;
         }

         const Node& getCorrespondingNode(const idVertex v) const
         {
#ifndef TTK_ENABLE_KAMIKAZE
            if (!isNode(v)) {
               std::cerr << "[FTR Graph]: " << v << " is not a node" << std::endl;
            }
#endif
            return nodes_[-vert2Graph_[v]];
         }

         const SuperArc& getCorrespondingArc(const idVertex v) const
         {
#ifndef TTK_ENABLE_KAMIKAZE
            if (!isArc(v)) {
               std::cerr << "[FTR Graph]: " << v << " is not an arc" << std::endl;
            }
#endif
            return arcs_[vert2Graph_[v]];
         }

         // Build structure
         // ---------------

         void addLeaf(const idVertex v)
         {
            leaves_.emplace_back(v);
         }

         void makeNode(const idVertex v)
         {
            nodes_.emplace_back(Node{v});
         }

         void makeArc(const idNode downId, const idNode upId)
         {
            arcs_.emplace_back(SuperArc{downId, upId});
         }

         // Initialize functions
         // --------------------

         void alloc() override;

         void init() override;

      };

   }
}

