/// \ingroup base
/// \class ttk::ftr::DynamicGraph
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-15
///
/// \brief TTK %fTRGraph dynamic graph tracking the evolution of level sets
///
/// This class deal with dynamic graph part of the algorithm, thracking the number of contour
/// on each vertex to deduce the Reeb graph.
/// This is done using an ST-tree.
///
/// \sa ttk::FTRGraph

#ifndef DYNAMICGRAPH_H
#define DYNAMICGRAPH_H

#include "FTRCommon.h"
#include "AtomicVector.h"

namespace ttk
{
   namespace ftr
   {
      template<typename ScalarType>
      struct DynGraphNode;

      using idRoot = int;

      template <typename ScalarType>
      class DynamicGraph : public Allocable
      {
        private:
         AtomicVector<DynGraphNode<ScalarType>> nodes_;
         idNode nbTrees_;

        public:
         DynamicGraph();
         virtual ~DynamicGraph();

         // Initialize functions
         // --------------------

         /// preallocate the array of tree
         /// representing the forset,
         /// \pre needs maxTreeForest_ to be set.
         void alloc() override;

         void init() override;

         // Dyn Graph functions
         // -------------------

         DynGraphNode<ScalarType>* find(const DynGraphNode<ScalarType>* const n) const
         {
            return n->findRoot();
         }

         idNode getNumberOfTree(void) const
         {
            return nbTrees_;
         }

         DynGraphNode<ScalarType>* addTree()
         {
            ++nbTrees_;
            const auto idNext = nodes_.getNext();
            return &nodes_[idNext];
         }

         void insertEdge(DynGraphNode<ScalarType>* const n1, DynGraphNode<ScalarType>* const n2,
                         const ScalarType w)
         {
            bool mergeTwoTree = n1->insertEdge(n2, w);
            if (mergeTwoTree) {
               --nbTrees_;
            }
         }

         void removeEdge(DynGraphNode<ScalarType>* const n)
         {
            n->removeEdge();
            ++nbTrees_;
         }

         // Debug

         void print(void);

         void test(void);

      };

      template<typename ScalarType>
      struct DynGraphNode {
         DynGraphNode* parent_;
         ScalarType    weight_;
         idNode        nbChilds_;
         idRoot        rootId_;

         explicit DynGraphNode() : parent_(nullptr), weight_(0), nbChilds_(0), rootId_(-1)
         {
         }

         DynGraphNode(const DynGraphNode& other)
         {
            operator=(other);
         }

         DynGraphNode& operator=(const DynGraphNode& other)
         {
            parent_   = other.parent_;
            weight_   = other.weight_;
            nbChilds_ = other.nbChilds_;
            rootId_   = other.rootId_;
            return *this;
         }

         // Accessor functions
         // ------------------

         idVertex getRootIndex(void) const
         {
            return rootId_;
         }

         void setRootIndex(const idVertex v)
         {
            rootId_ = v;
         }

         // Graph functions
         // ---------------

         /// Make this node the root of its tree
         // Various way to do that, test perfs ?
         void evert(void);

         /// Get representative node
         DynGraphNode* findRoot(void) const;

         /// Get representative node and kepp track of
         /// the node with the min weight on the path
         /// \ret tuple<root, node with min weight>
         std::tuple<DynGraphNode<ScalarType>*, DynGraphNode<ScalarType>*> findMinWeightRoot() const;

         /// Create a new edge between this node and the node n
         /// return true if we have merged two tree, false if it was just an intern operation
         bool insertEdge(DynGraphNode* const n, const ScalarType weight);

         void removeEdge(void);

      };
   }
}

#include "DynamicGraph_Template.h"

#endif /* end of include guard: DYNAMICGRAPH_H */
