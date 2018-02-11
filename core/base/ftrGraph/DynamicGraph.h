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

#include <set>
#include <vector>

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
         std::vector<DynGraphNode<ScalarType>> nodes_;

        public:
         DynamicGraph();
         virtual ~DynamicGraph();

         // Initialize functions
         // --------------------

         void setNumberOfNodes(const std::size_t nbNodes)
         {
            // we use the nbVerts_ as the number
            // of nodes to allocate
            setNumberOfVertices(nbNodes);
         }

         /// preallocate the array of tree
         /// representing the forset,
         /// \pre needs nbVerts_ to be set.
         void alloc() override;

         void init() override;

         // Dyn Graph functions
         // -------------------

         /// \brief get the node with the id: nid
         DynGraphNode<ScalarType>* getNode(const std::size_t nid)
         {
            return &nodes_[nid];
         }

         /// \brief get the id of the node: node
         std::size_t getNodeId(DynGraphNode<ScalarType>* node)
         {
            return std::distance(nodes_, node);
         }

         // check wether or not this node is connected to others
         bool isDisconnected(const DynGraphNode<ScalarType>* const node) const
         {
            return !node->hasParent();
         }

         // check wether or not this node is connected to others
         bool isDisconnected(const std::size_t nid)
         {
            return !getNode(nid)->hasParent();
         }

         /// \brief recover the root of a node
         DynGraphNode<ScalarType>* findRoot(const DynGraphNode<ScalarType>* const node)
         {
            return node->findRoot();
         }

         /// \brief recover the root of a node using its id
         DynGraphNode<ScalarType>* findRoot(const std::size_t nid)
         {
            return findRoot(getNode(nid));
         }

         /// \brief recover the root of several nodes once, using
         /// brace initializers style: findRoot({n1,n2})
         /// duplicate are removed by the std::set
         std::set<DynGraphNode<ScalarType>*> findRoot(
             std::initializer_list<DynGraphNode<ScalarType>*> nodes)
         {
            std::set<DynGraphNode<ScalarType>*> roots;
            for(auto* n : nodes) {
               roots.emplace(findRoot(n));
            }
            return roots;
         }

         /// \brief findRoot but using ids of the nodes
         std::set<DynGraphNode<ScalarType>*> findRoot(std::initializer_list<std::size_t> nodesIds)
         {
            std::set<DynGraphNode<ScalarType>*> roots;
            for(auto n : nodesIds) {
               roots.emplace(findRoot(n));
            }
            return roots;
         }

         /// \brief findRoot but using ids of the nodes in a vector
         template<typename type>
         std::set<DynGraphNode<ScalarType>*> findRoot(const std::vector<type>& nodesids)
         {
            std::set<DynGraphNode<ScalarType>*> roots;
            for(auto n : nodesids) {
               roots.emplace(findRoot(n));
            }
            return roots;
         }


         /// \brief count the number of connected component in a list of nodes
         std::size_t getNbCC(std::initializer_list<DynGraphNode<ScalarType>*> nodes)
         {
            return findRoot(nodes).size();
         }

         /// \brief count the number of connected component in a list of nodes using ids
         std::size_t getNbCC(std::initializer_list<std::size_t> nodesIds)
         {
            return findRoot(nodesIds).size();
         }


         /// \ret true if we have merged two tree, false if it was just an intern operation
         bool insertEdge(DynGraphNode<ScalarType>* const n1, DynGraphNode<ScalarType>* const n2,
                         const ScalarType w)
         {
           return n1->insertEdge(n2, w);
         }

         /// inert or replace existing edge between n1 and n2
         bool insertEdge(const std::size_t n1, const std::size_t n2, const ScalarType w)
         {
            return insertEdge(getNode(n1), getNode(n2), w);
         }

         /// remove the link btwn n and its parent
         void removeEdge(DynGraphNode<ScalarType>* const n)
         {
            n->removeEdge();
         }

         /// remove the link btwn n and its parent
         void removeEdge(const std::size_t nid)
         {
            removeEdge(getNode(nid));
         }

         /// remove the edge btwn n1 and n2
         /// \ret 0 if not an edge
         int removeEdge(DynGraphNode<ScalarType>* const n1, DynGraphNode<ScalarType>* const n2);

         /// remove the edge btwn n1 and n2
         /// \ret 0 if not an edge
         int removeEdge(const std::size_t nid1, const std::size_t nid2)
         {
            return removeEdge(getNode(nid1), getNode(nid2));
         }

         // Debug

         void print(void);

         void print(std::function<std::string(std::size_t)>);

         void test(void);
      };

      /// \brief class representing a node
      /// of a tree and the link to its parent if not the root
      template<typename ScalarType>
      struct DynGraphNode {
         DynGraphNode* parent_;
         ScalarType    weight_;
         idNode        nbChilds_;

         explicit DynGraphNode() : parent_(nullptr), weight_(0), nbChilds_(0)
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
            return *this;
         }

         // Accessor functions
         // ------------------

         ScalarType getWeight(void) const
         {
            return weight_;
         }

         idNode getNbChilds(void) const
         {
            return nbChilds_;
         }

         bool hasParent(void) const
         {
            return parent_;
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
         /// \ret true if we have merged two tree, false if it was just an intern operation
         bool insertEdge(DynGraphNode* const n, const ScalarType weight);

         /// Remove the link between this node and its parent, thus makeing a new root
         void removeEdge(void);

      };
   }
}

#include "DynamicGraph_Template.h"

#endif /* end of include guard: DYNAMICGRAPH_H */
