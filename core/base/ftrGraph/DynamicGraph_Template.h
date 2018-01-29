#ifndef DYNAMICGRAPH_TEMPLATE_H
#define DYNAMICGRAPH_TEMPLATE_H

#include "DynamicGraph.h"

namespace ttk
{
   namespace ftr
   {

      // DynamicGraph ----------------------------------

      template <typename ScalarType>
      DynamicGraph<ScalarType>::DynamicGraph(): nbTrees_(0)
      {
      }

      template <typename ScalarType>
      DynamicGraph<ScalarType>::~DynamicGraph()
      {
      }

      template <typename ScalarType>
      void DynamicGraph<ScalarType>::alloc()
      {
         nodes_.reserve(nbVerts_);
         nbTrees_ = 0;
      }

      template<typename ScalarType>
      void DynamicGraph<ScalarType>::init()
      {
      }

      template <typename ScalarType>
      void DynamicGraph<ScalarType>::print(void)
      {
         using namespace std;

         cout << "nb component : " << nbTrees_ << endl;
         for(const auto& node : nodes_) {
            cout << &node << " weight : " << node.weight_ << " parent : " << node.parent_;
            cout << " root : " << find(&node) << endl;
         }
      }

      template <typename ScalarType>
      void DynamicGraph<ScalarType>::test(void)
      {
         DynamicGraph testGraph;
         testGraph.setNumberOfVertices(4);
         testGraph.alloc();
         testGraph.init();

         DynGraphNode<ScalarType>* n1 = testGraph.addTree();
         DynGraphNode<ScalarType>* n2 = testGraph.addTree();

         if(testGraph.getNumberOfTree() != 2){
            std::cerr << "[Dynamic Graph]: bad initial number of CC in test" << std::endl;
         }

         testGraph.insertEdge(n1, n2, 0);

         // 1 -> 2

         if(testGraph.getNumberOfTree() != 1){
            std::cerr << "[Dynamic Graph]: bad number of CC after insert edge" << std::endl;
         }

         DynGraphNode<ScalarType>* n3 = testGraph.addTree();
         DynGraphNode<ScalarType>* n4 = testGraph.addTree();

         testGraph.insertEdge(n3, n4, 0);

         // 1 -> 2
         // 3 -> 4

         if(testGraph.getNumberOfTree() != 2){
            std::cerr << "[Dynamic Graph]: bad number of CC after insert edge 2" << std::endl;
         }

         testGraph.insertEdge(n2, n3, 0);

         // 1 -> 2 -> 3 -> 4

         if(testGraph.getNumberOfTree() != 1){
            std::cerr << "[Dynamic Graph]: bad number of CC after insert edge 3" << std::endl;
         }

         testGraph.removeEdge(n1);

         // 1
         // 2 -> 3 -> 4

         if(testGraph.getNumberOfTree() != 2){
            std::cerr << "[Dynamic Graph]: bad number of CC after remove edge" << std::endl;
         }

         // replace existing
         testGraph.insertEdge(n2, n4, 1);

         // 1
         // 3 -> 2 -> 4

         if(testGraph.getNumberOfTree() != 2){
            std::cerr << "[Dynamic Graph]: bad number of CC after replace edge" << std::endl;
         }

         testGraph.insertEdge(n1, n4, 0);

         // 3 -> 2 -> 4
         //      1 _/

         if(testGraph.getNumberOfTree() != 1){
            std::cerr << "[Dynamic Graph]: bad number of CC after insert edge 4" << std::endl;
         }

         testGraph.print();
      }

      // DynGraphNode ----------------------------------

      template <typename ScalarType>
      DynGraphNode<ScalarType>* DynGraphNode<ScalarType>::findRoot(void) const
      {
         DynGraphNode* curNode = const_cast<DynGraphNode<ScalarType>*>(this);
         while (curNode->parent_ != nullptr) {
            curNode = curNode->parent_;
         }
         return curNode;
      }

      template <typename ScalarType>
      std::tuple<DynGraphNode<ScalarType>*, DynGraphNode<ScalarType>*>
      DynGraphNode<ScalarType>::findMinWeightRoot() const
      {
         DynGraphNode* curNode = const_cast<DynGraphNode<ScalarType>*>(this);
         DynGraphNode* minNode = const_cast<DynGraphNode<ScalarType>*>(this);

         while (curNode->parent_) {
            if (curNode->weight_ < minNode->weight_) {
               minNode = curNode;
            }
            curNode = curNode->parent_;
         }
         return std::make_tuple(curNode, minNode);
      }

      template<typename ScalarType>
      bool DynGraphNode<ScalarType>::insertEdge(DynGraphNode<ScalarType>* const n, const ScalarType weight)
      {
         evert();
         auto nNodes = n->findMinWeightRoot();

         if (std::get<0>(nNodes) != this) {
            // The two nodes are in two different trees
            parent_ = n;
            weight_ = weight;
            n->nbChilds_++;
            return true;
         }

         // here the nodes are in the same tree

         if (weight > std::get<1>(nNodes)->weight_) {
            // We need replace this edge by the new one as the current weight is higher

            // add arc
            parent_ = n;
            weight_ = weight;
            n->nbChilds_++;

            // remove old
            std::get<1>(nNodes)->parent_->nbChilds_--;
            std::get<1>(nNodes)->parent_ = 0;
         }

         return false;
      }

      template <typename ScalarType>
      void DynGraphNode<ScalarType>::removeEdge(void)
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (!parent_) {
            std::cerr << "[FTR Graph]: DynGraph remove edge in root node" << std::endl;
         }
#endif

         parent_->nbChilds_--;
         parent_ = nullptr;
      }

      template <typename ScalarType>
      void DynGraphNode<ScalarType>::evert(void)
      {
         if(!parent_) return;

         nbChilds_++;

         DynGraphNode<ScalarType>* curNode      = this;
         ScalarType                curWeight    = curNode->weight_;
         DynGraphNode<ScalarType>* parentNode   = curNode->parent_;
         DynGraphNode<ScalarType>* gParentNode  = parentNode->parent_;

         // becom a new root
         curNode->parent_ = nullptr;

         // Reverse all the node until the root
         while (true) {
            parentNode->parent_ = curNode;
            parentNode->weight_ = curWeight;
            if (gParentNode) {
               curNode     = parentNode;
               curWeight   = curNode->weight_;
               parentNode  = gParentNode;
               gParentNode = gParentNode->parent_;
            } else {
               parentNode->nbChilds_--;
               break;
            }
         }
      }
   }
}

#endif /* end of include guard: DYNAMICGRAPH_TEMPLATE_H */
