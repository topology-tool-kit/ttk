#pragma once

#include "DynamicGraph.h"

#include <sstream>

namespace ttk {
  namespace ftr {

    // DynamicGraph ----------------------------------

    template <typename Type>
    DynamicGraph<Type>::DynamicGraph() {
    }

    template <typename Type>
    DynamicGraph<Type>::~DynamicGraph() {
    }

    template <typename Type>
    void DynamicGraph<Type>::alloc() {
      nodes_.resize(nbElmt_);
    }

    template <typename Type>
    void DynamicGraph<Type>::init() {
    }

    template <typename Type>
    int DynamicGraph<Type>::removeEdge(DynGraphNode<Type> *const n1,
                                       DynGraphNode<Type> *const n2) {
      if(n1->parent_ == n2) {
        removeEdge(n1);
        return 1;
      }

      if(n2->parent_ == n1) {
        removeEdge(n2);
        return 2;
      }

      return 0;
    }

    template <typename Type>
    std::string DynamicGraph<Type>::print(void) {

      std::stringstream res;

      for(const auto &node : nodes_) {
        if(1 or node.parent_) {
          res << "id: " << &node - &nodes_[0];
          if(node.parent_) {
            res << ", parent: " << node.parent_ - &nodes_[0];
          } else {
            res << ", parent: X";
          }
          res << " root: " << findRoot(&node) - &nodes_[0];
          res << " weight: " << (float)node.weight_;
          res << " cArc: " << node.corArc_;
          res << std::endl;
        }
      }
      return res.str();
    }

    template <typename Type>
    std::string DynamicGraph<Type>::print(
      const std::function<std::string(std::size_t)> &printFunction) {

      std::stringstream res;

      for(const auto &node : nodes_) {
        if(node.parent_) {
          res << "id: " << printFunction(&node - &nodes_[0])
              << " weight: " << (float)node.weight_;
          if(node.parent_) {
            res << ", parent: " << printFunction(node.parent_ - &nodes_[0]);
          } else {
            res << ", parent: X";
          }
          res << " root: " << printFunction(findRoot(&node) - &nodes_[0]);
        }
      }
      return res.str();
    }

    template <typename Type>
    std::string DynamicGraph<Type>::printNbCC(void) {

      std::stringstream res;
      std::vector<DynGraphNode<Type> *> roots;
      roots.reserve(nodes_.size());
      for(const auto &n : nodes_) {
        roots.emplace_back(n.findRoot());
      }
      std::sort(roots.begin(), roots.end());
      const auto it = std::unique(roots.begin(), roots.end());
      roots.erase(it, roots.end());
      res << "nb nodes " << nodes_.size() << std::endl;
      res << "nb cc " << roots.size() << std::endl;
      return res.str();
    }

    // DynGraphNode ----------------------------------

    template <typename Type>
    void DynGraphNode<Type>::evert(void) {
      if(!parent_)
        return;

      DynGraphNode<Type> *curNode = this;

      DynGraphNode<Type> *parentNode = curNode->parent_;
      Type parentWeight = curNode->weight_;

      DynGraphNode<Type> *gParentNode = parentNode->parent_;
      Type gParentWeight = parentNode->weight_;

      curNode->parent_ = nullptr;

      // Reverse all the node until the root
      while(true) {
        parentNode->parent_ = curNode;
        parentNode->weight_ = parentWeight;

        curNode = parentNode;
        parentNode = gParentNode;
        parentWeight = gParentWeight;

        if(gParentNode) {
          gParentWeight = gParentNode->weight_;
          gParentNode = gParentNode->parent_;
        } else {
          // keep same arc than the current root
          // if cur > this ?
          corArc_ = curNode->corArc_;
          break;
        }
      }
    }

    template <typename Type>
    DynGraphNode<Type> *DynGraphNode<Type>::findRoot(void) const {
      // the lastNode trick is used so we are sure to have a non null
      // return even if another thread is touching these nodes.
      DynGraphNode *curNode = const_cast<DynGraphNode<Type> *>(this);
      DynGraphNode *lastNode = curNode;
      while(curNode) {
        lastNode = curNode;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read
#endif
        curNode = curNode->parent_;
      }
      return lastNode;
    }

    template <typename Type>
    idSuperArc DynGraphNode<Type>::findRootArc(void) const {
      const auto root = findRoot();
      return root != nullptr ? root->corArc_ : -1;
    }

    template <typename Type>
    void DynGraphNode<Type>::setRootArc(const idSuperArc arcId) {
      corArc_ = arcId;
    }

    template <typename Type>
    std::tuple<DynGraphNode<Type> *, DynGraphNode<Type> *>
      DynGraphNode<Type>::findMinWeightRoot() const {

      DynGraphNode *minNode = const_cast<DynGraphNode<Type> *>(this);
      auto minW = minNode->weight_;
      DynGraphNode *curNode = minNode->parent_;

      if(!curNode)
        return std::make_tuple(curNode, minNode);

      while(curNode->parent_) {
        if(curNode->weight_ < minW) {
          minNode = curNode;
          minW = minNode->weight_;
        }
        curNode = curNode->parent_;
      }
      return std::make_tuple(curNode, minNode);
    }

    template <typename Type>
    bool DynGraphNode<Type>::insertEdge(DynGraphNode<Type> *const n,
                                        const Type weight,
                                        const idSuperArc corArc) {
      evert();
      auto nNodes = n->findMinWeightRoot();

      if(std::get<0>(nNodes) != this) {
        // The two nodes are in two different trees
        parent_ = n;
        weight_ = weight;
        n->corArc_ = corArc;
        return true;
      }

      // here the nodes are in the same tree

      if(weight > std::get<1>(nNodes)->weight_) {
        // We need replace the min edge by the new one as the current weight is
        // higher

        // add arc (Parsa like)
        parent_ = n;
        weight_ = weight;
        // corArc_ = corArc;

        // remove old
        std::get<1>(nNodes)->parent_ = 0;
        std::get<1>(nNodes)->corArc_ = corArc;
      } else {
        corArc_ = corArc;
      }

      return false;
    }

    template <typename Type>
    void DynGraphNode<Type>::removeEdge(void) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!parent_) {
        Debug dbg{};
        dbg.setDebugMsgPrefix("DynamicGraph");
        dbg.printErr("DynGraph remove edge in root node");
        return;
      }
#endif

      parent_ = nullptr;
    }

  } // namespace ftr
} // namespace ttk
