/// \ingroup base
/// \class ttk::DynamicTree
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \date 2019-09-15
///
/// \brief Implements the Dynamic Tree data-structure (or ST-Tree)
///
/// Adapted by Jules Vidal from the ttk::ftr::DynamicGraph class from
/// FTRGraph.h
///
/// \sa ttk::FTRGraph

#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace ttk {

  /// \brief class representing a node
  /// of a tree and the link to its parent if not the root
  struct DynTreeNode {
    DynTreeNode *parent_{nullptr};

    // Accessor functions
    // ------------------

    inline bool hasParent(void) const {
      return parent_ != nullptr;
    }

    // Graph functions
    // ---------------

    /// Make this node the root of its tree
    // Various way to do that, test perfs ?
    void evert(void);

    /// Get representative node
    DynTreeNode *findRoot(void) const;

    /// Create a new edge between this node and the node n
    /// @return true if we have merged two tree, false if it was just an intern
    /// operation
    bool insertEdge(DynTreeNode *const n);

    /// Remove the link between this node and its parent, thus makeing a new
    /// root
    inline void removeEdge(void) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!parent_) {
        std::cerr << "[DynamicTree]: DynTree remove edge in root node"
                  << std::endl;
        return;
      }
#endif
      parent_ = nullptr;
    }
  };

  // WARNING : this class allows to maintain the number of components
  // as it evolves, IT IS NOT PARALLEL SAFE
  class DynamicTree {
    std::vector<DynTreeNode> nodes_{};

  public:
    // Initialize functions
    // --------------------

    inline void alloc(const std::size_t nbNodes) {
      nodes_.resize(nbNodes);
    }

    /// \brief get the node with the id: nid
    inline DynTreeNode *getNode(const std::size_t nid) {
      return &nodes_[nid];
    }

    inline const DynTreeNode *getNode(const std::size_t nid) const {
      return &nodes_[nid];
    }

    /// \brief get the id of the node: node
    inline std::size_t getNodeId(DynTreeNode *node) {
      return node - &nodes_[0];
    }

    /// @return true if we have merged two tree, false if it was just an intern
    /// operation
    bool insertEdge(DynTreeNode *const n1, DynTreeNode *const n2) {
      return n1->insertEdge(n2);
    }

    /// inert or replace existing edge between n1 and n2
    inline bool insertEdge(const std::size_t n1, const std::size_t n2) {
      return getNode(n1)->insertEdge(getNode(n2));
    }

    /// remove the link btwn n and its parent
    inline void removeEdge(DynTreeNode *const n) {
      n->removeEdge();
    }

    /// remove the link btwn n and its parent
    inline void removeEdge(const std::size_t nid) {
      getNode(nid)->removeEdge();
    }

    /// remove the edge btwn n1 and n2
    /// @return 0 if not an edge
    int removeEdge(DynTreeNode *const n1, DynTreeNode *const n2) {
      if(n1->parent_ == n2) {
        n1->removeEdge();
        return 1;
      }
      if(n2->parent_ == n1) {
        n2->removeEdge();
        return 2;
      }
      return 0;
    }

    /// remove the edge btwn n1 and n2
    /// @return 0 if not an edge
    inline int removeEdge(const std::size_t nid1, const std::size_t nid2) {
      return removeEdge(getNode(nid1), getNode(nid2));
    }

    inline std::size_t getCCFromNode(const std::size_t n0) const {
      return getNode(n0)->findRoot() - &nodes_[0];
    }
    void retrieveNbCC(std::vector<size_t> &nbccIds) const {
      for(size_t nId = 0; nId < nodes_.size(); nId++) {
        if(nodes_[nId].parent_ == nullptr)
          nbccIds.emplace_back(nId);
      }
    }
    inline size_t getNbCC() const {
      return std::count_if(
        nodes_.begin(), nodes_.end(),
        [](const DynTreeNode &dtn) { return dtn.parent_ == nullptr; });
    }

    // Debug

    std::string print(void) const;
    std::string printNbCC(void) const;
  };

} // namespace ttk
