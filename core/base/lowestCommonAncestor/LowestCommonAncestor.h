/// \ingroup base
/// \class ttk::LowestCommonAncestor
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief Class to answer the lowest common ancestor requests of pairs of nodes
/// in a tree in constant time after a linear time preprocess.

#pragma once

// ttk includes
#include <Debug.h>
#include <RangeMinimumQuery.h>
// STL includes
#include <array>
#include <vector>

namespace ttk {

  class LowestCommonAncestor : virtual public Debug {
  public:
    class Node {
    public:
      inline void setAncestor(const int &id) {
        ancestor_ = id;
      }
      inline void addSuccessor(const int &id) {
        successor_.push_back(id);
      }
      inline int getNumberOfSuccessors() const {
        return static_cast<int>(successor_.size());
      }
      inline int getAncestorId() const {
        return ancestor_;
      }
      inline int getSuccessorId(const int &id) const {
        if(id < static_cast<int>(successor_.size())) {
          return successor_[id];
        } else {
          return -1;
        }
      }

    protected:
      int ancestor_{};
      std::vector<int> successor_{};
    };

    LowestCommonAncestor();

    /// Add a node in the tree
    /// \return Returns the id of the new node
    inline int addNode() {
      node_.push_back(Node());
      int id = static_cast<int>(node_.size());
      node_[id].setAncestor(id);
      return id;
    }

    inline void addNodes(const unsigned int &number) {
      node_.resize(node_.size() + number);
      for(unsigned int i = node_.size() - number; i < node_.size(); i++) {
        node_[i].setAncestor(i);
      }
    }

    /// Get the number of nodes in the tree
    inline unsigned int getNumberOfNodes() const {
      return node_.size();
    }

    /// \returns Returns a pointer to the id-th node
    inline Node &getNode(const unsigned int &id) {
      return node_[id];
    }

    /// Preprocess the tree structure to answer the query() calls in constant
    /// time. The preprocess takes linear time in number of nodes in the tree.
    int preprocess();

    /// Get the id of the lowest common ancestor of i and j.
    /// \pre preprocess() must have been called after the last change in the
    /// tree.
    inline int query(int i, int j) const {
      if(nodeFirstAppearence_[i] > nodeFirstAppearence_[j]) {
        std::swap(i, j);
      }
      return nodeOrder_[RMQuery(
        nodeFirstAppearence_[i], nodeFirstAppearence_[j])];
    }

  protected:
    int computeBlocs();
    int eulerianTransverse();
    int RMQuery(const int &i, const int &j) const;
    inline unsigned int min_pos_3(const std::array<int, 3> &triplet) const {
      if(triplet[0] < triplet[1]) {
        if(triplet[0] < triplet[2]) {
          return 0;
        } else {
          return 2;
        }
      } else {
        if(triplet[1] < triplet[2]) {
          return 1;
        } else {
          return 2;
        }
      }
    }

  protected:
    /* Tree structure */
    std::vector<Node> node_{};

    /* Eulerian Transverse */
    std::vector<int> nodeOrder_{};
    std::vector<int> nodeDepth_{};
    std::vector<int> nodeFirstAppearence_{};

    /* Range Minimum Query */
    int blocSize_{};
    // Boundaries of blocs
    std::vector<std::pair<int, int>> blocPartition_{};
    // Min values
    std::vector<int> blocMinimumValue_{};
    // RMQ of the blocMinimumValue_ vector
    RangeMinimumQuery<int> blocMinimumValueRMQ_{};
    // Positions of min values
    std::vector<int> blocMinimumPosition_{};
    // All queries for each possible bloc (positions)
    std::vector<std::vector<std::vector<int>>> normalizedBlocTable_{};
    // Corresponding normalized bloc for each bloc of nodeDepth_
    std::vector<int> blocToNormalizedBloc_{};
  };

} // namespace ttk
