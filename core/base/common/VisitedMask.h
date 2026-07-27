/// \ingroup base
/// \class ttk::VisitedMask
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date November 2021.

#pragma once

#include <DataTypes.h>

#include <algorithm>
#include <vector>

namespace ttk {

  /**
   * @brief Auto-cleaning re-usable graph propagations data structure
   *
   * This structure stores two vector references. The first vector
   * maps the nodes of a graph to a boolean (true if the node is
   * visited, false if not). The second reference points to a vector
   * that contains the visited nodes indices. Using RAII, the boolean
   * vector is cleaned using the indices vector when this structure
   * goes out of scope. See @ref
   * ttk::MorseSmaleComplex::getSaddleConnectors for an example use.
   *
   */
  struct VisitedMask {
    std::vector<bool> &isVisited_;
    std::vector<SimplexId> &visitedIds_;

    VisitedMask(std::vector<bool> &isVisited,
                std::vector<SimplexId> &visitedIds)
      : isVisited_{isVisited}, visitedIds_{visitedIds} {
    }

    ~VisitedMask() {
      // use RAII to clean & reset referenced vectors
      for(const auto id : this->visitedIds_) {
        this->isVisited_[id] = false;
      }
      // set size to 0 but keep allocated memory
      this->visitedIds_.clear();
    }

    inline void insert(const SimplexId id) {
      this->isVisited_[id] = true;
      this->visitedIds_.emplace_back(id);
    }

    inline bool remove(const SimplexId id) {
      if(!this->isVisited_[id]) {
        return false;
      }
      const auto it
        = std::find(this->visitedIds_.begin(), this->visitedIds_.end(), id);
      if(it != this->visitedIds_.end()) {
        this->visitedIds_.erase(it);
      }
      this->isVisited_[id] = false;
      return true;
    }
  };

} // namespace ttk
