#pragma once

#include <boost/heap/fibonacci_heap.hpp>
#include <vector>

namespace ttk {

  /// Superlevel Set Component Propagation
  template <typename IT>
  struct SSCPropagation {

    // union find members
    SSCPropagation<IT> *parent{this};

    // propagation data
    std::vector<IT> criticalPoints;
    boost::heap::fibonacci_heap<std::pair<IT, IT>> queue;
    IT segmentSize{0};
    std::vector<IT> segment;
    IT lastEncounteredCriticalPoint{-1};

    mutable IT nIterations{0};

    inline SSCPropagation *find() {
      if(this->parent == this)
        return this;
      else {
        auto tmp = this->parent->find();
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif // TTK_ENABLE_OPENMP
        this->parent = tmp;
        return this->parent;
      }
    }

    static inline SSCPropagation<IT> *unify(SSCPropagation<IT> *p0,
                                            SSCPropagation<IT> *p1,
                                            const IT *orderField) {
      SSCPropagation<IT> *parent = p0;
      SSCPropagation<IT> *child = p1;

      // determine parent and child based on persistence
      if(orderField[parent->criticalPoints[0]]
         < orderField[child->criticalPoints[0]]) {
        SSCPropagation<IT> *temp = parent;
        parent = child;
        child = temp;
      }

      // update union find tree
      // #pragma omp atomic write
      child->parent = parent;

      parent->segmentSize += child->segmentSize;

      // merge f. heaps
      parent->queue.merge(child->queue);

      return parent;
    }
  };
} // namespace ttk
