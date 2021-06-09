#pragma once

#include <boost/heap/fibonacci_heap.hpp>
#include <vector>

namespace ttk {
  namespace lts {

    /// Superlevel Set Component Propagation
    template <typename IT>
    struct Propagation {

      // union find parent (path compressed)
      Propagation<IT> *ufParent{this};

      // direct parent propagation
      Propagation<IT> *parent{this};

      // propagation data
      std::vector<IT> criticalPoints;
      boost::heap::fibonacci_heap<std::pair<IT, IT>> queue;
      IT segmentSize{0};
      std::vector<IT> segment;
      bool aborted{false};

      mutable IT nIterations{0};

      // find parent with path compression
      inline Propagation *find() {
        if(this->ufParent == this)
          return this;
        else {
          auto tmp = this->ufParent->find();
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif // TTK_ENABLE_OPENMP
          this->ufParent = tmp;
          return this->ufParent;
        }
      }

      static inline int unify(Propagation<IT> *p0, Propagation<IT> *p1) {
        if(p1 == nullptr)
          return 1;

        p0 = p0->find();
        p1 = p1->find();

        if(p0 == p1)
          return 1;

        // store parent without path compression
        p1->parent = p0;

        // update union find tree
        p1->ufParent = p0;

        p0->segmentSize += p1->segmentSize;

        // merge f. heaps
        p0->queue.merge(p1->queue);

        return 1;
      }
    };
  } // namespace lts
} // namespace ttk