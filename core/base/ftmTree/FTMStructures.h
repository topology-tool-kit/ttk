/// \ingroup base
//
/// \class ttk::FTMTree_MT
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date September 2016.
///
///\brief TTK structures for the contour tree

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <forward_list>
#include <iterator>
#include <memory>
#include <vector>

#include <boost/heap/fibonacci_heap.hpp>

#include "FTMAtomicVector.h"
#include "FTMDataTypes.h"

// todo remove
#include <iostream>

namespace ttk {
  namespace ftm {
    // Compute parameters (global)
    struct Params {
      TreeType treeType;
      bool segm = true;
      bool normalize = true;
      bool advStats = true;
      int samplingLvl = 0;
    };

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
    struct ActiveTask {
      float begin = -1;
      float end = -1;
      SimplexId origin = nullVertex;
    };
#endif

    // Scalar related containers (global)
    struct Scalars {
      SimplexId size;
      void *values;
      void *offsets;

      std::shared_ptr<std::vector<SimplexId>> sortedVertices, mirrorVertices;

      // Need vertices to be sorted : use mirrorVertices.

      bool isLower(SimplexId a, SimplexId b) const {
        return (*mirrorVertices)[a] < (*mirrorVertices)[b];
      }
      bool isEqLower(SimplexId a, SimplexId b) const {
        return (*mirrorVertices)[a] <= (*mirrorVertices)[b];
      }

      bool isHigher(SimplexId a, SimplexId b) const {
        return (*mirrorVertices)[a] > (*mirrorVertices)[b];
      }
      bool isEqHigher(SimplexId a, SimplexId b) const {
        return (*mirrorVertices)[a] >= (*mirrorVertices)[b];
      }

      Scalars()
        : size(0), values(nullptr), offsets(nullptr), sortedVertices(nullptr),
          mirrorVertices(nullptr) {
      }

      // Heavy
      Scalars(const Scalars &o)
        : size(o.size), values(o.values), offsets(o.offsets),
          sortedVertices(o.sortedVertices), mirrorVertices(o.mirrorVertices) {
        std::cout << "copy in depth, bad perfs" << std::endl;
      }

      // Sort
      template <typename type>
      void qsort(type arr[],
                 const long int begin,
                 const long int stop,
                 std::function<bool(type, type)> comp) const {
        if(begin >= stop)
          return;

#ifdef TTK_ENABLE_OPENMP
        static const long int MINSIZE = 10;
#endif

        long int left = begin - 1;
        long int right = stop + 1;
        const type pivot = arr[begin];

        while(1) {
          while(comp(pivot, arr[--right]))
            ;
          while(++left <= stop && !comp(pivot, arr[left]))
            ;

          if(left < right)
            swap_el<type>(arr, left, right);
          else
            break;
        }

        swap_el<type>(arr, begin, right);
#ifdef TTK_ENABLE_OPENMP
#pragma omp task untied if(right - begin > MINSIZE)
#endif
        qsort(arr, begin, right - 1, comp);
#ifdef TTK_ENABLE_OPENMP
#pragma omp task untied if(stop - right > MINSIZE)
#endif
        qsort(arr, right + 1, stop, comp);
      }

    private:
      template <typename type>
      static void swap_el(type arr[], const size_t a, const size_t b) {
        const type tmp = arr[a];
        arr[a] = arr[b];
        arr[b] = tmp;
      }
    };

    struct CurrentState {
      SimplexId vertex;
      boost::heap::fibonacci_heap<SimplexId, boost::heap::compare<VertCompFN>>
        propagation;

      CurrentState(SimplexId startVert, VertCompFN vertComp)
        : vertex(startVert), propagation(vertComp) {
      }

      CurrentState(VertCompFN vertComp)
        : vertex(nullVertex), propagation(vertComp) {
        // will need to use setStartVert before use
      }

      void setStartVert(const SimplexId v) {
        vertex = v;
      }

      SimplexId getNextMinVertex(void) {
        vertex = propagation.top();
        propagation.pop();
        return vertex;
      }

      void addNewVertex(const SimplexId v) {
        propagation.emplace(v);
      }

      void merge(CurrentState &other) {
        propagation.merge(other.propagation);
        vertex = propagation.top();
      }

      bool empty() {
        return propagation.empty();
      }

      // DEBUG ONLY
      bool find(SimplexId v) {
        return std::find(propagation.begin(), propagation.end(), v)
               != propagation.end();
      }
    };

    struct SharedData {
      SimplexId extrema;
      AtomicVector<CurrentState *> states;
      AtomicVector<idSuperArc> openedArcs;

      explicit SharedData(SimplexId e)
        : extrema(e), states(50), openedArcs(50) {
      }

      void addState(CurrentState *curState) {
        const idThread &thisTask = states.getNext();
        states[thisTask] = curState;
      }

      void addArc(const idSuperArc arc) {
        idSuperArc thisArc = openedArcs.getNext();
        openedArcs[thisArc] = arc;
      }

      void merge(SharedData &other) {
        for(auto *state : other.states) {
          addState(state);
        }

        for(auto &arc : other.openedArcs) {
          addArc(arc);
        }
      }

      void reserve(const size_t &s) {
        states.reserve(s);
        openedArcs.reserve(s);
      }
    };

    struct Comparison {
      VertCompFN vertLower, vertHigher;
    };

    using segm_it = std::vector<SimplexId>::iterator;
    using segm_rev_it = std::vector<SimplexId>::reverse_iterator;
    using segm_const_it = std::vector<SimplexId>::const_iterator;
    using segm_const_rev_it = std::vector<SimplexId>::const_reverse_iterator;

    // Segmentation data
    struct Region {
      // inverted in case of split tree
      segm_it segmentBegin;
      segm_it segmentEnd;
    };

  } // namespace ftm
} // namespace ttk

#endif /* end of include guard: STRUCTURES_H */
