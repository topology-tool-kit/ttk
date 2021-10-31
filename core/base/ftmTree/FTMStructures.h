/// \ingroup base
//
/// \class ttk::FTMTree_MT
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date September 2016.
///
///\brief TTK structures for the contour tree

#pragma once

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
      SimplexId size{};
      void *values{};
      const SimplexId *offsets{};

      // [0] -> vertex id of the global miminum
      // [size-1] -> vertex id of the global maximum
      std::vector<SimplexId> sortedVertices{};

      inline bool isLower(const SimplexId a, const SimplexId b) const {
        return this->offsets[a] < this->offsets[b];
      }
      inline bool isEqLower(const SimplexId a, const SimplexId b) const {
        return this->offsets[a] <= this->offsets[b];
      }
      inline bool isHigher(const SimplexId a, const SimplexId b) const {
        return this->offsets[a] > this->offsets[b];
      }
      inline bool isEqHigher(const SimplexId a, const SimplexId b) const {
        return this->offsets[a] >= this->offsets[b];
      }
    };

    struct CurrentState {
      SimplexId vertex;
      boost::heap::fibonacci_heap<SimplexId, boost::heap::compare<VertCompFN>>
        propagation;

      CurrentState(SimplexId startVert, VertCompFN &vertComp)
        : vertex(startVert), propagation(vertComp) {
      }

      CurrentState(VertCompFN &vertComp)
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
      FTMAtomicVector<CurrentState *> states;
      FTMAtomicVector<idSuperArc> openedArcs;

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
