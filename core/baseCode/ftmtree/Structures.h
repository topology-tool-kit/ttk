/// \ingroup baseCode
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

#include "AtomicVector.h"
#include "DataTypes.h"

// todo remove
#include<iostream>

namespace ttk
{
   // Compute parameters (global)
   struct Params {
      TreeType      treeType;
   };

   // Scalar related containers (global)
   struct Scalars {
      idVertex size;
      void*    values;

      // Actually, fields below are unused -----
      std::shared_ptr<std::vector<idVertex>> sosOffsets;
      std::shared_ptr<std::vector<idVertex>> sortedVertices, mirrorVertices;

      // Need vertices to be sorted : use mirrorVertices.

      bool isLower(idVertex a, idVertex b) const
      {
         return (*mirrorVertices)[a] < (*mirrorVertices)[b];
      }
      bool isEqLower(idVertex a, idVertex b) const
      {
         return (*mirrorVertices)[a] <= (*mirrorVertices)[b];
      }

      bool isHigher(idVertex a, idVertex b) const
      {
         return (*mirrorVertices)[a] > (*mirrorVertices)[b];
      }
      bool isEqHigher(idVertex a, idVertex b) const
      {
         return (*mirrorVertices)[a] >= (*mirrorVertices)[b];
      }

      Scalars() : sosOffsets(nullptr), sortedVertices(nullptr), mirrorVertices(nullptr)
      {
      }

      // Heavy
      Scalars(const Scalars& o)
          : sosOffsets(o.sosOffsets),
            sortedVertices(o.sortedVertices),
            mirrorVertices(o.mirrorVertices)
      {
          std::cout << "copy in depth" << std::endl;
      }

      // Sort
      template <typename type>
      void qsort(type arr[], const long int begin, const long int stop,
                 std::function<bool(type, type)> comp) const
      {
         if (begin >= stop)
            return;

         static const long int MINSIZE = 10;

         long int   left  = begin - 1;
         long int   right = stop + 1;
         const type pivot = arr[begin];

         while (1) {
            while (comp(pivot, arr[--right]))
               ;
            while (++left <= stop && !comp(pivot, arr[left]))
               ;

            if (left < right)
               swap_el<type>(arr, left, right);
            else
               break;
         }

         swap_el<type>(arr, begin, right);
#pragma omp task untied if (right - begin > MINSIZE)
         qsort(arr, begin, right - 1, comp);
#pragma omp task untied if (stop - right > MINSIZE)
         qsort(arr, right + 1, stop, comp);
      }

     private:
      template <typename type>
      inline void swap_el(type arr[], const size_t a, const size_t b) const
      {
         const type tmp = arr[a];
         arr[a]         = arr[b];
         arr[b]         = tmp;
        }
   };

   struct CurrentState {
      idVertex vertex;
      boost::heap::fibonacci_heap<idVertex, boost::heap::compare<VertCompFN>> propagation;

      CurrentState(idVertex startVert, VertCompFN vertComp)
          : vertex(startVert), propagation(vertComp)
      {
      }

      idVertex getNextMinVertex(void)
      {
          vertex = propagation.top();
          propagation.pop();
          return vertex;
      }

      void addNewVertex(const idVertex v)
      {
          propagation.emplace(v);
      }

      void merge(CurrentState& other)
      {
         propagation.merge(other.propagation);
         vertex = propagation.top();
      }

      bool empty()
      {
         return propagation.empty();
      }

      // DEBUG ONLY
      bool find(idVertex v)
      {
         return std::find(propagation.begin(), propagation.end(), v) != propagation.end();
      }
   };

   struct SharedData {
      idVertex                    extrema;
      AtomicVector<CurrentState*> states;
      AtomicVector<idSuperArc>    openedArcs;

      SharedData(idVertex e) : extrema(e), states(50), openedArcs(50)
      {
      }

      void addState(CurrentState* curState)
      {
         const idThread& thisTask = states.getNext();
         states[thisTask]         = curState;
      }

      void addArc(const idSuperArc arc)
      {
         idSuperArc thisArc  = openedArcs.getNext();
         openedArcs[thisArc] = arc;
      }

      void merge(SharedData& other)
      {
         for (auto* state : other.states) {
            addState(state);
         }

         for (auto& arc : other.openedArcs) {
            addArc(arc);
         }
      }

      void reserve(const size_t& s)
      {
          states.reserve(s);
          openedArcs.reserve(s);
      }
   };

   struct Comparison {
      VertCompFN vertLower, vertHigher;
   };

   using segm_it           = std::vector<idVertex>::iterator;
   using segm_rev_it       = std::vector<idVertex>::reverse_iterator;
   using segm_const_it     = std::vector<idVertex>::const_iterator;
   using segm_const_rev_it = std::vector<idVertex>::const_reverse_iterator;

   // Segmentation data
   struct Region {
      // inverted in case of split tree
      segm_it segmentBegin;
      segm_it segmentEnd;
   };
}

#endif /* end of include guard: STRUCTURES_H */
