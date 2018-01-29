/// \ingroup base
/// \class ttk::ftr::Propagation
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-15
///
/// \brief TTK %fTRGraph propagation management with Fibonacci heaps
///
/// This class deal with scalar related operations: store them, compare them, ...
///
/// \sa ttk::FTRGraph

#pragma once

#include "FTRCommon.h"

#include <boost/heap/fibonacci_heap.hpp>

namespace ttk
{
   namespace ftr
   {
      using VertCompFN = std::function<bool(idVertex, idVertex)>;

      class Propagation : public Allocable
      {
        private:
         idVertex curVert_;
         // priority queue
         boost::heap::fibonacci_heap<idVertex, boost::heap::compare<VertCompFN>> propagation_;

        public:
         Propagation(idVertex startVert, VertCompFN vertComp)
             : curVert_(startVert), propagation_(vertComp)
         {
         }

         idVertex getNextMinVertex(void)
         {
            curVert_ = propagation_.top();
            propagation_.pop();
            return curVert_;
         }

         void addNewVertex(const idVertex v)
         {
            propagation_.emplace(v);
         }

         void merge(Propagation& other)
         {
            propagation_.merge(other.propagation_);
            curVert_ = propagation_.top();
         }

         bool empty()
         {
            return propagation_.empty();
         }

         // DEBUG ONLY
         bool find(idVertex v)
         {
            return std::find(propagation_.begin(), propagation_.end(), v) != propagation_.end();
         }

         // Initialize functions
         // --------------------

         void alloc() override;

         void init() override;
      };
   }
}
