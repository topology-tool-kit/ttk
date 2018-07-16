/// \class ttk::ftr::Propagation
/// \ingroup base
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-15
///
/// \brief TTK %fTRGraph propagation management with Fibonacci heaps
///
/// This class deal with scalar related operations: store them, compare them, ...
///
/// \sa ttk::FTRGraph

#ifndef PROPAGATION_H
#define PROPAGATION_H


// local include
#include "FTRCommon.h"
#include "AtomicUF.h"

// base code includes
#include <Triangulation.h>

// library include
#include <boost/heap/fibonacci_heap.hpp>

namespace ttk
{
   namespace ftr
   {
      using VertCompFN = std::function<bool(idVertex, idVertex)>;

      class Propagation
      {
        private:
         // cache current simplex
         idVertex        curVert_;

         // comparison function
         VertCompFN comp_;
         // come from min/max leaf
         bool goUp_;

         // priority queue
         boost::heap::fibonacci_heap<idVertex, boost::heap::compare<VertCompFN>> propagation_;

         // representant (pos in array)
         AtomicUF id_;

        public:
         Propagation(idVertex startVert, VertCompFN vertComp, bool up)
             : curVert_{nullVertex},
               comp_{vertComp},
               goUp_{up},
               propagation_{vertComp},
               id_{this}
         {
            propagation_.emplace(startVert);
         }

         Propagation(const Propagation& other) = delete;

         idVertex getCurVertex(void) const
         {
            return curVert_;
         }

         AtomicUF* getId(void)
         {
            return id_.find();
         }

         idVertex nextVertex(void)
         {
            curVert_ = propagation_.top();
            removeDuplicates(curVert_);
            return curVert_;
         }

         idVertex getNextVertex(void) const
         {
#ifndef TTK_ENABLE_KAMIKAZE
            if (propagation_.empty()) {
               std::cerr << "[FTR]: Propagation get next on empty structure" << std::endl;
               return nullVertex;
            }
#endif
            return propagation_.top();
         }

         void removeDuplicates(const idVertex d) {
            while (!propagation_.empty() && propagation_.top() == d) {
               propagation_.pop();
            }
         }

         void removeBelow(const idVertex d, VertCompFN comp) {
            while (!propagation_.empty() && comp(propagation_.top(), d)) {
               propagation_.pop();
            }
         }

         void addNewVertex(const idVertex v)
         {
            propagation_.emplace(v);
         }

         void merge(Propagation& other)
         {
            if (&other == this) return;
            propagation_.merge(other.propagation_);
            AtomicUF::makeUnion(&id_, &other.id_);
            // TODO once after all the merge ?
            id_.find()->setPropagation(this);
         }

         bool empty() const
         {
            return propagation_.empty();
         }

         void clear(void)
         {
            propagation_.clear();
         }

         /// This comparison is reversed to the internal one
         /// because when we are growing, the natural order for the simplices
         /// is reversed to the order of the priority queue:
         /// We want the min to grow form lowest to highest
         bool compare(idVertex a, idVertex b) const
         {
            return comp_(b,a);
         }

         // true if propagation initiated form a min leaf
         bool goUp(void) const
         {
            return goUp_;
         }

         bool goDown(void) const
         {
            return !goUp_;
         }

         // DEBUG ONLY

         bool find(idVertex v) const
         {
            return std::find(propagation_.begin(), propagation_.end(), v) != propagation_.end();
         }

         std::size_t size() const
         {
            return propagation_.size();
         }

         std::string print() const
         {
            std::stringstream res;
            res << this << " localProp " << curVert_ << " : ";
#ifndef NDEBUG
            // only if perf are not important: copy
            boost::heap::fibonacci_heap<idVertex, boost::heap::compare<VertCompFN>> tmp(propagation_);
            while(!tmp.empty()){
               res << tmp.top() << " ";
               tmp.pop();
            }
#else
            res << "...";
#endif
            return res.str();
         }
      };
   }
}

#endif /* end of include guard: PROPAGATION_H */
