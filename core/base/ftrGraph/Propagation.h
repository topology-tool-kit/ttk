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

// C++ includes
#include <deque>

namespace ttk
{
   namespace ftr
   {
      class Propagation
      {
        private:
         // cache current simplex
         idVertex        curVert_;

         // comparison function
         VertCompFN comp_;
         // come from min/max leaf
         bool goUp_;

         // priority deque
         boost::heap::fibonacci_heap<idVertex, boost::heap::compare<VertCompFN>> propagation_;

         // representant (pos in array)
         AtomicUF id_;

         // Lazyness: do not impact the global DG
         // contains DG link between to node, lowest first
         // (see edge comparison from Parsa's paper)
         // along with the arc that visit this linkEdge
         std::deque<std::tuple<linkEdge, idSuperArc>> lazyAdd_, lazyDel_;

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

         // used to avoid thhe fibonacci heaps the local propagation is only
         // used to store comp and current vert during a sequential pass using
         // the reference alforithm and so it's a mock propagation used for
         // compliance with the rest of the code
         void setCurvert(const idVertex v)
         {
            curVert_ = v;
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
         bool compare(const idVertex a, const idVertex b) const
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

         void lazyAdd(const idEdge e0, const idEdge e1, const idSuperArc a)
         {
            lazyAdd_.emplace_back(std::make_tuple(std::make_pair(e0, e1), a));
         }

         void lazyDel(const idEdge e0, const idEdge e1)
         {
            // here the arc would be a non sense.
            lazyDel_.emplace_back(std::make_tuple(std::make_pair(e0, e1), nullSuperArc));
         }

         void sortLazyLists(LinkCompFN comp)
         {
            auto updateComp = [comp](const std::tuple<linkEdge, idSuperArc>& a,
                                     const std::tuple<linkEdge, idSuperArc>& b) {
               return comp(std::get<0>(a), std::get<0>(b));
            };
            std::sort(lazyAdd_.begin(), lazyAdd_.end(), updateComp);
            std::sort(lazyDel_.begin(), lazyDel_.end(), updateComp);
            // TODO maybe we could use only one list for both?
         }

         // return the head of lazyAdd / lazyDel (or a null link if empty)
         // would have used std::optional if possible
         std::tuple<linkEdge, idSuperArc> lazyAddNext()
         {
            if (lazyAdd_.empty()) {
               return  std::make_tuple(nullLink, nullSuperArc);
            } else {
               std::tuple<linkEdge, idSuperArc> add = lazyAdd_.front();
               lazyAdd_.pop_front();
               return add;
            }
         }

         std::tuple<linkEdge, idSuperArc> lazyDelNext()
         {
            if (lazyDel_.empty()) {
               return std::make_tuple(nullLink, nullSuperArc);
            } else {
               std::tuple<linkEdge, idSuperArc> del = lazyDel_.front();
               lazyDel_.pop_front();
               return del;
            }
         }

         bool lazyListsEmpty()
         {
            return lazyAdd_.empty() && lazyDel_.empty();
         }

         const decltype(lazyAdd_)& addEach() const
         {
            return lazyAdd_;
         }

         const decltype(lazyDel_)& delEach() const
         {
            return lazyDel_;
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
