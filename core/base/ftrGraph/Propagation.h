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

#include <Triangulation.h>
#include <UnionFind.h>

#include <boost/heap/fibonacci_heap.hpp>

namespace ttk
{
   namespace ftr
   {
      using VertCompFN = std::function<bool(idVertex, idVertex)>;

      class Propagation
      {
        private:
         // cached extrema: propagation_.top() will be the next one
         idVertex   curVert_;
         // comparison function
         VertCompFN comp_;
         // come from min/max leaf
         bool goUp_;
         // priority queue
         boost::heap::fibonacci_heap<idVertex, boost::heap::compare<VertCompFN>> propagation_;
         // representant
         UnionFind* rpz_;

        public:
         Propagation(idVertex startVert, VertCompFN vertComp, bool up, UnionFind* uf = nullptr)
             : curVert_(nullVertex), comp_(vertComp), goUp_(up), propagation_(vertComp)
         {
            rpz_ = (uf) ? uf : new UnionFind;
            propagation_.emplace(startVert);
         }

         idVertex getCurVertex(void) const
         {
            return curVert_;
         }

         UnionFind* getRpz(void) const
         {
            return rpz_->find();
         }

         idVertex nextVertex(void)
         {
            curVert_ = propagation_.top();
            propagation_.pop();
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
            while (propagation_.top() == d){
               propagation_.pop();
            }
         }

         void addNewVertex(const idVertex v)
         {
            propagation_.emplace(v);
         }

         void merge(Propagation& other)
         {
            propagation_.merge(other.propagation_);
            rpz_ = makeUnion(rpz_, other.rpz_);
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

#ifndef NDEBUG
         std::string print() const
         {
            std::stringstream res;
            boost::heap::fibonacci_heap<idVertex, boost::heap::compare<VertCompFN>> tmp(propagation_);
            res << "localProp " << curVert_ << " : ";
            while(!tmp.empty()){
               res << tmp.top() << " ";
               tmp.pop();
            }
            return res.str();
         }
#endif
      };
   }
}
