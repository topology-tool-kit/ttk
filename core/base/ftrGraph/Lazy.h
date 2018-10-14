
/// \ingroup base
/// \class ttk::ftr::Lazy
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-07-31
///
/// \brief TTK %DynamicGraph lazyness
///
/// This class maintain the information per arc
/// for the lazy evalutaion of each DG component.
///
/// \sa ttk::FTRGraph

#pragma once

// local includes
#include "DataTypesFTR.h"
#include "FTRCommon.h"

// c++ includes
#include <deque>
#include <vector>

namespace ttk
{
   namespace ftr
   {
      class Lazy : public Allocable
      {
        private:
         // When TTK_FTR_SINGLE_LIST is ON, we use the simple optimization
         // described in the original algorithm (parsa). In that
         // case we have a vector or size one. This is less efficient.
         // More over, the single list mode does not work in parallel,
         // as an update can impact other parts of the graph
         std::vector<std::deque<linkEdge>> lazyAdd_, lazyDel_;

        public:
         Lazy()
         {
         }

         void alloc() override
         {
#ifndef TTK_FTR_SINGLE_LIST
            lazyAdd_.resize(nbElmt_);
            lazyDel_.resize(nbElmt_);
#else
            lazyAdd_.resize(1);
            lazyDel_.resize(1);
#endif
         }

         void init() override
         {
            // nothing
         }

         void addEmplace(const idEdge e0, const idEdge e1, idSuperArc a)
         {
#ifdef TTK_FTR_SINGLE_LIST
            a = 0;
#endif
            lazyAdd_[a].emplace_back(std::make_pair(e0, e1));
         }

         void delEmplace(const idEdge e0, const idEdge e1, idSuperArc a)
         {
#ifdef TTK_FTR_SINGLE_LIST
            a = 0;
#endif
            // here the arc would be a non sense.
            lazyDel_[a].emplace_back(std::make_pair(e0, e1));
         }

         // return the head of lazyAdd / lazyDel (or a null link if empty)
         // would have used std::optional if possible
         linkEdge addGetNext(idSuperArc a)
         {
#ifdef TTK_FTR_SINGLE_LIST
            a = 0;
#endif
            if (lazyAdd_[a].empty()) {
               return nullLink;
            } else {
               linkEdge add = lazyAdd_[a].front();
               lazyAdd_[a].pop_front();
               return add;
            }
         }

         linkEdge delGetNext(idSuperArc a)
         {
#ifdef TTK_FTR_SINGLE_LIST
            a = 0;
#endif
            if (lazyDel_[a].empty()) {
               return nullLink;
            } else {
               linkEdge del = lazyDel_[a].front();
               lazyDel_[a].pop_front();
               return del;
            }
         }

         bool isEmpty(idSuperArc a)
         {
#ifdef TTK_FTR_SINGLE_LIST
            a = 0;
#endif
            return lazyAdd_[a].empty() && (lazyDel_[a].empty());
         }

         const decltype(lazyAdd_)& addEach() const
         {
            return lazyAdd_;
         }

         const decltype(lazyDel_)& delEach() const
         {
            return lazyDel_;
         }
      };

   }  // namespace ftr
}  // namespace ttk
