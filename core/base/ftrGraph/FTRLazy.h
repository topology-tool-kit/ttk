
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
#include "FTRCommon.h"
#include "FTRDataTypes.h"

// c++ includes
#include <set>
#include <vector>

namespace ttk {
  namespace ftr {
    class Lazy : public Allocable {
    private:
      std::vector<std::set<linkEdge>> lazyAdd_;

    public:
      Lazy() {
      }

      void alloc() override {
        lazyAdd_.resize(nbElmt_);
      }

      void init() override {
        // nothing
      }

      void addEmplace(const idEdge e0, const idEdge e1, const idSuperArc a) {
        lazyAdd_[a].emplace(std::make_pair(e0, e1));
      }

      void delEmplace(const idEdge e0, const idEdge e1, const idSuperArc a) {
        const auto p = std::make_pair(e0, e1);
        auto it = lazyAdd_[a].find(p);
        if(it != lazyAdd_[a].end()) {
          lazyAdd_[a].erase(it);
        }
      }

      // return the head of lazyAdd / lazyDel (or a null link if empty)
      // would have used std::optional if possible
      linkEdge addGetNext(const idSuperArc a) {
        if(lazyAdd_[a].empty()) {
          return nullLink;
        } else {
          linkEdge add = *lazyAdd_[a].begin();
          lazyAdd_[a].erase(lazyAdd_[a].begin());
          return add;
        }
      }

      bool isEmpty(const idSuperArc a) {
        return lazyAdd_[a].empty();
      }

      // const decltype(lazyAdd_)& addEach() const
      // {
      //    return lazyAdd_;
      // }

      // const decltype(lazyDel_)& delEach() const
      // {
      //    return lazyDel_;
      // }
    };

  } // namespace ftr
} // namespace ttk
