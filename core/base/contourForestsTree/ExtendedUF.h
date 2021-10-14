/// \ingroup base
/// \class ttk::ExtendedUnionFind
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.
///
/// \sa ttkContourForests.cpp %for a usage example.

#ifndef EXTENDEDUF_H
#define EXTENDEDUF_H

#include <vector>

#include "DeprecatedDataTypes.h"

namespace ttk {
  namespace cf {
    class ExtendedUnionFind {
    private:
      int rank_;
      ExtendedUnionFind *parent_;
      ufDataType data_;
      SimplexId origin_;

    public:
      inline ExtendedUnionFind(const SimplexId &origin) {
        rank_ = 0;
        parent_ = this;
        data_ = nullUfData;
        origin_ = origin;
      }

      inline ExtendedUnionFind(const ExtendedUnionFind &other) {
        rank_ = other.rank_;
        parent_ = this;
        data_ = other.data_;
        origin_ = other.origin_;
      }

      inline void setData(const ufDataType &d) {
        data_ = d;
      }

      inline void setOrigin(const SimplexId &origin) {
        origin_ = origin;
      }

      inline const ufDataType &getData(void) const {
        return data_;
      }

      inline const SimplexId &getOrigin(void) const {
        return origin_;
      }

      // heavy recursif
      inline ExtendedUnionFind *find() {
        if(parent_ == this)
          return this;
        else {
          parent_ = parent_->find();
          return parent_;
        }
      }

      inline int getRank() const {
        return rank_;
      }

      inline void setParent(ExtendedUnionFind *parent) {
        parent_ = parent;
      }

      inline void setRank(const int &rank) {
        rank_ = rank;
      }

      static inline ExtendedUnionFind *makeUnion(ExtendedUnionFind *uf0,
                                                 ExtendedUnionFind *uf1) {
        uf0 = uf0->find();
        uf1 = uf1->find();

        if(uf0 == uf1) {
          return uf0;
        } else if(uf0->getRank() > uf1->getRank()) {
          uf1->setParent(uf0);
          return uf0;
        } else if(uf0->getRank() < uf1->getRank()) {
          uf0->setParent(uf1);
          return uf1;
        } else {
          uf1->setParent(uf0);
          uf0->setRank(uf0->getRank() + 1);
          return uf0;
        }
      }

      static inline ExtendedUnionFind *
        makeUnion(std::vector<ExtendedUnionFind *> &sets) {
        ExtendedUnionFind *n = NULL;

        if(!sets.size())
          return NULL;

        if(sets.size() == 1)
          return sets[0];

        for(int i = 0; i < (int)sets.size() - 1; i++)
          n = makeUnion(sets[i], sets[i + 1]);

        return n;
      }

      inline bool operator<(const ExtendedUnionFind &other) const {
        return rank_ < other.rank_;
      }

      inline bool operator>(const ExtendedUnionFind &other) const {
        return rank_ > other.rank_;
      }
    };
  } // namespace cf
} // namespace ttk

#endif /* end of include guard: EXTENDEDUF_H */
