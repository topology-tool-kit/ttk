/// \ingroup base
/// \class ttk::UnionFind
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date July 2011.
///
/// \brief Union Find implementation for connectivity tracking.

#ifndef _UNION_FIND_H
#define _UNION_FIND_H

#include <Debug.h>

#include <algorithm>
#include <vector>

namespace ttk {

  class UnionFind : virtual public Debug {

  public:
    // 1) constructors, destructors, operators
    inline UnionFind();

    inline UnionFind(const UnionFind &other);

    inline bool operator<(const UnionFind &other) const {
      return rank_ < other.rank_;
    };

    inline bool operator>(const UnionFind &other) const {
      return rank_ > other.rank_;
    };

    // 2) functions
    inline UnionFind *find();

    inline int getRank() const {
      return rank_;
    };

    static inline UnionFind *makeUnion(UnionFind *uf0, UnionFind *uf1);

    static inline UnionFind *makeUnion(std::vector<UnionFind *> &sets);

    inline void setParent(UnionFind *parent) {
      parent_ = parent;
    };

    inline void setRank(const int &rank) {
      rank_ = rank;
    };

  protected:
    int rank_;
    UnionFind *parent_;
  };

  inline UnionFind::UnionFind() {

    rank_ = 0;
    parent_ = this;
  }

  inline UnionFind::UnionFind(const UnionFind &other) {

    rank_ = other.rank_;
    parent_ = this;
  }

  inline UnionFind *UnionFind::find() {

    if(parent_ == this)
      return this;
    else {
      parent_ = parent_->find();
      return parent_;
    }
  }

  static inline UnionFind *makeUnion(UnionFind *uf0, UnionFind *uf1) {

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

    return NULL;
  }

  static inline UnionFind *makeUnion(std::vector<UnionFind *> &sets) {

    UnionFind *n = NULL;

    if(!sets.size())
      return NULL;

    if(sets.size() == 1)
      return sets[0];

    for(int i = 0; i < (int)sets.size() - 1; i++)
      n = makeUnion(sets[i], sets[i + 1]);

    return n;
  }

} // namespace ttk

#endif
