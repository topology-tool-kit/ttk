/// \ingroup base
/// \class ttk::ftr::AtomicUF
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-07-16
///
///\brief Union find compliant with parallel find and maintaining
/// the current local propagation
/// inspired by ttk::UnionFind
///
/// \sa ftrGraph

#ifndef ATOMICUFFTR_H
#define ATOMICUFFTR_H

// c++ includes
#include <memory>
#include <vector>

namespace ttk {
  namespace ftr {

    class Propagation;

    class AtomicUF {
    private:
      unsigned rank_;
      AtomicUF *parent_;
      Propagation *prop_;

    public:
      inline explicit AtomicUF(Propagation *const p) : rank_(0), prop_{p} {
        parent_ = this;
      }

      // heavy recursif
      inline AtomicUF *find() {
        if(parent_ == this)
          return this;
        else {
          decltype(parent_) tmp = parent_->find();
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
          parent_ = tmp;

          return parent_;
        }
      }

      // Shared data get/set

      inline Propagation *getPropagation(void) {
        return prop_;
      }

      inline void setPropagation(Propagation *const p) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
        prop_ = p;
      }

      // UF get / set

      inline int getRank() const {
        return rank_;
      }

      inline void setRank(const int &rank) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
        rank_ = rank;
      }

      inline void setParent(AtomicUF *parent) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
        parent_ = parent;
      }

      static inline AtomicUF *makeUnion(AtomicUF *uf0, AtomicUF *uf1) {
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

      static inline AtomicUF *makeUnion(std::vector<AtomicUF *> &sets) {
        AtomicUF *n = NULL;

        if(!sets.size())
          return NULL;

        if(sets.size() == 1)
          return sets[0];

        for(int i = 0; i < (int)sets.size() - 1; i++) {
          n = makeUnion(sets[i], sets[i + 1]);
        }

        return n;
      }

      inline bool operator<(const AtomicUF &other) const {
        return rank_ < other.rank_;
      }

      inline bool operator>(const AtomicUF &other) const {
        return rank_ > other.rank_;
      }
    };

  } // namespace ftr
} // namespace ttk

#endif /* end of include guard: ATOMICUFFTR_H */
