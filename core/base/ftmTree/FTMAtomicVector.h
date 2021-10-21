/// \ingroup base
/// \class ttk::FTMAtomicVector
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2017-02-09
///
///\brief TTK processing package that manage a paralle vecrion of vector

#ifndef FTMATOMICVECTOR_H
#define FTMATOMICVECTOR_H

#ifdef TTK_ENABLE_OPENMP
#include <omp.h>
#endif // TTK_ENABLE_OPENMP
#include <iterator>
#include <vector>

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#include <typeinfo>
#endif

namespace ttk {
  template <typename type>
  class FTMAtomicVector : public std::vector<type> {
  private:
    std::size_t nextId;
    // for initialization
    const type defaultValue;

  public:
    FTMAtomicVector(const std::size_t initSize = 1, const type &dv = type{})
      : std::vector<type>(), nextId(0), defaultValue{dv} {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!initSize) {
        std::cout << "Caution, Atomic vector need a non-0 init size !"
                  << std::endl;
        std::vector<type>::resize(1, defaultValue);
      } else
#endif
      {
        std::vector<type>::resize(initSize, defaultValue);
      }
    }

    // copy constructor
    FTMAtomicVector(const FTMAtomicVector &other)
      : std::vector<type>(other), nextId(other.nextId) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!std::vector<type>::size()) {
        reserve(1);
      }
#endif
    }

    FTMAtomicVector(FTMAtomicVector &&other) noexcept = default;

    virtual ~FTMAtomicVector() = default;

    // ---
    // STL
    // ---

    // If we do not want to use default constructor for elements
    // pre-allocated
    void reserve(const std::size_t &newSize) {
      if(newSize > std::vector<type>::size()) {
#ifndef TTK_ENABLE_KAMIKAZE
#ifdef TTK_ENABLE_OPENMP
        if(omp_in_parallel()) {
          // WARNING: In parallel we do not want to make reserve as it can lead
          // to data race, we should not enter here
#pragma omp critical(AtomicUFReserve)
          { std::vector<type>::resize(newSize, defaultValue); }

        } else
#endif
#endif
        {
          std::vector<type>::resize(newSize, defaultValue);
        }
      }
    }

    void reset(const std::size_t &nId = 0) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
      nextId = nId;
    }

    void clear(void) {
      reset();

      // Remove old content
      std::size_t oldSize = std::vector<type>::size();
      std::vector<type>::clear();
      reserve(oldSize);
    }

    std::size_t getNext(void) {
      std::size_t resId;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
      resId = nextId++;

      if(nextId == std::vector<type>::size()) {
        reserve(std::vector<type>::size() * 2);
      }

      return resId;
    }

    std::size_t size(void) const {
      return nextId;
    }

    bool empty(void) const {
      return nextId == 0;
    }

    void push_back(const type &elmt) {
      const auto &curPos = getNext();
      (*this)[curPos] = elmt;
    }

    // --------
    // OPERATOR
    // --------

    FTMAtomicVector<type> &operator=(const FTMAtomicVector<type> &other) {
      if(&other != this) {
        std::vector<type>::operator=(other);
        nextId = other.nextId;
      }
      return *this;
    }

    // ---------
    // ITERATORS
    // ---------
    // allow foreach on the vector

    using iterator = typename std::vector<type>::iterator;
    using const_iterator = typename std::vector<type>::const_iterator;

    iterator end() {
      return this->begin() + nextId;
    }

    const_iterator cend() const {
      return this->cbegin() + nextId;
    }

    using riterator = typename std::vector<type>::reverse_iterator;
    using const_riterator = typename std::vector<type>::const_reverse_iterator;

    riterator rbegin() {
      return this->rend() - (nextId - 1);
    }

    const_riterator crbegin() const {
      return this->crend() - (nextId - 1);
    }
  };
} // namespace ttk

#endif /* end of include guard: FTMATOMICVECTOR_H */
