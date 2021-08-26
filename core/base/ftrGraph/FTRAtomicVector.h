/// \ingroup base
/// \class ttk::FTRAtomicVector
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-01-22
///
///\brief TTK processing package that manage a paralle version of vector *Same
/// as in FTM: Common ?*

#ifndef FTRATOMICVECTOR_H
#define FTRATOMICVECTOR_H

#include "BaseClass.h"

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
  class FTRAtomicVector : public std::vector<type> {
  private:
    std::size_t nextId;

  public:
    explicit FTRAtomicVector(const std::size_t initSize = 1)
      : std::vector<type>(), nextId(0) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!initSize) {
        std::cout << "Caution, Atomic vector need a non-0 init size !"
                  << std::endl;
        std::vector<type>::resize(initSize);
      } else
#endif
      {
        std::vector<type>::resize(initSize);
      }
    }

    // copy constructor
    FTRAtomicVector(const FTRAtomicVector &other)
      : std::vector<type>(other), nextId(other.nextId) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!std::vector<type>::size()) {
        reserve(1);
      }
#endif
    }

    // move constructor
    FTRAtomicVector(FTRAtomicVector &&other) noexcept = default;

    virtual ~FTRAtomicVector() = default;

    // ---
    // STL
    // ---

    void reserve(const std::size_t &newSize, const bool fromOther = false) {
      if(newSize > std::vector<type>::size()) {
#ifndef TTK_ENABLE_KAMIKAZE
#ifdef TTK_ENABLE_OPENMP
        if(omp_in_parallel()) {
          // WARNING: In parallel we do not want to make reserve as it can lead
          // to data race, we should not enter here
#pragma omp critical(AtomicUFReserve)
          {
            // if (fromOther)
            //    std::cout << " a function in the class ";
            // std::cout << "call RE-Reserve in FTRAtomicVector " << nextId;
            // std::cout << " ! Data Race may occurs ! " << typeid(this).name()
            // << std::endl;

            std::vector<type>::resize(newSize);
          }

        } else
#endif
#endif
        {
          std::vector<type>::resize(newSize);
        }
      }
      TTK_FORCE_USE(fromOther);
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
      reserve(oldSize, true);
    }

    std::size_t getNext(void) {
      std::size_t resId;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
      resId = nextId++;

      if(nextId == std::vector<type>::size()) {
        reserve(std::vector<type>::size() * 2, true);
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

    void emplace_back(const type &elmt) {
      // Not really constructed in place
      // here to ensure vector compatibility
      const auto &curPos = getNext();
      (*this)[curPos] = elmt;
    }

    void pop_back(void) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
      --nextId;
    }

    // --------
    // OPERATOR
    // --------

    FTRAtomicVector<type> &operator=(const FTRAtomicVector<type> &other) {
      if(&other != this) {
        std::vector<type>::operator=(other);
        nextId = other.nextId;
      }
      return *this;
    }

    FTRAtomicVector<type> &operator=(FTRAtomicVector<type> &&other) noexcept {
      nextId = std::move(other.nextId);
      std::vector<type>::operator=(std::move(other));
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

    const_iterator end() const {
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

    const_riterator rbegin() const {
      return this->rend() - (nextId - 1);
    }

    const_riterator crbegin() const {
      return this->crend() - (nextId - 1);
    }
  };
} // namespace ttk

#endif /* end of include guard: FTRATOMICVECTOR_H */
