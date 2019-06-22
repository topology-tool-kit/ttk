/// \ingroup base
/// \class ttk::AtomicVector
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-01-22
///
///\brief TTK processing package that manage a paralle version of vector *Same
/// as in FTM: Common ?*

#ifndef ATOMICVECTOR_H
#define ATOMICVECTOR_H

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
  class AtomicVector : public std::vector<type> {
  private:
    std::size_t nextId;

  public:
    explicit AtomicVector(const std::size_t initSize = 1)
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
    AtomicVector(const AtomicVector &other)
      : std::vector<type>(other), nextId(other.nextId) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!std::vector<type>::size()) {
        reserve(1);
      }
#endif
    }

    // move constructor
    AtomicVector(AtomicVector &&other) = default;

    virtual ~AtomicVector() = default;

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
            // std::cout << "call RE-Reserve in AtomicVector " << nextId;
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

    AtomicVector<type> &operator=(const AtomicVector<type> &other) {
      std::vector<type>::operator=(other);
      nextId = other.nextId;
      return *this;
    }

    AtomicVector<type> &operator=(AtomicVector<type> &&other) {
      std::vector<type>::operator=(std::move(other));
      nextId = std::move(other.nextId);
      return *this;
    }

    // ---------
    // ITERATORS
    // ---------
    // allow foreach on the vector

    typedef typename std::vector<type>::iterator iterator;
    typedef typename std::vector<type>::const_iterator const_iterator;

    iterator end() {
      return this->begin() + nextId;
    }

    const_iterator end() const {
      return this->begin() + nextId;
    }

    const_iterator cend() const {
      return this->cbegin() + nextId;
    }

    typedef typename std::vector<type>::reverse_iterator riterator;
    typedef typename std::vector<type>::const_reverse_iterator const_riterator;

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

#endif /* end of include guard: ATOMICVECTOR_H */
