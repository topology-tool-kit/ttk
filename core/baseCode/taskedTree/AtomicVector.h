/// \ingroup baseCode
/// \class ttk::AtomicVector
/// \author Charles Gueuent <charles.gueunet@lip6.fr>
/// \date 2017-02-09
///
///\brief TTK processing package that manage a paralle vecrion of vector

#ifndef ATOMICVECTOR_H
#define ATOMICVECTOR_H

#include <omp.h>
#include <iterator>
#include <vector>

#ifndef withKamikaze
# include <iostream>
# include <typeinfo>
#endif

template <typename type>
class AtomicVector : public std::vector<type>
{
  private:
   std::size_t nextId;

  public:
   AtomicVector(const std::size_t initSize = 1) : std::vector<type>(), nextId(0)
   {
#ifndef withKamikaze
      if (!initSize) {
          std::cout << "Caution, Atomic vector need a non-0 init size !" << std::endl;
          std::vector<type>::resize(initSize);
      } else
#endif
      {
          std::vector<type>::resize(initSize);
      }
   }

   // copy constructor
   AtomicVector(const AtomicVector &other) : std::vector<type>(other), nextId(other.nextId)
   {
#ifndef withKamikaze
      if (!std::vector<type>::size()) {
         reserve(1);
      }
#endif
   }

   virtual ~AtomicVector() = default;

   // ---
   // STL
   // ---

   void reserve(const std::size_t &newSize, const bool fromOther = false)
   {
      if (newSize > std::vector<type>::size()) {
#ifndef withKamikaze
# ifdef withOpenMP
         if (omp_in_parallel()) {
// In parallel we do not want to make reserve as it can lead to
// data race
# pragma omp critical(AtomicUFReserve)
            {
               if (fromOther)
                  std::cout << " a function in the class ";
               std::cout << "call RE-Reserve in AtomicVector " << nextId;
               std::cout << " ! Data Race may occurs ! " << typeid(this).name() << std::endl;

               std::vector<type>::resize(newSize);
            }

         } else
# endif
#endif
         {
            std::vector<type>::resize(newSize);
         }
      }
   }

   void reset(const std::size_t &nId = 0)
   {
#pragma omp atomic write
      nextId = nId;
   }

   void clear(void)
   {
      reset();

      // Remove old content
      std::size_t oldSize = std::vector<type>::size();
      std::vector<type>::clear();
      reserve(oldSize, true);
   }

   std::size_t getNext(void)
   {
      std::size_t resId;
#pragma omp atomic capture
      resId = nextId++;

      if (nextId == std::vector<type>::size()) {
         reserve(std::vector<type>::size() * 2, true);
      }

      return resId;
   }

   std::size_t size(void) const
   {
      return nextId;
   }

   void push_back(const type &elmt)
   {
       const auto& curPos = getNext();
       operator[](curPos) = elmt;
   }

   // ---------
   // ITERATORS
   // ---------
   // allow foreach on the vector

   typedef typename std::vector<type>::iterator iterator;
   typedef typename std::vector<type>::const_iterator const_iterator;

   iterator end()
   {
      return this->begin() + nextId;
   }

   const_iterator cend() const
   {
      return this->cbegin() + nextId;
   }
};

#endif /* end of include guard: ATOMICVECTOR_H */
