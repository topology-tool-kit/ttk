/// \ingroup base
//
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-01-22
///
///\brief TTK structures for the reeb graph

#pragma once

#include "FTRDataTypes.h"

#include <Debug.h>

#ifdef __APPLE__
#include <algorithm>
#include <numeric>
#else
#ifdef _WIN32
#include <algorithm>
#include <numeric>
#else
#ifdef __clang__
#include <algorithm>
#include <numeric>
#else
#include <parallel/algorithm>
#endif
#endif
#endif

#include <iostream>
#include <vector>

namespace ttk {
  namespace ftr {
    // Compute parameters (global)
    struct Params {
      bool singleSweep = false;
      bool segm = true;
      bool normalize = true;
      bool advStats = true;
      int samplingLvl = 0;

      idThread threadNumber = 1;
      int debugLevel = 1;

      void printSelf() {
        if(debugLevel) {
          std::cout << "[FTR Graph]: thread number: " << threadNumber
                    << std::endl;
          std::cout << "[FTR Graph]: debug lvl: " << debugLevel << std::endl;
          if(debugLevel > 2) {
            std::cout << "[FTR Graph]: segmentation: " << std::boolalpha << segm
                      << std::endl;
            std::cout << "[FTR Graph]: sampling level: " << samplingLvl
                      << std::endl;
          }
        }
      }
    };

    /// Force the class to use functions alloc and init
    /// alloc is used for all system allocation
    /// init is used to fill arrays that needs to be
    class Allocable : public Debug {
    protected:
      /// Allocation may depends on the number of vertices
      idVertex nbElmt_ = nullVertex;

    public:
      void setNumberOfElmt(const idVertex nbVerts) {
        nbElmt_ = nbVerts;
      }

      template <typename type>
      void fillVector(std::vector<type> &vect, const type &elmt) {
        const std::size_t nbIt = vect.size();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for(std::size_t i = 0; i < nbIt; i++) {
          vect[i] = elmt;
        }
      }

      virtual void alloc() = 0;

      virtual void init() = 0;
    };

    template <typename Iterator>
    void parallel_sort(Iterator begin, Iterator end) {
      // Sort the vertices array
#ifdef TTK_ENABLE_OPENMP
#ifdef _GLIBCXX_PARALLEL_FEATURES_H
      ::__gnu_parallel::sort(begin, end);
#else
      std::cout << "Caution, sequential sort" << std::endl;
      ::std::sort(begin, end);
#endif
#else
      ::std::sort(begin, end);
#endif
    }

    template <typename Iterator, typename el>
    void parallel_sort(Iterator begin,
                       Iterator end,
                       std::function<bool(el, el)> comp) {
      // Sort the vertices array
#ifdef TTK_ENABLE_OPENMP
#ifdef _GLIBCXX_PARALLEL_FEATURES_H
      ::__gnu_parallel::sort(begin, end, comp);
#else
      std::cout << "Caution, sequential sort" << std::endl;
      ::std::sort(begin, end, comp);
#endif
#else
      ::std::sort(begin, end, comp);
#endif
    }

    template <typename Iterator>
    void sort(Iterator begin, Iterator end) {
      ::std::sort(begin, end);
    }

    template <typename Iterator, typename el>
    void sort(Iterator begin, Iterator end, std::function<bool(el, el)> comp) {
      ::std::sort(begin, end, comp);
    }
  } // namespace ftr
} // namespace ttk
