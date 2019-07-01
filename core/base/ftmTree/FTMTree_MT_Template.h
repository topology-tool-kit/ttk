/// \ingroup base
/// \class ttk::FTMTree_MT
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the
/// sublevel set tree of scalar data and more
/// (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa ttkContourForests.cpp %for a usage example.

#ifndef FTMTREE_MT_TPL_H
#define FTMTREE_MT_TPL_H

#include <functional>
#if(defined(__GNUC__) && !defined(__clang__))
#include <parallel/algorithm>
#endif

#include "FTMTree_MT.h"

// ----
// Init
// ----

namespace ttk {
  namespace ftm {

    template <typename scalarType, typename idType>
    void ftm::FTMTree_MT::sortInput(void) {
      const auto &nbVertices = scalars_->size;

      auto *sortedVect = scalars_->sortedVertices.get();
      if(sortedVect == nullptr) {
        sortedVect = new std::vector<SimplexId>(0);
        scalars_->sortedVertices.reset(sortedVect);
      } else {
        sortedVect->clear();
      }

      auto indirect_sort = [&](const size_t &a, const size_t &b) {
        return isLower<scalarType, idType>(a, b);
      };

      sortedVect->resize(nbVertices, 0);
      std::iota(
        sortedVect->begin(), sortedVect->end(), static_cast<SimplexId>(0));

      // #pragma omp parallel
      // #pragma omp single
      //    scalars_->qsort<SimplexId>(sortedVect->data(), 0, scalars_->size -1,
      //    indirect_sort);

#ifdef TTK_ENABLE_OPENMP
#ifdef __clang__
      std::cout << "Caution, outside GCC, sequential sort" << std::endl;
      std::sort(sortedVect->begin(), sortedVect->end(), indirect_sort);
#else
      __gnu_parallel::sort(
        sortedVect->begin(), sortedVect->end(), indirect_sort);
#endif
#else
      std::sort(sortedVect->begin(), sortedVect->end(), indirect_sort);
#endif

      auto *mirrorVert = scalars_->mirrorVertices.get();
      if(mirrorVert == nullptr) {
        mirrorVert = new std::vector<SimplexId>(0);
        scalars_->mirrorVertices.reset(mirrorVert);
      } else {
        mirrorVert->clear();
      }

      scalars_->mirrorVertices->resize(nbVertices);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
      for(SimplexId i = 0; i < nbVertices; i++) {
        (*scalars_->mirrorVertices)[(*sortedVect)[i]] = i;
      }
    }

  } // namespace ftm
} // namespace ttk
// Process

#endif /* end of include guard: FTMTREE_MT_TPL_H */
