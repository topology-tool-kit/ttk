/// \ingroup baseCode
/// \class ttk::MergeTree
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
/// \sa vtkContourForests.cpp %for a usage example.

#ifndef MERGETREETEMPLATE_H
#define MERGETREETEMPLATE_H

#include <functional>

#include "MergeTree.h"

// ----
// Init
// ----
// {

template <typename scalarType>
void MergeTree::sortInput(void)
{
   const auto &nbVertices = scalars_->size;

   auto *sortedVect = scalars_->sortedVertices.get();
   if (sortedVect == nullptr) {
      sortedVect = new vector<idVertex>(0);
      scalars_->sortedVertices.reset(sortedVect);
   } else {
      sortedVect->clear();
   }

   auto indirect_sort = [&](const size_t &a, const size_t &b) { return isLower<scalarType>(a, b); };

   sortedVect->resize(nbVertices, 0);
   std::iota(sortedVect->begin(), sortedVect->end(), 0);

// #pragma omp parallel
// #pragma omp single
//    scalars_->qsort<idVertex>(sortedVect->data(), 0, scalars_->size -1, indirect_sort);

#ifdef withOpenMP
# ifdef __clang__
   cout << "Caution, outside GCC, sequential sort" << endl;
   std::sort(sortedVect->begin(), sortedVect->end(), indirect_sort);
# else
   __gnu_parallel::sort(sortedVect->begin(), sortedVect->end(), indirect_sort);
# endif
#else
   std::sort(sortedVect->begin(), sortedVect->end(), indirect_sort);
#endif

   auto *mirrorVert = scalars_->mirrorVertices.get();
   if (mirrorVert == nullptr) {
      mirrorVert = new vector<idVertex>(0);
      scalars_->mirrorVertices.reset(mirrorVert);
   } else {
      mirrorVert->clear();
   }

   scalars_->mirrorVertices->resize(nbVertices);

#pragma omp parallel for
   for (idVertex i = 0; i < nbVertices; i++) {
      (*scalars_->mirrorVertices)[(*sortedVect)[i]] = i;
   }
}
// }

// Process
// {
// }

#endif /* end of include guard: MERGETREETEMPLATE_H */

