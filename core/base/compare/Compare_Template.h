#ifndef COMPARE_TEMPLATE_H
#define COMPARE_TEMPLATE_H

#include "Compare.h"

namespace ttk
{
   template <typename Type>
   int Compare::computeVertDiff(void* const scalArray1, void* const scalArray2)
   {
      if (vertMapperM1toM2_.empty()) {
         std::cerr << "[Compare] comuteVertDiff needs vertMapper to be filled" << std::endl;
      }

      Type*          inputArray1 = (Type*)scalArray1;
      Type*          inputArray2 = (Type*)scalArray2;
      const idVertex nbVerts1    = mesh1_->getNumberOfVertices();
      const idVertex nbVerts2    = mesh2_->getNumberOfVertices();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for (idVertex i = 0; i < nbVerts1; ++i) {
         const idVertex posInM2 = vertMapperM1toM2_[i];
         if (posInM2 < 0 || posInM2 >= nbVerts2 || inputArray2[posInM2] != inputArray1[i]) {
            inputArray1[i] = 1;
         } else {
            inputArray1[i] = 0;
         }
      }

      return 0;
   }

   template <typename Type>
   int Compare::computeCellDiff(void* const scalArray1, void* const scalArray2)
   {
      if (cellMapperM1toM2_.empty()) {
         std::cerr << "[Compare] comuteCellDiff needs cellMapper to be filled" << std::endl;
      }

      Type*          inputArray1 = (Type*)scalArray1;
      Type*          inputArray2 = (Type*)scalArray2;
      const idCell nbCells1    = mesh1_->getNumberOfCells();
      const idCell nbCells2    = mesh2_->getNumberOfCells();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for (idCell i = 0; i < nbCells1; ++i) {
         const idCell posInM2 = cellMapperM1toM2_[i];
         if (posInM2 < 0 || posInM2 >= nbCells2 || inputArray2[posInM2] != inputArray1[i]) {
            inputArray1[i] = 1;
         } else {
            inputArray1[i] = 0;
         }
      }

      return 0;
   }
}

#endif /* end of include guard: COMPARE_TEMPLATE_H */
