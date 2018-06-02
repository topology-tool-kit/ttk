/// \ingroup base
/// \class ttk::Compare
/// \author Charles Gueunet
/// \date 2018-05-31
///
/// \brief TTK %compare processing package.
///
/// %Compare is a TTK processing package that
/// compare two data sets in terms of mesh and scalars
///
/// \sa ttk::Triangulation
/// \sa ttkCompare.cpp %for a usage example.

#ifndef COMPARE_H
#define COMPARE_H

#include <Triangulation.h>
#include <Wrapper.h>

#include "CompareTypes.h"

namespace ttk
{
   class Compare : public Debug
   {
     protected:
      Triangulation *       mesh1_, *mesh2_;
      std::vector<idVertex> vertMapperM1toM2_;
      std::vector<idCell>   cellMapperM1toM2_;

     public:
      Compare();

      ~Compare();

      int setupTriangulations(Triangulation *const m1, Triangulation *const m2)
      {
         mesh1_ = m1;
         mesh2_ = m2;
         return 0;
      }

      // return 0 if mesh are the same,
      // return 1 if vertices differs
      // return 2 if cells differs
      int computeMeshDiff(unsigned char* const vertArr, unsigned char* const cellArr);

      // return 0 if scalars are the same (use mesh1 to mesh2 mapper)
      // return 1 if there is a differemce
      // fill scalArray1 with 0 for identical scalars and 1 when differences
      template <typename Type>
      int computeVertDiff(void *const scalArray1, void *const scalArray2);

      template <typename Type>
      int computeCellDiff(void *const scalArray1, void *const scalArray2);

     private:
      // fill  vertMapperM1toM2_ and vertArr accordingly
      bool computeVertsDiff(unsigned char* const vertArr);

      // fill cellMapperM1toM2_ and cellArr accordingly
      // Nees vertMapperM1toM2_ to be filled
      bool computeCellDiff(unsigned char* const cellArr);
   };

}  // namespace ttk

#include "Compare_Template.h"

#endif /* end of include guard: COMPARE_H */
