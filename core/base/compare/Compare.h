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

#pragma once

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
      unsigned char *       diffVerts_, *diffCells_;

      bool hasDiffVerts_, hasDiffCells_;

     public:
      Compare();

      ~Compare();

      int setupTriangulations(Triangulation *const m1, Triangulation *const m2)
      {
         mesh1_ = m1;
         mesh2_ = m2;
         return 0;
      }

      void setVertsArray(void *arr)
      {
         diffVerts_ = (unsigned char *)arr;
      }

      void setCellsArray(void *arr)
      {
         diffCells_ = (unsigned char *)arr;
      }

      int computeMeshDiff(void);

      template <typename Type>
      int computeVertDiff(void *const scalArray1, void *const scalArray2);

      template <typename Type>
      int computeCellDiff(void *const scalArray1, void *const scalArray2);

     private:
      // fill  vertMapperM1toM2_ and diffVerts_ accordingly
      void computeVertsDiff(void);

      // fill cellMapperM1toM2_ and diffCells_ accordingly
      // Nees vertMapperM1toM2_ to be filled
      void computeCellDiff(void);
   };

}  // namespace ttk

#include "Compare_Template.h"
