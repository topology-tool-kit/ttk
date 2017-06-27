/// \ingroup baseCode
/// \class ttk::FTMTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date Dec 2016.
///
///\brief TTK processing package that efficiently computes the
/// contour tree of scalar data and more
/// (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa vtkFTMTree.cpp %for a usage example.

#ifndef FTMTREE_TPL_H
#define FTMTREE_TPL_H

#include "FTMTree.h"

// -------
// PROCESS
// -------

template <typename scalarType>
void FTMTree::build(void)
{
   // -----
   // INPUT
   // -----

   printParams();

#ifdef withOpenMP
   omp_set_num_threads(threadNumber_);
   omp_set_nested(1);
#endif

   // ----
   // INIT
   // ----

   setDebugLevel(debugLevel_);
   initNbScalars();

   // Alloc / reserve
   DebugTimer initTime;
   switch (params_->treeType) {
      case TreeType::Join:
         getJoinTree()->makeAlloc();
         break;
      case TreeType::Split:
         getSplitTree()->makeAlloc();
         break;
      case TreeType::Join_Split:
         getJoinTree()->makeAlloc();
         getSplitTree()->makeAlloc();
         break;
      case TreeType::Contour:
         getJoinTree()->makeAlloc();
         getSplitTree()->makeAlloc();
         makeAlloc();
         break;
      default:
         break;
   }
   printTime(initTime, "alloc step", -1, 3);

   DebugTimer startTime;

   // init values
   DebugTimer setTimer;
   switch (params_->treeType) {
      case TreeType::Join:
         getJoinTree()->makeInit();
         break;
      case TreeType::Split:
         getSplitTree()->makeInit();
         break;
      case TreeType::Join_Split:
         getJoinTree()->makeInit();
         getSplitTree()->makeInit();
         break;
      case TreeType::Contour:
         getJoinTree()->makeInit();
         getSplitTree()->makeInit();
         makeInit();
         break;
      default:
         break;
   }
   printTime(setTimer, "0 init");

   // for fast comparison
   // and regions / segmentation
   DebugTimer sortTime;
   initSoS();
   sortInput<scalarType>();
   printTime(sortTime, "1 sort step", -1, 1);

   // -----
   // BUILD
   // -----

   DebugTimer buildTime;
   FTMTree_CT::build(params_->treeType);
   printTime(buildTime, "9 build tree", -1, 1);

   printTime(startTime, "10 TOTAL ", -1, 1);
   // exit(0);

   // Build the list of regular vertices of the arc
   if (params_->segm) {
      switch (params_->treeType) {
         case TreeType::Join:
            getJoinTree()->buildSegmentation();
            getJoinTree()->finalizeSegmentation();
            break;
         case TreeType::Split:
            getSplitTree()->buildSegmentation();
            getSplitTree()->finalizeSegmentation();
            break;
         case TreeType::Join_Split:
            getJoinTree()->buildSegmentation();
            getSplitTree()->buildSegmentation();
            getJoinTree()->finalizeSegmentation();
            getSplitTree()->finalizeSegmentation();
            break;
         case TreeType::Contour:
            finalizeSegmentation();
            break;
         default:
            break;
      }
   }

   // Normalization
  if(params_->normalize) {
      switch (params_->treeType) {
         case TreeType::Join:
            getJoinTree()->normalizeIds();
            break;
         case TreeType::Split:
            getSplitTree()->normalizeIds();
            break;
         case TreeType::Join_Split:
            getJoinTree()->normalizeIds();
            getSplitTree()->normalizeIds();
            break;
         case TreeType::Contour:
            normalizeIds();
            break;
         default:
            break;
      }
  }

   if (debugLevel_ > 3) {
      switch (params_->treeType) {
         case TreeType::Join:
            jt_->printTree2();
            break;
         case TreeType::Split:
            st_->printTree2();
            break;
         case TreeType::Join_Split:
            jt_->printTree2();
            st_->printTree2();
            break;
         default:
            printTree2();
      }
   }
}

#endif /* end of include guard: FTMTREE_TPL_H */
