
/// \ingroup baseCode
/// \class ttk::ContourTree
/// \author Charles Gueuent <charles.gueunet@lip6.fr>
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
/// \sa vtkTaskedTree.cpp %for a usage example.

#ifndef TASKEDTREETEMPLATE_H
#define TASKEDTREETEMPLATE_H

#include "TaskedTree.h"

// -------
// PROCESS
// -------
// {

template <typename scalarType>
void TaskedTree::build(void)
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
   if (params_->treeType == TreeType::Join || params_->treeType == TreeType::Contour) {
      getJoinTree()->makeAlloc();
   }
   if (params_->treeType == TreeType::Split || params_->treeType == TreeType::Contour) {
      getSplitTree()->makeAlloc();
   }
   if (params_->treeType == TreeType::Contour) {
      makeAlloc();
   }
   printTime(initTime, "alloc step", -1, 3);

   DebugTimer startTime;

   // init values
   DebugTimer setTimer;
   if (params_->treeType == TreeType::Join || params_->treeType == TreeType::Contour) {
      getJoinTree()->makeInit();
   }
   if (params_->treeType == TreeType::Split || params_->treeType == TreeType::Contour) {
      getSplitTree()->makeInit();
   }
   if (params_->treeType == TreeType::Contour) {
      makeInit();
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
   ContourTree::build(params_->treeType);
   printTime(buildTime, "9 build tree", -1, 1);

   printTime(startTime, "10 TOTAL ", -1, 1);

   if(params_->treeType == TreeType::Join){
       getJoinTree()->buildSegmentation();
       getJoinTree()->finalizeSegmentation();
   }
   if(params_->treeType == TreeType::Split){
       getSplitTree()->buildSegmentation();
       getSplitTree()->finalizeSegmentation();
   }

   // exit(0);
}

// }

#endif /* end of include guard: TASKEDTREETEMPLATE_H */
