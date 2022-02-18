/// \ingroup base
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
/// \sa ttkFTMTree.cpp %for a usage example.

#ifndef FTMTREE_TPL_H
#define FTMTREE_TPL_H

#include "FTMTree.h"

// for std::isnan
#include <cmath>
// for std::numeric_limits
#include <limits>

// -------
// PROCESS
// -------

template <typename scalarType, class triangulationType>
void ttk::ftm::FTMTree::build(const triangulationType *mesh) {
  // -----
  // INPUT
  // -----

  printParams();

#ifdef TTK_ENABLE_OPENMP
  ParallelGuard pg{threadNumber_};
  omp_set_nested(1);
#ifdef TTK_ENABLE_OMP_PRIORITY
  if(omp_get_max_task_priority() < 5) {
    this->printWrn("OpenMP max priority is lower than 5");
  }
#endif
#endif

  // ----
  // INIT
  // ----

  setDebugLevel(debugLevel_);
  initNbScalars(mesh);

  // This section is aimed to prevent un-deterministic results if the data-set
  // have NaN values in it.
  // In this loop, we replace every NaN by a 0 value.
  // Recall: Equals values are distinguished using Simulation of Simplicity in
  // the FTM tree computation Note: Can we detect NaN using vtk ?
  if(::std::numeric_limits<scalarType>::has_quiet_NaN) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
    for(SimplexId i = 0; i < scalars_->size; i++) {
      if(::std::isnan((double)(((scalarType *)scalars_->values)[i]))) {
        ((scalarType *)scalars_->values)[i] = 0;
      }
    }
  }

  // Alloc / reserve
  Timer initTime;
  switch(params_->treeType) {
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
  printTime(initTime, "alloc", 3);

  Timer startTime;

  // init values
  Timer setTimer;
  switch(params_->treeType) {
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
  printTime(setTimer, "init", 3);

  // for fast comparison
  // and regions / segmentation
  Timer sortTime;
  sortInput<scalarType>();
  printTime(sortTime, "sort step", 3);

  // -----
  // BUILD
  // -----

  Timer buildTime;
  FTMTree_CT::build(mesh, params_->treeType);
  printTime(buildTime, "build tree", 3);

  printTime(startTime, "Total ", 1);

#ifdef PERF_TESTS
  exit(0);
#endif

  // Build the list of regular vertices of the arc
  if(params_->segm) {
    switch(params_->treeType) {
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
    switch(params_->treeType) {
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

  if(debugLevel_ > 4) {
    switch(params_->treeType) {
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
