#ifndef FTMTree_CT_Template_h_INCLUDED
#define FTMTree_CT_Template_h_INCLUDED

#include <FTMTree_CT.h>

namespace ttk {
  namespace ftm {

    template <class triangulationType>
    void FTMTree_CT::build(const triangulationType *mesh, TreeType tt) {
      Timer mergeTreesTime;

      const bool bothMT = tt == TreeType::Contour || tt == TreeType::Join_Split;

      initComp();

      if(bothMT) {
        // single leaf search for both tree
        // When executed from CT, both minima and maxima are extracted
        Timer precomputeTime;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif
        {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
          { leafSearch(mesh); }
        }
        printTime(precomputeTime, "leafSearch", 3);
      }

#ifdef TTK_ENABLE_OMP_PRIORITY
      {
        // Set priority
        if(st_->getNumberOfLeaves() < jt_->getNumberOfLeaves())
          st_->setPrior();
        else
          jt_->setPrior();
      }
#endif

      // JT & ST
      // clang-format off
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif
      {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
        {
          if(tt == TreeType::Join || bothMT) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task UNTIED() if(threadNumber_ > 1)
#endif
            jt_->build(mesh, tt == TreeType::Contour);
          }
          if(tt == TreeType::Split || bothMT) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task UNTIED() if(threadNumber_ > 1)
#endif
            st_->build(mesh, tt == TreeType::Contour);
          }
        }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
      }

      printTime(mergeTreesTime, "merge trees ", 3);

      // Combine
      if(tt == TreeType::Contour) {

        Timer combineFullTime;
        insertNodes();

        Timer combineTime;
        combine();
        printTime(combineTime, "combine trees", 4);
        printTime(combineFullTime, "combine full", 3);
      }
      // Debug
      if(debugLevel_ > 3) {
        std::string nbNodes;
        switch(tt) {
          case TreeType::Join:
            nbNodes = std::to_string(jt_->getNumberOfNodes());
            break;
          case TreeType::Split:
            nbNodes = std::to_string(st_->getNumberOfNodes());
            break;
          case TreeType::Join_Split:
            nbNodes
              = std::to_string(jt_->getNumberOfNodes() + st_->getNumberOfNodes());
            break;
          default:
            nbNodes = std::to_string(getNumberOfNodes());
        }
        this->printMsg({"- final number of nodes :", nbNodes});
      }
  }
// clang-format on
// clang format fail to use the right indentation level
// here, but it break the code if not disabled...

// ------------------------------------------------------------------------

template <class triangulationType>
int FTMTree_CT::leafSearch(const triangulationType *mesh) {
  const auto nbScalars = scalars_->size;
  const auto chunkSize = getChunkSize();
  const auto chunkNb = getChunkCount();

  // Extrema extract and launch tasks
  for(SimplexId chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId)
#endif
    {
      const SimplexId lowerBound = chunkId * chunkSize;
      const SimplexId upperBound
        = std::min(nbScalars, (chunkId + 1) * chunkSize);
      for(SimplexId v = lowerBound; v < upperBound; ++v) {
        const auto &neighNumb = mesh->getVertexNeighborNumber(v);
        valence upval = 0;
        valence downval = 0;

        for(valence n = 0; n < neighNumb; ++n) {
          SimplexId neigh;
          mesh->getVertexNeighbor(v, n, neigh);
          if(scalars_->isLower(neigh, v)) {
            ++downval;
          } else {
            ++upval;
          }
        }

        jt_->setValence(v, downval);
        st_->setValence(v, upval);

        if(!downval) {
          jt_->makeNode(v);
        }

        if(!upval) {
          st_->makeNode(v);
        }
      }
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
  return 0;
}

} // namespace ftm
} // namespace ttk
#endif // FTMTree_CT_Template_h_INCLUDED
