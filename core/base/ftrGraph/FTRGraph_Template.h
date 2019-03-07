#ifndef FTRGRAPH_TEMPLATE_H
#define FTRGRAPH_TEMPLATE_H

#include "FTRGraph.h"
#include "Propagation.h"
#include "Tasks.h"

#ifdef TTK_ENABLE_FTR_TASK_STATS
#include <iterator>
#endif
#ifdef TTK_ENABLE_FTR_VERT_STATS
#include <iterator>
#endif

#ifdef TTK_ENABLE_FTR_PRIORITY
#define OPTIONAL_PRIORITY(value) priority(value)
#else
#define OPTIONAL_PRIORITY(value)
#endif

#ifdef GPROFILE
#include <gperftools/profiler.h>
#endif

namespace ttk {
  namespace ftr {
    template <typename ScalarType>
    FTRGraph<ScalarType>::FTRGraph()
      : params_{}, scalars_{new Scalars<ScalarType>}, mesh_{} {
      // need a call to setupTriangulation later
    }

    template <typename ScalarType>
    FTRGraph<ScalarType>::FTRGraph(Triangulation *mesh)
      : params_{}, scalars_{new Scalars<ScalarType>}, mesh_{} {
      setupTriangulation(mesh);
    }

    template <typename ScalarType>
    FTRGraph<ScalarType>::~FTRGraph() {
      delete scalars_;
    }

    template <typename ScalarType>
    void FTRGraph<ScalarType>::build() {
      Timer t;

#ifndef TTK_ENABLE_KAMIKAZE
      if(!scalars_) {
        std::cerr << "[FTR Graph]: no scalars given" << std::endl;
        return;
      }
#endif
      // init some values

#ifdef TTK_ENABLE_OPENMP
      omp_set_num_threads(params_.threadNumber);
      omp_set_nested(1);
#ifdef TTK_ENABLE_FTR_PRIORITY
      if(omp_get_max_task_priority() < PriorityLevel::Max) {
        std::stringstream msg;
        msg << "[FTR Graph]: Warning, OpenMP max priority is lower than 5"
            << std::endl;
        dMsg(std::cerr, msg.str(), infoMsg);
      }
#endif
#endif

      params_.printSelf();

      // Precompute
      DebugTimer timeAlloc;
      alloc();
      printTime(timeAlloc, "[FTR Graph]: alloc time: ", infoMsg);

      DebugTimer timeInit;
      init();
      printTime(timeInit, "[FTR Graph]: init time: ", infoMsg);

      // std::cout << printMesh() << std::endl;
      // std::cout << mesh_.printEdges() << std::endl;

      DebugTimer finTime;
#ifdef GPROFILE
      std::cout << "Profiling enabled ..." << std::endl;
      ProfilerStart("ftr.log");
#endif

      DebugTimer timeSort;
      scalars_->sort();
      printTime(timeSort, "[FTR Graph]: sort time: ", infoMsg);

      DebugTimer timePreSortSimplices;
      mesh_.preSortEdges([&](const idVertex a, const idVertex b) {
        return scalars_->isLower(a, b);
      });
      mesh_.preSortTriangles([&](const idVertex a, const idVertex b) {
        return scalars_->isLower(a, b);
      });
      printTime(
        timePreSortSimplices, "[FTR Graph]: simplices sort time: ", infoMsg);

      // Build the graph

      DebugTimer timeBuild;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(params_.threadNumber)
#endif
      {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
        {
          // DebugTimer timeLeafSearch;
          // leafSearch();
          // printTime(timeLeafSearch, "[FTR Graph]: leaf search time: ",
          // timeMsg);

          DebugTimer timeCritSearch;
          criticalSearch();
          printTime(timeCritSearch, "[FTR Graph]: leaf search time ", timeMsg);

          DebugTimer timeSwipe;
          sweepFrowSeeds();
          // sweepSequential();
          printTime(timeSwipe, "[FTR Graph]: sweepFrowSeeds time: ", timeMsg);
        }
      }
      printTime(timeBuild, "[FTR Graph]: build time: ", timeMsg);

      // Debug print
#ifndef NDEBUG
      // std::cout << graph_.printVisit() << std::endl;
      // printGraph(4);
      // std::cout << "up" << std::endl;
      // std::cout << dynGraphs_.up.print() << std::endl;
      // std::cout << "down" << std::endl;
      // std::cout << dynGraphs_.down.print() << std::endl;
#endif

      // post-process
      DebugTimer postProcTime;
      graph_.mergeArcs<ScalarType>(scalars_);
      graph_.arcs2nodes<ScalarType>(scalars_);
      printTime(postProcTime, "[FTR Graph]: postProcess: ", advancedInfoMsg);

      // std::cout << "nb verts: " << mesh_.getNumberOfVertices() << std::endl;
      // std::cout << "nb triangle: " << mesh_.getNumberOfVertices() << std::endl;
      // std::cout << "nb nodes: " << graph_.getNumberOfNodes() << std::endl;
      // std::cout << "nb leaves: " << graph_.getNumberOfLeaves() << std::endl;
      // std::cout << "nb  arcs: " << graph_.getNumberOfVisibleArcs() << std::endl;

      printTime(finTime, "[FTR Graph]: *TOTAL* time: ", timeMsg);

#ifdef GPROFILE
      ProfilerStop();
#endif

      // list of regular vertices on each arc
      // explicit build: for sampling
      if(params_.samplingLvl) {
        graph_.buildArcSegmentation<ScalarType>(scalars_);
      }

      // Debug print
#ifndef NDEBUG
      // std::cout << graph_.printVisit() << std::endl;
      // printGraph(4);
      // std::cout << dynGraphs_.up.printNbCC() << std::endl;
#endif
      {
        std::stringstream msg;
        msg << "[FTR Graph]: " << graph_.getNumberOfVisibleArcs() << " arcs ("
            << graph_.getNumberOfArcs() << ")" << std::endl;
        dMsg(std::cout, msg.str(), advancedInfoMsg);
      }

#ifdef TTK_ENABLE_FTR_TASK_STATS
      std::cout << "propTimes_ :" << std::endl;
      copy(propTimes_.crbegin(), propTimes_.crend(),
           std::ostream_iterator<float>(std::cout, " "));
      std::cout << std::endl;
#endif

#ifdef TTK_ENABLE_FTR_VERT_STATS
      idVertex nbRevisit = graph_.getNbMultiTouch();
      std::cout << "revisit: " << nbRevisit << " / "
                << mesh_.getNumberOfVertices();
      std::cout << " = "
                << ((float)nbRevisit / mesh_.getNumberOfVertices()) * 100
                << std::endl;
      std::cout << "avoided: " << graph_.getNbAvoided() << std::endl;
#endif
    }

    // protected

    template <typename ScalarType>
    void FTRGraph<ScalarType>::leafSearch() {
      TaskChunk leafChunkParams(scalars_->getSize());
      leafChunkParams.grainSize = 10000;
      auto leafChunk = Tasks::getChunk(leafChunkParams);

      for(idPropagation leafChunkId = 0; leafChunkId < std::get<1>(leafChunk);
          ++leafChunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(leafChunkId)
#endif
        {
          const idVertex lowerBound
            = Tasks::getBegin(leafChunkId, std::get<0>(leafChunk));
          const idVertex upperBound = Tasks::getEnd(
            leafChunkId, std::get<0>(leafChunk), scalars_->getSize());

          for(idVertex v = lowerBound; v < upperBound; ++v) {
            const valence vNeighNumber = mesh_.getVertexNeighborNumber(v);
            bool isMax = false;
            bool isMin = true;

            for(valence n = 0; n < vNeighNumber; ++n) {
              idVertex neigh;
              mesh_.getVertexNeighbor(v, n, neigh);

              if(scalars_->isHigher(neigh, v)) {
                isMax = false;
              } else {
                isMin = false;
              }
            }

            // v is a minimum, add it to the leaves
            if(isMin || isMax) {
              graph_.addLeaf(v, isMin);
            }
          }
        } // end task
      }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
#ifdef TTK_ENABLE_FTR_TASK_STATS
      // Stats
      nbProp_ = graph_.getNumberOfLeaves();
      propTimes_.resize(nbProp_);
#endif
      {
        std::stringstream msg;
        msg << "[FTR Graph]: " << graph_.getNumberOfLeaves() << " leaves"
            << std::endl;
        dMsg(std::cout, msg.str(), infoMsg);
      }
    }

    template <typename ScalarType>
    void FTRGraph<ScalarType>::criticalSearch() {
      const bool addMin = true;
      const bool addMax = !params_.singleSweep;

      TaskChunk leafChunkParams(scalars_->getSize());
      leafChunkParams.grainSize = 10000;
      auto leafChunk = Tasks::getChunk(leafChunkParams);

      auto comp = [&](const idVertex a, const idVertex b) {
        return scalars_->isLower(a, b);
      };

      for(idPropagation leafChunkId = 0; leafChunkId < std::get<1>(leafChunk);
          ++leafChunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(leafChunkId)
#endif
        {
          const idVertex lowerBound
            = Tasks::getBegin(leafChunkId, std::get<0>(leafChunk));
          const idVertex upperBound = Tasks::getEnd(
            leafChunkId, std::get<0>(leafChunk), scalars_->getSize());

          // each task uses its local forests
          LocalForests localForests;
          localForests.up.setNumberOfNodes(60);
          localForests.up.init();
          localForests.up.alloc();
          localForests.down.setNumberOfNodes(60);
          localForests.down.init();
          localForests.down.alloc();

          for(idVertex v = lowerBound; v < upperBound; ++v) {
            std::tie(valences_.lower[v], valences_.upper[v])
              = getLinkNbCC(v, localForests, comp);

            // leaf cases
            if(addMin && valences_.lower[v] == 0) {
              graph_.addLeaf(v, true);
            }
            if(addMax && valences_.upper[v] == 0) {
              graph_.addLeaf(v, false);
            }
          }
        } // end task
      }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
#ifdef TTK_ENABLE_FTR_TASK_STATS
      // Stats
      nbProp_ = graph_.getNumberOfLeaves();
      propTimes_.resize(nbProp_);
#endif
      {
        std::stringstream msg;
        msg << "[FTR Graph]: " << graph_.getNumberOfLeaves() << " leaves"
            << std::endl;
        dMsg(std::cout, msg.str(), infoMsg);
      }
    }

    template <typename ScalarType>
    void FTRGraph<ScalarType>::sweepFrowSeeds() {
      const idNode nbSeed = graph_.getNumberOfLeaves();
      // used to interleave min and max
      // Note: useless if only start for min or max

      graph_.sortLeaves<ScalarType>(scalars_);

#ifdef TTK_ENABLE_FTR_TASK_STATS
      sweepStart_.reStart();
#endif
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskgroup
#endif
      {
        for(idNode i = 0; i < nbSeed; i++) {
          // alterneate min/max, string at the deepest
          idVertex l = (i % 2) ? i / 2 : nbSeed - 1 - i / 2;
          // increasing order, min first
          // idVertex l = i;
          // initialize structure
          const idVertex corLeaf = graph_.getLeaf(l);
          const bool fromMin = graph_.isLeafFromMin(l);
          Propagation *localPropagation = newPropagation(corLeaf, fromMin);
          const idSuperArc newArc
            = graph_.openArc(graph_.makeNode(corLeaf), localPropagation);
          // graph_.visit(corLeaf, newArc);
          // process
#ifdef TTK_ENABLE_OPENMP
#pragma omp task OPTIONAL_PRIORITY(PriorityLevel::Higher)
#endif
          growthFromSeed(corLeaf, localPropagation, newArc);
        }
      }
    }

    template <typename ScalarType>
    void FTRGraph<ScalarType>::sweepSequential() {
      growthSequential(0, mesh_.getNumberOfVertices());
    }

    template <typename ScalarType>
    void FTRGraph<ScalarType>::alloc() {
      mesh_.alloc();

      scalars_->setSize(mesh_.getNumberOfVertices());
      scalars_->alloc();

      graph_.setNumberOfElmt(mesh_.getNumberOfVertices());
      graph_.alloc();

      propagations_.setNumberOfElmt(mesh_.getNumberOfVertices());
      propagations_.alloc();

      dynGraphs_.up.setNumberOfElmt(mesh_.getNumberOfEdges());
      dynGraphs_.up.alloc();

      dynGraphs_.down.setNumberOfElmt(mesh_.getNumberOfEdges());
      dynGraphs_.down.alloc();

#ifndef TTK_DISABLE_FTR_LAZY
      lazy_.setNumberOfElmt(mesh_.getNumberOfVertices() * 2);
      lazy_.alloc();
#endif

      valences_.lower.resize(mesh_.getNumberOfVertices());
      valences_.upper.resize(mesh_.getNumberOfVertices());

#ifdef TTK_ENABLE_FTR_BFS
      bfsCells_.resize(mesh_.getNumberOfTriangles());
      bfsEdges_.resize(mesh_.getNumberOfEdges());
      bfsVerts_.resize(mesh_.getNumberOfVertices());
#endif
    }

    template <typename ScalarType>
    void FTRGraph<ScalarType>::init() {
      scalars_->removeNaN();
      scalars_->init();
      graph_.init();
      propagations_.init();
      dynGraphs_.up.init();
      dynGraphs_.down.init();

#ifndef TTK_DISABLE_FTR_LAZY
      lazy_.init();
#endif

#ifdef TTK_ENABLE_FTR_BFS
      fillVector<idCell>(bfsCells_, nullCell);
      fillVector<idEdge>(bfsEdges_, nullEdge);
      fillVector<idVertex>(bfsVerts_, nullVertex);
#endif

      // Stats
#ifdef TTK_ENABLE_FTR_TASK_STATS
      fillVector<float>(propTimes_, 0);
#endif
    }

  } // namespace ftr
} // namespace ttk

#endif /* end of include guard: FTRGRAPH_TEMPLATE_H */
