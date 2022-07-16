#pragma once

#include "FTRGraph.h"
#include "FTRPropagation.h"
#include "FTRTasks.h"

#ifdef TTK_ENABLE_FTR_TASK_STATS
#include <iterator>
#endif
#ifdef TTK_ENABLE_FTR_VERT_STATS
#include <iterator>
#endif

#ifdef TTK_ENABLE_OMP_PRIORITY
#define OPTIONAL_PRIORITY(value) priority(value)
#else
#define OPTIONAL_PRIORITY(value)
#endif

#ifdef GPROFILE
#include <gperftools/profiler.h>
#endif

namespace ttk {
  namespace ftr {
    template <typename ScalarType, typename triangulationType>
    FTRGraph<ScalarType, triangulationType>::FTRGraph() {
      this->setDebugMsgPrefix("FTRGraph");
    }

    template <typename ScalarType, typename triangulationType>
    FTRGraph<ScalarType, triangulationType>::FTRGraph(triangulationType *mesh) {
      this->setDebugMsgPrefix("FTRGraph");
      preconditionTriangulation(mesh);
    }

    template <typename ScalarType, typename triangulationType>
    void FTRGraph<ScalarType, triangulationType>::build() {
      Timer t;

      // init some values

#ifdef TTK_ENABLE_OPENMP
      ParallelGuard pg{params_.threadNumber};
      omp_set_nested(1);
#ifdef TTK_ENABLE_OMP_PRIORITY
      if(omp_get_max_task_priority() < PriorityLevel::Max) {
        this->printWrn("OpenMP max priority is lower than 5");
      }
#endif
#endif

      params_.printSelf();

      // Precompute
      Timer timeAlloc;
      alloc();
      this->printMsg(
        "alloc time: ", 1.0, timeAlloc.getElapsedTime(), this->threadNumber_);

      Timer timeInit;
      init();
      this->printMsg(
        "init time: ", 1.0, timeInit.getElapsedTime(), this->threadNumber_);

      // std::cout << printMesh() << std::endl;
      // std::cout << mesh_.printEdges() << std::endl;

      Timer finTime;
#ifdef GPROFILE
      std::cout << "Profiling enabled ..." << std::endl;
      ProfilerStart("ftr.log");
#endif

      Timer timeSort;
      this->printMsg(
        "sort time: ", 1.0, timeSort.getElapsedTime(), this->threadNumber_);

      Timer timePreSortSimplices;
      mesh_.preSortEdges([&](const idVertex a, const idVertex b) {
        return scalars_.isLower(a, b);
      });
      mesh_.preSortTriangles([&](const idVertex a, const idVertex b) {
        return scalars_.isLower(a, b);
      });
      this->printMsg("simplices sort time: ", 1.0,
                     timePreSortSimplices.getElapsedTime(),
                     this->threadNumber_);

      // Build the graph

      Timer timeBuild;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(params_.threadNumber)
#endif
      {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
        {
          Timer timeCritSearch;
          criticalSearch();
          this->printMsg("leaf search time ", 1.0,
                         timeCritSearch.getElapsedTime(), this->threadNumber_);

          Timer timeSwipe;
          sweepFrowSeeds();
          // sweepSequential();
          this->printMsg("sweepFrowSeeds time: ", 1.0,
                         timeSwipe.getElapsedTime(), this->threadNumber_);
        }
      }
      this->printMsg(
        "build time: ", 1.0, timeBuild.getElapsedTime(), this->threadNumber_);

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
      Timer postProcTime;
      graph_.mergeArcs<ScalarType>(&scalars_);
      graph_.arcs2nodes<ScalarType>(&scalars_);
      this->printMsg("postProcess: ", 1.0, postProcTime.getElapsedTime(),
                     this->threadNumber_);

      // std::cout << "nb verts: " << mesh_.getNumberOfVertices() << std::endl;
      // std::cout << "nb triangle: " << mesh_.getNumberOfVertices() <<
      // std::endl; std::cout << "nb nodes: " << graph_.getNumberOfNodes() <<
      // std::endl; std::cout << "nb leaves: " << graph_.getNumberOfLeaves() <<
      // std::endl; std::cout << "nb  arcs: " << graph_.getNumberOfVisibleArcs()
      // << std::endl;

      this->printMsg(
        "*TOTAL* time: ", 1.0, finTime.getElapsedTime(), this->threadNumber_);

#ifdef GPROFILE
      ProfilerStop();
#endif

      // list of regular vertices on each arc
      // explicit build: for sampling
      if(params_.samplingLvl) {
        graph_.buildArcSegmentation<ScalarType>(&scalars_);
      }

      // Debug print
#ifndef NDEBUG
      // std::cout << graph_.printVisit() << std::endl;
      // printGraph(4);
      // std::cout << dynGraphs_.up.printNbCC() << std::endl;
#endif
      this->printMsg(
        std::vector<std::vector<std::string>>{
          {"#Visible arcs", std::to_string(graph_.getNumberOfVisibleArcs())},
          {"#Arcs", std::to_string(graph_.getNumberOfArcs())},
        },
        debug::Priority::DETAIL);

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

    template <typename ScalarType, typename triangulationType>
    void FTRGraph<ScalarType, triangulationType>::criticalSearch() {
      const bool addMin = true;
      const bool addMax = !params_.singleSweep;

      ScalarFieldCriticalPoints critPoints;

      TaskChunk leafChunkParams(scalars_.getSize());
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
            leafChunkId, std::get<0>(leafChunk), scalars_.getSize());

          // each task uses its local forests
          for(idVertex v = lowerBound; v < upperBound; ++v) {
            bool lBoundary = false, uBoundary = false;
            critPoints.getNumberOfLowerUpperComponents(
              v, scalars_.getOffsets(), mesh_.getTriangulation(),
              valences_.lower[v], valences_.upper[v], lBoundary, uBoundary);

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
      this->printMsg(std::vector<std::vector<std::string>>{
        {"#Leaves", std::to_string(graph_.getNumberOfLeaves())}});
    }

    template <typename ScalarType, typename triangulationType>
    void FTRGraph<ScalarType, triangulationType>::sweepFrowSeeds() {
      const idNode nbSeed = graph_.getNumberOfLeaves();
      // used to interleave min and max
      // Note: useless if only start for min or max

      graph_.sortLeaves<ScalarType>(&scalars_);

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

    template <typename ScalarType, typename triangulationType>
    void FTRGraph<ScalarType, triangulationType>::sweepSequential() {
      growthSequential(0, mesh_.getNumberOfVertices());
    }

    template <typename ScalarType, typename triangulationType>
    void FTRGraph<ScalarType, triangulationType>::alloc() {
      mesh_.alloc();

      scalars_.setSize(mesh_.getNumberOfVertices());
      scalars_.alloc();

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
    }

    template <typename ScalarType, typename triangulationType>
    void FTRGraph<ScalarType, triangulationType>::init() {
      scalars_.removeNaN();
      scalars_.init();
      graph_.init();
      propagations_.init();
      dynGraphs_.up.init();
      dynGraphs_.down.init();

#ifndef TTK_DISABLE_FTR_LAZY
      lazy_.init();
#endif

      // Stats
#ifdef TTK_ENABLE_FTR_TASK_STATS
      fillVector<float>(propTimes_, 0);
#endif
    }

  } // namespace ftr
} // namespace ttk
