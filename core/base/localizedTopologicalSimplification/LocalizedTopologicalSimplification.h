/// \ingroup base
/// \class ttk::LocalizedTopologicalSimplification
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 09/23/2020
///
/// This module provides functions to compute a localized topological
/// simplification of a scalar array. This algorithm is an order of magnitude
/// faster than the default simplification procedure as LTS utlizes
/// shared-memory parallelism, and constraints computationally expensive
/// procedures onto the subsets of the domain that actually need to be
/// simplified.
///
/// \b Related \b publication: \n
/// 'Localized Topological Simplification of Scalar Data'
/// Jonas Lukasczyk, Christoph Garth, Ross Maciejewski and Julien Tierny.
/// IEEE Transactions on Visualization and Computer Graphics.
/// 2020.
///

#pragma once

// for parallel sort
#if(defined(__GNUC__) && !defined(__clang__))
#include <parallel/algorithm>
#endif

// for numerical perturbation
#include <boost/math/special_functions/next.hpp>

#include <Debug.h>
#include <Triangulation.h>

#include <SSCPropagation.h>

namespace ttk {

  inline std::string toFixed(const float &number, const int precision = 2) {
    std::stringstream vFraction;
    vFraction << std::fixed << std::setprecision(precision) << number;
    return vFraction.str();
  }

  template <typename IT>
  std::string
    toFixed(const IT &number0, const IT &number1, const int precision = 2) {
    return toFixed(((float)number0) / ((float)number1), precision);
  };

  class LocalizedTopologicalSimplification : virtual public Debug {

  public:
    LocalizedTopologicalSimplification() {
      this->setDebugMsgPrefix("LTS");
    };
    ~LocalizedTopologicalSimplification(){};

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    };

    template <typename IT>
    int allocateMemory(std::vector<IT> &segmentation,
                       std::vector<IT> &queueMask,
                       std::vector<IT> &localOrder,
                       std::vector<SSCPropagation<IT> *> &propagationMask,
                       std::vector<std::tuple<IT, IT, IT>> &sortedIndices,

                       const IT &nVertices) const {
      ttk::Timer timer;

      // -------------------------------------------------------------
      // allocate and init memory
      // -------------------------------------------------------------
      this->printMsg("Allocating Memory", 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      segmentation.resize(nVertices);
      queueMask.resize(nVertices);
      localOrder.resize(nVertices);
      propagationMask.resize(nVertices);
      sortedIndices.resize(nVertices);

      this->printMsg(
        "Allocating Memory", 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    };

    template <typename IT>
    int initializeMemory(IT *segmentation,
                         IT *queueMask,
                         IT *localOrder,
                         SSCPropagation<IT> **propagationMask,

                         const IT &nVertices) const {
      ttk::Timer timer;

      // -------------------------------------------------------------
      // allocate and init memory
      // -------------------------------------------------------------
      this->printMsg("Initializing Memory", 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT i = 0; i < nVertices; i++) {
        segmentation[i] = -1;
        queueMask[i] = -1;
        localOrder[i] = -1;
        propagationMask[i] = nullptr;
      }

      this->printMsg(
        "Initializing Memory", 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    };

    template <typename IT, class TT>
    int initializePropagations(
      std::vector<SSCPropagation<IT>> &propagations,
      IT *authorizationMask, // assumes preservation mask is initialized as -1
      IT *maximaBuffer,

      const IT *authorizedExtremaIndices,
      const IT &nAuthorizedExtremaIndices,
      const IT *inputOrder,
      const TT *triangulation) const {

      ttk::Timer timer;
      this->printMsg("Initializing Propagations", 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      const IT nVertices = triangulation->getNumberOfVertices();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT i = 0; i < nAuthorizedExtremaIndices; i++)
        authorizationMask[authorizedExtremaIndices[i]] = -2;

      IT writeIndex = 0;
// find discareded maxima
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT v = 0; v < nVertices; v++) {

        // if v needs to be preserved then skip
        if(authorizationMask[v] == -2)
          continue;

        // check if v has larger neighbors
        bool hasLargerNeighbor = false;
        const IT &vOrder = inputOrder[v];
        const IT nNeighbors = triangulation->getVertexNeighborNumber(v);
        for(IT n = 0; n < nNeighbors; n++) {
          IT u;
          triangulation->getVertexNeighbor(v, n, u);
          if(vOrder < inputOrder[u]) {
            hasLargerNeighbor = true;
            break;
          }
        }

        // if v has larger neighbors then v can not be maximum
        if(hasLargerNeighbor)
          continue;

        // get local write index for this thread
        IT localWriteIndex = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif // TTK_ENABLE_OPENMP
        localWriteIndex = writeIndex++;

// atomic write to prevent cache line conflicts
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif // TTK_ENABLE_OPENMP
        maximaBuffer[localWriteIndex] = v;
      }

      // init propagations
      propagations.resize(writeIndex);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT p = 0; p < writeIndex; p++) {
        propagations[p].criticalPoints.push_back(maximaBuffer[p]);
      }

      this->printMsg("Initializing Propagations (" + std::to_string(writeIndex)
                       + "|" + std::to_string(nVertices) + ")",
                     1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    };

    /// This is a simple superlevel set propagation procedure that just absorbs
    /// the largest neighbor of the current set until the propagation encounters
    /// a saddle. To gain additional speedup this procedure does NOT compute a
    /// list of visited critical points and does not check for an early escape.
    template <typename IT, typename TT>
    int computeSimplePropagation(
      SSCPropagation<IT> &propagation,
      SSCPropagation<IT> **propagationMask,
      IT *segmentation, // used here to also store registered larger vertices
      IT *queueMask, // used to mark vertices that have already been added to
                     // the queue by this thread

      const TT *triangulation,
      const IT *orderField) const {

      // pointer used to compare against representative
      auto *currentPropagation = &propagation;

      // frequently used propagation members
      IT &extremumIndex = currentPropagation->criticalPoints[0];
      auto *queue = &currentPropagation->queue;

      // add extremumIndex to queue
      queue->emplace(orderField[extremumIndex], extremumIndex);

      IT queueMaskLabel = extremumIndex;
      queueMask[extremumIndex] = queueMaskLabel;
      IT segmentId = extremumIndex;

      // grow segment until prop reaches a saddle and then decide if prop should
      // continue
      IT v = -1;
      while(!queue->empty()) {
        v = std::get<1>(queue->top());
        queue->pop();

        // continue if this thread has already seen this vertex
        if(propagationMask[v] != nullptr)
          continue;

        const IT &orderV = orderField[v];

        // add neighbors to queue AND check if v is a saddle
        bool isSaddle = false;
        const IT nNeighbors = triangulation->getVertexNeighborNumber(v);

        IT numberOfLargerNeighbors = 0;
        IT numberOfLargerNeighborsThisThreadVisited = 0;
        for(IT n = 0; n < nNeighbors; n++) {
          IT u;
          triangulation->getVertexNeighbor(v, n, u);

          const IT &orderU = orderField[u];

          // if larger neighbor
          if(orderU > orderV) {
            numberOfLargerNeighbors++;

            auto uPropagation = propagationMask[u];
            if(uPropagation == nullptr
               || currentPropagation != uPropagation->find())
              isSaddle = true;
            else
              numberOfLargerNeighborsThisThreadVisited++;
          } else if(queueMask[u] != queueMaskLabel) {
            queue->emplace(orderU, u);
            queueMask[u] = queueMaskLabel;
          }
        }

        // if v is a saddle we have to check if the current thread is the last
        // visitor
        if(isSaddle) {
          currentPropagation->lastEncounteredCriticalPoint = v;

          IT numberOfRegisteredLargerVertices = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif // TTK_ENABLE_OPENMP
          {
            segmentation[v] -= numberOfLargerNeighborsThisThreadVisited;
            numberOfRegisteredLargerVertices = segmentation[v];
          }

          // if this thread did not register the last remaining larger vertices
          // then terminate propagation
          if(numberOfRegisteredLargerVertices != -numberOfLargerNeighbors - 1)
            return 1;

          // merge propagations
          for(IT n = 0; n < nNeighbors; n++) {
            IT u;
            triangulation->getVertexNeighbor(v, n, u);

            if(propagationMask[u] != nullptr) {
              auto uPropagation = propagationMask[u]->find();
              if(uPropagation != currentPropagation) {
                currentPropagation = SSCPropagation<IT>::unify(
                  currentPropagation, uPropagation, orderField);
              }
            }
          }

          queue = &currentPropagation->queue;
          segmentId = v;
          queueMaskLabel = currentPropagation->criticalPoints[0];
        }

        // mark vertex as visited and continue
        propagationMask[v] = currentPropagation;
        segmentation[v] = segmentId;
        currentPropagation->segmentSize++;
      }

      this->printErr(
        "Simple propagations should never reach global minimum/maximum.");

      return 0;
    };

    template <typename IT, class TT>
    int computeSimplePropagations(std::vector<SSCPropagation<IT>> &propagations,
                                  SSCPropagation<IT> **propagationMask,
                                  IT *segmentation,
                                  IT *queueMask,

                                  const TT *triangulation,
                                  const IT *inputOrder) const {
      ttk::Timer timer;

      const IT nPropagations = propagations.size();
      this->printMsg(
        "Computing Propagations (" + std::to_string(nPropagations) + ")", 0, 0,
        this->threadNumber_, debug::LineMode::REPLACE);

      int status = 1;
// compute propagations
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT p = 0; p < nPropagations; p++) {
        int localStatus = this->computeSimplePropagation<IT, TT>(
          propagations[p], propagationMask, segmentation, queueMask,

          triangulation, inputOrder);

        if(!localStatus)
          status = 0;
      }

      if(!status)
        return 0;

      this->printMsg(
        "Computing Propagations (" + std::to_string(nPropagations) + ")", 1,
        timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    template <typename IT>
    int finalizePropagations(
      std::vector<SSCPropagation<IT> *> &parentPropagations,
      std::vector<SSCPropagation<IT>> &propagations,

      const IT nVertices) const {
      ttk::Timer timer;

      const IT nPropagations = propagations.size();

      this->printMsg(
        "Finalizing Propagations (" + std::to_string(nPropagations) + ")", 0,
        timer.getElapsedTime(), this->threadNumber_, debug::LineMode::REPLACE);

      IT nSegmentVertices = 0;
      IT nParentPropagations = 0;
      parentPropagations.clear();
      parentPropagations.resize(nPropagations);
      for(IT p = 0; p < nPropagations; p++) {
        auto *propagation = &propagations[p];
        if(propagation->parent == propagation) {
          nSegmentVertices = nSegmentVertices + propagation->segmentSize;
          parentPropagations[nParentPropagations++] = propagation;
        }
      }
      parentPropagations.resize(nParentPropagations);

      std::stringstream pFraction, vFraction;
      pFraction << std::fixed << std::setprecision(2)
                << ((float)nParentPropagations / (float)nPropagations);
      vFraction << std::fixed << std::setprecision(2)
                << ((float)nSegmentVertices / (float)nVertices);

      this->printMsg("Finalizing Propagations ("
                       + std::to_string(nParentPropagations) + "|"
                       + toFixed(nParentPropagations, nPropagations) + "|"
                       + toFixed(nSegmentVertices, nVertices) + ")",
                     1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    };

    template <typename IT, class TT>
    int computeSegment(IT *segmentation,
                       SSCPropagation<IT> *propagation,

                       const TT *triangulation) const {

      const IT &extremumIndex = propagation->criticalPoints[0];

      // collect segment
      auto &segment = propagation->segment;

      segment.resize(propagation->segmentSize);
      IT segmentIndex = 0;
      {
        std::vector<IT> queue(propagation->segmentSize);
        IT queueIndex = 0;

        // init queue
        {
          queue[queueIndex++] = extremumIndex;
          segmentation[extremumIndex] = -1000;
        }

        // flood fill by starting from extremum
        while(queueIndex > 0) {
          const IT v = queue[--queueIndex];

          segment[segmentIndex++] = v;

          IT nNeighbors = triangulation->getVertexNeighborNumber(v);
          for(IT n = 0; n < nNeighbors; n++) {
            IT u;
            triangulation->getVertexNeighbor(v, n, u);
            if(segmentation[u] >= 0) {
              segmentation[u] = -999; // mark as visited
              queue[queueIndex++] = u; // add to queue
            }
          }
        }
      }

      if(segmentIndex != propagation->segmentSize) {
        this->printErr("Segment size incorrect: " + std::to_string(segmentIndex)
                       + " " + std::to_string(propagation->segmentSize));
        return 0;
      }

      for(auto idx : propagation->segment)
        segmentation[idx] = extremumIndex;

      return 1;
    };

    template <typename IT, class TT>
    int computeSegments(IT *segmentation,
                        std::vector<SSCPropagation<IT> *> &propagations,

                        const TT *triangulation) const {

      const IT nPropagations = propagations.size();
      const IT nVertices = triangulation->getNumberOfVertices();

      ttk::Timer timer;
      this->printMsg(
        "Computing segments (" + std::to_string(nPropagations) + ")", 0, 0,
        this->threadNumber_, debug::LineMode::REPLACE);

      int status = 1;

// compute segments in parallel
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT p = 0; p < nPropagations; p++) {
        int localStatus
          = this->computeSegment<IT, TT>(segmentation, propagations[p],

                                         triangulation);
        if(!localStatus)
          status = 0;
      }
      if(!status)
        return 0;

      // print status
      if(this->debugLevel_ < 4) {
        this->printMsg(
          "Computing segments (" + std::to_string(nPropagations) + ")", 1,
          timer.getElapsedTime(), this->threadNumber_);
      } else {

        IT min = propagations[0]->segmentSize;
        IT max = min;
        IT avg = 0;

        for(IT p = 0; p < nPropagations; p++) {
          const auto propagation = propagations[p];
          if(min > propagation->segmentSize)
            min = propagation->segmentSize;
          if(max < propagation->segmentSize)
            max = propagation->segmentSize;
          avg += propagation->segmentSize;
        }

        avg /= nPropagations;

        this->printMsg("Computing segments (" + std::to_string(nPropagations)
                         + "|" + toFixed(min, nVertices) + "|"
                         + toFixed(avg, nVertices) + "|"
                         + toFixed(max, nVertices) + ")",
                       1, timer.getElapsedTime(), this->threadNumber_);
      }

      return 1;
    };

    template <typename IT, class TT>
    int computeLocalOrderOfSegmentIteration(
      IT *localOrder,
      IT *localVertexSequence,

      const bool &performSuperlevelSetPropagation,
      const TT *triangulation,
      const IT *segmentation,
      const IT &segmentId,
      const std::vector<IT> &boundary,
      const std::vector<IT> &segment,
      const IT &saddleIdx) const {
      // init priority queue
      std::priority_queue<std::pair<IT, IT>, std::vector<std::pair<IT, IT>>>
        queue;

      IT nSegmentVertices = segment.size();

      if(performSuperlevelSetPropagation) {
        // add saddle to queue
        queue.emplace(std::numeric_limits<IT>::max(), saddleIdx);
      } else {
        // invert order
        IT offset = -nSegmentVertices - 1; // ensures that inverted order < 0
        for(IT i = 0; i < nSegmentVertices; i++)
          localOrder[segment[i]] = offset - localOrder[segment[i]];

        // add all boundary vertices to queue
        for(const auto &v : boundary) {
          queue.emplace(localOrder[v], v);
          localOrder[v] = 0;
        }

        // also add saddle to queue
        queue.emplace(std::numeric_limits<IT>::min(), saddleIdx);
      }

      IT q = 0;
      while(!queue.empty()) {
        IT v = std::get<1>(queue.top());
        queue.pop();

        localVertexSequence[q++] = v;

        IT nNeighbors = triangulation->getVertexNeighborNumber(v);
        for(IT n = 0; n < nNeighbors; n++) {
          IT u;
          triangulation->getVertexNeighbor(v, n, u);

          // if u is in segment and has not already been added to the queue
          if(segmentation[u] == segmentId && localOrder[u] < 0) {
            queue.emplace(localOrder[u], u);
            localOrder[u] = 0;
          }
        }
      }

      if(performSuperlevelSetPropagation) {
        // skip first idx as it correponds to saddle
        for(IT i = 1; i <= nSegmentVertices; i++)
          localOrder[localVertexSequence[i]] = -i;
      } else {
        IT order = -nSegmentVertices;
        // skip last idx as it correponds to saddle
        for(IT i = 0; i < nSegmentVertices; i++)
          localOrder[localVertexSequence[i]] = order++;
      }

      return 1;
    };

    template <typename IT, class TT>
    int computeLocalOrderOfSegment(IT *localOrder,

                                   const SSCPropagation<IT> *propagation,
                                   const TT *triangulation,
                                   const IT *segmentation,
                                   const IT *inputOrder) const {
      // quick espace for small segments
      if(propagation->segmentSize == 1) {
        localOrder[propagation->segment[0]] = -2;
        return 1;
      }

      // init local order by input order
      // Note: force orders to be smaller than 0 for in situ updates in
      // propagation procedures
      IT nVertices = triangulation->getNumberOfVertices();
      for(const auto &v : propagation->segment)
        localOrder[v] = inputOrder[v] - nVertices;

      const IT &extremumIndex = propagation->criticalPoints[0];
      const IT &saddleIndex = propagation->lastEncounteredCriticalPoint;

      // vector to record boundary
      std::vector<IT> boundary;

      // make enough room for segment + saddle
      std::vector<IT> localVertexSequence(propagation->segmentSize + 1);

      int status = 1;
      bool containsResidualExtrema = true;
      bool performSuperlevelSetPropagation = true;
      while(containsResidualExtrema) {
        propagation->nIterations++;

        if(this->debugLevel_ > 3 && propagation->nIterations == 20)
          this->printWrn("Unable to converge with less than 20 iterations!");

        // execute superlevel set propagation
        status = this->computeLocalOrderOfSegmentIteration<IT, TT>(
          localOrder, localVertexSequence.data(),

          performSuperlevelSetPropagation, triangulation, segmentation,
          extremumIndex, boundary, propagation->segment, saddleIndex);
        if(!status)
          return 0;

        performSuperlevelSetPropagation = !performSuperlevelSetPropagation;

        IT boundaryWriteIdx = 0;
        IT nResidualMaxima = 0;
        IT nResidualMinima = 0;

        for(const auto &v : propagation->segment) {
          bool isOnSegmentBoundary = false;
          bool hasSmallerNeighbor = false;
          bool hasLargerNeighbor = false;

          const auto &vOrder = localOrder[v];

          IT nNeighbors = triangulation->getVertexNeighborNumber(v);
          for(IT n = 0; n < nNeighbors; n++) {
            IT u;
            triangulation->getVertexNeighbor(v, n, u);

            if(u == saddleIndex) {
              hasLargerNeighbor = true;
              continue;
            }

            // if u is not inside segment -> v is on segment boundary
            if(segmentation[u] != extremumIndex) {
              isOnSegmentBoundary = true;
            } else {
              if(localOrder[u] > vOrder)
                hasLargerNeighbor = true;
              else
                hasSmallerNeighbor = true;
            }
          }

          if(!hasLargerNeighbor)
            nResidualMaxima++;

          if(isOnSegmentBoundary) {
            localVertexSequence[boundaryWriteIdx++] = v;
          } else if(!hasSmallerNeighbor) {
            nResidualMinima++;
          }
        }

        containsResidualExtrema = nResidualMinima > 0 || nResidualMaxima > 0;

        if(containsResidualExtrema && boundary.size() == 0) {
          boundary.resize(boundaryWriteIdx);
          for(IT i = 0; i < boundaryWriteIdx; i++) {
            boundary[i] = localVertexSequence[i];
          }
        }
      }

      return 1;
    };

    template <typename IT, class TT>
    int computeLocalOrderOfSegments(
      IT *localOrder,

      const TT *triangulation,
      const IT *segmentation,
      const IT *inputOrder,
      const std::vector<SSCPropagation<IT> *> &propagations) const {
      ttk::Timer timer;

      const IT nPropagations = propagations.size();
      this->printMsg("Computing Local Order of Segments ("
                       + std::to_string(nPropagations) + ")",
                     0, 0, this->threadNumber_, debug::LineMode::REPLACE);

      int status = 1;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT p = 0; p < nPropagations; p++) {
        int localStatus = this->computeLocalOrderOfSegment<IT>(
          localOrder,

          propagations[p], triangulation, segmentation, inputOrder);
        if(!localStatus)
          status = 0;
      }
      if(!status)
        return 0;

// enforce that saddles have the highest local order
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT p = 0; p < nPropagations; p++) {
        localOrder[propagations[p]->lastEncounteredCriticalPoint]
          = std::numeric_limits<IT>::max();
      }

      if(this->debugLevel_ < 4) {
        this->printMsg("Computing Local Order of Segments ("
                         + std::to_string(nPropagations) + ")",
                       1, timer.getElapsedTime(), this->threadNumber_);
      } else {
        IT min = propagations[0]->nIterations;
        IT max = min;
        IT avg = 0;

        for(IT p = 0; p < nPropagations; p++) {
          const auto propagation = propagations[p];
          if(min > propagation->nIterations)
            min = propagation->nIterations;
          if(max < propagation->nIterations)
            max = propagation->nIterations;
          avg += propagation->nIterations;
        }

        avg /= nPropagations;

        this->printMsg("Computing Local Order of Segments ("
                         + std::to_string(nPropagations) + "|"
                         + std::to_string(min) + "|" + std::to_string(avg) + "|"
                         + std::to_string(max) + ")",
                       1, timer.getElapsedTime(), this->threadNumber_);
      }

      return 1;
    };

    template <typename IT>
    int
      flattenOrder(IT *outputOrder,

                   const std::vector<SSCPropagation<IT> *> &parentPropagations,
                   const IT &nVertices) const {
      ttk::Timer timer;
      this->printMsg("Flattening Order Field", 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      const IT nParentPropagations = parentPropagations.size();

// flatten segemnt to order of last encountered saddles
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT p = 0; p < nParentPropagations; p++) {
        const auto *propagation = parentPropagations[p];
        const auto &order
          = outputOrder[propagation->lastEncounteredCriticalPoint];
        for(const auto &v : propagation->segment)
          outputOrder[v] = order;
      }

      this->printMsg("Flattening Order Field", 1, timer.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    };

    template <typename DT, typename IT>
    int flattenScalars(DT *outputScalars,

                       const std::vector<SSCPropagation<IT>> &propagationsMin,
                       const std::vector<SSCPropagation<IT>> &propagationsMax,
                       const DT *inputScalars,
                       const IT &nVertices) const {
      ttk::Timer timer;
      this->printMsg("Flattening Scalar Field", 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT v = 0; v < nVertices; v++)
        outputScalars[v] = inputScalars[v];

      std::vector<const std::vector<SSCPropagation<IT>> *> propagationsPair
        = {&propagationsMin, &propagationsMax};

      for(const auto propagations : propagationsPair) {
        const IT nPropagations = propagations->size();

// flatten scalars to saddle value
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT p = 0; p < nPropagations; p++) {
          const auto &propagation = (*propagations)[p];
          if(propagation.parent != &propagation)
            continue;

          const auto &saddleScalar
            = outputScalars[propagation.lastEncounteredCriticalPoint];
          for(const auto &v : propagation.segment)
            outputScalars[v] = saddleScalar;
        }
      }

      this->printMsg("Flattening Scalar Field", 1, timer.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    };

    template <typename DT, typename IT>
    int computeGlobalOrder(IT *outputOrder,
                           std::vector<std::tuple<DT, IT, IT>> &sortedIndices,

                           const DT *rank1,
                           const IT *rank2) const {
      ttk::Timer timer;

      const IT nVertices = sortedIndices.size();

// init tuples
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT i = 0; i < nVertices; i++) {
        auto &t = sortedIndices[i];
        std::get<0>(t) = rank1[i];
        std::get<1>(t) = rank2[i];
        std::get<2>(t) = i;
      }

      this->printMsg("Computing Global Order", 0.2, timer.getElapsedTime(),
                     this->threadNumber_, debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#ifdef __clang__
      this->printWrn("Caution, outside GCC, sequential sort");
      std::sort(sortedIndices.begin(), sortedIndices.end());
#else
      omp_set_num_threads(this->threadNumber_);
      __gnu_parallel::sort(sortedIndices.begin(), sortedIndices.end());
      omp_set_num_threads(1);
#endif
#else
      this->printWrn("Caution, outside GCC, sequential sort");
      std::sort(sortedIndices.begin(), sortedIndices.end());
#endif

      this->printMsg("Computing Global Order", 0.8, timer.getElapsedTime(),
                     this->threadNumber_, debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT i = 0; i < nVertices; i++)
        outputOrder[std::get<2>(sortedIndices[i])] = i;

      this->printMsg("Computing Global Order", 1, timer.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    };

    template <typename DT, typename IT>
    int computeNumericalPerturbation(
      DT *outputScalars,
      const std::vector<std::tuple<IT, IT, IT>> &sortedIndices) const {
      ttk::Timer timer;
      this->printMsg("Applying numerical perturbation", 0, 0,
                     this->threadNumber_, debug::LineMode::REPLACE);

      const IT nVertices = sortedIndices.size();
      for(IT i = 1; i < nVertices; i++) {
        const IT &v0 = std::get<2>(sortedIndices[i - 1]);
        const IT &v1 = std::get<2>(sortedIndices[i]);
        if(outputScalars[v0] >= outputScalars[v1])
          outputScalars[v1] = boost::math::float_next(outputScalars[v0]);
      }

      this->printMsg("Applying numerical perturbation", 1,
                     timer.getElapsedTime(), this->threadNumber_);

      return 1;
    };

    template <typename IT, class TT>
    int detectAndRemoveUnauthorizedMaxima(
      IT *outputOrder,
      IT *segmentation,
      IT *queueMask,
      IT *localOrder,
      SSCPropagation<IT> **propagationMask,
      std::vector<SSCPropagation<IT>> &propagations,
      std::vector<std::tuple<IT, IT, IT>> &sortedIndices,

      const TT *triangulation,
      const IT *authorizedExtremaIndices,
      const IT &nAuthorizedExtremaIndices) const {
      const IT nVertices = triangulation->getNumberOfVertices();

      int status = 0;

      // init propagations
      status = this->initializeMemory<IT>(segmentation, queueMask, localOrder,
                                          propagationMask,

                                          nVertices);
      if(!status)
        return 0;

      // init propagations
      status = this->initializePropagations<IT, TT>(
        propagations,
        queueMask, // use as authorization mask (will be overriden by subsequent
                   // procedures)
        localOrder, // use as maxima buffer (will be overriden by subsequent
                    // procedures)

        authorizedExtremaIndices, nAuthorizedExtremaIndices, outputOrder,
        triangulation);
      if(!status)
        return 0;

      // compute propagations
      status = this->computeSimplePropagations<IT>(
        propagations, propagationMask, segmentation, queueMask,

        triangulation, outputOrder);
      if(!status)
        return 0;

      // finalize master propagations
      std::vector<SSCPropagation<IT> *> parentPropagations;
      status = this->finalizePropagations<IT>(parentPropagations, propagations,

                                              nVertices);
      if(!status)
        return 0;

      // compute segments
      status = this->computeSegments<IT, TT>(segmentation, parentPropagations,

                                             triangulation);
      if(!status)
        return 0;

      // compute local order of segments
      status = this->computeLocalOrderOfSegments<IT, TT>(
        localOrder,

        triangulation, segmentation, outputOrder, parentPropagations);
      if(!status)
        return 0;

      // flatten order
      status = this->flattenOrder<IT>(outputOrder,

                                      parentPropagations, nVertices);
      if(!status)
        return 0;

      // compute global offsets
      status = this->computeGlobalOrder<IT, IT>(outputOrder, sortedIndices,

                                                outputOrder, localOrder);
      if(!status)
        return 0;

      return 1;
    };

    template <typename IT>
    int invertArray(IT *outputOrder,

                    const IT *inputOrder,
                    const IT &nVertices) const {
      ttk::Timer timer;
      this->printMsg(
        "Inverting Array", 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT v = 0; v < nVertices; v++)
        outputOrder[v] = nVertices - inputOrder[v];

      this->printMsg(
        "Inverting Array", 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    template <typename DT, typename IT, class TT>
    int removeUnauthorizedExtrema(DT *outputScalars,
                                  IT *outputOrder,

                                  const TT *triangulation,
                                  const DT *inputScalars,
                                  const IT *inputOrder,
                                  const IT *authorizedExtremaIndices,
                                  const IT &nAuthorizedExtremaIndices,
                                  const bool &addPerturbation) const {
      ttk::Timer globalTimer;

      IT nVertices = triangulation->getNumberOfVertices();

      // Allocating Memory
      int status = 0;
      std::vector<IT> segmentation;
      std::vector<IT> queueMask;
      std::vector<IT> localOrder;
      std::vector<SSCPropagation<IT> *> propagationMask;
      std::vector<SSCPropagation<IT>> propagationsMax;
      std::vector<SSCPropagation<IT>> propagationsMin;
      std::vector<std::tuple<IT, IT, IT>> sortedIndices;

      this->allocateMemory<IT>(segmentation, queueMask, localOrder,
                               propagationMask, sortedIndices,

                               nVertices);

      {
        this->printMsg("----------- [Removing unauthorized minima]",
                       ttk::debug::Separator::L2);

        // init output order as input order
        status = this->invertArray(outputOrder,

                                   inputOrder, nVertices);
        if(!status)
          return 0;

        status = this->detectAndRemoveUnauthorizedMaxima<IT, TT>(
          outputOrder, segmentation.data(), queueMask.data(), localOrder.data(),
          propagationMask.data(), propagationsMin, sortedIndices,

          triangulation, authorizedExtremaIndices, nAuthorizedExtremaIndices);
        if(!status)
          return 0;
      }

      // Maxima
      {
        this->printMsg("----------- [Removing unauthorized maxima]",
                       ttk::debug::Separator::L2);

        // invert outputOrder
        status = this->invertArray(outputOrder,

                                   outputOrder, nVertices);
        if(!status)
          return 0;

        status = this->detectAndRemoveUnauthorizedMaxima<IT, TT>(
          outputOrder, segmentation.data(), queueMask.data(), localOrder.data(),
          propagationMask.data(), propagationsMax, sortedIndices,

          triangulation, authorizedExtremaIndices, nAuthorizedExtremaIndices);
        if(!status)
          return 0;
      }

      status = this->flattenScalars<DT, IT>(outputScalars,

                                            propagationsMin, propagationsMax,
                                            inputScalars, nVertices);
      if(!status)
        return 0;

      // optionally add perturbation
      if(addPerturbation) {
        this->printMsg(debug::Separator::L2);
        status = this->computeNumericalPerturbation<DT, IT>(outputScalars,

                                                            sortedIndices);
        if(!status)
          return 0;
      }

      this->printMsg(debug::Separator::L2);
      this->printMsg(
        "Complete", 1, globalTimer.getElapsedTime(), this->threadNumber_);

      this->printMsg(debug::Separator::L1);

      return 1;
    };

  }; // class
} // namespace ttk
