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

// for numerical perturbation
#include <boost/math/special_functions/next.hpp>

#include <Debug.h>
#include <Triangulation.h>

#include <Propagation.h>

namespace ttk {
  namespace lts {

    inline std::string toFixed(const float &number, const int precision = 2) {
      std::stringstream vFraction;
      vFraction << std::fixed << std::setprecision(precision) << number;
      return vFraction.str();
    }

    template <typename IT>
    std::string
      toFixed(const IT &number0, const IT &number1, const int precision = 2) {
      return toFixed(((float)number0) / ((float)number1), precision);
    }

    class LocalizedTopologicalSimplification : virtual public Debug {

    public:
      enum class PAIR_TYPE {
        EXTREMUM_SADDLE = 0,
        MINIMUM_SADDLE = 1,
        MAXIMUM_SADDLE = 2
      };

      LocalizedTopologicalSimplification() {
        this->setDebugMsgPrefix("LTS");
      }
      ~LocalizedTopologicalSimplification() = default;

      int preconditionTriangulation(
        ttk::AbstractTriangulation *triangulation) const {
        return triangulation->preconditionVertexNeighbors();
      }

      /// This method allocates all temporary memory required for LTS
      /// procedures.
      template <typename IT>
      int allocateMemory(std::vector<IT> &segmentation,
                         std::vector<IT> &queueMask,
                         std::vector<IT> &localOrder,
                         std::vector<Propagation<IT> *> &propagationMask,
                         std::vector<std::tuple<IT, IT, IT>> &sortedIndices,

                         const IT &nVertices) const {
        ttk::Timer timer;

        constexpr char msg[] = "Allocating Memory";
        this->printMsg(
          msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

        segmentation.resize(nVertices);
        queueMask.resize(nVertices);
        localOrder.resize(nVertices);
        propagationMask.resize(nVertices);
        sortedIndices.resize(nVertices);

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

        return 1;
      }

      /// This method initializes all temporary memory for LTS procedures
      /// (assumes memory is already allocated).
      template <typename IT>
      int initializeMemory(IT *segmentation,
                           IT *queueMask,
                           IT *localOrder,
                           Propagation<IT> **propagationMask,

                           const IT &nVertices) const {
        ttk::Timer timer;

        constexpr char msg[] = "Initializing Memory";
        this->printMsg(
          msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT i = 0; i < nVertices; i++) {
          segmentation[i] = -1;
          queueMask[i] = -1;
          localOrder[i] = -1;
          propagationMask[i] = nullptr;
        }

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

        return 1;
      }

      /// This method iterates over an order array and detects all maxima, for
      /// which it will then initialize a propagation data structure. One can
      /// also pass a whitelist of authorized maxima that will be skipped. The
      /// maximaBuffer parameter needs to be big enough to fit all possible
      /// maxima indices (the buffer is explicitly exposed for memory reuse).
      /// Note, propagations are sorted in ascending order of their
      /// corresponding maxima.
      template <typename IT, class TT>
      int initializePropagations(
        std::vector<Propagation<IT>> &propagations,
        IT *authorizationMask, // assumes preservation mask is initialized as -1
        IT *maximaBuffer,

        const IT *authorizedExtremaIndices,
        const IT &nAuthorizedExtremaIndices,
        const IT *order,
        const TT *triangulation) const {

        ttk::Timer timer;
        this->printMsg("Initializing Propagations", 0, 0, this->threadNumber_,
                       debug::LineMode::REPLACE);

        const IT nVertices = triangulation->getNumberOfVertices();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads( \
  this->threadNumber_) if(nAuthorizedExtremaIndices > 1000)
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
          const IT &vOrder = order[v];
          const IT nNeighbors = triangulation->getVertexNeighborNumber(v);
          for(IT n = 0; n < nNeighbors; n++) {
            IT u;
            triangulation->getVertexNeighbor(v, n, u);
            if(vOrder < order[u]) {
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

        // sort maxima in ascending order
        std::sort(maximaBuffer, maximaBuffer + writeIndex,
                  [=](const IT &a, const IT &b) -> bool {
                    return order[a] < order[b];
                  });

        // if there are no authorized maxima then always skip most persistent
        // prop
        if(nAuthorizedExtremaIndices < 1)
          writeIndex--;

        // init propagations
        propagations.resize(writeIndex);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT p = 0; p < writeIndex; p++) {
          propagations[p].criticalPoints.push_back(maximaBuffer[p]);
        }

        this->printMsg("Initializing Propagations ("
                         + std::to_string(writeIndex) + "|"
                         + std::to_string(nVertices) + ")",
                       1, timer.getElapsedTime(), this->threadNumber_);

        return 1;
      }

      /// This is a simple superlevel set propagation procedure that just
      /// absorbs the largest neighbor of the current set until the propagation
      /// encounters a saddle. To gain additional speedup this procedure keeps
      /// track of vertices that have already been added to the Fibonacci heap
      /// via the queueMask (duplicate entries slow down the heap).
      template <typename IT, typename TT>
      int computeSimplePropagation(
        Propagation<IT> &propagation,
        Propagation<IT> **propagationMask,
        IT *segmentation, // used here to also store registered larger vertices
        IT *queueMask, // used to mark vertices that have already been added to
                       // the queue by this thread

        const TT *triangulation,
        const IT *order) const {

        // pointer used to compare against representative
        auto *currentPropagation = &propagation;

        // add extremumIndex (stored in segmentId) to queue
        IT segmentId = currentPropagation->criticalPoints[0];
        auto queue = &currentPropagation->queue;
        queue->emplace(order[segmentId], segmentId);

        queueMask[segmentId] = segmentId;

        // grow segment until prop reaches a saddle and then decide if prop
        // should continue
        IT v = -1;
        while(!queue->empty()) {
          v = std::get<1>(queue->top());
          queue->pop();

          // continue if this thread has already seen this vertex
          if(propagationMask[v] != nullptr)
            continue;

          const IT &orderV = order[v];

          // add neighbors to queue AND check if v is a saddle
          bool isSaddle = false;
          const IT nNeighbors = triangulation->getVertexNeighborNumber(v);

          IT numberOfLargerNeighbors = 0;
          IT numberOfLargerNeighborsThisThreadVisited = 0;
          for(IT n = 0; n < nNeighbors; n++) {
            IT u;
            triangulation->getVertexNeighbor(v, n, u);

            const IT &orderU = order[u];

            // if larger neighbor
            if(orderU > orderV) {
              numberOfLargerNeighbors++;

              auto uPropagation = propagationMask[u];
              if(uPropagation == nullptr
                 || currentPropagation != uPropagation->find())
                isSaddle = true;
              else
                numberOfLargerNeighborsThisThreadVisited++;
            } else if(queueMask[u] != segmentId) {
              queue->emplace(orderU, u);
              queueMask[u] = segmentId;
            }
          }

          // if v is a saddle we have to check if the current thread is the last
          // visitor
          if(isSaddle) {
            currentPropagation->criticalPoints.push_back(v);

            IT numberOfRegisteredLargerVertices = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif // TTK_ENABLE_OPENMP
            {
              segmentation[v] -= numberOfLargerNeighborsThisThreadVisited;
              numberOfRegisteredLargerVertices = segmentation[v];
            }

            // if this thread did not register the last remaining larger
            // vertices then terminate propagation
            if(numberOfRegisteredLargerVertices != -numberOfLargerNeighbors - 1)
              return 1;

            // get most dominant propagation
            std::vector<Propagation<IT> *> neighborPropagations(
              nNeighbors, nullptr);
            for(IT n = 0; n < nNeighbors; n++) {
              IT u;
              triangulation->getVertexNeighbor(v, n, u);
              if(propagationMask[u] != nullptr) {
                auto &neighborPropagation = neighborPropagations[n];
                neighborPropagation = propagationMask[u]->find();
                if(order[neighborPropagation->criticalPoints[0]]
                   > order[currentPropagation->criticalPoints[0]])
                  currentPropagation = neighborPropagation;
              }
            }

            // merge propagations
            for(IT n = 0; n < nNeighbors; n++) {
              Propagation<IT>::unify(
                currentPropagation, neighborPropagations[n]);
            }

            queue = &currentPropagation->queue;
            segmentId = currentPropagation->criticalPoints[0];
          }

          // mark vertex as visited and continue
          propagationMask[v] = currentPropagation;
          segmentation[v] = segmentId;
          currentPropagation->segmentSize++;
        }

        this->printErr(
          "Simple propagations should never reach global minimum/maximum.");

        return 0;
      }

      /// Basically the same as the simple propagation procedure, except that a
      /// propagation keeps track of the persistence of the computed
      /// propagation. As soon as the persistence strictly exceeds the given
      /// threshold the propagation is aborted.
      template <typename IT, typename DT, typename TT>
      int computePersistenceSensitivePropagation(
        Propagation<IT> &propagation,
        Propagation<IT> **propagationMask,
        IT *segmentation, // used here to also store registered larger vertices
        IT *queueMask, // used to mark vertices that have already been added to
                       // the queue by this thread

        const TT *triangulation,
        const IT *order,
        const DT *scalars,
        const DT persistenceThreshold) const {

        // pointer used to compare against representative
        auto *currentPropagation = &propagation;

        // add extremumIndex (stored in segmentId) to queue
        IT segmentId = currentPropagation->criticalPoints[0];
        auto queue = &currentPropagation->queue;
        queue->emplace(order[segmentId], segmentId);

        DT s0 = scalars[segmentId];

        queueMask[segmentId] = segmentId;

        // grow segment until prop reaches a saddle and then decide if prop
        // should continue
        IT v = -1;
        while(!queue->empty()) {
          v = std::get<1>(queue->top());
          queue->pop();

          // continue if this thread has already seen this vertex
          if(propagationMask[v] != nullptr)
            continue;

          // if propagation exceeds persistence threshold stop
          const DT &s1 = scalars[v];
          const DT &sd = s0 < s1 ? s1 - s0 : s0 - s1;
          if(sd > persistenceThreshold) {
            currentPropagation->aborted = true;
            return 1;
          }

          const IT &orderV = order[v];

          // add neighbors to queue AND check if v is a saddle
          bool isSaddle = false;
          const IT nNeighbors = triangulation->getVertexNeighborNumber(v);

          IT numberOfLargerNeighbors = 0;
          IT numberOfLargerNeighborsThisThreadVisited = 0;
          for(IT n = 0; n < nNeighbors; n++) {
            IT u;
            triangulation->getVertexNeighbor(v, n, u);

            const IT &orderU = order[u];

            // if larger neighbor
            if(orderU > orderV) {
              numberOfLargerNeighbors++;

              auto uPropagation = propagationMask[u];
              if(uPropagation == nullptr
                 || currentPropagation != uPropagation->find())
                isSaddle = true;
              else
                numberOfLargerNeighborsThisThreadVisited++;
            } else if(queueMask[u] != segmentId) {
              queue->emplace(orderU, u);
              queueMask[u] = segmentId;
            }
          }

          // if v is a saddle we have to check if the current thread is the last
          // visitor
          if(isSaddle) {
            currentPropagation->criticalPoints.push_back(v);

            IT numberOfRegisteredLargerVertices = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif // TTK_ENABLE_OPENMP
            {
              segmentation[v] -= numberOfLargerNeighborsThisThreadVisited;
              numberOfRegisteredLargerVertices = segmentation[v];
            }

            // if this thread did not register the last remaining larger
            // vertices then terminate propagation
            if(numberOfRegisteredLargerVertices != -numberOfLargerNeighbors - 1)
              return 1;

            // get most dominant propagation
            std::vector<Propagation<IT> *> neighborPropagations(
              nNeighbors, nullptr);
            for(IT n = 0; n < nNeighbors; n++) {
              IT u;
              triangulation->getVertexNeighbor(v, n, u);
              if(propagationMask[u] != nullptr) {
                auto &neighborPropagation = neighborPropagations[n];
                neighborPropagation = propagationMask[u]->find();
                if(order[neighborPropagation->criticalPoints[0]]
                   > order[currentPropagation->criticalPoints[0]])
                  currentPropagation = neighborPropagation;
              }
            }

            // merge propagations
            for(IT n = 0; n < nNeighbors; n++) {
              Propagation<IT>::unify(
                currentPropagation, neighborPropagations[n]);
            }

            queue = &currentPropagation->queue;
            segmentId = currentPropagation->criticalPoints[0];
            s0 = scalars[segmentId];
          }

          // mark vertex as visited and continue
          propagationMask[v] = currentPropagation;
          segmentation[v] = segmentId;
          currentPropagation->segmentSize++;
        }

        return 0;
      }

      /// This method computes (optionally in parallel) a list of simple
      /// propagations.
      template <typename IT, class TT>
      int computeSimplePropagations(std::vector<Propagation<IT>> &propagations,
                                    Propagation<IT> **propagationMask,
                                    IT *segmentation,
                                    IT *queueMask,

                                    const TT *triangulation,
                                    const IT *inputOrder) const {
        ttk::Timer timer;
        const IT nPropagations = propagations.size();
        const std::string msg
          = "Computing Propagations (" + std::to_string(nPropagations) + ")";
        this->printMsg(
          msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

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

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

        return 1;
      }

      /// This method computes (optionally in parallel) a list of
      /// persistence-sensitive propagations.
      template <typename IT, typename DT, class TT>
      int computePersistenceSensitivePropagations(
        std::vector<Propagation<IT>> &propagations,
        Propagation<IT> **propagationMask,
        IT *segmentation,
        IT *queueMask,

        const TT *triangulation,
        const IT *order,
        const DT *scalars,
        const DT persistenceThreshold) const {
        ttk::Timer timer;
        const IT nPropagations = propagations.size();
        const std::string msg = "Computing Persistence-Sensitive Propagations ("
                                + std::to_string(nPropagations) + ")";
        this->printMsg(
          msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

        int status = 1;
// compute propagations
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT p = 0; p < nPropagations; p++) {
          int localStatus
            = this->computePersistenceSensitivePropagation<IT, DT, TT>(
              propagations[p], propagationMask, segmentation, queueMask,

              triangulation, order, scalars, persistenceThreshold);

          if(!localStatus)
            status = 0;
        }

        if(!status)
          return 0;

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

        return 1;
      }

      /// This method identifies from a set of propagations so-called parent
      /// propagations, which are those propagations that either terminated at a
      /// saddle who has larger unvisited neighbors, or that were merged into a
      /// propagation that was aborted (during persistence-sensitive
      /// simplification).
      template <typename IT>
      int
        finalizePropagations(std::vector<Propagation<IT> *> &parentPropagations,
                             std::vector<Propagation<IT>> &propagations,

                             const IT nVertices) const {
        ttk::Timer timer;

        const IT nPropagations = propagations.size();

        this->printMsg(
          "Finalizing Propagations (" + std::to_string(nPropagations) + ")", 0,
          timer.getElapsedTime(), this->threadNumber_,
          debug::LineMode::REPLACE);

        IT nSegmentVertices = 0;
        IT nParentPropagations = 0;
        parentPropagations.resize(nPropagations);
        for(IT p = 0; p < nPropagations; p++) {
          auto *propagation = &propagations[p];
          if(!propagation->aborted
             && (propagation->parent == propagation
                 || propagation->parent->aborted)) {
            nSegmentVertices = nSegmentVertices + propagation->segmentSize;
            parentPropagations[nParentPropagations++] = propagation;
          }
        }
        parentPropagations.resize(nParentPropagations);

        this->printMsg("Finalizing Propagations ("
                         + std::to_string(nParentPropagations) + "|"
                         + toFixed(nParentPropagations, nPropagations) + "|"
                         + toFixed(nSegmentVertices, nVertices) + ")",
                       1, timer.getElapsedTime(), this->threadNumber_);

        return 1;
      }

      /// This method computes the domain segment of a given propagation. To
      /// this end, it stores a list of all segment vertices on the propagation
      /// data structure, and adds labels to the segmentation array.
      template <typename IT, class TT>
      int computeSegment(IT *segmentation,
                         Propagation<IT> *propagation,

                         const IT *order,
                         const TT *triangulation) const {

        const IT &extremumIndex = propagation->criticalPoints[0];
        const IT &saddleIndex = propagation->criticalPoints.back();
        const IT &saddleOrder = order[saddleIndex];

        // collect segment
        auto &segment = propagation->segment;

        segment.resize(propagation->segmentSize);
        IT segmentIndex = 0;
        if(propagation->segmentSize > 0) {
          std::vector<IT> queue(propagation->segmentSize);
          IT queueIndex = 0;

          // init queue
          {
            queue[queueIndex++] = extremumIndex;
            segmentation[extremumIndex] = -1000; // mark as visited
          }

          // flood fill by starting from extremum
          while(queueIndex > 0) {
            const IT v = queue[--queueIndex];

            segment[segmentIndex++] = v;

            IT nNeighbors = triangulation->getVertexNeighborNumber(v);
            for(IT n = 0; n < nNeighbors; n++) {
              IT u;
              triangulation->getVertexNeighbor(v, n, u);
              auto &s = segmentation[u];
              if(s >= 0 && order[u] > saddleOrder) {
                s = -1000; // mark as visited
                queue[queueIndex++] = u; // add to queue
              }
            }
          }
        }

        if(segmentIndex != propagation->segmentSize) {
          this->printErr("Segment size incorrect: "
                         + std::to_string(segmentIndex) + " "
                         + std::to_string(propagation->segmentSize));
          return 0;
        }

        for(auto idx : propagation->segment)
          segmentation[idx] = extremumIndex;

        return 1;
      }

      /// This method computes the segments of a given list of propagations.
      template <typename IT, class TT>
      int computeSegments(IT *segmentation,
                          std::vector<Propagation<IT> *> &propagations,

                          const IT *order,
                          const TT *triangulation) const {

        const IT nPropagations = propagations.size();
        const IT nVertices = triangulation->getNumberOfVertices();

        ttk::Timer timer;
        const std::string msg
          = "Computing Segments (" + std::to_string(nPropagations) + ")";
        this->printMsg(
          msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

        int status = 1;

// compute segments in parallel
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT p = 0; p < nPropagations; p++) {
          int localStatus
            = this->computeSegment<IT, TT>(segmentation, propagations[p],

                                           order, triangulation);
          if(!localStatus)
            status = 0;
        }
        if(!status)
          return 0;

        // print status
        if(this->debugLevel_ < 4 || nPropagations == 0) {
          this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
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

          this->printMsg("Computing Segments (" + std::to_string(nPropagations)
                           + "|" + toFixed(min, nVertices) + "|"
                           + toFixed(avg, nVertices) + "|"
                           + toFixed(max, nVertices) + ")",
                         1, timer.getElapsedTime(), this->threadNumber_);
        }

        return 1;
      }

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
      }

      template <typename IT, class TT>
      int computeLocalOrderOfSegment(IT *localOrder,

                                     const Propagation<IT> *propagation,
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
        const IT &saddleIndex = propagation->criticalPoints.back();

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
      }

      template <typename IT, class TT>
      int computeLocalOrderOfSegments(
        IT *localOrder,

        const TT *triangulation,
        const IT *segmentation,
        const IT *inputOrder,
        const std::vector<Propagation<IT> *> &propagations) const {
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
          localOrder[propagations[p]->criticalPoints.back()]
            = std::numeric_limits<IT>::max();
        }

        if(this->debugLevel_ < 4 || nPropagations == 0) {
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
                           + std::to_string(min) + "|" + std::to_string(avg)
                           + "|" + std::to_string(max) + ")",
                         1, timer.getElapsedTime(), this->threadNumber_);
        }

        return 1;
      }

      template <typename IT>
      int flattenOrder(
        IT *outputOrder,
        const std::vector<Propagation<IT> *> &parentPropagations) const {

        ttk::Timer timer;
        this->printMsg("Flattening Order Array", 0, 0, this->threadNumber_,
                       debug::LineMode::REPLACE);

        const IT nParentPropagations = parentPropagations.size();

// flatten segemnt to order of last encountered saddles
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT p = 0; p < nParentPropagations; p++) {
          const auto *propagation = parentPropagations[p];
          const auto &saddleOrder
            = outputOrder[propagation->criticalPoints.back()];
          for(const auto &v : propagation->segment)
            outputOrder[v] = saddleOrder;
        }

        this->printMsg("Flattening Order Array", 1, timer.getElapsedTime(),
                       this->threadNumber_);

        return 1;
      }

      template <typename DT, typename IT>
      int flattenScalars(DT *scalars,
                         const std::vector<Propagation<IT>> &propagationsA,
                         const std::vector<Propagation<IT>> &propagationsB
                         = {}) const {

        ttk::Timer timer;
        this->printMsg("Flattening Scalar Array", 0, 0, this->threadNumber_,
                       debug::LineMode::REPLACE);

        std::vector<const std::vector<Propagation<IT>> *> propagationsPair
          = {&propagationsA, &propagationsB};

        for(const auto propagations : propagationsPair) {
          const IT nPropagations = propagations->size();

// flatten scalars to saddle value
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
          for(IT p = 0; p < nPropagations; p++) {
            const auto propagation = &((*propagations)[p]);

            if(!propagation->aborted
               && (propagation->parent == propagation
                   || propagation->parent->aborted)) {
              const auto &saddleScalar
                = scalars[propagation->criticalPoints.back()];
              for(const auto &v : propagation->segment)
                scalars[v] = saddleScalar;
            }
          }
        }

        this->printMsg("Flattening Scalar Array", 1, timer.getElapsedTime(),
                       this->threadNumber_);

        return 1;
      }

      template <typename IT>
      int computeGlobalOrder(
        IT *order,
        const IT *localOrder,
        std::vector<std::tuple<IT, IT, IT>> &sortedIndices) const {
        ttk::Timer timer;

        const IT nVertices = sortedIndices.size();

// init tuples
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT i = 0; i < nVertices; i++) {
          auto &t = sortedIndices[i];
          std::get<0>(t) = order[i];
          std::get<1>(t) = localOrder[i];
          std::get<2>(t) = i;
        }

        this->printMsg("Computing Global Order", 0.2, timer.getElapsedTime(),
                       this->threadNumber_, debug::LineMode::REPLACE);

        TTK_PSORT(
          this->threadNumber_, sortedIndices.begin(), sortedIndices.end());

        this->printMsg("Computing Global Order", 0.8, timer.getElapsedTime(),
                       this->threadNumber_, debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT i = 0; i < nVertices; i++)
          order[std::get<2>(sortedIndices[i])] = i;

        this->printMsg("Computing Global Order", 1, timer.getElapsedTime(),
                       this->threadNumber_);

        return 1;
      }

      template <typename DT, typename IT>
      int computeNumericalPerturbation(
        DT *scalars,
        const std::vector<std::tuple<IT, IT, IT>> &sortedIndices,
        const bool descending = false) const {
        ttk::Timer timer;
        this->printMsg("Applying numerical perturbation", 0, 0,
                       this->threadNumber_, debug::LineMode::REPLACE);

        const IT nVertices = sortedIndices.size();

        if(descending)
          for(IT i = 1; i < nVertices; i++) {
            const IT &v0 = std::get<2>(sortedIndices[i - 1]);
            const IT &v1 = std::get<2>(sortedIndices[i]);
            if(scalars[v0] >= scalars[v1])
              scalars[v1] = boost::math::float_next(scalars[v0]);
          }
        else
          for(IT i = nVertices - 1; i > 0; i--) {
            const IT &v0 = std::get<2>(sortedIndices[i]);
            const IT &v1 = std::get<2>(sortedIndices[i - 1]);
            if(scalars[v0] >= scalars[v1])
              scalars[v1] = boost::math::float_next(scalars[v0]);
          }

        this->printMsg("Applying numerical perturbation", 1,
                       timer.getElapsedTime(), this->threadNumber_);

        return 1;
      }

      template <typename IT, class TT>
      int detectAndRemoveUnauthorizedMaxima(
        IT *order,
        IT *segmentation,
        IT *queueMask,
        IT *localOrder,
        Propagation<IT> **propagationMask,
        std::vector<Propagation<IT>> &propagations,
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
          queueMask, // use as authorization mask (will be overriden by
                     // subsequent procedures)
          localOrder, // use as maxima buffer (will be overriden by subsequent
                      // procedures)

          authorizedExtremaIndices, nAuthorizedExtremaIndices, order,
          triangulation);
        if(!status)
          return 0;

        // compute propagations
        status = this->computeSimplePropagations<IT, TT>(
          propagations, propagationMask, segmentation, queueMask,

          triangulation, order);
        if(!status)
          return 0;

        // finalize master propagations
        std::vector<Propagation<IT> *> parentPropagations;
        status
          = this->finalizePropagations<IT>(parentPropagations, propagations,

                                           nVertices);
        if(!status)
          return 0;

        // compute segments
        status = this->computeSegments<IT, TT>(segmentation, parentPropagations,

                                               order, triangulation);
        if(!status)
          return 0;

        // compute local order of segments
        status = this->computeLocalOrderOfSegments<IT, TT>(
          localOrder,

          triangulation, segmentation, order, parentPropagations);
        if(!status)
          return 0;

        // flatten order
        status = this->flattenOrder<IT>(order, parentPropagations);
        if(!status)
          return 0;

        // compute global offsets
        status = this->computeGlobalOrder<IT>(order, localOrder, sortedIndices);
        if(!status)
          return 0;

        return 1;
      }

      template <typename IT, typename DT, class TT>
      int detectAndRemoveNonPersistentMaxima(
        DT *scalars,
        IT *order,
        IT *segmentation,
        IT *queueMask,
        IT *localOrder,
        Propagation<IT> **propagationMask,
        std::vector<Propagation<IT>> &propagations,
        std::vector<std::tuple<IT, IT, IT>> &sortedIndices,

        const TT *triangulation,
        const DT persistenceThreshold) const {

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
          queueMask, // will be ignored since there are no authorized maxima
          localOrder, // use as maxima buffer (will be overriden by subsequent
                      // procedures)

          nullptr, 0, order, triangulation);
        if(!status)
          return 0;

        // compute propagations
        status = this->computePersistenceSensitivePropagations<IT, DT, TT>(
          propagations, propagationMask, segmentation, queueMask,

          triangulation, order, scalars, persistenceThreshold);
        if(!status)
          return 0;

        // finalize master propagations
        std::vector<Propagation<IT> *> parentPropagations;
        status
          = this->finalizePropagations<IT>(parentPropagations, propagations,

                                           nVertices);
        if(!status)
          return 0;

        // compute segments
        status = this->computeSegments<IT, TT>(segmentation, parentPropagations,

                                               order, triangulation);
        if(!status)
          return 0;

        // compute local order of segments
        status = this->computeLocalOrderOfSegments<IT, TT>(
          localOrder,

          triangulation, segmentation, order, parentPropagations);
        if(!status)
          return 0;

        // flatten order
        status = this->flattenOrder<IT>(order, parentPropagations);
        if(!status)
          return 0;

        // compute global offsets
        status = this->computeGlobalOrder<IT>(order, localOrder, sortedIndices);
        if(!status)
          return 0;

        status = this->flattenScalars<DT, IT>(scalars, propagations);
        if(!status)
          return 0;

        return 1;
      }

      template <typename IT>
      int invertOrder(IT *outputOrder, const IT &nVertices) const {
        ttk::Timer timer;
        this->printMsg("Inverting Order", 0, 0, this->threadNumber_,
                       debug::LineMode::REPLACE);

        const auto nVerticesM1 = nVertices - 1;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT v = 0; v < nVertices; v++)
          outputOrder[v] = nVerticesM1 - outputOrder[v];

        this->printMsg(
          "Inverting Order", 1, timer.getElapsedTime(), this->threadNumber_);

        return 1;
      }

      template <typename DT, typename IT, class TT>
      int removeUnauthorizedExtrema(DT *scalars,
                                    IT *order,

                                    const TT *triangulation,
                                    const IT *authorizedExtremaIndices,
                                    const IT &nAuthorizedExtremaIndices,
                                    const bool &computePerturbation) const {
        ttk::Timer globalTimer;

        IT nVertices = triangulation->getNumberOfVertices();

        // Allocating Memory
        int status = 0;
        std::vector<IT> segmentation;
        std::vector<IT> queueMask;
        std::vector<IT> localOrder;
        std::vector<Propagation<IT> *> propagationMask;
        std::vector<Propagation<IT>> propagationsMax;
        std::vector<Propagation<IT>> propagationsMin;
        std::vector<std::tuple<IT, IT, IT>> sortedIndices;

        this->allocateMemory<IT>(segmentation, queueMask, localOrder,
                                 propagationMask, sortedIndices,

                                 nVertices);

        // Maxima
        {
          this->printMsg("----------- [Removing Unauthorized Maxima]",
                         ttk::debug::Separator::L2);

          status = this->detectAndRemoveUnauthorizedMaxima<IT, TT>(
            order, segmentation.data(), queueMask.data(), localOrder.data(),
            propagationMask.data(), propagationsMax, sortedIndices,

            triangulation, authorizedExtremaIndices, nAuthorizedExtremaIndices);
          if(!status)
            return 0;
        }

        // Minima
        {
          this->printMsg("----------- [Removing Unauthorized Minima]",
                         ttk::debug::Separator::L2);

          // invert order
          if(!this->invertOrder(order, nVertices))
            return 0;

          status = this->detectAndRemoveUnauthorizedMaxima<IT, TT>(
            order, segmentation.data(), queueMask.data(), localOrder.data(),
            propagationMask.data(), propagationsMin, sortedIndices,

            triangulation, authorizedExtremaIndices, nAuthorizedExtremaIndices);
          if(!status)
            return 0;

          // revert order
          if(!this->invertOrder(order, nVertices))
            return 0;
        }

        // flatten scalars
        status = this->flattenScalars<DT, IT>(
          scalars, propagationsMax, propagationsMin);
        if(!status)
          return 0;

        // optionally compute perturbation
        if(computePerturbation) {
          this->printMsg(debug::Separator::L2);
          status = this->computeNumericalPerturbation<DT, IT>(
            scalars, sortedIndices);
          if(!status)
            return 0;
        }

        this->printMsg(debug::Separator::L2);
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime(), this->threadNumber_);

        this->printMsg(debug::Separator::L1);

        return 1;
      }

      template <typename DT, typename IT, class TT>
      int removeNonPersistentExtrema(DT *scalars,
                                     IT *order,

                                     const TT *triangulation,
                                     const DT persistenceThreshold,
                                     const bool &computePerturbation,
                                     const PAIR_TYPE &pairType
                                     = PAIR_TYPE::EXTREMUM_SADDLE) const {
        ttk::Timer globalTimer;

        IT nVertices = triangulation->getNumberOfVertices();

        // Allocating Memory
        int status = 0;
        std::vector<IT> segmentation;
        std::vector<IT> queueMask;
        std::vector<IT> localOrder;
        std::vector<Propagation<IT> *> propagationMask;
        std::vector<Propagation<IT>> propagationsMax;
        std::vector<Propagation<IT>> propagationsMin;
        std::vector<std::tuple<IT, IT, IT>> sortedIndices;

        this->allocateMemory<IT>(segmentation, queueMask, localOrder,
                                 propagationMask, sortedIndices,

                                 nVertices);

        // Maxima
        if(pairType == PAIR_TYPE::EXTREMUM_SADDLE
           || pairType == PAIR_TYPE::MAXIMUM_SADDLE) {
          this->printMsg("----------- [Removing Non-Persistent Maxima]",
                         ttk::debug::Separator::L2);

          status = this->detectAndRemoveNonPersistentMaxima<IT, DT, TT>(
            scalars, order, segmentation.data(), queueMask.data(),
            localOrder.data(), propagationMask.data(), propagationsMax,
            sortedIndices,

            triangulation, persistenceThreshold);
          if(!status)
            return 0;
        }

        // Minima
        if(pairType == PAIR_TYPE::EXTREMUM_SADDLE
           || pairType == PAIR_TYPE::MINIMUM_SADDLE) {
          this->printMsg("----------- [Removing Non-Persistent Minima]",
                         ttk::debug::Separator::L2);

          if(!this->invertOrder(order, nVertices))
            return 0;

          status = this->detectAndRemoveNonPersistentMaxima<IT, DT, TT>(
            scalars, order, segmentation.data(), queueMask.data(),
            localOrder.data(), propagationMask.data(), propagationsMin,
            sortedIndices,

            triangulation, persistenceThreshold);
          if(!status)
            return 0;

          if(!this->invertOrder(order, nVertices))
            return 0;
        }

        // optionally compute perturbation
        if(computePerturbation) {
          this->printMsg(debug::Separator::L2);
          status = this->computeNumericalPerturbation<DT, IT>(
            scalars, sortedIndices, pairType == PAIR_TYPE::MAXIMUM_SADDLE);
          if(!status)
            return 0;
        }

        this->printMsg(debug::Separator::L2);
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime(), this->threadNumber_);

        this->printMsg(debug::Separator::L1);

        return 1;
      }

    }; // class
  } // namespace lts
} // namespace ttk
