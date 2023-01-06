/// \ingroup base
/// \class ttk::IntegralLines
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date March 2016
/// \date MPI implementation: December 2022
///
/// \brief TTK processing package for the computation of edge-based integral
/// lines of the gradient of an input scalar field defined on a PL manifold.
///
/// Given a list of sources, the package produces forward or backward integral
/// lines along the edges of the input triangulation.
///
/// \sa ttkIntegralLines.cpp %for a usage example.

#pragma once

// base code includes
#include <ArrayLinkedList.h>
#include <Geometry.h>
#include <ScalarFieldCriticalPoints.h>
#include <Triangulation.h>
// std includes
#include <algorithm>
#include <iterator>
#include <limits>
#include <unordered_set>
#define TABULAR_SIZE 50
#ifdef TTK_ENABLE_MPI
#define IS_ELEMENT_TO_PROCESS 0
#define IS_MESSAGE_SIZE 1
#endif

namespace ttk {

  struct IntegralLine {
    std::vector<ttk::SimplexId> *trajectory;
    std::vector<double> *distanceFromSeed;
    std::vector<ttk::SimplexId> *localVertexIdentifier;
    ttk::SimplexId seedIdentifier;
  };

#ifdef TTK_ENABLE_MPI
  /* Counters for keeping tracks of number of integral lines
   * left to compute. Each counter is aligned on a cache line
   * to prevent false-sharing
   */
  static int finishedElement_ __attribute__((aligned(64)));
  static int addedElement_ __attribute__((aligned(64)));

  /*
   * For each integral line continuing on another process, we send one layer of
   * ghost cells and we keep one layers for this process, meaning that vertices
   * Id1 belongs to this process and Id2 belongs to a neighboring
   * process.
   */
  struct ElementToBeSent {
    ttk::SimplexId Id1;
    ttk::SimplexId Id2;
    double DistanceFromSeed1;
    double DistanceFromSeed2;
    ttk::SimplexId LocalVertexIdentifier1;
    ttk::SimplexId LocalVertexIdentifier2;
    ttk::SimplexId SeedIdentifier;
  };
#endif

  enum Direction { Forward = 0, Backward };

  class IntegralLines : virtual public Debug {

  public:
    IntegralLines();
    ~IntegralLines() override;

    template <class triangulationType = ttk::AbstractTriangulation>
    int execute(triangulationType *triangulation);

    /**
     * @brief Computes the integral line starting at the vertex of global id
     * seedIdentifier.
     *
     * @tparam triangulationType
     * @param triangulation
     * @param integralLine integral line to compute
     * @param offsets Order array of the scalar array
     */
    template <class triangulationType = ttk::AbstractTriangulation>
    void computeIntegralLine(const triangulationType *triangulation,
                             ttk::IntegralLine integralLine,
                             const ttk::SimplexId *offsets) const;

    /**
     * @brief Create an OpenMP task that contains the computation of nbElement
     * integral lines.
     *
     * @tparam triangulationType
     * @param triangulation
     * @param chunkIntegralLine integral lines to compute within one task
     * @param offsets Order array of the scalar array
     * @param nbElement number of integral lines in chunkIntegralLine
     */
    template <class triangulationType>
    void createTask(const triangulationType *triangulation,
                    std::vector<ttk::IntegralLine> &chunkIntegralLine,
                    const ttk::SimplexId *offsets,
                    int nbElement) const;
    /**
     * @brief Initializes the three attributes of an integral line: the global
     * id of its seed, its trajectory, and the distances of its points with
     * regards to its seed. Then stores the pointers to those objects in
     * chunkIntegralLine to use it for task creation.
     *
     * @tparam triangulationType
     * @param triangulation
     * @param chunkIntegralLine integral lines to compute within one task
     * @param startingIndex index of the first seed of seeds to be used to
     * create an integral line and add it to chunkIntegralLine
     * @param nbElement number of integral lines in chunkIntegralLine
     * @param seeds starting points of the integral lines to be computed
     */
    template <class triangulationType>
    void prepareForTask(const triangulationType *triangulation,
                        std::vector<ttk::IntegralLine> &chunkIntegralLine,
                        int startingIndex,
                        int nbElement,
                        std::vector<SimplexId> *seeds) const;

    inline void setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
    }

    inline void setSeedNumber(const SimplexId &seedNumber) {
      seedNumber_ = seedNumber;
    }

    inline void setDirection(int direction) {
      direction_ = direction;
    }

#ifdef TTK_ENABLE_MPI

    /**
     * @brief Checks if an integral line should be continued on another process
     * or not. If so, encapsulate the necessary data in a struct and stores it
     * in toSend_.
     *
     * @tparam triangulationType
     * @param triangulation
     * @param integralLine integral line to inspect
     * @param isMax boolean that stores whether the computation should be
     * furthered locally
     */
    template <class triangulationType>
    void storeToSendIfNecessary(const triangulationType *triangulation,
                                ttk::IntegralLine integralLine,
                                bool &isMax) const;

    /**
     * @brief Extract the data of element to initialize the three attributes of
     * an integral line and stores their pointers the in chunk vectors at index.
     * When chunk vectors are full, the task is created and index is
     * reinitialized to 0.
     *
     * @tparam triangulationType
     * @param triangulation
     * @param element received integral line to unpack
     * @param chunkIntegralLine integral lines to compute within one task
     * @param index current number of integral lines in chunks
     * @param taskSize number of integral lines that should be computed in one
     * task
     * @param offsets Order array of the scalar array
     */
    template <class triangulationType>
    void receiveElement(const triangulationType *triangulation,
                        ElementToBeSent &element,
                        std::vector<ttk::IntegralLine> &chunkIntegralLine,
                        int &index,
                        int taskSize,
                        const ttk::SimplexId *offsets);

    inline void
      setToSend(std::vector<std::vector<std::vector<ElementToBeSent>>> *send) {
      toSend_ = send;
    }

    inline void setGlobalElementCounter(int counter) {
      globalElementCounter_ = counter;
    }

    inline void setNeighbors(const std::vector<int> &neighbors) {
      neighbors_ = neighbors;
      neighborNumber_ = neighbors_.size();
      int idx = 0;
      for(int neighbor : neighbors) {
        neighborsToId_[neighbor] = idx;
        idx++;
      }
    }
    /**
     * @brief Create a MPI type for sending integral lines from
     * one process to another
     *
     */
    void createMessageType() {
      ttk::SimplexId id = 0;
      MPI_Datatype types[]
        = {getMPIType(id), getMPIType(id), MPI_DOUBLE,    MPI_DOUBLE,
           getMPIType(id), getMPIType(id), getMPIType(id)};
      int lengths[] = {1, 1, 1, 1, 1, 1, 1};
      const long int mpi_offsets[]
        = {offsetof(ElementToBeSent, Id1),
           offsetof(ElementToBeSent, Id2),
           offsetof(ElementToBeSent, DistanceFromSeed1),
           offsetof(ElementToBeSent, DistanceFromSeed2),
           offsetof(ElementToBeSent, LocalVertexIdentifier1),
           offsetof(ElementToBeSent, LocalVertexIdentifier2),
           offsetof(ElementToBeSent, SeedIdentifier)};
      MPI_Type_create_struct(
        7, lengths, mpi_offsets, types, &(this->MessageType));
      MPI_Type_commit(&(this->MessageType));
    }
#endif

    inline void setOutputLocalVertexIdentifiers(
      std::vector<ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
        *localVertexIdentifiers) {
      outputLocalVertexIdentifiers_ = localVertexIdentifiers;
    }

    int preconditionTriangulation(ttk::AbstractTriangulation *triangulation) {
      int status = triangulation->preconditionVertexNeighbors();
      status += triangulation->preconditionEdges();
      status += triangulation->preconditionEdgeStars();
      status += triangulation->preconditionVertexEdges();
      // For critical points computation
      status += triangulation->preconditionVertexStars();
      status += triangulation->preconditionBoundaryVertices();
      this->scalarFieldCriticalPoints_.setDomainDimension(2);
      return status;
    }

    inline void setInputScalarField(void *data) {
      inputScalarField_ = data;
    }

    /**
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p data buffer prior to any
     * computation (the VTK wrapper already includes a mechanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    inline void setInputOffsets(const SimplexId *const data) {
      inputOffsets_ = data;
    }

    inline void
      setVertexIdentifierScalarField(std::vector<SimplexId> *const data) {
      vertexIdentifierScalarField_ = data;
    }

    inline void setOutputTrajectories(
      std::vector<ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
        *trajectories) {
      outputTrajectories_ = trajectories;
    }

    inline void setOutputDistancesFromSeed(
      std::vector<ArrayLinkedList<std::vector<double>, TABULAR_SIZE>>
        *distancesFromSeed) {
      outputDistancesFromSeed_ = distancesFromSeed;
    }

    inline void setOutputSeedIdentifiers(
      std::vector<ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
        *seedIdentifiers) {
      outputSeedIdentifiers_ = seedIdentifiers;
    }

    inline void setChunkSize(int size) {
      chunkSize_ = size;
    }

    inline void buildScalarFieldCriticalPoints() {
      this->scalarFieldCriticalPoints_ = ttk::ScalarFieldCriticalPoints();
    }

  protected:
    ttk::SimplexId vertexNumber_;
    ttk::SimplexId seedNumber_;
    ttk::SimplexId chunkSize_;
    ttk::SimplexId direction_;
    void *inputScalarField_;
    const ttk::SimplexId *inputOffsets_;
    std::vector<ttk::SimplexId> *vertexIdentifierScalarField_;
    std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
      *outputTrajectories_;
    std::vector<ttk::ArrayLinkedList<std::vector<double>, TABULAR_SIZE>>
      *outputDistancesFromSeed_;
    std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
      *outputSeedIdentifiers_;
    std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
      *outputLocalVertexIdentifiers_;
    ttk::ScalarFieldCriticalPoints scalarFieldCriticalPoints_;

#ifdef TTK_ENABLE_MPI
    const int *vertexRankArray_{nullptr};
    std::vector<std::vector<std::vector<ElementToBeSent>>> *toSend_{nullptr};
    int neighborNumber_;
    std::unordered_map<int, int> neighborsToId_;
    std::vector<int> neighbors_;
    SimplexId keepWorking_;
    SimplexId globalElementCounter_;
    MPI_Datatype MessageType;
#endif
  };
} // namespace ttk

#ifdef TTK_ENABLE_MPI

template <class triangulationType>
void ttk::IntegralLines::receiveElement(
  const triangulationType *triangulation,
  ElementToBeSent &element,
  std::vector<ttk::IntegralLine> &chunkIntegralLine,
  int &index,
  int taskSize,
  const ttk::SimplexId *offsets) {
  ttk::IntegralLine integralLine
    = ttk::IntegralLine{nullptr, nullptr, nullptr, element.SeedIdentifier};
  // Create integral line object on this process
  int threadNum{0};
#if TTK_ENABLE_OPENMP
  threadNum = omp_get_thread_num();
#endif
  integralLine.trajectory = outputTrajectories_->at(threadNum).addArrayElement(
    std::vector<ttk::SimplexId>(
      {triangulation->getVertexLocalId(element.Id1),
       triangulation->getVertexLocalId(element.Id2)}));
  integralLine.distanceFromSeed
    = outputDistancesFromSeed_->at(threadNum).addArrayElement(
      std::vector<double>(
        {element.DistanceFromSeed1, element.DistanceFromSeed2}));
  integralLine.localVertexIdentifier
    = outputLocalVertexIdentifiers_->at(threadNum).addArrayElement(
      std::vector<ttk::SimplexId>(
        {element.LocalVertexIdentifier1, element.LocalVertexIdentifier2}));
  outputSeedIdentifiers_->at(threadNum).addArrayElement(element.SeedIdentifier);

  // Add to chunks for task granularity
  chunkIntegralLine[index].trajectory = integralLine.trajectory;
  chunkIntegralLine[index].seedIdentifier = integralLine.seedIdentifier;
  chunkIntegralLine[index].distanceFromSeed = integralLine.distanceFromSeed;
  chunkIntegralLine[index].localVertexIdentifier
    = integralLine.localVertexIdentifier;

  // Start task if necessary
  if(index == taskSize - 1) {
    this->createTask<triangulationType>(
      triangulation, chunkIntegralLine, offsets, taskSize);
    index = 0;
  } else {
    index++;
  }
}

template <class triangulationType>
void ttk::IntegralLines::storeToSendIfNecessary(
  const triangulationType *triangulation,
  ttk::IntegralLine integralLine,
  bool &isMax) const {
#ifdef TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    int size = integralLine.trajectory->size();
    if(size > 1) {
      int rankArray = vertexRankArray_[integralLine.trajectory->back()];
      if(rankArray != ttk::MPIrank_) {
        ElementToBeSent element
          = ElementToBeSent{-1, -1, 0, 0, -1, -1, integralLine.seedIdentifier};
        element.Id2
          = triangulation->getVertexGlobalId(integralLine.trajectory->back());
        element.Id1 = triangulation->getVertexGlobalId(
          integralLine.trajectory->at(size - 2));
        element.DistanceFromSeed2 = integralLine.distanceFromSeed->back();
        element.DistanceFromSeed1 = integralLine.distanceFromSeed->at(size - 2);
        element.LocalVertexIdentifier2
          = integralLine.localVertexIdentifier->back();
        element.LocalVertexIdentifier1
          = integralLine.localVertexIdentifier->at(size - 2);
#if TTK_ENABLE_OPENMP
        toSend_
          ->at(neighborsToId_.find(rankArray)->second)[omp_get_thread_num()]
          .push_back(element);
#else
        toSend_->at(neighborsToId_.find(rankArray)->second)[0].push_back(
          element);
#endif
        isMax = true;
      }
    }
  }
#endif
}
#endif

template <class triangulationType>
void ttk::IntegralLines::computeIntegralLine(
  const triangulationType *triangulation,
  ttk::IntegralLine integralLine,
  const SimplexId *offsets) const {
  double distance = integralLine.distanceFromSeed->back();
  ttk::SimplexId v = integralLine.trajectory->back();
  float p0[3];
  float p1[3];
  triangulation->getVertexPoint(v, p0[0], p0[1], p0[2]);
  bool isMax{};
  std::vector<std::vector<ttk::SimplexId>> *components;
  while(!isMax) {
    std::vector<std::vector<ttk::SimplexId>> upperComponents;
    std::vector<std::vector<ttk::SimplexId>> lowerComponents;
    // Computation of the critical type
    char criticalType
      = this->scalarFieldCriticalPoints_.getCriticalType<triangulationType>(
        v, offsets, triangulation, &upperComponents, &lowerComponents);
    // End of integral line if an appropriate maxima is reached
    if((criticalType == (char)CriticalType::Local_maximum
        && direction_ == static_cast<int>(Direction::Forward))
       || (criticalType == (char)CriticalType::Local_minimum
           && direction_ != static_cast<int>(Direction::Forward))) {
      isMax = true;
#ifdef TTK_ENABLE_MPI
      if(ttk::isRunningWithMPI() && vertexRankArray_[v] == ttk::MPIrank_) {
#if TTK_ENABLE_OPENMP
#pragma omp atomic update seq_cst
#endif
        finishedElement_++;
      }
#endif
    } else {
      // If the integral line continues, the computation continues using lower
      // and upper components
      components = &(lowerComponents);
      if(criticalType == (char)CriticalType::Saddle1
         || criticalType == (char)CriticalType::Saddle2
         || criticalType == (char)CriticalType::Degenerate) {
        if(direction_ == static_cast<int>(Direction::Forward)) {
          components = &(upperComponents);
        }
      } else {
        // In case the vertex is not a saddle point, all neighbor vertices
        // are used for the computation
        components->clear();
        ttk::SimplexId neighborNumber
          = triangulation->getVertexNeighborNumber(v);
        components->push_back(std::vector<ttk::SimplexId>());
        ttk::SimplexId id;
        for(ttk::SimplexId i = 0; i < neighborNumber; i++) {
          triangulation->getVertexNeighbor(v, i, id);
          components->at(0).push_back(id);
        }
      }
      // For each connected components, the max (or the min) is computed
      // and a task is created to further the computation of the integral line
      // except for one component that keeps computing within the current task
      ttk::SimplexId numberOfComponents = components->size();
      for(int i = 0; i < numberOfComponents; i++) {
        SimplexId vnext = -1;
        ttk::SimplexId fnext = offsets[v];
        ttk::SimplexId elementInComponentNumber = components->at(i).size();
        for(SimplexId k = 0; k < elementInComponentNumber; ++k) {
          if(direction_ == static_cast<int>(Direction::Forward)) {
            if(fnext < offsets[components->at(i)[k]]) {
              vnext = components->at(i)[k];
              fnext = offsets[components->at(i)[k]];
            }
          } else {
            if(fnext > offsets[components->at(i)[k]]) {
              vnext = components->at(i)[k];
              fnext = offsets[components->at(i)[k]];
            }
          }
        }
        if(i == numberOfComponents - 1) {
          triangulation->getVertexPoint(vnext, p1[0], p1[1], p1[2]);
          distance += Geometry::distance(p0, p1, 3);
          integralLine.trajectory->push_back(vnext);
          p0[0] = p1[0];
          p0[1] = p1[1];
          p0[2] = p1[2];
          integralLine.distanceFromSeed->push_back(distance);
          integralLine.localVertexIdentifier->push_back(
            integralLine.localVertexIdentifier->back() + 1);
          v = vnext;
        } else {
          // The task is created
#ifdef TTK_ENABLE_MPI
          ttk::SimplexId seedIdentifier = triangulation->getVertexGlobalId(v);
#else
          ttk::SimplexId seedIdentifier = v;
#endif
          ttk::IntegralLine integralLineFork
            = ttk::IntegralLine{nullptr, nullptr, nullptr, seedIdentifier};
          triangulation->getVertexPoint(vnext, p1[0], p1[1], p1[2]);
          double distanceFork = Geometry::distance(p0, p1, 3);
#ifdef TTK_ENABLE_MPI
#if TTK_ENABLE_OPENMP
#pragma omp atomic update seq_cst
#endif
          addedElement_++;
#endif
          int threadNum{0};
#if TTK_ENABLE_OPENMP
          threadNum = omp_get_thread_num();
#endif
          integralLineFork.trajectory
            = outputTrajectories_->at(threadNum).addArrayElement(
              std::vector<ttk::SimplexId>({v, vnext}));
          integralLineFork.distanceFromSeed
            = outputDistancesFromSeed_->at(threadNum).addArrayElement(
              std::vector<double>({0, distanceFork}));
          integralLineFork.localVertexIdentifier
            = outputLocalVertexIdentifiers_->at(threadNum).addArrayElement(
              std::vector<ttk::SimplexId>(
                {integralLine.localVertexIdentifier->back() + 1,
                 integralLine.localVertexIdentifier->back() + 2}));
          outputSeedIdentifiers_->at(threadNum).addArrayElement(seedIdentifier);
#if TTK_ENABLE_OPENMP
#pragma omp task firstprivate(integralLineFork)
          {
#endif
#ifdef TTK_ENABLE_MPI
            bool hasBeenSent = false;
            this->storeToSendIfNecessary<triangulationType>(
              triangulation, integralLineFork, hasBeenSent);
            if(!hasBeenSent) {
#endif
              this->computeIntegralLine<triangulationType>(
                triangulation, integralLineFork, offsets);
#ifdef TTK_ENABLE_MPI
            }
#endif
#if TTK_ENABLE_OPENMP
          }
#endif
        }
      }
    }
#ifdef TTK_ENABLE_MPI
    this->storeToSendIfNecessary<triangulationType>(
      triangulation, integralLine, isMax);
#endif
  }
}

template <class triangulationType>
void ttk::IntegralLines::prepareForTask(
#ifdef TTK_ENABLE_MPI
  const triangulationType *triangulation,
#else
  const triangulationType *ttkNotUsed(triangulation),
#endif
  std::vector<ttk::IntegralLine> &chunkIntegralLine,
  int startingIndex,
  int nbElement,
  std::vector<SimplexId> *seeds) const {

  for(SimplexId j = 0; j < nbElement; j++) {
    SimplexId v{seeds->at(j + startingIndex)};
#ifdef TTK_ENABLE_MPI
    chunkIntegralLine[j].seedIdentifier = triangulation->getVertexGlobalId(v);
#else
    chunkIntegralLine[j].seedIdentifier = v;
#endif
    int threadNum{0};
#if TTK_ENABLE_OPENMP
    threadNum = omp_get_thread_num();
#endif
    chunkIntegralLine[j].trajectory
      = outputTrajectories_->at(threadNum).addArrayElement(
        std::vector<ttk::SimplexId>(1, v));
    chunkIntegralLine[j].distanceFromSeed
      = outputDistancesFromSeed_->at(threadNum).addArrayElement(
        std::vector<double>(1, 0));
    chunkIntegralLine[j].localVertexIdentifier
      = outputLocalVertexIdentifiers_->at(threadNum).addArrayElement(
        std::vector<ttk::SimplexId>(1, 0));
    outputSeedIdentifiers_->at(threadNum).addArrayElement(
      chunkIntegralLine[j].seedIdentifier);
  }
}

template <class triangulationType>
void ttk::IntegralLines::createTask(
  const triangulationType *triangulation,
  std::vector<ttk::IntegralLine> &chunkIntegralLine,
  const ttk::SimplexId *offsets,
  int nbElement) const {
#if TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkIntegralLine)
  {
#endif
    for(int j = 0; j < nbElement; j++) {
      this->computeIntegralLine<triangulationType>(
        triangulation, chunkIntegralLine[j], offsets);
    }
#if TTK_ENABLE_OPENMP
  }
#endif
}

template <class triangulationType>
int ttk::IntegralLines::execute(triangulationType *triangulation) {

#ifdef TTK_ENABLE_MPI
  keepWorking_ = 1;
  finishedElement_ = 0;
  addedElement_ = 0;
#endif
  const SimplexId *offsets = inputOffsets_;
  std::vector<SimplexId> *seeds = vertexIdentifierScalarField_;
  Timer t;

  std::vector<ttk::IntegralLine> chunkIntegralLine(chunkSize_);
  int taskNumber = (int)seedNumber_ / chunkSize_;
#ifdef TTK_ENABLE_OPENMP
#ifdef TTK_ENABLE_MPI
#pragma omp parallel shared(finishedElement_, toSend_, addedElement_) \
  num_threads(threadNumber_)
  {
#else
#pragma omp parallel num_threads(threadNumber_)
  {
#endif // TTK_ENABLE_MPI
#pragma omp master
    {
#endif // TTK_ENABLE_OPENMP
      for(SimplexId i = 0; i < taskNumber; ++i) {
        this->prepareForTask<triangulationType>(
          triangulation, chunkIntegralLine, i * chunkSize_, chunkSize_, seeds);
        this->createTask<triangulationType>(
          triangulation, chunkIntegralLine, offsets, chunkSize_);
      }
      int rest = seedNumber_ % chunkSize_;
      if(rest > 0) {
        this->prepareForTask<triangulationType>(
          triangulation, chunkIntegralLine, taskNumber * chunkSize_, rest,
          seeds);
        this->createTask<triangulationType>(
          triangulation, chunkIntegralLine, offsets, rest);
      }
#ifdef TTK_ENABLE_OPENMP
    }
  }
#endif
#ifdef TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    int i;
    int finishedElementReceived = 0;
    std::vector<int> sendMessageSize(neighborNumber_);
    std::vector<int> recvMessageSize(neighborNumber_);
    std::vector<std::vector<ElementToBeSent>> send_buf(neighborNumber_);
    std::vector<std::vector<ElementToBeSent>> recv_buf(neighborNumber_);
    for(i = 0; i < neighborNumber_; i++) {
      send_buf.reserve((int)seedNumber_ * 0.005);
      recv_buf.reserve((int)seedNumber_ * 0.005);
    }
    std::vector<MPI_Request> requests(2 * neighborNumber_, MPI_REQUEST_NULL);
    std::vector<MPI_Status> statuses(4 * neighborNumber_);
    int taskSize;
    int index;
    int totalMessageSize;
    while(keepWorking_) {
      finishedElement_ -= addedElement_;
      addedElement_ = 0;
      // Exchange of the number of integral lines finished on all processes
      MPI_Allreduce(&finishedElement_, &finishedElementReceived, 1, MPI_INTEGER,
                    MPI_SUM, ttk::MPIcomm_);
      finishedElement_ = 0;
      // Update the number of integral lines left to compute
      globalElementCounter_ -= finishedElementReceived;
      // Stop working in case there are no more computation to be done
      if(globalElementCounter_ == 0) {
        keepWorking_ = 0;
      }
      if(keepWorking_) {
        totalMessageSize = 0;
        // Preparation of the send buffers and exchange of the size of messages
        // to be sent
        for(i = 0; i < neighborNumber_; i++) {
          for(int j = 0; j < threadNumber_; j++) {
            send_buf[i].insert(send_buf[i].end(), toSend_->at(i)[j].begin(),
                               toSend_->at(i)[j].end());
            toSend_->at(i)[j].clear();
          }
          sendMessageSize[i] = (int)send_buf[i].size();
          MPI_Isend(&sendMessageSize[i], 1, MPI_INTEGER, neighbors_.at(i),
                    IS_MESSAGE_SIZE, ttk::MPIcomm_, &requests[2 * i]);
          MPI_Irecv(&recvMessageSize[i], 1, MPI_INTEGER, neighbors_.at(i),
                    IS_MESSAGE_SIZE, ttk::MPIcomm_, &requests[2 * i + 1]);
        }
        MPI_Waitall(2 * neighborNumber_, requests.data(), MPI_STATUSES_IGNORE);
        // Exchange of the data
        for(i = 0; i < neighborNumber_; i++) {
          if(recv_buf[i].size() < (size_t)recvMessageSize[i]) {
            recv_buf[i].resize(recvMessageSize[i]);
          }
          if(recvMessageSize[i] > 0) {
            MPI_Irecv(recv_buf[i].data(), recvMessageSize[i], this->MessageType,
                      neighbors_.at(i), IS_ELEMENT_TO_PROCESS, ttk::MPIcomm_,
                      &requests[2 * i]);
            totalMessageSize += recvMessageSize[i];
          }

          if(sendMessageSize[i] > 0) {
            MPI_Isend(send_buf[i].data(), sendMessageSize[i], this->MessageType,
                      neighbors_.at(i), IS_ELEMENT_TO_PROCESS, ttk::MPIcomm_,
                      &requests[2 * i + 1]);
          }
        }
        MPI_Waitall(2 * neighborNumber_, requests.data(), MPI_STATUSES_IGNORE);
        for(i = 0; i < neighborNumber_; i++) {
          send_buf[i].clear();
        }
        // Extraction of the received data and creation of the tasks
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel shared(finishedElement_, toSend_) \
  num_threads(threadNumber_)
        {
#pragma omp master
          {
#endif // TTK_ENABLE_OPENMP
            index = 0;
            taskSize = std::min(
              (ttk::SimplexId)std::max(totalMessageSize / (threadNumber_ * 100),
                                       std::min(totalMessageSize, 50)),
              chunkSize_);
            chunkIntegralLine.resize(taskSize);
            for(i = 0; i < neighborNumber_; i++) {
              for(int j = 0; j < recvMessageSize[i]; j++) {
                this->receiveElement<triangulationType>(
                  triangulation, recv_buf[i][j], chunkIntegralLine, index,
                  taskSize, offsets);
              }
            }
            if(index > 0) {
              this->createTask<triangulationType>(
                triangulation, chunkIntegralLine, offsets, index);
            }
#ifdef TTK_ENABLE_OPENMP
          }
        }
#endif // TTK_ENABLE_OPENMP
      }
    }
  }
#endif // TTK_ENABLE_MPI
  {
    std::stringstream msg;
    msg << "Processed " << vertexNumber_ << " points";
    this->printMsg(msg.str(), 1, t.getElapsedTime(), threadNumber_);
  }
  return 0;
}
