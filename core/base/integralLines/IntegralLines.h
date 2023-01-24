/// \ingroup base
/// \class ttk::IntegralLines
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date March 2016
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
#include <Geometry.h>
#include <Triangulation.h>

// std includes
#include <limits>
#include <unordered_set>

namespace ttk {

  struct IntegralLine {
    std::vector<ttk::SimplexId> *trajectory;
    std::vector<double> *distanceFromSeed;
    std::vector<ttk::SimplexId> *localVertexIdentifier;
    ttk::SimplexId seedIdentifier;
    ttk::SimplexId forkIdentifier = -1;
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
    ttk::SimplexId ForkIdentifier;
  };
#endif

  enum Direction { Forward = 0, Backward };

  class IntegralLines : virtual public Debug {

  public:
    IntegralLines();
    ~IntegralLines() override;

    template <class triangulationType>
    inline float getDistance(const triangulationType *triangulation,
                             const SimplexId &a,
                             const SimplexId &b) const {
      float p0[3];
      triangulation->getVertexPoint(a, p0[0], p0[1], p0[2]);
      float p1[3];
      triangulation->getVertexPoint(b, p1[0], p1[1], p1[2]);

      return Geometry::distance(p0, p1, 3);
    }

    template <typename dataType, class triangulationType>
    inline float getGradient(const triangulationType *triangulation,
                             const SimplexId &a,
                             const SimplexId &b,
                             dataType *scalars) const {
      return std::fabs(static_cast<float>(scalars[b] - scalars[a]))
             / getDistance<triangulationType>(triangulation, a, b);
    }

    template <typename dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int execute(const triangulationType *) const;

    template <typename dataType,
              class Compare,
              class triangulationType = ttk::AbstractTriangulation>
    int execute(Compare, const triangulationType *) const;

    inline void setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
    }

    inline void setSeedNumber(const SimplexId &seedNumber) {
      seedNumber_ = seedNumber;
    }

    inline void setDirection(int direction) {
      direction_ = direction;
    }

    void findNextVertex(ttk::SimplexId &vnext,
                        ttk::SimplexId &fnext,
                        std::vector<ttk::SimplexId> &component,
                        const SimplexId *offsets) const;

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
        = {getMPIType(id), getMPIType(id), MPI_DOUBLE,     MPI_DOUBLE,
           getMPIType(id), getMPIType(id), getMPIType(id), getMPIType(id)};
      int lengths[] = {1, 1, 1, 1, 1, 1, 1, 1};
      const long int mpi_offsets[]
        = {offsetof(ElementToBeSent, Id1),
           offsetof(ElementToBeSent, Id2),
           offsetof(ElementToBeSent, DistanceFromSeed1),
           offsetof(ElementToBeSent, DistanceFromSeed2),
           offsetof(ElementToBeSent, LocalVertexIdentifier1),
           offsetof(ElementToBeSent, LocalVertexIdentifier2),
           offsetof(ElementToBeSent, SeedIdentifier),
           offsetof(ElementToBeSent, ForkIdentifier)};
      MPI_Type_create_struct(
        8, lengths, mpi_offsets, types, &(this->MessageType));
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

    inline void setVertexIdentifierScalarField(SimplexId *const data) {
      vertexIdentifierScalarField_ = data;
    }

    inline void
      setOutputTrajectories(std::vector<std::vector<SimplexId>> *trajectories) {
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

    inline void setOutputForkIdentifiers(
      std::vector<ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
        *forkIdentifiers) {
      outputForkIdentifiers_ = forkIdentifiers;
    }

    inline void setChunkSize(int size) {
      chunkSize_ = size;
    }

    inline void buildScalarFieldCriticalPoints() {
      this->scalarFieldCriticalPoints_ = ttk::ScalarFieldCriticalPoints();
    }

  protected:
    SimplexId vertexNumber_;
    SimplexId seedNumber_;
    int direction_;
    void *inputScalarField_;
    const ttk::SimplexId *inputOffsets_;
    std::vector<ttk::SimplexId> *vertexIdentifierScalarField_;
    std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
      *outputTrajectories_;
    std::vector<ttk::ArrayLinkedList<std::vector<double>, TABULAR_SIZE>>
      *outputDistancesFromSeed_;
    std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
      *outputSeedIdentifiers_;
    std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
      *outputForkIdentifiers_;
    std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
      *outputLocalVertexIdentifiers_;
    ttk::ScalarFieldCriticalPoints scalarFieldCriticalPoints_;

#ifdef TTK_ENABLE_MPI
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
  ttk::IntegralLine integralLine = ttk::IntegralLine{
    nullptr, nullptr, nullptr, element.SeedIdentifier, element.ForkIdentifier};
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
  outputForkIdentifiers_->at(threadNum).addArrayElement(element.ForkIdentifier);

  // Add to chunks for task granularity
  chunkIntegralLine[index].trajectory = integralLine.trajectory;
  chunkIntegralLine[index].seedIdentifier = integralLine.seedIdentifier;
  chunkIntegralLine[index].forkIdentifier = integralLine.forkIdentifier;
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
      int rankArray
        = triangulation->getVertexRank(integralLine.trajectory->back());
      if(rankArray != ttk::MPIrank_) {
        ElementToBeSent element = ElementToBeSent{-1,
                                                  -1,
                                                  0,
                                                  0,
                                                  -1,
                                                  -1,
                                                  integralLine.seedIdentifier,
                                                  integralLine.forkIdentifier};
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

void ttk::IntegralLines::findNextVertex(ttk::SimplexId &vnext,
                                        ttk::SimplexId &fnext,
                                        std::vector<ttk::SimplexId> &component,
                                        const SimplexId *offsets) const {
  ttk::SimplexId elementInComponentNumber = component.size();
  for(ttk::SimplexId k = 0; k < elementInComponentNumber; ++k) {
    if(direction_ == static_cast<int>(Direction::Forward)) {
      if(fnext < offsets[component[k]]) {
        vnext = component[k];
        fnext = offsets[component[k]];
      }
    } else {
      if(fnext > offsets[component[k]]) {
        vnext = component[k];
        fnext = offsets[component[k]];
      }
    }
  }
}

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
      if(ttk::isRunningWithMPI()
         && triangulation->getVertexRank(v) == ttk::MPIrank_) {
#if TTK_ENABLE_OPENMP
#pragma omp atomic update seq_cst
#endif
        finishedElement_++;
      } else {
        this->storeToSendIfNecessary<triangulationType>(
          triangulation, integralLine, isMax);
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
        // For each connected components, the max (or the min) is computed
        // and a task is created to further the computation of the integral line
        ttk::SimplexId numberOfComponents = components->size();
#ifdef TTK_ENABLE_MPI
#if TTK_ENABLE_OPENMP
#pragma omp atomic update seq_cst
#endif
        finishedElement_++;
#if TTK_ENABLE_OPENMP
#pragma omp atomic update seq_cst
#endif
        addedElement_ += numberOfComponents;
#endif
        isMax = true;
        for(int i = 0; i < numberOfComponents; i++) {
          SimplexId vnext = -1;
          ttk::SimplexId fnext = offsets[v];
          this->findNextVertex(vnext, fnext, components->at(i), offsets);
          // The task is created
#ifdef TTK_ENABLE_MPI
          ttk::SimplexId forkIdentifier
            = triangulation->getVertexGlobalId(vnext);
#else
          ttk::SimplexId forkIdentifier = vnext;
#endif
          ttk::IntegralLine integralLineFork
            = ttk::IntegralLine{nullptr, nullptr, nullptr,
                                integralLine.seedIdentifier, forkIdentifier};
          triangulation->getVertexPoint(vnext, p1[0], p1[1], p1[2]);
          double distanceFork = Geometry::distance(p0, p1, 3);
          int threadNum{0};
#if TTK_ENABLE_OPENMP
          threadNum = omp_get_thread_num();
#endif
          integralLineFork.trajectory
            = outputTrajectories_->at(threadNum).addArrayElement(
              std::vector<ttk::SimplexId>({v, vnext}));
          integralLineFork.distanceFromSeed
            = outputDistancesFromSeed_->at(threadNum).addArrayElement(
              std::vector<double>({distance, distance + distanceFork}));
          integralLineFork.localVertexIdentifier
            = outputLocalVertexIdentifiers_->at(threadNum).addArrayElement(
              std::vector<ttk::SimplexId>(
                {integralLine.localVertexIdentifier->back(),
                 integralLine.localVertexIdentifier->back() + 1}));
          outputSeedIdentifiers_->at(threadNum).addArrayElement(
            integralLine.seedIdentifier);
          outputForkIdentifiers_->at(threadNum).addArrayElement(forkIdentifier);

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
        SimplexId vnext = -1;
        ttk::SimplexId fnext = offsets[v];
        this->findNextVertex(vnext, fnext, components->at(0), offsets);
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
#ifdef TTK_ENABLE_MPI
        this->storeToSendIfNecessary<triangulationType>(
          triangulation, integralLine, isMax);
#endif
      }
    }
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
    outputForkIdentifiers_->at(threadNum).addArrayElement(
      chunkIntegralLine[j].forkIdentifier);
  }

  {
    std::stringstream msg;
    msg << "Processed " << vertexNumber_ << " points";
    this->printMsg(msg.str(), 1, t.getElapsedTime(), 1);
  }

  return 0;
}

template <typename dataType, class Compare, class triangulationType>
int ttk::IntegralLines::execute(Compare cmp,
                                const triangulationType *triangulation) const {
  const auto offsets = inputOffsets_;
  SimplexId *identifiers
    = static_cast<SimplexId *>(vertexIdentifierScalarField_);
  dataType *scalars = static_cast<dataType *>(inputScalarField_);
  std::vector<std::vector<SimplexId>> *trajectories = outputTrajectories_;

  Timer t;

  // get the seeds
  std::unordered_set<SimplexId> isSeed;
  for(SimplexId k = 0; k < seedNumber_; ++k)
    isSeed.insert(identifiers[k]);
  std::vector<SimplexId> seeds(isSeed.begin(), isSeed.end());
  isSeed.clear();

  trajectories->resize(seeds.size());
  for(SimplexId i = 0; i < (SimplexId)seeds.size(); ++i) {
    SimplexId v{seeds[i]};
    (*trajectories)[i].push_back(v);

    bool isMax{};
    while(!isMax) {
      SimplexId vnext{-1};
      float fnext = std::numeric_limits<float>::min();
      SimplexId neighborNumber = triangulation->getVertexNeighborNumber(v);
      for(SimplexId k = 0; k < neighborNumber; ++k) {
        SimplexId n;
        triangulation->getVertexNeighbor(v, k, n);

        if((direction_ == static_cast<int>(Direction::Forward))
           xor (offsets[n] < offsets[v])) {
          const float f
            = getGradient<dataType, triangulationType>(v, n, scalars);
          if(f > fnext) {
            vnext = n;
            fnext = f;
          }
        }
      }

      if(vnext == -1)
        isMax = true;
      else {
        v = vnext;
        (*trajectories)[i].push_back(v);

        if(cmp(v))
          isMax = true;
      }
    }
  }

  {
    std::stringstream msg;
    msg << "Processed " << vertexNumber_ << " points";
    this->printMsg(msg.str(), 1, t.getElapsedTime(), 1);
  }

  return 0;
}
