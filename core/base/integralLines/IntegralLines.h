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
#include <numeric>
#include <string>
#include <unordered_set>

#define INTEGRAL_LINE_TABULAR_SIZE 50
#ifdef TTK_ENABLE_MPI
#define INTEGRAL_LINE_IS_ELEMENT_TO_PROCESS 0
#define INTEGRAL_LINE_IS_MESSAGE_SIZE 1
#endif

namespace ttk {
  namespace intgl {

    /**
     * @brief Struct containing the data of an integral line.
     * trajectories: vector of identifiers of each vertex the integral line
     * passes on distanceFromSeed: distance of each vertex for the integral line
     * from the seed localVertexIdentifiers: vector of local vertex ids. It
     * represents the number of the vertex in the integral line. seedIdentifier:
     * identifier of the seed of the integral line forkIdentifier: identifier of
     * the last fork the integral line encountered
     */
    struct IntegralLine {
      std::vector<ttk::SimplexId> trajectory;
      std::vector<double> distanceFromSeed;
      std::vector<ttk::SimplexId> localVertexIdentifier;
      ttk::SimplexId seedIdentifier;
      ttk::SimplexId forkIdentifier = -1;
    };

#ifdef TTK_ENABLE_MPI
    /**
     * @brief Struct used for sorting ghost data during the generation of global
     * ids. A GhostElementsToSort object is created for each segment of integral
     * line that ends or begins with a ghost vertex. minLocalVertexId,
     * seedIdentifier and forkIdentifier are used to identify this segment of
     * integral line on different processes.
     */
    struct GhostElementsToSort {
      ttk::SimplexId ownedGlobalId;
      ttk::SimplexId minLocalVertexId;
      ttk::SimplexId seedIdentifier;
      ttk::SimplexId forkIdentifier{-1};
      ttk::SimplexId globalEdgeId{-1};
      ttk::SimplexId ghostVertexLocalId;
      ttk::SimplexId ghostEdgeLocalId;
    };

    inline bool operator<(const GhostElementsToSort &left,
                          const GhostElementsToSort &right) {
      if(left.seedIdentifier != right.seedIdentifier) {
        return left.seedIdentifier < right.seedIdentifier;
      }
      if(left.forkIdentifier != right.forkIdentifier) {
        return left.forkIdentifier < right.forkIdentifier;
      }
      return left.minLocalVertexId < right.minLocalVertexId;
    };

    /* Counters for keeping tracks of number of integral lines
     * left to compute. Each counter is aligned on a cache line
     * to prevent false-sharing
     */
    static int finishedElement_ __attribute__((aligned(64)));
    static int addedElement_ __attribute__((aligned(64)));

    /*
     * For each integral line continuing on another process, we send one layer
     * of ghost cells and we keep one layers for this process, meaning that
     * vertices Id1 belongs to this process and Id2 belongs to a neighboring
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
  } // namespace intgl

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
                             ttk::intgl::IntegralLine *integralLine,
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
                    std::vector<ttk::intgl::IntegralLine *> &chunkIntegralLine,
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
    void
      prepareForTask(const triangulationType *triangulation,
                     std::vector<ttk::intgl::IntegralLine *> &chunkIntegralLine,
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

    /**
     * @brief Finds the vertex of highest or lowest offsets (depending on the
     * direction of the integral line) in the component vector.
     *
     * @param vnext identifier of the vertex
     * @param fnext offset of the vertex
     * @param component list of identifiers to do the search in
     * @param offsets Order array of the scalar array
     */
    void findNextVertex(ttk::SimplexId &vnext,
                        ttk::SimplexId &fnext,
                        std::vector<ttk::SimplexId> &component,
                        const SimplexId *offsets) const;

#ifdef TTK_ENABLE_MPI

    /**
     * @brief Constructs the global identifiers for vertices and edges.
     *
     * @tparam triangulationType
     * @param globalVertexId vector of global identifiers for vertices
     * @param globalCellId vector of global identifiers for cells
     * @param integralLines linked list of integralLines
     * @param triangulation Triangulation
     * @return int
     */
    template <typename triangulationType>
    int getGlobalIdentifiers(
      std::vector<ttk::SimplexId> &globalVertexId,
      std::vector<ttk::SimplexId> &globalCellId,
      const std::vector<ttk::ArrayLinkedList<ttk::intgl::IntegralLine,
                                             INTEGRAL_LINE_TABULAR_SIZE>>
        &integralLines,
      const triangulationType *triangulation);

    /**
     * @brief Sorts ghosts and exchanges their global identifiers between
     * processes.
     *
     * @param unmatchedGhosts Ghosts vertices and edges for which the global
     * identifier is to be determined
     * @param globalVertexId vector of global identifiers for vertices
     * @param globalCellId vector of global identifiers for cells
     * @return int 0 for success
     */
    int exchangeGhosts(
      const std::vector<std::vector<ttk::intgl::GhostElementsToSort>>
        &unmatchedGhosts,
      std::vector<ttk::SimplexId> &globalVertexId,
      std::vector<ttk::SimplexId> &globalCellId);

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
                                ttk::intgl::IntegralLine *integralLine,
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
    void
      receiveElement(const triangulationType *triangulation,
                     intgl::ElementToBeSent &element,
                     std::vector<ttk::intgl::IntegralLine *> &chunkIntegralLine,
                     int &index,
                     int taskSize,
                     const ttk::SimplexId *offsets);

    inline void setToSend(
      std::vector<std::vector<std::vector<intgl::ElementToBeSent>>> *send) {
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
        = {offsetof(intgl::ElementToBeSent, Id1),
           offsetof(intgl::ElementToBeSent, Id2),
           offsetof(intgl::ElementToBeSent, DistanceFromSeed1),
           offsetof(intgl::ElementToBeSent, DistanceFromSeed2),
           offsetof(intgl::ElementToBeSent, LocalVertexIdentifier1),
           offsetof(intgl::ElementToBeSent, LocalVertexIdentifier2),
           offsetof(intgl::ElementToBeSent, SeedIdentifier),
           offsetof(intgl::ElementToBeSent, ForkIdentifier)};
      MPI_Type_create_struct(
        8, lengths, mpi_offsets, types, &(this->MessageType));
      MPI_Type_commit(&(this->MessageType));
    }
#endif

    int preconditionTriangulation(ttk::AbstractTriangulation *triangulation) {
      int status = triangulation->preconditionVertexNeighbors();
      // For critical points computation
      status += triangulation->preconditionVertexStars();
      status += triangulation->preconditionBoundaryVertices();
      this->scalarFieldCriticalPoints_.preconditionTriangulation(triangulation);
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

    inline void setOutputIntegralLines(
      std::vector<ttk::ArrayLinkedList<ttk::intgl::IntegralLine,
                                       INTEGRAL_LINE_TABULAR_SIZE>>
        *integralLines) {
      this->outputIntegralLines_ = integralLines;
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
    std::vector<ttk::ArrayLinkedList<ttk::intgl::IntegralLine,
                                     INTEGRAL_LINE_TABULAR_SIZE>>
      *outputIntegralLines_;
    ttk::ScalarFieldCriticalPoints scalarFieldCriticalPoints_;
    bool EnableForking{false};

#ifdef TTK_ENABLE_MPI
    std::vector<std::vector<std::vector<intgl::ElementToBeSent>>> *toSend_{
      nullptr};
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
  intgl::ElementToBeSent &element,
  std::vector<ttk::intgl::IntegralLine *> &chunkIntegralLine,
  int &index,
  int taskSize,
  const ttk::SimplexId *offsets) {

  // Create integral line object on this process
  int threadNum{0};
#ifdef TTK_ENABLE_OPENMP4
  threadNum = omp_get_thread_num();
#endif // TTK_ENABLE_OPENMP4

  ttk::intgl::IntegralLine *integralLine
    = outputIntegralLines_->at(threadNum).addArrayElement(
      ttk::intgl::IntegralLine{
        std::vector<ttk::SimplexId>(
          {triangulation->getVertexLocalId(element.Id1),
           triangulation->getVertexLocalId(element.Id2)}),
        std::vector<double>(
          {element.DistanceFromSeed1, element.DistanceFromSeed2}),
        std::vector<ttk::SimplexId>(
          {element.LocalVertexIdentifier1, element.LocalVertexIdentifier2}),
        element.SeedIdentifier, element.ForkIdentifier});

  // Add to chunks for task granularity
  chunkIntegralLine[index] = integralLine;

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
  ttk::intgl::IntegralLine *integralLine,
  bool &isMax) const {
#ifdef TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    int size = integralLine->trajectory.size();
    if(size > 1) {
      int rankArray
        = triangulation->getVertexRank(integralLine->trajectory.back());
      if(rankArray != ttk::MPIrank_) {
        intgl::ElementToBeSent element
          = intgl::ElementToBeSent{-1,
                                   -1,
                                   0,
                                   0,
                                   -1,
                                   -1,
                                   integralLine->seedIdentifier,
                                   integralLine->forkIdentifier};
        element.Id2
          = triangulation->getVertexGlobalId(integralLine->trajectory.back());
        element.Id1 = triangulation->getVertexGlobalId(
          integralLine->trajectory.at(size - 2));
        element.DistanceFromSeed2 = integralLine->distanceFromSeed.back();
        element.DistanceFromSeed1 = integralLine->distanceFromSeed.at(size - 2);
        element.LocalVertexIdentifier2
          = integralLine->localVertexIdentifier.back();
        element.LocalVertexIdentifier1
          = integralLine->localVertexIdentifier.at(size - 2);
#ifdef TTK_ENABLE_OPENMP4
        toSend_
          ->at(neighborsToId_.find(rankArray)->second)[omp_get_thread_num()]
          .push_back(element);
#else
        toSend_->at(neighborsToId_.find(rankArray)->second)[0].push_back(
          element);
#endif // TTK_ENABLE_OPENMP4
        isMax = true;
      }
    }
  }
#endif
}
#endif

inline void
  ttk::IntegralLines::findNextVertex(ttk::SimplexId &vnext,
                                     ttk::SimplexId &fnext,
                                     std::vector<ttk::SimplexId> &component,
                                     const SimplexId *offsets) const {
  ttk::SimplexId const elementInComponentNumber = component.size();
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
  ttk::intgl::IntegralLine *integralLine,
  const SimplexId *offsets) const {
  double distance = integralLine->distanceFromSeed.back();
  ttk::SimplexId v = integralLine->trajectory.back();
  float p0[3];
  float p1[3];
  triangulation->getVertexPoint(v, p0[0], p0[1], p0[2]);
  bool isMax{};
  std::vector<std::vector<ttk::SimplexId>> *components;
  while(!isMax) {
    std::vector<std::vector<ttk::SimplexId>> upperComponents;
    std::vector<std::vector<ttk::SimplexId>> lowerComponents;
    // Computation of the critical type
    char const criticalType
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
#ifdef TTK_ENABLE_OPENMP4
#pragma omp atomic update seq_cst
#endif // TTK_ENABLE_OPENMP4
        ttk::intgl::finishedElement_++;
      } else {
        this->storeToSendIfNecessary<triangulationType>(
          triangulation, integralLine, isMax);
      }
#endif
    } else {
      // If the integral line continues, the computation continues using lower
      // and upper components
      components = &(lowerComponents);
      if(direction_ == static_cast<int>(Direction::Forward)) {
        components = &(upperComponents);
      }
      if((criticalType == (char)CriticalType::Saddle1
          || criticalType == (char)CriticalType::Saddle2
          || criticalType == (char)CriticalType::Degenerate)
         && EnableForking) {
        // For each connected components, the max (or the min) is computed
        // and a task is created to further the computation of the integral line
        ttk::SimplexId const numberOfComponents = components->size();
#ifdef TTK_ENABLE_MPI
#ifdef TTK_ENABLE_OPENMP4
#pragma omp atomic update seq_cst
#endif // TTK_ENABLE_OPENMP4
        ttk::intgl::finishedElement_++;
#ifdef TTK_ENABLE_OPENMP4
#pragma omp atomic update seq_cst
#endif // TTK_ENABLE_OPENMP4
        ttk::intgl::addedElement_ += numberOfComponents;
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
          ttk::SimplexId const forkIdentifier = vnext;
#endif
          triangulation->getVertexPoint(vnext, p1[0], p1[1], p1[2]);
          double const distanceFork = Geometry::distance(p0, p1, 3);
          int threadNum{0};
#ifdef TTK_ENABLE_OPENMP4
          threadNum = omp_get_thread_num();
#endif // TTK_ENABLE_OPENMP4
          ttk::intgl::IntegralLine *integralLineFork
            = outputIntegralLines_->at(threadNum).addArrayElement(
              ttk::intgl::IntegralLine{
                std::vector<ttk::SimplexId>({v, vnext}),
                std::vector<double>({distance, distance + distanceFork}),
                std::vector<ttk::SimplexId>(
                  {integralLine->localVertexIdentifier.back(),
                   integralLine->localVertexIdentifier.back() + 1}),
                integralLine->seedIdentifier, forkIdentifier});

#ifdef TTK_ENABLE_OPENMP4
#pragma omp task firstprivate(integralLineFork)
          {
#endif // TTK_ENABLE_OPENMP4
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
#ifdef TTK_ENABLE_OPENMP4
          }
#endif // TTK_ENABLE_OPENMP4
        }
      } else {
        // In case the vertex is not a saddle point, all neighbor vertices
        // are used for the computation
        components->clear();
        ttk::SimplexId const neighborNumber
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
        integralLine->trajectory.push_back(vnext);
        p0[0] = p1[0];
        p0[1] = p1[1];
        p0[2] = p1[2];
        integralLine->distanceFromSeed.push_back(distance);
        integralLine->localVertexIdentifier.push_back(
          integralLine->localVertexIdentifier.back() + 1);
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
  std::vector<ttk::intgl::IntegralLine *> &chunkIntegralLine,
  int startingIndex,
  int nbElement,
  std::vector<SimplexId> *seeds) const {

  for(SimplexId j = 0; j < nbElement; j++) {
    SimplexId const v{seeds->at(j + startingIndex)};
    int seedIdentifier;
#ifdef TTK_ENABLE_MPI
    seedIdentifier = triangulation->getVertexGlobalId(v);
#else
    seedIdentifier = v;
#endif
    int threadNum{0};
#ifdef TTK_ENABLE_OPENMP4
    threadNum = omp_get_thread_num();
#endif // TTK_ENABLE_OPENMP4
    chunkIntegralLine[j] = outputIntegralLines_->at(threadNum).addArrayElement(
      ttk::intgl::IntegralLine{
        std::vector<ttk::SimplexId>(1, v), std::vector<double>(1, 0),
        std::vector<ttk::SimplexId>(1, 0), seedIdentifier, -1});
  }
}

template <class triangulationType>
void ttk::IntegralLines::createTask(
  const triangulationType *triangulation,
  std::vector<ttk::intgl::IntegralLine *> &chunkIntegralLine,
  const ttk::SimplexId *offsets,
  int nbElement) const {
#ifdef TTK_ENABLE_OPENMP4
#pragma omp task firstprivate(chunkIntegralLine)
  {
#endif // TTK_ENABLE_OPENMP4
    for(int j = 0; j < nbElement; j++) {
      this->computeIntegralLine<triangulationType>(
        triangulation, chunkIntegralLine[j], offsets);
    }
#ifdef TTK_ENABLE_OPENMP4
  }
#endif // TTK_ENABLE_OPENMP4
}

template <class triangulationType>
int ttk::IntegralLines::execute(triangulationType *triangulation) {

#ifdef TTK_ENABLE_MPI
  keepWorking_ = 1;
  ttk::intgl::finishedElement_ = 0;
  ttk::intgl::addedElement_ = 0;
#endif
  const SimplexId *offsets = inputOffsets_;
  std::vector<SimplexId> *seeds = vertexIdentifierScalarField_;
  Timer t;

  std::vector<ttk::intgl::IntegralLine *> chunkIntegralLine(chunkSize_);
  int const taskNumber = (int)seedNumber_ / chunkSize_;
#ifdef TTK_ENABLE_OPENMP4
#ifdef TTK_ENABLE_MPI
#pragma omp parallel shared(                                        \
  ttk::intgl::finishedElement_, toSend_, ttk::intgl::addedElement_) \
  num_threads(threadNumber_)
  {
#else
#pragma omp parallel num_threads(threadNumber_)
  {
#endif // TTK_ENABLE_MPI
#pragma omp master
    {
#endif // TTK_ENABLE_OPENMP4
      for(SimplexId i = 0; i < taskNumber; ++i) {
        this->prepareForTask<triangulationType>(
          triangulation, chunkIntegralLine, i * chunkSize_, chunkSize_, seeds);
        this->createTask<triangulationType>(
          triangulation, chunkIntegralLine, offsets, chunkSize_);
      }
      int const rest = seedNumber_ % chunkSize_;
      if(rest > 0) {
        this->prepareForTask<triangulationType>(
          triangulation, chunkIntegralLine, taskNumber * chunkSize_, rest,
          seeds);
        this->createTask<triangulationType>(
          triangulation, chunkIntegralLine, offsets, rest);
      }
#ifdef TTK_ENABLE_OPENMP4
    }
  }
#endif
#ifdef TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    int i;
    int finishedElementReceived = 0;
    std::vector<int> sendMessageSize(neighborNumber_);
    std::vector<int> recvMessageSize(neighborNumber_);
    std::vector<std::vector<intgl::ElementToBeSent>> send_buf(neighborNumber_);
    std::vector<std::vector<intgl::ElementToBeSent>> recv_buf(neighborNumber_);
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
      ttk::intgl::finishedElement_ -= ttk::intgl::addedElement_;
      ttk::intgl::addedElement_ = 0;
      // Exchange of the number of integral lines finished on all processes
      MPI_Allreduce(&ttk::intgl::finishedElement_, &finishedElementReceived, 1,
                    MPI_INTEGER, MPI_SUM, ttk::MPIcomm_);
      ttk::intgl::finishedElement_ = 0;
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
                    INTEGRAL_LINE_IS_MESSAGE_SIZE, ttk::MPIcomm_,
                    &requests[2 * i]);
          MPI_Irecv(&recvMessageSize[i], 1, MPI_INTEGER, neighbors_.at(i),
                    INTEGRAL_LINE_IS_MESSAGE_SIZE, ttk::MPIcomm_,
                    &requests[2 * i + 1]);
        }
        MPI_Waitall(2 * neighborNumber_, requests.data(), MPI_STATUSES_IGNORE);
        // Exchange of the data
        for(i = 0; i < neighborNumber_; i++) {
          if(recv_buf[i].size() < (size_t)recvMessageSize[i]) {
            recv_buf[i].resize(recvMessageSize[i]);
          }
          if(recvMessageSize[i] > 0) {
            MPI_Irecv(recv_buf[i].data(), recvMessageSize[i], this->MessageType,
                      neighbors_.at(i), INTEGRAL_LINE_IS_ELEMENT_TO_PROCESS,
                      ttk::MPIcomm_, &requests[2 * i]);
            totalMessageSize += recvMessageSize[i];
          }

          if(sendMessageSize[i] > 0) {
            MPI_Isend(send_buf[i].data(), sendMessageSize[i], this->MessageType,
                      neighbors_.at(i), INTEGRAL_LINE_IS_ELEMENT_TO_PROCESS,
                      ttk::MPIcomm_, &requests[2 * i + 1]);
          }
        }
        MPI_Waitall(2 * neighborNumber_, requests.data(), MPI_STATUSES_IGNORE);
        for(i = 0; i < neighborNumber_; i++) {
          send_buf[i].clear();
        }
        // Extraction of the received data and creation of the tasks
#ifdef TTK_ENABLE_OPENMP4
#pragma omp parallel shared(ttk::intgl::finishedElement_, toSend_) \
  num_threads(threadNumber_)
        {
#pragma omp master
          {
#endif // TTK_ENABLE_OPENMP4
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
#ifdef TTK_ENABLE_OPENMP4
          }
        }
#endif // TTK_ENABLE_OPENMP4
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

#ifdef TTK_ENABLE_MPI
template <typename triangulationType>
int ttk::IntegralLines::getGlobalIdentifiers(
  std::vector<ttk::SimplexId> &globalVertexId,
  std::vector<ttk::SimplexId> &globalCellId,
  const std::vector<
    ttk::ArrayLinkedList<ttk::intgl::IntegralLine, INTEGRAL_LINE_TABULAR_SIZE>>
    &integralLines,
  const triangulationType *triangulation) {
  ttk::SimplexId outputVertexNumber = 0;
  ttk::SimplexId outputCellNumber = 0;
  ttk::SimplexId realVertexNumber = 0;
  ttk::SimplexId realCellNumber = 0;
  ttk::SimplexId intervalSize;
  // Counts vertices and edges number (with and without ghosts)
#ifdef TTK_ENABLE_OPENMP4
#pragma omp parallel for reduction(+:outputVertexNumber,outputCellNumber,realCellNumber,realVertexNumber) schedule(static,1) private(intervalSize)
#endif
  for(int thread = 0; thread < threadNumber_; thread++) {
    std::list<std::array<ttk::intgl::IntegralLine,
                         INTEGRAL_LINE_TABULAR_SIZE>>::const_iterator
      integralLine
      = integralLines[thread].list_.begin();
    while(integralLine != integralLines[thread].list_.end()) {
      for(int i = 0; i < INTEGRAL_LINE_TABULAR_SIZE; i++) {
        if(integralLine->at(i).trajectory.size() > 0) {
          intervalSize = static_cast<ttk::SimplexId>(
            integralLine->at(i).trajectory.size());
          outputVertexNumber += intervalSize;
          outputCellNumber += intervalSize - 1;
          if(triangulation->getVertexRank(integralLine->at(i).trajectory.at(0))
             != ttk::MPIrank_) {
            intervalSize--;
          }
          if(integralLine->at(i).trajectory.size() > 1) {
            realCellNumber += intervalSize - 1;
            if(triangulation->getVertexRank(
                 integralLine->at(i).trajectory.back())
               != ttk::MPIrank_) {
              intervalSize--;
            }
          }
          realVertexNumber += intervalSize;
        } else {
          break;
        }
      }
      integralLine++;
    }
  }

  ttk::SimplexId vertIndex;
  ttk::SimplexId cellIndex;

  // Perform exclusive prefix sum to find local offset for vertices and
  // cells
  MPI_Exscan(&realVertexNumber, &vertIndex, 1, ttk::getMPIType(vertIndex),
             MPI_SUM, ttk::MPIcomm_);
  MPI_Exscan(&realCellNumber, &cellIndex, 1, ttk::getMPIType(cellIndex),
             MPI_SUM, ttk::MPIcomm_);

  // Rank 0 received garbage values, it is replaced by the correct offset
  // (always 0)
  if(ttk::MPIrank_ == 0) {
    vertIndex = 0;
    cellIndex = 0;
  }

  // Generate Global identifers and package ghost data for the process that owns
  // the ghost data
  ttk::SimplexId startCellId = 0;
  ttk::SimplexId startVertexId = 0;
  globalVertexId.resize(outputVertexNumber);
  globalCellId.resize(outputCellNumber);
  std::vector<std::vector<ttk::intgl::GhostElementsToSort>> unmatchedGhosts(
    neighborNumber_);
  std::vector<std::vector<ttk::SimplexId>> unmatchedGhostsVertexLocalId(
    neighborNumber_);
  std::vector<std::vector<ttk::SimplexId>> unmatchedGhostsEdgeLocalId(
    neighborNumber_);
  for(int thread = 0; thread < threadNumber_; thread++) {
    std::list<std::array<ttk::intgl::IntegralLine,
                         INTEGRAL_LINE_TABULAR_SIZE>>::const_iterator
      integralLine
      = integralLines[thread].list_.begin();
    while(integralLine != integralLines[thread].list_.end()) {
      for(int i = 0; i < INTEGRAL_LINE_TABULAR_SIZE; i++) {
        if(integralLine->at(i).trajectory.size() > 0) {
          if(triangulation->getVertexRank(integralLine->at(i).trajectory.at(0))
             != ttk::MPIrank_) {
            globalVertexId.at(startVertexId) = -1;
          } else {
            globalVertexId.at(startVertexId) = vertIndex;
            vertIndex++;
          }
          startVertexId++;
          if(integralLine->at(i).trajectory.size() > 1) {
            if(triangulation->getVertexRank(
                 integralLine->at(i).trajectory.at(0))
               != ttk::MPIrank_) {
              globalCellId.at(startCellId) = -1;
              unmatchedGhosts
                .at(neighborsToId_[triangulation->getVertexRank(
                  integralLine->at(i).trajectory.at(0))])
                .push_back(ttk::intgl::GhostElementsToSort{
                  vertIndex, integralLine->at(i).localVertexIdentifier.at(0),
                  integralLine->at(i).seedIdentifier,
                  integralLine->at(i).forkIdentifier, -1, startVertexId - 1,
                  startCellId});
            } else {
              globalCellId.at(startCellId) = cellIndex;
              cellIndex++;
            }
            startCellId++;
            if(integralLine->at(i).trajectory.size() > 2) {
              std::iota(globalCellId.begin() + startCellId,
                        globalCellId.begin() + startCellId
                          + integralLine->at(i).trajectory.size() - 2,
                        cellIndex);
              std::iota(globalVertexId.begin() + startVertexId,
                        globalVertexId.begin() + startVertexId
                          + integralLine->at(i).trajectory.size() - 2,
                        vertIndex);
              startCellId += integralLine->at(i).trajectory.size() - 2;
              startVertexId += integralLine->at(i).trajectory.size() - 2;
              vertIndex += integralLine->at(i).trajectory.size() - 2;
              cellIndex += integralLine->at(i).trajectory.size() - 2;
            }
            if(triangulation->getVertexRank(
                 integralLine->at(i).trajectory.back())
               != ttk::MPIrank_) {
              globalVertexId.at(startVertexId) = -1;
              unmatchedGhosts
                .at(neighborsToId_[triangulation->getVertexRank(
                  integralLine->at(i).trajectory.back())])
                .push_back(ttk::intgl::GhostElementsToSort{
                  globalVertexId.at(startVertexId - 1),
                  integralLine->at(i).localVertexIdentifier.at(
                    integralLine->at(i).localVertexIdentifier.size() - 2),
                  integralLine->at(i).seedIdentifier,
                  integralLine->at(i).forkIdentifier,
                  globalCellId.at(startCellId - 1), startVertexId,
                  startCellId - 1});
            } else {
              globalVertexId.at(startVertexId) = vertIndex;
              vertIndex++;
            }
            startVertexId++;
          }
        } else {
          break;
        }
      }
      integralLine++;
    }
  }
  // unmatchedGhosts contains the ghost data for each process.
  // For two processes i and j that are neighbors, their respective
  // unmatchedGhosts[k_i] and unmatchedGhosts[k_j] have the same
  // number of elements and each element represents a segment of
  // integral line. This means that if unmatchedGhosts[k_i] and
  // unmatchedGhosts[k_j] are sorted using the same comparator,
  // element l of unmatchedGhosts[k_i] and element l of unmatchedGhosts[k_j]
  // represent the same segment of integral line. Therefore, once
  // unmatchedGhosts vectors are sorted for each process, all that
  // is left to do is exchange the global identifiers of vertices
  // and edges, as the receiving process already knows to which local
  // vertices in its domain the ghosts corresponds to.
  if(ttk::MPIsize_ > 1) {
    for(int i = 0; i < neighborNumber_; i++) {
      TTK_PSORT(threadNumber_, unmatchedGhosts.at(i).begin(),
                unmatchedGhosts.at(i).end());
    }
    this->exchangeGhosts(unmatchedGhosts, globalVertexId, globalCellId);
  }

  return 0;
}

inline int ttk::IntegralLines::exchangeGhosts(
  const std::vector<std::vector<ttk::intgl::GhostElementsToSort>>
    &unmatchedGhosts,
  std::vector<ttk::SimplexId> &globalVertexId,
  std::vector<ttk::SimplexId> &globalCellId) {
  ttk::SimplexId id{0};
  std::vector<std::vector<ttk::SimplexId>> globalIdsToSend(neighborNumber_);
  std::vector<std::vector<ttk::SimplexId>> globalIdsToReceive(neighborNumber_);
  for(int i = 0; i < neighborNumber_; i++) {
    globalIdsToSend[i].resize(2 * unmatchedGhosts[i].size());
    globalIdsToReceive[i].resize(2 * unmatchedGhosts[i].size());
  }
  // The sending buffer is prepared
#ifdef TTK_ENABLE_OPENMP4
#pragma omp parallel for
#endif
  for(int i = 0; i < neighborNumber_; i++) {
    for(size_t j = 0; j < unmatchedGhosts[i].size(); j++) {
      globalIdsToSend[i][2 * j] = unmatchedGhosts[i][j].ownedGlobalId;
      globalIdsToSend[i][2 * j + 1] = unmatchedGhosts[i][j].globalEdgeId;
    }
  }
  // Each process sends and receives global identifiers of ghosts
  // to its neighbors
  for(int neighbor : neighbors_) {
    MPI_Sendrecv(globalIdsToSend[neighborsToId_[neighbor]].data(),
                 globalIdsToSend[neighborsToId_[neighbor]].size(),
                 ttk::getMPIType(id), neighbor, ttk::MPIrank_,
                 globalIdsToReceive[neighborsToId_[neighbor]].data(),
                 globalIdsToSend[neighborsToId_[neighbor]].size(),
                 ttk::getMPIType(id), neighbor, neighbor, ttk::MPIcomm_,
                 MPI_STATUS_IGNORE);
  }
  // The global identifiers of ghosts are inserted in the globalVertexId
  // and globalCellId vectors
#ifdef TTK_ENABLE_OPENMP4
#pragma omp parallel for
#endif
  for(int i = 0; i < neighborNumber_; i++) {
    for(size_t j = 0; j < unmatchedGhosts[i].size(); j++) {
      globalVertexId.at(unmatchedGhosts[i][j].ghostVertexLocalId)
        = globalIdsToReceive[i][2 * j];
      if(globalCellId.at(unmatchedGhosts[i][j].ghostEdgeLocalId) == -1) {
        globalCellId.at(unmatchedGhosts[i][j].ghostEdgeLocalId)
          = globalIdsToReceive[i][j * 2 + 1];
      }
    }
  }
  return 0;
}

#endif // TTK_ENABLE_MPI
