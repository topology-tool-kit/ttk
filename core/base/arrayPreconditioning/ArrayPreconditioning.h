/// \ingroup base
/// \class ttk::ArrayPreconditioning
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date 2022.
///
/// This module defines the %ArrayPreconditioning class that generates order
/// arrays from a selection of scalar field arrays.
/// In distributed, using GlobalOrder set to True, this module will compute a
/// global order, otherwise each process will locally compute its order.
///
/// \b Online \b examples:\n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mpiExample/">
///   MPI example</a> \n
#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <psort.h>
#include <vector>

namespace ttk {

#ifdef TTK_ENABLE_MPI

#ifdef TTK_ENABLE_MPI_RANK_ID_INT
  using RankId = int;
#else
  using RankId = char;
#endif

  namespace globalOrder {

    template <typename datatype>
    struct vertexToSort {
      ttk::SimplexId globalId;
      datatype value;
      ttk::RankId rank;
    };

    struct sortedVertex {
      ttk::SimplexId globalId;
      ttk::SimplexId order;
    };

    template <typename datatype>
    bool comp(const vertexToSort<datatype> a, const vertexToSort<datatype> b) {
      return (b.value > a.value)
             || (a.value == b.value && a.globalId < b.globalId);
    };

    template <typename datatype>
    bool oppositeComp(const vertexToSort<datatype> a,
                      const vertexToSort<datatype> b) {
      return (a.value > b.value)
             || (a.value == b.value && a.globalId > b.globalId);
    }
  } // namespace globalOrder

#endif // TTK_ENABLE_MPI
  /**
   * The ArrayPreconditioning class provides methods to generate order arrays
   * from a selection of scalar field arrays.
   */
  class ArrayPreconditioning : virtual public Debug {

  public:
    ArrayPreconditioning();

    int preconditionTriangulation(AbstractTriangulation *triangulation) {
      // Pre-condition functions.
      if(triangulation) {
#ifdef TTK_ENABLE_MPI
        triangulation->preconditionExchangeGhostVertices();
#endif // TTK_ENABLE_MPI
      }
      return 0;
    }

    void setGlobalOrder(bool order) {
      this->GlobalOrder = order;
    }

#ifdef TTK_ENABLE_MPI
    /**
     * This function computes the post-processing steps of the
     * algorithm: computing the global order of a vertex using the offset
     * of the current process and the local order of the vertex, sending and
     * receiving the order of vertices than the process has computed or owns and
     * placing the received global order in the orderArray array. When more than
     * ChunkSize elements are being handled on one process at once, the
     * post-processing is done in chunks, that starts at the index begin and
     * ends at index end. In the case where some processes will send more
     * messages than others, only the send/recv step will be computed.
     */
    template <typename DT, typename triangulationType>
    int
      postprocessing(const triangulationType *triangulation,
                     ttk::SimplexId *orderArray,
                     const ttk::SimplexId begin,
                     const ttk::SimplexId end,
                     const ttk::SimplexId orderOffset,
                     std::vector<globalOrder::vertexToSort<DT>> &verticesToSort,
                     MPI_Datatype &MPI_sortedVertexType,
                     MPI_Datatype &MPI_SimplexId,
                     const ttk::SimplexId itNumber,
                     const ttk::SimplexId localNbChunk) const {
      ttk::RankId rank{0};
      std::vector<std::vector<globalOrder::sortedVertex>> verticesSorted(
        ttk::MPIsize_, std::vector<globalOrder::sortedVertex>());
      if(itNumber < localNbChunk + 1) {
        // Compute the order and prepare the vector to send
#ifdef TTK_ENABLE_OPENMP
        std::vector<std::vector<std::vector<globalOrder::sortedVertex>>>
          verticesSortedThread(
            this->threadNumber_,
            std::vector<std::vector<globalOrder::sortedVertex>>(
              ttk::MPIsize_, std::vector<globalOrder::sortedVertex>()));
#pragma omp parallel firstprivate(rank) num_threads(threadNumber_)
        {
          int threadNumber = omp_get_thread_num();
#pragma omp for schedule(static)
          for(ttk::SimplexId i = begin; i < end; i++) {
            rank = verticesToSort.at(i).rank;
            if(rank == ttk::MPIrank_) {
              orderArray[triangulation->getVertexLocalId(
                verticesToSort.at(i).globalId)]
                = orderOffset + i;
            } else {
              verticesSortedThread.at(threadNumber)
                .at(rank)
                .push_back(globalOrder::sortedVertex{
                  verticesToSort.at(i).globalId, orderOffset + i});
            }
          }
        }
        // For a smaller memory footprint
        verticesToSort.resize(begin + 1);
        verticesToSort.shrink_to_fit();
        // Concatenate the vector produced by each thread
#pragma omp parallel for schedule(static, 1)
        for(int j = 0; j < ttk::MPIsize_; j++) {
          for(int i = 0; i < this->threadNumber_; i++) {
            if(j != ttk::MPIrank_) {
              verticesSorted.at(j).insert(
                verticesSorted.at(j).end(),
                verticesSortedThread.at(i).at(j).begin(),
                verticesSortedThread.at(i).at(j).end());
            }
          }
        }
        verticesSortedThread.clear();
        verticesSortedThread.shrink_to_fit();
#else
        for(ttk::SimplexId i = begin; i < end; i++) {
          rank = verticesToSort.at(i).rank;
          if(rank == ttk::MPIrank_) {
            orderArray[triangulation->getVertexLocalId(
              verticesToSort.at(i).globalId)]
              = orderOffset + i;
          } else {
            verticesSorted.at(rank).push_back(globalOrder::sortedVertex{
              verticesToSort.at(i).globalId, orderOffset + i});
          }
        }
#endif // TTK_ENABLE_OPENMP
      }
      // Send and receive the global orders. The use of MPI_Waitsome,
      // Isend and Irecv enables the computation to overlap communications.
      std::vector<MPI_Request> sendRequests(ttk::MPIsize_ - 1);
      std::vector<MPI_Request> recvRequests(ttk::MPIsize_ - 1);
      std::vector<MPI_Status> sendStatus(ttk::MPIsize_ - 1);
      std::vector<MPI_Status> recvStatus(ttk::MPIsize_ - 1);
      std::vector<ttk::SimplexId> sendMessageSize(ttk::MPIsize_, 0);
      std::vector<ttk::SimplexId> recvMessageSize(ttk::MPIsize_, 0);
      std::vector<int> recvCompleted(ttk::MPIsize_ - 1, 0);
      std::vector<int> sendCompleted(ttk::MPIsize_ - 1, 0);
      int sendPerformedCount = 0;
      int recvPerformedCount = 0;
      int sendPerformedCountTotal = 0;
      int recvPerformedCountTotal = 0;
      int count = 0;
      for(int i = 0; i < ttk::MPIsize_; i++) {
        // Send size of verticesToSort.at(i)
        if(i != ttk::MPIrank_) {
          sendMessageSize[i] = verticesSorted.at(i).size();
          MPI_Isend(&sendMessageSize[i], 1, MPI_SimplexId, i, 0, ttk::MPIcomm_,
                    &sendRequests[count]);
          MPI_Irecv(&recvMessageSize[i], 1, MPI_SimplexId, i, 0, ttk::MPIcomm_,
                    &recvRequests[count]);
          count++;
        }
      }
      std::vector<std::vector<globalOrder::sortedVertex>> recvVerticesSorted(
        ttk::MPIsize_, std::vector<globalOrder::sortedVertex>());
      std::vector<MPI_Request> sendRequestsData(ttk::MPIsize_ - 1);
      std::vector<MPI_Request> recvRequestsData(ttk::MPIsize_ - 1);
      std::vector<MPI_Status> recvStatusData(ttk::MPIsize_ - 1);
      int recvCount = 0;
      int sendCount = 0;
      int r;
      while((sendPerformedCountTotal < ttk::MPIsize_ - 1
             || recvPerformedCountTotal < ttk::MPIsize_ - 1)) {
        if(sendPerformedCountTotal < ttk::MPIsize_ - 1) {
          MPI_Waitsome(ttk::MPIsize_ - 1, sendRequests.data(),
                       &sendPerformedCount, sendCompleted.data(),
                       sendStatus.data());
          if(sendPerformedCount > 0) {
            for(int i = 0; i < sendPerformedCount; i++) {
              r = sendCompleted[i];
              if(ttk::MPIrank_ <= sendCompleted[i]) {
                r++;
              }
              if((sendMessageSize[r] > 0)) {
                MPI_Isend(verticesSorted.at(r).data(), sendMessageSize[r],
                          MPI_sortedVertexType, r, 1, ttk::MPIcomm_,
                          &sendRequestsData[sendCount]);
                sendCount++;
              }
            }
            sendPerformedCountTotal += sendPerformedCount;
          }
        }
        if(recvPerformedCountTotal < ttk::MPIsize_ - 1) {
          MPI_Waitsome(ttk::MPIsize_ - 1, recvRequests.data(),
                       &recvPerformedCount, recvCompleted.data(),
                       recvStatus.data());
          if(recvPerformedCount > 0) {
            for(int i = 0; i < recvPerformedCount; i++) {
              r = recvStatus[i].MPI_SOURCE;
              if((recvMessageSize[r] > 0)) {
                recvVerticesSorted.at(r).resize(recvMessageSize[r]);
                MPI_Irecv(recvVerticesSorted.at(r).data(), recvMessageSize[r],
                          MPI_sortedVertexType, r, 1, ttk::MPIcomm_,
                          &recvRequestsData[recvCount]);

                recvCount++;
              }
            }
            recvPerformedCountTotal += recvPerformedCount;
          }
        }
      }
      recvPerformedCountTotal = 0;
      while(recvPerformedCountTotal < recvCount) {
        MPI_Waitsome(recvCount, recvRequestsData.data(), &recvPerformedCount,
                     recvCompleted.data(), recvStatusData.data());
        if(recvPerformedCount > 0) {
          for(int i = 0; i < recvPerformedCount; i++) {
            r = recvStatusData[i].MPI_SOURCE;
#pragma omp parallel for schedule(static)
            for(int j = 0; j < recvMessageSize[r]; j++) {
              orderArray[triangulation->getVertexLocalId(
                recvVerticesSorted.at(r).at(j).globalId)]
                = recvVerticesSorted.at(r).at(j).order;
            }
          }
          recvPerformedCountTotal += recvPerformedCount;
        }
      }
      MPI_Waitall(sendCount, sendRequestsData.data(), MPI_STATUSES_IGNORE);
      return 0;
    }

#endif // TTK_ENABLE_MPI

    template <typename DT, typename triangulationType>
    int processScalarArray(const triangulationType *triangulation,
                           ttk::SimplexId *orderArray,
                           const DT *scalarArray,
                           const size_t nVerts) const { // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator
      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(nVerts)},
      });
      this->printMsg(ttk::debug::Separator::L1);

// -----------------------------------------------------------------------
// Computing order Array
// -----------------------------------------------------------------------
#ifdef TTK_ENABLE_MPI
      if(ttk::isRunningWithMPI()) {
        ttk::SimplexId id = 0;
        MPI_Datatype MPI_SimplexId = getMPIType(id);
        // Prepare the vector to sort.
        std::vector<globalOrder::vertexToSort<DT>> verticesToSort;
        verticesToSort.reserve(nVerts);
#ifdef TTK_ENABLE_OPENMP
#pragma omp declare reduction (merge : std::vector<globalOrder::vertexToSort<DT>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp parallel for reduction(merge : verticesToSort) schedule(static)
#endif
        for(size_t i = 0; i < nVerts; i++) {
          if(triangulation->getVertexRank(i) == ttk::MPIrank_) {
            verticesToSort.emplace_back(globalOrder::vertexToSort<DT>{
              triangulation->getVertexGlobalId(i), scalarArray[i],
              (ttk::RankId)ttk::MPIrank_});
          }
        }

        MPI_Datatype MPI_vertexToSortType;

        /*
         *  WARNING: the struct is sent as an array of char, as experiments show
         * that using MPI's built-in struct management yields poor performance
         * when used with a templated struct.
         */
        MPI_Type_contiguous(sizeof(globalOrder::vertexToSort<DT>), MPI_CHAR,
                            &MPI_vertexToSortType);
        MPI_Type_commit(&MPI_vertexToSortType);

        std::vector<ttk::SimplexId> vertexDistribution(ttk::MPIsize_);
        ttk::SimplexId localVertexNumber = verticesToSort.size();
        MPI_Allgather(&localVertexNumber, 1, MPI_SimplexId,
                      vertexDistribution.data(), 1, MPI_SimplexId,
                      ttk::MPIcomm_);
        // The distributed sort will sort using the scalar field first,
        // and the global id for disambiguation.
        p_sort::parallel_sort<globalOrder::vertexToSort<DT>>(
          verticesToSort, globalOrder::comp<DT>, globalOrder::oppositeComp<DT>,
          vertexDistribution, MPI_vertexToSortType, MPI_SimplexId,
          threadNumber_);

        // Now, the computation for the global order starts, using the sorted
        // vector.

        MPI_Datatype types[] = {MPI_SimplexId, MPI_SimplexId};
        int lengths[] = {1, 1};
        MPI_Datatype MPI_sortedVertexType;
        const long int mpi_offsets[]
          = {offsetof(globalOrder::sortedVertex, globalId),
             offsetof(globalOrder::sortedVertex, order)};
        MPI_Type_create_struct(
          2, lengths, mpi_offsets, types, &MPI_sortedVertexType);
        MPI_Type_commit(&MPI_sortedVertexType);

        ttk::SimplexId verticesToSortSize = verticesToSort.size();
        // Compute the order of the first element of the current process
        ttk::SimplexId orderOffset
          = std::accumulate(vertexDistribution.begin(),
                            vertexDistribution.begin() + ttk::MPIrank_, 0);
        // nbChunk, rest and nbChunkTotal are used to compute the
        // post-processing bit by bit.
        ttk::SimplexId nbChunk = std::floor(verticesToSortSize / ChunkSize);
        ttk::SimplexId rest = verticesToSortSize % ChunkSize;
        ttk::SimplexId nbChunkTotal{nbChunk};
        if(rest > 0) {
          nbChunkTotal++;
        }
        MPI_Allreduce(MPI_IN_PLACE, &nbChunkTotal, 1, MPI_SimplexId, MPI_MAX,
                      ttk::MPIcomm_);
        ttk::SimplexId beginningChunk = nbChunk;

        // The end of the vector is computed first to resize and delete easily
        // elements of verticesToSort once they will no longer be used in the
        // rest of the algorithm.
        if(rest > 0) {
          beginningChunk++;
          postprocessing<DT, triangulationType>(
            triangulation, orderArray, nbChunk * ChunkSize,
            nbChunk * ChunkSize + rest, orderOffset, verticesToSort,
            MPI_sortedVertexType, MPI_SimplexId, nbChunk, nbChunk);
        }
        for(int chunk = nbChunk - 1; chunk >= 0; chunk--) {
          postprocessing<DT, triangulationType>(
            triangulation, orderArray, chunk * ChunkSize,
            (chunk + 1) * ChunkSize, orderOffset, verticesToSort,
            MPI_sortedVertexType, MPI_SimplexId, chunk, nbChunk);
        }

        // In case some processes still have messages to send after the current
        // process has sent all of his.
        for(int chunk = beginningChunk; chunk < nbChunkTotal; chunk++) {
          postprocessing<DT, triangulationType>(
            triangulation, orderArray, nbChunk * ChunkSize + rest,
            nbChunk * ChunkSize + rest, orderOffset, verticesToSort,
            MPI_sortedVertexType, MPI_SimplexId, chunk, nbChunk);
        }

        // Exchange of the order array computed for the ghost vertices.
        ttk::exchangeGhostVertices<ttk::SimplexId, triangulationType>(
          orderArray, triangulation, ttk::MPIcomm_, 1);
      }
#else
      this->printMsg("MPI not enabled!");
      TTK_FORCE_USE(orderArray);
      TTK_FORCE_USE(scalarArray);
      TTK_FORCE_USE(triangulation);
      return 0;
#endif // TTK_ENABLE_MPI

      // ---------------------------------------------------------------------
      // print global performance
      // ---------------------------------------------------------------------
      {
        this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime() // global progress, time
        );
        this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
      }

      return 1; // return success
    }

  protected:
    bool GlobalOrder{false};
    // This value has been chosen for systems of 128 Gb of memory per computing
    // node. For systems with much smaller memory, it may be inadequate and
    // require a smaller value.
    int ChunkSize{1000000000};
  }; // ArrayPreconditioning class

} // namespace ttk
