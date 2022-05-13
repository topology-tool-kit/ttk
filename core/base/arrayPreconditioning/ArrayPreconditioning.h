/// \ingroup base
/// \class ttk::ArrayPreconditioning
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date 2022.
///
/// This module defines the %ArrayPreconditioning class that generates order
/// arrays from a selection of scalar field arrays.
///

#pragma once

// ttk common includes
#include <Debug.h>

#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
namespace ttk {

  /**
   * The ArrayPreconditioning class provides methods to generate order arrays
   * from a selection of scalar field arrays.
   */
  class ArrayPreconditioning : virtual public Debug {

  public:
    ArrayPreconditioning();

    template <typename DT, typename IT>
    struct value {
      DT scalar;
      IT globalId;

      value(DT _scalar, IT _globalId) : scalar(_scalar), globalId(_globalId) {
      }
    };

    // creates a vector of value structs based on the given pointers
    // only takes values of vertices which mainly belong to the current rank
    // ( so only vertices which are no ghost cells)
    template <typename DT, typename IT>
    void populateVector(std::vector<value<DT, IT>> &valuesToSortVector,
                        std::vector<IT> &gidsToGetVector,
                        std::unordered_map<IT, IT> &gidToLidMap,
                        const size_t nVerts,
                        const DT *const scalars,
                        const IT *const globalIds,
                        const int *const rankArray,
                        const int rank) const {
      for(size_t i = 0; i < nVerts; i++) {
        IT globalId = globalIds[i];
        gidToLidMap[globalId] = i;
        if(rankArray[i] == rank) {
          valuesToSortVector.emplace_back(scalars[i], globalId);
        } else {
          gidsToGetVector.push_back(globalId);
        }
      }
    }

    // orders an value vector first by their scalar value and then by global id
    template <typename DT, typename IT>
    void sortVerticesDistributed(std::vector<value<DT, IT>> &values,
                                 const int nThreads) const {

      TTK_FORCE_USE(nThreads);

      TTK_PSORT(
        nThreads, values.begin(), values.end(),
        [](value<DT, IT> v1, value<DT, IT> v2) {
          return (v1.scalar < v2.scalar)
                 || (v1.scalar == v2.scalar && v1.globalId < v2.globalId);
        });
    }

    // send the highest burstSize values and decrease the vector by that amount
    // check if there are actually that many elements in the vector
    template <typename DT, typename IT>
    void returnVectorForBurstsize(std::vector<value<DT, IT>> &outVector,
                                  std::vector<value<DT, IT>> &values,
                                  size_t burstSize) const {

      if(burstSize > values.size()) {
        outVector.resize(values.size(), {0, 0});
        outVector.assign(values.begin(), values.end());
        values.clear();
      } else {
        outVector.resize(burstSize, {0, 0});
        outVector.assign(values.end() - burstSize, values.end());
        values.erase(values.end() - burstSize, values.end());
      }
    }

    /**
     * @brief Sort vertices according to scalars disambiguated by offsets
     *
     * @param[in] nInts number of long ints, as globalid / order pairs
     * @param[in] orderedValuesForRank array of size nInts, the actual pairs,
     * [i] = gId, [i+1] = order
     * @param[in] gidToLidMap map which maps scalar values to a defined order
     * @param[out] order array of size nVerts, computed order of vertices, this
     * procedure doesn't fill it completely
     * @param[in] nThreads number of parallel threads
     */
    template <typename IT>
    void buildArrayForReceivedData(const size_t nInts,
                                   const IT *const orderedValuesForRank,
                                   std::unordered_map<IT, IT> &gidToLidMap,
                                   SimplexId *const order,
                                   const int nThreads) const {

      TTK_FORCE_USE(nThreads);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nInts; i += 2) {
        order[gidToLidMap[orderedValuesForRank[i]]]
          = orderedValuesForRank[i + 1];
      }
    }

#ifdef TTK_ENABLE_MPI
    // this method gets the current maximum value and works with it
    // if this leads to the vector of one rank being empty,
    // it sends the ordering to the rank and requests new values
    template <typename DT, typename IT>
    void getMax(int intTag,
                int structTag,
                IT &currentOrder,
                int burstSize,
                MPI_Datatype MPI_IT,
                std::vector<value<DT, IT>> &finalValues,
                std::vector<std::vector<value<DT, IT>>> &unsortedReceivedValues,
                std::vector<std::vector<IT>> &orderResendValues,
                std::vector<IT> &orderedValuesForRank,
                std::vector<value<DT, IT>> &sortingValues) const {
      // take the current maximum scalar over all ranks
      int rankIdOfMaxScalar = -1;
      DT maxScalar = std::numeric_limits<DT>::lowest();
      IT maxGId = -1;
      for(size_t i = 0; i < unsortedReceivedValues.size(); i++) {
        if(unsortedReceivedValues[i].size() > 0) {
          const auto &v = unsortedReceivedValues[i].back();
          if(v.scalar == maxScalar ? v.globalId > maxGId
                                   : v.scalar > maxScalar) {
            maxScalar = v.scalar;
            maxGId = v.globalId;
            rankIdOfMaxScalar = i;
          }
        }
      }
      if(rankIdOfMaxScalar == -1) {

        this->printMsg("FinalValues.size: "
                       + std::to_string(finalValues.size()));
        this->printErr("All vectors are empty, but out final vector is "
                       "not complete yet. Either something went wrong or "
                       "some rank didn't send their values yet.");
      }
      // move the struct from the unsortedReceivedValues subvector to the
      // finalValues vector to get an ordering
      value<DT, IT> currentValue
        = unsortedReceivedValues[rankIdOfMaxScalar].back();
      // we send the globalId and the the order via one send command,
      // therefore we need to check two concurrent values at the same time
      // later on
      orderResendValues[rankIdOfMaxScalar].push_back(currentValue.globalId);
      orderResendValues[rankIdOfMaxScalar].push_back(currentOrder);
      currentOrder--;
      finalValues.push_back(currentValue);
      unsortedReceivedValues[rankIdOfMaxScalar].pop_back();
      if(unsortedReceivedValues[rankIdOfMaxScalar].size() == 0) {
        if(rankIdOfMaxScalar == 0) {
          // append the ordered values to the correct vector
          orderedValuesForRank.insert(
            orderedValuesForRank.end(),
            orderResendValues[rankIdOfMaxScalar].begin(),
            orderResendValues[rankIdOfMaxScalar].end());
          orderResendValues[rankIdOfMaxScalar].clear();

          if(sortingValues.size() > 0) {
            std::vector<value<DT, IT>> ownValues;
            this->returnVectorForBurstsize<DT, IT>(
              ownValues, sortingValues, burstSize);

            unsortedReceivedValues[rankIdOfMaxScalar] = ownValues;
          } else {
            this->printMsg("We are done with rank 0");
          }
        } else {
          // receive more values from rank, send ordering to the rank
          // send to the finished rank that we want more
          MPI_Send(orderResendValues[rankIdOfMaxScalar].data(),
                   orderResendValues[rankIdOfMaxScalar].size(), MPI_IT,
                   rankIdOfMaxScalar, intTag * rankIdOfMaxScalar,
                   MPI_COMM_WORLD);
          orderResendValues[rankIdOfMaxScalar].clear();
          // we receive more vals, if the size is still zero afterwards, we are
          // finished with this rank
          this->ReceiveAndAddToVector(
            rankIdOfMaxScalar, structTag, unsortedReceivedValues);
          if(unsortedReceivedValues[rankIdOfMaxScalar].size() == 0) {
            this->printMsg("We are done with rank "
                           + std::to_string(rankIdOfMaxScalar));
          }
        }
      }
    }
#endif

#ifdef TTK_ENABLE_MPI
    template <typename DT, typename IT>
    void ReceiveAndAddToVector(
      int rankFrom,
      int structTag,
      std::vector<std::vector<value<DT, IT>>> &unsortedReceivedValues) const {
      std::vector<value<DT, IT>> &receivedValues
        = unsortedReceivedValues[rankFrom];
      // be prepared to receive burstsize of elements, resize after receiving to
      // the correct size
      int amount;
      MPI_Status status;
      MPI_Probe(rankFrom, structTag * rankFrom, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_CHAR, &amount);
      receivedValues.resize(amount / sizeof(value<DT, IT>), {0, 0});
      MPI_Recv(receivedValues.data(), amount, MPI_CHAR, rankFrom,
               structTag * rankFrom, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
#endif

    template <typename DT, typename IT>
    int processScalarArray(ttk::SimplexId *orderArray,
                           const DT *scalarArray,
                           const IT *globalIds,
                           const int *rankArray,
                           const size_t nVerts,
                           const int burstSize) const { // start global timer
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
      {
        int numProcs;
        int rank;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int intTag = 101;
        int structTag = 102;
        if(rank == 0)
          this->printMsg("Global Point Ids and RankArray exist, "
                         "therefore we are in distributed mode!");
        this->printMsg("#Ranks " + std::to_string(numProcs) + ", this is rank "
                       + std::to_string(rank));
        MPI_Datatype MPI_IT = ttk::getMPIType(static_cast<IT>(0));

        this->printMsg("#Points in Rank " + std::to_string(rank) + ": "
                       + std::to_string(nVerts));
        ttk::Timer fillAndSortTimer;
        std::vector<value<DT, IT>> sortingValues;
        std::vector<IT> gidsToGetVector;
        std::unordered_map<IT, IT> gidToLidMap;
        this->populateVector<DT, IT>(sortingValues, gidsToGetVector,
                                     gidToLidMap, nVerts, scalarArray,
                                     globalIds, rankArray, rank);

        // sort the scalar array distributed first by the scalar value itself,
        // then by the global id
        this->sortVerticesDistributed<DT, IT>(
          sortingValues, this->threadNumber_);
        this->printMsg("#Unique Points in Rank " + std::to_string(rank) + ": "
                       + std::to_string(sortingValues.size()));
        this->printMsg("#Ghostpoints in Rank " + std::to_string(rank) + ": "
                       + std::to_string(gidsToGetVector.size()));
        // when all are done sorting, rank 0 requests the highest values and
        // merges them
        if(rank == 0) {
          this->printMsg(
            "Filling vector and sorting for each rank done, starting merge.", 1,
            fillAndSortTimer.getElapsedTime());
        }

        std::vector<IT> orderedValuesForRank;
        std::vector<value<DT, IT>> finalValues;
        ttk::Timer mergeTimer;
        IT localSize = sortingValues.size();
        IT totalSize;
        // get the complete size  of the dataset by summing up the local sizes
        this->printMsg("Localsize: " + std::to_string(localSize));
        MPI_Reduce(
          &localSize, &totalSize, 1, MPI_IT, MPI_SUM, 0, MPI_COMM_WORLD);
        if(rank == 0) {
          this->printMsg("Total amount of distributed points: "
                         + std::to_string(totalSize));
          this->printMsg("Rank 0 starts merging");
          finalValues.reserve(totalSize);
          IT currentOrder = totalSize - 1;
          std::vector<std::vector<value<DT, IT>>> unsortedReceivedValues;
          unsortedReceivedValues.resize(numProcs);
          std::vector<std::vector<IT>> orderResendValues;
          orderResendValues.resize(numProcs);
          for(int i = 0; i < numProcs; i++) {
            orderResendValues[i].reserve(burstSize);
          }

          // receive the first batch of values
          for(int i = 0; i < numProcs; i++) {
            if(i == 0) {
              std::vector<value<DT, IT>> ownValues;
              this->returnVectorForBurstsize<DT, IT>(
                ownValues, sortingValues, burstSize);
              unsortedReceivedValues[i] = ownValues;
            } else {
              this->ReceiveAndAddToVector<DT, IT>(
                i, structTag, unsortedReceivedValues);
            }
          }
          while(finalValues.size() < totalSize) {
            this->getMax<DT, IT>(intTag, structTag, currentOrder, burstSize,
                                 MPI_IT, finalValues, unsortedReceivedValues,
                                 orderResendValues, orderedValuesForRank,
                                 sortingValues);
          }

          this->printMsg("Finished with sorting, max value is "
                         + std::to_string(finalValues[0].scalar)
                         + ", min value is "
                         + std::to_string(finalValues.back().scalar));
        } else { // other Ranks
          // send the next burstsize values and then wait for an answer from the
          // root rank
          while(sortingValues.size() > 0) {
            std::vector<value<DT, IT>> sendValues;
            this->returnVectorForBurstsize<DT, IT>(
              sendValues, sortingValues, burstSize);
            int size = sendValues.size();
            MPI_Send(sendValues.data(), size * sizeof(value<DT, IT>), MPI_CHAR,
                     0, structTag * rank, MPI_COMM_WORLD);
            std::vector<IT> receivedValues;

            // be prepared to receive burstsize of elements, resize after
            // receiving to the correct size
            receivedValues.resize(burstSize * 2);
            MPI_Status status;
            int amount;
            MPI_Recv(receivedValues.data(), burstSize * 2, MPI_IT, 0,
                     intTag * rank, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_IT, &amount);

            receivedValues.resize(amount);
            orderedValuesForRank.insert(orderedValuesForRank.end(),
                                        receivedValues.begin(),
                                        receivedValues.end());
          }
          // afterwards send once a message of length 0 to root to show that we
          // are done
          MPI_Send(sortingValues.data(), 0, MPI_CHAR, 0, structTag * rank,
                   MPI_COMM_WORLD);
        }

        // all ranks do the following
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0) {
          this->printMsg(
            "Merging done and results sent to ranks, ranks are getting order "
            "for ghost cells and constructing order array.",
            1, mergeTimer.getElapsedTime());
        }
        this->printMsg("#Orders received for Rank " + std::to_string(rank)
                       + ": "
                       + std::to_string(orderedValuesForRank.size() / 2));

        ttk::Timer orderTimer;
        this->buildArrayForReceivedData<IT>(
          orderedValuesForRank.size(), orderedValuesForRank.data(), gidToLidMap,
          orderArray, this->threadNumber_);
        this->printMsg("Built own values, getting Ghost Cell values");

        // we receive the values at the ghostcells through the abstract
        // exchangeGhostCells method
        ttk::exchangeGhostCells<ttk::SimplexId, IT>(orderArray, rankArray,
                                                    globalIds, gidToLidMap,
                                                    nVerts, MPI_COMM_WORLD);
      }
#else
      this->printMsg("MPI not enabled!");
      TTK_FORCE_USE(orderArray);
      TTK_FORCE_USE(scalarArray);
      TTK_FORCE_USE(globalIds);
      TTK_FORCE_USE(rankArray);
      TTK_FORCE_USE(burstSize);
      return 0;
#endif

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

  }; // ArrayPreconditioning class

} // namespace ttk
