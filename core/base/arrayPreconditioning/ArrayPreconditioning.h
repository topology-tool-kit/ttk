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

#ifdef TTK_ENABLE_MPI
#include <mpi.h>
#endif

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
#ifdef TTK_ENABLE_MPI
    MPI_Datatype getMPIType(float val) const {
      TTK_FORCE_USE(val);
      return MPI_FLOAT;
    };
    MPI_Datatype getMPIType(int val) const {
      TTK_FORCE_USE(val);
      return MPI_INT;
    };
    MPI_Datatype getMPIType(unsigned int val) const {
      TTK_FORCE_USE(val);
      return MPI_UNSIGNED;
    };
    MPI_Datatype getMPIType(double val) const {
      TTK_FORCE_USE(val);
      return MPI_DOUBLE;
    };
    MPI_Datatype getMPIType(long double val) const {
      TTK_FORCE_USE(val);
      return MPI_LONG_DOUBLE;
    };
    MPI_Datatype getMPIType(long val) const {
      TTK_FORCE_USE(val);
      return MPI_LONG;
    };
    MPI_Datatype getMPIType(unsigned long val) const {
      TTK_FORCE_USE(val);
      return MPI_UNSIGNED_LONG;
    };
    MPI_Datatype getMPIType(long long val) const {
      TTK_FORCE_USE(val);
      return MPI_LONG_LONG;
    };
    MPI_Datatype getMPIType(unsigned long long val) const {
      TTK_FORCE_USE(val);
      return MPI_UNSIGNED_LONG_LONG;
    };
#endif
    template <typename DT, typename IT>
    struct value {
      DT scalar;
      IT globalId;
      IT localId;
      IT ordering = -1;

      value(DT _scalar, IT _globalId, IT _localId)
        : scalar(_scalar), globalId(_globalId), localId(_localId) {
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
                        const char *const ghostCells) const {
      for(size_t i = 0; i < nVerts; i++) {
        IT globalId = globalIds[i];
        IT localId = i;
        gidToLidMap[globalId] = localId;
        if((int)ghostCells[i] == 0) {
          float scalarValue = scalars[i];
          valuesToSortVector.emplace_back(scalarValue, globalId, localId);
        } else {
          gidsToGetVector.push_back(globalId);
        }
      }
    }

    // orders an value vector first by their scalar value and then by global id
    template <typename DT, typename IT>
    void sortVerticesDistributed(std::vector<value<DT, IT>> &values) const {
      std::sort(
        values.begin(), values.end(), [](value<DT, IT> v1, value<DT, IT> v2) {
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
        outVector.assign(values.begin(), values.end());
        values.clear();
      } else {
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

    // gets all the neighbors of a rank based on the rankarray of the vertices
    // belonging to this rank
    inline void getNeighbors(std::unordered_set<int> &neighbors,
                             const size_t nVertices,
                             const int ownRank,
                             const int *const rankArray) const {
      for(size_t i = 0; i < nVertices; i++) {
        int r = rankArray[i];
        if(r != ownRank)
          neighbors.emplace(r);
      }
    }

    // returns a vector of gid vectors, one vector for each neighbor
    // this shows us where we need to get the order for these gids from
    template <typename IT>
    void getGIdForNeighbors(std::vector<std::vector<IT>> &outVector,
                            const size_t nGIds,
                            const size_t nNeighbors,
                            std::vector<IT> &gidsToGetVector,
                            std::unordered_set<int> &neighbors,
                            std::unordered_map<IT, IT> &gidToLidMap,
                            const int *const rankArray) const {
      outVector.resize(nNeighbors);
      for(size_t i = 0; i < nGIds; i++) {
        IT gId = gidsToGetVector[i];
        IT lId = gidToLidMap[gId];
        int rankForGId = rankArray[lId];
        auto distance
          = std::distance(neighbors.begin(), neighbors.find(rankForGId));
        outVector[distance].push_back(gId);
      }
    }

    template <typename IT>
    void getOrderForGIds(std::vector<IT> &outVector,
                         const size_t nGIds,
                         const IT *const gIds,
                         std::unordered_map<IT, IT> &gidToLidMap,
                         const SimplexId *const order,
                         const int nThreads) const {

      TTK_FORCE_USE(nThreads);
      outVector.resize(nGIds * 2);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nGIds; i++) {
        IT gId = gIds[i];
        IT lId = gidToLidMap[gId];
        SimplexId orderForThisGId = order[lId];
        outVector[2 * i] = gId;
        outVector[2 * i + 1] = orderForThisGId;
      }
    }
#ifdef TTK_ENABLE_MPI
    template <typename DT, typename IT>
    void ReceiveAndAddToVector(
      MPI_Datatype mpi_values,
      int rankFrom,
      int structTag,
      int intTag,
      std::vector<std::vector<value<DT, IT>>> &unsortedReceivedValues) const {
      std::vector<value<DT, IT>> receivedValues;
      // be prepared to receive burstsize of elements, resize after receiving to
      // the correct size
      int amount;
      MPI_Recv(&amount, 1, MPI_INT, rankFrom, intTag, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      receivedValues.resize(amount, {0, 0, 0});
      MPI_Recv(receivedValues.data(), amount, mpi_values, rankFrom, structTag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      unsortedReceivedValues[rankFrom] = receivedValues;
    }
#endif
    template <typename DT, typename IT>
    int processScalarArray(ttk::SimplexId *orderArray,
                           const DT *scalarArray,
                           const IT *globalIds,
                           const int *rankArray,
                           const char *ghostCells,
                           const size_t nVertices,
                           const int burstSize) const { // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(nVertices)},
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
        int intTag = 100;
        int structTag = 101;
        int boolTag = 102;
        if(rank == 0)
          this->printMsg("Global Point Ids, Ghost Cells and RankArray exist, "
                         "therefore we are in distributed mode!");
        this->printMsg("#Ranks " + std::to_string(numProcs) + ", this is rank "
                       + std::to_string(rank));
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Datatype mpi_values;
        const int nitems = 4;
        int blocklengths[4] = {1, 1, 1, 1};
        MPI_Datatype MPI_DT = this->getMPIType(static_cast<DT>(0));
        MPI_Datatype MPI_IT = this->getMPIType(static_cast<IT>(0));
        MPI_Datatype types[4] = {MPI_DT, MPI_IT, MPI_IT, MPI_IT};
        MPI_Aint offsets[4];
        typedef value<DT, IT> value_DT_IT;
        offsets[0] = offsetof(value_DT_IT, scalar);
        offsets[1] = offsetof(value_DT_IT, globalId);
        offsets[2] = offsetof(value_DT_IT, localId);
        offsets[3] = offsetof(value_DT_IT, ordering);
        MPI_Type_create_struct(
          nitems, blocklengths, offsets, types, &mpi_values);
        MPI_Type_commit(&mpi_values);

        this->printMsg("#Points in Rank " + std::to_string(rank) + ": "
                       + std::to_string(nVertices));
        ttk::Timer fillAndSortTimer;
        std::vector<value_DT_IT> sortingValues;
        std::vector<IT> gidsToGetVector;
        std::unordered_map<IT, IT> gidToLidMap;
        this->populateVector<DT, IT>(sortingValues, gidsToGetVector,
                                     gidToLidMap, nVertices, scalarArray,
                                     globalIds, ghostCells);

        // sort the scalar array distributed first by the scalar value itself,
        // then by the global id
        this->sortVerticesDistributed<DT, IT>(sortingValues);
        this->printMsg("#Unique Points in Rank " + std::to_string(rank) + ": "
                       + std::to_string(sortingValues.size()));
        this->printMsg("#Ghostpoints in Rank " + std::to_string(rank) + ": "
                       + std::to_string(gidsToGetVector.size()));
        // when all are done sorting, rank 0 requests the highest values and
        // merges them
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0) {
          this->printMsg(
            "Filling vector and sorting for each rank done, starting merge.", 1,
            fillAndSortTimer.getElapsedTime());
        }

        std::vector<IT> orderedValuesForRank;
        std::vector<value_DT_IT> finalValues;
        ttk::Timer mergeTimer;
        size_t totalSize = 0;
        if(rank == 0) {
          this->printMsg("Rank 0 starts merging");
          // get the nVertices from each rank, add them to get the complete size
          // of the dataset
          for(int i = 0; i < numProcs; i++) {
            if(i == 0) {
              totalSize += sortingValues.size();
            } else {
              IT receivedSize;
              MPI_Recv(&receivedSize, 1, MPI_IT, i, intTag, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
              totalSize += receivedSize;
            }
          }
          this->printMsg("Total amount of distributed points: "
                         + std::to_string(totalSize));
          IT currentOrder = totalSize - 1;
          std::vector<std::vector<value_DT_IT>> unsortedReceivedValues;
          unsortedReceivedValues.resize(numProcs);
          std::vector<std::vector<IT>> orderResendValues;
          orderResendValues.resize(numProcs);

          // receive the first batch of values
          for(int i = 0; i < numProcs; i++) {
            if(i == 0) {
              std::vector<value_DT_IT> ownValues;
              this->returnVectorForBurstsize<DT, IT>(
                ownValues, sortingValues, burstSize);
              unsortedReceivedValues[i] = ownValues;
            } else {
              this->ReceiveAndAddToVector<DT, IT>(
                mpi_values, i, structTag, intTag, unsortedReceivedValues);
            }
          }

          while(finalValues.size() < totalSize) {
            // take the current maximum scalar over all ranks
            int rankIdOfMaxScalar = -1;
            float maxScalar = -FLT_MAX;
            IT maxGId = -1;
            for(int i = 0; i < numProcs; i++) {
              if(unsortedReceivedValues[i].size() > 0) {
                int thisId = i;
                float thisScalar = unsortedReceivedValues[i].back().scalar;
                IT thisGId = unsortedReceivedValues[i].back().globalId;
                if(thisScalar > maxScalar
                   || (thisScalar == maxScalar && thisGId > maxGId)) {
                  maxScalar = thisScalar;
                  maxGId = thisGId;
                  rankIdOfMaxScalar = thisId;
                }
              }
            }
            if(rankIdOfMaxScalar == -1) {
              this->printMsg("FinalValues.size: "
                             + std::to_string(finalValues.size()));
              this->printMsg("All vectors are empty, but out final vector is "
                             "not complete yet. Either something went wrong or "
                             "some rank didn't send their values yet.");
              return 0;
            }

            // move the struct from the unsortedReceivedValues subvector to the
            // finalValues vector to get an ordering
            value_DT_IT currentValue
              = unsortedReceivedValues[rankIdOfMaxScalar].back();
            // we send the globalId and the the order via one send command,
            // therefore we need to check two concurrent values at the same time
            // later on
            orderResendValues[rankIdOfMaxScalar].push_back(
              currentValue.globalId);
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
                  std::vector<value_DT_IT> ownValues;
                  this->returnVectorForBurstsize<DT, IT>(
                    ownValues, sortingValues, burstSize);
                  unsortedReceivedValues[rankIdOfMaxScalar] = ownValues;
                } else {
                  this->printMsg("We are done with rank 0!");
                }
              } else {
                // receive more values from rank, send ordering to the rank
                // send to the finished rank that we want more
                MPI_Send(orderResendValues[rankIdOfMaxScalar].data(),
                         orderResendValues[rankIdOfMaxScalar].size(), MPI_IT,
                         rankIdOfMaxScalar, intTag, MPI_COMM_WORLD);
                orderResendValues[rankIdOfMaxScalar].clear();
                // check if there are more values to be received. If so, receive
                // them
                bool moreVals = false;
                MPI_Recv(&moreVals, 1, MPI_CXX_BOOL, rankIdOfMaxScalar, boolTag,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(moreVals) {
                  this->ReceiveAndAddToVector(mpi_values, rankIdOfMaxScalar,
                                              structTag, intTag,
                                              unsortedReceivedValues);
                } else {
                  this->printMsg("We are done with rank "
                                 + std::to_string(rankIdOfMaxScalar));
                }
              }
            }
          }

          this->printMsg("Finished with sorting, max value is "
                         + std::to_string(finalValues[0].scalar)
                         + ", min value is "
                         + std::to_string(finalValues.back().scalar));
        } else {
          IT nValues = sortingValues.size();
          MPI_Send(&nValues, 1, MPI_IT, 0, intTag, MPI_COMM_WORLD);

          // send the next burstsize values and then wait for an answer from the
          // root rank
          while(sortingValues.size() > 0) {
            std::vector<value_DT_IT> sendValues;
            this->returnVectorForBurstsize<DT, IT>(
              sendValues, sortingValues, burstSize);
            int size = sendValues.size();
            MPI_Send(&size, 1, MPI_INT, 0, intTag, MPI_COMM_WORLD);
            MPI_Send(sendValues.data(), size, mpi_values, 0, structTag,
                     MPI_COMM_WORLD);
            std::vector<IT> receivedValues;

            // be prepared to receive burstsize of elements, resize after
            // receiving to the correct size
            receivedValues.resize(burstSize * 2);
            MPI_Status status;
            int amount;
            MPI_Recv(receivedValues.data(), burstSize * 2, MPI_IT, 0, intTag,
                     MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_IT, &amount);
            receivedValues.resize(amount);
            orderedValuesForRank.insert(orderedValuesForRank.end(),
                                        receivedValues.begin(),
                                        receivedValues.end());

            // afterwards send to root if there are still values to be sent
            bool moreVals = sortingValues.size() > 0;
            MPI_Send(&moreVals, 1, MPI_CXX_BOOL, 0, boolTag, MPI_COMM_WORLD);
          }
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
                       + ": " + std::to_string(orderedValuesForRank.size()));

        ttk::Timer orderTimer;

        this->buildArrayForReceivedData<IT>(
          orderedValuesForRank.size(), orderedValuesForRank.data(), gidToLidMap,
          orderArray, this->threadNumber_);

        // we still need to get the data for the ghostcells from their ranks
        std::unordered_set<int> neighbors;
        this->getNeighbors(neighbors, nVertices, rank, rankArray);
        this->printMsg("Rank " + std::to_string(rank)
                       + " #Neighbors: " + std::to_string(neighbors.size()));
        std::vector<std::vector<IT>> gIdForNeighbors;
        this->getGIdForNeighbors<IT>(gIdForNeighbors, gidsToGetVector.size(),
                                     neighbors.size(), gidsToGetVector,
                                     neighbors, gidToLidMap, rankArray);
        MPI_Request req;
        std::vector<size_t> sizesToSend;
        sizesToSend.resize(gIdForNeighbors.size());
        for(size_t i = 0; i < gIdForNeighbors.size(); i++) {
          auto iter = neighbors.begin();
          std::advance(iter, i);
          int rankToSend = *iter;
          std::vector<IT> gIdToSend = gIdForNeighbors[i];
          sizesToSend[i] = gIdToSend.size();
          // we need to send to all neighbors how many we will send them,
          // and then send the values to neighbors where the amount > 0
          // first we Isend everything we need to the neighbors
          MPI_Isend(&sizesToSend[i], 1, MPI_IT, rankToSend, intTag,
                    MPI_COMM_WORLD, &req);
          MPI_Request_free(&req);
          if(sizesToSend[i] > 0) {
            MPI_Isend(gIdForNeighbors[i].data(), sizesToSend[i], MPI_IT,
                      rankToSend, intTag, MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
          }
        }

        // then we blockingly receive from our neighbors
        std::vector<std::vector<IT>> ordersToSend;
        ordersToSend.resize(neighbors.size());
        for(size_t i = 0; i < neighbors.size(); i++) {
          auto iter = neighbors.begin();
          std::advance(iter, i);
          int rankToRecv = *iter;
          size_t sizeToRecv;
          MPI_Recv(&sizeToRecv, 1, MPI_IT, rankToRecv, intTag, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);

          if(sizeToRecv > 0) {
            std::vector<IT> gIdToRecv;
            gIdToRecv.resize(sizeToRecv);
            MPI_Recv(gIdToRecv.data(), sizeToRecv, MPI_IT, rankToRecv, intTag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // prepare what they need
            this->getOrderForGIds<IT>(ordersToSend[i], gIdToRecv.size(),
                                      gIdToRecv.data(), gidToLidMap, orderArray,
                                      this->threadNumber_);

            // Isend it to them
            MPI_Isend(ordersToSend[i].data(), ordersToSend[i].size(), MPI_IT,
                      rankToRecv, intTag, MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
          }
        }

        // and blockingly receive the order from our neighbors
        for(size_t i = 0; i < gIdForNeighbors.size(); i++) {
          auto iter = neighbors.begin();
          std::advance(iter, i);
          int rankToRecvOrder = *iter;
          auto sizeToRecvOrder = gIdForNeighbors[i].size() * 2;
          std::vector<IT> ordersToRecv;
          ordersToRecv.resize(sizeToRecvOrder);
          if(sizeToRecvOrder > 0) {
            MPI_Recv(ordersToRecv.data(), sizeToRecvOrder, MPI_IT,
                     rankToRecvOrder, intTag, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            // and add it to our orderarray
            this->buildArrayForReceivedData<IT>(
              ordersToRecv.size(), ordersToRecv.data(), gidToLidMap, orderArray,
              this->threadNumber_);
          }
        }

        // at the end, free up the MPI Datatype
        MPI_Type_free(&mpi_values);
      }
#else
      this->printMsg("MPI not enabled!");
      TTK_FORCE_USE(orderArray);
      TTK_FORCE_USE(scalarArray);
      TTK_FORCE_USE(globalIds);
      TTK_FORCE_USE(rankArray);
      TTK_FORCE_USE(ghostCells);
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
