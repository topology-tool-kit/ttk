/// \ingroup base
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date April 2022
///
/// \brief Utilities for MPI implementation.

#pragma once

#include <BaseClass.h>
#include <Timer.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#if TTK_ENABLE_MPI
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

namespace ttk {

  inline MPI_Datatype getMPIType(const float ttkNotUsed(val)) {
    return MPI_FLOAT;
  };
  inline MPI_Datatype getMPIType(const int ttkNotUsed(val)) {
    return MPI_INT;
  };
  inline MPI_Datatype getMPIType(const unsigned int ttkNotUsed(val)) {
    return MPI_UNSIGNED;
  };
  inline MPI_Datatype getMPIType(const double ttkNotUsed(val)) {
    return MPI_DOUBLE;
  };
  inline MPI_Datatype getMPIType(const long double ttkNotUsed(val)) {
    return MPI_LONG_DOUBLE;
  };
  inline MPI_Datatype getMPIType(const long ttkNotUsed(val)) {
    return MPI_LONG;
  };
  inline MPI_Datatype getMPIType(const unsigned long ttkNotUsed(val)) {
    return MPI_UNSIGNED_LONG;
  };
  inline MPI_Datatype getMPIType(const long long ttkNotUsed(val)) {
    return MPI_LONG_LONG;
  };
  inline MPI_Datatype getMPIType(const unsigned long long ttkNotUsed(val)) {
    return MPI_UNSIGNED_LONG_LONG;
  };

  inline bool isRunningWithMPI() {
    return ttk::MPIsize_ > 1;
  };

  inline int startMPITimer(Timer &t, int rank, int size) {
    if(size > 0) {
      MPI_Barrier(MPI_COMM_WORLD);
      if(rank == 0) {
        t.reStart();
      }
    }
    return 0;
  };

  inline double endMPITimer(Timer &t, int rank, int size) {
    double elapsedTime = 0;
    if(size > 0) {
      MPI_Barrier(MPI_COMM_WORLD);
      if(rank == 0) {
        elapsedTime = t.getElapsedTime();
      }
    }

    return elapsedTime;
  };

  /**
   * @brief Gather vectors on a specific rank
   *
   * @param[out] dst Container storing vectors copied from other ranks
   * @param[in] src Data to send to \ref destRank and to store in \ref dst
   * @param[in] destRank Destination process identifier
   * @return 0 in case of success
   */
  template <typename T>
  int gatherVectors(std::vector<std::vector<T>> &dst,
                    const std::vector<T> &src,
                    const int destRank) {

    if(!ttk::isRunningWithMPI()) {
      return -1;
    }

    // src sizes (gathered on rank 0)
    std::vector<unsigned long> vecSizes{};

    if(ttk::MPIrank_ == destRank) {
      vecSizes.resize(MPIsize_);
      dst.resize(MPIsize_);
    }

    const unsigned long localSize = src.size();
    // gather src sizes on destRank
    MPI_Gather(&localSize, 1, MPI_UNSIGNED_LONG, vecSizes.data(), 1,
               MPI_UNSIGNED_LONG, destRank, MPI_COMM_WORLD);

    if(ttk::MPIrank_ == destRank) {
      // allocate dst with vecSizes
      for(int i = 0; i < ttk::MPIsize_; ++i) {
        if(i == destRank) {
          continue;
        }
        dst[i].resize(vecSizes[i]);
      }

      for(int i = 0; i < ttk::MPIsize_; ++i) {
        if(i == destRank) {
          continue;
        }
        // receive src content from other ranks
        MPI_Recv(dst[i].data(), dst[i].size(), ttk::getMPIType(src[0]), i,
                 MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      dst[destRank] = std::move(src);

    } else {
      // send src content to destRank
      MPI_Send(src.data(), src.size(), ttk::getMPIType(src[0]), destRank, 0,
               MPI_COMM_WORLD);
    }

    return 0;
  }

  /**
   * @brief Request all ghost cell scalar data from one rank from their owning
   * ranks and get the data
   *
   * @param[out] scalarArray the scalar array which we want to fill and which is
   * filled on the other ranks
   * @param[in] rankArray the owner array for the scalar data
   * @param[in] globalIds the global id array for the scalar data
   * @param[in] gidToLidMap a map which translates global ids to local,
   * rank-based ids
   * @param[in] neighbors a set of rank neighbors of rankToSend, to prevent
   * unnecessary communication
   * @param[in] rankToSend Destination process identifier
   * @param[in] nVerts number of vertices in the arrays
   * @param[in] communicator the communicator over which the ranks are connected
   * (most likely MPI_COMM_WORLD)
   * @return 0 in case of success
   */
  template <typename DT, typename IT>
  int getGhostCellScalars(DT *scalarArray,
                          const int *const rankArray,
                          const IT *const globalIds,
                          const std::unordered_map<IT, IT> &gidToLidMap,
                          const std::unordered_set<int> &neighbors,
                          const int rankToSend,
                          const IT nVerts,
                          MPI_Comm communicator) {
    if(!ttk::isRunningWithMPI()) {
      return -1;
    }
    MPI_Datatype MPI_DT = getMPIType(static_cast<DT>(0));
    MPI_Datatype MPI_IT = getMPIType(static_cast<IT>(0));
    // we need unique tags for each rankToSend, otherwise messages might become
    // entangled
    int tagMultiplier = rankToSend + 1;
    int amountTag = 101 * tagMultiplier;
    int idsTag = 102 * tagMultiplier;
    int valuesTag = 103 * tagMultiplier;
    if(rankToSend == ttk::MPIrank_) {
      // initialize the inner vectors with size 0
      std::vector<std::vector<IT>> rankVectors(
        ttk::MPIsize_, std::vector<IT>(0));
      // aggregate the needed ids

      for(IT i = 0; i < nVerts; i++) {
        if(ttk::MPIrank_ != rankArray[i]) {
          rankVectors[rankArray[i]].push_back(globalIds[i]);
        }
      }
      // send the amount of ids and the needed ids themselves
      for(int r = 0; r < ttk::MPIsize_; r++) {
        if(ttk::MPIrank_ != r && neighbors.find(r) != neighbors.end()) {
          IT nValues = rankVectors[r].size();
          MPI_Send(&nValues, 1, MPI_IT, r, amountTag, communicator);
          if(nValues > 0) {
            MPI_Send(
              rankVectors[r].data(), nValues, MPI_IT, r, idsTag, communicator);
          }
        }
      }
      // receive the scalar values
      for(int r = 0; r < ttk::MPIsize_; r++) {
        if(ttk::MPIrank_ != r && neighbors.find(r) != neighbors.end()) {
          IT nValues = rankVectors[r].size();
          std::vector<DT> receivedValues(nValues);
          if(nValues > 0) {
            MPI_Recv(receivedValues.data(), nValues, MPI_DT, r, valuesTag,
                     communicator, MPI_STATUS_IGNORE);
            for(IT i = 0; i < nValues; i++) {
              DT receivedVal = receivedValues[i];
              IT globalId = rankVectors[r][i];
              IT localId = gidToLidMap.at(globalId);
              scalarArray[localId] = receivedVal;
            }
          }
        }
      }

    } else { // owner ranks
      // if rankToSend is not the neighbor of the current rank, we do not need
      // to do anything
      if(neighbors.find(rankToSend) != neighbors.end()) {
        // receive the amount of ids and the needed ids themselves
        IT nValues;

        MPI_Recv(&nValues, 1, MPI_IT, rankToSend, amountTag, communicator,
                 MPI_STATUS_IGNORE);

        if(nValues > 0) {
          std::vector<IT> receivedIds(nValues);
          MPI_Recv(receivedIds.data(), nValues, MPI_IT, rankToSend, idsTag,
                   communicator, MPI_STATUS_IGNORE);

          // assemble the scalar values
          std::vector<DT> valuesToSend(nValues);
          for(IT i = 0; i < nValues; i++) {
            IT globalId = receivedIds[i];
            IT localId = gidToLidMap.at(globalId);
            valuesToSend[i] = scalarArray[localId];
          }

          // send the scalar values
          MPI_Send(valuesToSend.data(), nValues, MPI_DT, rankToSend, valuesTag,
                   communicator);
        }
      }
    }

    return 0;
  }

  /**
   * @brief get the neighbors of a rank by traversing the rankArray
   *
   * @param[out] neighbors a set containing the ranks which are neighbors of
   * this rank
   * @param[in] rankArray the owner array for the scalar data
   * @param[in] nVerts the number of vertices in rankArray
   * @return 0 in case of success
   */
  template <typename IT>
  int getNeighbors(std::unordered_set<int> &neighbors,
                   const int *const rankArray,
                   const IT nVerts) {
    for(IT i = 0; i < nVerts; i++) {
      if(rankArray[i] != ttk::MPIrank_) {
        neighbors.emplace(rankArray[i]);
      }
    }
    return 0;
  }

  /**
   * @brief exchange all ghost cell information by calling getGhostCellScalars
   * for every rank
   *
   * @param[out] scalarArray the scalar array which we want to fill and which is
   * filled on the other ranks
   * @param[in] rankArray the owner array for the scalar data
   * @param[in] globalIds the global id array for the scalar data
   * @param[in] gidToLidMap a map which translates global ids to local,
   * rank-based ids
   * @param[in] nVerts number of vertices in the arrays
   * @param[in] communicator the communicator over which the ranks are connected
   * (most likely MPI_COMM_WORLD)
   * @return 0 in case of success
   */
  template <typename DT, typename IT>
  int exchangeGhostCells(DT *scalarArray,
                         const int *const rankArray,
                         const IT *const globalIds,
                         const std::unordered_map<IT, IT> gidToLidMap,
                         const IT nVerts,
                         MPI_Comm communicator) {
    if(!ttk::isRunningWithMPI()) {
      return -1;
    }
    std::unordered_set<int> neighbors;
    getNeighbors<IT>(neighbors, rankArray, nVerts);
    for(int r = 0; r < ttk::MPIsize_; r++) {
      getGhostCellScalars<DT, IT>(scalarArray, rankArray, globalIds,
                                  gidToLidMap, neighbors, r, nVerts,
                                  communicator);
      MPI_Barrier(communicator);
    }
    return 0;
  }

  void inline produceRankArray(std::vector<int> &rankArray,
                               long int *globalIds,
                               unsigned char *ghostCells,
                               int nVertices) {
    MPI_Datatype MIT = ttk::getMPIType(static_cast<ttk::SimplexId>(0));
    MPI_Comm ttkGhostCellPreconditioningComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &ttkGhostCellPreconditioningComm);
    std::vector<ttk::SimplexId> currentRankUnknownIds;
    std::vector<std::vector<ttk::SimplexId>> allUnknownIds(ttk::MPIsize_);
    std::unordered_set<ttk::SimplexId> gIdSet;
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdToLocalMap;

    for(int i = 0; i < nVertices; i++) {
      int ghostCellVal = ghostCells[i];
      ttk::SimplexId globalId = globalIds[i];
      if(ghostCellVal == 0) {
        // if the ghost cell value is 0, then this vertex mainly belongs to
        // this rank
        rankArray[i] = ttk::MPIrank_;
        gIdSet.insert(globalId);
      } else {
        // otherwise the vertex belongs to another rank and we need to find
        // out to which one this needs to be done by broadcasting the global
        // id and hoping for some other rank to answer
        currentRankUnknownIds.push_back(globalId);
        gIdToLocalMap[globalId] = i;
      }
    }

    allUnknownIds[ttk::MPIrank_] = currentRankUnknownIds;
    ttk::SimplexId sizeOfCurrentRank;
    // first each rank gets the information which rank needs which globalid
    for(int r = 0; r < ttk::MPIsize_; r++) {
      if(r == ttk::MPIrank_)
        sizeOfCurrentRank = currentRankUnknownIds.size();
      MPI_Bcast(&sizeOfCurrentRank, 1, MIT, r, ttkGhostCellPreconditioningComm);
      allUnknownIds[r].resize(sizeOfCurrentRank);
      MPI_Bcast(allUnknownIds[r].data(), sizeOfCurrentRank, MIT, r,
                ttkGhostCellPreconditioningComm);
    }
    // then we check if the needed globalid values are present in the local
    // globalid map if so, we send the rank value to the requesting rank
    std::vector<ttk::SimplexId> gIdsToSend;
    for(int r = 0; r < ttk::MPIsize_; r++) {
      if(r != ttk::MPIrank_) {
        // send the needed values to r
        gIdsToSend.clear();
        for(ttk::SimplexId gId : allUnknownIds[r]) {
          if(gIdSet.count(gId)) {
            // add the value to the vector which will be sent
            gIdsToSend.push_back(gId);
          }
        }
        // send whole vector of data
        MPI_Send(gIdsToSend.data(), gIdsToSend.size(), MIT, r, 101,
                 ttkGhostCellPreconditioningComm);
      } else {
        // receive a variable amount of values from different ranks
        size_t i = 0;
        std::vector<ttk::SimplexId> receivedGlobals;
        while(i < allUnknownIds[ttk::MPIrank_].size()) {
          receivedGlobals.resize(allUnknownIds[ttk::MPIrank_].size());
          MPI_Status status;
          int amount;
          MPI_Recv(receivedGlobals.data(), allUnknownIds[ttk::MPIrank_].size(),
                   MIT, MPI_ANY_SOURCE, MPI_ANY_TAG,
                   ttkGhostCellPreconditioningComm, &status);
          int sourceRank = status.MPI_SOURCE;
          MPI_Get_count(&status, MIT, &amount);
          receivedGlobals.resize(amount);
          for(ttk::SimplexId receivedGlobal : receivedGlobals) {
            ttk::SimplexId localVal = gIdToLocalMap[receivedGlobal];
            rankArray[localVal] = sourceRank;
            i++;
          }
        }
      }
    }
    // free the communicator once we are done with everything MPI
    MPI_Comm_free(&ttkGhostCellPreconditioningComm);
  }
} // namespace ttk
#endif
