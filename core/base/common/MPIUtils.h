/// \ingroup base
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date April 2022
///
/// \brief Utilities for MPI implementation.

#pragma once

#include <BaseClass.h>
#include <Timer.h>
#include <array>
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
      MPI_Barrier(ttk::MPIcomm_);
      if(rank == 0) {
        t.reStart();
      }
    }
    return 0;
  };

  inline double endMPITimer(Timer &t, int rank, int size) {
    double elapsedTime = 0;
    if(size > 0) {
      MPI_Barrier(ttk::MPIcomm_);
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
               MPI_UNSIGNED_LONG, destRank, ttk::MPIcomm_);

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
                 MPI_ANY_TAG, ttk::MPIcomm_, MPI_STATUS_IGNORE);
      }
      dst[destRank] = std::move(src);

    } else {
      // send src content to destRank
      MPI_Send(src.data(), src.size(), ttk::getMPIType(src[0]), destRank, 0,
               ttk::MPIcomm_);
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
   * (most likely ttk::MPIcomm_)
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
                   const IT nVerts,
                   MPI_Comm communicator) {
    for(IT i = 0; i < nVerts; i++) {
      if(rankArray[i] != ttk::MPIrank_) {
        neighbors.emplace(rankArray[i]);
      }
    }
    std::vector<int> sendVector(neighbors.begin(), neighbors.end());
    int localSize = neighbors.size();
    int sizes[ttk::MPIsize_];
    int displacements[ttk::MPIsize_];
    MPI_Gather(&localSize, 1, MPI_INT, sizes, 1, MPI_INT, 0, communicator);
    int totalSize = 0;
    if(ttk::MPIrank_ == 0) {
      for(int i = 0; i < ttk::MPIsize_; i++) {
        totalSize += sizes[i];
        if(i == 0) {
          displacements[i] = 0;
        } else {
          displacements[i] = displacements[i - 1] + sizes[i - 1];
        }
      }
    }
    std::vector<int> rootVector(totalSize);

    MPI_Gatherv(sendVector.data(), sendVector.size(), MPI_INT,
                rootVector.data(), sizes, displacements, MPI_INT, 0,
                communicator);
    std::vector<int> scatterVector;

    if(ttk::MPIrank_ == 0) {
      // now we transform this 1d vector in a correct vector of sets which we
      // will transform back to scatter it
      std::vector<std::unordered_set<int>> setsFromRanks(ttk::MPIsize_);
      auto begin = rootVector.begin();
      auto end = rootVector.begin();
      for(int i = 0; i < ttk::MPIsize_; i++) {
        end = begin + sizes[i];
        std::unordered_set<int> s(begin, end);
        setsFromRanks[i] = s;
        begin = end;
        // std::cout << "R" << std::to_string(i) << " nr needs something from "
        // << std::to_string(setsFromRanks[i].size()) << " neighbors." <<
        // std::endl;
      }
      // now we need to check for each rank if they are a neighbor of any other
      // rank. If so, we need to add those to the neighbors
      for(int i = 0; i < ttk::MPIsize_; i++) {
        for(int j = 0; i < ttk::MPIsize_; i++) {
          if(setsFromRanks[j].find(i) != setsFromRanks[j].end()) {
            setsFromRanks[i].emplace(j);
          }
        }
      }
      // now we transform this vector of sets back into a 1d vector to scatter
      // it
      for(int i = 0; i < ttk::MPIsize_; i++) {
        sizes[i] = setsFromRanks[i].size();
        if(i == 0) {
          displacements[i] = 0;
        } else {
          displacements[i] = displacements[i - 1] + sizes[i - 1];
        }
        scatterVector.insert(scatterVector.end(), setsFromRanks[i].begin(),
                             setsFromRanks[i].end());
      }
    }

    // scatter first the size and then the vector itself
    int receivedSize;
    MPI_Scatter(sizes, 1, MPI_INT, &receivedSize, 1, MPI_INT, 0, communicator);

    // and then the actual neighbors
    std::vector<int> receivedNeighbors(receivedSize);
    MPI_Scatterv(scatterVector.data(), sizes, displacements, MPI_INT,
                 receivedNeighbors.data(), receivedSize, MPI_INT, 0,
                 communicator);
    // then we turn the vector back into a set
    std::unordered_set<int> finalSet(
      receivedNeighbors.begin(), receivedNeighbors.end());
    neighbors = finalSet;

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
   * (most likely ttk::MPIcomm_)
   * @return 0 in case of success
   */
  template <typename DT, typename IT>
  int exchangeGhostCells(DT *scalarArray,
                         const int *const rankArray,
                         const IT *const globalIds,
                         const std::unordered_map<IT, IT> gidToLidMap,
                         const IT nVerts,
                         MPI_Comm communicator = MPI_COMM_WORLD) {
    if(!ttk::isRunningWithMPI()) {
      return -1;
    }
    std::unordered_set<int> neighbors;
    getNeighbors<IT>(neighbors, rankArray, nVerts, communicator);
    for(int r = 0; r < ttk::MPIsize_; r++) {
      getGhostCellScalars<DT, IT>(scalarArray, rankArray, globalIds,
                                  gidToLidMap, neighbors, r, nVerts,
                                  communicator);
      MPI_Barrier(communicator);
    }
    return 0;
  }

  // returns true if bounding boxes intersect, false if not
  bool inline checkForIntersection(double *myBB, double *theirBB) {
    return !(
      myBB[0] > theirBB[1] // my left side is right of their right side
      || myBB[1] < theirBB[0] // my right side is left of their left side
      || myBB[2] > theirBB[3] // my bottom side is above their top side
      || myBB[3] < theirBB[2] // my top side is under their bottom side
      || myBB[4] > theirBB[5] // my front side is behind their back side
      || myBB[5] < theirBB[4] // my back side is in front of their front side
    );
  }

  /**
   * @brief produce the RankArray array, that stores rank ownership information
   *
   * @param[out] rankArray the owner array for the scalar data
   * @param[in] globalIds the global id array for the scalar data
   * @param[in] ghostCells the ghost array for the scalar data
   * @param[in] nVertices number of vertices in the arrays
   */

  void inline produceRankArray(std::vector<int> &rankArray,
                               long int *globalIds,
                               unsigned char *ghostCells,
                               int nVertices,
                               double *boundingBox) {
    std::vector<std::array<double, 6>> rankBoundingBoxes(
      ttk::MPIsize_, std::array<double, 6>({}));
    std::copy(
      boundingBox, boundingBox + 6, rankBoundingBoxes[ttk::MPIrank_].begin());
    for(int r = 0; r < ttk::MPIsize_; r++) {
      MPI_Bcast(rankBoundingBoxes[r].data(), 6, MPI_DOUBLE, r, ttk::MPIcomm_);
    }

    double epsilon = 0.00001;
    // inflate our own bounding box by epsilon
    for(int i = 0; i < 6; i++) {
      if(i % 2 == 0)
        boundingBox[i] -= epsilon;
      if(i % 2 == 1)
        boundingBox[i] += epsilon;
    }
    std::vector<int> neighbors;
    for(int i = 0; i < ttk::MPIsize_; i++) {
      if(i != ttk::MPIrank_) {
        double *theirBoundingBox = rankBoundingBoxes[i].data();
        if(checkForIntersection(boundingBox, theirBoundingBox)) {
          neighbors.push_back(i);
        }
      }
    }
    MPI_Datatype MIT = ttk::getMPIType(static_cast<ttk::SimplexId>(0));
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

    ttk::SimplexId sizeOfCurrentRank = currentRankUnknownIds.size();
    std::vector<ttk::SimplexId> gIdsToSend;
    std::vector<ttk::SimplexId> receivedGlobals;
    receivedGlobals.resize(sizeOfCurrentRank);
    ttk::SimplexId sizeOfNeighbor;
    std::vector<ttk::SimplexId> neighborUnknownIds;
    for(int neighbor : neighbors) {
      // we first send the size and then all needed ids to the neighbor
      MPI_Sendrecv(&sizeOfCurrentRank, 1, MIT, neighbor, ttk::MPIrank_,
                   &sizeOfNeighbor, 1, MIT, neighbor, neighbor, ttk::MPIcomm_,
                   MPI_STATUS_IGNORE);
      neighborUnknownIds.resize(sizeOfNeighbor);
      gIdsToSend.reserve(sizeOfNeighbor);

      MPI_Sendrecv(currentRankUnknownIds.data(), sizeOfCurrentRank, MIT,
                   neighbor, ttk::MPIrank_, neighborUnknownIds.data(),
                   sizeOfNeighbor, MIT, neighbor, neighbor, ttk::MPIcomm_,
                   MPI_STATUS_IGNORE);

      // then we check if the needed globalid values are present in the local
      // globalid set if so, we send the rank value to the requesting rank
      for(ttk::SimplexId gId : neighborUnknownIds) {
        if(gIdSet.count(gId)) {
          // add the value to the vector which will be sent
          gIdsToSend.push_back(gId);
        }
      }
      MPI_Status status;
      int amount;

      MPI_Sendrecv(gIdsToSend.data(), gIdsToSend.size(), MIT, neighbor,
                   ttk::MPIrank_, receivedGlobals.data(),
                   currentRankUnknownIds.size(), MIT, neighbor, neighbor,
                   ttk::MPIcomm_, &status);

      MPI_Get_count(&status, MIT, &amount);
      receivedGlobals.resize(amount);

      for(ttk::SimplexId receivedGlobal : receivedGlobals) {
        ttk::SimplexId localVal = gIdToLocalMap[receivedGlobal];
        rankArray[localVal] = neighbor;
      }
      // cleanup
      gIdsToSend.clear();
      receivedGlobals.resize(sizeOfCurrentRank);
    }
  }
} // namespace ttk
#endif
