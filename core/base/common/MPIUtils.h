/// \ingroup base
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date April 2022
///
/// \brief Utilities for MPI implementation.

#pragma once

#include <BaseClass.h>
#include <Timer.h>
#include <iostream>

#include <algorithm>
#include <array>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef TTK_ENABLE_MPI

// disable the MPI C++ API
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
  inline MPI_Datatype getMPIType(const unsigned char ttkNotUsed(val)) {
    return MPI_UNSIGNED_CHAR;
  };

  template <typename DT, typename IT>
  struct value {
    DT scalar;
    IT globalId;

    value(DT _scalar, IT _globalId) : scalar(_scalar), globalId(_globalId) {
    }
  };

  inline bool isRunningWithMPI() {
    return ttk::MPIsize_ > 1;
  };

  inline bool hasInitializedMPI() {
    return ttk::MPIsize_ > 0;
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
   * @brief Request all ghost cell scalar data from one rank from their owning
   * ranks and get the data
   *
   * @param[out] scalarArray the scalar array which we want to fill and which is
   * filled on the other ranks
   * @param[in] triangulation the triangulation for the data
   * @param[in] rankToSend Destination process identifier
   * @param[in] communicator the communicator over which the ranks are connected
   * (most likely ttk::MPIcomm_)
   * @return 0 in case of success
   */
  template <typename DT, typename triangulationType>
  int getGhostCellScalars(DT *scalarArray,
                          const triangulationType *triangulation,
                          const int rankToSend,
                          MPI_Comm communicator,
                          const int dimensionNumber) {
    const std::vector<int> &neighbors = triangulation->getNeighborRanks();
    if(!ttk::isRunningWithMPI()) {
      return -1;
    }
    MPI_Datatype MPI_DT = getMPIType(static_cast<DT>(0));
    // we need unique tags for each rankToSend, otherwise messages might become
    // entangled
    int tagMultiplier = rankToSend + 1;
    int valuesTag = 103 * tagMultiplier;
    if(rankToSend == ttk::MPIrank_) {
      int neighborNumber = neighbors.size();
      const auto &ghostCellsPerOwner = triangulation->getGhostCellsPerOwner();
      // receive the scalar values
      for(int r = 0; r < neighborNumber; r++) {
        ttk::SimplexId nValues = ghostCellsPerOwner[neighbors[r]].size();
        std::vector<DT> receivedValues(nValues * dimensionNumber);
        if(nValues > 0) {
          MPI_Recv(receivedValues.data(), nValues * dimensionNumber, MPI_DT,
                   neighbors.at(r), valuesTag, communicator, MPI_STATUS_IGNORE);

          for(ttk::SimplexId i = 0; i < nValues; i++) {
            for(int j = 0; j < dimensionNumber; j++) {
              DT receivedVal = receivedValues[i * dimensionNumber + j];
              ttk::SimplexId globalId = ghostCellsPerOwner[neighbors[r]][i];
              ttk::SimplexId localId = triangulation->getCellLocalId(globalId);
              scalarArray[localId * dimensionNumber + j] = receivedVal;
            }
          }
        }
      }
    } else { // owner ranks
      // if rankToSend is not the neighbor of the current rank, we do not need
      // to do anything
      if(std::find(neighbors.begin(), neighbors.end(), rankToSend)
         != neighbors.end()) {
        // get the needed globalids from the triangulation
        const auto &ghostCellsForThisRank
          = triangulation->getRemoteGhostCells()[rankToSend];
        ttk::SimplexId nValues = ghostCellsForThisRank.size();
        if(nValues > 0) {
          // assemble the scalar values
          std::vector<DT> valuesToSend(nValues * dimensionNumber);
          for(ttk::SimplexId i = 0; i < nValues; i++) {
            for(int j = 0; j < dimensionNumber; j++) {
              ttk::SimplexId globalId = ghostCellsForThisRank[i];
              ttk::SimplexId localId = triangulation->getCellLocalId(globalId);
              valuesToSend[i * dimensionNumber + j]
                = scalarArray[localId * dimensionNumber + j];
            }
          }

          // send the scalar values
          MPI_Send(valuesToSend.data(), nValues * dimensionNumber, MPI_DT,
                   rankToSend, valuesTag, communicator);
        }
      }
    }

    return 0;
  }

  template <typename DT, typename triangulationType>
  int getGhostVertexScalars(DT *scalarArray,
                            const triangulationType *triangulation,
                            const int rankToSend,
                            MPI_Comm communicator,
                            const int dimensionNumber) {
    const std::vector<int> &neighbors = triangulation->getNeighborRanks();
    if(!ttk::isRunningWithMPI()) {
      return -1;
    }
    MPI_Datatype MPI_DT = getMPIType(static_cast<DT>(0));
    // we need unique tags for each rankToSend, otherwise messages might become
    // entangled
    int tagMultiplier = rankToSend + 1;
    int valuesTag = 103 * tagMultiplier;
    if(rankToSend == ttk::MPIrank_) {
      int neighborNumber = neighbors.size();
      const auto &ghostVerticesPerOwner
        = triangulation->getGhostVerticesPerOwner();
      // receive the scalar values
      for(int r = 0; r < neighborNumber; r++) {
        ttk::SimplexId nValues = ghostVerticesPerOwner[neighbors[r]].size();
        std::vector<DT> receivedValues(nValues * dimensionNumber);
        if(nValues > 0) {
          MPI_Recv(receivedValues.data(), nValues * dimensionNumber, MPI_DT,
                   neighbors.at(r), valuesTag, communicator, MPI_STATUS_IGNORE);
          for(ttk::SimplexId i = 0; i < nValues; i++) {
            for(int j = 0; j < dimensionNumber; j++) {
              DT receivedVal = receivedValues[i * dimensionNumber + j];
              ttk::SimplexId globalId = ghostVerticesPerOwner[neighbors[r]][i];
              ttk::SimplexId localId
                = triangulation->getVertexLocalId(globalId);
              scalarArray[localId * dimensionNumber + j] = receivedVal;
            }
          }
        }
      }
    } else { // owner ranks
      // if rankToSend is not the neighbor of the current rank, we do not need
      // to do anything
      if(std::find(neighbors.begin(), neighbors.end(), rankToSend)
         != neighbors.end()) {
        // get the needed globalids from the triangulation
        const auto &ghostVerticesForThisRank
          = triangulation->getRemoteGhostVertices()[rankToSend];
        ttk::SimplexId nValues = ghostVerticesForThisRank.size();
        if(nValues > 0) {
          // assemble the scalar values
          std::vector<DT> valuesToSend(nValues * dimensionNumber);
          for(ttk::SimplexId i = 0; i < nValues; i++) {
            for(int j = 0; j < dimensionNumber; j++) {
              ttk::SimplexId globalId = ghostVerticesForThisRank[i];
              ttk::SimplexId localId
                = triangulation->getVertexLocalId(globalId);
              valuesToSend[i * dimensionNumber + j]
                = scalarArray[localId * dimensionNumber + j];
            }
          }

          // send the scalar values
          MPI_Send(valuesToSend.data(), nValues * dimensionNumber, MPI_DT,
                   rankToSend, valuesTag, communicator);
        }
      }
    }

    return 0;
  }

  /**
   * @brief Request all ghost cell scalar data from one rank from their owning
   * ranks and get the data
   * this method is for usage without a triangulation, if a triangulation is
   * available, use getGhostCellScalars(), it is more performant when used
   * multiple times
   *
   * @param[out] scalarArray the scalar array which we want to fill and which is
   * filled on the other ranks
   * @param[in] getVertexRank lambda to get rank from vertex
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
  template <typename DT,
            typename IT,
            typename GVGID,
            typename GVR,
            typename GVLID>
  int getGhostDataScalarsWithoutTriangulation(DT *scalarArray,
                                              const GVR &getVertexRank,
                                              const GVGID &getVertexGlobalId,
                                              const GVLID &getVertexLocalId,
                                              const std::vector<int> &neighbors,
                                              const int rankToSend,
                                              const IT nVerts,
                                              MPI_Comm communicator,
                                              const int dimensionNumber) {
    int neighborNumber = neighbors.size();
    if(!ttk::isRunningWithMPI()) {
      return -1;
    }
    MPI_Datatype MPI_DT = getMPIType(static_cast<DT>(0));
    MPI_Datatype MPI_IT = getMPIType(static_cast<IT>(0));
    using globalIdType = decltype(getVertexGlobalId(0));
    MPI_Datatype MPI_GIT = getMPIType(static_cast<globalIdType>(0));
    // we need unique tags for each rankToSend, otherwise messages might become
    // entangled
    int tagMultiplier = rankToSend + 1;
    int amountTag = 101 * tagMultiplier;
    int idsTag = 102 * tagMultiplier;
    int valuesTag = 103 * tagMultiplier;
    if(rankToSend == ttk::MPIrank_) {
      // initialize the inner vectors with size 0
      std::vector<std::vector<globalIdType>> rankVectors(
        ttk::MPIsize_, std::vector<globalIdType>(0));
      // aggregate the needed ids

      for(IT i = 0; i < nVerts; i++) {
        if(ttk::MPIrank_ != getVertexRank(i)) {
          rankVectors[getVertexRank(i)].push_back(getVertexGlobalId(i));
        }
      }
      // send the amount of ids and the needed ids themselves
      for(int r = 0; r < neighborNumber; r++) {
        IT nValues = rankVectors[neighbors[r]].size();
        MPI_Send(&nValues, 1, MPI_IT, neighbors[r], amountTag, communicator);
        if(nValues > 0) {
          MPI_Send(rankVectors[neighbors[r]].data(), nValues, MPI_GIT,
                   neighbors[r], idsTag, communicator);
        }
      }
      // receive the scalar values
      for(int r = 0; r < neighborNumber; r++) {
        IT nValues = rankVectors[neighbors[r]].size();
        std::vector<DT> receivedValues(nValues * dimensionNumber);
        if(nValues > 0) {
          MPI_Recv(receivedValues.data(), nValues * dimensionNumber, MPI_DT,
                   neighbors[r], valuesTag, communicator, MPI_STATUS_IGNORE);
          for(IT i = 0; i < nValues; i++) {
            for(int j = 0; j < dimensionNumber; j++) {
              DT receivedVal = receivedValues[i * dimensionNumber + j];
              const auto globalId = rankVectors[neighbors[r]][i];
              IT localId = getVertexLocalId(globalId);
              scalarArray[localId * dimensionNumber + j] = receivedVal;
            }
          }
        }
      }
    } else { // owner ranks
      // if rankToSend is not the neighbor of the current rank, we do not need
      // to do anything
      if(std::find(neighbors.begin(), neighbors.end(), rankToSend)
         != neighbors.end()) {
        // receive the amount of ids and the needed ids themselves
        IT nValues;

        MPI_Recv(&nValues, 1, MPI_IT, rankToSend, amountTag, communicator,
                 MPI_STATUS_IGNORE);

        if(nValues > 0) {
          std::vector<globalIdType> receivedIds(nValues);
          MPI_Recv(receivedIds.data(), nValues, MPI_GIT, rankToSend, idsTag,
                   communicator, MPI_STATUS_IGNORE);

          // assemble the scalar values
          std::vector<DT> valuesToSend(nValues * dimensionNumber);
          for(IT i = 0; i < nValues; i++) {
            for(int j = 0; j < dimensionNumber; j++) {
              IT globalId = receivedIds[i];
              IT localId = getVertexLocalId(globalId);
              valuesToSend[i * dimensionNumber + j]
                = scalarArray[localId * dimensionNumber + j];
            }
          }

          // send the scalar values
          MPI_Send(valuesToSend.data(), nValues * dimensionNumber, MPI_DT,
                   rankToSend, valuesTag, communicator);
        }
      }
    }

    return 0;
  }

  /**
   * @brief get the neighbors of a rank
   *
   * @param[out] neighbors a set containing the ranks which are neighbors of
   * this rank
   * @param[in] getVertexRank lambda to get rank from vertex
   * @param[in] nVerts the number of vertices
   * @return 0 in case of success
   */
  template <typename IT, typename GVR>
  int preconditionNeighborsUsingRankArray(std::vector<int> &neighbors,
                                          const GVR &getVertexRank,
                                          const IT nVerts,
                                          MPI_Comm communicator) {
    std::unordered_set<int> neighborSet{};
    for(IT i = 0; i < nVerts; i++) {
      if(getVertexRank(i) != ttk::MPIrank_) {
        neighborSet.emplace(getVertexRank(i));
      }
    }
    std::vector<int> sendVector(neighborSet.begin(), neighborSet.end());
    int localSize = neighborSet.size();
    std::vector<int> sizes(ttk::MPIsize_);
    std::vector<int> displacements(ttk::MPIsize_);
    MPI_Gather(
      &localSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, communicator);
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
                rootVector.data(), sizes.data(), displacements.data(), MPI_INT,
                0, communicator);
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
      }
      // now we need to check for each rank if they are a neighbor of any other
      // rank. If so, we need to add those to the neighbors
      for(int i = 0; i < ttk::MPIsize_; i++) {
        for(int j = 0; j < ttk::MPIsize_; j++) {
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
    MPI_Scatter(
      sizes.data(), 1, MPI_INT, &receivedSize, 1, MPI_INT, 0, communicator);

    // and then the actual neighbors
    std::vector<int> receivedNeighbors(receivedSize);
    MPI_Scatterv(scatterVector.data(), sizes.data(), displacements.data(),
                 MPI_INT, receivedNeighbors.data(), receivedSize, MPI_INT, 0,
                 communicator);
    // then we turn the vector back into a set
    std::unordered_set<int> finalSet(
      receivedNeighbors.begin(), receivedNeighbors.end());
    neighborSet = finalSet;

    // We copy the set as a vector
    neighbors.clear();
    for(int neighbor : neighborSet) {
      neighbors.push_back(neighbor);
    }
    std::sort(neighbors.begin(), neighbors.end());

    return 0;
  }

  /**
   * @brief exchange all ghost cell information by calling getGhostCellScalars
   * for every rank
   *
   * @param[out] scalarArray the scalar array which we want to fill and which is
   * filled on the other ranks
   * @param[in] triangulation the triangulation for the data
   * @param[in] communicator the communicator over which the ranks are connected
   * (most likely ttk::MPIcomm_)
   * @return 0 in case of success
   */
  template <typename DT, typename triangulationType>
  int exchangeGhostCells(DT *scalarArray,
                         const triangulationType *triangulation,
                         MPI_Comm communicator,
                         const int dimensionNumber = 1) {
    if(!ttk::isRunningWithMPI()) {
      return -1;
    }
    if(!triangulation->hasPreconditionedDistributedCells()) {
      return -1;
    }
    for(int r = 0; r < ttk::MPIsize_; r++) {
      getGhostCellScalars<DT, triangulationType>(
        scalarArray, triangulation, r, communicator, dimensionNumber);
      MPI_Barrier(communicator);
    }
    return 0;
  }

  template <typename DT, typename triangulationType>
  int exchangeGhostVertices(DT *scalarArray,
                            const triangulationType *triangulation,
                            MPI_Comm communicator,
                            const int dimensionNumber = 1) {
    if(!ttk::isRunningWithMPI()) {
      return -1;
    }
    if(!triangulation->hasPreconditionedDistributedVertices()) {
      return -1;
    }
    for(int r = 0; r < ttk::MPIsize_; r++) {
      getGhostVertexScalars<DT, triangulationType>(
        scalarArray, triangulation, r, communicator, dimensionNumber);
      MPI_Barrier(communicator);
    }
    return 0;
  }

  /**
     * @brief exchange all ghost cell information by calling getGhostCellScalars
     * for every rank
     * this method is for usage without a triangulation, if a triangulation is
     available,
     * use exchangeGhostCells(), it is more performant when used multiple
     times

     * @param[out] scalarArray the scalar array which we want to fill and which
     is
     * filled on the other ranks
     * @param[in] getVertexRank lambda to get rank from vertex
     * @param[in] globalIds the global id array for the scalar data
     * @param[in] gidToLidMap a map which translates global ids to local,
     * rank-based ids
     * @param[in] nVerts number of vertices in the arrays
     * @param[in] communicator the communicator over which the ranks are
     connected
     * (most likely ttk::MPIcomm_)
     * @return 0 in case of success
     */
  template <typename DT,
            typename IT,
            typename GVGID,
            typename GVR,
            typename GVLID>
  int exchangeGhostDataWithoutTriangulation(DT *scalarArray,
                                            const GVR &getVertexRank,
                                            const GVGID &getVertexGlobalId,
                                            const GVLID &getVertexLocalId,
                                            const IT nVerts,
                                            MPI_Comm communicator,
                                            const std::vector<int> &neighbors,
                                            const int dimensionNumber = 1) {
    if(!ttk::isRunningWithMPI()) {
      return -1;
    }
    for(int r = 0; r < ttk::MPIsize_; r++) {
      getGhostDataScalarsWithoutTriangulation(
        scalarArray, getVertexRank, getVertexGlobalId, getVertexLocalId,
        neighbors, r, nVerts, communicator, dimensionNumber);
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

  void inline preconditionNeighborsUsingBoundingBox(
    double *boundingBox, std::vector<int> &neighbors) {

    std::vector<std::array<double, 6>> rankBoundingBoxes(ttk::MPIsize_);
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

    for(int i = 0; i < ttk::MPIsize_; i++) {
      if(i != ttk::MPIrank_) {
        double *theirBoundingBox = rankBoundingBoxes[i].data();
        if(checkForIntersection(boundingBox, theirBoundingBox)) {
          neighbors.push_back(i);
        }
      }
    }
  }

  /**
   * @brief produce the RankArray array, that stores rank ownership information
   *
   * @param[out] rankArray the owner array for the scalar data
   * @param[in] globalIds the global id array for the scalar data
   * @param[in] ghostCells the ghost array for the scalar data
   * @param[in] nVertices number of vertices in the arrays
   * @param[in] neighbors vector of neighboring ranks
   */
  void inline produceRankArray(std::vector<int> &rankArray,
                               const LongSimplexId *globalIds,
                               const unsigned char *ghostCells,
                               int nVertices,
                               double *boundingBox,
                               std::vector<int> &neighbors) {
    if(neighbors.empty()) {
      ttk::preconditionNeighborsUsingBoundingBox(boundingBox, neighbors);
    }
    MPI_Datatype MIT = ttk::getMPIType(ttk::SimplexId{});
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

  /**
   * @brief send the highest burstSize values and decrease the vector by that
   * amount, check if there are actually that many elements in the vector
   *
   * @param[out] outVector the vector which is returned
   * @param[in] values the data vector from which we take the elements
   * @param[in] burstSize the amount of elements extracted
   */
  template <typename DT, typename IT>
  void returnVectorForBurstsize(std::vector<value<DT, IT>> &outVector,
                                std::vector<value<DT, IT>> &values,
                                size_t burstSize) {

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
   * @brief receive elements from a rank and add them to the correct subvector
   *
   * @param[in] rankFrom the rank from which we receive
   * @param[in] structTag an MPI tag representing messages in which structs are
   * sent
   * @param[out] unsortedReceivedValues a vector of vectors representing
   * received values for each rank
   */
  template <typename DT, typename IT>
  void ReceiveAndAddToVector(
    int rankFrom,
    int structTag,
    std::vector<std::vector<value<DT, IT>>> &unsortedReceivedValues) {
    std::vector<value<DT, IT>> &receivedValues
      = unsortedReceivedValues[rankFrom];
    // probe to get the correct size to receive
    int amount;
    MPI_Status status;
    MPI_Probe(rankFrom, structTag * rankFrom, ttk::MPIcomm_, &status);
    MPI_Get_count(&status, MPI_CHAR, &amount);
    receivedValues.resize(amount / sizeof(value<DT, IT>), {0, 0});
    MPI_Recv(receivedValues.data(), amount, MPI_CHAR, rankFrom,
             structTag * rankFrom, ttk::MPIcomm_, MPI_STATUS_IGNORE);
  }

  /**
   * @brief get the current maximum value and works with it, if this leads to
   * the vector of one rank being empty, send the ordering to the rank and
   * request new values
   *
   * @param[in] intTag an MPI tag representing messages in which integers are
   * sent
   * @param[in] structTag an MPI tag representing messages in which structs are
   * sent
   * @param[in, out] currentOrder integer keeping track of the current order,
   * counting down
   * @param[in] burstSize number of values sent in one communication step
   * @param[in] MPI_IT MPI datatype representing integers
   * @param[out] processedValueCounter a counter keeping track of the number
   * of values already processed
   * @param[in, out] unsortedReceivedValues a vector of vectors representing
   * received values for each rank
   * @param[in, out] orderResendValues a vector of vectors of integers
   * representing globalIds and their orders for each rank
   * @param[in, out] orderedValuesForRank a vector representing integers
   * representing globalIds and their orders for this rank
   * @param[in, out] sortingValues vector of value structs with global ids
   * belonging to this rank
   */
  template <typename DT, typename IT>
  void getMax(int intTag,
              int structTag,
              IT &currentOrder,
              int burstSize,
              MPI_Datatype MPI_IT,
              IT &processedValueCounter,
              std::vector<std::vector<value<DT, IT>>> &unsortedReceivedValues,
              std::vector<std::vector<IT>> &orderResendValues,
              std::vector<IT> &orderedValuesForRank,
              std::vector<value<DT, IT>> &sortingValues) {
    // take the current maximum scalar over all ranks
    int rankIdOfMaxScalar = -1;
    DT maxScalar = std::numeric_limits<DT>::lowest();
    IT maxGId = -1;
    for(size_t i = 0; i < unsortedReceivedValues.size(); i++) {
      if(unsortedReceivedValues[i].size() > 0) {
        const auto &v = unsortedReceivedValues[i].back();
        if(v.scalar == maxScalar ? v.globalId > maxGId : v.scalar > maxScalar) {
          maxScalar = v.scalar;
          maxGId = v.globalId;
          rankIdOfMaxScalar = i;
        }
      }
    }

    value<DT, IT> currentValue
      = unsortedReceivedValues[rankIdOfMaxScalar].back();
    // we send the globalId and the the order via one send command,
    // therefore we need to check two concurrent values at the same time
    // later on
    orderResendValues[rankIdOfMaxScalar].push_back(currentValue.globalId);
    orderResendValues[rankIdOfMaxScalar].push_back(currentOrder);
    currentOrder--;
    processedValueCounter++;
    unsortedReceivedValues[rankIdOfMaxScalar].pop_back();
    if(unsortedReceivedValues[rankIdOfMaxScalar].size() == 0) {
      if(rankIdOfMaxScalar == 0) {
        // append the ordered values to the correct vector
        orderedValuesForRank.insert(
          orderedValuesForRank.end(),
          orderResendValues[rankIdOfMaxScalar].begin(),
          orderResendValues[rankIdOfMaxScalar].end());
        orderResendValues[rankIdOfMaxScalar].clear();

        if(!sortingValues.empty()) {
          std::vector<value<DT, IT>> ownValues;
          returnVectorForBurstsize<DT, IT>(ownValues, sortingValues, burstSize);

          unsortedReceivedValues[rankIdOfMaxScalar] = ownValues;
        }
      } else {
        // receive more values from rank, send ordering to the rank
        // send to the finished rank that we want more
        MPI_Send(orderResendValues[rankIdOfMaxScalar].data(),
                 orderResendValues[rankIdOfMaxScalar].size(), MPI_IT,
                 rankIdOfMaxScalar, intTag * rankIdOfMaxScalar, ttk::MPIcomm_);
        orderResendValues[rankIdOfMaxScalar].clear();
        // we receive more vals, if the size is still zero afterwards, we are
        // finished with this rank
        ReceiveAndAddToVector(
          rankIdOfMaxScalar, structTag, unsortedReceivedValues);
      }
    }
  }

  /**
   * @brief fill the order array based on global ids and their order values
   *
   * @param[in] nInts number of long ints, as globalid / order pairs
   * @param[in] orderedValuesForRank array of size nInts, the actual pairs,
   * [i] = gId, [i+1] = order
   * @param[in] gidToLidMap map which maps scalar values to a defined order
   * @param[out] order array of size nVerts, computed order of vertices, this
   * procedure doesn't fill it completely, ghostcells are missing
   * @param[in] nThreads number of parallel threads
   */
  template <typename IT, typename GVLID>
  void buildArrayForReceivedData(const size_t nInts,
                                 const IT *const orderedValuesForRank,
                                 const GVLID &getVertexLocalId,
                                 SimplexId *const order,
                                 const int nThreads) {

    TTK_FORCE_USE(nThreads);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nInts; i += 2) {
      order[getVertexLocalId(orderedValuesForRank[i])]
        = orderedValuesForRank[i + 1];
    }
  }

  /**
   * @brief creates a vector of value structs based on the given pointers, only
   * takes values of vertices which belong to the current rank, so no ghost
   * cells
   *
   * @param[out] valuesToSortVector vector of value structs with global ids
   * belonging to this rank
   * @param[out] gidToLidMap map mapping global ids to local rank ids
   * @param[in] nVerts number of vertices
   * @param[in] scalars the scalar data array
   * @param[in] globalIds the global id array for the scalar data
   * @param[in] getVertexRank lambda to get rank from vertex
   */
  template <typename DT, typename IT, typename GVGID, typename GVR>
  void populateVector(std::vector<value<DT, IT>> &valuesToSortVector,
                      const size_t nVerts,
                      const DT *const scalars,
                      const GVGID &getVertexGlobalId,
                      const GVR &getVertexRank) {
    for(size_t i = 0; i < nVerts; i++) {
      IT globalId = getVertexGlobalId(i);
      if(getVertexRank(i) == ttk::MPIrank_) {
        valuesToSortVector.emplace_back(scalars[i], globalId);
      }
    }
  }

  /**
   * @brief Sort vertices according to scalars disambiguated by global ids
   *
   * @param[in, out] values vector of size nVerts, structs with a scalar value
   * and a global id
   * @param[in] nThreads number of parallel threads
   */
  template <typename DT, typename IT>
  void sortVerticesDistributed(std::vector<value<DT, IT>> &values,
                               const int nThreads) {

    TTK_FORCE_USE(nThreads);

    TTK_PSORT(nThreads, values.begin(), values.end(),
              [](value<DT, IT> v1, value<DT, IT> v2) {
                return (v1.scalar < v2.scalar)
                       || (v1.scalar == v2.scalar && v1.globalId < v2.globalId);
              });
  }

  /**
   * @brief produce the orderArray, that defines a unique ordering based on the
   * scalar values and the global ids
   *
   * @param[out] orderArray the order array for the scalar data
   * @param[in] scalarArray the scalar data array
   * @param[in] globalIds the global id array for the scalar data
   * @param[in] getVertexRank lambda to get rank from vertex
   * @param[in] nVerts number of vertices in the arrays
   * @param[in] burstSize number of values sent in one communication step
   */
  template <typename DT, typename GVGID, typename GVR, typename GVLID>
  void produceOrdering(SimplexId *orderArray,
                       const DT *scalarArray,
                       const GVGID &getVertexGlobalId,
                       const GVR &getVertexRank,
                       const GVLID &getVertexLocalId,
                       const size_t nVerts,
                       const int burstSize,
                       std::vector<int> &neighbors) {
    int intTag = 101;
    int structTag = 102;
    if(neighbors.empty()) {
      ttk::preconditionNeighborsUsingRankArray(
        neighbors, getVertexRank, nVerts, ttk::MPIcomm_);
    }
    MPI_Barrier(ttk::MPIcomm_);

    using IT = decltype(getVertexGlobalId(0));
    MPI_Datatype MPI_IT = ttk::getMPIType(static_cast<IT>(0));

    ttk::Timer fillAndSortTimer;
    std::vector<value<DT, IT>> sortingValues;
    populateVector(
      sortingValues, nVerts, scalarArray, getVertexGlobalId, getVertexRank);

    // sort the scalar array distributed first by the scalar value itself,
    // then by the global id
    sortVerticesDistributed<DT, IT>(sortingValues, ttk::globalThreadNumber_);

    // when all are done sorting, rank 0 requests the highest values and
    // merges them
    MPI_Barrier(ttk::MPIcomm_);

    std::vector<IT> orderedValuesForRank;
    IT processedValueCounter = 0;
    ttk::Timer mergeTimer;
    IT localSize = sortingValues.size();
    IT totalSize;
    // get the complete size  of the dataset by summing up the local sizes
    MPI_Reduce(&localSize, &totalSize, 1, MPI_IT, MPI_SUM, 0, ttk::MPIcomm_);
    if(ttk::MPIrank_ == 0) {
      IT currentOrder = totalSize - 1;
      std::vector<std::vector<value<DT, IT>>> unsortedReceivedValues;
      unsortedReceivedValues.resize(ttk::MPIsize_);
      std::vector<std::vector<IT>> orderResendValues;
      orderResendValues.resize(ttk::MPIsize_);
      for(int i = 0; i < ttk::MPIsize_; i++) {
        orderResendValues[i].reserve(burstSize);
      }

      // receive the first batch of values
      for(int i = 0; i < ttk::MPIsize_; i++) {
        if(i == 0) {
          std::vector<value<DT, IT>> ownValues;
          returnVectorForBurstsize<DT, IT>(ownValues, sortingValues, burstSize);
          unsortedReceivedValues[i] = ownValues;
        } else {
          ReceiveAndAddToVector<DT, IT>(i, structTag, unsortedReceivedValues);
        }
      }
      while(processedValueCounter < totalSize) {
        getMax<DT, IT>(intTag, structTag, currentOrder, burstSize, MPI_IT,
                       processedValueCounter, unsortedReceivedValues,
                       orderResendValues, orderedValuesForRank, sortingValues);
      }

    } else { // other Ranks
      // send the next burstsize values and then wait for an answer from the
      // root rank
      while(!sortingValues.empty()) {
        std::vector<value<DT, IT>> sendValues;
        returnVectorForBurstsize<DT, IT>(sendValues, sortingValues, burstSize);
        int size = sendValues.size();
        MPI_Send(sendValues.data(), size * sizeof(value<DT, IT>), MPI_CHAR, 0,
                 structTag * ttk::MPIrank_, ttk::MPIcomm_);
        std::vector<IT> receivedValues;

        // be prepared to receive burstsize of elements, resize after
        // receiving to the correct size
        receivedValues.resize(burstSize * 2);
        MPI_Status status;
        int amount;
        MPI_Recv(receivedValues.data(), burstSize * 2, MPI_IT, 0,
                 intTag * ttk::MPIrank_, ttk::MPIcomm_, &status);
        MPI_Get_count(&status, MPI_IT, &amount);

        receivedValues.resize(amount);
        orderedValuesForRank.insert(orderedValuesForRank.end(),
                                    receivedValues.begin(),
                                    receivedValues.end());
      }
      // afterwards send once a message of length 0 to root to show that we
      // are done
      MPI_Send(sortingValues.data(), 0, MPI_CHAR, 0, structTag * ttk::MPIrank_,
               ttk::MPIcomm_);
    }

    // all ranks do the following
    MPI_Barrier(ttk::MPIcomm_);

    ttk::Timer orderTimer;
    buildArrayForReceivedData<IT>(orderedValuesForRank.size(),
                                  orderedValuesForRank.data(), getVertexLocalId,
                                  orderArray, ttk::globalThreadNumber_);

    // we receive the values at the ghostcells through the abstract
    // exchangeGhostCells method
    ttk::exchangeGhostDataWithoutTriangulation<ttk::SimplexId, IT>(
      orderArray, getVertexRank, getVertexGlobalId, getVertexLocalId, nVerts,
      ttk::MPIcomm_, neighbors);
  }

  /**
   * @brief sends the content of a vector to a neighbor.
   *
   * @tparam dataType data type of the vector
   * @param sendBuffer vector to send
   * @param messageType MPI data type of the vector
   * @param neighbor rank of the process to send the vector to
   */

  template <typename dataType>
  void sendVector(std::vector<dataType> &sendBuffer,
                  MPI_Datatype messageType,
                  int neighbor) {
    ttk::SimplexId dataSize = sendBuffer.size();
    MPI_Send(&dataSize, 1, getMPIType(dataSize), neighbor, ttk::MPIrank_,
             ttk::MPIcomm_);
    MPI_Send(sendBuffer.data(), dataSize, messageType, neighbor, ttk::MPIrank_,
             ttk::MPIcomm_);
  }

  /**
   * @brief receives the content of a vector from a neighbor.
   *
   * @tparam dataType data type of the vector
   * @param receiveBuffer vector to receive
   * @param messageType MPI data type of the vector
   * @param neighbor rank of the process to receive the vector from
   */

  template <typename dataType>
  void recvVector(std::vector<dataType> &receiveBuffer,
                  ttk::SimplexId &recvMessageSize,
                  MPI_Datatype messageType,
                  int neighbor) {
    MPI_Recv(&recvMessageSize, 1, getMPIType(recvMessageSize), neighbor,
             neighbor, ttk::MPIcomm_, MPI_STATUS_IGNORE);
    receiveBuffer.resize(recvMessageSize);
    MPI_Recv(receiveBuffer.data(), recvMessageSize, messageType, neighbor,
             neighbor, ttk::MPIcomm_, MPI_STATUS_IGNORE);
  }

} // namespace ttk

#endif // TTK_ENABLE_MPI
