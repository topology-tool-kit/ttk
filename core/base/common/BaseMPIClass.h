/// \ingroup base
/// \class ttk::BaseMPIClass
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date April 2022
///
/// \brief Base Class and utilities for MPI implementation.

#pragma once
#include <BaseClass.h>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#if TTK_ENABLE_MPI
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

namespace ttk {
  COMMON_EXPORTS extern int MPIrank_;

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
    int flag_i;
    MPI_Initialized(&flag_i);
    return flag_i;
  }

  template <typename IT>
  void getGhostCellInfo(std::vector<std::queue<IT>> &rankQueues,
                        const int *const rankArray,
                        const int rank,
                        const IT nVerts) {
    for(IT i = 0; i < nVerts; i++) {
      if(rank != rankArray[i]) {
        rankQueues[rankArray[i]].push(i);
      }
    }
  }

  // this is a collective operation, therefore all processeses
  // in the communicator need to call it for it to work correctly
  // otherwise it will freeze because it calls MPI_Send from rankToSend
  // and expects an MPI_Recv from all other ranks
  template <typename DT, typename IT>
  void getGhostCellScalars(DT *scalarArray,
                           const int *const rankArray,
                           const IT *const globalIds,
                           const std::unordered_map<IT, IT> gidToLidMap,
                           const std::unordered_set<int> neighbors,
                           const int rankToSend,
                           const IT nVerts,
                           const MPI_Comm communicator) {
    MPI_Datatype MPI_DT = getMPIType(static_cast<DT>(0));
    MPI_Datatype MPI_IT = getMPIType(static_cast<IT>(0));
    int rank;
    int nRanks;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &nRanks);
    // we need unique tags for each rankToSend, otherwise messages might become
    // entangled
    int tagMultiplier = rankToSend + 1;
    int amountTag = 101 * tagMultiplier;
    int idsTag = 102 * tagMultiplier;
    int valuesTag = 103 * tagMultiplier;
    if(rankToSend == rank) {
      // initialize the inner vectors with size 0
      std::vector<std::vector<IT>> rankVectors(nRanks, std::vector<IT>(0));
      // aggregate the needed ids
      for(IT i = 0; i < nVerts; i++) {
        if(rank != rankArray[i]) {
          rankVectors[rankArray[i]].push_back(globalIds[i]);
        }
      }

      // send the amount of ids and the needed ids themselves
      for(int r = 0; r < nRanks; r++) {
        if(rankToSend != r && neighbors.find(r) != neighbors.end()) {
          IT nValues = rankVectors[r].size();
          MPI_Send(&nValues, 1, MPI_IT, r, amountTag, communicator);
          if(nValues > 0) {
            MPI_Send(
              rankVectors[r].data(), nValues, MPI_IT, r, idsTag, communicator);
          }
        }
      }

      // receive the scalar values
      for(int r = 0; r < nRanks; r++) {
        if(rankToSend != r && neighbors.find(r) != neighbors.end()) {
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

    } else {
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
  }

  // this gets the neighbors for each rank by traversing the rankArray
  template <typename IT>
  std::unordered_set<int> getNeighbors(const int *const rankArray,
                                       const IT nVerts,
                                       const int ownRank) {
    std::unordered_set<int> neighbors;
    for(IT i = 0; i < nVerts; i++) {
      if(rankArray[i] != ownRank) {
        neighbors.emplace(rankArray[i]);
      }
    }
    return neighbors;
  }

  // this method exchanges all ghost cell information
  // it does so by calling getGhostCellScalars for every rank in the
  // communicator
  template <typename DT, typename IT>
  void exchangeGhostCells(DT *scalarArray,
                          const int *const rankArray,
                          const IT *const globalIds,
                          const std::unordered_map<IT, IT> gidToLidMap,
                          const IT nVerts,
                          const MPI_Comm communicator) {
    int nRanks;
    int rank;
    MPI_Comm_size(communicator, &nRanks);
    MPI_Comm_rank(communicator, &rank);
    auto neighbors = getNeighbors<IT>(rankArray, nVerts, rank);
    for(int r = 0; r < nRanks; r++) {
      getGhostCellScalars<DT, IT>(scalarArray, rankArray, globalIds,
                                  gidToLidMap, neighbors, r, nVerts,
                                  communicator);
      MPI_Barrier(communicator);
    }
  }

  class BaseMPIClass : public BaseClass {

  public:
    BaseMPIClass();
    virtual ~BaseMPIClass() = default;
  };
} // namespace ttk

#endif