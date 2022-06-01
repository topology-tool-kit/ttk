/// \ingroup base
/// \class ttk::PathCompression
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date 2022.
///
/// This module defines the %PathCompression class that computes
/// for each vertex of a triangulation the vertices to which ascending and
/// descending manifolds it belongs.
///
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <limits.h>
#include <map>
#include <stdint.h>
#include <unordered_map>

namespace ttk {

  /**
   * The PathCompression class provides methods to compute for
   * each vertex of a triangulation the vertices to which ascending and
   * descending manifolds it belongs.
   */
  class PathCompression : virtual public Debug {

  public:
    struct globalIdOwner {
      ttk::SimplexId globalId;
      int ownerRank;
      ttk::SimplexId ascendingTarget = -1;
      ttk::SimplexId descendingTarget = -1;

      globalIdOwner(ttk::SimplexId _globalId = -1, int _owner = -1)
        : globalId(_globalId), ownerRank(_owner) {
      }

      bool operator<(const globalIdOwner &other) const {
        return globalId < other.globalId;
      }

      bool operator==(const globalIdOwner &other) const {
        return globalId == other.globalId;
      }
    };
    PathCompression();

    int preconditionTriangulation(AbstractTriangulation *triangulation) {
      // Pre-condition functions.
      if(triangulation) {
        triangulation->preconditionVertexNeighbors();
      }
      return 0;
    }

    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computeCompression(ttk::SimplexId *descendingManifold,
                           ttk::SimplexId *ascendingManifold,
                           const dataType *inputData,
                           const triangulationType *triangulation,
                           const int *rankArray = nullptr,
                           const ttk::SimplexId *globalIds = nullptr) const {
      // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
      });
      this->printMsg(ttk::debug::Separator::L1);

      // -----------------------------------------------------------------------
      // Compute Vertex Minima
      // -----------------------------------------------------------------------
      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing compression",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_);

        ttk::SimplexId nVertices = triangulation->getNumberOfVertices();

        std::vector<ttk::SimplexId> previousDesc(nVertices);
        std::vector<ttk::SimplexId> currentDesc(nVertices);
        std::vector<ttk::SimplexId> previousAsc(nVertices);
        std::vector<ttk::SimplexId> currentAsc(nVertices);
        bool useMPI = false;
        std::vector<globalIdOwner> foreignVertices;

        TTK_FORCE_USE(useMPI);
        TTK_FORCE_USE(foreignVertices);
        TTK_FORCE_USE(rankArray);
        TTK_FORCE_USE(globalIds);

#if TTK_ENABLE_MPI
        if(ttk::isRunningWithMPI() && rankArray != nullptr
           && globalIds != nullptr) {
          useMPI = true;
        }

#endif

        // for the first step we initialize each vertex with the id of their
        // largest / smallest neighbor. Afterwards we only compare the arrays
        for(ttk::SimplexId i = 0; i < nVertices; i++) {
          int nNeighbors = triangulation->getVertexNeighborNumber(i);
          ttk::SimplexId neighborId;
          dataType smallest = inputData[i];
          dataType largest = inputData[i];
          // if there is no larger / smaller neighbor, the vertex points to
          // itself and is therefore a maximum / minimum we do not need to check
          // for equality, because we use the order array
          previousDesc[i] = i;
          previousAsc[i] = i;
          // if the vertex belongs to ourselves, we don't need to strictly point
          // to ourselves, but to the largest neighbor
#ifdef TTK_ENABLE_MPI
          if(useMPI) {
            if(rankArray[i] == ttk::MPIrank_) {
              for(int j = 0; j < nNeighbors; j++) {
                triangulation->getVertexNeighbor(i, j, neighborId);
                // and for the largest neighbor to get to the ascending manifold
                if(inputData[neighborId] > largest) {
                  previousAsc[i] = neighborId;
                  largest = inputData[neighborId];
                }
                // we're checking for the smallest neighbor to get the
                // descending manifold
                if(inputData[neighborId] < smallest) {
                  previousDesc[i] = neighborId;
                  smallest = inputData[neighborId];
                }
              }
            } else {
              globalIdOwner GIO = {globalIds[i], rankArray[i]};
              foreignVertices.push_back(GIO);
            }
          }
#else
          for(int j = 0; j < nNeighbors; j++) {
            triangulation->getVertexNeighbor(i, j, neighborId);
            // and for the largest neighbor to get to the ascending manifold
            if(inputData[neighborId] > largest) {
              previousAsc[i] = neighborId;
              largest = inputData[neighborId];
            }
            // we're checking for the smallest neighbor to get the descending
            // manifold
            if(inputData[neighborId] < smallest) {
              previousDesc[i] = neighborId;
              smallest = inputData[neighborId];
            }
          }
#endif
        }

        this->printMsg("Initializing values done");
        // now we swap between the two arrays until nothing changes anymore ergo
        // all paths are finished
        int step = 0;
        bool same = false;
        while(!same) {
          this->printMsg("Running Step " + std::to_string(step));
          same = true;
          if(step % 2 == 0) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
            for(ttk::SimplexId i = 0; i < nVertices; i++) {
              int nextDesc = previousDesc[previousDesc[i]];
              int nextAsc = previousAsc[previousAsc[i]];
              if(nextDesc != currentDesc[i]) {
                currentDesc[i] = nextDesc;
                same = false;
              }
              if(nextAsc != currentAsc[i]) {
                currentAsc[i] = nextAsc;
                same = false;
              }
            }
          } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
            for(ttk::SimplexId i = 0; i < nVertices; i++) {
              int nextDesc = currentDesc[currentDesc[i]];
              int nextAsc = currentAsc[currentAsc[i]];
              if(nextDesc != previousDesc[i]) {
                previousDesc[i] = nextDesc;
                same = false;
              }
              if(nextAsc != previousAsc[i]) {
                previousAsc[i] = nextAsc;
                same = false;
              }
            }
          }

          step++;
        }

#ifdef TTK_ENABLE_MPI
        // now we need to transform local ids into global ids to correctly work
        // over all ranks
        if(useMPI) {
          std::unordered_map<SimplexId, SimplexId> gIdTolIdMap;
          triangulation->getVertexGlobalIdMap(gIdTolIdMap);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
          for(ttk::SimplexId i = 0; i < nVertices; i++) {
            currentDesc[i] = globalIds[currentDesc[i]];
            currentAsc[i] = globalIds[currentAsc[i]];
          }
          this->printMsg("Rank " + std::to_string(ttk::MPIrank_)
                           + ", finished own values in Step "
                           + std::to_string(step),
                         1, localTimer.getElapsedTime());

          // now we need to request the values we still need from other ranks
          // R0 builds up out transferance map over ranks
          MPI_Barrier(MPI_COMM_WORLD);
          // MPI_Reduce to get the size of everything we need
          int localSize = foreignVertices.size();
          this->printMsg("Localsize: " + std::to_string(localSize));
          int totalSize;
          MPI_Reduce(
            &localSize, &totalSize, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
          MPI_Bcast(&totalSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
          this->printMsg("Rank " + std::to_string(ttk::MPIrank_)
                         + " got the totalsize " + std::to_string(totalSize));
          std::vector<globalIdOwner> edgesWithTargets(totalSize);
          int sizes[ttk::MPIsize_];
          int displacements[ttk::MPIsize_];
          std::vector<globalIdOwner> sendValues;
          int receivedSize;
          std::vector<globalIdOwner> receivedIds;
          std::vector<globalIdOwner> allValuesFromRanks;
          if(ttk::MPIrank_ == 0) {
            std::vector<globalIdOwner> edges(totalSize);
            // construct the set with the globalids not owned by R0

            // first we use MPI_Gather to get the size of each rank to populate
            // the displacements
            MPI_Gather(
              &localSize, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
            for(int i = 0; i < ttk::MPIsize_; i++) {
              sizes[i] = sizes[i] * sizeof(globalIdOwner);
            }
            displacements[0] = 0;
            // build our displacements
            for(int i = 1; i < ttk::MPIsize_; i++) {
              displacements[i] = displacements[i - 1] + sizes[i - 1];
            }

            // then we use MPI_Gatherv to get the data that is needed
            // R0 gets all the needed global ids and their original owners,
            // which rank actually needs those gids doesn't matter, because in
            // the end we send the finished edges back
            MPI_Gatherv(foreignVertices.data(),
                        foreignVertices.size() * sizeof(globalIdOwner),
                        MPI_CHAR, edges.data(), sizes, displacements, MPI_CHAR,
                        0, MPI_COMM_WORLD);
            this->printMsg("R0 received " + std::to_string(edges.size())
                           + " ids which are needed");

            // we have all the gids which are needed by _some_ rank and the
            // owner of them, now we have to request from the owners to which
            // vertices these gids are pointing to build our map
            std::vector<std::vector<globalIdOwner>> valuesFromRanks;
            valuesFromRanks.resize(ttk::MPIsize_);
            for(globalIdOwner currentId : edges) {
              valuesFromRanks[currentId.ownerRank].push_back(currentId);
            }
            this->printMsg("R0 reordered Ids");

            for(int i = 0; i < ttk::MPIsize_; i++) {
              sizes[i] = valuesFromRanks[i].size() * sizeof(globalIdOwner);
            }
            displacements[0] = 0;
            // build our displacements
            for(int i = 1; i < ttk::MPIsize_; i++) {
              displacements[i] = displacements[i - 1] + sizes[i - 1];
            }
            // we turn our vector of vectors into a 1D vector to send it via
            // scatter
            for(auto &&v : valuesFromRanks) {
              allValuesFromRanks.insert(
                allValuesFromRanks.end(), v.begin(), v.end());
            }
          } else { // the other ranks
            // first send the number of ids this rank needs to the root, then
            // the ids themselves the NULL attributes are only relevant for the
            // root rank
            MPI_Gather(
              &localSize, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gatherv(foreignVertices.data(),
                        foreignVertices.size() * sizeof(globalIdOwner),
                        MPI_CHAR, NULL, NULL, NULL, MPI_CHAR, 0,
                        MPI_COMM_WORLD);
          }

          // we need to receive the results to which the gids are pointing from
          // the ranks and build our map
          // we broadcast sizes and displacements and gather all the results
          // with allgatherv
          MPI_Bcast(sizes, ttk::MPIsize_, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(displacements, ttk::MPIsize_, MPI_INT, 0, MPI_COMM_WORLD);
          this->printMsg("R" + std::to_string(ttk::MPIrank_)
                         + " received sizes and displacements from R0");

          receivedSize = sizes[ttk::MPIrank_] / sizeof(globalIdOwner);
          this->printMsg("R" + std::to_string(ttk::MPIrank_) + " owns "
                         + std::to_string(receivedSize) + " ids");

          // and then the actual gids
          receivedIds.resize(receivedSize);
          MPI_Scatterv(allValuesFromRanks.data(), sizes, displacements,
                       MPI_CHAR, receivedIds.data(),
                       receivedSize * sizeof(globalIdOwner), MPI_CHAR, 0,
                       MPI_COMM_WORLD);
          this->printMsg("R" + std::to_string(ttk::MPIrank_) + " received ids");

          // now we need to find to where the gids point and send the values
          // back to R0
          sendValues.resize(receivedSize);
          for(ttk::SimplexId i = 0; i < receivedSize; i++) {
            globalIdOwner currentVal = receivedIds[i];
            ttk::SimplexId lId = gIdTolIdMap.at(currentVal.globalId);
            currentVal.ascendingTarget = currentAsc[lId];
            currentVal.descendingTarget = currentDesc[lId];
            sendValues[i] = currentVal;
          }
          this->printMsg("R" + std::to_string(ttk::MPIrank_)
                         + " is done with their values");

          MPI_Allgatherv(sendValues.data(),
                         sendValues.size() * sizeof(globalIdOwner), MPI_CHAR,
                         edgesWithTargets.data(), sizes, displacements,
                         MPI_CHAR, MPI_COMM_WORLD);
          this->printMsg("R" + std::to_string(ttk::MPIrank_)
                         + " got the results");

          // now each rank has a vector consisting of gIds, the ranks to which
          // they belong and the ascending / descending target we have all the
          // information on R0 which we need to resolve any manifolds stretching
          // over multiple ranks
          std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdToAscendingMap;
          std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdToDescendingMap;

          for(size_t i = 0; i < edgesWithTargets.size(); i++) {
            globalIdOwner currentVal = edgesWithTargets[i];
            gIdToAscendingMap.insert(
              std::make_pair(currentVal.globalId, currentVal.ascendingTarget));
            gIdToDescendingMap.insert(
              std::make_pair(currentVal.globalId, currentVal.descendingTarget));
          }

          // now we need to check for graphs in the map and iteratively compress
          // them
          bool changed = true;
          while(changed) {
            changed = false;
            for(auto &it : gIdToAscendingMap) {
              if(gIdToAscendingMap.count(it.second) && (it.first != it.second)
                 && (gIdToAscendingMap[it.first]
                     != gIdToAscendingMap[it.second])) {
                gIdToAscendingMap[it.first] = gIdToAscendingMap[it.second];
                changed = true;
              }
            }

            for(auto &it : gIdToDescendingMap) {
              if(gIdToDescendingMap.count(it.second) && (it.first != it.second)
                 && (gIdToDescendingMap[it.first]
                     != gIdToDescendingMap[it.second])) {
                gIdToDescendingMap[it.first] = gIdToDescendingMap[it.second];
                changed = true;
              }
            }
          }

          // now each rank simply needs to walk over their vertices and replace
          // ones aiming to ghostcells with the correct ones from the map
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
          for(ttk::SimplexId i = 0; i < nVertices; i++) {
            ttk::SimplexId descVal = currentDesc[i];
            ttk::SimplexId gid = globalIds[i];
            ttk::SimplexId ascVal = currentAsc[i];
            if(gIdToDescendingMap.count(descVal)) {
              currentDesc[i] = gIdToDescendingMap[descVal];
            } else if(gIdToDescendingMap.count(gid)) {
              currentDesc[i] = gIdToDescendingMap[gid];
            }
            if(gIdToAscendingMap.count(ascVal)) {
              currentAsc[i] = gIdToAscendingMap[ascVal];
            } else if(gIdToAscendingMap.count(gid)) {
              currentAsc[i] = gIdToAscendingMap[gid];
            }
          }
        }
#endif // TTK_ENABLE_MPI

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(ttk::SimplexId i = 0; i < nVertices; i++) {
          descendingManifold[i] = currentDesc[i];
          ascendingManifold[i] = currentAsc[i];
        }

        // print the progress of the current subprocedure with elapsed time
        this->printMsg("Computed Compression",
                       1, // progress
                       localTimer.getElapsedTime(), this->threadNumber_);
      }

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

  }; // PathCompression class

} // namespace ttk
