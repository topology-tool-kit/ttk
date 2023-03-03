///
/// \ingroup base
/// \class ttk::Identifiers
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date October 2022.
///
/// This module defines the %Identifiers class that computes for each vertex of
/// a triangulation the average scalar value of itself and its direct neighbors.
///
///

#pragma once

// ttk common includes
#include <Debug.h>
#ifdef TTK_ENABLE_MPI
#include <Geometry.h>
#include <KDTree.h>
#include <map>
#endif

namespace ttk {

#ifdef TTK_ENABLE_MPI
  struct Point {
    double x;
    double y;
    double z;
    ttk::SimplexId localId;
  };

  struct Response {
    ttk::SimplexId id;
    ttk::LongSimplexId globalId;
  };

#endif

  /**
   * The Identifiers class provides methods to compute for each vertex and
   * cell of a data set a global id. It uses RankArray and outdated global
   * identifiers if they are defined to accelerate computation for data sets
   * of PolyData type.
   */
  class Identifiers : virtual public Debug {

  public:
    Identifiers();

    ~Identifiers() override = default;

  protected:
    ttk::SimplexId vertexNumber_{};
    ttk::SimplexId cellNumber_{};
    ttk::LongSimplexId *vertexIdentifiers_;
    ttk::LongSimplexId *cellIdentifiers_;
#ifdef TTK_ENABLE_MPI
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> *vertGtoL_;
    std::vector<int> *neighbors_;
    std::map<int, int> neighborToId_;
    int neighborNumber_;
    double *bounds_;
    int *dims_;
    double *spacing_;
    MPI_Datatype mpiIdType_;
    MPI_Datatype mpiResponseType_;
    MPI_Datatype mpiPointType_;
    int dimension_{};
    int hasSentData_{0};
    int *vertexRankArray_{nullptr};
    int *cellRankArray_{nullptr};
    unsigned char *vertGhost_{nullptr};
    unsigned char *cellGhost_{nullptr};
    std::vector<std::vector<ttk::SimplexId>> pointsToCells_;
    float *pointSet_;
    ttk::LongSimplexId *connectivity_;
    ttk::LongSimplexId *outdatedGlobalPointIds_{nullptr};
    ttk::LongSimplexId *outdatedGlobalCellIds_{nullptr};
    std::map<ttk::LongSimplexId, ttk::SimplexId> vertOutdatedGtoL_;
    std::map<ttk::LongSimplexId, ttk::SimplexId> cellOutdatedGtoL_;
    ttk::KDTree<float, std::array<float, 3>> kdt_;
#endif

  public:
    void setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
    }

    void setCellNumber(const SimplexId &cellNumber) {
      cellNumber_ = cellNumber;
    }

    void setVertexIdentifiers(ttk::LongSimplexId *vertexIdentifiers) {
      this->vertexIdentifiers_ = vertexIdentifiers;
    }

    void setCellIdentifiers(ttk::LongSimplexId *cellIdentifiers) {
      this->cellIdentifiers_ = cellIdentifiers;
    }

#ifdef TTK_ENABLE_MPI

    inline void setDomainDimension(const int &dimension) {
      dimension_ = dimension;
    }

    inline void setPointsToCells(
      const std::vector<std::vector<ttk::SimplexId>> &pointsToCells) {
      pointsToCells_ = pointsToCells;
    }

    void setVertexRankArray(int *vertexRankArray) {
      this->vertexRankArray_ = vertexRankArray;
    }

    void setCellRankArray(int *cellRankArray) {
      this->cellRankArray_ = cellRankArray;
    }

    void setVertGhost(unsigned char *vertGhost) {
      this->vertGhost_ = vertGhost;
    }

    void setCellGhost(unsigned char *cellGhost) {
      this->cellGhost_ = cellGhost;
    }

    void setBounds(double *bounds) {
      this->bounds_ = bounds;
    }

    void setVertGtoL(
      std::unordered_map<ttk::SimplexId, ttk::SimplexId> *vertGtoL) {
      this->vertGtoL_ = vertGtoL;
    }

    void setSpacing(double *spacing) {
      this->spacing_ = spacing;
    }

    void setDims(int *dims) {
      this->dims_ = dims;
    }

    void setPointSet(float *pointSet) {
      pointSet_ = pointSet;
    }

    void setConnectivity(ttk::LongSimplexId *connectivity) {
      connectivity_ = connectivity;
    }

    void setOutdatedGlobalPointIds(ttk::LongSimplexId *outdatedGlobalPointIds) {
      outdatedGlobalPointIds_ = outdatedGlobalPointIds;
    }

    void buildKDTree() {
      kdt_ = ttk::KDTree<float, std::array<float, 3>>(false, 2);
      kdt_.build(pointSet_, vertexNumber_, 3);
    }

    void setOutdatedGlobalCellIds(ttk::LongSimplexId *outdatedGlobalCellIds) {
      outdatedGlobalCellIds_ = outdatedGlobalCellIds;
    }

    void initializeMPITypes() {
      ttk::SimplexId id{-1};
      ttk::LongSimplexId longId{-1};
      // Initialize id type
      mpiIdType_ = getMPIType(id);

      // Initialize Point Type
      MPI_Datatype types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, mpiIdType_};
      int lengths[] = {1, 1, 1, 1};
      const long int mpi_offsets[]
        = {offsetof(Point, x), offsetof(Point, y), offsetof(Point, z),
           offsetof(Point, localId)};
      MPI_Type_create_struct(4, lengths, mpi_offsets, types, &mpiPointType_);
      MPI_Type_commit(&mpiPointType_);

      // Initialize Response Type
      MPI_Datatype typesResponse[] = {mpiIdType_, getMPIType(longId)};
      int lengthsResponse[] = {1, 1};
      const long int mpi_offsetsResponse[]
        = {offsetof(Response, id), offsetof(Response, globalId)};
      MPI_Type_create_struct(2, lengthsResponse, mpi_offsetsResponse,
                             typesResponse, &mpiResponseType_);
      MPI_Type_commit(&mpiResponseType_);
    }

    void inline findPoint(ttk::SimplexId &id, float x, float y, float z) {
      std::array<float, 3> coordinates = {x, y, z};
      std::vector<ttk::KDTree<float, std::array<float, 3>> *> neighbours;
      std::vector<float> costs;
      kdt_.getKClosest(1, coordinates, neighbours, costs);
      id = neighbours[0]->id_;
    }

    void initializeNeighbors(double *boundingBox,
                             std::vector<int> &neighborRanks) {
      if(neighborRanks.empty()) {
        preconditionNeighborsUsingBoundingBox(boundingBox, neighborRanks);
      }
      neighbors_ = &neighborRanks;
      neighborToId_.clear();
      neighborNumber_ = neighbors_->size();
      for(int i = 0; i < neighborNumber_; i++) {
        neighborToId_[neighbors_->at(i)] = i;
      }
    }

    /**
     * @brief Finds the local point closest to the point coordinates using a
     * kd-tree, if the coordinates of the received point are within the bounds
     * of the process. The global identifier is then stored in the vector
     * locatedSimplices for each thread. The data is then copied to a simple
     * vector to send it back to the neighbor.
     *
     * @param locatedSimplices for each thread, stores the global id and local
     * id of the located points.
     * @param receivedPoints points to be located received from a neighbor
     * @param recvMessageSize size of the receivedPoints vector
     * @param send_buf send buffer
     */

    void inline locatePoints(
#ifdef TTK_ENABLE_OPENMP
      std::vector<std::vector<Response>> &locatedSimplices,
#else
      std::vector<std::vector<Response>> &ttkNotUsed(locatedSimplices),
#endif
      std::vector<Point> &receivedPoints,
      ttk::SimplexId &recvMessageSize,
      std::vector<Response> &send_buf) {
      ttk::LongSimplexId globalId{-1};
      ttk::SimplexId id{-1};
      send_buf.clear();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) firstprivate(globalId, id) \
  shared(locatedSimplices)
#endif
      for(int n = 0; n < recvMessageSize; n++) {
        if(bounds_[0] <= receivedPoints[n].x
           && bounds_[1] >= receivedPoints[n].x
           && bounds_[2] <= receivedPoints[n].y
           && bounds_[3] >= receivedPoints[n].y
           && bounds_[4] <= receivedPoints[n].z
           && bounds_[5] >= receivedPoints[n].z) {
          this->findPoint(
            id, receivedPoints[n].x, receivedPoints[n].y, receivedPoints[n].z);
          if((vertexRankArray_ != nullptr
              && vertexRankArray_[id] == ttk::MPIrank_)
             || (vertGhost_ != nullptr && vertGhost_[id] == 0)) {
            globalId = vertexIdentifiers_[id];
            if(globalId >= 0) {
#ifdef TTK_ENABLE_OPENMP
              locatedSimplices[omp_get_thread_num()].push_back(
                Response{receivedPoints[n].localId, globalId});
#else
              send_buf.push_back(Response{receivedPoints[n].localId, globalId});
#endif
            }
          }
        }
      }
#ifdef TTK_ENABLE_OPENMP
      for(int n = 0; n < threadNumber_; n++) {
        send_buf.insert(send_buf.end(), locatedSimplices[n].begin(),
                        locatedSimplices[n].end());
        locatedSimplices[n].clear();
      }
#endif
    }
    /**
     * @brief Identifies the local point corresponding to the received point
     * using its outdated global identifier. The valid global identifier is then
     * stored in the vector locatedSimplices for each thread. The data is then
     * copied to a simple vector to send it back to the neighbor.
     *
     * @param locatedSimplices for each thread, stores the global id and local
     * id of the located points.
     * @param receivedOutdatedGlobalIds outdated global ids to be identified
     * @param recvMessageSize size of the receivedOutdatedGlobalIds vector
     * @param send_buf send buffer
     */

    void inline identifyPoints(
#ifdef TTK_ENABLE_OPENMP
      std::vector<std::vector<Response>> &locatedSimplices,
#else
      std::vector<std::vector<Response>> &ttkNotUsed(locatedSimplices),
#endif
      std::vector<ttk::SimplexId> &receivedOutdatedGlobalIds,
      ttk::SimplexId &recvMessageSize,
      std::vector<Response> &send_buf) {
      send_buf.clear();
      ttk::SimplexId globalId{-1};
      std::map<ttk::LongSimplexId, ttk::SimplexId>::iterator search;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(int n = 0; n < recvMessageSize; n++) {
        search = vertOutdatedGtoL_.find(receivedOutdatedGlobalIds[n]);
        if(search != vertOutdatedGtoL_.end()) {
          globalId = vertexIdentifiers_[search->second];
          if(globalId >= 0) {
#ifdef TTK_ENABLE_OPENMP
            locatedSimplices[omp_get_thread_num()].push_back(
              Response{receivedOutdatedGlobalIds[n], globalId});
#else
            send_buf.push_back(
              Response{receivedOutdatedGlobalIds[n], globalId});
#endif
          }
        }
      }
#ifdef TTK_ENABLE_OPENMP
      for(int n = 0; n < threadNumber_; n++) {
        send_buf.insert(send_buf.end(), locatedSimplices[n].begin(),
                        locatedSimplices[n].end());
        locatedSimplices[n].clear();
      }
#endif
    }
    /**
     * @brief Finds the local cell for which the vertices have the same global
     * ids as the vertices of the received cell. The global identifier is then
     * stored in the vector locatedSimplices for each thread. The data is then
     * copied to a simple vector to send it back to the neighbor.
     *
     * @param locatedSimplices for each thread, stores the global id and local
     * id of the located cells.
     * @param receivedCells cells to be located received from a neighbor
     * @param recvMessageSize size of the receivedCells vector
     * @param send_buf send buffer
     */

    void inline locateCells(
#ifdef TTK_ENABLE_OPENMP
      std::vector<std::vector<Response>> &locatedSimplices,
#else
      std::vector<std::vector<Response>> &ttkNotUsed(locatedSimplices),
#endif
      std::vector<ttk::SimplexId> &receivedCells,
      ttk::SimplexId &recvMessageSize,
      std::vector<Response> &send_buf) {
      int id{-1};
      send_buf.clear();
      std::unordered_map<ttk::SimplexId, ttk::SimplexId>::iterator search;
      std::vector<ttk::SimplexId> localPointIds;
      localPointIds.reserve(dimension_ + 1);
      size_t expectedSize = static_cast<size_t>(dimension_) + 1;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  firstprivate(search, localPointIds)
#endif
      for(int n = 0; n < recvMessageSize; n += dimension_ + 2) {
        localPointIds.clear();
        for(int k = 1; k < dimension_ + 2; k++) {
          search = vertGtoL_->find(receivedCells[n + k]);
          if(search != vertGtoL_->end()) {
            localPointIds.push_back(search->second);
          } else {
            break;
          }
        }
        if(localPointIds.size() == expectedSize) {
          bool foundIt = false;
          int k = 0;
          int l;
          int m = 0;
          while(!foundIt && m < dimension_ + 1) {
            int size = pointsToCells_[localPointIds[m]].size();
            k = 0;
            while(!foundIt && k < size) {
              l = 0;
              while(l < dimension_ + 1) {
                id = connectivity_[pointsToCells_[localPointIds[m]][k]
                                     * (dimension_ + 1)
                                   + l];
                auto it = find(localPointIds.begin(), localPointIds.end(), id);
                if(it == localPointIds.end()) {
                  break;
                }
                l++;
              }
              if(l == dimension_ + 1) {
                foundIt = true;
#ifdef TTK_ENABLE_OPENMP
                locatedSimplices[omp_get_thread_num()].push_back(Response{
                  receivedCells[n],
                  cellIdentifiers_[pointsToCells_[localPointIds[m]][k]]});
#else
                send_buf.push_back(Response{
                  receivedCells[n],
                  cellIdentifiers_[pointsToCells_[localPointIds[m]][k]]});
#endif
              }
              k++;
            }
            m++;
          }
        }
      }
#ifdef TTK_ENABLE_OPENMP
      for(int n = 0; n < threadNumber_; n++) {
        send_buf.insert(send_buf.end(), locatedSimplices[n].begin(),
                        locatedSimplices[n].end());
        locatedSimplices[n].clear();
      }
#endif
    }
    /**
     * @brief Identifies the local cell corresponding to the received cell using
     * its outdated global identifier. The valid global identifier is then
     * stored in the vector locatedSimplices for each thread. The data is then
     * copied to a simple vector to send it back to the neighbor.
     *
     * @param locatedSimplices for each thread, stores the global id and local
     * id of the located points.
     * @param receivedOutdatedGlobalIds outdated global ids to be identified
     * @param recvMessageSize size of the receivedOutdatedGlobalIds vector
     * @param send_buf send buffer
     */

    void inline identifyCells(
#ifdef TTK_ENABLE_OPENMP
      std::vector<std::vector<Response>> &locatedSimplices,
#else
      std::vector<std::vector<Response>> &ttkNotUsed(locatedSimplices),
#endif
      std::vector<ttk::SimplexId> &receivedOutdatedGlobalIds,
      ttk::SimplexId &recvMessageSize,
      std::vector<Response> &send_buf) {
      ttk::LongSimplexId globalId{-1};
      std::map<ttk::LongSimplexId, ttk::SimplexId>::iterator search;
      send_buf.clear();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(int n = 0; n < recvMessageSize; n++) {
        search = cellOutdatedGtoL_.find(receivedOutdatedGlobalIds[n]);
        if(search != cellOutdatedGtoL_.end()) {
          globalId = cellIdentifiers_[search->second];
          if(globalId >= 0) {
#ifdef TTK_ENABLE_OPENMP
            locatedSimplices[omp_get_thread_num()].push_back(
              Response{receivedOutdatedGlobalIds[n], globalId});
#else
            send_buf.push_back(
              Response{receivedOutdatedGlobalIds[n], globalId});
#endif
          }
        }
      }
#ifdef TTK_ENABLE_OPENMP
      for(int n = 0; n < threadNumber_; n++) {
        send_buf.insert(send_buf.end(), locatedSimplices[n].begin(),
                        locatedSimplices[n].end());
        locatedSimplices[n].clear();
      }
#endif
    }

    template <typename dataType>
    void sendToAllNeighbors(std::vector<dataType> &vectorToSend,
                            MPI_Datatype messageType) const {
      for(int j = 0; j < neighborNumber_; j++) {
        sendVector<dataType>(vectorToSend, messageType, neighbors_->at(j));
      }
    }

    template <typename dataType>
    void sendToAllNeighborsUsingRankArray(
      std::vector<std::vector<dataType>> &vectorToSend,
      MPI_Datatype messageType) {

      for(int j = 0; j < neighborNumber_; j++) {
        sendVector<dataType>(vectorToSend[neighborToId_[neighbors_->at(j)]],
                             messageType, neighbors_->at(j));
      }
    }

    /**
     * @brief Computes the number of vertices owned by the current process and
     * stores ghost points. An exclusive prefix sum is performed to compute the
     * offset of each process and the value is used to generate the global ids.
     * Then global ids are generated.
     *
     * @param vertGhostCoordinatesPerRank similar to `vertGhostCoordinates`,
     * useful when RankArray is defined for vertices.
     * vertGhostCoordinatesPerRank[i] stores the vertices of process
     * neighbors_[i]
     * @param vertGhostCoordinates stores a Point struct for each ghost point
     * @param vertGhostGlobalIdsPerRank similar to `vertGhostGlobalIds`, useful
     * when RankArray and outdatedGlobalPointIds are defined for points.
     * vertGhostGlobalIdsPerRank[i] stores the vertex data of process
     * neighbors_[i]
     * @param vertGhostGlobalIds stores the global id of a ghost point, in case
     * outdatedGlobalPointIds_ is defined.
     */
    void generateGlobalIds(
      std::vector<std::vector<Point>> &vertGhostCoordinatesPerRank,
      std::vector<Point> &vertGhostCoordinates,
      std::vector<std::vector<ttk::SimplexId>> &vertGhostGlobalIdsPerRank,
      std::vector<ttk::SimplexId> &vertGhostGlobalIds) {
      ttk::SimplexId realVertexNumber = vertexNumber_;
      ttk::SimplexId realCellNumber = cellNumber_;
      float p[3];

      // Computes the number of vertices owned by the current process
      // If the vertex is not owned, it will be added to the vector of
      // ghosts.

      if(vertexRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
          if(vertexRankArray_[i] != ttk::MPIrank_) {
            realVertexNumber--;
            if(outdatedGlobalPointIds_ == nullptr) {
              p[0] = pointSet_[i * 3];
              p[1] = pointSet_[i * 3 + 1];
              p[2] = pointSet_[i * 3 + 2];
              vertGhostCoordinatesPerRank[neighborToId_[vertexRankArray_[i]]]
                .push_back(Point{p[0], p[1], p[2], i});
            } else {
              vertGhostGlobalIdsPerRank[neighborToId_[vertexRankArray_[i]]]
                .push_back(
                  static_cast<ttk::SimplexId>(outdatedGlobalPointIds_[i]));
            }
          }
        }
      } else {
        for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
          if(vertGhost_[i] != 0) {
            realVertexNumber--;
            if(outdatedGlobalPointIds_ == nullptr) {

              p[0] = pointSet_[i * 3];
              p[1] = pointSet_[i * 3 + 1];
              p[2] = pointSet_[i * 3 + 2];
              vertGhostCoordinates.push_back(Point{p[0], p[1], p[2], i});
            } else {
              vertGhostGlobalIds.push_back(
                static_cast<ttk::SimplexId>(outdatedGlobalPointIds_[i]));
            }
          }
        }
      }

      // The number of cells owned by the process is computed
      if(cellRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellRankArray_[i] != ttk::MPIrank_) {
            realCellNumber--;
          }
        }
      } else {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellGhost_[i] != 0) {
            realCellNumber--;
          }
        }
      }

      ttk::SimplexId vertIndex;
      ttk::SimplexId cellIndex;

      // Perform exclusive prefix sum to find local offset for vertices and
      // cells
      MPI_Exscan(
        &realVertexNumber, &vertIndex, 1, mpiIdType_, MPI_SUM, ttk::MPIcomm_);
      MPI_Exscan(
        &realCellNumber, &cellIndex, 1, mpiIdType_, MPI_SUM, ttk::MPIcomm_);

      // Rank 0 received garbage values, it is replaced by the correct offset
      // (always 0)
      if(ttk::MPIrank_ == 0) {
        vertIndex = 0;
        cellIndex = 0;
      }

      // Generate global ids for vertices
      if(vertexRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
          if(vertexRankArray_[i] == ttk::MPIrank_) {
            vertexIdentifiers_[i] = vertIndex;
            (*vertGtoL_)[vertIndex] = i;
            vertIndex++;
          }
          if(outdatedGlobalPointIds_ != nullptr) {
            vertOutdatedGtoL_[outdatedGlobalPointIds_[i]] = i;
          }
        }
      } else {
        for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
          if(vertGhost_[i] == 0) {
            vertexIdentifiers_[i] = vertIndex;
            (*vertGtoL_)[vertIndex] = i;
            vertIndex++;
          }
          if(outdatedGlobalPointIds_ != nullptr) {
            vertOutdatedGtoL_[outdatedGlobalPointIds_[i]] = i;
          }
        }
      }

      // Generate global ids for cells
      if(cellRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellRankArray_[i] == ttk::MPIrank_) {
            cellIdentifiers_[i] = cellIndex;
            cellIndex++;
          }
          if(outdatedGlobalCellIds_ != nullptr) {
            cellOutdatedGtoL_[outdatedGlobalCellIds_[i]] = i;
          }
        }
      } else {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellGhost_[i] == 0) {
            cellIdentifiers_[i] = cellIndex;
            cellIndex++;
          }
          if(outdatedGlobalCellIds_ != nullptr) {
            cellOutdatedGtoL_[outdatedGlobalCellIds_[i]] = i;
          }
        }
      }
    }

    /**
     * @brief Generates global ids for the PolyData data set type
     *
     * @return int: 1 for success
     */

    int executePolyData() {
      vertGtoL_->clear();
      std::vector<Point> vertGhostCoordinates;
      std::vector<ttk::SimplexId> vertGhostGlobalIds;
      std::vector<ttk::SimplexId> cellGhostGlobalVertexIds;
      std::vector<ttk::SimplexId> cellGhostGlobalIds;
      std::vector<std::vector<Point>> vertGhostCoordinatesPerRank;
      std::vector<std::vector<ttk::SimplexId>> vertGhostGlobalIdsPerRank;
      std::vector<std::vector<ttk::SimplexId>> cellGhostGlobalVertexIdsPerRank;
      std::vector<std::vector<ttk::SimplexId>> cellGhostGlobalIdsPerRank;
      // In case RankArray is defined for vertices, a vector needs to be resized
      // In case outdated global point ids are defined, then it is
      // vertGhostGlobalIdsPerRank that need to be resized. This vector will
      // store the outdated global ids of ghost points and their local id.
      // Otherwise, vertGhostCoordinatesPerRank needs to be resized. This vector
      // will store the coordinates of a ghost point and its local id.
      if(vertexRankArray_ != nullptr) {
        if(outdatedGlobalPointIds_ == nullptr) {
          vertGhostCoordinatesPerRank.resize(neighborNumber_);
        } else {
          vertGhostGlobalIdsPerRank.resize(neighborNumber_);
        }
      }

      // Similar as above, but with cells
      if(cellRankArray_ != nullptr) {
        if(outdatedGlobalCellIds_ == nullptr) {
          cellGhostGlobalVertexIdsPerRank.resize(neighborNumber_);
        } else {
          cellGhostGlobalIdsPerRank.resize(neighborNumber_);
        }
      }

      this->generateGlobalIds(vertGhostCoordinatesPerRank, vertGhostCoordinates,
                              vertGhostGlobalIdsPerRank, vertGhostGlobalIds);

      // Start exchange of data to identify the global ids of ghost simplices
      // records whether the process has sent its ghost simplices
      hasSentData_ = 0;
      ttk::SimplexId recvMessageSize = 0;
      std::vector<Point> receivedPoints;
      std::vector<ttk::SimplexId> receivedIds;
      std::vector<Response> receivedResponse;
      std::vector<Response> send_buf;
      std::vector<std::vector<Response>> locatedSimplices(threadNumber_);
      // For each neighbor, a process will receive the ghost points of all its
      // neighbors or, if it is its turn, will receive the ghost points of all
      // its neighbors
      for(int i = 0; i < neighborNumber_ + 1; i++) {
        if((i == neighborNumber_ && !hasSentData_)
           || (!hasSentData_ && ttk::MPIrank_ < neighbors_->at(i))) {
          // it is the turn of the current process to send its ghosts
          if(outdatedGlobalPointIds_ == nullptr) {
            if(vertexRankArray_ == nullptr) {
              sendToAllNeighbors<Point>(vertGhostCoordinates, mpiPointType_);
            } else {
              sendToAllNeighborsUsingRankArray<Point>(
                vertGhostCoordinatesPerRank, mpiPointType_);
            }
          } else {
            if(vertexRankArray_ == nullptr) {
              sendToAllNeighbors<ttk::SimplexId>(
                vertGhostGlobalIds, mpiIdType_);
            } else {
              sendToAllNeighborsUsingRankArray<ttk::SimplexId>(
                vertGhostGlobalIdsPerRank, mpiIdType_);
            }
          }
          // The current process receives the responses of all the processes
          // and associates its ghosts with the right global identifier
          for(int j = 0; j < neighborNumber_; j++) {
            recvVector<Response>(receivedResponse, recvMessageSize,
                                 mpiResponseType_, neighbors_->at(j));
            if(outdatedGlobalPointIds_ == nullptr) {
              for(int n = 0; n < recvMessageSize; n++) {
                vertexIdentifiers_[receivedResponse[n].id]
                  = receivedResponse[n].globalId;
                (*vertGtoL_)[receivedResponse[n].globalId]
                  = receivedResponse[n].id;
              }
            } else {
              for(int n = 0; n < recvMessageSize; n++) {
                vertexIdentifiers_[vertOutdatedGtoL_[receivedResponse[n].id]]
                  = receivedResponse[n].globalId;
                (*vertGtoL_)[receivedResponse[n].globalId]
                  = receivedResponse[n].id;
              }
            }
          }
          hasSentData_ = 1;
        } else {
          // It is not the turn of the current process to send its data.
          if(outdatedGlobalPointIds_ == nullptr) {
            recvVector<Point>(receivedPoints, recvMessageSize, mpiPointType_,
                              neighbors_->at(i - hasSentData_));
            // Point coordinates are matched to a local point using a kd-tree
            // and its global id is added to the vector send_buf
            locatePoints(
              locatedSimplices, receivedPoints, recvMessageSize, send_buf);
          } else {
            recvVector<ttk::SimplexId>(receivedIds, recvMessageSize, mpiIdType_,
                                       neighbors_->at(i - hasSentData_));
            // Outdated global ids are matched to a local point
            // and its global id is added to the vector send_buf
            identifyPoints(
              locatedSimplices, receivedIds, recvMessageSize, send_buf);
          }
          // The points founds are sent back to the neighbor with their global
          // ids
          sendVector<Response>(
            send_buf, mpiResponseType_, neighbors_->at(i - hasSentData_));
        }
      }
      // Start of computation of ghost information for cells
      int id{-1};
      // If outdated global ids exist, for each ghost cell is added its local id
      // and its outdated global id
      if(cellRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellRankArray_[i] != ttk::MPIrank_) {
            if(outdatedGlobalCellIds_ == nullptr) {
              cellGhostGlobalVertexIdsPerRank[neighborToId_[cellRankArray_[i]]]
                .push_back(i);
              for(int k = 0; k < dimension_ + 1; k++) {
                id = connectivity_[i * (dimension_ + 1) + k];
                cellGhostGlobalVertexIdsPerRank
                  [neighborToId_[cellRankArray_[i]]]
                    .push_back(vertexIdentifiers_[id]);
              }
            } else {
              cellGhostGlobalIdsPerRank[neighborToId_[cellRankArray_[i]]]
                .push_back(
                  static_cast<ttk::SimplexId>(outdatedGlobalCellIds_[i]));
            }
          }
        }
      } else {
        // If no outdated global ids exist, for each ghost cell is added its
        // local id and the global id of all of its vertices
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellGhost_[i] != 0) {
            if(outdatedGlobalCellIds_ == nullptr) {
              cellGhostGlobalVertexIds.push_back(i);
              for(int k = 0; k < dimension_ + 1; k++) {
                id = connectivity_[i * (dimension_ + 1) + k];
                cellGhostGlobalVertexIds.push_back(vertexIdentifiers_[id]);
              }
            } else {
              cellGhostGlobalIds.push_back(
                static_cast<ttk::SimplexId>(outdatedGlobalCellIds_[i]));
            }
          }
        }
      }
      hasSentData_ = 0;

      // Exchange cells similarly to what is done for vertices
      for(int i = 0; i < neighborNumber_ + 1; i++) {
        if((i == neighborNumber_ && !hasSentData_)
           || (!hasSentData_ && ttk::MPIrank_ < neighbors_->at(i))) {
          // For each neighbor, a process will receive the ghost cells of all
          // its neighbors or, if it is its turn, will receive the ghost cells
          // of all its neighbors
          if(outdatedGlobalCellIds_ == nullptr) {
            if(cellRankArray_ == nullptr) {
              sendToAllNeighbors<ttk::SimplexId>(
                cellGhostGlobalVertexIds, mpiIdType_);
            } else {
              sendToAllNeighborsUsingRankArray<ttk::SimplexId>(
                cellGhostGlobalVertexIdsPerRank, mpiIdType_);
            }
          } else {
            if(cellRankArray_ == nullptr) {
              sendToAllNeighbors<ttk::SimplexId>(
                cellGhostGlobalIds, mpiIdType_);
            } else {
              sendToAllNeighborsUsingRankArray<ttk::SimplexId>(
                cellGhostGlobalIdsPerRank, mpiIdType_);
            }
          }
          // The current process receives the responses of all the processes
          // and associates its ghosts with the right global identifier
          for(int j = 0; j < neighborNumber_; j++) {
            recvVector<Response>(receivedResponse, recvMessageSize,
                                 mpiResponseType_, neighbors_->at(j));
            if(outdatedGlobalCellIds_ == nullptr) {
              for(int n = 0; n < recvMessageSize; n++) {
                cellIdentifiers_[receivedResponse[n].id]
                  = receivedResponse[n].globalId;
              }
            } else {
              for(int n = 0; n < recvMessageSize; n++) {
                cellIdentifiers_[cellOutdatedGtoL_[receivedResponse[n].id]]
                  = receivedResponse[n].globalId;
              }
            }
          }
          hasSentData_ = 1;
        } else {
          // It is not the turn of the current process to send its data.
          recvVector<ttk::SimplexId>(receivedIds, recvMessageSize, mpiIdType_,
                                     neighbors_->at(i - hasSentData_));
          if(outdatedGlobalCellIds_ == nullptr) {
            // Point global ids are matched to a local cell and the global
            // id of the cell is added to the vector send_buf
            locateCells(
              locatedSimplices, receivedIds, recvMessageSize, send_buf);
          } else {
            // Outdated global ids are matched to a local cell
            // and its global id is added to the vector send_buf
            identifyCells(
              locatedSimplices, receivedIds, recvMessageSize, send_buf);
          }
          // The cells founds are sent back to the neighbor with their global
          // ids
          sendVector<Response>(
            send_buf, mpiResponseType_, neighbors_->at(i - hasSentData_));
        }
      }

      return 1; // return success
    }

    /**
     * @brief Generates global ids for the ImageData data set type
     *
     * @return int: 1 for success
     */
    int executeImageData() {

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      // Reorganize bounds to only execute Allreduce twice
      double tempBounds[6] = {
        bounds_[0], bounds_[2], bounds_[4], bounds_[1], bounds_[3], bounds_[5]};
      double tempGlobalBounds[6];
      // Compute and send to all processes the lower bounds of the data set
      MPI_Allreduce(
        tempBounds, tempGlobalBounds, 3, MPI_DOUBLE, MPI_MIN, ttk::MPIcomm_);

      // Compute and send to all processes the higher bounds of the data set
      MPI_Allreduce(tempBounds + 3, tempGlobalBounds + 3, 3, MPI_DOUBLE,
                    MPI_MAX, ttk::MPIcomm_);
      // Global bounds
      double globalBounds[6]
        = {tempGlobalBounds[0], tempGlobalBounds[3], tempGlobalBounds[1],
           tempGlobalBounds[4], tempGlobalBounds[2], tempGlobalBounds[5]};
      // Compute global width and height of the data set
      int width
        = static_cast<int>((globalBounds[1] - globalBounds[0]) / spacing_[0])
          + 1;
      int height
        = static_cast<int>((globalBounds[3] - globalBounds[2]) / spacing_[1])
          + 1;

      // Compute offset of the current process for each direction
      int offsetWidth = static_cast<int>(
        std::round((bounds_[0] - globalBounds[0]) / spacing_[0]));
      int offsetHeight = static_cast<int>(
        std::round((bounds_[2] - globalBounds[2]) / spacing_[1]));
      int offsetLength = static_cast<int>(
        std::round((bounds_[4] - globalBounds[4]) / spacing_[2]));
      // Generate global ids for vertices
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(int k = 0; k < dims_[2]; k++) {
        for(int j = 0; j < dims_[1]; j++) {
          for(int i = 0; i < dims_[0]; i++) {
            vertexIdentifiers_[k * dims_[0] * dims_[1] + j * dims_[0] + i]
              = i + offsetWidth + (j + offsetHeight) * width
                + (k + offsetLength) * width * height;
          }
        }
      }

      // Generate global ids for cells
      dims_[0] -= 1;
      dims_[1] -= 1;
      dims_[2] -= 1;
      width -= 1;
      height -= 1;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(int k = 0; k < dims_[2]; k++) {
        for(int j = 0; j < dims_[1]; j++) {
          for(int i = 0; i < dims_[0]; i++) {
            cellIdentifiers_[k * dims_[0] * dims_[1] + j * dims_[0] + i]
              = i + offsetWidth + (j + offsetHeight) * width
                + (k + offsetLength) * width * height;
          }
        }
      }

      return 1; // return success
    }

#endif

    /**
     * @brief Generates global ids for all data set type in sequential
     *
     * @return int
     */
    int executeSequential() {
      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(SimplexId i = 0; i < vertexNumber_; i++) {
        // avoid any processing if the abort signal is sent
        vertexIdentifiers_[i] = i;
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(SimplexId i = 0; i < cellNumber_; i++) {
        // avoid any processing if the abort signal is sent
        cellIdentifiers_[i] = i;
      }

      return 1; // return success
    }

  }; // Identifiers class

} // namespace ttk
