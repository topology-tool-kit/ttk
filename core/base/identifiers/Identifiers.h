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
    ttk::SimplexId globalId;
  };

#endif

  /**
   * The Identifiers class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class Identifiers : virtual public Debug {

  public:
    Identifiers();

    void setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
    }

    void setCellNumber(const SimplexId &cellNumber) {
      cellNumber_ = cellNumber;
    }

    void setVertexIdentifiers(std::vector<ttk::SimplexId> *vertexIdentifiers) {
      this->vertexIdentifiers_ = vertexIdentifiers;
    }

    void setCellIdentifiers(std::vector<ttk::SimplexId> *cellIdentifiers) {
      this->cellIdentifiers_ = cellIdentifiers;
    }

#ifdef TTK_ENABLE_MPI
    inline void setDomainDimension(const int &dimension) {
      dimension_ = dimension;
    }

    inline void setPointsToCells(
      const std::vector<std::vector<ttk::SimplexId>> pointsToCells) {
      pointsToCells_ = pointsToCells;
    }

    void setVertRankArray(ttk::SimplexId *vertRankArray) {
      this->vertRankArray_ = vertRankArray;
    }

    void setCellRankArray(ttk::SimplexId *cellRankArray) {
      this->cellRankArray_ = cellRankArray;
    }

    void setVertGhost(unsigned char *vertGhost) {
      this->vertGhost_ = vertGhost;
    }

    void setCellGhost(unsigned char *cellGhost) {
      this->cellGhost_ = cellGhost;
    }

    void setNbPoints(int nbPoints) {
      this->nbPoints_ = nbPoints;
    }

    void setBounds(double *bounds) {
      this->bounds_ = bounds;
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
      kdt_ = KDTree<float>(false, 2);
      kdt_.build(pointSet_, vertexNumber_, dimension_);
    }

    void setOutdatedGlobalCellIds(ttk::LongSimplexId *outdatedGlobalCellIds) {
      outdatedGlobalCellIds_ = outdatedGlobalCellIds;
    }

    void initializeMPITypes() {
      ttk::SimplexId id{-1};
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
      MPI_Datatype typesResponse[] = {mpiIdType_, mpiIdType_};
      int lengthsResponse[] = {1, 1};
      const long int mpi_offsetsResponse[]
        = {offsetof(Response, id), offsetof(Response, globalId)};
      MPI_Type_create_struct(2, lengthsResponse, mpi_offsetsResponse,
                             typesResponse, &mpiResponseType_);
      MPI_Type_commit(&mpiResponseType_);
    }

    void inline findPoint(ttk::SimplexId &id, float x, float y, float z) {
      std::vector<float> coordinates = {x, y, z};
      std::vector<KDTree<float> *> neighbours;
      std::vector<float> costs;
      kdt_.getKClosest(1, coordinates, neighbours, costs);
      id = neighbours[0]->id_;
    }

    void initializeNeighbors(double *boundingBox) {
      getNeighborsUsingBoundingBox(boundingBox, neighbors_);
      neighborNumber_ = neighbors_.size();
      for(int i = 0; i < neighborNumber_; i++) {
        neighborToId_[neighbors_[i]] = i;
      }
    }

    void inline locatePoints(
      std::vector<std::vector<Response>> &locatedSimplices,
      std::vector<Point> &receivedPoints,
      ttk::SimplexId &recvMessageSize,
      std::vector<Response> &send_buf) {
      ttk::SimplexId globalId{-1};
      ttk::SimplexId id{-1};
#pragma omp parallel for num_threads(threadNumber_) firstprivate(globalId, id) \
  shared(locatedSimplices)
      for(int n = 0; n < recvMessageSize; n++) {
        if(bounds_[0] <= receivedPoints[n].x
           && bounds_[1] >= receivedPoints[n].x
           && bounds_[2] <= receivedPoints[n].y
           && bounds_[3] >= receivedPoints[n].y
           && bounds_[4] <= receivedPoints[n].z
           && bounds_[5] >= receivedPoints[n].z) {
          this->findPoint(
            id, receivedPoints[n].x, receivedPoints[n].y, receivedPoints[n].z);
          if((vertGhost_ != nullptr && vertGhost_[id] == 0)
             || (vertRankArray_ != nullptr
                 && vertRankArray_[id] == ttk::MPIrank_)) {
            globalId = vertexIdentifiers_->at(id);
            if(globalId >= 0) {
              locatedSimplices[omp_get_thread_num()].push_back(
                Response{receivedPoints[n].localId, globalId});
            }
          }
        }
      }
      send_buf.clear();
      for(int n = 0; n < threadNumber_; n++) {
        send_buf.insert(send_buf.end(), locatedSimplices[n].begin(),
                        locatedSimplices[n].end());
        locatedSimplices[n].clear();
      }
    }

    void inline identifyPoints(
      std::vector<std::vector<Response>> &locatedSimplices,
      std::vector<ttk::SimplexId> &receivedOutdatedGlobalIds,
      ttk::SimplexId &recvMessageSize,
      std::vector<Response> &send_buf) {
      ttk::SimplexId globalId{-1};
      std::map<ttk::LongSimplexId, ttk::SimplexId>::iterator search;
#pragma omp parallel for num_threads(threadNumber_)
      for(int n = 0; n < recvMessageSize; n++) {
        search = vertOutdatedGtoL_.find(receivedOutdatedGlobalIds[n]);
        if(search != vertOutdatedGtoL_.end()) {
          globalId = vertexIdentifiers_->at(search->second);
          if(globalId >= 0) {
            locatedSimplices[omp_get_thread_num()].push_back(
              Response{receivedOutdatedGlobalIds[n], globalId});
          }
        }
      }
      send_buf.clear();
      for(int n = 0; n < threadNumber_; n++) {
        send_buf.insert(send_buf.end(), locatedSimplices[n].begin(),
                        locatedSimplices[n].end());
        locatedSimplices[n].clear();
      }
    }

    void inline locateCells(
      std::vector<std::vector<Response>> &locatedSimplices,
      std::vector<ttk::SimplexId> &receivedCells,
      ttk::SimplexId &recvMessageSize,
      std::vector<Response> &send_buf) {
      send_buf.clear();
      int id{-1};
      std::map<ttk::SimplexId, ttk::SimplexId>::iterator search;
      std::vector<ttk::SimplexId> localPointIds;
      localPointIds.reserve(nbPoints_);
#pragma omp parallel for num_threads(threadNumber_) firstprivate(localPointIds)
      for(int n = 0; n < recvMessageSize; n += nbPoints_ + 1) {
        localPointIds.clear();
        for(int k = 1; k < nbPoints_ + 1; k++) {
          search = vertGtoL_.find(receivedCells[n + k]);
          if(search != vertGtoL_.end()) {
            localPointIds.push_back(search->second);
          } else {
            break;
          }
        }
        if(localPointIds.size() == static_cast<size_t>(nbPoints_)) {
          bool foundIt = false;
          int k = 0;
          int l;
          int m = 0;
          while(!foundIt && m < nbPoints_) {
            int size = pointsToCells_[localPointIds[m]].size();
            k = 0;
            while(!foundIt && k < size) {
              l = 0;
              while(l < nbPoints_) {
                id = connectivity_[pointsToCells_[localPointIds[m]][k]
                                     * nbPoints_
                                   + l];
                auto it = find(localPointIds.begin(), localPointIds.end(), id);
                if(it == localPointIds.end()) {
                  break;
                }
                l++;
              }
              if(l == nbPoints_) {
                foundIt = true;
                locatedSimplices[omp_get_thread_num()].push_back(Response{
                  receivedCells[n],
                  cellIdentifiers_->at(pointsToCells_[localPointIds[m]][k])});
              }
              k++;
            }
            m++;
          }
        }
      }
      for(int n = 0; n < threadNumber_; n++) {
        send_buf.insert(send_buf.end(), locatedSimplices[n].begin(),
                        locatedSimplices[n].end());
        locatedSimplices[n].clear();
      }
    }

    void inline identifyCells(
      std::vector<std::vector<Response>> &locatedSimplices,
      std::vector<ttk::SimplexId> &receivedOutdatedGlobalIds,
      ttk::SimplexId &recvMessageSize,
      std::vector<Response> &send_buf) {
      ttk::SimplexId globalId{-1};
      std::map<ttk::LongSimplexId, ttk::SimplexId>::iterator search;
      send_buf.clear();
#pragma omp parallel for num_threads(threadNumber_)
      for(int n = 0; n < recvMessageSize; n++) {
        search = cellOutdatedGtoL_.find(receivedOutdatedGlobalIds[n]);
        if(search != cellOutdatedGtoL_.end()) {
          globalId = cellIdentifiers_->at(search->second);
          if(globalId >= 0) {
            locatedSimplices[omp_get_thread_num()].push_back(
              Response{receivedOutdatedGlobalIds[n],
                       static_cast<ttk::SimplexId>(globalId)});
          }
        }
      }
      for(int n = 0; n < threadNumber_; n++) {
        send_buf.insert(send_buf.end(), locatedSimplices[n].begin(),
                        locatedSimplices[n].end());
        locatedSimplices[n].clear();
      }
    }

    template <typename dataType>
    void SendRecvVector(std::vector<dataType> &vectorToSend,
                        std::vector<dataType> &receiveBuffer,
                        ttk::SimplexId &recvMessageSize,
                        MPI_Datatype messageType,
                        int neighbor) const {
      ttk::SimplexId dataSize = vectorToSend.size();
      receiveBuffer.clear();
      MPI_Sendrecv(&dataSize, 1, getMPIType(dataSize), neighbor, ttk::MPIrank_,
                   &recvMessageSize, 1, getMPIType(dataSize), neighbor,
                   neighbor, ttk::MPIcomm_, MPI_STATUS_IGNORE);
      receiveBuffer.resize(recvMessageSize);
      MPI_Sendrecv(vectorToSend.data(), dataSize, messageType, neighbor,
                   ttk::MPIrank_, receiveBuffer.data(), recvMessageSize,
                   messageType, neighbor, neighbor, ttk::MPIcomm_,
                   MPI_STATUS_IGNORE);
    }

    template <typename dataType>
    void sendVector(std::vector<dataType> &vectorToSend,
                    MPI_Datatype messageType,
                    int neighbor) const {
      ttk::SimplexId dataSize = vectorToSend.size();
      MPI_Send(&dataSize, 1, getMPIType(dataSize), neighbor, ttk::MPIrank_,
               ttk::MPIcomm_);
      MPI_Send(vectorToSend.data(), dataSize, messageType, neighbor,
               ttk::MPIrank_, ttk::MPIcomm_);
    }

    template <typename dataType>
    void sendToAllNeighbors(std::vector<dataType> &vectorToSend,
                            MPI_Datatype messageType) const {

      for(int j = 0; j < neighborNumber_; j++) {
        sendVector<dataType>(vectorToSend, messageType, neighbors_[j]);
      }
    }

    template <typename dataType>
    void sendToAllNeighborsUsingRankArray(
      std::vector<std::vector<dataType>> &vectorToSend,
      MPI_Datatype messageType) {

      for(int j = 0; j < neighborNumber_; j++) {
        sendVector<dataType>(vectorToSend[neighborToId_[neighbors_[j]]],
                             messageType, neighbors_[j]);
      }
    }

    template <typename dataType>
    void recvVector(std::vector<dataType> &receiveBuffer,
                    ttk::SimplexId &recvMessageSize,
                    MPI_Datatype messageType,
                    int neighbor) const {
      MPI_Recv(&recvMessageSize, 1, getMPIType(recvMessageSize), neighbor,
               neighbor, ttk::MPIcomm_, MPI_STATUS_IGNORE);
      receiveBuffer.resize(recvMessageSize);
      MPI_Recv(receiveBuffer.data(), recvMessageSize, messageType, neighbor,
               neighbor, ttk::MPIcomm_, MPI_STATUS_IGNORE);
    }

    void generateGlobalIds(
      std::vector<std::vector<Point>> &vertGhostCoordinatesPerRank,
      std::vector<Point> &vertGhostCoordinates,
      std::vector<std::vector<ttk::SimplexId>> &vertGhostGlobalIdsPerRank,
      std::vector<ttk::SimplexId> &vertGhostGlobalIds) {
      ttk::SimplexId realVertexNumber = vertexNumber_;
      ttk::SimplexId realCellNumber = cellNumber_;
      float p[3];
      if(vertRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
          if(vertRankArray_[i] != ttk::MPIrank_) {
            realVertexNumber--;
            if(outdatedGlobalPointIds_ == nullptr) {
              p[0] = pointSet_[i * 3];
              p[1] = pointSet_[i * 3 + 1];
              p[2] = pointSet_[i * 3 + 2];
              vertGhostCoordinatesPerRank[neighborToId_[vertRankArray_[i]]]
                .push_back(Point{p[0], p[1], p[2], i});
            } else {
              vertGhostGlobalIdsPerRank[neighborToId_[vertRankArray_[i]]]
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

      // Perform exclusive prefix sum
      MPI_Exscan(
        &realVertexNumber, &vertIndex, 1, mpiIdType_, MPI_SUM, ttk::MPIcomm_);
      MPI_Exscan(
        &realCellNumber, &cellIndex, 1, mpiIdType_, MPI_SUM, ttk::MPIcomm_);

      if(ttk::MPIrank_ == 0) {
        vertIndex = 0;
        cellIndex = 0;
      }
      if(vertRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
          if(vertRankArray_[i] == ttk::MPIrank_) {
            vertexIdentifiers_->at(i) = vertIndex;
            vertGtoL_[vertIndex] = i;
            vertIndex++;
          }
          if(outdatedGlobalPointIds_ != nullptr) {
            vertOutdatedGtoL_[outdatedGlobalPointIds_[i]] = i;
          }
        }
      } else {
        for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
          if(vertGhost_[i] == 0) {
            vertexIdentifiers_->at(i) = vertIndex;
            vertGtoL_[vertIndex] = i;
            vertIndex++;
          }
          if(outdatedGlobalPointIds_ != nullptr) {
            vertOutdatedGtoL_[outdatedGlobalPointIds_[i]] = i;
          }
        }
      }
      if(cellRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellRankArray_[i] == ttk::MPIrank_) {
            cellIdentifiers_->at(i) = cellIndex;
            cellIndex++;
          }
          if(outdatedGlobalCellIds_ != nullptr) {
            cellOutdatedGtoL_[outdatedGlobalCellIds_[i]] = i;
          }
        }
      } else {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellGhost_[i] == 0) {
            cellIdentifiers_->at(i) = cellIndex;
            cellIndex++;
          }
          if(outdatedGlobalCellIds_ != nullptr) {
            cellOutdatedGtoL_[outdatedGlobalCellIds_[i]] = i;
          }
        }
      }
    }

    int executePolyData() {

      std::vector<Point> vertGhostCoordinates;
      std::vector<ttk::SimplexId> vertGhostGlobalIds;
      std::vector<ttk::SimplexId> cellGhostGlobalVertexIds;
      std::vector<ttk::SimplexId> cellGhostGlobalIds;
      std::vector<std::vector<Point>> vertGhostCoordinatesPerRank;
      std::vector<std::vector<ttk::SimplexId>> vertGhostGlobalIdsPerRank;
      std::vector<std::vector<ttk::SimplexId>> cellGhostGlobalVertexIdsPerRank;
      std::vector<std::vector<ttk::SimplexId>> cellGhostGlobalIdsPerRank;
      if(vertRankArray_ != nullptr) {
        vertGhostCoordinatesPerRank.resize(neighborNumber_);
        vertGhostGlobalIdsPerRank.resize(neighborNumber_);
      }
      if(cellRankArray_ != nullptr) {
        cellGhostGlobalVertexIdsPerRank.resize(neighborNumber_);
        cellGhostGlobalIdsPerRank.resize(neighborNumber_);
      }

      this->generateGlobalIds(vertGhostCoordinatesPerRank, vertGhostCoordinates,
                              vertGhostGlobalIdsPerRank, vertGhostGlobalIds);

      ttk::SimplexId recvMessageSize{0};
      std::vector<Point> receivedPoints;
      std::vector<ttk::SimplexId> receivedIds;
      std::vector<Response> receivedResponse;
      std::vector<Response> send_buf;
      std::vector<std::vector<Response>> locatedSimplices(threadNumber_);
      for(int i = 0; i < neighborNumber_ + 1; i++) {
        if((i == neighborNumber_ && !hasSentData_)
           || (ttk::MPIrank_ < neighbors_[i] && !hasSentData_)) {
          if(outdatedGlobalPointIds_ == nullptr) {
            if(vertRankArray_ == nullptr) {
              sendToAllNeighbors<Point>(vertGhostCoordinates, mpiPointType_);
            } else {
              sendToAllNeighborsUsingRankArray<Point>(
                vertGhostCoordinatesPerRank, mpiPointType_);
            }
          } else {
            if(vertRankArray_ == nullptr) {
              sendToAllNeighbors<ttk::SimplexId>(
                vertGhostGlobalIds, mpiIdType_);
            } else {
              sendToAllNeighborsUsingRankArray<ttk::SimplexId>(
                vertGhostGlobalIdsPerRank, mpiIdType_);
            }
          }
          for(int j = 0; j < neighborNumber_; j++) {
            recvVector<Response>(receivedResponse, recvMessageSize,
                                 mpiResponseType_, neighbors_[j]);
            if(outdatedGlobalPointIds_ == nullptr) {
              for(int n = 0; n < recvMessageSize; n++) {
                vertexIdentifiers_->at(receivedResponse[n].id)
                  = receivedResponse[n].globalId;
                vertGtoL_[receivedResponse[n].globalId]
                  = receivedResponse[n].id;
              }
            } else {
              for(int n = 0; n < recvMessageSize; n++) {
                vertexIdentifiers_->at(
                  vertOutdatedGtoL_[receivedResponse[n].id])
                  = receivedResponse[n].globalId;
                vertGtoL_[receivedResponse[n].globalId]
                  = receivedResponse[n].id;
              }
            }
          }
          hasSentData_ = 1;
        } else {
          if(outdatedGlobalPointIds_ == nullptr) {
            recvVector<Point>(receivedPoints, recvMessageSize, mpiPointType_,
                              neighbors_[i - hasSentData_]);
            locatePoints(
              locatedSimplices, receivedPoints, recvMessageSize, send_buf);
          } else {
            recvVector<ttk::SimplexId>(receivedIds, recvMessageSize, mpiIdType_,
                                       neighbors_[i - hasSentData_]);
            identifyPoints(
              locatedSimplices, receivedIds, recvMessageSize, send_buf);
          }
          sendVector<Response>(
            send_buf, mpiResponseType_, neighbors_[i - hasSentData_]);
        }
      }
      int id{-1};
      if(cellRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellRankArray_[i] != ttk::MPIrank_) {
            if(outdatedGlobalCellIds_ == nullptr) {
              cellGhostGlobalVertexIdsPerRank[neighborToId_[cellRankArray_[i]]]
                .push_back(i);
              for(int k = 0; k < nbPoints_; k++) {
                id = connectivity_[i * nbPoints_ + k];
                cellGhostGlobalVertexIdsPerRank
                  [neighborToId_[cellRankArray_[i]]]
                    .push_back(vertexIdentifiers_->at(id));
              }
            } else {
              cellGhostGlobalIdsPerRank[neighborToId_[cellRankArray_[i]]]
                .push_back(
                  static_cast<ttk::SimplexId>(outdatedGlobalCellIds_[i]));
            }
          }
        }
      } else {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellGhost_[i] != 0) {
            if(outdatedGlobalCellIds_ == nullptr) {
              cellGhostGlobalVertexIds.push_back(i);
              for(int k = 0; k < nbPoints_; k++) {
                id = connectivity_[i * nbPoints_ + k];
                cellGhostGlobalVertexIds.push_back(vertexIdentifiers_->at(id));
              }
            } else {
              cellGhostGlobalIds.push_back(
                static_cast<ttk::SimplexId>(outdatedGlobalCellIds_[i]));
            }
          }
        }
      }
      hasSentData_ = 0;

      // Exchange cells
      for(int i = 0; i < neighborNumber_ + 1; i++) {
        if((i == neighborNumber_ && !hasSentData_)
           || (ttk::MPIrank_ < neighbors_[i] && !hasSentData_)) {
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
          for(int j = 0; j < neighborNumber_; j++) {
            recvVector<Response>(receivedResponse, recvMessageSize,
                                 mpiResponseType_, neighbors_[j]);
            if(outdatedGlobalCellIds_ == nullptr) {
              for(int n = 0; n < recvMessageSize; n++) {
                cellIdentifiers_->at(receivedResponse[n].id)
                  = receivedResponse[n].globalId;
              }
            } else {
              for(int n = 0; n < recvMessageSize; n++) {
                cellIdentifiers_->at(cellOutdatedGtoL_[receivedResponse[n].id])
                  = receivedResponse[n].globalId;
              }
            }
          }
          hasSentData_ = 1;
        } else {
          recvVector<ttk::SimplexId>(receivedIds, recvMessageSize, mpiIdType_,
                                     neighbors_[i - hasSentData_]);
          if(outdatedGlobalCellIds_ == nullptr) {
            locateCells(
              locatedSimplices, receivedIds, recvMessageSize, send_buf);
          } else {
            identifyCells(
              locatedSimplices, receivedIds, recvMessageSize, send_buf);
          }
          sendVector<Response>(
            send_buf, mpiResponseType_, neighbors_[i - hasSentData_]);
        }
      }

      return 1; // return success
    }

    int executeImageData() {

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      double tempBounds[6] = {
        bounds_[0], bounds_[2], bounds_[4], bounds_[1], bounds_[3], bounds_[5]};
      double tempGlobalBounds[6];
      MPI_Allreduce(
        tempBounds, tempGlobalBounds, 3, MPI_DOUBLE, MPI_MIN, ttk::MPIcomm_);
      MPI_Allreduce(tempBounds + 3, tempGlobalBounds + 3, 3, MPI_DOUBLE,
                    MPI_MAX, ttk::MPIcomm_);
      double globalBounds[6]
        = {tempGlobalBounds[0], tempGlobalBounds[3], tempGlobalBounds[1],
           tempGlobalBounds[4], tempGlobalBounds[2], tempGlobalBounds[5]};

      int width
        = static_cast<int>((globalBounds[1] - globalBounds[0]) / spacing_[0])
          + 1;
      int height
        = static_cast<int>((globalBounds[3] - globalBounds[2]) / spacing_[1])
          + 1;

      int offsetWidth
        = static_cast<int>((bounds_[0] - globalBounds[0]) / spacing_[0]);
      int offsetHeight
        = static_cast<int>((bounds_[2] - globalBounds[2]) / spacing_[1]);
      int offsetLength
        = static_cast<int>((bounds_[4] - globalBounds[4]) / spacing_[2]);

      // Generate global ids for vertices
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(int k = 0; k < dims_[2]; k++) {
        for(int j = 0; j < dims_[1]; j++) {
          for(int i = 0; i < dims_[0]; i++) {
            vertexIdentifiers_->at(k * dims_[0] * dims_[1] + j * dims_[0] + i)
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
            cellIdentifiers_->at(k * dims_[0] * dims_[1] + j * dims_[0] + i)
              = i + offsetWidth + (j + offsetHeight) * width
                + (k + offsetLength) * width * height;
          }
        }
      }

      return 1; // return success
    }

#endif

    int execute() {
      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(SimplexId i = 0; i < vertexNumber_; i++) {
        // avoid any processing if the abort signal is sent
        vertexIdentifiers_->at(i) = i;
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(SimplexId i = 0; i < cellNumber_; i++) {
        // avoid any processing if the abort signal is sent
        cellIdentifiers_->at(i) = i;
      }

      return 1; // return success
    }

  protected:
    ttk::SimplexId vertexNumber_{};
    ttk::SimplexId cellNumber_{};
    std::vector<ttk::SimplexId> *vertexIdentifiers_;
    std::vector<ttk::SimplexId> *cellIdentifiers_;
#ifdef TTK_ENABLE_MPI
    int nbPoints_{0};
    std::map<ttk::SimplexId, ttk::SimplexId> vertGtoL_;
    std::vector<int> neighbors_;
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
    ttk::SimplexId *vertRankArray_{nullptr};
    ttk::SimplexId *cellRankArray_{nullptr};
    unsigned char *vertGhost_{nullptr};
    unsigned char *cellGhost_{nullptr};
    std::vector<std::vector<ttk::SimplexId>> pointsToCells_;
    float *pointSet_;
    ttk::LongSimplexId *connectivity_;
    ttk::LongSimplexId *outdatedGlobalPointIds_{nullptr};
    ttk::LongSimplexId *outdatedGlobalCellIds_{nullptr};
    std::map<ttk::LongSimplexId, ttk::SimplexId> vertOutdatedGtoL_;
    std::map<ttk::LongSimplexId, ttk::SimplexId> cellOutdatedGtoL_;
    KDTree<float> kdt_;

#endif
  }; // Identifiers class

} // namespace ttk