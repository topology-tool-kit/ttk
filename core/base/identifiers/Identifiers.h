/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::Identifiers
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %Identifiers class that computes for each vertex of
/// a triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'Identifiers'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Geometry.h>
#include <Triangulation.h>
#include <map>

namespace ttk {

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

  /**
   * The Identifiers class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class Identifiers : virtual public Debug {

  public:
    Identifiers();

    int setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    int setCellNumber(const SimplexId &cellNumber) {
      cellNumber_ = cellNumber;
      return 0;
    }

    inline void setDomainDimension(const int &dimension) {
      dimension_ = dimension;
    }

  protected:
    int nbPoints_{0};
    std::map<ttk::SimplexId, ttk::SimplexId> vertGtoL_;
    std::vector<int> neighbors_;
    std::map<int, int> neighborToId_;
    int neighborNumber_;
    double *bounds_;
    ttk::SimplexId vertexNumber_{};
    ttk::SimplexId cellNumber_{};
    MPI_Datatype mpiIdType_;
    MPI_Datatype mpiResponseType_;
    MPI_Datatype mpiPointType_;
    int dimension_{};
    ttk::SimplexId *vertRankArray_{nullptr};
    ttk::SimplexId *cellRankArray_{nullptr};
    unsigned char *vertGhost_{nullptr};
    unsigned char *cellGhost_{nullptr};
    std::vector<ttk::SimplexId> *vertexIdentifiers_;
    std::vector<ttk::SimplexId> *cellIdentifiers_;
    std::vector<std::vector<ttk::SimplexId>> pointsToCells_;

    void initializeNeighbors(double *boundingBox) {
      getNeighborsUsingBoundingBox(boundingBox, neighbors_);
      neighborNumber_ = neighbors_.size();
      for(int i = 0; i < neighborNumber_; i++) {
        neighborToId_[neighbors_[i]] = i;
      }
    }
    template <class triangulationType>
    void initializePointsToCells(triangulationType *triangulation) {
      pointsToCells_.resize(vertexNumber_);
      int id{-1};
      for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
        int starNumber = triangulation->getVertexStarNumber(i);
        for(int j = 0; j < starNumber; j++) {
          triangulation->getVertexStar(i, j, id);
          if(cellGhost_[id] == 0) {
            pointsToCells_[i].push_back(id);
          }
        }
      }
    }

  public:
    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(ttk::AbstractTriangulation *triangulation) {
      setDomainDimension(triangulation->getDimensionality());
      setVertexNumber(triangulation->getNumberOfVertices());
      return triangulation->preconditionVertexStars();
      ;
    }

    /**
     * TODO 3: Implmentation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */
    template <class triangulationType = ttk::AbstractTriangulation>
    int execute(triangulationType *triangulation) {
      // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      std::vector<Point> vertGhostCoordinates;
      std::vector<ttk::SimplexId> cellGhostGlobalVertexIds;
      std::vector<ttk::SimplexId> cellGhostLocalIds;
      std::vector<std::vector<Point>> vertGhostCoordinatesPerRank;
      std::vector<std::vector<ttk::SimplexId>> cellGhostLocalIdsPerRank;
      std::vector<std::vector<ttk::SimplexId>> cellGhostGlobalVertexIdsPerRank;

      if(vertRankArray_ != nullptr) {
        vertGhostCoordinatesPerRank.resize(neighborNumber_);
      }
      if(cellRankArray_ != nullptr) {
        cellGhostGlobalVertexIdsPerRank.resize(neighborNumber_);
        cellGhostLocalIdsPerRank.resize(neighborNumber_);
      }

      this->generateGlobalIds<triangulationType>(
        triangulation, vertGhostCoordinatesPerRank, vertGhostCoordinates);

      ttk::SimplexId recvMessageSize{0};
      std::vector<Point> receivedPoints;
      std::vector<ttk::SimplexId> receivedCells;
      std::vector<Response> receivedResponse;
      std::vector<Response> locatedSimplices;

      for(int i = 0; i < neighborNumber_; i++) {
        if(vertRankArray_ == nullptr) {
          this->exchangeAndLocatePoints<triangulationType>(
            locatedSimplices, vertGhostCoordinates, receivedPoints,
            receivedResponse, neighbors_[i], recvMessageSize, triangulation);
          int count = 0;
          for(int n = 0; n < vertGhostCoordinates.size(); n++) {
            if(vertGhostCoordinates[n - count].localId
               == receivedResponse[count].id) {
              vertGhostCoordinates.erase(vertGhostCoordinates.begin() + n
                                         - count);
              count++;
            }
          }
        } else {
          this->exchangeAndLocatePoints<triangulationType>(
            locatedSimplices, vertGhostCoordinatesPerRank[i], receivedPoints,
            receivedResponse, neighbors_[i], recvMessageSize, triangulation);
        }
      }
      printMsg("POINTS DONE, START CELLS");

      initializePointsToCells<triangulationType>(triangulation);

      int id{-1};

      if(cellRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellRankArray_[i] != ttk::MPIrank_) {
            int vertexNumber = triangulation->getCellVertexNumber(i);
            cellGhostLocalIdsPerRank[neighborToId_[cellRankArray_[i]]]
              .push_back(i);
            for(int k = 0; k < vertexNumber; k++) {
              triangulation->getCellVertex(i, k, id);
              cellGhostGlobalVertexIdsPerRank[neighborToId_[cellRankArray_[i]]]
                .push_back(vertexIdentifiers_->at(id));
            }
          }
        }
      } else {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellGhost_[i] != 0) {
            int vertexNumber = triangulation->getCellVertexNumber(i);
            cellGhostLocalIds.push_back(i);
            for(int k = 0; k < vertexNumber; k++) {
              triangulation->getCellVertex(i, k, id);
              cellGhostGlobalVertexIds.push_back(vertexIdentifiers_->at(id));
            }
          }
        }
      }
      MPI_Barrier(ttk::MPIcomm_);
      printMsg("Start communication phase");

      for(int i = 0; i < neighborNumber_; i++) {
        if(cellRankArray_ == nullptr) {
          this->exchangeAndLocateCells<triangulationType>(
            locatedSimplices, cellGhostGlobalVertexIds, cellGhostLocalIds,
            receivedCells, receivedResponse, neighbors_[i], recvMessageSize,
            triangulation);

          for(int n = 0; n < recvMessageSize; n++) {
            cellIdentifiers_->at(
              cellGhostLocalIds[receivedResponse[n].id / nbPoints_ - n])
              = receivedResponse[n].globalId;
            cellGhostLocalIds.erase(cellGhostLocalIds.begin()
                                    + receivedResponse[n].id / nbPoints_ - n);
            cellGhostGlobalVertexIds.erase(
              cellGhostGlobalVertexIds.begin() + receivedResponse[n].id
                - n * nbPoints_,
              cellGhostGlobalVertexIds.begin() + receivedResponse[n].id
                - n * nbPoints_ + nbPoints_);
          }
        } else {
          this->exchangeAndLocateCells<triangulationType>(
            locatedSimplices, cellGhostGlobalVertexIdsPerRank[i],
            cellGhostLocalIdsPerRank[i], receivedCells, receivedResponse,
            neighbors_[i], recvMessageSize, triangulation);
          for(int n = 0; n < recvMessageSize; n++) {
            cellIdentifiers_->at(
              cellGhostLocalIdsPerRank[i][receivedResponse[n].id / nbPoints_])
              = receivedResponse[n].globalId;
          }
        }
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

    void setVertexIdentifiers(std::vector<ttk::SimplexId> *vertexIdentifiers) {
      this->vertexIdentifiers_ = vertexIdentifiers;
    }

    void setCellIdentifiers(std::vector<ttk::SimplexId> *cellIdentifiers) {
      this->cellIdentifiers_ = cellIdentifiers;
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

    template <class triangulationType>
    void inline findPoint(ttk::SimplexId &id,
                          float x,
                          float y,
                          float z,
                          triangulationType *triangulation) {
      float pointToFind[3] = {x, y, z};
      float dist{0};
      float p[3];
      triangulation->getVertexPoint(0, p[0], p[1], p[2]);
      float minDist = Geometry::distance(pointToFind, p);
      ttk::SimplexId indexMin = 0;
      for(int i = 1; i < vertexNumber_; i++) {
        triangulation->getVertexPoint(i, p[0], p[1], p[2]);
        dist = Geometry::distance(pointToFind, p, dimension_);
        if(dist < minDist) {
          minDist = dist;
          indexMin = i;
        }
      }
      id = indexMin;
    }

    template <typename triangulationType>
    void inline exchangeAndLocatePoints(
      std::vector<Response> &locatedSimplices,
      std::vector<Point> &simplicesCoordinates,
      std::vector<Point> &receivedPoints,
      std::vector<Response> &receivedResponse,
      int neighbor,
      int &recvMessageSize,
      triangulationType *triangulation) {
      locatedSimplices.clear();
      this->SendRecvVector<Point>(simplicesCoordinates, receivedPoints,
                                  recvMessageSize, mpiPointType_, neighbor);
      ttk::SimplexId globalId{-1};
      ttk::SimplexId id{-1};
      for(int n = 0; n < recvMessageSize; n++) {
        if(bounds_[0] <= receivedPoints[n].x
           && bounds_[1] >= receivedPoints[n].x
           && bounds_[2] <= receivedPoints[n].y
           && bounds_[3] >= receivedPoints[n].y
           && bounds_[4] <= receivedPoints[n].z
           && bounds_[5] >= receivedPoints[n].z) {
          this->findPoint<triangulationType>(
            id, receivedPoints[n].x, receivedPoints[n].y, receivedPoints[n].z,
            triangulation);
          if(vertGhost_[id] == 0) {
            globalId = vertexIdentifiers_->at(id);
            if(globalId >= 0) {
              locatedSimplices.push_back(
                Response{receivedPoints[n].localId, globalId});
            }
          }
        }
      }
      this->SendRecvVector<Response>(locatedSimplices, receivedResponse,
                                     recvMessageSize, mpiResponseType_,
                                     neighbor);
      for(int n = 0; n < recvMessageSize; n++) {
        vertexIdentifiers_->at(receivedResponse[n].id)
          = receivedResponse[n].globalId;
        vertGtoL_[receivedResponse[n].globalId] = receivedResponse[n].id;
      }
    }

    template <typename triangulationType>
    void inline exchangeAndLocateCells(
      std::vector<Response> &locatedSimplices,
      std::vector<ttk::SimplexId> &cellGhostGlobalVertexIds,
      std::vector<ttk::SimplexId> &cellGhostLocalIds,
      std::vector<ttk::SimplexId> &receivedCells,
      std::vector<Response> &receivedResponse,
      int neighbor,
      int &recvMessageSize,
      triangulationType *triangulation) {
      int id{-1};
      std::map<ttk::SimplexId, ttk::SimplexId>::iterator search;
      std::vector<ttk::SimplexId> localPointIds;
      localPointIds.reserve(nbPoints_);
      locatedSimplices.clear();
      this->SendRecvVector<ttk::SimplexId>(cellGhostGlobalVertexIds,
                                           receivedCells, recvMessageSize,
                                           mpiIdType_, neighbor);
      for(int n = 0; n < recvMessageSize; n += nbPoints_) {
        localPointIds.clear();
        for(int k = 0; k < nbPoints_; k++) {
          search = vertGtoL_.find(receivedCells[n + k]);
          if(search != vertGtoL_.end()) {
            localPointIds.push_back(search->second);
          } else {
            break;
          }
        }
        if(localPointIds.size() != static_cast<size_t>(nbPoints_)) {
          break;
        }

        bool foundIt = false;
        int k = 0;
        int l;
        int m = 0;
        while(!foundIt && m < nbPoints_) {
          int size = pointsToCells_[localPointIds[m]].size();
          k = 0;
          while(!foundIt && k < size) {
            int vertexNumber = triangulation->getCellVertexNumber(
              pointsToCells_[localPointIds[m]][k]);
            l = 0;
            while(l < vertexNumber) {
              triangulation->getCellVertex(
                pointsToCells_[localPointIds[m]][k], l, id);
              auto it = find(localPointIds.begin(), localPointIds.end(), id);
              if(it == localPointIds.end()) {
                break;
              }
              l++;
            }
            if(l == vertexNumber) {
              foundIt = true;
              locatedSimplices.push_back(Response{
                n, cellIdentifiers_->at(pointsToCells_[localPointIds[m]][k])});
            }
            k++;
          }
          m++;
        }
      }
      if(recvMessageSize > 0 && cellGhostGlobalVertexIds.size() > 0)
        this->SendRecvVector<Response>(locatedSimplices, receivedResponse,
                                       recvMessageSize, mpiResponseType_,
                                       neighbor);
    }

    template <typename dataType>
    void SendRecvVector(std::vector<dataType> &vectorToSend,
                        std::vector<dataType> &receiveBuffer,
                        int &recvMessageSize,
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

    template <typename triangulationType>
    void generateGlobalIds(
      triangulationType *triangulation,
      std::vector<std::vector<Point>> &vertGhostCoordinatesPerRank,
      std::vector<Point> &vertGhostCoordinates) {
      ttk::SimplexId realVertexNumber = vertexNumber_;
      ttk::SimplexId realCellNumber = cellNumber_;
      float p[3];
      if(vertRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
          if(vertRankArray_[i] != ttk::MPIrank_) {
            realVertexNumber--;
            triangulation->getVertexPoint(i, p[0], p[1], p[2]);
            vertGhostCoordinatesPerRank[neighborToId_[vertRankArray_[i]]]
              .push_back(Point{p[0], p[1], p[2], i});
          }
        }
      } else {
        for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
          if(vertGhost_[i] != 0) {
            realVertexNumber--;
            triangulation->getVertexPoint(i, p[0], p[1], p[2]);
            vertGhostCoordinates.push_back(Point{p[0], p[1], p[2], i});
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
        }
      } else {
        for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
          if(vertGhost_[i] == 0) {
            vertexIdentifiers_->at(i) = vertIndex;
            vertGtoL_[vertIndex] = i;

            vertIndex++;
          }
        }
      }
      if(cellRankArray_ != nullptr) {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellRankArray_[i] == ttk::MPIrank_) {
            cellIdentifiers_->at(i) = cellIndex;
            cellIndex++;
          }
        }
      } else {
        for(ttk::SimplexId i = 0; i < cellNumber_; i++) {
          if(cellGhost_[i] == 0) {
            cellIdentifiers_->at(i) = cellIndex;
            cellIndex++;
          }
        }
      }
    }

  }; // Identifiers class

} // namespace ttk
