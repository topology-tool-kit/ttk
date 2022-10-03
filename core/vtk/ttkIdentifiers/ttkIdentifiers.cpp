#include <ttkIdentifiers.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <map>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIdentifiers);

ttkIdentifiers::ttkIdentifiers() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  setDebugMsgPrefix("Identifiers");

  //   vtkWarningMacro("`TTK Identifiers' is now deprecated. Please use "
  //                   "`Generate Ids' instead.");
}

ttkIdentifiers::~ttkIdentifiers() = default;

int ttkIdentifiers::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkIdentifiers::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

void ttkIdentifiers::createMPIPointType(MPI_Datatype *mpiPointType) {
  ttk::SimplexId id{-1};
  MPI_Datatype types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, getMPIType(id)};
  int lengths[] = {1, 1, 1, 1};
  const long int mpi_offsets[] = {offsetof(Point, x), offsetof(Point, y),
                                  offsetof(Point, z), offsetof(Point, localId)};
  MPI_Type_create_struct(4, lengths, mpi_offsets, types, mpiPointType);
  MPI_Type_commit(mpiPointType);
}

void ttkIdentifiers::createMPIResponseType(MPI_Datatype *mpiResponseType) {
  ttk::SimplexId id = 0;
  MPI_Datatype types[] = {getMPIType(id), getMPIType(id)};
  int lengths[] = {1, 1};
  const long int mpi_offsets[]
    = {offsetof(Response, id), offsetof(Response, globalId)};
  MPI_Type_create_struct(2, lengths, mpi_offsets, types, mpiResponseType);
  MPI_Type_commit(mpiResponseType);
}

template <typename dataType>
void ttkIdentifiers::SendRecvVector(std::vector<dataType> &vectorToSend,
                                    std::vector<dataType> &receiveBuffer,
                                    int &recvMessageSize,
                                    MPI_Datatype &messageType,
                                    int neighbor) {
  ttk::SimplexId dataSize = vectorToSend.size();
  receiveBuffer.clear();
  MPI_Sendrecv(&dataSize, 1, getMPIType(dataSize), neighbor, ttk::MPIrank_,
               &recvMessageSize, 1, getMPIType(dataSize), neighbor, neighbor,
               ttk::MPIcomm_, MPI_STATUS_IGNORE);
  receiveBuffer.resize(recvMessageSize);
  MPI_Sendrecv(vectorToSend.data(), dataSize, messageType, neighbor,
               ttk::MPIrank_, receiveBuffer.data(), recvMessageSize,
               messageType, neighbor, neighbor, ttk::MPIcomm_,
               MPI_STATUS_IGNORE);
}

template <typename triangulationType>
void ttkIdentifiers::exchangeAndLocatePoints(
  std::vector<Response> &locatedSimplices,
  std::vector<Point> &simplicesCoordinates,
  std::vector<Point> &receivedPoints,
  std::vector<Response> &receivedResponse,
  int neighbor,
  MPI_Datatype mpiPointType,
  MPI_Datatype mpiResponseType,
  int recvMessageSize,
  vtkDataSet *input,
  double *bounds,
  vtkIntArray *vertexIdentifiers,
  std::map<ttk::SimplexId, ttk::SimplexId> &vertGtoL,
  triangulationType *triangulation) {
  locatedSimplices.clear();
  this->SendRecvVector<Point>(simplicesCoordinates, receivedPoints,
                              recvMessageSize, mpiPointType, neighbor);
  ttk::SimplexId globalId{-1};
  ttk::SimplexId id{-1};
  for(int n = 0; n < recvMessageSize; n++) {
    if(bounds[0] <= receivedPoints[n].x && bounds[1] >= receivedPoints[n].x
       && bounds[2] <= receivedPoints[n].y && bounds[3] >= receivedPoints[n].y
       && bounds[4] <= receivedPoints[n].z
       && bounds[5] >= receivedPoints[n].z) {
      this->findPoint<triangulationType>(id, receivedPoints[n].x,
                                         receivedPoints[n].y,
                                         receivedPoints[n].z, triangulation);
      printMsg("id: " + std::to_string(id));
      if(id >= 0) {
        globalId = vertexIdentifiers->GetTuple1(id);
        if(globalId >= 0) {
          locatedSimplices.push_back(
            Response{receivedPoints[n].localId, globalId});
        }
      }
    }
  }
  this->SendRecvVector<Response>(locatedSimplices, receivedResponse,
                                 recvMessageSize, mpiResponseType, neighbor);
  for(int n = 0; n < recvMessageSize; n++) {
    vertexIdentifiers->SetTuple1(
      receivedResponse[n].id, receivedResponse[n].globalId);
    vertGtoL[receivedResponse[n].globalId] = receivedResponse[n].id;
  }
}

int ttkIdentifiers::RequestData(vtkInformation *ttkNotUsed(request),
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);
  // vtkImageData* data = vtkImageData::SafeDownCast(input);
  // data->ComputeBounds();
  // double* bounds = data->GetBounds();
  //   cout << "bounds: " << bounds[0] << "," << bounds[1] << "  " << bounds[2]
  //        << "," << bounds[3] << " " << bounds[4] << "," <<
  //        bounds[5] << endl;
  // double* point = data->GetPoint(0);
  // printMsg("Point: "+std::to_string(point[0])+" "+std::to_string(point[1])+"
  // "+std::to_string(point[2])); data->SetExtent(wholeExtent); data->Update();

  Timer t;

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifiers
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vtkSmartPointer<ttkSimplexIdTypeArray> cellIdentifiers
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();

  vertexIdentifiers->SetName(VertexFieldName.data());
  vertexIdentifiers->SetNumberOfComponents(1);
  vertexIdentifiers->SetNumberOfTuples(input->GetNumberOfPoints());
  vertexIdentifiers->Fill(-1);

  cellIdentifiers->SetName(CellFieldName.data());
  cellIdentifiers->SetNumberOfComponents(1);
  cellIdentifiers->SetNumberOfTuples(input->GetNumberOfCells());
  cellIdentifiers->Fill(-1);

  SimplexId vertexNumber = input->GetNumberOfPoints();
  SimplexId cellNumber = input->GetNumberOfCells();

  SimplexId realVertexNumber = vertexNumber;
  SimplexId realCellNumber = cellNumber;

  vtkDataArray *vertRankArray = input->GetPointData()->GetArray("RankArray");
  vtkDataArray *cellRankArray = input->GetCellData()->GetArray("RankArray");

  vtkDataArray *vertGhost = input->GetPointData()->GetArray("vtkGhostType");
  vtkDataArray *cellGhost = input->GetCellData()->GetArray("vtkGhostType");
  std::map<ttk::SimplexId, ttk::SimplexId> vertGtoL;

  std::vector<Point> vertGhostCoordinates;
  std::vector<std::vector<Point>> vertGhostCoordinatesPerRank;
  std::vector<ttk::SimplexId> cellGhostGlobalVertexIds;
  std::vector<ttk::SimplexId> cellGhostLocalIds;
  std::vector<std::vector<ttk::SimplexId>> cellGhostLocalIdsPerRank;
  std::vector<std::vector<ttk::SimplexId>> cellGhostGlobalVertexIdsPerRank;
  std::vector<int> neighbors;
  std::map<int, int> neighborToId;
  double *boundingBox = input->GetBounds();
  getNeighborsUsingBoundingBox(boundingBox, neighbors);
  int neighborNumber = neighbors.size();
  for(int i = 0; i < neighborNumber; i++) {
    neighborToId[neighbors[i]] = i;
  }
  double bounds[6];
  double p[3];
  vtkIdList *points = vtkIdList::New();
  ttk::SimplexId nbPoints;
  MPI_Datatype mpiIdType = getMPIType(nbPoints);
  input->GetCellPoints(0, points);
  nbPoints = points->GetNumberOfIds();
  input->GetBounds(bounds);
  if(vertRankArray != nullptr) {
    for(SimplexId i = 0; i < vertexNumber; i++) {
      if(vertRankArray->GetTuple1(i) != ttk::MPIrank_) {
        realVertexNumber--;
        input->GetPoint(i, p);
        vertGhostCoordinatesPerRank[neighborToId[vertRankArray->GetTuple1(i)]]
          .push_back(Point{p[0], p[1], p[2], i});
      }
    }
  } else {
    for(SimplexId i = 0; i < vertexNumber; i++) {
      if(vertGhost->GetTuple1(i) != 0) {
        realVertexNumber--;
        input->GetPoint(i, p);
        vertGhostCoordinates.push_back(Point{p[0], p[1], p[2], i});
      }
    }
  }
  if(cellRankArray != nullptr) {
    for(SimplexId i = 0; i < cellNumber; i++) {
      if(cellRankArray->GetTuple1(i) != ttk::MPIrank_) {
        realCellNumber--;
      }
    }
  } else {
    for(SimplexId i = 0; i < cellNumber; i++) {
      if(cellGhost->GetTuple1(i) != 0) {
        realCellNumber--;
      }
    }
  }

  ttk::SimplexId vertIndex;
  ttk::SimplexId cellIndex;

  // Perform exclusive prefix sum
  MPI_Exscan(
    &realVertexNumber, &vertIndex, 1, mpiIdType, MPI_SUM, ttk::MPIcomm_);
  MPI_Exscan(&realCellNumber, &cellIndex, 1, mpiIdType, MPI_SUM, ttk::MPIcomm_);

  if(ttk::MPIrank_ == 0) {
    vertIndex = 0;
    cellIndex = 0;
  }

  if(vertRankArray != nullptr) {
    for(SimplexId i = 0; i < vertexNumber; i++) {
      if(vertRankArray->GetTuple1(i) == ttk::MPIrank_) {
        vertexIdentifiers->SetTuple1(i, vertIndex);
        vertGtoL[vertIndex] = i;
        vertIndex++;
      }
    }
  } else {
    for(SimplexId i = 0; i < vertexNumber; i++) {
      if(vertGhost->GetTuple1(i) == 0) {
        vertexIdentifiers->SetTuple1(i, vertIndex);
        vertGtoL[vertIndex] = i;

        vertIndex++;
      }
    }
  }
  if(cellRankArray != nullptr) {
    for(SimplexId i = 0; i < cellNumber; i++) {
      if(cellRankArray->GetTuple1(i) == ttk::MPIrank_) {
        cellIdentifiers->SetTuple1(i, cellIndex);
        cellIndex++;
      }
    }
  } else {
    for(SimplexId i = 0; i < cellNumber; i++) {
      if(cellGhost->GetTuple1(i) == 0) {
        cellIdentifiers->SetTuple1(i, cellIndex);
        cellIndex++;
      }
    }
  }

  // Handle Ghost exchange
  MPI_Datatype mpiPointType;
  createMPIPointType(&mpiPointType);
  MPI_Datatype mpiResponseType;
  createMPIResponseType(&mpiResponseType);
  ttk::SimplexId recvMessageSize{0};
  std::vector<Point> receivedPoints;
  std::vector<ttk::SimplexId> receivedCells;
  std::vector<Response> receivedResponse;
  std::vector<Response> locatedSimplices;

  for(int i = 0; i < neighborNumber; i++) {
    if(vertRankArray == nullptr) {
      ttkTemplateMacro(
        triangulation->getType(),
        this->exchangeAndLocatePoints(
          locatedSimplices, vertGhostCoordinates, receivedPoints,
          receivedResponse, neighbors[i], mpiPointType, mpiResponseType,
          recvMessageSize, input, bounds, vertexIdentifiers, vertGtoL,
          (TTK_TT *)triangulation->getData()));
      int count = 0;
      for(int n = 0; n < vertGhostCoordinates.size(); n++) {
        if(vertGhostCoordinates[n - count].localId
           == receivedResponse[count].id) {
          vertGhostCoordinates.erase(vertGhostCoordinates.begin() + n - count);
          count++;
        }
      }
    } else {
      ttkTemplateMacro(
        triangulation->getType(),
        this->exchangeAndLocatePoints(
          locatedSimplices, vertGhostCoordinatesPerRank[i], receivedPoints,
          receivedResponse, neighbors[i], mpiPointType, mpiResponseType,
          recvMessageSize, input, bounds, vertexIdentifiers, vertGtoL,
          (TTK_TT *)triangulation->getData()));
    }
  }
  printMsg("POINTS DONE, START CELLS");

  std::vector<std::vector<ttk::SimplexId>> pointsToCells(vertexNumber);
  vtkIdList *cells = vtkIdList::New();
  for(ttk::SimplexId i = 0; i < vertexNumber; i++) {
    input->GetPointCells(i, cells);
    for(int j = 0; j < cells->GetNumberOfIds(); j++) {
      if(cellGhost->GetTuple1(cells->GetId(j)) == 0) {
        pointsToCells[i].push_back(cells->GetId(j));
      }
    }
  }
  if(cellRankArray != nullptr) {
    for(SimplexId i = 0; i < cellNumber; i++) {
      if(cellRankArray->GetTuple1(i) != ttk::MPIrank_) {
        input->GetCellPoints(i, points);
        cellGhostLocalIdsPerRank[neighborToId[cellRankArray->GetTuple1(i)]]
          .push_back(i);
        for(int k = 0; k < nbPoints; k++) {
          cellGhostGlobalVertexIdsPerRank[neighborToId[cellRankArray->GetTuple1(
                                            i)]]
            .push_back(vertexIdentifiers->GetTuple1(points->GetId(k)));
        }
      }
    }
  } else {
    for(SimplexId i = 0; i < cellNumber; i++) {
      if(cellGhost->GetTuple1(i) != 0) {
        input->GetCellPoints(i, points);
        cellGhostLocalIds.push_back(i);
        for(int k = 0; k < nbPoints; k++) {
          cellGhostGlobalVertexIds.push_back(
            vertexIdentifiers->GetTuple1(points->GetId(k)));
        }
      }
    }
  }
  MPI_Barrier(ttk::MPIcomm_);
  printMsg("Start communication phase");

  for(int i = 0; i < neighborNumber; i++) {
    vtkIdList *point = vtkIdList::New();
    std::map<ttk::SimplexId, ttk::SimplexId>::iterator search;
    std::vector<ttk::SimplexId> localPointIds;
    localPointIds.reserve(nbPoints);
    locatedSimplices.clear();
    this->SendRecvVector<ttk::SimplexId>(cellGhostGlobalVertexIds,
                                         receivedCells, recvMessageSize,
                                         mpiIdType, neighbors[i]);
    printMsg("size globalVertexIds: "
             + std::to_string(cellGhostGlobalVertexIds.size()));
    printMsg("size receivedCells: " + std::to_string(receivedCells.size()));
    for(int n = 0; n < recvMessageSize; n += nbPoints) {
      localPointIds.clear();
      for(int k = 0; k < nbPoints; k++) {
        search = vertGtoL.find(receivedCells[n + k]);

        if(search != vertGtoL.end()) {
          localPointIds.push_back(search->second);
        } else {
          break;
        }
      }
      if(localPointIds.size() != static_cast<size_t>(nbPoints)) {
        break;
      }

      bool foundIt = false;
      int k = 0;
      int l;
      int m = 0;
      while(!foundIt && m < nbPoints) {
        int size = pointsToCells[localPointIds[m]].size();
        k = 0;
        while(!foundIt && k < size) {
          input->GetCellPoints(pointsToCells[localPointIds[m]][k], point);
          l = 0;
          while(l < nbPoints) {
            auto it = find(
              localPointIds.begin(), localPointIds.end(), point->GetId(l));
            if(it == localPointIds.end()) {
              break;
            }
            l++;
          }
          if(l == nbPoints) {
            foundIt = true;
            locatedSimplices.push_back(Response{
              n,
              cellIdentifiers->GetTuple1(pointsToCells[localPointIds[m]][k])});
          }
          k++;
        }
        m++;
      }
    }
    this->SendRecvVector<Response>(locatedSimplices, receivedResponse,
                                   recvMessageSize, mpiResponseType,
                                   neighbors[i]);
    printMsg("size: " + std::to_string(locatedSimplices.size()));

    for(int n = 0; n < recvMessageSize; n++) {
      cellIdentifiers->SetTuple1(
        cellGhostLocalIds[receivedResponse[n].id / nbPoints - n],
        receivedResponse[n].globalId);
      cellGhostLocalIds.erase(cellGhostLocalIds.begin()
                              + receivedResponse[n].id / nbPoints - n);
      cellGhostGlobalVertexIds.erase(cellGhostGlobalVertexIds.begin()
                                       + receivedResponse[n].id - n * nbPoints,
                                     cellGhostGlobalVertexIds.begin()
                                       + receivedResponse[n].id - n * nbPoints
                                       + nbPoints);
    }
  }
  MPI_Barrier(ttk::MPIcomm_);
  printMsg("cellGhostGlobalVertexIds size: "
           + std::to_string(cellGhostGlobalVertexIds.size()));

  // #ifdef TTK_ENABLE_OPENMP
  // #pragma omp parallel for num_threads(threadNumber_)
  // #endif
  //   for(SimplexId i = 0; i < vertexNumber; i++) {
  //     // avoid any processing if the abort signal is sent
  //     vertexIdentifiers->SetTuple1(i, i);
  //   }

  // #ifdef TTK_ENABLE_OPENMP
  // #pragma omp parallel for num_threads(threadNumber_)
  // #endif
  //   for(SimplexId i = 0; i < cellNumber; i++) {
  //     // avoid any processing if the abort signal is sent
  //     cellIdentifiers->SetTuple1(i, i);
  //   }

  output->GetPointData()->AddArray(vertexIdentifiers);
  output->GetCellData()->AddArray(cellIdentifiers);

  printMsg("Processed " + std::to_string(vertexNumber) + " vertices and "
             + std::to_string(cellNumber) + " cells",
           1, t.getElapsedTime(), threadNumber_);

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
