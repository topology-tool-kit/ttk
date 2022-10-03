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

int ttkIdentifiers::RequestData(vtkInformation *ttkNotUsed(request),
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);
  // int exactExtent;
  // vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  // exactExtent =
  // inInfo->Get(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT());
  //  printMsg("exact extent: "+std::to_string(exactExtent));
  // int wholeExtent[6];
  // inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent);
  // int ghost_levels =
  // inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  // int piece =
  // inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  // printMsg("ghost level: "+std::to_string(ghost_levels));
  // printMsg("hasAnyGhostCells: "+std::to_string(input->HasAnyGhostCells()));
  // double* origin;
  // vtkImageData* data = vtkImageData::SafeDownCast(input);
  // data->GetExtent(wholeExtent);
  // data->ComputeBounds();
  // double* bounds = data->GetBounds();
  // origin = data->GetOrigin();
  // printMsg("number of points: "+std::to_string(data->GetNumberOfPoints()));
  // if (ttk::MPIrank_ == 0)
  // {
  //   cout << "Whole extent: " << wholeExtent[0] << "," << wholeExtent[1] << "
  //   " << wholeExtent[2]         << "," << wholeExtent[3] << " " <<
  //   wholeExtent[4] << "," <<
  //        wholeExtent[5] << endl;
  // }
  //   cout << "bounds: " << bounds[0] << "," << bounds[1] << "  " << bounds[2]
  //        << "," << bounds[3] << " " << bounds[4] << "," <<
  //        bounds[5] << endl;
  // double* point = data->GetPoint(0);
  // printMsg("Point: "+std::to_string(point[0])+" "+std::to_string(point[1])+"
  // "+std::to_string(point[2])); data->SetExtent(wholeExtent); data->Update();
  // printMsg("number of points: "+std::to_string(data->GetNumberOfPoints()));
  // printMsg("piece number: "+std::to_string(piece));
  // double* point;
  // int point_ijk[3];
  // data->GetDimensions(point_ijk);
  // //data->ComputeStructuredCoordinates(point, point_ijk, pcoords);
  // printMsg("x: "+std::to_string(origin[0])+" y: "+std::to_string(origin[1])+"
  // z: "+std::to_string(origin[2]));
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
  std::vector<ttk::SimplexId> cellGhostGlobalVertexIds;
  std::vector<ttk::SimplexId> cellGhostLocalIds;
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
        vertGhostCoordinates.push_back(Point{p[0], p[1], p[2], i});
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
  ttk::SimplexId comm_buf[2] = {realVertexNumber, realCellNumber};
  std::vector<ttk::SimplexId> recv_offset(2 * ttk::MPIsize_);
  std::vector<ttk::SimplexId> offset_count(2 * ttk::MPIsize_);
  MPI_Gather(
    comm_buf, 2, mpiIdType, recv_offset.data(), 2, mpiIdType, 0, ttk::MPIcomm_);
  offset_count[0] = 0;
  offset_count[1] = 0;
  for(ttk::SimplexId i = 1; i < ttk::MPIsize_; i++) {
    offset_count[2 * i] = offset_count[2 * (i - 1)] + recv_offset[2 * (i - 1)];
    offset_count[2 * i + 1] = offset_count[2 * i - 1] + recv_offset[2 * i - 1];
  }

  MPI_Scatter(offset_count.data(), 2, mpiIdType, comm_buf, 2, mpiIdType, 0,
              ttk::MPIcomm_);

  ttk::SimplexId vertIndex = comm_buf[0];
  ttk::SimplexId cellIndex = comm_buf[1];

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
  ttk::SimplexId vertGhostNumber;
  ttk::SimplexId cellGhostNumber;
  double *boundingBox = input->GetBounds();
  std::vector<int> neighbors;
  getNeighborsUsingBoundingBox(boundingBox, neighbors);
  int neighborNumber = neighbors.size();
  ttk::SimplexId recvMessageSize;
  std::vector<Point> receivedPoints;
  std::vector<ttk::SimplexId> receivedCells;
  std::vector<Response> receivedResponse;
  std::vector<Response> locatedSimplices;
  ttk::SimplexId id = -1;
  ttk::SimplexId globalId = -1;

  for(int i = 0; i < neighborNumber; i++) {
    locatedSimplices.clear();
    this->SendRecvVector<Point>(vertGhostCoordinates, receivedPoints,
                                recvMessageSize, mpiPointType, neighbors[i]);
    for(int n = 0; n < recvMessageSize; n++) {
      if(bounds[0] <= receivedPoints[n].x && bounds[1] >= receivedPoints[n].x
         && bounds[2] <= receivedPoints[n].y && bounds[3] >= receivedPoints[n].y
         && bounds[4] <= receivedPoints[n].z
         && bounds[5] >= receivedPoints[n].z) {
        id = input->FindPoint(
          receivedPoints[n].x, receivedPoints[n].y, receivedPoints[n].z);

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
                                   recvMessageSize, mpiResponseType,
                                   neighbors[i]);
    for(int n = 0; n < recvMessageSize; n++) {
      vertexIdentifiers->SetTuple1(
        receivedResponse[n].id, receivedResponse[n].globalId);
      vertGtoL[receivedResponse[n].globalId] = receivedResponse[n].id;
    }

    int count = 0;
    for(int n = 0; n < vertGhostCoordinates.size(); n++) {
      if(vertGhostCoordinates[n - count].localId
         == receivedResponse[count].id) {
        vertGhostCoordinates.erase(vertGhostCoordinates.begin() + n - count);
        count++;
      }
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
        cellGhostLocalIds.push_back(i);
        for(int k = 0; k < nbPoints; k++) {
          cellGhostGlobalVertexIds.push_back(
            vertexIdentifiers->GetTuple1(points->GetId(k)));
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
  std::map<ttk::SimplexId, ttk::SimplexId>::iterator search;
  globalId = -1;
  std::vector<ttk::SimplexId> localPointIds;
  localPointIds.reserve(nbPoints);
  for(int i = 0; i < neighborNumber; i++) {
    locatedSimplices.clear();
    this->SendRecvVector<ttk::SimplexId>(cellGhostGlobalVertexIds,
                                         receivedCells, recvMessageSize,
                                         mpiIdType, neighbors[i]);

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
          input->GetCellPoints(pointsToCells[localPointIds[m]][k], points);
          l = 0;
          while(l < nbPoints) {
            auto it = find(
              localPointIds.begin(), localPointIds.end(), points->GetId(l));
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
