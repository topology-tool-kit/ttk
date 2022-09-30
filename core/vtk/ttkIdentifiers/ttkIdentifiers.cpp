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

void ttkIdentifiers::createMPIResponseCellType(MPI_Datatype *mpiResponseType) {
  ttk::SimplexId id = 0;
  MPI_Datatype types[] = {getMPIType(id), getMPIType(id), getMPIType(id)};
  int lengths[] = {1, 1, 1};
  const long int mpi_offsets[]
    = {offsetof(ResponseCell, localId), offsetof(ResponseCell, vectorId),
       offsetof(ResponseCell, globalId)};
  MPI_Type_create_struct(3, lengths, mpi_offsets, types, mpiResponseType);
  MPI_Type_commit(mpiResponseType);
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
  double bounds[6];
  double p[3];
  vtkIdList *points = vtkIdList::New();
  ttk::SimplexId nbPoints;
  input->GetCellPoints(0, points);
  nbPoints = points->GetNumberOfIds();
  // printMsg("nbPoints: "+std::to_string(nbPoints));
  input->GetBounds(bounds);
  printMsg("bounds: " + std::to_string(bounds[0]) + ", "
           + std::to_string(bounds[1]) + ", " + std::to_string(bounds[2]) + ", "
           + std::to_string(bounds[3]) + ", " + std::to_string(bounds[4]) + ", "
           + std::to_string(bounds[5]));
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
  MPI_Gather(comm_buf, 2, getMPIType(realVertexNumber), recv_offset.data(), 2,
             getMPIType(realVertexNumber), 0, ttk::MPIcomm_);
  offset_count[0] = 0;
  offset_count[1] = 0;
  for(ttk::SimplexId i = 1; i < ttk::MPIsize_; i++) {
    offset_count[2 * i] = offset_count[2 * (i - 1)] + recv_offset[2 * (i - 1)];
    offset_count[2 * i + 1] = offset_count[2 * i - 1] + recv_offset[2 * i - 1];
  }

  MPI_Scatter(offset_count.data(), 2, getMPIType(realVertexNumber), comm_buf, 2,
              getMPIType(realVertexNumber), 0, ttk::MPIcomm_);

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
        if(ttk::MPIrank_ == 1) {
          if(vertIndex == 40 || vertIndex == 55 || vertIndex == 65
             || vertIndex == 66) {
            printMsg("local index: " + std::to_string(i)
                     + " for global index: " + std::to_string(vertIndex));
          }
        }

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
  MPI_Datatype mpiResponseCellType;
  createMPIResponseCellType(&mpiResponseCellType);
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
  std::vector<ResponseCell> receivedResponseCell;
  std::vector<Response> locatedPoints;
  std::vector<ResponseCell> locatedCells;
  ttk::SimplexId id = -1;
  ttk::SimplexId globalId = -1;

  for(int i = 0; i < neighborNumber; i++) {
    receivedPoints.clear();
    locatedPoints.clear();
    receivedResponse.clear();
    vertGhostNumber = vertGhostCoordinates.size();
    MPI_Sendrecv(&vertGhostNumber, 1, getMPIType(vertGhostNumber), neighbors[i],
                 ttk::MPIrank_, &recvMessageSize, 1,
                 getMPIType(vertGhostNumber), neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);
    receivedPoints.resize(recvMessageSize);
    MPI_Sendrecv(vertGhostCoordinates.data(), vertGhostNumber, mpiPointType,
                 neighbors[i], ttk::MPIrank_, receivedPoints.data(),
                 recvMessageSize, mpiPointType, neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);
    for(int n = 0; n < recvMessageSize; n++) {
      if(bounds[0] <= receivedPoints[n].x && bounds[1] >= receivedPoints[n].x
         && bounds[2] <= receivedPoints[n].y && bounds[3] >= receivedPoints[n].y
         && bounds[4] <= receivedPoints[n].z
         && bounds[5] >= receivedPoints[n].z) {
        id = input->FindPoint(
          receivedPoints[n].x, receivedPoints[n].y, receivedPoints[n].z);

        if(id >= 0) {
          globalId = vertexIdentifiers->GetTuple1(id);
          if(globalId == 4641) {
            printMsg("PROBLEM, local id: " + std::to_string(id)
                     + " with coord: x:" + std::to_string(receivedPoints[n].x)
                     + " y: " + std::to_string(receivedPoints[n].y)
                     + " z:" + std::to_string(receivedPoints[n].z));
          }
          if(globalId >= 0) {
            locatedPoints.push_back(
              Response{receivedPoints[n].localId, globalId});
          }
        }
      }
    }
    ttk::SimplexId locatedPointNumber = locatedPoints.size();
    MPI_Sendrecv(&locatedPointNumber, 1, getMPIType(locatedPointNumber),
                 neighbors[i], ttk::MPIrank_, &recvMessageSize, 1,
                 getMPIType(locatedPointNumber), neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);

    receivedResponse.resize(recvMessageSize);

    MPI_Sendrecv(locatedPoints.data(), locatedPointNumber, mpiResponseType,
                 neighbors[i], ttk::MPIrank_, receivedResponse.data(),
                 recvMessageSize, mpiResponseType, neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);

    for(int n = 0; n < recvMessageSize; n++) {
      vertexIdentifiers->SetTuple1(
        receivedResponse[n].id, receivedResponse[n].globalId);
      vertGtoL[receivedResponse[n].globalId] = receivedResponse[n].id;
      if(ttk::MPIrank_ == 1) {
        if(receivedResponse[n].globalId == 40
           || receivedResponse[n].globalId == 55
           || receivedResponse[n].globalId == 65
           || receivedResponse[n].globalId == 66) {
          printMsg("local index: " + std::to_string(receivedResponse[n].id)
                   + " for global index: "
                   + std::to_string(receivedResponse[n].globalId));
        }
      }
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
  // printMsg("vertGhostSize: " + std::to_string(vertGhostCoordinates.size()));
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
  printMsg("PointsToCell done");
  if(cellRankArray != nullptr) {
    for(SimplexId i = 0; i < cellNumber; i++) {
      if(cellRankArray->GetTuple1(i) != ttk::MPIrank_) {
        input->GetCellPoints(i, points);
        cellGhostGlobalVertexIds.push_back(i);
        for(int k = 0; k < nbPoints; k++) {
          cellGhostGlobalVertexIds.push_back(
            vertexIdentifiers->GetTuple1(points->GetId(k)));
          if(i == 100 && ttk::MPIrank_ == 3) {
            printMsg(
              "Global id points: "
              + std::to_string(vertexIdentifiers->GetTuple1(points->GetId(k))));
            printMsg("Local id points: " + std::to_string(points->GetId(k)));
          }
        }
      }
    }
  } else {
    for(SimplexId i = 0; i < cellNumber; i++) {
      if(cellGhost->GetTuple1(i) != 0) {
        input->GetCellPoints(i, points);
        cellGhostGlobalVertexIds.push_back(i);
        for(int k = 0; k < nbPoints; k++) {
          cellGhostGlobalVertexIds.push_back(
            vertexIdentifiers->GetTuple1(points->GetId(k)));
          if(i == 100 && ttk::MPIrank_ == 3) {
            printMsg(
              "Global id points: "
              + std::to_string(vertexIdentifiers->GetTuple1(points->GetId(k))));
            printMsg("Local id points: " + std::to_string(points->GetId(k)));
          }
          // printMsg("point to send
          // :"+std::to_string(cellGhostGlobalVertexIds.back()));
        }
      }
    }
  }
  MPI_Barrier(ttk::MPIcomm_);
  printMsg("Start communication phase");
  std::map<ttk::SimplexId, ttk::SimplexId>::iterator search;
  // printMsg("number of cells: "+std::to_string(input->GetNumberOfCells()));
  globalId = -1;
  // MPI_Barrier(ttk::MPIcomm_);
  std::vector<ttk::SimplexId> localPointIds;
  localPointIds.reserve(nbPoints);
  for(int i = 0; i < neighborNumber; i++) {
    std::string s = "";
    for(int h = 0; h < cellGhostGlobalVertexIds.size(); h += nbPoints + 1) {
      for(int g = 0; g < nbPoints + 1; g++) {
        s += std::to_string(cellGhostGlobalVertexIds[h + g]) + ", ";
      }
      s += "\n";
    }
    if(ttk::MPIrank_ == 3)
      printMsg(s);
    printMsg("Round: " + std::to_string(i)
             + ", neighbor: " + std::to_string(neighbors[i]));
    receivedCells.clear();
    receivedResponseCell.clear();
    locatedCells.clear();
    cellGhostNumber = cellGhostGlobalVertexIds.size();
    printMsg("cellGhostNumber: " + std::to_string(cellGhostNumber));
    MPI_Sendrecv(&cellGhostNumber, 1, getMPIType(cellGhostNumber), neighbors[i],
                 ttk::MPIrank_, &recvMessageSize, 1,
                 getMPIType(cellGhostNumber), neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);
    // MPI_Barrier(ttk::MPIcomm_);
    // printMsg("recvMessageSize for cellGhostGlobalVertexIds: "
    //            + std::to_string(recvMessageSize));
    receivedCells.resize(recvMessageSize);
    // printMsg("Message size sent for ghost cells");
    MPI_Sendrecv(cellGhostGlobalVertexIds.data(), cellGhostNumber,
                 getMPIType(realVertexNumber), neighbors[i], ttk::MPIrank_,
                 receivedCells.data(), recvMessageSize,
                 getMPIType(realVertexNumber), neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);
    // printMsg("Data sent");
    for(int n = 0; n < recvMessageSize; n += (nbPoints + 1)) {
      int localCellId = receivedCells[n];
      // if(localCellId == 100 && neighbors[i] == 3) {
      //   printMsg("DAAAAAAA");
      // }
      localPointIds.clear();
      for(int k = 1; k < nbPoints + 1; k++) {
        search = vertGtoL.find(receivedCells[n + k]);
        if(localCellId == 100 && neighbors[i] == 3) {
          printMsg("received global id[" + std::to_string(receivedCells[n + k])
                   + "]: " + std::to_string(receivedCells[n + k]));
        }

        if(search != vertGtoL.end()) {
          localPointIds.push_back(search->second);
        } else {
          // if(localCellId == 100 && neighbors[i] == 3) {
          //   printMsg("WE GOT A PROBLEM1");
          // }
          break;
        }
      }
      if(localPointIds.size() != static_cast<size_t>(nbPoints)) {
        // if(localCellId == 12000 && neighbors[i] == 0) {
        //   printMsg("WE GOT A PROBLEM2");
        // }
        break;
        // printMsg("Da0");
      }
      // else {
      //   printMsg("OK");
      // }

      bool foundIt = false;
      int k = 0;
      int l;
      int m = 0;
      while(!foundIt && m < nbPoints) {
        int size = pointsToCells[localPointIds[m]].size();
        k = 0;
        while(!foundIt && k < size) {
          // printMsg("local cell id here:
          // "+std::to_string(pointsToCells[localPointIds[m]][k]));
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
            locatedCells.push_back(ResponseCell{
              localCellId, n,
              cellIdentifiers->GetTuple1(pointsToCells[localPointIds[m]][k])});
            // if(localCellId == 100 && neighbors[i] == 3) {
            // //   printMsg("Found it");
            //   printMsg("Global cell id:
            //   "+std::to_string(locatedCells.back().globalId));
            // }
          } else {
            if(localCellId == 100 && neighbors[i] == 3) {
              // //   printMsg("Not found: ");
              //    printMsg("local cell id here: "
              //             +
              //             std::to_string(pointsToCells[localPointIds[m]][k]));
              if(pointsToCells[localPointIds[m]][k] == 44) {
                for(int f = 0; f < localPointIds.size(); f++) {
                  printMsg("localPointIds[" + std::to_string(f)
                           + "]: " + std::to_string(localPointIds[f]));
                }
              }
            }
            //   for(int f = 0; f < points->GetNumberOfIds(); f++) {
            //     printMsg("points[" + std::to_string(f)
            //              + "]: " + std::to_string(points->GetId(f)));
            //   }
            //   }
            //   for(int f = 0; f < localPointIds.size(); f++) {
            //     printMsg("localPointIds[" + std::to_string(f)
            //              + "]: " + std::to_string(localPointIds[f]));
            //   }

            //   printMsg("local cell id: " + std::to_string(localCellId));
            //   for(int f = 0; f < points->GetNumberOfIds(); f++) {
            //     printMsg("points[" + std::to_string(f)
            //              + "]: " + std::to_string(points->GetId(f)));
            //   }
            // }
            // if(localCellId == 100 && neighbors[i] == 3) {
            //  printMsg("Didn't find it");
            //}
          }
          k++;
        }
        m++;
      }
    }
    // printMsg("Cells located");
    ttk::SimplexId locatedPointNumber = locatedCells.size();
    MPI_Sendrecv(&locatedPointNumber, 1, getMPIType(locatedPointNumber),
                 neighbors[i], ttk::MPIrank_, &recvMessageSize, 1,
                 getMPIType(locatedPointNumber), neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);

    receivedResponseCell.resize(recvMessageSize);

    MPI_Sendrecv(locatedCells.data(), locatedPointNumber, mpiResponseCellType,
                 neighbors[i], ttk::MPIrank_, receivedResponseCell.data(),
                 recvMessageSize, mpiResponseCellType, neighbors[i],
                 neighbors[i], ttk::MPIcomm_, MPI_STATUS_IGNORE);
    // printMsg("recvMessageSize for located cells: " +
    // std::to_string(recvMessageSize));
    for(int n = 0; n < recvMessageSize; n++) {
      cellIdentifiers->SetTuple1(
        receivedResponseCell[n].localId, receivedResponseCell[n].globalId);
      // if(receivedResponseCell[n].id == 100 && ttk::MPIrank_ == 3) {
      //   printMsg("Here's your associated global id:
      //   "+std::to_string(receivedResponseCell[n].globalId));
      // }
      cellGhostGlobalVertexIds.erase(
        cellGhostGlobalVertexIds.begin() + receivedResponseCell[n].vectorId
          - n * (nbPoints + 1),
        cellGhostGlobalVertexIds.begin() + receivedResponseCell[n].vectorId
          - n * (nbPoints + 1) + nbPoints + 1);
    }

    // printMsg("About to delete");
    // int count = 0;
    // for(int n = 0; n < cellGhostGlobalVertexIds.size(); n++) {
    //      if(receivedResponse[count].id == 100 && ttk::MPIrank_ == 3) {
    //        printMsg("WE GOT A PROBLEM3");
    //      }
    //   if(cellGhostGlobalVertexIds[n - count*(nbPoints + 1)] ==
    //   receivedResponse[count].id) {

    //     count++;
    //   }
    // }
  }
  MPI_Barrier(ttk::MPIcomm_);
  printMsg("cellGhostGlobalVertexIds size: "
           + std::to_string(cellGhostGlobalVertexIds.size()));
  // / printMsg("first : "+std::to_string(cellGhostGlobalVertexIds[0].x)+"
  // / "+std::to_string(cellGhostGlobalVertexIds[0].y)+"
  // / "+std::to_string(cellGhostGlobalVertexIds[0].z));

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
