#include <ttkIdentifiers.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <map>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

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
  MPI_Datatype types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  int lengths[] = {1, 1, 1};
  const long int mpi_offsets[]
    = {offsetof(Point, x), offsetof(Point, y), offsetof(Point, z)};
  MPI_Type_create_struct(3, lengths, mpi_offsets, types, mpiPointType);
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

int ttkIdentifiers::RequestData(vtkInformation *ttkNotUsed(request),
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

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

  printMsg("DA0");

  std::vector<Point> vertGhostCoordinates;
  std::vector<ttk::SimplexId> cellGhostGlobalVertexIds;
  double bounds[6];
  double p[3];
  vtkIdList *points = vtkIdList::New();
  ttk::SimplexId nbPoints;
  input->GetCellPoints(0, points);
  nbPoints = points->GetNumberOfIds();
  input->GetCellBounds(0, bounds);
  if(vertRankArray != nullptr) {
    for(SimplexId i = 0; i < vertexNumber; i++) {
      if(vertRankArray->GetTuple1(i) != ttk::MPIrank_) {
        realVertexNumber--;
        input->GetPoint(i, p);
        vertGhostCoordinates.push_back(Point{p[0], p[1], p[2]});
      }
    }
  } else {
    for(SimplexId i = 0; i < vertexNumber; i++) {
      if(vertGhost->GetTuple1(i) != 0) {
        realVertexNumber--;
        input->GetPoint(i, p);
        vertGhostCoordinates.push_back(Point{p[0], p[1], p[2]});
        // printMsg("x: "+std::to_string(vertGhostCoordinates.back().x)+" y:
        // "+std::to_string(vertGhostCoordinates.back().y)+" z:
        // "+std::to_string(vertGhostCoordinates.back().z));
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
  std::vector<Response> locatedPoints;
  ttk::SimplexId id = -1;
  ttk::SimplexId globalId = -1;

  for(int i = 0; i < neighborNumber; i++) {
    // printMsg("Round: "+std::to_string(i));
    receivedPoints.clear();
    locatedPoints.clear();
    receivedResponse.clear();
    vertGhostNumber = vertGhostCoordinates.size();
    printMsg("vertGhostnumber: " + std::to_string(vertGhostNumber));
    MPI_Sendrecv(&vertGhostNumber, 1, getMPIType(vertGhostNumber), neighbors[i],
                 ttk::MPIrank_, &recvMessageSize, 1,
                 getMPIType(vertGhostNumber), neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);
    printMsg("recvMessageSize: " + std::to_string(recvMessageSize));
    receivedPoints.resize(recvMessageSize);
    printMsg("Message size sent for ghost points");
    MPI_Sendrecv(vertGhostCoordinates.data(), vertGhostNumber, mpiPointType,
                 neighbors[i], ttk::MPIrank_, receivedPoints.data(),
                 recvMessageSize, mpiPointType, neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);
    printMsg("Message sent for ghost points");
    for(int n = 0; n < recvMessageSize; n++) {
      id = input->FindPoint(
        receivedPoints[n].x, receivedPoints[n].y, receivedPoints[n].z);
      printMsg("x: " + std::to_string(receivedPoints[n].x)
               + " y: " + std::to_string(receivedPoints[n].y)
               + " z: " + std::to_string(receivedPoints[n].z));
      if(id >= 0) {
        globalId = vertexIdentifiers->GetTuple1(id);
        if(globalId >= 0) {
          locatedPoints.push_back(Response{n, globalId});
          // printMsg("located points: "+std::to_string(id)+", global:
          // "+std::to_string(globalId));
        }
      }
    }
    // printMsg("Points located");
    ttk::SimplexId locatedPointNumber = locatedPoints.size();
    MPI_Sendrecv(&locatedPointNumber, 1, getMPIType(locatedPointNumber),
                 neighbors[i], ttk::MPIrank_, &recvMessageSize, 1,
                 getMPIType(locatedPointNumber), neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);
    // printMsg("Message size sent for located points");
    receivedResponse.resize(recvMessageSize);

    MPI_Sendrecv(locatedPoints.data(), locatedPointNumber, mpiResponseType,
                 neighbors[i], ttk::MPIrank_, receivedResponse.data(),
                 recvMessageSize, mpiResponseType, neighbors[i], neighbors[i],
                 ttk::MPIcomm_, MPI_STATUS_IGNORE);
    // printMsg("Message sent for located points");
    // MPI_Barrier(ttk::MPIcomm_);
    for(int n = 0; n < recvMessageSize; n++) {
      id = input->FindPoint(vertGhostCoordinates[receivedResponse[n].id].x,
                            vertGhostCoordinates[receivedResponse[n].id].y,
                            vertGhostCoordinates[receivedResponse[n].id].z);
      if(id >= 0) {
        vertexIdentifiers->SetTuple1(id, receivedResponse[n].globalId);
        vertGtoL[receivedResponse[n].globalId] = id;
      }
    }
    // printMsg("recvMessageSize of located Points:
    // "+std::to_string(recvMessageSize)+" size of vertGhost:
    // "+std::to_string(vertGhostCoordinates.size())+" from: "
    // +std::to_string(neighbors[i]));
    for(int n = 0; n < recvMessageSize; n++) {
      vertGhostCoordinates.erase(vertGhostCoordinates.begin()
                                 + receivedResponse[n].id - n);
    }
  }
  printMsg("vertGhostSize: " + std::to_string(vertGhostCoordinates.size()));
  printMsg("POINTS DONE, START CELLS");
  // if (cellRankArray != nullptr){
  //   for(SimplexId i = 0; i < cellNumber; i++) {
  //     if (cellRankArray->GetTuple1(i) != ttk::MPIrank_){
  //       input->GetCellPoints(i, points);
  //       for (int k = 0; k < nbPoints; k++){
  //         cellGhostGlobalVertexIds.push_back(vertexIdentifiers->GetTuple1(points->GetId(k)));
  //         //printMsg("point to send :
  //         "+std::to_string(cellGhostGlobalVertexIds.back()));
  //       }
  //     }
  //   }
  // }
  // else {
  //   for(SimplexId i = 0; i < cellNumber; i++) {
  //     if (cellGhost->GetTuple1(i) != 0){
  //       input->GetCellPoints(i, points);
  //       for (int k = 0; k < nbPoints; k++){
  //         cellGhostGlobalVertexIds.push_back(vertexIdentifiers->GetTuple1(points->GetId(k)));
  //         //printMsg("point to send :
  //         "+std::to_string(cellGhostGlobalVertexIds.back()));
  //       }
  //     }
  //   }
  // }
  // vtkIdList* intersectionCells = vtkIdList::New();
  // vtkIdList* cells = vtkIdList::New();
  // std::map<ttk::SimplexId,ttk::SimplexId>::iterator search;
  // printMsg("number of cells: "+std::to_string(input->GetNumberOfCells()));
  // for(int i = 0; i < neighborNumber; i++) {
  //   printMsg("Round: "+std::to_string(i));
  //   receivedCells.clear();
  //   receivedResponse.clear();
  //   locatedPoints.clear();
  //   cellGhostNumber = cellGhostGlobalVertexIds.size();
  //   printMsg("cellGhostNumber: "+std::to_string(cellGhostNumber));
  //   MPI_Sendrecv(&cellGhostNumber, 1, getMPIType(cellGhostNumber),
  //   neighbors[i], ttk::MPIrank_,
  //             &recvMessageSize, 1, getMPIType(cellGhostNumber), neighbors[i],
  //             neighbors[i], ttk::MPIcomm_, MPI_STATUS_IGNORE);
  //   printMsg("recvMessageSize: "+std::to_string(recvMessageSize));
  //   receivedCells.resize(recvMessageSize);
  //   printMsg("Message size sent for ghost cells");
  //   MPI_Sendrecv(cellGhostGlobalVertexIds.data(), cellGhostNumber,
  //   getMPIType(realVertexNumber), neighbors[i], ttk::MPIrank_,
  //             receivedCells.data(), recvMessageSize,
  //             getMPIType(realVertexNumber), neighbors[i], neighbors[i],
  //             ttk::MPIcomm_, MPI_STATUS_IGNORE);
  //   printMsg("Message sent for ghost cells");
  //   for (int n = 0; n < recvMessageSize/nbPoints; n++){
  //     //center[0] = receivedPoints[n].x; center[1] = receivedPoints[n].y;
  //     center[2] = receivedPoints[n].z;
  //     //printMsg("point to find :
  //     "+std::to_string(receivedCells[n*nbPoints]));
  //     //printMsg("map result :
  //     "+std::to_string(vertGtoL[receivedCells[n*nbPoints]])); search =
  //     vertGtoL.find(receivedCells[n*nbPoints]); if (search !=
  //     vertGtoL.end()){
  //       //printMsg("Da0");
  //       input->GetPointCells(vertGtoL[receivedCells[n*nbPoints]],
  //       intersectionCells); int k = 1; while (k < nbPoints &&
  //       intersectionCells->GetNumberOfIds() > 1){
  //         search = vertGtoL.find(receivedCells[n*nbPoints+k]);
  //         if (search != vertGtoL.end()){
  //           input->GetPointCells(vertGtoL[receivedCells[n*nbPoints+k]],
  //           cells); intersectionCells->IntersectWith(cells); k++;
  //         } else {
  //           break;
  //         }
  //       }
  //       if (intersectionCells->GetNumberOfIds() == 1){
  //         printMsg("Daaaaa");
  //         locatedPoints.push_back(Response{n*nbPoints,
  //         cellIdentifiers->GetTuple1(intersectionCells->GetId(0))});
  //       }
  //       else {
  //         printMsg("Size of intersection:
  //         "+std::to_string(intersectionCells->GetNumberOfIds()));
  //       }
  //     }
  //   }
  //   //printMsg("cells located");
  //   ttk::SimplexId locatedPointNumber = locatedPoints.size();
  //   MPI_Sendrecv(&locatedPointNumber, 1, getMPIType(locatedPointNumber),
  //   neighbors[i], ttk::MPIrank_,
  //             &recvMessageSize, 1, getMPIType(locatedPointNumber),
  //             neighbors[i], neighbors[i], ttk::MPIcomm_, MPI_STATUS_IGNORE);
  //   //printMsg("Message size sent for located cells");
  //   receivedResponse.resize(recvMessageSize);

  //   MPI_Sendrecv(locatedPoints.data(), locatedPointNumber , mpiResponseType,
  //   neighbors[i], ttk::MPIrank_,
  //            receivedResponse.data(), recvMessageSize, mpiResponseType,
  //            neighbors[i], neighbors[i], ttk::MPIcomm_,
  //             MPI_STATUS_IGNORE);
  //   //printMsg("Message sent for located cells");
  //   //MPI_Barrier(ttk::MPIcomm_);
  //   for (int n = 0; n < recvMessageSize; n++){
  //     //center[0] = receivedPoints[n].x; center[1] = receivedPoints[n].y;
  //     center[2] = receivedPoints[n].z; search =
  //     vertGtoL.find(cellGhostGlobalVertexIds[receivedResponse[n].id*nbPoints]);
  //     if (search != vertGtoL.end()){
  //       input->GetPointCells(vertGtoL[cellGhostGlobalVertexIds[receivedResponse[n].id*nbPoints]],
  //       intersectionCells); int k = 1; while (k < nbPoints &&
  //       intersectionCells->GetNumberOfIds() > 1){
  //         search =
  //         vertGtoL.find(cellGhostGlobalVertexIds[receivedResponse[n].id*nbPoints+k]);
  //         if (search != vertGtoL.end()){
  //           input->GetPointCells(vertGtoL[cellGhostGlobalVertexIds[receivedResponse[n].id*nbPoints+k]],
  //           cells); intersectionCells->IntersectWith(cells); k++;
  //         } else {
  //           break;
  //         }
  //       }
  //       if (intersectionCells->GetNumberOfIds() == 1){
  //         cellIdentifiers->SetTuple1(intersectionCells->GetId(0),
  //         receivedResponse[n].globalId);
  //       }
  //     }
  //   }
  //   printMsg("recvMessageSize of located cells:
  //   "+std::to_string(recvMessageSize)+" size of vertGhost:
  //   "+std::to_string(cellGhostGlobalVertexIds.size())+" from: "
  //   +std::to_string(neighbors[i])); for (int n = 0; n < recvMessageSize;
  //   n++){
  //     cellGhostGlobalVertexIds.erase(cellGhostGlobalVertexIds.begin()+receivedResponse[n].id-n);
  //   }
  // }
  // MPI_Barrier(ttk::MPIcomm_);
  // printMsg("cellGhostGlobalVertexIds size:
  // "+std::to_string(cellGhostGlobalVertexIds.size()));
  /// printMsg("first : "+std::to_string(cellGhostGlobalVertexIds[0].x)+"
  /// "+std::to_string(cellGhostGlobalVertexIds[0].y)+"
  /// "+std::to_string(cellGhostGlobalVertexIds[0].z));

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
