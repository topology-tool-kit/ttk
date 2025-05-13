#include <ttkTrajectoryStatistics.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkTable.h> 
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <algorithm> // std::sort, std::unique
#include <vector> 
#include <map>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkTrajectoryStatistics);


ttkTrajectoryStatistics::ttkTrajectoryStatistics() {
  this->setDebugMsgPrefix("TrajectoryStatistics");
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}


int ttkTrajectoryStatistics::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }

  if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),
              "vtkDataSet");
    return 1;
  }

  return 0;
}


int ttkTrajectoryStatistics::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }

  if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }

  return 0;
}


int ttkTrajectoryStatistics::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  vtkUnstructuredGrid *inputGrid = vtkUnstructuredGrid::GetData(inputVector[0]);
  if (!inputGrid) {
    this->printErr("Invalid input : vtkUnstructuredGrid missing");
    return 0;
  }

  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[1]);
  if (!inputDataSet) {
    this->printErr("Invalid input : vtkDataSet missing");
    return 0;
  }


  vtkIntArray *compIdArray = vtkIntArray::SafeDownCast(
            inputGrid->GetCellData()->GetArray("ConnectedComponentId"));
  vtkIntArray *timeArray = vtkIntArray::SafeDownCast(
            inputGrid->GetPointData()->GetArray("TimeStep"));
  vtkIntArray *vertexGlobalIdArray = vtkIntArray::SafeDownCast(
            inputGrid->GetPointData()->GetArray("VertexGlobalId"));
  vtkIntArray *compLength = vtkIntArray::SafeDownCast(
            inputGrid->GetCellData()->GetArray("ComponentLength"));

  if (!compIdArray || !timeArray || !vertexGlobalIdArray || !compLength){
    this->printErr("Missing data in input vtu");
    return 0;
  }


  this->printMsg("Extraction ConnectedComponentId && TimeStep && VertexGlobalId done");

  vtkIdType numCells = inputGrid->GetNumberOfCells();

  if (numCells == 0)
    return 1;

  // map avec : {id de la traj, [cellsID, ...]}
  std::map<int, std::vector<vtkIdType>> groupTraj;

  for (vtkIdType cellId = 0; cellId < numCells; ++cellId){
    int trajId = compIdArray->GetValue(cellId); // ConnectedComponentId de la traj
    groupTraj[trajId].push_back(cellId);
  }

  this->printMsg("Traj Id Array created");

  const size_t numTraj = groupTraj.size();

  std::vector<std::vector<int>>    trajTime(numTraj); // TimeStep
  std::vector<std::vector<double>> trajX(numTraj), trajY(numTraj), trajZ(numTraj); // coord
  std::vector<std::vector<int>>    trajVertexId(numTraj);  // VertexGlobalId

  vtkNew<vtkIdList> cellPointsIds;
  size_t trajIndex = 0;
  for (const auto &trajEntry : groupTraj){
    
    const std::vector<vtkIdType> &cellIds = trajEntry.second;
    std::vector<vtkIdType> pointsIds; 
    pointsIds.reserve(cellIds.size()*2);


    // recup tt les points de la traj
    for (vtkIdType cellId : cellIds){
        cellPointsIds->Reset();
        inputGrid->GetCellPoints(cellId, cellPointsIds);
        vtkIdType numPts = cellPointsIds->GetNumberOfIds();
        for (vtkIdType i=0; i<numPts; ++i){
            pointsIds.push_back(cellPointsIds->GetId(i));
        }   
    }


    std::sort(pointsIds.begin(), pointsIds.end());
    pointsIds.erase(std::unique(pointsIds.begin(), pointsIds.end()), 
                    pointsIds.end()); // doublon
    std::sort(pointsIds.begin(), pointsIds.end(), [&](vtkIdType a, vtkIdType b){
        return timeArray->GetValue(a) < timeArray->GetValue(b);
    }); // tri TimeStep

    const size_t numPoints = pointsIds.size();
    trajTime[trajIndex].reserve(numPoints);
    trajX[trajIndex].reserve(numPoints);
    trajY[trajIndex].reserve(numPoints);
    trajZ[trajIndex].reserve(numPoints);
    trajVertexId[trajIndex].reserve(numPoints); 

    for (vtkIdType pointId : pointsIds){
        int t = timeArray->GetValue(pointId);
        double coords[3];
        int globalId = vertexGlobalIdArray->GetValue(pointId);
        inputGrid->GetPoint(pointId, coords);
        trajTime[trajIndex].push_back(t);
        trajX[trajIndex].push_back(coords[0]);
        trajY[trajIndex].push_back(coords[1]);
        trajZ[trajIndex].push_back(coords[2]);
        trajVertexId[trajIndex].push_back(globalId);
    }
    ++trajIndex;
  }

  this->printMsg("Pre-fetching done, "+ std::to_string(trajIndex) + " unique trajectory find");

  //Appel core/base
  std::vector<int> startFrames(numTraj), endFrames(numTraj), durations(numTraj);
  std::vector<double> VX(numTraj), VY(numTraj), surfMin(numTraj), surfMax(numTraj), surfMoy(numTraj);
   

  std::vector<vtkDataArray *> inputScalarFieldsRaw;
  std::vector<vtkDataArray *> inputScalarFields;
  const auto pointData= inputDataSet->GetPointData();

  if (!pointData){
    this->printErr("scalarArray missing");
    return 0;
  }

  int numberOfInputFields = pointData->GetNumberOfArrays();
  this->printMsg("numFrame = " + std::to_string(numberOfInputFields));
  
  vtkDataArray *firstScalarField = pointData->GetArray(0);

  for(int i = 0; i < numberOfInputFields; ++i) {
    vtkDataArray *currentScalarField = pointData->GetArray(i);
    if(currentScalarField == nullptr
       || currentScalarField->GetName() == nullptr) {
      continue;
    }
    std::string const sfname{currentScalarField->GetName()};
    if(sfname.rfind("_Order") == (sfname.size() - 6)) {
      continue;
    }
    if(firstScalarField->GetDataType() != currentScalarField->GetDataType()) {
      this->printErr("Inconsistent field data type or size between fields `"
                     + std::string{firstScalarField->GetName()} + "' and `"
                     + sfname + "'");
      return -1;
    }
    inputScalarFieldsRaw.push_back(currentScalarField);
  }

  std::sort(inputScalarFieldsRaw.begin(), inputScalarFieldsRaw.end(),
            [](vtkDataArray *a, vtkDataArray *b) {
              std::string s1 = a->GetName();
              std::string s2 = b->GetName();
              return std::lexicographical_compare(
                s1.begin(), s1.end(), s2.begin(), s2.end());
            });

  numberOfInputFields = inputScalarFieldsRaw.size();
  this->printMsg("New number of frame = " + std::to_string(numberOfInputFields));
  for(int i = 0; i < numberOfInputFields ; i++) {
    vtkDataArray *currentScalarField = inputScalarFieldsRaw[i];
    // Print scalar field names:
    // std::cout << currentScalarField->GetName() << std::endl;
    inputScalarFields.push_back(currentScalarField);
  }

  const int nFields = static_cast<int>(inputScalarFields.size());
  if(nFields == 0) {
    this->printErr("No scalar fields selected after sampling.");
    return 0;
  }

  vtkIdType nPts = inputScalarFields[0]->GetNumberOfTuples();
  this->printMsg("nombre de points : " + std::to_string(nPts));


  for(int f = 1; f < nFields; ++f) {
    if(inputScalarFields[f]->GetNumberOfTuples() != nPts) {
        this->printErr("Scalar fields have inconsistent number of points.");
        return 0;
    }
  }

  std::vector<std::vector<double>> vertexScalars(
    nPts, std::vector<double>(nFields)
  );

  for(int f = 0; f < nFields; ++f) {
    vtkDataArray *fieldArr = inputScalarFields[f];
    for(vtkIdType pid = 0; pid < nPts; ++pid) {
      vertexScalars[pid][f] = fieldArr->GetTuple1(pid);
    }
  }

  this->printMsg("Scalars recup");

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);

  int status = 0;

  status = this->execute(
                    trajTime, 
                    trajX,
                    trajY,
                    trajZ,
                    trajVertexId,
                    vertexScalars,
                    startFrames,
                    endFrames,
                    durations,
                    VX,
                    VY,
                    surfMin, 
                    surfMax,
                    surfMoy,
                    triangulation->getData()
                    );
  
  if (status != 1)
    return 0;

  // OUTPUT 

  vtkTable *outputTable = vtkTable::GetData(outputVector, 0);
  if (!outputTable){
    this->printErr("output vtkTable");
    return 0;
  }

  vtkUnstructuredGrid *outputGrid = vtkUnstructuredGrid::GetData(outputVector, 1);
  if(!outputGrid) {
    this->printErr("Null output grid.");
    return 0;
  }


  //Colonnes
  
  vtkSmartPointer<vtkIntArray> colStartFrame = vtkSmartPointer<vtkIntArray>::New();
  colStartFrame->SetName("StartFrame");
  colStartFrame->SetNumberOfTuples(numTraj);
  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
    colStartFrame->SetValue(i, startFrames[i]);
  }

  vtkSmartPointer<vtkIntArray> colEndFrame = vtkSmartPointer<vtkIntArray>::New();
  colEndFrame->SetName("EndFrame");
  colEndFrame->SetNumberOfTuples(numTraj);
  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
    colEndFrame->SetValue(i, endFrames[i]);
  }

  vtkSmartPointer<vtkIntArray> colDuration = vtkSmartPointer<vtkIntArray>::New();
  colDuration->SetName("Duration");
  colDuration->SetNumberOfTuples(numTraj);
  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
    colDuration->SetValue(i, durations[i]);
  }

  vtkSmartPointer<vtkDoubleArray> colVX = vtkSmartPointer<vtkDoubleArray>::New();
  colVX->SetName("VX");
  colVX->SetNumberOfTuples(numTraj);
  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
    colVX->SetValue(i, VX[i]);
  }

  vtkSmartPointer<vtkDoubleArray> colVY = vtkSmartPointer<vtkDoubleArray>::New();
  colVY->SetName("VY");
  colVY->SetNumberOfTuples(numTraj);
  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
    colVY->SetValue(i, VY[i]);
  }

  outputTable->AddColumn(colStartFrame);
  outputTable->AddColumn(colEndFrame);
  outputTable->AddColumn(colDuration);
  outputTable->AddColumn(colVX);
  outputTable->AddColumn(colVY);

  // NEW VTU 
  
  vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
  newPoints->SetDataType(inputGrid->GetPoints()->GetDataType());
  outputGrid->SetPoints(newPoints);
  outputGrid->Allocate(numCells);  
 
  std::vector<vtkIdType> oldToNewPointId(inputGrid->GetNumberOfPoints(), -1);
  vtkNew<vtkIdList> cellPointIds;
  vtkSmartPointer<vtkIntArray> newTrajIdArray = vtkSmartPointer<vtkIntArray>::New();
  newTrajIdArray->SetName("NewTrajectoryId");
  newTrajIdArray->SetNumberOfTuples(numCells);

  vtkIdType outCellId = 0;
  size_t trajIndexx = 0;
  for(const auto &trajEntry : groupTraj) {
    // trajEntry.first = ancien ConnectedComponentId, trajEntry.second = liste de cellIds
    for(vtkIdType cellId : trajEntry.second) {
        cellPointIds->Reset();
        inputGrid->GetCellPoints(cellId, cellPointIds);
        vtkIdType n = cellPointIds->GetNumberOfIds();
        std::vector<vtkIdType> newPtIds;
        newPtIds.reserve(n);
        for(vtkIdType i = 0; i < n; ++i) {
            vtkIdType oldPid = cellPointIds->GetId(i);
            // Si le point n'a pas encore été ajouté, on l'ajoute
            if(oldToNewPointId[oldPid] < 0) {
                double coord[3];
                inputGrid->GetPoint(oldPid, coord);
                vtkIdType newPid = newPoints->InsertNextPoint(coord);
                oldToNewPointId[oldPid] = newPid;
            }
        newPtIds.push_back(oldToNewPointId[oldPid]);
        }
        // Ajout de la cellule (même type VTK que l’originale) avec les nouveaux IDs de points
        outputGrid->InsertNextCell(inputGrid->GetCellType(cellId), n, newPtIds.data());
        // Assigner l'ID de trajectoire nouveau à cette cellule
        newTrajIdArray->SetValue(outCellId++, static_cast<int>(trajIndexx));
    }
    ++trajIndexx;
  }

  outputGrid->GetCellData()->AddArray(newTrajIdArray);

  this->printMsg("Fin TrajectoryStatistic");

  return 1;
}
