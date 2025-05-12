#include <ttkTrajectoryStatistics.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkTable.h> 
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>

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
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}


int ttkTrajectoryStatistics::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}


int ttkTrajectoryStatistics::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
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


  vtkIntArray *compIdArray = vtkIntArray::SafeDownCast(
            inputGrid->GetCellData()->GetArray("ConnectedComponentId"));
  vtkIntArray *timeArray = vtkIntArray::SafeDownCast(
            inputGrid->GetPointData()->GetArray("TimeStep"));

  if (!compIdArray){
    this->printErr("ConnectrdComponentId missing");
    return 0;
  }

  if (!timeArray){
    this->printErr(" TimeStep missing");
    return 0;
  }

  this->printMsg("Extraction ConnectedComponentId && TimeStep done");

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

  std::vector<std::vector<int>> trajTime(numTraj); // TimeStep
  std::vector<std::vector<double>> trajX(numTraj), trajY(numTraj), trajZ(numTraj); // coord

  vtkNew<vtkIdList> cellPointsIds;
  size_t trajIndex = 0;
  for (const auto &trajEntry : groupTraj){
    
    const std::vector<vtkIdType> &cellIds = trajEntry.second;
    std::vector<vtkIdType> pointsIds; 
    pointsIds.reserve(cellIds.size()*2);

    this->printMsg("Allocation pointIds done");

    // recup tt les points de la traj
    for (vtkIdType cellId : cellIds){
        cellPointsIds->Reset();
        this->printMsg("Getting cell : "+ std::to_string(cellId));
        inputGrid->GetCellPoints(cellId, cellPointsIds);
        this->printMsg("Got it");
        vtkIdType numPts = cellPointsIds->GetNumberOfIds();
        for (vtkIdType i=0; i<numPts; ++i){
            this->printMsg("pushing it");
            pointsIds.push_back(cellPointsIds->GetId(i));
        }   
    }

    this->printMsg("All Points in pointIds (with double)");

    std::sort(pointsIds.begin(), pointsIds.end());
    pointsIds.erase(std::unique(pointsIds.begin(), pointsIds.end()), 
                    pointsIds.end()); // doublon
    std::sort(pointsIds.begin(), pointsIds.end(), [&](vtkIdType a, vtkIdType b){
        return timeArray->GetValue(a) < timeArray->GetValue(b);
    }); // tri TimeStep

    this->printMsg("All unique points extrac from traj : " + std::to_string(trajIndex));

    const size_t numPoints = pointsIds.size();
    trajTime[trajIndex].reserve(numPoints);
    trajX[trajIndex].reserve(numPoints);
    trajY[trajIndex].reserve(numPoints);
    trajZ[trajIndex].reserve(numPoints);

    for (vtkIdType pointId : pointsIds){
        int t = timeArray->GetValue(pointId);
        double coords[3];
        inputGrid->GetPoint(pointId, coords);
        trajTime[trajIndex].push_back(t);
        trajX[trajIndex].push_back(coords[0]);
        trajY[trajIndex].push_back(coords[1]);
        trajZ[trajIndex].push_back(coords[2]);
    }
    ++trajIndex;
  }

  //Appel core/base
  std::vector<int> startFrames(numTraj), endFrames(numTraj), durations(numTraj);
  std::vector<double> lengths(numTraj);
  std::vector<double> moyVelX(numTraj), moyVelY(numTraj);

  //int status = this->execute(...);
  

  // Construction vtkTable
  
  vtkTable *outputTable = vtkTable::GetData(outputVector, 0);
  if (!outputTable){
    this->printErr("output vtkTable");
    return 0;
  }

  //Colonnes
  
  //vtkSmartPointer<vtkIntArray> colStartFrame = vtkSmartPointer<vtkIntArray>::New();
  //colStartFrame->SetName("StartFrame");
  //colStartFrame->SetNumerOfTuples(numTraj);
  //for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
  //  colStartFrame->SetValue(i, startFrames[i]);
  //}

  //vtkSmartPointer<vtkIntArray> colEndFrame = vtkSmartPointer<vtkIntArray>::New();
  //colEndFrame->SetName("EndFrame");
  //colEndFrame->SetNumberOfTuples(numTraj);
  //for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
  //  colEndFrame->SetValue(i, endFrames[i]);
  //}

  //vtkSmartPointer<vtkIntArray> colDuration = vtkSmartPointer<vtkIntArray>::New();
  //colDuration->SetName("Duration");
  //colDuration->SetNumberOfTuples(numTraj);
  //for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
  //  colDuration->SetValue(i, durations[i]);
  //}

  //vtkSmartPointer<vtkDoubleArray> colLength = 
  //                          vtkSmartPointer<vtkDoubleArray>::New();
  //colLength->SetName("Length");
  //colLength->SetNumberOfTuples(numTraj);
  //for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
  //  colLength->SetValue(i, lengths[i]);
  //}

  //vtkSmartPointer<vtkDoubleArray> colMoyVelX = 
  //                   vtkSmartPointer<vtkDoubleArray>::New();
  //colMoyVelX->SetName("moyVelX");
  //colMoyVelX->SetNumberOfTuples(numTraj);
  //for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
  //  colMoyVelX->SetValue(i, moyVelX[i]);
  //}
  
  //vtkSmartPointer<vtkDoubleArray> colMoyVelY = 
  //                vtkSmartPointer<vtkDoubleArray>::New();
  //colMoyVelY->SetName("MeanVelY");
  //colMoyVelY->SetNumberOfTuples(numTraj);
  //for(vtkIdType i = 0; i < static_cast<vtkIdType>(numTraj); ++i) {
  //  colMoyVelY->SetValue(i, meanVelY[i]);
  //}
  
  //outputTable->AddColumn(colStartFrame);
  //outputTable->AddColumn(colEndFrame);
  //outputTable->AddColumn(colDuration);
  //outputTable->AddColumn(colLength);
  //outputTable->AddColumn(colMoyVelX);
  //outputTable->AddColumn(colMoyVelY);
  return 1;
}
