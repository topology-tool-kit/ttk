#include <ttkTrajectoryStatistics.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCharArray.h>
#include <vtkTable.h> 
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkGradientFilter.h>
#include <vtkLine.h>

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
  this->SetNumberOfOutputPorts(5);
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


  if (port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 1);
    return 1;
  }
  
  if (port == 3) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }

  if (port == 4) {
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

/* #######################################################
 * 
 *                      TRAJ DATA 
 *
 * #######################################################
*/
  vtkIntArray *compIdArray = vtkIntArray::SafeDownCast(
            inputGrid->GetCellData()->GetArray("ConnectedComponentId"));
  vtkIntArray *timeArray = vtkIntArray::SafeDownCast(
            inputGrid->GetPointData()->GetArray("TimeStep"));
  vtkIntArray *vertexGlobalIdArray = vtkIntArray::SafeDownCast(
            inputGrid->GetPointData()->GetArray("VertexGlobalId"));
  vtkIntArray *compLength = vtkIntArray::SafeDownCast(
            inputGrid->GetCellData()->GetArray("ComponentLength"));
  vtkDoubleArray *persistanceArray = vtkDoubleArray::SafeDownCast(
            inputGrid->GetPointData()->GetArray("InstantPersistence"));

  if (!compIdArray || !timeArray || !vertexGlobalIdArray || !compLength){
    this->printErr("Missing data in input vtu");
    return 0;
  }

  vtkIdType numCells = inputGrid->GetNumberOfCells();

  if (numCells == 0)
    return 1;

  // map avec : {id de la traj, [cellsID, ...]}
  std::map<int, std::vector<vtkIdType>> groupTraj;

  for (vtkIdType cellId = 0; cellId < numCells; ++cellId){
    int trajId = compIdArray->GetValue(cellId); // ConnectedComponentId de la traj
    groupTraj[trajId].push_back(cellId);
  }

  const size_t numTraj = groupTraj.size();

  std::vector<std::vector<int>>    trajTime(numTraj); // TimeStep
  std::vector<std::vector<double>> trajX(numTraj), trajY(numTraj), instantPersistance(numTraj); 
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

    trajVertexId[trajIndex].reserve(numPoints); 

    for (vtkIdType pointId : pointsIds){
        int t = timeArray->GetValue(pointId);
        double persistance = persistanceArray->GetValue(pointId); 
        int globalId = vertexGlobalIdArray->GetValue(pointId);
        trajTime[trajIndex].push_back(t);
        trajVertexId[trajIndex].push_back(globalId);
        instantPersistance[trajIndex].push_back(persistance);
        double coords[3];
        inputGrid->GetPoint(pointId, coords);
        trajX[trajIndex].push_back(coords[0]);
        trajY[trajIndex].push_back(coords[1]);
    }
    ++trajIndex;
  }


/* #######################################################
 * 
 *               SCALAR DATASET  
 *
 * #######################################################
*/
  std::vector<vtkDataArray *> inputScalarFieldsRaw;
  std::vector<vtkDataArray *> inputScalarFields;
  const auto pointData= inputDataSet->GetPointData();

  if (!pointData){
    this->printErr("scalarArray missing");
    return 0;
  }

  int numberOfInputFields = pointData->GetNumberOfArrays();
  
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

  for(int i = 0; i < numberOfInputFields ; i++) {
    vtkDataArray *currentScalarField = inputScalarFieldsRaw[i];
    // Print scalar field names:
    // std::cout << currentScalarField->GetName() << std::endl;
    inputScalarFields.push_back(currentScalarField);
  }

  int const fieldNumber = inputScalarFields.size();
  std::vector<void *> inputFields(fieldNumber);
  for(int i = 0; i < fieldNumber; i++) {
    inputFields[i] = ttkUtils::GetVoidPointer(inputScalarFields[i]);
  }
  this->setInputScalars(inputFields);
  this->setFiltreX(filtreX);
  this->setFiltreY(filtreY);
  this->setCosCol(cosCol);
  this->setMaxRadus(maxRadus);
  this->setMinFrameDist(minFrameDist);
  this->setMaxFrameDist(maxFrameDist);
  this->setSpatialScale(1/(spatialScale*1000));
  this->setInterFrame(interFrame*std::pow(10.0,-6.0));
  this->setConvertDur(convertDur);
  this->setMinVx(minVx);
  this->setCoordCratere(coordCratere);
  this->setThreshCratereAngle(threshCratereAngle);
  this->setMaxX(maxX);
  this->setMaxY(maxY);
  this->setMinY(minY);
  this->setMinX(minX);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);
  
  const int nPts = inputScalarFields[0]->GetNumberOfTuples();

/* #######################################################
 * 
 *               OUTPUT VECTOR   
 *
 * #######################################################
*/

  //############## Correct trajectory ####################

  std::vector<std::vector<double>> newTraj(numTraj);
  std::vector<std::vector<double>> finalTraj;
  std::vector<FuseRecord> fuseRecords;


  this->correctTrajectory(trajTime, trajX, trajY, finalTraj,newTraj, fuseRecords); 

  std::vector<int>  durations(numTraj);
  std::vector<double> VX(numTraj), VY(numTraj), 
                      surfMin(numTraj), surfMax(numTraj), surfMean(numTraj);
  std::vector<std::vector<ttk::SimplexId>> allVertexDebris(numTraj);
  std::vector<ttk::SimplexId> excludedCriticalPoints;

  std::vector<std::vector<double>> gradientNorms;
  if(!computeAllGradientMagnitudes(inputDataSet,
                                   inputScalarFieldsRaw,
                                   gradientNorms)) {
    this->printErr("Impossible de calculer gradient magnitudes.");
    return 0;
  }

  int status = 0;

  ttkVtkTemplateMacro(inputScalarFields[0]->GetDataType(), triangulation->getType(),
      (status = this->execute<VTK_TT, TTK_TT>(
                        trajTime, 
                        trajVertexId,
                        durations,
                        VX,
                        VY,
                        surfMin, 
                        surfMax,
                        surfMean,
                        allVertexDebris,
                        excludedCriticalPoints,
                        frameSurface,
                        errSurf,
                        gradientNorms,
                        finalTraj,
                        (TTK_TT *)triangulation->getData()
                        )));
  
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

  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 2);
  if(!outputDataSet) {
    this->printErr("Null output grid2.");
    return 0;
  }

  vtkUnstructuredGrid *outputTraj = vtkUnstructuredGrid::GetData(outputVector, 3);
  if(!outputTraj) {
    this->printErr("Null output gridTraj.");
    return 0;
  }

  vtkUnstructuredGrid *outputMerge = vtkUnstructuredGrid::GetData(outputVector, 4);
  if(!outputMerge) {
    this->printErr("Null output gridMerge.");
    return 0;
  }

  //Colonnes

  const int numMerge = finalTraj.size();
  
  vtkSmartPointer<vtkIntArray> colStartFrame = vtkSmartPointer<vtkIntArray>::New();
  colStartFrame->SetName("StartFrame");
  colStartFrame->SetNumberOfTuples(numMerge);
  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numMerge); ++i) {
    colStartFrame->SetValue(i, finalTraj[i][4]);
  }

  vtkSmartPointer<vtkIntArray> colEndFrame = vtkSmartPointer<vtkIntArray>::New();
  colEndFrame->SetName("EndFrame");
  colEndFrame->SetNumberOfTuples(numMerge);
  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numMerge); ++i) {
    colEndFrame->SetValue(i, finalTraj[i][5]);
  }

  vtkSmartPointer<vtkIntArray> colDuration = vtkSmartPointer<vtkIntArray>::New();
  colDuration->SetName("Duration");
  colDuration->SetNumberOfTuples(numMerge);
  double dur_mean = 0;
  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numMerge); ++i) {
    colDuration->SetValue(i, durations[i]);
    dur_mean += durations[i];
  }
  dur_mean = dur_mean/numMerge;
  this->printMsg("DUREE MOYENNE = "+std::to_string(dur_mean));

  vtkSmartPointer<vtkDoubleArray> colVX = vtkSmartPointer<vtkDoubleArray>::New();
  colVX->SetName("VX");
  colVX->SetNumberOfTuples(numMerge);
  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numMerge); ++i) {
    colVX->SetValue(i, VX[i]);
  }

  vtkSmartPointer<vtkDoubleArray> colVY = vtkSmartPointer<vtkDoubleArray>::New();
  colVY->SetName("VY");
  colVY->SetNumberOfTuples(numMerge);
  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numMerge); ++i) {
    colVY->SetValue(i, VY[i]);
  }


  std::vector<std::vector<double>> surfaceFinal(finalTraj.size());
  std::vector<double> surfFinalMean(finalTraj.size());
  for (int i=0; i<numTraj; i++){
    if (newTraj[i][4] != -1)
     surfaceFinal[newTraj[i][4]].push_back(surfMean[i]);
  }
  for (int i=0; i<finalTraj.size(); i++){
    double sum = 0;
    for (int j=0;j<surfaceFinal[i].size(); j++){
        sum += surfaceFinal[i][j];
    }
    surfFinalMean[i] = sum/surfaceFinal[i].size();
  }

/*  vtkSmartPointer<vtkDoubleArray> colSurfMin = vtkSmartPointer<vtkDoubleArray>::New();
  colSurfMin->SetName("SurfaceMin");
  colSurfMin->SetNumberOfTuples(numTraj);
  for(vtkIdType i = 0; i < numTraj; ++i) {
    colSurfMin->SetValue(i, surfMin[i]);
  }
  vtkSmartPointer<vtkDoubleArray> colSurfMax = vtkSmartPointer<vtkDoubleArray>::New();
  colSurfMax->SetName("SurfaceMax");
  colSurfMax->SetNumberOfTuples(numTraj);
  for(vtkIdType i = 0; i < numTraj; ++i) {
    colSurfMax->SetValue(i, surfMax[i]);
  }
*/
  vtkSmartPointer<vtkDoubleArray> colSurfMean = vtkSmartPointer<vtkDoubleArray>::New();
  colSurfMean->SetName("SurfaceMean");
  colSurfMean->SetNumberOfTuples(finalTraj.size());
  double surf_mean = 0;
  for(vtkIdType i = 0; i < finalTraj.size(); ++i) {
    colSurfMean->SetValue(i, surfFinalMean[i]);
    surf_mean += surfFinalMean[i];
  }
  surf_mean = surf_mean/numMerge;
  this->printMsg("SURFACE MOYENNE =" +std::to_string(surf_mean));
  
  outputTable->AddColumn(colStartFrame);
  outputTable->AddColumn(colEndFrame);
  outputTable->AddColumn(colDuration);
  outputTable->AddColumn(colVX);
  outputTable->AddColumn(colVY);
  //outputTable->AddColumn(colSurfMin);
  //outputTable->AddColumn(colSurfMax);
  outputTable->AddColumn(colSurfMean);


  vtkSmartPointer<vtkDataSet> surfOutput = vtkSmartPointer<vtkDataSet>::Take(inputDataSet->NewInstance());
  surfOutput->ShallowCopy(inputDataSet);
  outputDataSet->ShallowCopy(surfOutput);

  vtkSmartPointer<vtkCharArray> surfaceVertexArray = vtkSmartPointer<vtkCharArray>::New();
  surfaceVertexArray->SetName("SurfaceVertex");

  vtkSmartPointer<vtkIntArray> surfaceTrajId = vtkSmartPointer<vtkIntArray>::New();
  surfaceTrajId->SetName("FinalTrajId");


  vtkIdType numPoints = inputDataSet->GetNumberOfPoints();
  surfaceVertexArray->SetNumberOfTuples(numPoints);
  surfaceVertexArray->FillComponent(0, 0); // tous les points à 0 (false)

  surfaceTrajId->SetNumberOfTuples(numPoints);
  surfaceTrajId->FillComponent(0,-1);

  for (int i=0; i<allVertexDebris.size(); i++){
      std::vector<ttk::SimplexId> &trajSurface = allVertexDebris[i];
      //if (allVertexDebris[i].size() != 0)
        //this->printMsg("surface trouvé -> " + std::to_string(allVertexDebris[i].size()) + " traj Id = " + std::to_string(newTraj[i][4]));
      int finalId = -1;
      if (newTraj[i][4] != -1)
        finalId = newTraj[i][4]; 
      auto [minIt, maxIt] = std::minmax_element(instantPersistance[i].begin(), instantPersistance[i].end());
      double persisMean
          = std::accumulate(
              instantPersistance[i].begin(),
              instantPersistance[i].end(),
              0.0
            )
          / instantPersistance[i].size();
      //if (allVertexDebris[i].size() != 0)
       // this->printMsg("persis = min -> " + std::to_string(*minIt) + " max -> " + std::to_string(*maxIt) + " mean -> " + std::to_string(persisMean));
      for(const auto vertexId : trajSurface) {
        if(vertexId >= 0 && vertexId < numPoints){
            surfaceTrajId->SetValue(vertexId, finalId);
            if (surfaceVertexArray->GetValue(vertexId) == 0) 
                surfaceVertexArray->SetValue(vertexId, 1); // true
            else
                surfaceVertexArray->SetValue(vertexId, 2); // doublon

        }
      }
  }

  for (const auto vertexId: excludedCriticalPoints){
    surfaceVertexArray->SetValue(vertexId, 3); // diverge        
  } 

  outputDataSet->GetPointData()->AddArray(surfaceVertexArray);
  outputDataSet->GetPointData()->AddArray(surfaceTrajId);
  
    
  vtkSmartPointer<vtkPoints> linearPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> linearLines  = vtkSmartPointer<vtkCellArray>::New();
  linearPoints->SetNumberOfPoints(2*newTraj.size());

  vtkSmartPointer<vtkIntArray> linearFinalId = vtkSmartPointer<vtkIntArray>::New();
  linearFinalId->SetName("linearFinalId");
  linearFinalId->SetNumberOfTuples(newTraj.size());
  
  vtkSmartPointer<vtkIntArray> trajId = vtkSmartPointer<vtkIntArray>::New();
  trajId->SetName("trajId");
  trajId->SetNumberOfTuples(newTraj.size());

  for (int i=0; i<newTraj.size(); i++){
    double x,y;

    const auto &coef = newTraj[i];
    const int startFrame = trajTime[i].front(); 
    const int endFrame = trajTime[i].back();
    
    x =  coef[0]*startFrame + coef[2];
    y =  coef[1]*startFrame + coef[3];
    linearPoints->SetPoint(2*i, x, y, startFrame);
    
    x =  coef[0]*endFrame + coef[2];
    y =  coef[1]*endFrame + coef[3];
    linearPoints->SetPoint(2*i+1, x, y, endFrame);

    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, 2*i + 0);
    line->GetPointIds()->SetId(1, 2*i + 1);
    linearLines->InsertNextCell(line);
    linearFinalId->SetValue(i , coef[4]);
    trajId->SetValue(i,i);
  }

  outputGrid->SetPoints(linearPoints);
  outputGrid->SetCells(VTK_LINE, linearLines);
  outputGrid->GetCellData()->AddArray(linearFinalId);
  outputGrid->GetCellData()->AddArray(trajId);

  vtkSmartPointer<vtkPoints> fusionPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> fusionLines  = vtkSmartPointer<vtkCellArray>::New();
  fusionPoints->SetNumberOfPoints(2*fuseRecords.size());
 
  vtkSmartPointer<vtkIntArray> fusionFinalId = vtkSmartPointer<vtkIntArray>::New();
  fusionFinalId->SetName("fusionFinalId");
  fusionFinalId->SetNumberOfTuples(fuseRecords.size()); 

  // merge 
  int count =0;
  for (FuseRecord &f : fuseRecords){
    double x,y;

    const auto &coef = newTraj[f.i];
    const int startFrame = f.startFrame; 
    const int endFrame = f.endFrame;
    
    x =  coef[0]*endFrame + coef[2];
    y =  coef[1]*endFrame + coef[3];
    fusionPoints->SetPoint(2*count, x, y, endFrame);
    const auto &coef2 = newTraj[f.j];
    x =  coef2[0]*startFrame + coef2[2];
    y =  coef2[1]*startFrame + coef2[3];
    fusionPoints->SetPoint(2*count+1, x, y, startFrame);

    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, 2*count + 0);
    line->GetPointIds()->SetId(1, 2*count + 1);
    fusionLines->InsertNextCell(line);
    fusionFinalId->SetValue(count , coef[4]);
    count++;
  }

  outputTraj->SetPoints(fusionPoints);
  outputTraj->SetCells(VTK_LINE, fusionLines);
  outputTraj->GetCellData()->AddArray(fusionFinalId);

  vtkSmartPointer<vtkPoints> mergePoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> mergeLines  = vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkIntArray> newMergeIdArray = vtkSmartPointer<vtkIntArray>::New();
  newMergeIdArray->SetName("MergeId");
  newMergeIdArray->SetNumberOfTuples(numMerge);

  //finaltraj
  mergePoints->SetNumberOfPoints(2*finalTraj.size());
  for (int i=0; i<finalTraj.size(); i++){
    double x, y; 
    double startFrame; 
    const auto &coef = finalTraj[i];
    if (!extendTraj)
        startFrame = finalTraj[i][4];
    else 
        startFrame = 0;
    const double endFrame = finalTraj[i][5];
    
    x =  coef[0]*startFrame + coef[2];
    y =  coef[1]*startFrame + coef[3];
    mergePoints->SetPoint(2*i,x,y, startFrame);
    
    x =  coef[0]*endFrame+ coef[2];
    y =  coef[1]*endFrame + coef[3];
    mergePoints->SetPoint(2*i+1, x, y, endFrame);

    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, 2*i + 0);
    line->GetPointIds()->SetId(1, 2*i + 1);
    mergeLines->InsertNextCell(line);
    newMergeIdArray->SetValue(i , i);

  }

  outputMerge->SetPoints(mergePoints);
  outputMerge->SetCells(VTK_LINE, mergeLines);
  outputMerge->GetCellData()->AddArray(newMergeIdArray);


  this->printMsg("Fin TrajectoryStatistic");

  return 1;
}

int ttkTrajectoryStatistics::computeAllGradientMagnitudes(
  vtkDataSet *inputDataSet,
  const std::vector<vtkDataArray *> &inputScalarFields,
  std::vector<std::vector<double>> &gradientNorms
) {
  if(!inputDataSet) {
    this->printErr("inputDataSet is nullptr.");
    return 0;
  }

  const size_t nFields = inputScalarFields.size();
  if(nFields == 0) {
    this->printErr("pas de scalar fields fournis.");
    return 0;
  }

  vtkIdType nPts = inputDataSet->GetNumberOfPoints();
  if(nPts <= 0) {
    this->printErr("maillage vide (nPts <= 0).");
    return 0;
  }

  gradientNorms.clear();
  gradientNorms.resize(nFields);
  for(size_t f = 0; f < nFields; ++f) {
    gradientNorms[f].assign(static_cast<size_t>(nPts), 0.0);
  }

  for(size_t f = 0; f < nFields; f++) {
    vtkDataArray *currScalar = inputScalarFields[f];
    if(!currScalar || !currScalar->GetName()) {
      this->printErr("scalar array invalide en frame " + std::to_string(f));
      return 0;
    }
    const char *scalarName = currScalar->GetName();

    vtkSmartPointer<vtkGradientFilter> gradFilter = vtkSmartPointer<vtkGradientFilter>::New();
    gradFilter->SetInputData(inputDataSet);
    gradFilter->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);


    gradFilter->Update();

    vtkDataSet *gradOutput = gradFilter->GetOutput();
    if (!gradOutput){
        this->printErr("grad output missing");
    }

    vtkDataArray *gradArray = gradOutput->GetPointData()->GetArray("Gradients");
    if(!gradArray) {
      this->printErr(
            "recup gradArray echec"
      );
      return 0;
    }

    for(vtkIdType pid = 0; pid < nPts; ++pid) {
      double gx = gradArray->GetComponent(pid, 0);
      double gy = gradArray->GetComponent(pid, 1);
      double gz = gradArray->GetComponent(pid, 2);
      gradientNorms[f][static_cast<size_t>(pid)] = std::sqrt(gx * gx + gy * gy + gz * gz);
    }
  }

  return 1;
}

