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
  this->SetNumberOfOutputPorts(4);
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

  this->printMsg("Traj Id Array created");

  const size_t numTraj = groupTraj.size();

  std::vector<std::vector<int>>    trajTime(numTraj); // TimeStep
  std::vector<std::vector<double>> trajX(numTraj), trajY(numTraj), trajZ(numTraj); 
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
        int globalId = vertexGlobalIdArray->GetValue(pointId);
        trajTime[trajIndex].push_back(t);
        trajVertexId[trajIndex].push_back(globalId);
        double coords[3];
        inputGrid->GetPoint(pointId, coords);
        trajX[trajIndex].push_back(coords[0]);
        trajY[trajIndex].push_back(coords[1]);
        trajZ[trajIndex].push_back(coords[2]);
    }
    ++trajIndex;
  }

  this->printMsg("Pre-fetching done, "+ std::to_string(trajIndex) + " unique trajectory find");



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


  std::vector<std::vector<double>> gradientNorms;
  if(!computeAllGradientMagnitudes(inputDataSet,
                                   inputScalarFieldsRaw,
                                   gradientNorms)) {
    this->printErr("Impossible de calculer gradient magnitudes.");
    return 0;
  }

  numberOfInputFields = inputScalarFieldsRaw.size();
  this->printMsg("New number of frame = " + std::to_string(numberOfInputFields));


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

  this->printMsg("Scalars recup");




  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);
  
  int surfMethods = 1;
  if (gradThresh && !addThresh)
    surfMethods = 2;
  else if (addThresh && !gradThresh)
    surfMethods = 1;
  else {
    this->printErr("0 ou + que 2 méthodes choisis pour les surfaces");
    return 0;
  }
   
  const int nPts = inputScalarFields[0]->GetNumberOfTuples();
  this->printMsg("nPts = " + std::to_string(nPts));


/* #######################################################
 * 
 *               OUTPUT VECTOR   
 *
 * #######################################################
*/

  std::vector<int> startFrames(numTraj), endFrames(numTraj), durations(numTraj);
  std::vector<double> VX(numTraj), VY(numTraj), 
                      surfMin(numTraj), surfMax(numTraj), surfMean(numTraj);
  std::vector<std::vector<ttk::SimplexId>> allVertexDebris(numTraj);
  std::vector<ttk::SimplexId> excludedCriticalPoints;
  std::vector<std::vector<double>> newTraj(numTraj);

  int status = 0;

  ttkVtkTemplateMacro(inputScalarFields[0]->GetDataType(), triangulation->getType(),
      (status = this->execute<VTK_TT, TTK_TT>(
                        trajTime, 
                        trajX,
                        trajY,
                        trajZ,
                        trajVertexId,
                        startFrames,
                        endFrames,
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
                        surfMethods,
                        newTraj,
                        gradientNorms,
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
  if(!outputGrid) {
    this->printErr("Null output grid2.");
    return 0;
  }

  vtkUnstructuredGrid *outputTraj = vtkUnstructuredGrid::GetData(outputVector, 3);
  if(!outputTraj) {
    this->printErr("Null output gridTraj.");
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

  vtkSmartPointer<vtkDoubleArray> colSurfMin = vtkSmartPointer<vtkDoubleArray>::New();
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
  vtkSmartPointer<vtkDoubleArray> colSurfMean = vtkSmartPointer<vtkDoubleArray>::New();
  colSurfMean->SetName("SurfaceMean");
  colSurfMean->SetNumberOfTuples(numTraj);
  for(vtkIdType i = 0; i < numTraj; ++i) {
    colSurfMean->SetValue(i, surfMean[i]);
  }
  
  outputTable->AddColumn(colStartFrame);
  outputTable->AddColumn(colEndFrame);
  outputTable->AddColumn(colDuration);
  outputTable->AddColumn(colVX);
  outputTable->AddColumn(colVY);
  outputTable->AddColumn(colSurfMin);
  outputTable->AddColumn(colSurfMax);
  outputTable->AddColumn(colSurfMean);



  
  outputGrid->ShallowCopy(inputGrid);
  
  vtkSmartPointer<vtkIntArray> newTrajIdArray = vtkSmartPointer<vtkIntArray>::New();
  newTrajIdArray->SetName("NewTrajectoryId");
  newTrajIdArray->SetNumberOfTuples(inputGrid->GetNumberOfCells());

  vtkIdType outCellId = 0;
  size_t trajIndexx = 0;
  for(const auto &trajEntry : groupTraj) {
    for(vtkIdType cellId : trajEntry.second) {
      newTrajIdArray->SetValue(cellId, static_cast<int>(trajIndexx));
    }
    ++trajIndexx;
  }

  outputGrid->GetCellData()->AddArray(newTrajIdArray); 


  vtkSmartPointer<vtkDataSet> copy = vtkSmartPointer<vtkDataSet>::Take(inputDataSet->NewInstance());
  copy->ShallowCopy(inputDataSet); 
  outputDataSet->ShallowCopy(copy);


  vtkIdType numPoints = inputDataSet->GetNumberOfPoints();

  vtkSmartPointer<vtkIntArray> surfaceVertexArray = vtkSmartPointer<vtkIntArray>::New();
  surfaceVertexArray->SetName("SurfaceVertex");

  surfaceVertexArray->SetNumberOfTuples(numPoints);
  surfaceVertexArray->FillComponent(0, 0); // tous les points à 0 (false)

  vtkSmartPointer<vtkDoubleArray> gradArray = vtkSmartPointer<vtkDoubleArray>::New();
  gradArray->SetName("gradArray");

  gradArray->SetNumberOfTuples(numPoints);
  gradArray->FillComponent(0, 0); 
  
  vtkSmartPointer<vtkDoubleArray> trajV = vtkSmartPointer<vtkDoubleArray>::New();
  trajV->SetName("trajV");

  trajV->SetNumberOfTuples(numPoints);
  trajV->FillComponent(0, 0);


  for (size_t i = 0; i < excludedCriticalPoints.size(); i++){
    surfaceVertexArray->SetValue(excludedCriticalPoints[i], 2);
  }
  for(size_t trajId = 0; trajId < allVertexDebris.size(); ++trajId) {
    const auto &trajSurfaces = allVertexDebris[trajId];
    for(const auto vertexId : trajSurfaces) {
        if(vertexId >= 0 && vertexId < numPoints){
          int doublon = surfaceVertexArray->GetValue(vertexId);
          if (doublon == 1){
            surfaceVertexArray->SetValue(vertexId, 3);
            trajV->SetValue(vertexId, trajId);
          }
          else  {
            surfaceVertexArray->SetValue(vertexId, 1); // true
            trajV->SetValue(vertexId, trajId);
          }
        }
    }
  }

  for (int i =0; i<numPoints; i++){
    gradArray->SetValue(i, gradientNorms[0][i]);
  }
  outputDataSet->GetPointData()->AddArray(surfaceVertexArray);
  outputDataSet->GetPointData()->AddArray(gradArray);
  outputDataSet->GetPointData()->AddArray(trajV);




  auto nTraj = static_cast<vtkIdType>(newTraj.size());
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines  = vtkSmartPointer<vtkCellArray>::New();

  points->SetNumberOfPoints(2*nTraj);
  vtkSmartPointer<vtkIntArray> TrajIdArray = vtkSmartPointer<vtkIntArray>::New();
  TrajIdArray->SetName("NewTrajectoryId");
  TrajIdArray->SetNumberOfTuples(nTraj);
  

  for (vtkIdType i = 0; i<nTraj; ++i){
 
    double t0 = static_cast<double>(startFrames[i]);
    double t1 = static_cast<double>(endFrames[i]);

    const auto &coef = newTraj[i];
    double x0 = coef[0] * t0 + coef[2];
    double y0 = coef[1] * t0 + coef[3];
    double z0 = t0 *10;
    points->SetPoint(2*i + 0, x0, y0, z0);
    double x1 = coef[0] * t1 + coef[2];
    double y1 = coef[1] * t1 + coef[3];
    double z1 = t1 * 10;
    points->SetPoint(2*i + 1, x1, y1, z1);

    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, 2*i + 0);
    line->GetPointIds()->SetId(1, 2*i + 1);
    lines->InsertNextCell(line);


    TrajIdArray->SetValue(i, i);
  }

  outputTraj->SetPoints(points);
  outputTraj->SetCells(VTK_LINE, lines);
  outputTraj->GetCellData()->AddArray(TrajIdArray);







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
  this->printMsg("nFields = " + std::to_string(nFields));
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

  // Pour chaque frame (= chacun des inputScalarFields[f])
  for(size_t f = 0; f < nFields; f++) {
    vtkDataArray *currScalar = inputScalarFields[f];
    if(!currScalar || !currScalar->GetName()) {
      this->printErr("scalar array invalide en frame " + std::to_string(f));
      return 0;
    }
    const char *scalarName = currScalar->GetName();

    // 1) Instanciation de vtkGradientFilter
    vtkSmartPointer<vtkGradientFilter> gradFilter = vtkSmartPointer<vtkGradientFilter>::New();
    gradFilter->SetInputData(inputDataSet);
    gradFilter->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);


    gradFilter->Update();

    // 3) Récupérer la sortie contenant l’array vectoriel "Gradients_<scalarName>"
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

    // 4) Calculer la magnitude au niveau de chaque sommet<
    for(vtkIdType pid = 0; pid < nPts; ++pid) {
      double gx = gradArray->GetComponent(pid, 0);
      double gy = gradArray->GetComponent(pid, 1);
      double gz = gradArray->GetComponent(pid, 2);
      gradientNorms[f][static_cast<size_t>(pid)] = std::sqrt(gx * gx + gy * gy + gz * gz);
    }
  }

  return 1;
}

