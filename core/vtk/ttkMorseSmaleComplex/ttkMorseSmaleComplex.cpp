#include                  <ttkMorseSmaleComplex.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkMorseSmaleComplex)

  ttkMorseSmaleComplex::ttkMorseSmaleComplex():
    ScalarField{},
    InputOffsetScalarFieldName{"OutputOffsetScalarField"},
    UseInputOffsetScalarField{},
    IterationThreshold{-1},
    ReverseSaddleMaximumConnection{true},
    ReverseSaddleSaddleConnection{true},
    ComputeAscendingSeparatrices1{true},
    ComputeDescendingSeparatrices1{true},
    ComputeSaddleConnectors{true},
    ComputeAscendingSeparatrices2{true},
    ComputeDescendingSeparatrices2{true},
    ComputeAscendingSegmentation{true},
    ComputeDescendingSegmentation{true},
    ComputeFinalSegmentation{true},
    ScalarFieldId{},
    OffsetFieldId{-1},
    ReturnSaddleConnectors{false},
    SaddleConnectorsPersistenceThreshold{0},

    triangulation_{},
    defaultOffsets_{},
    hasUpdatedMesh_{}
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(4);
}

ttkMorseSmaleComplex::~ttkMorseSmaleComplex(){
  if(defaultOffsets_)
    defaultOffsets_->Delete();
}

int ttkMorseSmaleComplex::FillInputPortInformation(int port, 
    vtkInformation *info){

  switch(port){
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
      break;
  }

  return 1;
}

int ttkMorseSmaleComplex::FillOutputPortInformation(int port, 
    vtkInformation *info){

  switch(port){
    case 0:
    case 1:
    case 2:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;

    case 3:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
      break;
  }

  return 1;
}

int ttkMorseSmaleComplex::setupTriangulation(vtkDataSet* input){
  hasUpdatedMesh_=false;

  triangulation_=ttkTriangulation::getTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_){
    cerr << "[ttkMorseSmaleComplex] Error : ttkTriangulation::getTriangulation() is null." << endl;
    return -1;
  }
#endif

  triangulation_->setWrapper(this);
  // setupTriangulation() is called first to select the correct algorithm (2D or 3D)
  morseSmaleComplex_.setupTriangulation(triangulation_);
  morseSmaleComplex_.setWrapper(this);

  if(triangulation_->isEmpty() or ttkTriangulation::hasChangedConnectivity(triangulation_, input, this)){
    hasUpdatedMesh_=true;
    Modified();
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_->isEmpty()){
    cerr << "[ttkMorseSmaleComplex] Error : ttkTriangulation allocation problem." << endl;
    return -1;
  }
#endif

  return 0;
}

vtkDataArray* ttkMorseSmaleComplex::getScalars(vtkDataSet* input){
  vtkDataArray* inputScalars{};

  vtkPointData* pointData=input->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData){
    cerr << "[ttkMorseSmaleComplex] Error : input has no point data." << endl;
    return inputScalars;
  }
#endif

  if(ScalarField.length()){
    inputScalars=pointData->GetArray(ScalarField.data());
  }
  else{
    inputScalars=pointData->GetArray(ScalarFieldId);
    if(inputScalars)
      ScalarField=inputScalars->GetName();
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars){
    cerr 
      << "[ttkMorseSmaleComplex] Error : input scalar field pointer is null." 
      << endl;
    return inputScalars;
  }
#endif

  return inputScalars;
}

vtkDataArray* ttkMorseSmaleComplex::getOffsets(vtkDataSet* input){
  vtkDataArray* inputOffsets{};

  if(OffsetFieldId != -1){
    inputOffsets = input->GetPointData()->GetArray(OffsetFieldId);
    if(inputOffsets){
      InputOffsetScalarFieldName = inputOffsets->GetName();
      UseInputOffsetScalarField = true;
    }
  }

  if(UseInputOffsetScalarField and InputOffsetScalarFieldName.length()){
    inputOffsets=
      input->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  }
  else{
    if(hasUpdatedMesh_ and defaultOffsets_){
      defaultOffsets_->Delete();
      defaultOffsets_=nullptr;
    }

    if(!defaultOffsets_){
      const int numberOfVertices=input->GetNumberOfPoints();

      defaultOffsets_=vtkIntArray::New();
      defaultOffsets_->SetNumberOfComponents(1);
      defaultOffsets_->SetNumberOfTuples(numberOfVertices);
      defaultOffsets_->SetName("OffsetsScalarField");
      for(int i=0; i<numberOfVertices; ++i)
        defaultOffsets_->SetTuple1(i,i);
    }

    inputOffsets=defaultOffsets_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOffsets){
    cerr 
      << "[ttkMorseSmaleComplex] Error : wrong input offset scalar field." 
      << endl;
    return inputOffsets;
  }
#endif

  return inputOffsets;
}

int ttkMorseSmaleComplex::doIt(vector<vtkDataSet *> &inputs,
    vector<vtkDataSet *> &outputs){
  Memory m;
  int ret{};

  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *outputCriticalPoints = 
    vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  vtkUnstructuredGrid *outputSeparatrices1 = 
    vtkUnstructuredGrid::SafeDownCast(outputs[1]);
  vtkUnstructuredGrid *outputSeparatrices2 = 
    vtkUnstructuredGrid::SafeDownCast(outputs[2]);
  vtkDataSet *outputMorseComplexes = outputs[3];

  ret=setupTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret){
    cerr << "[ttkMorseSmaleComplex] Error : wrong triangulation." << endl;
    return -1;
  }
#endif

  vtkDataArray* inputScalars=getScalars(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars){
    cerr << "[ttkMorseSmaleComplex] Error : wrong scalars." << endl;
    return -1;
  }
#endif

  vtkDataArray* inputOffsets=getOffsets(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOffsets){
    cerr << "[ttkMorseSmaleComplex] Error : wrong offsets." << endl;
    return -1;
  }
#endif

  {
    stringstream msg;
    msg << "[ttkMorseSmaleComplex] Launching computation on field `"
      << inputScalars->GetName() << "'..." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // critical points
  int criticalPoints_numberOfPoints{};
  vector<float> criticalPoints_points;
  vector<int> criticalPoints_points_cellDimensions;
  vector<int> criticalPoints_points_cellIds;
  vector<char> criticalPoints_points_isOnBoundary;
  vector<int> criticalPoints_points_PLVertexIdentifiers;
  vector<int> criticalPoints_points_manifoldSize;

  // 1-separatrices
  int separatrices1_numberOfPoints{};
  vector<float> separatrices1_points;
  vector<char> separatrices1_points_smoothingMask;
  vector<int> separatrices1_points_cellDimensions;
  vector<int> separatrices1_points_cellIds;
  int separatrices1_numberOfCells{};
  vector<int> separatrices1_cells;
  vector<int> separatrices1_cells_sourceIds;
  vector<int> separatrices1_cells_destinationIds;
  vector<int> separatrices1_cells_separatrixIds;
  vector<char> separatrices1_cells_separatrixTypes;
  vector<char> separatrices1_cells_isOnBoundary;

  // 2-separatrices
  int separatrices2_numberOfPoints{};
  vector<float> separatrices2_points;
  int separatrices2_numberOfCells{};
  vector<int> separatrices2_cells;
  vector<int> separatrices2_cells_sourceIds;
  vector<int> separatrices2_cells_separatrixIds;
  vector<char> separatrices2_cells_separatrixTypes;
  vector<char> separatrices2_cells_isOnBoundary;

  const int dimensionality=triangulation_->getCellVertexNumber(0)-1;

  // morse complexes
  const int numberOfVertices=triangulation_->getNumberOfVertices();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfVertices){
    cerr << "[ttkMorseSmaleComplex] Error : input has no vertices." << endl;
    return -1;
  }
#endif

  vtkSmartPointer<vtkIntArray> 
    ascendingManifold=vtkSmartPointer<vtkIntArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!ascendingManifold){
    cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation problem." << 
      endl;
    return -1;
  }
#endif
  ascendingManifold->SetNumberOfComponents(1);
  ascendingManifold->SetNumberOfTuples(numberOfVertices);
  ascendingManifold->SetName("AscendingManifold");

  vtkSmartPointer<vtkIntArray> 
    descendingManifold=vtkSmartPointer<vtkIntArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!descendingManifold){
    cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation problem." << 
      endl;
    return -1;
  }
#endif
  descendingManifold->SetNumberOfComponents(1);
  descendingManifold->SetNumberOfTuples(numberOfVertices);
  descendingManifold->SetName("DescendingManifold");

  vtkSmartPointer<vtkIntArray> 
    morseSmaleManifold=vtkSmartPointer<vtkIntArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!morseSmaleManifold){
    cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation problem." << 
      endl;
    return -1;
  }
#endif
  morseSmaleManifold->SetNumberOfComponents(1);
  morseSmaleManifold->SetNumberOfTuples(numberOfVertices);
  morseSmaleManifold->SetName("MorseSmaleManifold");

  morseSmaleComplex_.setIterationThreshold(IterationThreshold);

  morseSmaleComplex_.setReverseSaddleMaximumConnection(
      ReverseSaddleMaximumConnection);

  morseSmaleComplex_.setReverseSaddleSaddleConnection(
      ReverseSaddleSaddleConnection);

  morseSmaleComplex_.setComputeCriticalPoints(
      ComputeCriticalPoints);

  morseSmaleComplex_.setComputeAscendingSeparatrices1(
      ComputeAscendingSeparatrices1);

  morseSmaleComplex_.setComputeDescendingSeparatrices1(
      ComputeDescendingSeparatrices1);
  morseSmaleComplex_.setComputeSaddleConnectors(ComputeSaddleConnectors);

  morseSmaleComplex_.setComputeAscendingSeparatrices2(
      ComputeAscendingSeparatrices2);

  morseSmaleComplex_.setComputeDescendingSeparatrices2(
      ComputeDescendingSeparatrices2);

  morseSmaleComplex_.setComputeAscendingSegmentation(
      ComputeAscendingSegmentation);

  morseSmaleComplex_.setComputeDescendingSegmentation(
      ComputeDescendingSegmentation);
  morseSmaleComplex_.setComputeFinalSegmentation(ComputeFinalSegmentation);

  morseSmaleComplex_.setReturnSaddleConnectors(
      ReturnSaddleConnectors);
  morseSmaleComplex_.setSaddleConnectorsPersistenceThreshold(
      SaddleConnectorsPersistenceThreshold);

  morseSmaleComplex_.setPrioritizeSpeedOverMemory(
      PrioritizeSpeedOverMemory);

  morseSmaleComplex_.setInputScalarField(inputScalars->GetVoidPointer(0));
  morseSmaleComplex_.setInputOffsets(inputOffsets->GetVoidPointer(0));

  morseSmaleComplex_.setOutputMorseComplexes(
      ascendingManifold->GetVoidPointer(0),
      descendingManifold->GetVoidPointer(0),
      morseSmaleManifold->GetVoidPointer(0));

  switch(inputScalars->GetDataType()){
#ifndef TTK_ENABLE_KAMIKAZE
    vtkTemplateMacro({
      // critical points
      vector<VTK_TT> criticalPoints_points_cellScalars;

      // 1-separatrices
      vector<VTK_TT> separatrices1_cells_separatrixFunctionMaxima;
      vector<VTK_TT> separatrices1_cells_separatrixFunctionMinima;
      vector<VTK_TT> separatrices1_cells_separatrixFunctionDiffs;

      // 2-separatrices
      vector<VTK_TT> separatrices2_cells_separatrixFunctionMaxima;
      vector<VTK_TT> separatrices2_cells_separatrixFunctionMinima;
      vector<VTK_TT> separatrices2_cells_separatrixFunctionDiffs;

      morseSmaleComplex_.setOutputCriticalPoints(
        &criticalPoints_numberOfPoints,
        &criticalPoints_points,
        &criticalPoints_points_cellDimensions,
        &criticalPoints_points_cellIds,
        &criticalPoints_points_cellScalars,
        &criticalPoints_points_isOnBoundary,
        &criticalPoints_points_PLVertexIdentifiers,
        &criticalPoints_points_manifoldSize);


      morseSmaleComplex_.setOutputSeparatrices1(
        &separatrices1_numberOfPoints,
        &separatrices1_points,
        &separatrices1_points_smoothingMask,
        &separatrices1_points_cellDimensions,
        &separatrices1_points_cellIds,
        &separatrices1_numberOfCells,
        &separatrices1_cells,
        &separatrices1_cells_sourceIds,
        &separatrices1_cells_destinationIds,
        &separatrices1_cells_separatrixIds,
        &separatrices1_cells_separatrixTypes,
        &separatrices1_cells_separatrixFunctionMaxima,
        &separatrices1_cells_separatrixFunctionMinima,
        &separatrices1_cells_separatrixFunctionDiffs,
        &separatrices1_cells_isOnBoundary);


      morseSmaleComplex_.setOutputSeparatrices2(
        &separatrices2_numberOfPoints,
        &separatrices2_points,
        &separatrices2_numberOfCells,
        &separatrices2_cells,
        &separatrices2_cells_sourceIds,
        &separatrices2_cells_separatrixIds,
        &separatrices2_cells_separatrixTypes,
        &separatrices2_cells_separatrixFunctionMaxima,
        &separatrices2_cells_separatrixFunctionMinima,
        &separatrices2_cells_separatrixFunctionDiffs,
        &separatrices2_cells_isOnBoundary);

      ret = morseSmaleComplex_.execute<VTK_TT>();
      if (ret) {
        cerr
          << "[ttkMorseSmaleComplex] Error : MorseSmaleComplex.execute() "
          << "error code : " << ret << endl;
        return -1;
      }

      // critical points
      {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        if (!points) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkPoints allocation "
            << "problem." << endl;
          return -1;
        }

        vtkSmartPointer<vtkIntArray> cellDimensions =
          vtkSmartPointer<vtkIntArray>::New();
        if (!cellDimensions) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
            << "problem." << endl;
          return -1;
        }
        cellDimensions->SetNumberOfComponents(1);
        cellDimensions->SetName("CellDimension");

        vtkSmartPointer<vtkIntArray> cellIds =
          vtkSmartPointer<vtkIntArray>::New();
        if (!cellIds) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
            << "problem." << endl;
          return -1;
        }
        cellIds->SetNumberOfComponents(1);
        cellIds->SetName("CellId");

        vtkDataArray* cellScalars = inputScalars->NewInstance();
        if (!cellScalars) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
            << "problem." << endl;
          return -1;
        }
        cellScalars->SetNumberOfComponents(1);
        cellScalars->SetName(ScalarField.data());

        vtkSmartPointer<vtkCharArray> isOnBoundary =
          vtkSmartPointer<vtkCharArray>::New();
        if (!isOnBoundary) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkCharArray allocation "
            << "problem." << endl;
          return -1;
        }
        isOnBoundary->SetNumberOfComponents(1);
        isOnBoundary->SetName("IsOnBoundary");

        vtkSmartPointer<vtkIntArray> PLVertexIdentifiers =
          vtkSmartPointer<vtkIntArray>::New();
        if (!PLVertexIdentifiers) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
            << "problem." << endl;
          return -1;
        }
        PLVertexIdentifiers->SetNumberOfComponents(1);
        PLVertexIdentifiers->SetName("VertexIdentifier");

        vtkSmartPointer<vtkIntArray> manifoldSizeScalars = 
          vtkSmartPointer<vtkIntArray>::New();
        if (!manifoldSizeScalars) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
            << "problem." << endl;
          return -1;
        }
        manifoldSizeScalars->SetNumberOfComponents(1);
        manifoldSizeScalars->SetName("ManifoldSize");

        for (int i = 0; i<criticalPoints_numberOfPoints; ++i) {
          points->InsertNextPoint(criticalPoints_points[3 * i],
            criticalPoints_points[3 * i + 1],
            criticalPoints_points[3 * i + 2]);

          cellDimensions->InsertNextTuple1(
            criticalPoints_points_cellDimensions[i]);
          cellIds->InsertNextTuple1(criticalPoints_points_cellIds[i]);

          cellScalars->InsertNextTuple1(
            criticalPoints_points_cellScalars[i]);

          isOnBoundary->InsertNextTuple1(
            criticalPoints_points_isOnBoundary[i]);

          PLVertexIdentifiers->InsertNextTuple1(
            criticalPoints_points_PLVertexIdentifiers[i]);

          if (ComputeAscendingSegmentation and ComputeDescendingSegmentation)
            manifoldSizeScalars->InsertNextTuple1(
              criticalPoints_points_manifoldSize[i]);
          else
            manifoldSizeScalars->InsertNextTuple1(-1);
        }
        outputCriticalPoints->SetPoints(points);

        vtkPointData* pointData = outputCriticalPoints->GetPointData();
        if (!pointData) {
          cerr << "[ttkMorseSmaleComplex] Error : outputCriticalPoints has "
            << "no point data." << endl;
          return -1;
        }

        pointData->AddArray(cellDimensions);
        pointData->AddArray(cellIds);
        pointData->AddArray(cellScalars);
        pointData->AddArray(isOnBoundary);
        pointData->AddArray(PLVertexIdentifiers);
        pointData->AddArray(manifoldSizeScalars);
      }

      // 1-separatrices
      if (ComputeAscendingSeparatrices1 or ComputeDescendingSeparatrices1 or
        ComputeSaddleConnectors) {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        if (!points) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkPoints allocation "
            << "problem." << endl;
          return -1;
        }
        vtkSmartPointer<vtkCharArray>
          smoothingMask = vtkSmartPointer<vtkCharArray>::New();
        if (!smoothingMask) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkCharArray allocation "
            << "problem." << endl;
          return -1;
        }
        smoothingMask->SetNumberOfComponents(1);
        smoothingMask->SetName("SmoothingMask");

        vtkSmartPointer<vtkIntArray>
          cellDimensions = vtkSmartPointer<vtkIntArray>::New();
        if (!cellDimensions) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
            << "problem." << endl;
          return -1;
        }
        cellDimensions->SetNumberOfComponents(1);
        cellDimensions->SetName("CellDimension");

        vtkSmartPointer<vtkIntArray>
          cellIds = vtkSmartPointer<vtkIntArray>::New();
        if (!cellIds) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
            << "problem." << endl;
          return -1;
        }
        cellIds->SetNumberOfComponents(1);
        cellIds->SetName("CellId");

        vtkSmartPointer<vtkIntArray>
          sourceIds = vtkSmartPointer<vtkIntArray>::New();
        if (!sourceIds) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
            << "problem." << endl;
          return -1;
        }
        sourceIds->SetNumberOfComponents(1);
        sourceIds->SetName("SourceId");

        vtkSmartPointer<vtkIntArray>
          destinationIds = vtkSmartPointer<vtkIntArray>::New();
        if (!destinationIds) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
            << "problem." << endl;
          return -1;
        }
        destinationIds->SetNumberOfComponents(1);
        destinationIds->SetName("DestinationId");

        vtkSmartPointer<vtkIntArray>
          separatrixIds = vtkSmartPointer<vtkIntArray>::New();
        if (!separatrixIds) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
            << "problem." << endl;
          return -1;
        }
        separatrixIds->SetNumberOfComponents(1);
        separatrixIds->SetName("SeparatrixId");

        vtkSmartPointer<vtkCharArray>
          separatrixTypes = vtkSmartPointer<vtkCharArray>::New();
        if (!separatrixTypes) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkCharArray allocation "
            << "problem." << endl;
          return -1;
        }
        separatrixTypes->SetNumberOfComponents(1);
        separatrixTypes->SetName("CriticalPointIndex");

        vtkDataArray* separatrixFunctionMaxima = inputScalars->NewInstance();
        if (!separatrixFunctionMaxima) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
            << "problem." << endl;
          return -1;
        }
        separatrixFunctionMaxima->SetNumberOfComponents(1);
        separatrixFunctionMaxima->SetName("SeparatrixFunctionMaximum");

        vtkDataArray* separatrixFunctionMinima = inputScalars->NewInstance();
        if (!separatrixFunctionMinima) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
            << "problem." << endl;
          return -1;
        }
        separatrixFunctionMinima->SetNumberOfComponents(1);
        separatrixFunctionMinima->SetName("SeparatrixFunctionMinimum");

        vtkDataArray* separatrixFunctionDiffs = inputScalars->NewInstance();
        if (!separatrixFunctionDiffs) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
            << "problem." << endl;
          return -1;
        }
        separatrixFunctionDiffs->SetNumberOfComponents(1);
        separatrixFunctionDiffs->SetName("SeparatrixFunctionDifference");

        vtkSmartPointer<vtkCharArray>
          isOnBoundary = vtkSmartPointer<vtkCharArray>::New();
        if (!isOnBoundary) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkCharArray allocation "
            << "problem." << endl;
          return -1;
        }
        isOnBoundary->SetNumberOfComponents(1);
        isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");

        for (int i = 0; i<separatrices1_numberOfPoints; ++i) {
          points->InsertNextPoint(separatrices1_points[3 * i],
            separatrices1_points[3 * i + 1],
            separatrices1_points[3 * i + 2]);

          
          smoothingMask->InsertNextTuple1(
            separatrices1_points_smoothingMask[i]);
          cellDimensions->InsertNextTuple1(
            separatrices1_points_cellDimensions[i]);
            cellIds->InsertNextTuple1(separatrices1_points_cellIds[i]);
        }
        outputSeparatrices1->SetPoints(points);

        outputSeparatrices1->Allocate(separatrices1_numberOfCells);
        int ptr{};
        for (int i = 0; i<separatrices1_numberOfCells; ++i) {
          vtkIdType line[2];
          line[0] = separatrices1_cells[ptr + 1];
          line[1] = separatrices1_cells[ptr + 2];

          outputSeparatrices1->InsertNextCell(VTK_LINE, 2, line);

          sourceIds->InsertNextTuple1(separatrices1_cells_sourceIds[i]);

          destinationIds->InsertNextTuple1(
            separatrices1_cells_destinationIds[i]);

          separatrixIds->InsertNextTuple1(separatrices1_cells_separatrixIds[i]);

          
          separatrixTypes->InsertNextTuple1(
            separatrices1_cells_separatrixTypes[i]);

          separatrixFunctionMaxima->InsertNextTuple1(
            separatrices1_cells_separatrixFunctionMaxima[i]);

          separatrixFunctionMinima->InsertNextTuple1(
            separatrices1_cells_separatrixFunctionMinima[i]);

          separatrixFunctionDiffs->InsertNextTuple1(
            separatrices1_cells_separatrixFunctionDiffs[i]);

          isOnBoundary->InsertNextTuple1(separatrices1_cells_isOnBoundary[i]);

          ptr += (separatrices1_cells[ptr] + 1);
        }

        vtkPointData* pointData = outputSeparatrices1->GetPointData();
        if (!pointData) {
          cerr << "[ttkMorseSmaleComplex] Error : outputSeparatrices1 has "
            << "no point data." << endl;
          return -1;
        }

        pointData->AddArray(smoothingMask);
        pointData->AddArray(cellDimensions);
        pointData->AddArray(cellIds);

        vtkCellData* cellData = outputSeparatrices1->GetCellData();
        if (!cellData) {
          cerr << "[ttkMorseSmaleComplex] Error : outputSeparatrices1 has "
            << "no cell data." << endl;
          return -1;
        }

        cellData->AddArray(sourceIds);
        cellData->AddArray(destinationIds);
        cellData->AddArray(separatrixIds);
        cellData->AddArray(separatrixTypes);
        cellData->AddArray(separatrixFunctionMaxima);
        cellData->AddArray(separatrixFunctionMinima);
        cellData->AddArray(separatrixFunctionDiffs);
        cellData->AddArray(isOnBoundary);
      }

      // 2-separatrices
      if (dimensionality == 3 and (ComputeAscendingSeparatrices2 or 
        ComputeDescendingSeparatrices2)) {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        if (!points) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkPoints allocation problem."
            << endl;
          return -1;
        }

        vtkSmartPointer<vtkIntArray> sourceIds = 
          vtkSmartPointer<vtkIntArray>::New();
        if (!sourceIds) {
          cerr 
            << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation problem."
            << endl;
          return -1;
        }
        sourceIds->SetNumberOfComponents(1);
        sourceIds->SetName("SourceId");

        vtkSmartPointer<vtkIntArray> separatrixIds = 
    vtkSmartPointer<vtkIntArray>::New();
        if (!separatrixIds) {
          cerr 
            << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation problem."
            << endl;
          return -1;
        }
        separatrixIds->SetNumberOfComponents(1);
        separatrixIds->SetName("SeparatrixId");

        vtkSmartPointer<vtkCharArray> separatrixTypes = 
    vtkSmartPointer<vtkCharArray>::New();
        if (!separatrixTypes) {
          cerr << 
            "[ttkMorseSmaleComplex] Error : vtkCharArray allocation  problem."
            << endl;
          return -1;
        }
        separatrixTypes->SetNumberOfComponents(1);
        separatrixTypes->SetName("CriticalPointIndex");

        vtkDataArray* separatrixFunctionMaxima = inputScalars->NewInstance();
        if (!separatrixFunctionMaxima) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
            << "problem." << endl;
          return -1;
        }
        separatrixFunctionMaxima->SetNumberOfComponents(1);
        separatrixFunctionMaxima->SetName("SeparatrixFunctionMaximum");

        vtkDataArray* separatrixFunctionMinima = inputScalars->NewInstance();
        if (!separatrixFunctionMinima) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
            << "problem." << endl;
          return -1;
        }
        separatrixFunctionMinima->SetNumberOfComponents(1);
        separatrixFunctionMinima->SetName("SeparatrixFunctionMinimum");

        vtkDataArray* separatrixFunctionDiffs = inputScalars->NewInstance();
        if (!separatrixFunctionDiffs) {
          cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
            << "problem." << endl;
          return -1;
        }
        separatrixFunctionDiffs->SetNumberOfComponents(1);
        separatrixFunctionDiffs->SetName("SeparatrixFunctionDifference");

        vtkSmartPointer<vtkCharArray> isOnBoundary = 
          vtkSmartPointer<vtkCharArray>::New();
        if (!isOnBoundary) {
          cerr << 
            "[ttkMorseSmaleComplex] Error : vtkCharArray allocation problem."
            << endl;
          return -1;
        }
        isOnBoundary->SetNumberOfComponents(1);
        isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");

        for (int i = 0; i<separatrices2_numberOfPoints; ++i) {
          points->InsertNextPoint(separatrices2_points[3 * i],
            separatrices2_points[3 * i + 1],
            separatrices2_points[3 * i + 2]);
        }
        outputSeparatrices2->SetPoints(points);

        outputSeparatrices2->Allocate(separatrices2_numberOfCells);
        int ptr{};
        for (int i = 0; i<separatrices2_numberOfCells; ++i) {
          const int vertexNumber = separatrices2_cells[ptr];

          if (vertexNumber == 3) {
            vtkIdType triangle[3];
            triangle[0] = separatrices2_cells[ptr + 1];
            triangle[1] = separatrices2_cells[ptr + 2];
            triangle[2] = separatrices2_cells[ptr + 3];

            outputSeparatrices2->InsertNextCell(
              VTK_TRIANGLE, vertexNumber, triangle);
          }
          else {
            vtkIdType ids[16];
            for (int j = 1; j <= vertexNumber; ++j)
              ids[j - 1] = separatrices2_cells[ptr + j];

            outputSeparatrices2->InsertNextCell(VTK_POLYGON, vertexNumber, ids);
          }

          sourceIds->InsertNextTuple1(separatrices2_cells_sourceIds[i]);
          separatrixIds->InsertNextTuple1(separatrices2_cells_separatrixIds[i]);
          
          separatrixTypes->InsertNextTuple1(
            separatrices2_cells_separatrixTypes[i]);
          separatrixFunctionMaxima->InsertNextTuple1(
            separatrices2_cells_separatrixFunctionMaxima[i]);

          separatrixFunctionMinima->InsertNextTuple1(
            separatrices2_cells_separatrixFunctionMinima[i]);

          separatrixFunctionDiffs->InsertNextTuple1(
            separatrices2_cells_separatrixFunctionDiffs[i]);
          isOnBoundary->InsertNextTuple1(separatrices2_cells_isOnBoundary[i]);

          ptr += (separatrices2_cells[ptr] + 1);
        }

        vtkCellData* cellData = outputSeparatrices2->GetCellData();
        if (!cellData) {
          cerr << "[ttkMorseSmaleComplex] Error : "
            << "outputSeparatrices2 has no cell data." << endl;
          return -1;
        }

        cellData->AddArray(sourceIds);
        cellData->AddArray(separatrixIds);
        cellData->AddArray(separatrixTypes);
        cellData->AddArray(separatrixFunctionMaxima);
        cellData->AddArray(separatrixFunctionMinima);
        cellData->AddArray(separatrixFunctionDiffs);
        cellData->AddArray(isOnBoundary);
      }
    });
  #else
    vtkTemplateMacro({
      // critical points
      vector<VTK_TT> criticalPoints_points_cellScalars;

      // 1-separatrices
      vector<VTK_TT> separatrices1_cells_separatrixFunctionMaxima;
      vector<VTK_TT> separatrices1_cells_separatrixFunctionMinima;
      vector<VTK_TT> separatrices1_cells_separatrixFunctionDiffs;

      // 2-separatrices
      vector<VTK_TT> separatrices2_cells_separatrixFunctionMaxima;
      vector<VTK_TT> separatrices2_cells_separatrixFunctionMinima;
      vector<VTK_TT> separatrices2_cells_separatrixFunctionDiffs;

      morseSmaleComplex_.setOutputCriticalPoints(
        &criticalPoints_numberOfPoints,
        &criticalPoints_points,
        &criticalPoints_points_cellDimensions,
        &criticalPoints_points_cellIds,
        &criticalPoints_points_cellScalars,
        &criticalPoints_points_isOnBoundary,
        &criticalPoints_points_PLVertexIdentifiers,
        &criticalPoints_points_manifoldSize);


      morseSmaleComplex_.setOutputSeparatrices1(
        &separatrices1_numberOfPoints,
        &separatrices1_points,
        &separatrices1_points_smoothingMask,
        &separatrices1_points_cellDimensions,
        &separatrices1_points_cellIds,
        &separatrices1_numberOfCells,
        &separatrices1_cells,
        &separatrices1_cells_sourceIds,
        &separatrices1_cells_destinationIds,
        &separatrices1_cells_separatrixIds,
        &separatrices1_cells_separatrixTypes,
        &separatrices1_cells_separatrixFunctionMaxima,
        &separatrices1_cells_separatrixFunctionMinima,
        &separatrices1_cells_separatrixFunctionDiffs,
        &separatrices1_cells_isOnBoundary);


      morseSmaleComplex_.setOutputSeparatrices2(
        &separatrices2_numberOfPoints,
        &separatrices2_points,
        &separatrices2_numberOfCells,
        &separatrices2_cells,
        &separatrices2_cells_sourceIds,
        &separatrices2_cells_separatrixIds,
        &separatrices2_cells_separatrixTypes,
        &separatrices2_cells_separatrixFunctionMaxima,
        &separatrices2_cells_separatrixFunctionMinima,
        &separatrices2_cells_separatrixFunctionDiffs,
        &separatrices2_cells_isOnBoundary);

      morseSmaleComplex_.execute<VTK_TT>();

      // critical points
      {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        vtkSmartPointer<vtkIntArray> cellDimensions =
          vtkSmartPointer<vtkIntArray>::New();

        cellDimensions->SetNumberOfComponents(1);
        cellDimensions->SetName("CellDimension");

        vtkSmartPointer<vtkIntArray> cellIds =
          vtkSmartPointer<vtkIntArray>::New();

        cellIds->SetNumberOfComponents(1);
        cellIds->SetName("CellId");

        vtkDataArray* cellScalars = inputScalars->NewInstance();

        cellScalars->SetNumberOfComponents(1);
        cellScalars->SetName(ScalarField.data());

        vtkSmartPointer<vtkCharArray> isOnBoundary =
          vtkSmartPointer<vtkCharArray>::New();

        isOnBoundary->SetNumberOfComponents(1);
        isOnBoundary->SetName("IsOnBoundary");

        vtkSmartPointer<vtkIntArray> PLVertexIdentifiers =
          vtkSmartPointer<vtkIntArray>::New();

        PLVertexIdentifiers->SetNumberOfComponents(1);
        PLVertexIdentifiers->SetName("VertexIdentifier");

        vtkSmartPointer<vtkIntArray> manifoldSizeScalars = 
          vtkSmartPointer<vtkIntArray>::New();

        manifoldSizeScalars->SetNumberOfComponents(1);
        manifoldSizeScalars->SetName("ManifoldSize");

        for (int i = 0; i<criticalPoints_numberOfPoints; ++i) {
          points->InsertNextPoint(criticalPoints_points[3 * i],
            criticalPoints_points[3 * i + 1],
            criticalPoints_points[3 * i + 2]);

          cellDimensions->InsertNextTuple1(
            criticalPoints_points_cellDimensions[i]);
          cellIds->InsertNextTuple1(criticalPoints_points_cellIds[i]);

          cellScalars->InsertNextTuple1(
            criticalPoints_points_cellScalars[i]);

          isOnBoundary->InsertNextTuple1(
            criticalPoints_points_isOnBoundary[i]);

          PLVertexIdentifiers->InsertNextTuple1(
            criticalPoints_points_PLVertexIdentifiers[i]);

          if (ComputeAscendingSegmentation and ComputeDescendingSegmentation)
            manifoldSizeScalars->InsertNextTuple1(
              criticalPoints_points_manifoldSize[i]);
          else
            manifoldSizeScalars->InsertNextTuple1(-1);
        }
        outputCriticalPoints->SetPoints(points);

        vtkPointData* pointData = outputCriticalPoints->GetPointData();

        pointData->AddArray(cellDimensions);
        pointData->AddArray(cellIds);
        pointData->AddArray(cellScalars);
        pointData->AddArray(isOnBoundary);
        pointData->AddArray(PLVertexIdentifiers);
        pointData->AddArray(manifoldSizeScalars);
      }

      // 1-separatrices
      if (ComputeAscendingSeparatrices1 or ComputeDescendingSeparatrices1 or
        ComputeSaddleConnectors) {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        vtkSmartPointer<vtkCharArray>
          smoothingMask = vtkSmartPointer<vtkCharArray>::New();

        smoothingMask->SetNumberOfComponents(1);
        smoothingMask->SetName("SmoothingMask");

        vtkSmartPointer<vtkIntArray>
          cellDimensions = vtkSmartPointer<vtkIntArray>::New();

        cellDimensions->SetNumberOfComponents(1);
        cellDimensions->SetName("CellDimension");

        vtkSmartPointer<vtkIntArray>
          cellIds = vtkSmartPointer<vtkIntArray>::New();

        cellIds->SetNumberOfComponents(1);
        cellIds->SetName("CellId");

        vtkSmartPointer<vtkIntArray>
          sourceIds = vtkSmartPointer<vtkIntArray>::New();

        sourceIds->SetNumberOfComponents(1);
        sourceIds->SetName("SourceId");

        vtkSmartPointer<vtkIntArray>
          destinationIds = vtkSmartPointer<vtkIntArray>::New();

        destinationIds->SetNumberOfComponents(1);
        destinationIds->SetName("DestinationId");

        vtkSmartPointer<vtkIntArray>
          separatrixIds = vtkSmartPointer<vtkIntArray>::New();

        separatrixIds->SetNumberOfComponents(1);
        separatrixIds->SetName("SeparatrixId");

        vtkSmartPointer<vtkCharArray>
          separatrixTypes = vtkSmartPointer<vtkCharArray>::New();

        separatrixTypes->SetNumberOfComponents(1);
        separatrixTypes->SetName("CriticalPointIndex");

        vtkDataArray* separatrixFunctionMaxima = inputScalars->NewInstance();

        separatrixFunctionMaxima->SetNumberOfComponents(1);
        separatrixFunctionMaxima->SetName("SeparatrixFunctionMaximum");

        vtkDataArray* separatrixFunctionMinima = inputScalars->NewInstance();

        separatrixFunctionMinima->SetNumberOfComponents(1);
        separatrixFunctionMinima->SetName("SeparatrixFunctionMinimum");

        vtkDataArray* separatrixFunctionDiffs = inputScalars->NewInstance();

        separatrixFunctionDiffs->SetNumberOfComponents(1);
        separatrixFunctionDiffs->SetName("SeparatrixFunctionDifference");

        vtkSmartPointer<vtkCharArray>
          isOnBoundary = vtkSmartPointer<vtkCharArray>::New();

        isOnBoundary->SetNumberOfComponents(1);
        isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");

        for (int i = 0; i<separatrices1_numberOfPoints; ++i) {
          points->InsertNextPoint(separatrices1_points[3 * i],
            separatrices1_points[3 * i + 1],
            separatrices1_points[3 * i + 2]);

          
          smoothingMask->InsertNextTuple1(
            separatrices1_points_smoothingMask[i]);
          
          cellDimensions->InsertNextTuple1(
            separatrices1_points_cellDimensions[i]);
          cellIds->InsertNextTuple1(separatrices1_points_cellIds[i]);
        }
        outputSeparatrices1->SetPoints(points);

        outputSeparatrices1->Allocate(separatrices1_numberOfCells);
        int ptr{};
        for (int i = 0; i<separatrices1_numberOfCells; ++i) {
          vtkIdType line[2];
          line[0] = separatrices1_cells[ptr + 1];
          line[1] = separatrices1_cells[ptr + 2];

          outputSeparatrices1->InsertNextCell(VTK_LINE, 2, line);

          sourceIds->InsertNextTuple1(separatrices1_cells_sourceIds[i]);

          destinationIds->InsertNextTuple1(
            separatrices1_cells_destinationIds[i]);

          separatrixIds->InsertNextTuple1(separatrices1_cells_separatrixIds[i]);

          separatrixTypes->InsertNextTuple1(
            separatrices1_cells_separatrixTypes[i]);

          separatrixFunctionMaxima->InsertNextTuple1(
            separatrices1_cells_separatrixFunctionMaxima[i]);

          separatrixFunctionMinima->InsertNextTuple1(
            separatrices1_cells_separatrixFunctionMinima[i]);

          separatrixFunctionDiffs->InsertNextTuple1(
            separatrices1_cells_separatrixFunctionDiffs[i]);

          isOnBoundary->InsertNextTuple1(separatrices1_cells_isOnBoundary[i]);

          ptr += (separatrices1_cells[ptr] + 1);
        }

        vtkPointData* pointData = outputSeparatrices1->GetPointData();

        pointData->AddArray(smoothingMask);
        pointData->AddArray(cellDimensions);
        pointData->AddArray(cellIds);

        vtkCellData* cellData = outputSeparatrices1->GetCellData();

        cellData->AddArray(sourceIds);
        cellData->AddArray(destinationIds);
        cellData->AddArray(separatrixIds);
        cellData->AddArray(separatrixTypes);
        cellData->AddArray(separatrixFunctionMaxima);
        cellData->AddArray(separatrixFunctionMinima);
        cellData->AddArray(separatrixFunctionDiffs);
        cellData->AddArray(isOnBoundary);
      }

      // 2-separatrices
      if (dimensionality == 3 and (ComputeAscendingSeparatrices2 or 
        ComputeDescendingSeparatrices2)) {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        vtkSmartPointer<vtkIntArray> sourceIds = 
          vtkSmartPointer<vtkIntArray>::New();

        sourceIds->SetNumberOfComponents(1);
        sourceIds->SetName("SourceId");

        vtkSmartPointer<vtkIntArray> separatrixIds = 
          vtkSmartPointer<vtkIntArray>::New();

        separatrixIds->SetNumberOfComponents(1);
        separatrixIds->SetName("SeparatrixId");

        vtkSmartPointer<vtkCharArray> separatrixTypes = 
          vtkSmartPointer<vtkCharArray>::New();

        separatrixTypes->SetNumberOfComponents(1);
        separatrixTypes->SetName("CriticalPointIndex");

        vtkDataArray* separatrixFunctionMaxima = inputScalars->NewInstance();

        separatrixFunctionMaxima->SetNumberOfComponents(1);
        separatrixFunctionMaxima->SetName("SeparatrixFunctionMaximum");

        vtkDataArray* separatrixFunctionMinima = inputScalars->NewInstance();

        separatrixFunctionMinima->SetNumberOfComponents(1);
        separatrixFunctionMinima->SetName("SeparatrixFunctionMinimum");

        vtkDataArray* separatrixFunctionDiffs = inputScalars->NewInstance();

        separatrixFunctionDiffs->SetNumberOfComponents(1);
        separatrixFunctionDiffs->SetName("SeparatrixFunctionDifference");

        vtkSmartPointer<vtkCharArray> isOnBoundary = 
          vtkSmartPointer<vtkCharArray>::New();

        isOnBoundary->SetNumberOfComponents(1);
        isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");

        for (int i = 0; i<separatrices2_numberOfPoints; ++i) {
          points->InsertNextPoint(separatrices2_points[3 * i],
            separatrices2_points[3 * i + 1],
            separatrices2_points[3 * i + 2]);
        }
        outputSeparatrices2->SetPoints(points);

        outputSeparatrices2->Allocate(separatrices2_numberOfCells);
        int ptr{};
        for (int i = 0; i<separatrices2_numberOfCells; ++i) {
          const int vertexNumber = separatrices2_cells[ptr];

          if (vertexNumber == 3) {
            vtkIdType triangle[3];
            triangle[0] = separatrices2_cells[ptr + 1];
            triangle[1] = separatrices2_cells[ptr + 2];
            triangle[2] = separatrices2_cells[ptr + 3];

            outputSeparatrices2->InsertNextCell(
              VTK_TRIANGLE, vertexNumber, triangle);
          }
          else {
            vtkIdType ids[16];
            for (int j = 1; j <= vertexNumber; ++j)
              ids[j - 1] = separatrices2_cells[ptr + j];

            outputSeparatrices2->InsertNextCell(VTK_POLYGON, vertexNumber, ids);
          }

          sourceIds->InsertNextTuple1(separatrices2_cells_sourceIds[i]);
          separatrixIds->InsertNextTuple1(separatrices2_cells_separatrixIds[i]);
          
          separatrixTypes->InsertNextTuple1(
            separatrices2_cells_separatrixTypes[i]);
          separatrixFunctionMaxima->InsertNextTuple1(
            separatrices2_cells_separatrixFunctionMaxima[i]);

          separatrixFunctionMinima->InsertNextTuple1(
            separatrices2_cells_separatrixFunctionMinima[i]);

          separatrixFunctionDiffs->InsertNextTuple1(
            separatrices2_cells_separatrixFunctionDiffs[i]);
          isOnBoundary->InsertNextTuple1(separatrices2_cells_isOnBoundary[i]);

          ptr += (separatrices2_cells[ptr] + 1);
        }

        vtkCellData* cellData = outputSeparatrices2->GetCellData();

        cellData->AddArray(sourceIds);
        cellData->AddArray(separatrixIds);
        cellData->AddArray(separatrixTypes);
        cellData->AddArray(separatrixFunctionMaxima);
        cellData->AddArray(separatrixFunctionMinima);
        cellData->AddArray(separatrixFunctionDiffs);
        cellData->AddArray(isOnBoundary);
      }
    });
#endif
  }

  // morse complexes
  if(ComputeAscendingSegmentation or ComputeDescendingSegmentation){
    outputMorseComplexes->ShallowCopy(input);

    vtkPointData* pointData=outputMorseComplexes->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData){
      cerr 
        << "[ttkMorseSmaleComplex] Error : outputMorseComplexes has no point "
        << "data." << endl;
      return -1;
    }
#endif
    pointData->AddArray(descendingManifold);
    pointData->AddArray(ascendingManifold);
    pointData->AddArray(morseSmaleManifold);
  }

  {
    stringstream msg;
    msg << "[ttkMorseSmaleComplex] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
