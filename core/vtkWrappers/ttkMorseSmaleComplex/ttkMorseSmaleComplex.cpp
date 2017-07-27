#include                  <ttkMorseSmaleComplex.h>

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
#ifndef withKamikaze
  if(!triangulation_){
    cerr << "[vtkContourForests] Error : ttkTriangulation::getTriangulation() is null." << endl;
    return -1;
  }
#endif

  triangulation_->setWrapper(this);
  morseSmaleComplex_.setWrapper(this);
  morseSmaleComplex_.setupTriangulation(triangulation_);

  if(triangulation_->isEmpty() or ttkTriangulation::hasChangedConnectivity(triangulation_, input, this)){
    hasUpdatedMesh_=true;
    Modified();
  }

#ifndef withKamikaze
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

#ifndef withKamikaze
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

#ifndef withKamikaze
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

#ifndef withKamikaze
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
#ifndef withKamikaze
  if(ret){
    cerr << "[ttkMorseSmaleComplex] Error : wrong triangulation." << endl;
    return -1;
  }
#endif

  vtkDataArray* inputScalars=getScalars(input);
#ifndef withKamikaze
  if(!inputScalars){
    cerr << "[ttkMorseSmaleComplex] Error : wrong scalars." << endl;
    return -2;
  }
#endif

  vtkDataArray* inputOffsets=getOffsets(input);
#ifndef withKamikaze
  if(!inputOffsets){
    cerr << "[ttkMorseSmaleComplex] Error : wrong offsets." << endl;
    return -3;
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
#ifndef withKamikaze
  if(!numberOfVertices){
    cerr << "[ttkMorseSmaleComplex] Error : input has no vertices." << endl;
    return -4;
  }
#endif

  vtkSmartPointer<vtkIntArray> 
    ascendingManifold=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
  if(!ascendingManifold){
    cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation problem." << 
      endl;
    return -5;
  }
#endif
  ascendingManifold->SetNumberOfComponents(1);
  ascendingManifold->SetNumberOfTuples(numberOfVertices);
  ascendingManifold->SetName("AscendingManifold");

  vtkSmartPointer<vtkIntArray> 
    descendingManifold=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
  if(!descendingManifold){
    cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation problem." << 
      endl;
    return -6;
  }
#endif
  descendingManifold->SetNumberOfComponents(1);
  descendingManifold->SetNumberOfTuples(numberOfVertices);
  descendingManifold->SetName("DescendingManifold");

  vtkSmartPointer<vtkIntArray> 
    morseSmaleManifold=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
  if(!morseSmaleManifold){
    cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation problem." << 
      endl;
    return -7;
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

  morseSmaleComplex_.setInputScalarField(inputScalars->GetVoidPointer(0));
  morseSmaleComplex_.setInputOffsets(inputOffsets->GetVoidPointer(0));

  morseSmaleComplex_.setOutputMorseComplexes(
      ascendingManifold->GetVoidPointer(0),
      descendingManifold->GetVoidPointer(0),
      morseSmaleManifold->GetVoidPointer(0));

  switch(inputScalars->GetDataType()){
    vtkTemplateMacro(({
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

          ret=morseSmaleComplex_.execute<VTK_TT>();
#ifndef withKamikaze
          if(ret){
            cerr 
              << "[ttkMorseSmaleComplex] Error : MorseSmaleComplex.execute() "
              << "error code : " << ret << endl;
            return -8;
          }
#endif

          // critical points
          {
            vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
#ifndef withKamikaze
            if(!points){
              cerr << "[ttkMorseSmaleComplex] Error : vtkPoints allocation "
                << "problem." << endl;
              return -9;
            }
#endif

            vtkSmartPointer<vtkIntArray> cellDimensions = 
              vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!cellDimensions){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
                << "problem." << endl;
              return -10;
            }
#endif
            cellDimensions->SetNumberOfComponents(1);
            cellDimensions->SetName("CellDimension");

            vtkSmartPointer<vtkIntArray> cellIds=
              vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!cellIds){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
                << "problem." << endl;
              return -11;
            }
#endif
            cellIds->SetNumberOfComponents(1);
            cellIds->SetName("CellId");

            vtkDataArray* cellScalars=inputScalars->NewInstance();
#ifndef withKamikaze
            if(!cellScalars){
              cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
                << "problem." << endl;
              return -12;
            }
#endif
            cellScalars->SetNumberOfComponents(1);
            cellScalars->SetName(ScalarField.data());

            vtkSmartPointer<vtkCharArray> isOnBoundary=
              vtkSmartPointer<vtkCharArray>::New();
#ifndef withKamikaze
            if(!isOnBoundary){
              cerr << "[ttkMorseSmaleComplex] Error : vtkCharArray allocation "
                << "problem." << endl;
              return -13;
            }
#endif
            isOnBoundary->SetNumberOfComponents(1);
            isOnBoundary->SetName("IsOnBoundary");

            vtkSmartPointer<vtkIntArray> PLVertexIdentifiers=
              vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!PLVertexIdentifiers){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
                << "problem." << endl;
              return -14;
            }
#endif
            PLVertexIdentifiers->SetNumberOfComponents(1);
            PLVertexIdentifiers->SetName("VertexIdentifier");

            vtkSmartPointer<vtkIntArray> manifoldSizeScalars=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!manifoldSizeScalars){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
                << "problem." << endl;
              return -15;
            }
#endif
            manifoldSizeScalars->SetNumberOfComponents(1);
            manifoldSizeScalars->SetName("ManifoldSize");

            for(int i=0; i<criticalPoints_numberOfPoints; ++i){
              points->InsertNextPoint(criticalPoints_points[3*i],
                  criticalPoints_points[3*i+1],
                  criticalPoints_points[3*i+2]);

              cellDimensions->InsertNextTuple1(
                  criticalPoints_points_cellDimensions[i]);
              cellIds->InsertNextTuple1(criticalPoints_points_cellIds[i]);

              cellScalars->InsertNextTuple1(
                  criticalPoints_points_cellScalars[i]);

              isOnBoundary->InsertNextTuple1(
                  criticalPoints_points_isOnBoundary[i]);

              PLVertexIdentifiers->InsertNextTuple1(
                  criticalPoints_points_PLVertexIdentifiers[i]);

              if(ComputeAscendingSegmentation and ComputeDescendingSegmentation)
                manifoldSizeScalars->InsertNextTuple1(
                    criticalPoints_points_manifoldSize[i]);
              else
                manifoldSizeScalars->InsertNextTuple1(-1);
            }
            outputCriticalPoints->SetPoints(points);

            vtkPointData* pointData=outputCriticalPoints->GetPointData();
#ifndef withKamikaze
            if(!pointData){
              cerr << "[ttkMorseSmaleComplex] Error : outputCriticalPoints has "
                << "no point data." << endl;
              return -16;
            }
#endif

            pointData->AddArray(cellDimensions);
            pointData->AddArray(cellIds);
            pointData->AddArray(cellScalars);
            pointData->AddArray(isOnBoundary);
            pointData->AddArray(PLVertexIdentifiers);
            pointData->AddArray(manifoldSizeScalars);
          }

          // 1-separatrices
          if(ComputeAscendingSeparatrices1 or ComputeDescendingSeparatrices1 or 
              ComputeSaddleConnectors){
            vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
#ifndef withKamikaze
            if(!points){
              cerr << "[ttkMorseSmaleComplex] Error : vtkPoints allocation "
                << "problem." << endl;
              return -17;
            }
#endif
            vtkSmartPointer<vtkCharArray>
              smoothingMask=vtkSmartPointer<vtkCharArray>::New();
#ifndef withKamikaze
            if(!smoothingMask){
              cerr << "[ttkMorseSmaleComplex] Error : vtkCharArray allocation "
                << "problem." << endl;
              return -18;
            }
#endif
            smoothingMask->SetNumberOfComponents(1);
            smoothingMask->SetName("SmoothingMask");

            vtkSmartPointer<vtkIntArray> 
              cellDimensions=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!cellDimensions){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
                << "problem." << endl;
              return -19;
            }
#endif
            cellDimensions->SetNumberOfComponents(1);
            cellDimensions->SetName("CellDimension");

            vtkSmartPointer<vtkIntArray> 
              cellIds=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!cellIds){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
                << "problem." << endl;
              return -20;
            }
#endif
            cellIds->SetNumberOfComponents(1);
            cellIds->SetName("CellId");

            vtkSmartPointer<vtkIntArray> 
              sourceIds=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!sourceIds){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
                << "problem." << endl;
              return -21;
            }
#endif
            sourceIds->SetNumberOfComponents(1);
            sourceIds->SetName("SourceId");

            vtkSmartPointer<vtkIntArray> 
              destinationIds=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!destinationIds){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
                << "problem." << endl;
              return -22;
            }
#endif
            destinationIds->SetNumberOfComponents(1);
            destinationIds->SetName("DestinationId");

            vtkSmartPointer<vtkIntArray> 
              separatrixIds=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!separatrixIds){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
                << "problem." << endl;
              return -23;
            }
#endif
            separatrixIds->SetNumberOfComponents(1);
            separatrixIds->SetName("SeparatrixId");

            vtkSmartPointer<vtkCharArray> 
              separatrixTypes=vtkSmartPointer<vtkCharArray>::New();
#ifndef withKamikaze
            if(!separatrixTypes){
              cerr << "[ttkMorseSmaleComplex] Error : vtkCharArray allocation "
                << "problem." << endl;
              return -24;
            }
#endif
            separatrixTypes->SetNumberOfComponents(1);
            separatrixTypes->SetName("CriticalPointIndex");

            vtkDataArray* separatrixFunctionMaxima=inputScalars->NewInstance();
#ifndef withKamikaze
            if(!separatrixFunctionMaxima){
              cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
                << "problem." << endl;
              return -25;
            }
#endif
            separatrixFunctionMaxima->SetNumberOfComponents(1);
            separatrixFunctionMaxima->SetName("SeparatrixFunctionMaximum");

            vtkDataArray* separatrixFunctionMinima=inputScalars->NewInstance();
#ifndef withKamikaze
            if(!separatrixFunctionMinima){
              cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
                << "problem." << endl;
              return -26;
            }
#endif
            separatrixFunctionMinima->SetNumberOfComponents(1);
            separatrixFunctionMinima->SetName("SeparatrixFunctionMinimum");

            vtkDataArray* separatrixFunctionDiffs=inputScalars->NewInstance();
#ifndef withKamikaze
            if(!separatrixFunctionDiffs){
              cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
                << "problem." << endl;
              return -27;
            }
#endif
            separatrixFunctionDiffs->SetNumberOfComponents(1);
            separatrixFunctionDiffs->SetName("SeparatrixFunctionDifference");

            vtkSmartPointer<vtkCharArray> 
              isOnBoundary=vtkSmartPointer<vtkCharArray>::New();
#ifndef withKamikaze
            if(!isOnBoundary){
              cerr << "[ttkMorseSmaleComplex] Error : vtkCharArray allocation "
                << "problem." << endl;
              return -28;
            }
#endif
            isOnBoundary->SetNumberOfComponents(1);
            isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");

            for(int i=0; i<separatrices1_numberOfPoints; ++i){
              points->InsertNextPoint(separatrices1_points[3*i],
                  separatrices1_points[3*i+1],
                  separatrices1_points[3*i+2]);

              smoothingMask->InsertNextTuple1(separatrices1_points_smoothingMask[i]);
              cellDimensions->InsertNextTuple1(separatrices1_points_cellDimensions[i]);
              cellIds->InsertNextTuple1(separatrices1_points_cellIds[i]);
            }
            outputSeparatrices1->SetPoints(points);

            outputSeparatrices1->Allocate(separatrices1_numberOfCells);
            int ptr{};
            for(int i=0; i<separatrices1_numberOfCells; ++i){
              vtkIdType line[2];
              line[0]=separatrices1_cells[ptr+1];
              line[1]=separatrices1_cells[ptr+2];

              outputSeparatrices1->InsertNextCell(VTK_LINE, 2, line);

              sourceIds->InsertNextTuple1(separatrices1_cells_sourceIds[i]);

              destinationIds->InsertNextTuple1(separatrices1_cells_destinationIds[i]);

              separatrixIds->InsertNextTuple1(separatrices1_cells_separatrixIds[i]);

              separatrixTypes->InsertNextTuple1(separatrices1_cells_separatrixTypes[i]);

              separatrixFunctionMaxima->InsertNextTuple1(
                  separatrices1_cells_separatrixFunctionMaxima[i]);

              separatrixFunctionMinima->InsertNextTuple1(
                  separatrices1_cells_separatrixFunctionMinima[i]);

              separatrixFunctionDiffs->InsertNextTuple1(
                  separatrices1_cells_separatrixFunctionDiffs[i]);

              isOnBoundary->InsertNextTuple1(separatrices1_cells_isOnBoundary[i]);

              ptr+=(separatrices1_cells[ptr]+1);
            }

            vtkPointData* pointData=outputSeparatrices1->GetPointData();
#ifndef withKamikaze
            if(!pointData){
              cerr << "[ttkMorseSmaleComplex] Error : outputSeparatrices1 has "
                << "no point data." << endl;
              return -29;
            }
#endif

            pointData->AddArray(smoothingMask);
            pointData->AddArray(cellDimensions);
            pointData->AddArray(cellIds);

            vtkCellData* cellData=outputSeparatrices1->GetCellData();
#ifndef withKamikaze
            if(!cellData){
              cerr << "[ttkMorseSmaleComplex] Error : outputSeparatrices1 has "
                << "no cell data." << endl;
              return -30;
            }
#endif

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
          if(dimensionality==3 and (ComputeAscendingSeparatrices2 or ComputeDescendingSeparatrices2)){
            vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
#ifndef withKamikaze
            if(!points){
              cerr << "[ttkMorseSmaleComplex] Error : vtkPoints allocation problem." 
                << endl;
              return -32;
            }
#endif

            vtkSmartPointer<vtkIntArray> sourceIds=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!sourceIds){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation problem." 
                << endl;
              return -32;
            }
#endif
            sourceIds->SetNumberOfComponents(1);
            sourceIds->SetName("SourceId");

            vtkSmartPointer<vtkIntArray> separatrixIds=vtkSmartPointer<vtkIntArray>::New();
#ifndef withKamikaze
            if(!separatrixIds){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation problem." 
                << endl;
              return -33;
            }
#endif
            separatrixIds->SetNumberOfComponents(1);
            separatrixIds->SetName("SeparatrixId");

            vtkSmartPointer<vtkCharArray> separatrixTypes=vtkSmartPointer<vtkCharArray>::New();
#ifndef withKamikaze
            if(!separatrixTypes){
              cerr << "[ttkMorseSmaleComplex] Error : vtkCharArray allocation problem." 
                << endl;
              return -34;
            }
#endif
            separatrixTypes->SetNumberOfComponents(1);
            separatrixTypes->SetName("CriticalPointIndex");

            vtkDataArray* separatrixFunctionMaxima=inputScalars->NewInstance();
#ifndef withKamikaze
            if(!separatrixFunctionMaxima){
              cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
                << "problem." << endl;
              return -35;
            }
#endif
            separatrixFunctionMaxima->SetNumberOfComponents(1);
            separatrixFunctionMaxima->SetName("SeparatrixFunctionMaximum");

            vtkDataArray* separatrixFunctionMinima=inputScalars->NewInstance();
#ifndef withKamikaze
            if(!separatrixFunctionMinima){
              cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
                << "problem." << endl;
              return -36;
            }
#endif
            separatrixFunctionMinima->SetNumberOfComponents(1);
            separatrixFunctionMinima->SetName("SeparatrixFunctionMinimum");

            vtkDataArray* separatrixFunctionDiffs=inputScalars->NewInstance();
#ifndef withKamikaze
            if(!separatrixFunctionDiffs){
              cerr << "[ttkMorseSmaleComplex] Error : vtkDataArray allocation "
                << "problem." << endl;
              return -37;
            }
#endif
            separatrixFunctionDiffs->SetNumberOfComponents(1);
            separatrixFunctionDiffs->SetName("SeparatrixFunctionDifference");

            vtkSmartPointer<vtkCharArray> isOnBoundary=vtkSmartPointer<vtkCharArray>::New();
#ifndef withKamikaze
            if(!isOnBoundary){
              cerr << "[ttkMorseSmaleComplex] Error : vtkCharArray allocation problem." 
                << endl;
              return -38;
            }
#endif
            isOnBoundary->SetNumberOfComponents(1);
            isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");

            for(int i=0; i<separatrices2_numberOfPoints; ++i){
              points->InsertNextPoint(separatrices2_points[3*i],
                  separatrices2_points[3*i+1],
                  separatrices2_points[3*i+2]);
            }
            outputSeparatrices2->SetPoints(points);

            outputSeparatrices2->Allocate(separatrices2_numberOfCells);
            int ptr{};
            for(int i=0; i<separatrices2_numberOfCells; ++i){
              const int vertexNumber=separatrices2_cells[ptr];

              if(vertexNumber==3){
                vtkIdType triangle[3];
                triangle[0]=separatrices2_cells[ptr+1];
                triangle[1]=separatrices2_cells[ptr+2];
                triangle[2]=separatrices2_cells[ptr+3];

                outputSeparatrices2->InsertNextCell(VTK_TRIANGLE, vertexNumber, triangle);
              }
              else{
                vtkIdType ids[16];
                for(int j=1; j<=vertexNumber; ++j)
                  ids[j-1]=separatrices2_cells[ptr+j];

                outputSeparatrices2->InsertNextCell(VTK_POLYGON, vertexNumber, ids);
              }

              sourceIds->InsertNextTuple1(separatrices2_cells_sourceIds[i]);
              separatrixIds->InsertNextTuple1(separatrices2_cells_separatrixIds[i]);
              separatrixTypes->InsertNextTuple1(separatrices2_cells_separatrixTypes[i]);
              separatrixFunctionMaxima->InsertNextTuple1(
                  separatrices2_cells_separatrixFunctionMaxima[i]);

              separatrixFunctionMinima->InsertNextTuple1(
                  separatrices2_cells_separatrixFunctionMinima[i]);

              separatrixFunctionDiffs->InsertNextTuple1(
                  separatrices2_cells_separatrixFunctionDiffs[i]);
              isOnBoundary->InsertNextTuple1(separatrices2_cells_isOnBoundary[i]);

              ptr+=(separatrices2_cells[ptr]+1);
            }

            vtkCellData* cellData=outputSeparatrices2->GetCellData();
#ifndef withKamikaze
            if(!cellData){
              cerr << "[ttkMorseSmaleComplex] Error : "
                << "outputSeparatrices2 has no cell data." << endl;
              return -39;
            }
#endif

            cellData->AddArray(sourceIds);
            cellData->AddArray(separatrixIds);
            cellData->AddArray(separatrixTypes);
            cellData->AddArray(separatrixFunctionMaxima);
            cellData->AddArray(separatrixFunctionMinima);
            cellData->AddArray(separatrixFunctionDiffs);
            cellData->AddArray(isOnBoundary);
          }
    }));
  }

  // morse complexes
  outputMorseComplexes->ShallowCopy(input);
  if(ComputeAscendingSegmentation or ComputeDescendingSegmentation){
    vtkPointData* pointData=outputMorseComplexes->GetPointData();
#ifndef withKamikaze
    if(!pointData){
      cerr 
        << "[ttkMorseSmaleComplex] Error : outputMorseComplexes has no point "
        << "data." << endl;
      return -40;
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
