#include<ttkDiscreteGradient.h>

vtkStandardNewMacro(ttkDiscreteGradient)

  ttkDiscreteGradient::ttkDiscreteGradient():
    UseAllCores{},
    ThreadNumber{},
    ScalarField{},
    InputOffsetScalarFieldName{},
    UseInputOffsetScalarField{},
    ReverseSaddleMaximumConnection{},
    ReverseSaddleSaddleConnection{},
    AllowSecondPass{},
    AllowThirdPass{},
    ComputeGradientGlyphs{},
    ScalarFieldId{},
    OffsetFieldId{},

    triangulation_{},
    inputScalars_{},
    offsets_{},
    inputOffsets_{},
    hasUpdatedMesh_{}
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(2);
}

ttkDiscreteGradient::~ttkDiscreteGradient(){
  if(offsets_)
    offsets_->Delete();
}

int ttkDiscreteGradient::FillInputPortInformation(int port, vtkInformation* info)
{
  switch(port){
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
      break;
  }

  return 1;
}

int ttkDiscreteGradient::FillOutputPortInformation(int port, vtkInformation* info){
  switch(port){
    case 0:
    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
  }

  return 1;
}

int ttkDiscreteGradient::setupTriangulation(vtkDataSet* input){
  triangulation_=ttkTriangulation::getTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_){
    cerr << "[ttkDiscreteGradient] Error : ttkTriangulation::getTriangulation() is null." << endl;
    return -1;
  }
#endif

  hasUpdatedMesh_=ttkTriangulation::hasChangedConnectivity(triangulation_, input, this);

  triangulation_->setWrapper(this);
  discreteGradient_.setWrapper(this);
  discreteGradient_.setupTriangulation(triangulation_);
  Modified();

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_->isEmpty()){
    cerr << "[vtkIntegralLines] Error : vtkTriangulation allocation problem." << endl;
    return -1;
  }
#endif
  return 0;
}

int ttkDiscreteGradient::getScalars(vtkDataSet* input){
  vtkPointData* pointData=input->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData){
    cerr << "[ttkDiscreteGradient] Error : input has no point data." << endl;
    return -1;
  }

  if(!ScalarField.length()){
    cerr << "[ttkDiscreteGradient] Error : scalar field has no name." << endl;
    return -2;
  }
#endif

  if(ScalarField.length()){
    inputScalars_=pointData->GetArray(ScalarField.data());
  }
  else{
    inputScalars_=pointData->GetArray(ScalarFieldId);
    if(inputScalars_)
      ScalarField=inputScalars_->GetName();
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars_){
    cerr << "[ttkDiscreteGradient] Error : input scalar field pointer is null." << endl;
    return -3;
  }
#endif

  return 0;
}

int ttkDiscreteGradient::getOffsets(vtkDataSet* input){
  if(OffsetFieldId != -1){
    inputOffsets_=input->GetPointData()->GetArray(OffsetFieldId);
    if(inputOffsets_){
      InputOffsetScalarFieldName=inputOffsets_->GetName();
      UseInputOffsetScalarField=true;
    }
  }

  if(UseInputOffsetScalarField and InputOffsetScalarFieldName.length())
    inputOffsets_=input->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  else{
    if(hasUpdatedMesh_ and offsets_){
      offsets_->Delete();
      offsets_=nullptr;
    }

    if(!offsets_){
      const int numberOfVertices=input->GetNumberOfPoints();

      offsets_=vtkIntArray::New();
      offsets_->SetNumberOfComponents(1);
      offsets_->SetNumberOfTuples(numberOfVertices);
      offsets_->SetName("OffsetsScalarField");
      for(int i=0; i<numberOfVertices; ++i)
        offsets_->SetTuple1(i,i);
    }

    inputOffsets_=offsets_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOffsets_){
    cerr << "[ttkDiscreteGradient] Error : wrong input offset scalar field." << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkDiscreteGradient::doIt(vector<vtkDataSet *> &inputs,
    vector<vtkDataSet *> &outputs){
  vtkDataSet* input=inputs[0];
  vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  vtkUnstructuredGrid* outputGradientGlyphs=vtkUnstructuredGrid::SafeDownCast(outputs[1]);

  int ret{};

  ret=setupTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret){
    cerr << "[ttkDiscreteGradient] Error : wrong triangulation." << endl;
    return -1;
  }
#endif

  ret=getScalars(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret){
    cerr << "[ttkDiscreteGradient] Error : wrong scalars." << endl;
    return -2;
  }
#endif

  ret=getOffsets(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret){
    cerr << "[ttkDiscreteGradient] Error : wrong offsets." << endl;
    return -3;
  }
#endif

  // critical points
  int criticalPoints_numberOfPoints{};
  vector<float> criticalPoints_points;
  vector<int> criticalPoints_points_cellDimensions;
  vector<int> criticalPoints_points_cellIds;
  vector<char> criticalPoints_points_isOnBoundary;
  vector<int> criticalPoints_points_PLVertexIdentifiers;
  vector<int> criticalPoints_points_manifoldSize;

  // gradient pairs
  int gradientGlyphs_numberOfPoints{};
  vector<float> gradientGlyphs_points;
  vector<int> gradientGlyphs_points_pairOrigins;
  int gradientGlyphs_numberOfCells{};
  vector<int> gradientGlyphs_cells;
  vector<int> gradientGlyphs_cells_pairTypes;

  // baseCode processing
  discreteGradient_.setWrapper(this);
  discreteGradient_.setIterationThreshold(IterationThreshold);
  discreteGradient_.setReverseSaddleMaximumConnection(ReverseSaddleMaximumConnection);
  discreteGradient_.setReverseSaddleSaddleConnection(ReverseSaddleSaddleConnection);
  discreteGradient_.setInputScalarField(inputScalars_->GetVoidPointer(0));
  discreteGradient_.setInputOffsets(inputOffsets_->GetVoidPointer(0));

  discreteGradient_.setOutputGradientGlyphs(&gradientGlyphs_numberOfPoints,
      &gradientGlyphs_points,
      &gradientGlyphs_points_pairOrigins,
      &gradientGlyphs_numberOfCells,
      &gradientGlyphs_cells,
      &gradientGlyphs_cells_pairTypes);

  const int dimensionality=triangulation_->getDimensionality();

  switch(inputScalars_->GetDataType()){
    vtkTemplateMacro(({
          vector<VTK_TT> criticalPoints_points_cellScalars;

          discreteGradient_.setOutputCriticalPoints(&criticalPoints_numberOfPoints,
              &criticalPoints_points,
              &criticalPoints_points_cellDimensions,
              &criticalPoints_points_cellIds,
              &criticalPoints_points_cellScalars,
              &criticalPoints_points_isOnBoundary,
              &criticalPoints_points_PLVertexIdentifiers,
              &criticalPoints_points_manifoldSize);

          ret=discreteGradient_.buildGradient<VTK_TT>();
#ifndef TTK_ENABLE_KAMIKAZE
          if(ret){
          cerr << "[ttkDiscreteGradient] Error : DiscreteGradient.buildGradient() error code : " << ret << endl;
          return -8;
          }
#endif

          if(AllowSecondPass){
            ret=discreteGradient_.buildGradient2<VTK_TT>();
#ifndef TTK_ENABLE_KAMIKAZE
            if(ret){
              cerr << "[ttkDiscreteGradient] Error : DiscreteGradient.buildGradient2() error code : " << ret << endl;
              return -9;
            }
#endif
          }

          if(dimensionality==3 and AllowThirdPass){
            ret=discreteGradient_.buildGradient3<VTK_TT>();
#ifndef TTK_ENABLE_KAMIKAZE
            if(ret){
              cerr << "[ttkDiscreteGradient] Error : DiscreteGradient.buildGradient2() error code : " << ret << endl;
              return -10;
            }
#endif
          }

          ret=discreteGradient_.reverseGradient<VTK_TT>();
#ifndef TTK_ENABLE_KAMIKAZE
          if(ret){
            cerr << "[ttkDiscreteGradient] Error : DiscreteGradient.reverseGradient() error code : " << ret << endl;
            return -11;
          }
#endif

          // critical points
          {
            discreteGradient_.setCriticalPoints<VTK_TT>();

            vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
#ifndef TTK_ENABLE_KAMIKAZE
            if(!points){
              cerr << "[ttkDiscreteGradient] Error : vtkPoints allocation problem." << endl;
              return -12;
            }
#endif

            vtkSmartPointer<vtkIntArray> cellDimensions=vtkSmartPointer<vtkIntArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
            if(!cellDimensions){
              cerr << "[ttkDiscreteGradient] Error : vtkIntArray allocation problem." << endl;
              return -13;
            }
#endif
            cellDimensions->SetNumberOfComponents(1);
            cellDimensions->SetName("CellDimension");

            vtkSmartPointer<vtkIntArray> cellIds=vtkSmartPointer<vtkIntArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
            if(!cellIds){
              cerr << "[ttkDiscreteGradient] Error : vtkIntArray allocation problem." << endl;
              return -14;
            }
#endif
            cellIds->SetNumberOfComponents(1);
            cellIds->SetName("CellId");

            vtkDataArray* cellScalars=inputScalars_->NewInstance();
#ifndef TTK_ENABLE_KAMIKAZE
            if(!cellScalars){
              cerr << "[ttkDiscreteGradient] Error : vtkDataArray allocation problem." << endl;
              return -15;
            }
#endif
            cellScalars->SetNumberOfComponents(1);
            cellScalars->SetName(ScalarField.data());

            vtkSmartPointer<vtkCharArray> isOnBoundary=vtkSmartPointer<vtkCharArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
            if(!isOnBoundary){
              cerr << "[vtkMorseSmaleComplex] Error : vtkCharArray allocation problem." << endl;
              return -16;
            }
#endif
            isOnBoundary->SetNumberOfComponents(1);
            isOnBoundary->SetName("IsOnBoundary");

            vtkSmartPointer<vtkIntArray> PLVertexIdentifiers=
              vtkSmartPointer<vtkIntArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
            if(!PLVertexIdentifiers){
              cerr << "[ttkMorseSmaleComplex] Error : vtkIntArray allocation "
                << "problem." << endl;
              return -10;
            }
#endif
            PLVertexIdentifiers->SetNumberOfComponents(1);
            PLVertexIdentifiers->SetName("VertexIdentifier");

            for(int i=0; i<criticalPoints_numberOfPoints; ++i){
              points->InsertNextPoint(criticalPoints_points[3*i],
                  criticalPoints_points[3*i+1],
                  criticalPoints_points[3*i+2]);

              cellDimensions->InsertNextTuple1(criticalPoints_points_cellDimensions[i]);
              cellIds->InsertNextTuple1(criticalPoints_points_cellIds[i]);
              cellScalars->InsertNextTuple1(criticalPoints_points_cellScalars[i]);
              isOnBoundary->InsertNextTuple1(criticalPoints_points_isOnBoundary[i]);
              PLVertexIdentifiers->InsertNextTuple1(criticalPoints_points_PLVertexIdentifiers[i]);
            }
            outputCriticalPoints->SetPoints(points);

            vtkPointData* pointData=outputCriticalPoints->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
            if(!pointData){
              cerr << "[ttkDiscreteGradient] Error : outputCriticalPoints has no point data." << endl;
              return -17;
            }
#endif

            pointData->AddArray(cellDimensions);
            pointData->AddArray(cellIds);
            pointData->AddArray(cellScalars);
            pointData->AddArray(isOnBoundary);
            pointData->AddArray(PLVertexIdentifiers);
          }

    }));
  }

  // gradient glyphs
  if(ComputeGradientGlyphs){
    discreteGradient_.setGradientGlyphs();

    vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!points){
      cerr << "[ttkDiscreteGradient] Error : vtkPoints allocation problem." << endl;
      return -18;
    }
#endif

    vtkSmartPointer<vtkIntArray> pairOrigins=vtkSmartPointer<vtkIntArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pairOrigins){
      cerr << "[ttkDiscreteGradient] Error : vtkIntArray allocation problem." << endl;
      return -19;
    }
#endif
    pairOrigins->SetNumberOfComponents(1);
    pairOrigins->SetName("PairOrigin");

    vtkSmartPointer<vtkIntArray> pairTypes=vtkSmartPointer<vtkIntArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pairTypes){
      cerr << "[ttkDiscreteGradient] Error : vtkIntArray allocation problem." << endl;
      return -20;
    }
#endif
    pairTypes->SetNumberOfComponents(1);
    pairTypes->SetName("PairType");

    for(int i=0; i<gradientGlyphs_numberOfPoints; ++i){
      points->InsertNextPoint(gradientGlyphs_points[3*i],
          gradientGlyphs_points[3*i+1],
          gradientGlyphs_points[3*i+2]);

      pairOrigins->InsertNextTuple1(gradientGlyphs_points_pairOrigins[i]);
    }
    outputGradientGlyphs->SetPoints(points);

    outputGradientGlyphs->Allocate(gradientGlyphs_numberOfCells);
    int ptr{};
    for(int i=0; i<gradientGlyphs_numberOfCells; ++i){
      vtkIdType line[2];
      line[0]=gradientGlyphs_cells[ptr+1];
      line[1]=gradientGlyphs_cells[ptr+2];

      outputGradientGlyphs->InsertNextCell(VTK_LINE, 2, line);

      pairTypes->InsertNextTuple1(gradientGlyphs_cells_pairTypes[i]);

      ptr+=(gradientGlyphs_cells[ptr]+1);
    }

    vtkPointData* pointData=outputGradientGlyphs->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData){
      cerr << "[ttkDiscreteGradient] Error : outputGradientGlyphs has no point data." << endl;
      return -21;
    }
#endif

    pointData->AddArray(pairOrigins);

    vtkCellData* cellData=outputGradientGlyphs->GetCellData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!cellData){
      cerr << "[ttkDiscreteGradient] Error : outputGradientGlyphs has no cell data." << endl;
      return -22;
    }
#endif

    cellData->AddArray(pairTypes);
  }

  return 0;
}

