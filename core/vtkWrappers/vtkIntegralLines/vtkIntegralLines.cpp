#include<vtkIntegralLines.h>

vtkStandardNewMacro(vtkIntegralLines)

  vtkIntegralLines::vtkIntegralLines():
    hasUpdatedMesh_{false},
    inputScalars_{nullptr},
    offsets_{nullptr},
    inputOffsets_{nullptr},
    identifiers_{nullptr}
{
  SetNumberOfInputPorts(2);
  triangulation_ = NULL;
  
  OffsetScalarFieldName = "OutputOffsetScalarField";
}

vtkIntegralLines::~vtkIntegralLines(){
  if(offsets_)
    offsets_->Delete();
}

int vtkIntegralLines::FillInputPortInformation(int port, vtkInformation* info){
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");

  return 1;
}

int vtkIntegralLines::FillOutputPortInformation(int port, vtkInformation* info){
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");

  return 1;
}

int vtkIntegralLines::getTriangulation(vtkDataSet* input){
 
  triangulation_ = vtkTriangulation::getTriangulation(input);
  
  if(!triangulation_)
    return -1;
  
  triangulation_->setWrapper(this);
  integralLines_.setWrapper(this);
  
  integralLines_.setupTriangulation(triangulation_);
  Modified();
  hasUpdatedMesh_ = true;
  
#ifndef withKamikaze
  // allocation problem
  if(triangulation_->isEmpty()){
    cerr << "[vtkIntegralLines] Error : vtkTriangulation allocation problem." << endl;
    return -1;
  }
#endif

  return 0;
}

int vtkIntegralLines::getScalars(vtkDataSet* input){
  vtkPointData* pointData=input->GetPointData();

#ifndef withKamikaze
  if(!pointData){
    cerr << "[vtkIntegralLines] Error : input has no point data." << endl;
    return -1;
  }

  if(!ScalarField.length()){
    cerr << "[vtkIntegralLines] Error : scalar field has no name." << endl;
    return -2;
  }
#endif

  inputScalars_=pointData->GetArray(ScalarField.data());

#ifndef withKamikaze
  if(!inputScalars_){
    cerr << "[vtkIntegralLines] Error : input scalar field pointer is null." << endl;
    return -3;
  }
#endif

  return 0;
}

int vtkIntegralLines::getOffsets(vtkDataSet* input){
  if(UseOffsetScalarField and OffsetScalarFieldName.length()){
    inputOffsets_=input->GetPointData()->GetArray(OffsetScalarFieldName.data());
  }
  else{
    if(hasUpdatedMesh_ and offsets_){
      offsets_->Delete();
      offsets_=nullptr;
      hasUpdatedMesh_ = false;
    }

    if(!offsets_){
      const int numberOfPoints=input->GetNumberOfPoints();

      offsets_=vtkIntArray::New();
      offsets_->SetNumberOfComponents(1);
      offsets_->SetNumberOfTuples(numberOfPoints);
      offsets_->SetName("OffsetScalarField");

      for(int i=0; i<numberOfPoints; ++i)
        offsets_->SetTuple1(i,i);
    }

    inputOffsets_=offsets_;
  }

#ifndef withKamikaze
  // allocation problem
  if(!inputOffsets_){
    cerr << "[vtkIntegralLines] Error : wrong input offset scalar field." << endl;
    return -1;
  }
#endif

  return 0;
}

int vtkIntegralLines::getIdentifiers(vtkPointSet* input){
  if(VertexIdentifierScalarFieldName.length())
    identifiers_=input->GetPointData()->GetArray(VertexIdentifierScalarFieldName.data());

#ifndef withKamikaze
  // allocation problem
  if(!identifiers_){
    cerr << "[vtkIntegralLines] Error : wrong input vertex identifier scalar field." << endl;
    return -1;
  }
#endif

  return 0;
}

int vtkIntegralLines::getTrajectories(vtkDataSet* input,
    vector<vector<int>>& trajectories,
    vtkUnstructuredGrid* output){
  vtkSmartPointer<vtkUnstructuredGrid> ug=vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> pts=vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkFloatArray> dist=vtkSmartPointer<vtkFloatArray>::New();
  dist->SetNumberOfComponents(1);
  dist->SetName("DistanceFromSeed");

  // here, copy the original scalars
  int numberOfArrays=input->GetPointData()->GetNumberOfArrays();

  vector<vtkDataArray*> scalarArrays;
  for(int k=0; k<numberOfArrays; ++k){
    auto a=input->GetPointData()->GetArray(k);

    if(a->GetNumberOfComponents()==1)
      scalarArrays.push_back(a);
  }
  // not efficient, implicit conversion to double
  vector<vtkSmartPointer<vtkDoubleArray>> inputScalars(scalarArrays.size());
  for(unsigned int k=0; k<scalarArrays.size(); ++k){
    inputScalars[k]=vtkSmartPointer<vtkDoubleArray>::New();
    inputScalars[k]->SetNumberOfComponents(1);
    inputScalars[k]->SetName(scalarArrays[k]->GetName());
  }

  float p0[3];
  float p1[3];
  vtkIdType ids[2];
  for(unsigned int i=0; i<trajectories.size(); ++i){
    if(trajectories[i].size()){
      int vertex=trajectories[i][0];
      // init
      triangulation_->getVertexPoint(vertex,p0[0],p0[1],p0[2]);
      ids[0]=pts->InsertNextPoint(p0);
      // distanceScalars
      float distanceFromSeed{};
      dist->InsertNextTuple1(distanceFromSeed);
      // inputScalars
      for(unsigned int k=0; k<scalarArrays.size(); ++k)
        inputScalars[k]->InsertNextTuple1(scalarArrays[k]->GetTuple1(vertex));

      for(unsigned int j=1; j<trajectories[i].size(); ++j){
        vertex=trajectories[i][j];
        triangulation_->getVertexPoint(vertex,p1[0],p1[1],p1[2]);
        ids[1]=pts->InsertNextPoint(p1);
        // distanceScalars
        distanceFromSeed+=Geometry::distance(p0,p1,3);
        dist->InsertNextTuple1(distanceFromSeed);
        // inputScalars
        for(unsigned int k=0; k<scalarArrays.size(); ++k)
          inputScalars[k]->InsertNextTuple1(scalarArrays[k]->GetTuple1(vertex));

        ug->InsertNextCell(VTK_LINE,2,ids);

        // iteration
        ids[0]=ids[1];
        p0[0]=p1[0];
        p0[1]=p1[1];
        p0[2]=p1[2];
      }
    }
  }
  ug->SetPoints(pts);
  ug->GetPointData()->AddArray(dist);
  for(unsigned int k=0; k<scalarArrays.size(); ++k)
    ug->GetPointData()->AddArray(inputScalars[k]);

  output->ShallowCopy(ug);

  return 0;
}

int vtkIntegralLines::doIt(vector<vtkDataSet *> &inputs,
  vector<vtkDataSet *> &outputs){
//                            vtkPointSet* seeds, vtkUnstructuredGrid* output){
  
  Memory m;
  
  vtkDataSet *domain = inputs[0];
  vtkPointSet *seeds = vtkPointSet::SafeDownCast(inputs[1]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  
  int ret{};

  ret=getTriangulation(domain);
#ifndef withKamikaze
  // triangulation problem
  if(ret){
    cerr << "[vtkIntegralLines] Error : wrong triangulation." << endl;
    return -1;
  }
#endif

  ret=getScalars(domain);
#ifndef withKamikaze
  // field problem
  if(ret){
    cerr << "[vtkIntegralLines] Error : wrong scalar field." << endl;
    return -2;
  }
#endif

  ret=getOffsets(domain);
#ifndef withKamikaze
  // field problem
  if(ret){
    cerr << "[vtkIntegralLines] Error : wrong offsets." << endl;
    return -3;
  }
#endif

  ret=getIdentifiers(seeds);
#ifndef withKamikaze
  // field problem
  if(ret){
    cerr << "[vtkIntegralLines] Error : wrong identifiers." << endl;
    return -3;
  }
#endif

  const int numberOfPointsInDomain=domain->GetNumberOfPoints();
#ifndef withKamikaze
  // no points.
  if(numberOfPointsInDomain<=0){
    cerr << "[vtkIntegralLines] Error : domain has no points." << endl;
    return -4;
  }
#endif

  const int numberOfPointsInSeeds=seeds->GetNumberOfPoints();
#ifndef withKamikaze
  // no points.
  if(numberOfPointsInSeeds<=0){
    cerr << "[vtkIntegralLines] Error : seeds have no points." << endl;
    return -5;
  }
#endif

  vector<vector<int>> trajectories;

  integralLines_.setVertexNumber(numberOfPointsInDomain);
  integralLines_.setSeedNumber(numberOfPointsInSeeds);
  integralLines_.setDirection(Direction);
  integralLines_.setInputScalarField(inputScalars_->GetVoidPointer(0));
  integralLines_.setInputOffsets(inputOffsets_->GetVoidPointer(0));
  
  integralLines_.setVertexIdentifierScalarField(
    identifiers_->GetVoidPointer(0));
  integralLines_.setOutputTrajectories(&trajectories);
  
  switch(inputScalars_->GetDataType()){
    vtkTemplateMacro(({ret=integralLines_.execute<VTK_TT>();}));
  }
#ifndef withKamikaze
  // something wrong in baseCode
  if(ret){
    cerr << "[vtkIntegralLines] IntegralLines.execute() error code : " << ret << endl;
    return -6;
  }
#endif

  // make the vtk trajectories
  ret=getTrajectories(domain, trajectories, output);
#ifndef withKamikaze
  // trajectories problem
  if(ret){
    cerr << "[vtkIntegralLines] Error : wrong trajectories." << endl;
    return -7;
  }
#endif

  {
    stringstream msg;
    msg << "[vtkIntegralLines] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}