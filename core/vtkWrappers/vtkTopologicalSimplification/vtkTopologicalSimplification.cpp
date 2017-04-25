#include                  <vtkTopologicalSimplification.h>

vtkStandardNewMacro(vtkTopologicalSimplification)

  vtkTopologicalSimplification::vtkTopologicalSimplification():
    hasUpdatedMesh_{false},
    identifiers_{},
    inputScalars_{},
    offsets_{},
    inputOffsets_{}
{
  SetNumberOfInputPorts(2);
  triangulation_ = NULL;

  ScalarFieldId = 0;
  OutputOffsetScalarFieldName = "OutputOffsetScalarField";
  VertexIdentifierScalarField = "VertexIdentifier";
  ConsiderIdentifierAsBlackList = false;
  InputOffsetScalarFieldName = "OutputOffsetScalarField";
}

vtkTopologicalSimplification::~vtkTopologicalSimplification(){
  if(offsets_)
    offsets_->Delete();
}

int vtkTopologicalSimplification::FillInputPortInformation(int port, 
  vtkInformation *info){
  
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");

  return 1;
}

int vtkTopologicalSimplification::getTriangulation(vtkDataSet* input){
  
  triangulation_ = vtkTriangulation::getTriangulation(input);
  
  if(!triangulation_)
    return -1;
  
  triangulation_->setWrapper(this);
  topologicalSimplification_.setWrapper(this);
  topologicalSimplification_.setupTriangulation(triangulation_);
  Modified();
  hasUpdatedMesh_ = true;
  
#ifndef withKamikaze
  if(triangulation_->isEmpty()){
    cerr << "[vtkTopologicalSimplification] Error : vtkTriangulation allocation problem." << endl;
    return -1;
  }
#endif

  return 0;
}

int vtkTopologicalSimplification::getScalars(vtkDataSet* input){
  vtkPointData* pointData=input->GetPointData();

#ifndef withKamikaze
  if(!pointData){
    cerr << "[vtkTopologicalSimplification] Error : input has no point data." << endl;
    return -1;
  }
#endif

  if(ScalarField.length()){
    inputScalars_=pointData->GetArray(ScalarField.data());
  }
  else{
    inputScalars_=pointData->GetArray(ScalarFieldId);
    if(inputScalars_)
      ScalarField = inputScalars_->GetName();
  }

#ifndef withKamikaze
  if(!inputScalars_){
    cerr << "[vtkTopologicalSimplification] Error : input scalar field pointer is null." << endl;
    return -3;
  }
#endif

  return 0;
}

int vtkTopologicalSimplification::getIdentifiers(vtkPointSet* input){
  if(VertexIdentifierScalarField.length())
    identifiers_=input->GetPointData()->GetArray(VertexIdentifierScalarField.data());

#ifndef withKamikaze
  if(!identifiers_){
    cerr << "[vtkTopologicalSimplification] Error : wrong vertex identifier scalar field." << endl;
    return -1;
  }
#endif

  return 0;
}

int vtkTopologicalSimplification::getOffsets(vtkDataSet* input){
  if(UseInputOffsetScalarField and InputOffsetScalarFieldName.length())
    inputOffsets_=input->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  else{
    if(hasUpdatedMesh_ and offsets_){
      offsets_->Delete();
      offsets_=nullptr;
      hasUpdatedMesh_ = false;
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

#ifndef withKamikaze
  if(!inputOffsets_){
    cerr << "[vtkTopologicalSimplification] Error : wrong input offset scalar field." << endl;
    return -1;
  }
#endif

  return 0;
}

int vtkTopologicalSimplification::doIt(vector<vtkDataSet *> &inputs,
  vector<vtkDataSet *> &outputs){
  
  Memory m;
  
  vtkDataSet *domain = inputs[0];
  vtkPointSet *constraints = vtkPointSet::SafeDownCast(inputs[1]);
  vtkDataSet *output = outputs[0];
  
  int ret{};

  ret=getTriangulation(domain);
#ifndef withKamikaze
  if(ret){
    cerr << "[vtkTopologicalSimplification] Error : wrong triangulation." << endl;
    return -1;
  }
#endif

  ret=getScalars(domain);
#ifndef withKamikaze
  if(ret){
    cerr << "[vtkTopologicalSimplification] Error : wrong scalars." << endl;
    return -2;
  }
#endif

  ret=getIdentifiers(constraints);
#ifndef withKamikaze
  if(ret){
    cerr << "[vtkTopologicalSimplification] Error : wrong identifiers." << endl;
    return -3;
  }
#endif

  ret=getOffsets(domain);
#ifndef withKamikaze
  if(ret){
    cerr << "[vtkTopologicalSimplification] Error : wrong offsets." << endl;
    return -4;
  }
#endif

  const int numberOfVertices=domain->GetNumberOfPoints();
#ifndef withKamikaze
  if(numberOfVertices<=0){
    cerr << "[vtkTopologicalSimplification] Error : domain has no points." << endl;
    return -5;
  }
#endif

#ifndef withKamikaze
  if(OutputOffsetScalarFieldName.length()<=0){
    cerr << "[vtkTopologicalSimplification] Error : output offset scalar field has no name." << endl;
    return -6;
  }
#endif

  vtkSmartPointer<vtkIntArray> outputOffsets=vtkSmartPointer<vtkIntArray>::New();
  if(outputOffsets){
    outputOffsets->SetNumberOfComponents(1);
    outputOffsets->SetNumberOfTuples(numberOfVertices);
    outputOffsets->SetName(OutputOffsetScalarFieldName.data());
  }
#ifndef withKamikaze
  else{
    cerr << "[vtkTopologicalSimplification] Error : vtkIntArray allocation problem." << endl;
    return -7;
  }
#endif

  vtkDataArray* outputScalars{};
  switch(inputScalars_->GetDataType()){
    case VTK_DOUBLE:
      outputScalars=vtkDoubleArray::New();
      break;

    case VTK_FLOAT:
      outputScalars=vtkFloatArray::New();
      break;

    case VTK_INT:
      outputScalars=vtkIntArray::New();
      break;

    case VTK_SHORT:
      outputScalars=vtkShortArray::New();
      break;

    case VTK_UNSIGNED_SHORT:
      outputScalars=vtkUnsignedShortArray::New();
      break;

    case VTK_CHAR:
      outputScalars=vtkCharArray::New();
      break;

    case VTK_UNSIGNED_CHAR:
      outputScalars=vtkUnsignedCharArray::New();
      break;

#ifndef withKamikaze
    default:
      cerr << "[vtkTopologicalSimplification] Error : Unsupported data type." << endl;
      return -8;
#endif
  }
  if(outputScalars){
    outputScalars->SetNumberOfTuples(numberOfVertices);
    outputScalars->SetName(inputScalars_->GetName());
  }
#ifndef withKamikaze
  else{
    cerr << "[vtkTopologicalSimplification] Error : vtkDataArray allocation problem." << endl;
    return -9;
  }
#endif

  topologicalSimplification_.setVertexNumber(domain->GetNumberOfPoints());
  topologicalSimplification_.setConstraintNumber(
    constraints->GetNumberOfPoints());
  topologicalSimplification_.setInputScalarFieldPointer(
    inputScalars_->GetVoidPointer(0));
  topologicalSimplification_.setVertexIdentifierScalarFieldPointer(
    identifiers_->GetVoidPointer(0));
  topologicalSimplification_.setConsiderIdentifierAsBlackList(
    ConsiderIdentifierAsBlackList);
  topologicalSimplification_.setAddPerturbation(AddPerturbation);
  
  topologicalSimplification_.setInputOffsetScalarFieldPointer(
    inputOffsets_->GetVoidPointer(0));
  
  topologicalSimplification_.setOutputScalarFieldPointer(
    outputScalars->GetVoidPointer(0));
  
  topologicalSimplification_.setOutputOffsetScalarFieldPointer(
    outputOffsets->GetVoidPointer(0));

  switch(inputScalars_->GetDataType()){
    vtkTemplateMacro(({ret=topologicalSimplification_.execute<VTK_TT>();}));
  }
#ifndef withKamikaze
  // something wrong in baseCode
  if(ret){
    cerr << "[vtkTopologicalSimplification] TopologicalSimplification.execute() error code : " << ret << endl;
    return -10;
  }
#endif

  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(outputOffsets);
  output->GetPointData()->AddArray(outputScalars);
  outputScalars->Delete();

  {
    stringstream msg;
    msg << "[vtkTopologicalSimplification] Memory usage: " << 
m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
