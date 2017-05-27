#include                  <ttkPersistenceDiagram.h>

vtkStandardNewMacro(ttkPersistenceDiagram)

  ttkPersistenceDiagram::ttkPersistenceDiagram():
    UseAllCores{},
    inputScalars_{},
    CTPersistenceDiagram_{vtkUnstructuredGrid::New()},
    offsets_{},
    inputOffsets_{},
    varyingMesh_{}
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);

  ScalarFieldId = 0;
  OffsetFieldId = -1;
  ComputeSaddleConnectors = false;
  InputOffsetScalarFieldName = "OutputOffsetScalarField";
  
  triangulation_ = NULL;
}

ttkPersistenceDiagram::~ttkPersistenceDiagram(){
  if(CTPersistenceDiagram_)
    CTPersistenceDiagram_->Delete();
  if(offsets_)
    offsets_->Delete();
}

int ttkPersistenceDiagram::FillOutputPortInformation(int port, 
  vtkInformation *info){
  switch (port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
  }

  return 1;
}

int ttkPersistenceDiagram::getScalars(vtkDataSet* input){
  vtkPointData* pointData=input->GetPointData();

#ifndef withKamikaze
  if(!pointData){
    cerr << "[ttkPersistenceDiagram] Error : input has no point data." << endl;
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
    cerr 
      << "[ttkPersistenceDiagram] Error : input scalar field pointer is null." 
      << endl;
    return -3;
  }
#endif

  stringstream msg;
  msg << "[ttkPersistenceDiagram] Starting computation on field `"
    << inputScalars_->GetName() << "'..." << endl;
  dMsg(cout, msg.str(), infoMsg);

  return 0;
}

int ttkPersistenceDiagram::getTriangulation(vtkDataSet* input){
  varyingMesh_=false;

  triangulation_ = ttkTriangulation::getTriangulation(input);
  if(!triangulation_)
    return 0;
  
  triangulation_->setWrapper(this);
  persistenceDiagram_.setupTriangulation(triangulation_);
  
  if(triangulation_->isEmpty() 
    or ttkTriangulation::hasChangedConnectivity(triangulation_, input, this)){
    Modified();
    varyingMesh_=true;
  }

  return 0;
}

int ttkPersistenceDiagram::getOffsets(vtkDataSet* input){

  if(OffsetFieldId != -1){
    inputOffsets_ = input->GetPointData()->GetArray(OffsetFieldId);
    if(inputOffsets_){
      InputOffsetScalarFieldName = inputOffsets_->GetName();
      UseInputOffsetScalarField = true;
    }
  }

  if(UseInputOffsetScalarField and InputOffsetScalarFieldName.length())
    inputOffsets_=
      input->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  else{
    if(varyingMesh_ and offsets_){
      offsets_->Delete();
      offsets_=nullptr;
    }

    if(!offsets_){
      const int numberOfVertices=input->GetNumberOfPoints();

      offsets_=vtkIntArray::New();
      offsets_->SetNumberOfComponents(1);
      offsets_->SetNumberOfTuples(numberOfVertices);
      offsets_->SetName("OffsetScalarField");
      for(int i=0; i<numberOfVertices; ++i)
        offsets_->SetTuple1(i,i);
    }

    inputOffsets_=offsets_;
  }

#ifndef withKamikaze
  if(!inputOffsets_){
    cerr << "[ttkPersistenceDiagram] Error : wrong input offset scalar field." 
      << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkPersistenceDiagram::doIt(vector<vtkDataSet *> &inputs,
  vector<vtkDataSet *> &outputs){
  
  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *outputCTPersistenceDiagram = 
    vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  
  Memory m;
  
  int ret{};

  ret=getScalars(input);
#ifndef withKamikaze
  if(ret){
    cerr << "[ttkPersistenceDiagram] Error : wrong scalars." << endl;
    return -1;
  }
#endif

  ret=getTriangulation(input);
#ifndef withKamikaze
  if(ret){
    cerr << "[ttkPersistenceDiagram] Error : wrong triangulation." << endl;
    return -2;
  }
#endif

  ret=getOffsets(input);
#ifndef withKamikaze
  if(ret){
    cerr << "[ttkPersistenceDiagram] Error : wrong offsets." << endl;
    return -3;
  }
#endif

  vector<tuple<Cell,Cell>> dmt_pairs;
  persistenceDiagram_.setDMTPairs(&dmt_pairs);

  persistenceDiagram_.setWrapper(this);
  persistenceDiagram_.setInputScalars(inputScalars_->GetVoidPointer(0));
  persistenceDiagram_.setInputOffsets(inputOffsets_->GetVoidPointer(0));
  persistenceDiagram_.setComputeSaddleConnectors(ComputeSaddleConnectors);
  switch(inputScalars_->GetDataType()){
    vtkTemplateMacro(({
          vector<tuple<idVertex,NodeType,idVertex,NodeType,VTK_TT,idVertex>> 
            CTDiagram;

          persistenceDiagram_.setOutputCTDiagram(&CTDiagram);
          ret=persistenceDiagram_.execute<VTK_TT>();
#ifndef withKamikaze
          if(ret){
          cerr << "[ttkPersistenceDiagram] PersistenceDiagram.execute() "
            << "error code : " << ret << endl;
          return -4;
          }
#endif

          ret=getPersistenceDiagram<VTK_TT>(TreeType::Contour, CTDiagram);
#ifndef withKamikaze
          if(ret){
            cerr << "[ttkPersistenceDiagram] Error : "
              << "build of contour tree persistence diagram has failed." 
              << endl;
            return -5;
          }
#endif
    }));
  }

  outputCTPersistenceDiagram->ShallowCopy(CTPersistenceDiagram_);

  {
    stringstream msg;
    msg << "[ttkPersistenceDiagram] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}