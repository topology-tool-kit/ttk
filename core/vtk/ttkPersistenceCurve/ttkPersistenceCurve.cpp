#include                  <ttkPersistenceCurve.h>

using namespace ftm;

vtkStandardNewMacro(ttkPersistenceCurve)

  ttkPersistenceCurve::ttkPersistenceCurve():
    UseAllCores{},
    inputScalars_{},
    JTPersistenceCurve_{vtkTable::New()},
    MSCPersistenceCurve_{vtkTable::New()},
    STPersistenceCurve_{vtkTable::New()},
    CTPersistenceCurve_{vtkTable::New()},
    offsets_{},
    inputOffsets_{},
    varyingMesh_{}
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(4);
  
  triangulation_ = NULL;
  ScalarFieldId = 0;
  InputOffsetScalarFieldName = "OutputOffsetScalarField";

  inputTriangulation_ = vtkSmartPointer<ttkTriangulationFilter>::New();
}

ttkPersistenceCurve::~ttkPersistenceCurve(){
  if(JTPersistenceCurve_)
    JTPersistenceCurve_->Delete();
  if(MSCPersistenceCurve_)
    MSCPersistenceCurve_->Delete();
  if(STPersistenceCurve_)
    STPersistenceCurve_->Delete();
  if(CTPersistenceCurve_)
    CTPersistenceCurve_->Delete();
  if(offsets_)
    offsets_->Delete();
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkPersistenceCurve::needsToAbort(){
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkPersistenceCurve::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkPersistenceCurve] " << progress*100
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkPersistenceCurve::FillOutputPortInformation(int port, vtkInformation* info){
  switch (port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
      break;

    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
      break;

    case 2:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
      break;

    case 3:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
      break;
  }

  return 1;
}

int ttkPersistenceCurve::getScalars(vtkDataSet* input){
  vtkPointData* pointData=input->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData){
    cerr << "[ttkPersistenceCurve] Error : input has no point data." << endl;
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

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars_){
    cerr << "[ttkPersistenceCurve] Error : input scalar field pointer is null." 
      << endl;
    return -2;
  }
#endif

  return 0;
}

int ttkPersistenceCurve::getTriangulation(vtkDataSet* input){
  
  varyingMesh_=false;

  triangulation_ = ttkTriangulation::getTriangulation(input);
  
  if(!triangulation_)
    return 0;
  
  triangulation_->setWrapper(this);
  persistenceCurve_.setupTriangulation(triangulation_);
  
  if(triangulation_->isEmpty() 
    or ttkTriangulation::hasChangedConnectivity(triangulation_, input, this)){
    Modified();
    varyingMesh_=true;
  }

  return 0;
}

int ttkPersistenceCurve::getOffsets(vtkDataSet* input){
  if(UseInputOffsetScalarField and InputOffsetScalarFieldName.length())
    inputOffsets_ =
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

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOffsets_){
    cerr << "[ttkPersistenceCurve] Error : wrong input offset scalar field." 
      << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkPersistenceCurve::doIt(vtkDataSet *input,
    vtkTable* outputJTPersistenceCurve,
    vtkTable* outputMSCPersistenceCurve,
    vtkTable* outputSTPersistenceCurve,
    vtkTable* outputCTPersistenceCurve){
  int ret{};

  ret=getScalars(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret){
    cerr << "[ttkPersistenceCurve] Error : wrong scalars." << endl;
    return -1;
  }
#endif

  ret=getTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret){
    cerr << "[ttkPersistenceCurve] Error : wrong triangulation." << endl;
    return -2;
  }
#endif

  ret=getOffsets(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret){
    cerr << "[ttkPersistenceCurve] Error : wrong offsets." << endl;
    return -3;
  }
#endif

  persistenceCurve_.setWrapper(this);
  persistenceCurve_.setInputScalars(inputScalars_->GetVoidPointer(0));
  persistenceCurve_.setInputOffsets(inputOffsets_->GetVoidPointer(0));
  persistenceCurve_.setComputeSaddleConnectors(ComputeSaddleConnectors);
  switch(inputScalars_->GetDataType()){
    vtkTemplateMacro(({
          vector<pair<VTK_TT, idVertex>> JTPlot;
          vector<pair<VTK_TT, idVertex>> STPlot;
          vector<pair<VTK_TT, idVertex>> MSCPlot;
          vector<pair<VTK_TT, idVertex>> CTPlot;

          persistenceCurve_.setOutputJTPlot(&JTPlot);
          persistenceCurve_.setOutputMSCPlot(&MSCPlot);
          persistenceCurve_.setOutputSTPlot(&STPlot);
          persistenceCurve_.setOutputCTPlot(&CTPlot);
          ret=persistenceCurve_.execute<VTK_TT>();

#ifndef TTK_ENABLE_KAMIKAZE
          if(ret){
          cerr 
            << "[ttkPersistenceCurve] PersistenceCurve.execute() error code : " 
            << ret << endl;
          return -4;
          }
#endif

          ret=getPersistenceCurve<vtkDoubleArray,VTK_TT>(TreeType::Join, JTPlot);
#ifndef TTK_ENABLE_KAMIKAZE
          if(ret){
          cerr << "[ttkPersistenceCurve] Error :"
            << " build of join tree persistence curve has failed." << endl;
          return -5;
          }
#endif

          ret=getMSCPersistenceCurve<vtkDoubleArray,VTK_TT>(MSCPlot);
#ifndef TTK_ENABLE_KAMIKAZE
          if(ret){
          cerr << "[ttkPersistenceCurve] Error : "
            << "build of saddle-saddle persistence curve has failed." << endl;
          return -6;
          }
#endif

          ret=getPersistenceCurve<vtkDoubleArray,VTK_TT>(TreeType::Split, STPlot);
#ifndef TTK_ENABLE_KAMIKAZE
          if(ret){
          cerr << "[ttkPersistenceCurve] Error : "
            << "build of split tree persistence curve has failed." << endl;
          return -7;
          }
#endif

          ret=getPersistenceCurve<vtkDoubleArray,VTK_TT>(TreeType::Contour, CTPlot);
#ifndef TTK_ENABLE_KAMIKAZE
          if(ret){
          cerr << "[ttkPersistenceCurve] Error : "
            << "build of contour tree persistence curve has failed." << endl;
          return -8;
          }
#endif
    }));
  }

  outputJTPersistenceCurve->ShallowCopy(JTPersistenceCurve_);
  outputMSCPersistenceCurve->ShallowCopy(MSCPersistenceCurve_);
  outputSTPersistenceCurve->ShallowCopy(STPersistenceCurve_);
  outputCTPersistenceCurve->ShallowCopy(CTPersistenceCurve_);

  return 0;
}

int ttkPersistenceCurve::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  
  inputTriangulation_->SetInputData(input);
  inputTriangulation_->Update();

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkTable* outputJTPersistenceCurve = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  outInfo = outputVector->GetInformationObject(1);
  vtkTable* outputMSCPersistenceCurve = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  outInfo = outputVector->GetInformationObject(2);
  vtkTable* outputSTPersistenceCurve = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  outInfo = outputVector->GetInformationObject(3);
  vtkTable* outputCTPersistenceCurve = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  doIt(inputTriangulation_->GetOutput(),
      outputJTPersistenceCurve,
      outputMSCPersistenceCurve,
      outputSTPersistenceCurve,
      outputCTPersistenceCurve);

  {
    stringstream msg;
    msg << "[ttkPersistenceCurve] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
