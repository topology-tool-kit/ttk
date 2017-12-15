#include <ttkFieldSelector.h>

vtkStandardNewMacro(ttkFieldSelector)

  // transmit abort signals
  bool ttkFieldSelector::needsToAbort(){
    return GetAbortExecute();
  }

// transmit progress status
int ttkFieldSelector::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkFieldSelector] " << progress*100 
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkFieldSelector::doIt(vtkDataSet* input, vtkDataSet* output){
  Memory m;

  output->ShallowCopy(input);

  vtkPointData* inputPointData=input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputPointData){
    cerr << "[ttkFieldSelector] Error: input has no point data." << endl;
    return -1;
  }
#endif

  vtkSmartPointer<vtkPointData> outputPointData=vtkSmartPointer<vtkPointData>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputPointData){
    cerr << "[ttkFieldSelector] Error: vtkPointData memory allocation problem." << endl;
    return -1;
  }
#endif

  for(auto& scalar : ScalarFields){
    vtkDataArray* arr=inputPointData->GetArray(scalar.data());
    if(arr) outputPointData->AddArray(arr);
  }

  if(ScalarFields.size())
    output->GetPointData()->ShallowCopy(outputPointData);

  ScalarFields.clear();

  {
    stringstream msg;
    msg << "[ttkFieldSelector] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkFieldSelector::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector){

  Memory m;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkFieldSelector] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
