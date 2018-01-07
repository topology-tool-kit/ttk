#include <ttkPointDataSelector.h>

vtkStandardNewMacro(ttkPointDataSelector)

  // transmit abort signals
  bool ttkPointDataSelector::needsToAbort(){
    return GetAbortExecute();
  }

// transmit progress status
int ttkPointDataSelector::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkPointDataSelector] " << progress*100 
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkPointDataSelector::doIt(vtkDataSet* input, vtkDataSet* output){
  Memory m;

  output->ShallowCopy(input);

  vtkPointData* inputPointData=input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputPointData){
    cerr << "[ttkPointDataSelector] Error: input has no point data." << endl;
    return -1;
  }
#endif

  vtkSmartPointer<vtkPointData> outputPointData=vtkSmartPointer<vtkPointData>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputPointData){
    cerr << "[ttkPointDataSelector] Error: vtkPointData memory allocation problem." << endl;
    return -1;
  }
#endif

  for(auto& scalar : ScalarFields){
    if(scalar.length()>0){
      vtkDataArray* arr=inputPointData->GetArray(scalar.data());
      if(arr) outputPointData->AddArray(arr);
    }
  }

  output->GetPointData()->ShallowCopy(outputPointData);

  {
    stringstream msg;
    msg << "[ttkPointDataSelector] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkPointDataSelector::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector){
  Memory m;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkPointDataSelector] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
