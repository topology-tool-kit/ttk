#include <ttkTableDataSelector.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkTableDataSelector)

  // transmit abort signals
  bool ttkTableDataSelector::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status
int ttkTableDataSelector::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkTableDataSelector] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkTableDataSelector::doIt(vtkTable *input, vtkTable *output) {
  Memory m;

  output->ShallowCopy(input);

  vtkFieldData *inputRowData = input->GetRowData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputRowData) {
    cerr << "[ttkTableDataSelector] Error: input has no row data." << endl;
    return -1;
  }
#endif

  vtkSmartPointer<vtkFieldData> outputRowData
    = vtkSmartPointer<vtkFieldData>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputRowData) {
    cerr
      << "[ttkTableDataSelector] Error: vtkFieldData memory allocation problem."
      << endl;
    return -1;
  }
#endif

  for(auto &scalar : ScalarFields) {
    if(scalar.length() > 0) {
      vtkDataArray *arr = inputRowData->GetArray(scalar.data());
      if(arr)
        outputRowData->AddArray(arr);
    }
  }

  output->GetRowData()->ShallowCopy(outputRowData);

  {
    stringstream msg;
    msg << "[ttkTableDataSelector] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkTableDataSelector::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {
  Memory m;

  vtkTable *input = vtkTable::GetData(inputVector[0]);
  vtkTable *output = vtkTable::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkTableDataSelector] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
