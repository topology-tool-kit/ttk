#include <ttkTableDataSelector.h>

#include <regex>

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

  if(AvailableCols.empty()) {
    // loading from states
    FillAvailableCols(input);
  }

  for(auto &col : SelectedCols) {
    // check valid col
    if(col.empty()) {
      continue;
    }
    // check bounds in the range
    ptrdiff_t pos = find(AvailableCols.begin(), AvailableCols.end(), col)
                    - AvailableCols.begin();
    if(pos < RangeId[0] || pos > RangeId[1]) {
      continue;
    }
    // filter by regex
    if(!regex_match(col, regex(RegexpString))) {
      continue;
    }
    // add the attay
    vtkDataArray *arr = inputRowData->GetArray(col.c_str());
    if(arr)
      outputRowData->AddArray(arr);
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

int ttkTableDataSelector::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkTable *input = vtkTable::GetData(inputVector[0]);
  FillAvailableCols(input);
  return vtkTableAlgorithm::RequestInformation(
    request, inputVector, outputVector);
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

void ttkTableDataSelector::FillAvailableCols(vtkTable *input) {
  int nbColumns = input->GetNumberOfColumns();
  AvailableCols.clear();
  AvailableCols.resize(nbColumns);
  for(int i = 0; i < nbColumns; ++i) {
    AvailableCols[i] = input->GetColumnName(i);
  }
}
