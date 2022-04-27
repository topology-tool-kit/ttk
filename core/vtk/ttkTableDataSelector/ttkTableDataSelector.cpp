#include <ttkMacros.h>
#include <ttkTableDataSelector.h>
#include <ttkUtils.h>

#include <vtkObjectFactory.h>

#include <regex>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkTableDataSelector);

int ttkTableDataSelector::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkTableDataSelector::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkTableDataSelector::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkTable *input = vtkTable::GetData(inputVector[0]);
  FillAvailableCols(input);
  return ttkAlgorithm::RequestInformation(request, inputVector, outputVector);
}

int ttkTableDataSelector::RequestData(vtkInformation *ttkNotUsed(request),
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {
  vtkTable *input = vtkTable::GetData(inputVector[0]);
  vtkTable *output = vtkTable::GetData(outputVector);

  Timer timer;

  output->ShallowCopy(input);

  vtkDataSetAttributes *inputRowData = input->GetRowData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputRowData) {
    this->printErr("Input has no row data.");
    return -1;
  }
#endif

  vtkSmartPointer<vtkDataSetAttributes> outputRowData
    = vtkSmartPointer<vtkDataSetAttributes>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputRowData) {
    this->printErr("vtkFieldData memory allocation problem.");
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
    // add the array
    vtkDataArray *arr = inputRowData->GetArray(col.c_str());
    if(arr)
      outputRowData->AddArray(arr);
  }

  output->GetRowData()->ShallowCopy(outputRowData);

  {
    this->printMsg(ttk::debug::Separator::L1);
    this->printMsg("Complete", 1, timer.getElapsedTime());
    this->printMsg(ttk::debug::Separator::L1);
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
