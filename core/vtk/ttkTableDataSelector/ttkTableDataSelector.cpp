#include <ttkTableDataSelector.h>

#include <vtkDataSetAttributes.h>
#include <vtkFieldData.h>

#include <regex>

vtkStandardNewMacro(ttkTableDataSelector);

ttkTableDataSelector::ttkTableDataSelector() {
  this->setDebugMsgPrefix("TableDataSelector");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkTableDataSelector::~ttkTableDataSelector() {
}

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

int ttkTableDataSelector::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  ttk::Timer t;

  this->printMsg("Selecting Columns", 0, 0, ttk::debug::LineMode::REPLACE);

  auto input = vtkTable::GetData(inputVector[0]);
  auto output = vtkTable::GetData(outputVector);

  output->ShallowCopy(input);

  vtkFieldData *inputRowData = input->GetRowData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputRowData) {
    this->printErr("Input has no row data.");
    return -1;
  }
#endif

  auto outputRowData = vtkSmartPointer<vtkFieldData>::New();

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
    if(!std::regex_match(col, std::regex(RegexpString))) {
      continue;
    }
    // add the attay
    vtkDataArray *arr = inputRowData->GetArray(col.c_str());
    if(arr)
      outputRowData->AddArray(arr);
  }

  output->GetRowData()->ShallowCopy(outputRowData);

  this->printMsg("Selecting Columns", 1, t.getElapsedTime());

  return 1;
}

int ttkTableDataSelector::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkTable *input = vtkTable::GetData(inputVector[0]);
  FillAvailableCols(input);
  return ttkAlgorithm::RequestInformation(request, inputVector, outputVector);
}

void ttkTableDataSelector::FillAvailableCols(vtkTable *input) {
  int nbColumns = input->GetNumberOfColumns();
  AvailableCols.clear();
  AvailableCols.resize(nbColumns);
  for(int i = 0; i < nbColumns; ++i) {
    AvailableCols[i] = input->GetColumnName(i);
  }
}