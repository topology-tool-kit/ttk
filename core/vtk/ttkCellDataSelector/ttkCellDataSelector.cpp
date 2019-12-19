#include <regex>

#include <ttkCellDataSelector.h>

#include <vtkPointData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCellDataSelector)

  // transmit abort signals
  bool ttkCellDataSelector::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status
int ttkCellDataSelector::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkCellDataSelector] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkCellDataSelector::doIt(vtkDataSet *input, vtkDataSet *output) {
  Memory m;

  if((SelectedFieldName.size() == 1) && (RenameSelected)) {
    output->DeepCopy(input);
  } else {
    output->ShallowCopy(input);
  }

  vtkCellData *inputCellData = input->GetCellData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputCellData) {
    cerr << "[ttkCellDataSelector] Error: input has no cell data." << endl;
    return -1;
  }
#endif

  vtkSmartPointer<vtkCellData> outputCellData
    = vtkSmartPointer<vtkCellData>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputCellData) {
    cerr
      << "[ttkCellDataSelector] Error: vtkCellData memory allocation problem."
      << endl;
    return -1;
  }
#endif

  if(AvailableFields.empty()) {
    // when loading states
    FillAvailableFields(input);
  }

  try {
    for(auto &scalar : SelectedFields) {
      // valid array
      if(scalar.empty()) {
        continue;
      }
      // check bounds in the range
      ptrdiff_t pos
        = find(AvailableFields.begin(), AvailableFields.end(), scalar)
          - AvailableFields.begin();
      if(pos < RangeId[0] || pos > RangeId[1]) {
        continue;
      }
      // filter by regex
      if(!regex_match(scalar, regex(RegexpString))) {
        continue;
      }
      // add the attay
      vtkDataArray *arr = inputCellData->GetArray(scalar.data());
      if(arr) {

        if(RenameSelected) {
          if(SelectedFieldName.size() != 1 && RangeId[1] - RangeId[0] != 0) {
            vtkErrorMacro("Can't rename more than one field.");
            return 0;
          }

          if(localFieldCopy_) {
            localFieldCopy_->Delete();
            localFieldCopy_ = nullptr;
          }

          localFieldCopy_ = arr->NewInstance();

          if(localFieldCopy_) {
            localFieldCopy_->DeepCopy(arr);
            localFieldCopy_->SetName(SelectedFieldName.data());
            arr = localFieldCopy_;
          }
        }

        outputCellData->AddArray(arr);
      }
    }
  } catch(std::regex_error &) {
    vtkWarningMacro("[ttkCellDataSelector]: Bad regexp.");
  }

  output->GetCellData()->ShallowCopy(outputCellData);

  {
    stringstream msg;
    msg << "[ttkCellDataSelector] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkCellDataSelector::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  FillAvailableFields(input);
  return vtkDataSetAlgorithm::RequestInformation(
    request, inputVector, outputVector);
}

int ttkCellDataSelector::RequestData(vtkInformation *request,
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {

  Memory m;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkCellDataSelector] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}

void ttkCellDataSelector::FillAvailableFields(vtkDataSet *input) {
  int nbScalars = input->GetCellData()->GetNumberOfArrays();
  AvailableFields.clear();
  AvailableFields.resize(nbScalars);
  for(int i = 0; i < nbScalars; ++i) {
    AvailableFields[i] = input->GetCellData()->GetArrayName(i);
  }
}
