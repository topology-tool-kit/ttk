#include <algorithm>
#include <cstddef>
#include <regex>
#include <ttkPointDataSelector.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPointDataSelector);

ttkPointDataSelector::ttkPointDataSelector() {
  this->setDebugMsgPrefix("PointDataSelector");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkPointDataSelector::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPointDataSelector::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkPointDataSelector::doIt(vtkDataSet *input, vtkDataSet *output) {
  output->ShallowCopy(input);

  vtkPointData *inputPointData = input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputPointData) {
    cerr << "[ttkPointDataSelector] Error: input has no point data." << endl;
    return -1;
  }
#endif

  vtkSmartPointer<vtkPointData> outputPointData
    = vtkSmartPointer<vtkPointData>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputPointData) {
    cerr
      << "[ttkPointDataSelector] Error: vtkPointData memory allocation problem."
      << endl;
    return -1;
  }
#endif

  if(AvailableFields.empty()) {
    // when loading from statefiles
    // or vtk script
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
      // retrieve array if match
      if(!regex_match(scalar, regex(RegexpString))) {
        continue;
      }
      // Add the array
      vtkDataArray *arr = inputPointData->GetArray(scalar.c_str());
      if(arr) {

        if(RenameSelected) {
          if(SelectedFields.size() != 1 && RangeId[1] - RangeId[0] != 0) {
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

        outputPointData->AddArray(arr);
      }
    }
  } catch(std::regex_error &) {
    vtkWarningMacro("[ttkPointDataSelector]: Bad regexp.");
  }

  output->GetPointData()->ShallowCopy(outputPointData);

  return 0;
}

int ttkPointDataSelector::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  FillAvailableFields(input);
  return ttkAlgorithm::RequestInformation(request, inputVector, outputVector);
}

int ttkPointDataSelector::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  doIt(input, output);

  return 1;
}

void ttkPointDataSelector::FillAvailableFields(vtkDataSet *input) {
  int nbScalars = input->GetPointData()->GetNumberOfArrays();
  AvailableFields.clear();
  AvailableFields.resize(nbScalars);
  for(int i = 0; i < nbScalars; ++i) {
    AvailableFields[i] = input->GetPointData()->GetArrayName(i);
  }
}
