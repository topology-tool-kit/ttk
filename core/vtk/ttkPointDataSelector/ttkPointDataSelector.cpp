#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <algorithm>
#include <regex>
#include <ttkPointDataSelector.h>

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

void ttkPointDataSelector::FillAvailableFields(vtkDataSet *input) {
  int nbScalars = input->GetPointData()->GetNumberOfArrays();
  AvailableFields.clear();
  AvailableFields.resize(nbScalars);
  for(int i = 0; i < nbScalars; ++i) {
    AvailableFields[i] = input->GetPointData()->GetArrayName(i);
  }
}

int ttkPointDataSelector::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  FillAvailableFields(input);
  return ttkAlgorithm::RequestInformation(request, inputVector, outputVector);
}

int ttkPointDataSelector::RequestData(vtkInformation *ttkNotUsed(request),
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  output->ShallowCopy(input);

  vtkPointData *inputPointData = input->GetPointData();
  vtkNew<vtkPointData> outputPointData{};

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputPointData) {
    this->printErr("Input has no point data.");
    return -1;
  }

  if(!outputPointData) {
    this->printErr("vtkPointData memory allocation problem.");
    return -1;
  }
#endif

  if(AvailableFields.empty()) {
    // when loading from statefiles or vtk script
    FillAvailableFields(input);
  }

  try {
    for(auto &scalar : SelectedFields) {
      // valid array
      if(scalar.empty()) {
        continue;
      }
      // check bounds in the range
      std::ptrdiff_t pos
        = std::find(AvailableFields.begin(), AvailableFields.end(), scalar)
          - AvailableFields.begin();
      if(pos < RangeId[0] || pos > RangeId[1]) {
        continue;
      }
      // retrieve array if match
      if(!std::regex_match(scalar, std::regex(RegexpString))) {
        continue;
      }
      // Add the array
      vtkDataArray *arr = inputPointData->GetArray(scalar.c_str());
      if(arr) {

        if(RenameSelected) {
          if(SelectedFields.size() != 1 && RangeId[1] - RangeId[0] != 0) {
            this->printErr("Can't rename more than one field.");
            return 0;
          }

          vtkSmartPointer<vtkDataArray> localFieldCopy{arr->NewInstance()};
          if(localFieldCopy) {
            localFieldCopy->DeepCopy(arr);
            localFieldCopy->SetName(SelectedFieldName.data());
            arr = localFieldCopy;
          }
        }

        outputPointData->AddArray(arr);
      }
    }
  } catch(std::regex_error &) {
    this->printErr("Bad regexp.");
  }

  output->GetPointData()->ShallowCopy(outputPointData);

  return 1;
}
