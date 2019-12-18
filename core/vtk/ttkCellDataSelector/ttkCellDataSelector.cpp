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

  if((ScalarFields.size() == 1) && (RenameSelected)) {
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

  try {
    const int lastField = std::min((int)ScalarFields.size(), RangeId[1] + 1);
    for(int i = RangeId[0]; i < lastField; i++) {
      auto &scalar = ScalarFields[i];
      if(scalar.length() > 0 && regex_match(scalar, regex(RegexpString))) {
        vtkDataArray *arr = inputCellData->GetArray(scalar.data());
        if(arr) {

          if(RenameSelected) {
            if(ScalarFields.size() != 1 && RangeId[1] - RangeId[0] != 0) {
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
  NbScalars = input->GetCellData()->GetNumberOfArrays();
  return vtkDataSetAlgorithm::RequestInformation(request, inputVector, outputVector);
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
