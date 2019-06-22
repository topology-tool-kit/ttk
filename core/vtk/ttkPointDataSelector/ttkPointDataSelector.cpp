#include <regex>
#include <ttkPointDataSelector.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPointDataSelector)

  // transmit abort signals
  bool ttkPointDataSelector::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status
int ttkPointDataSelector::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkPointDataSelector] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkPointDataSelector::doIt(vtkDataSet *input, vtkDataSet *output) {
  Memory m;

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

  try {
    for(auto &scalar : ScalarFields) {
      if(scalar.length() > 0 && regex_match(scalar, regex(RegexpString))) {
        vtkDataArray *arr = inputPointData->GetArray(scalar.data());
        if(arr) {

          if((ScalarFields.size() == 1) && (RenameSelected)) {

            if(localFieldCopy_) {
              localFieldCopy_->Delete();
              localFieldCopy_ = NULL;
            }

            switch(arr->GetDataType()) {
              case VTK_CHAR:
                localFieldCopy_ = vtkCharArray::New();
                break;

              case VTK_DOUBLE:
                localFieldCopy_ = vtkDoubleArray::New();
                break;

              case VTK_FLOAT:
                localFieldCopy_ = vtkFloatArray::New();
                break;

              case VTK_INT:
                localFieldCopy_ = vtkIntArray::New();
                break;

              case VTK_ID_TYPE:
                localFieldCopy_ = vtkIdTypeArray::New();
                break;

              case VTK_UNSIGNED_SHORT:
                localFieldCopy_ = vtkUnsignedShortArray::New();
                break;

              default: {
                stringstream msg;
                msg << "[ttkPointDataSelector] Unsupported data type :("
                    << endl;
                dMsg(cerr, msg.str(), fatalMsg);
              } break;
            }

            if(localFieldCopy_) {
              localFieldCopy_->DeepCopy(arr);
              localFieldCopy_->SetName(SelectedFieldName.data());
              arr = localFieldCopy_;
            }
          }

          outputPointData->AddArray(arr);
        }
      }
    }
  } catch(std::regex_error &) {
    vtkWarningMacro("[ttkPointDataSelector]: Bad regexp.");
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
                                      vtkInformationVector *outputVector) {
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
