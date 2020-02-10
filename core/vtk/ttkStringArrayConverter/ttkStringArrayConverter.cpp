#include <ttkStringArrayConverter.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkDelimitedTextReader.h>
#include <vtkFieldData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkTable.h>

#include <ttkUtils.h>

vtkStandardNewMacro(ttkStringArrayConverter);

ttkStringArrayConverter::ttkStringArrayConverter() {
  this->setDebugMsgPrefix("StringArrayConverter");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkStringArrayConverter::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkStringArrayConverter::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkStringArrayConverter::RequestData(vtkInformation *request,
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {
  ttk::Timer timer;

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

  if(input == nullptr || output == nullptr) {
    this->printErr("Null input data, aborting");
    return 0;
  }

  // point data
  const auto pd = input->GetPointData();
  // string array
  const auto sa = pd->GetAbstractArray(this->InputStringArray.data());

  if(sa == nullptr || !sa->IsA("vtkStringArray")) {
    this->printErr("Cannot find any string array with the name "
                   + this->InputStringArray);
    return 0;
  }

  const auto nvalues = sa->GetNumberOfTuples();
  const auto saIt = sa->NewIterator();

  for(vtkIdType i = 0; i < nvalues; ++i) {
  }

  // print stats
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete (#rows: )", 1, timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
