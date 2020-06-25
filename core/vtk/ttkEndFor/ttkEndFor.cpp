#include <ttkEndFor.h>

#include <vtkCompositeDataPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>

vtkStandardNewMacro(ttkEndFor);

ttkEndFor::ttkEndFor() {
  this->setDebugMsgPrefix("EndFor");

  SetNumberOfInputPorts(2);
  SetNumberOfOutputPorts(1);
}

ttkEndFor::~ttkEndFor(){};

int ttkEndFor::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject", 1);
    return 1;
  }
  return 0;
}

int ttkEndFor::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkEndFor::RequestUpdateExtent(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  // Request next index for data input
  inputVector[0]->GetInformationObject(0)->Set(
    vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->Iteration);

  return 1;
}

int ttkEndFor::RequestData(vtkInformation *request,
                           vtkInformationVector **inputVector,
                           vtkInformationVector *outputVector) {
  // Get Data Input
  auto inputData = vtkDataObject::GetData(inputVector[0]);

  // Get For Input
  auto inputFor = vtkDataObject::GetData(inputVector[1]);

  // Get iteration information from For
  auto iterationInformation = vtkDoubleArray::SafeDownCast(
    inputFor->GetFieldData()->GetAbstractArray("_ttk_IterationInfo"));
  if(!iterationInformation) {
    this->printErr(
      "Unable to retrieve iteration information from ForEach head");
    return 0;
  }
  this->Iteration = iterationInformation->GetValue(0) + 1;
  int nIterations = iterationInformation->GetValue(1);

  // Print status
  this->printMsg("Iteration ( " + std::to_string(this->Iteration - 1) + " / "
                   + std::to_string(nIterations - 1) + " ) complete ",
                 ttk::debug::Separator::BACKSLASH);

  if(this->Iteration < nIterations) {
    // Request Next Element
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
  } else {
    // Stop iterations
    request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
    this->Iteration = 0;

    // Copy Input to Output
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    auto output = outInfo->Get(vtkDataObject::DATA_OBJECT());
    output->ShallowCopy(inputData);
    output->GetFieldData()->RemoveArray("_ttk_IterationInfo");
  }

  return 1;
}