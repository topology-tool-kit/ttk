#include <ttkEndFor.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>

vtkStandardNewMacro(ttkEndFor);

ttkEndFor::ttkEndFor() {
  this->setDebugMsgPrefix("EndFor");

  SetNumberOfInputPorts(2);
  SetNumberOfOutputPorts(1);
}

ttkEndFor::~ttkEndFor(){};

int ttkEndFor::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  else
    return 0;
  return 1;
}

int ttkEndFor::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;
  return 1;
}

int ttkEndFor::RequestInformation(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  // Reset index
  this->nextIndex = 0;
  return 1;
}

int ttkEndFor::RequestUpdateExtent(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  // Request next index
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  inInfo->Set(
    vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->nextIndex);

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
  auto iterationInformationFor = vtkDoubleArray::SafeDownCast(
    inputFor->GetFieldData()->GetAbstractArray("_ttk_IterationInfo"));
  if(!iterationInformationFor) {
    this->printErr(
      "Unable to retrieve iteration information from ForEach head");
    return 0;
  }
  this->nextIndex = iterationInformationFor->GetValue(0) + 1;
  int nIteration = iterationInformationFor->GetValue(1);

  if(this->nextIndex < nIteration) {
    // Request Next Element
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);

    // Print status
    this->printMsg("Next Iteration: ( " + std::to_string(this->nextIndex + 1)
                     + " / " + std::to_string(nIteration) + " ) ",
                   ttk::debug::Separator::BACKSLASH);
  } else {
    // Stop iterations
    request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
    this->nextIndex = 0;

    // Print status
    this->printMsg("Iterations complete ", ttk::debug::Separator::BACKSLASH);

    // Copy Input to Output
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    auto output = outInfo->Get(vtkDataObject::DATA_OBJECT());
    output->ShallowCopy(inputData);
    output->GetFieldData()->RemoveArray("_ttk_IterationInfo");
  }

  return 1;
}