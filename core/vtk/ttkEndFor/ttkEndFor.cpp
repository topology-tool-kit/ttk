#include <ttkEndFor.h>

#include <vtkCompositeDataPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>

#include <ttkForEach.h>

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

int removeFieldDataRecursively(vtkDataObject *object) {
  object->GetFieldData()->RemoveArray("_ttk_IterationInfo");
  if(object->IsA("vtkMultiBlockDataSet")) {
    auto objectAsMB = static_cast<vtkMultiBlockDataSet *>(object);
    for(size_t i = 0; i < objectAsMB->GetNumberOfBlocks(); i++)
      removeFieldDataRecursively(objectAsMB->GetBlock(i));
  }
  return 1;
}

int ttkEndFor::RequestData(vtkInformation *request,
                           vtkInformationVector **inputVector,
                           vtkInformationVector *outputVector) {

  // find for each head
  ttkForEach *forEach = nullptr;
  {
    auto inputAlgorithm = this->GetInputAlgorithm(1, 0);
    while(inputAlgorithm && !inputAlgorithm->IsA("ttkForEach")) {
      inputAlgorithm = inputAlgorithm->GetInputAlgorithm();
    }
    forEach = ttkForEach::SafeDownCast(inputAlgorithm);
  }

  if(!forEach) {
    this->printErr("Second input not connected to a ttkForEach filter.");
    return 0;
  }

  // get iteration info
  int i = forEach->GetIterationIdx() - 1;
  int n = forEach->GetIterationNumber();

  bool isRepeatedIteration = this->LastIterationIdx == i && i > 0;
  this->LastIterationIdx = i;

  if(isRepeatedIteration)
    this->printMsg("For Loop Modified -> Restarting Iterations",
                   ttk::debug::Separator::BACKSLASH);
  else
    this->printMsg("Iteration ( " + std::to_string(i) + " / "
                     + std::to_string(n - 1) + " ) complete ",
                   ttk::debug::Separator::BACKSLASH);

  if(i >= n - 1 && !isRepeatedIteration) {
    // if this is the last iteration
    auto input = vtkDataObject::GetData(inputVector[0]);
    auto output = vtkDataObject::GetData(outputVector);
    output->ShallowCopy(input);
    removeFieldDataRecursively(output);
    request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
  } else {
    // if this is an intermediate iteration
    forEach->Modified();
    this->GetInputAlgorithm(0, 0)->Update();

    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
  }

  return 1;
}