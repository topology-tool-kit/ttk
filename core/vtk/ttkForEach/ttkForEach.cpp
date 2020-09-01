#include <ttkForEach.h>

#include <vtkDataSet.h>
#include <vtkInformation.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTable.h>

vtkStandardNewMacro(ttkForEach);

ttkForEach::ttkForEach() {
  this->setDebugMsgPrefix("ForEach");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkForEach::~ttkForEach(){};

int ttkForEach::RequestInformation(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  // These values need to exists to automatically enable temporal streaming
  double dummy[2] = {-99999, -99998};
  outputVector->GetInformationObject(0)->Set(
    vtkStreamingDemandDrivenPipeline::TIME_STEPS(), dummy, 2);
  outputVector->GetInformationObject(0)->Set(
    vtkStreamingDemandDrivenPipeline::TIME_RANGE(), dummy, 2);

  return 1;
}

int addRecursivelyToFieldData(vtkDataObject *object,
                              vtkSmartPointer<vtkDataArray> array) {
  object->GetFieldData()->AddArray(array);
  if(object->IsA("vtkMultiBlockDataSet")) {
    auto objectAsMB = (vtkMultiBlockDataSet *)object;
    for(size_t i = 0; i < objectAsMB->GetNumberOfBlocks(); i++)
      addRecursivelyToFieldData(objectAsMB->GetBlock(i), array);
  }
  return 1;
}

int ttkForEach::RequestData(vtkInformation *request,
                            vtkInformationVector **inputVector,
                            vtkInformationVector *outputVector) {
  // retrieve iterationIndex from pipeline
  auto iterationIndex = outputVector->GetInformationObject(0)->Get(
    vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

  // Get Input and Output
  auto input = vtkDataObject::GetData(inputVector[0]);

  // Determine Mode
  std::string modeStrings[5] = {"B", "R", "G", "V", "A"};

  // Iteration info
  auto iterationInformation = vtkSmartPointer<vtkDoubleArray>::New();
  iterationInformation->SetName("_ttk_IterationInfo");
  iterationInformation->SetNumberOfComponents(2);
  iterationInformation->SetNumberOfTuples(1);
  iterationInformation->SetValue(0, iterationIndex);
  iterationInformation->SetValue(1, 0);

  auto mode = this->GetExtractionMode();
  if(mode < 0) {
    if(input->IsA("vtkMultiBlockDataSet"))
      mode = 0;
    else if(input->IsA("vtkTable"))
      mode = 1;
    else {
      this->printErr("Unable to automatically determine iteration mode.");
      return 0;
    }
  }

  // Get Iteration Bounds
  if(mode == 0) {
    if(!input->IsA("vtkMultiBlockDataSet")) {
      this->printErr("Block iteration requires 'vtkMultiBlockDataSet' input.");
      return 0;
    }
    iterationInformation->SetValue(
      1, ((vtkMultiBlockDataSet *)input)->GetNumberOfBlocks());
  } else if(mode == 1) {
    if(!input->IsA("vtkTable")) {
      this->printErr("Row iteration requires 'vtkTable' input.");
      return 0;
    }
    iterationInformation->SetValue(1, ((vtkTable *)input)->GetNumberOfRows());
  } else if(mode == 3) {
    auto inputArray = this->GetInputArrayToProcess(0, inputVector);
    if(!inputArray) {
      this->printErr("Unable to retrieve input array.");
      return 0;
    }
    iterationInformation->SetValue(1, inputArray->GetNumberOfTuples());
  } else if(mode == 4) {
    auto inputFD
      = input->GetAttributesAsFieldData(this->GetArrayAttributeType());
    if(inputFD == nullptr) {
      this->printErr("Unable to retrieve attribute type.");
      return 0;
    } else {
      iterationInformation->SetValue(1, inputFD->GetNumberOfArrays());
    }
  } else {
    this->printErr("Unsupported mode");
    return 0;
  }

  this->printMsg(
    "[" + modeStrings[mode] + "] Iteration: ( "
      + std::to_string((int)(iterationIndex)) + " / "
      + std::to_string((int)(iterationInformation->GetValue(1)) - 1) + " ) ",
    ttk::debug::Separator::SLASH);

  this->SetExpressionString(std::to_string((int)iterationIndex));
  this->SetExtractUniqueValues(false);
  this->SetOutputArrayName("IterationArray");

  if(!ttkExtract::RequestData(request, inputVector, outputVector))
    return 0;

  // Add iteration info to output
  auto output = vtkDataObject::GetData(outputVector);
  addRecursivelyToFieldData(output, iterationInformation);

  return 1;
}
