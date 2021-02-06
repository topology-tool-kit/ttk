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

int addRecursivelyToFieldData(vtkDataObject *object,
                              const vtkSmartPointer<vtkDataArray> &array) {
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
  // Get Input and Output
  auto input = vtkDataObject::GetData(inputVector[0]);

  if(this->LastInput != input || this->IterationIdx >= this->IterationNumber) {
    this->LastInput = input;
    this->IterationIdx = 0;
  }

  // Determine Mode
  auto mode = this->GetExtractionMode();
  if(mode == EXTRACTION_MODE::AUTO) {
    if(input->IsA("vtkMultiBlockDataSet"))
      mode = EXTRACTION_MODE::BLOCKS;
    else if(input->IsA("vtkTable"))
      mode = EXTRACTION_MODE::ROWS;
    else {
      this->printErr("Unable to automatically determine iteration mode.");
      return 0;
    }
  }

  // Get Iteration Bounds
  switch(mode) {
    case EXTRACTION_MODE::BLOCKS: {
      auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast(input);
      if(!inputAsMB) {
        this->printErr(
          "Block iteration requires 'vtkMultiBlockDataSet' input.");
        return 0;
      }
      this->IterationNumber = inputAsMB->GetNumberOfBlocks();
      break;
    }

    case EXTRACTION_MODE::BLOCK_TUPLES: {
      auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast(input);

      if(!inputAsMB || inputAsMB->GetNumberOfBlocks() < 1) {
        this->printErr(
          "Block Tuple iteration requires 'vtkMultiBlockDataSet' input that "
          "contains at least one 'vtkMultiBlockDataSet'.");
        return 0;
      }
      auto firstComponent
        = vtkMultiBlockDataSet::SafeDownCast(inputAsMB->GetBlock(0));
      if(!firstComponent) {
        this->printErr(
          "Block Tuple iteration requires 'vtkMultiBlockDataSet' input that "
          "contains at least one 'vtkMultiBlockDataSet'.");
        return 0;
      }
      this->IterationNumber = firstComponent->GetNumberOfBlocks();
      break;
    }

    case EXTRACTION_MODE::ROWS: {
      auto inputAsT = vtkTable::SafeDownCast(input);
      if(!inputAsT) {
        this->printErr("Row iteration requires 'vtkTable' input.");
        return 0;
      }
      this->IterationNumber = inputAsT->GetNumberOfRows();
      break;
    }

    case EXTRACTION_MODE::ARRAY_VALUES: {
      auto inputArray = this->GetInputArrayToProcess(0, inputVector);
      if(!inputArray) {
        this->printErr("Unable to retrieve input array.");
        return 0;
      }
      this->IterationNumber = inputArray->GetNumberOfTuples();
      break;
    }

    case EXTRACTION_MODE::ARRAYS: {
      auto inputFD
        = input->GetAttributesAsFieldData(this->GetArrayAttributeType());
      if(inputFD) {
        this->IterationNumber = inputFD->GetNumberOfArrays();
      } else {
        this->printErr("Unable to retrieve attribute type.");
        return 0;
      }
      break;
    }

    default: {
      this->printErr("Unsupported mode");
      return 0;
    }
  }

  // Iteration info
  auto iterationInformation = vtkSmartPointer<vtkDoubleArray>::New();
  iterationInformation->SetName("_ttk_IterationInfo");
  iterationInformation->SetNumberOfComponents(2);
  iterationInformation->SetNumberOfTuples(1);
  iterationInformation->SetValue(0, this->IterationIdx);
  iterationInformation->SetValue(1, this->IterationNumber);

  std::string modeStrings[6] = {"B", "R", "G", "V", "A", "BT"};
  this->printMsg("[" + modeStrings[static_cast<int>(mode)] + "] Iteration: ( "
                   + std::to_string(this->IterationIdx) + " / "
                   + std::to_string(this->IterationNumber - 1) + " ) ",
                 ttk::debug::Separator::SLASH);

  this->SetExpressionString(std::to_string(this->IterationIdx));
  this->SetExtractUniqueValues(false);
  this->SetOutputArrayName("IterationArray");

  if(!ttkExtract::RequestData(request, inputVector, outputVector))
    return 0;

  // Add iteration info to output
  auto output = vtkDataObject::GetData(outputVector);
  addRecursivelyToFieldData(output, iterationInformation);

  this->IterationIdx++;

  return 1;
}
