#include <ttkBlockAggregator.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkSmartPointer.h>

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkMultiBlockDataSet.h>

vtkStandardNewMacro(ttkBlockAggregator);

ttkBlockAggregator::ttkBlockAggregator() {
  this->setDebugMsgPrefix("BlockAggregator");

  this->Reset();

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkBlockAggregator::~ttkBlockAggregator() {
}

int ttkBlockAggregator::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

int ttkBlockAggregator::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkBlockAggregator::Reset() {
  this->AggregatedMultiBlockDataSet
    = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  return 1;
}

int copyObjects(vtkDataObject *source, vtkDataObject *copy) {
  if(source->IsA("vtkMultiBlockDataSet")) {
    auto sourceAsMB = vtkMultiBlockDataSet::SafeDownCast(source);
    auto copyAsMB = vtkMultiBlockDataSet::SafeDownCast(copy);

    if(sourceAsMB == nullptr || copyAsMB == nullptr) {
      return 0;
    }

    const auto sourceFD = sourceAsMB->GetFieldData();
    auto copyFD = copyAsMB->GetFieldData();

    if(sourceFD == nullptr || copyFD == nullptr) {
      return 0;
    }

    copyFD->ShallowCopy(sourceFD);

    for(size_t i = 0; i < sourceAsMB->GetNumberOfBlocks(); i++) {
      auto block = sourceAsMB->GetBlock(i);
      auto blockCopy
        = vtkSmartPointer<vtkDataObject>::Take(block->NewInstance());

      copyObjects(block, blockCopy);
      copyAsMB->SetBlock(i, blockCopy);
    }
  } else {
    copy->ShallowCopy(source);
  }

  return 1;
}

int ttkBlockAggregator::AggregateBlock(vtkDataObject *dataObject) {
  ttk::Timer t;
  size_t nBlocks = this->AggregatedMultiBlockDataSet->GetNumberOfBlocks();
  this->printMsg("Adding object add index " + std::to_string(nBlocks), 0,
                 ttk::debug::LineMode::REPLACE);

  auto copy = vtkSmartPointer<vtkDataObject>::Take(dataObject->NewInstance());
  copyObjects(dataObject, copy);

  this->AggregatedMultiBlockDataSet->SetBlock(nBlocks, copy);

  this->printMsg("Adding object at block index " + std::to_string(nBlocks), 1,
                 t.getElapsedTime());

  return 1;
}

int ttkBlockAggregator::RequestData(vtkInformation *ttkNotUsed(request),
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  // Get iteration information
  double iterationIndex = 0;
  this->SetInputArrayToProcess(0, 0, 0, 2, "_ttk_IterationInfo");
  auto iterationInformation = vtkDoubleArray::SafeDownCast(
    this->GetInputArrayToProcess(0, inputVector));
  if(iterationInformation) {
    iterationIndex = iterationInformation->GetValue(0);
    this->AggregatedMultiBlockDataSet->GetFieldData()->AddArray(
      iterationInformation);
  }

  // Check if AggregatedMultiBlockDataSet needs to be reset
  if(!iterationInformation || this->GetForceReset() || iterationIndex == 0)
    this->Reset();

  // Add all inputs
  size_t nInputs = inputVector[0]->GetNumberOfInformationObjects();
  for(size_t i = 0; i < nInputs; i++) {
    auto input = vtkDataObject::GetData(inputVector[0], i);

    if(this->GetFlattenInput() && input->IsA("vtkMultiBlockDataSet")) {
      auto inputAsMB = (vtkMultiBlockDataSet *)input;
      auto nBlocks = inputAsMB->GetNumberOfBlocks();
      for(size_t j = 0; j < nBlocks; j++)
        this->AggregateBlock(inputAsMB->GetBlock(j));
    } else
      this->AggregateBlock(input);
  }

  // Prepare output
  auto output = vtkMultiBlockDataSet::GetData(outputVector);
  output->ShallowCopy(this->AggregatedMultiBlockDataSet);

  return 1;
}
