#include <ttkBlockAggregator.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkInformationVector.h>

vtkStandardNewMacro(ttkBlockAggregator);

ttkBlockAggregator::ttkBlockAggregator() {
  this->setDebugMsgPrefix("BlockAggregator");

  this->Reset();

  this->SetNumberOfInputPorts(5);
  this->SetNumberOfOutputPorts(1);
}
ttkBlockAggregator::~ttkBlockAggregator() {
}

int ttkBlockAggregator::Reset() {
  this->AggregatedMultiBlockDataSet
    = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  return 1;
}

int ttkBlockAggregator::AggregateBlock(vtkDataObject *dataObject,
                                       bool useShallowCopy) {
  ttk::Timer t;
  size_t nBlocks = this->AggregatedMultiBlockDataSet->GetNumberOfBlocks();
  this->printMsg("Adding object add index " + std::to_string(nBlocks), 0,
                 ttk::debug::LineMode::REPLACE);

  auto temp = vtkSmartPointer<vtkDataObject>::Take(dataObject->NewInstance());
  useShallowCopy ? temp->ShallowCopy(dataObject) : temp->DeepCopy(dataObject);

  this->AggregatedMultiBlockDataSet->SetBlock(nBlocks, temp);
  this->printMsg("Adding object at block index " + std::to_string(nBlocks), 1,
                 t.getElapsedTime());

  return 1;
}

int ttkBlockAggregator::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  if(port > 0)
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  return 1;
}

int ttkBlockAggregator::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
}

int ttkBlockAggregator::RequestData(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  // Get Input
  auto firstInput = vtkDataObject::GetData(inputVector[0]);

  // Get iteration information
  double iteration = 0;
  auto iterationInformation = vtkDoubleArray::SafeDownCast(
    firstInput->GetFieldData()->GetAbstractArray("_ttk_IterationInfo"));
  bool useIterations = iterationInformation != nullptr;
  if(useIterations) {
    iteration = iterationInformation->GetValue(0);
    this->AggregatedMultiBlockDataSet->GetFieldData()->AddArray(
      iterationInformation);
  }

  // Check if AggregatedMultiBlockDataSet needs to be reset
  if(!useIterations || this->GetForceReset() || iteration == 0)
    this->Reset();

  // useShallowCopy only if there are no iterations
  bool useShallowCopy = !useIterations;

  // Iterate over input ports
  for(size_t i = 0; i < 5; i++) {
    vtkInformationVector *inVector = inputVector[i];
    if(inVector->GetNumberOfInformationObjects() != 1)
      continue;

    vtkInformation *inInfo = inVector->GetInformationObject(0);
    auto input = inInfo->Get(vtkDataObject::DATA_OBJECT());
    if(input) {
      auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast(input);

      if(inputAsMB && this->GetFlattenInput()) {
        auto nBlocks = inputAsMB->GetNumberOfBlocks();
        for(size_t j = 0; j < nBlocks; j++)
          this->AggregateBlock(inputAsMB->GetBlock(j), useShallowCopy);
      } else
        this->AggregateBlock(input, useShallowCopy);
    }
  }

  // Get Output
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto outMB = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  outMB->ShallowCopy(this->AggregatedMultiBlockDataSet);

  return 1;
}