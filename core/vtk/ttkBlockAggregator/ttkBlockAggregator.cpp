#include <ttkBlockAggregator.h>

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkStreamingDemandDrivenPipeline.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkBlockAggregator)

  int ttkBlockAggregator::AggregateBlock(vtkDataObject *dataObject,
                                         bool useShallowCopy) {

  size_t nBlocks = this->AggregatedMultiBlockDataSet->GetNumberOfBlocks();

  auto temp0 = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  temp0->SetBlock(0, dataObject);

  auto temp1 = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  if(useShallowCopy)
    temp1->ShallowCopy(temp0);
  else
    temp1->DeepCopy(temp0);

  this->AggregatedMultiBlockDataSet->SetBlock(nBlocks, temp1->GetBlock(0));

  // Print status
  stringstream msg;
  msg << "[ttkBlockAggregator] Added block at index " << nBlocks << endl;
  dMsg(cout, msg.str(), infoMsg);

  return 1;
}

int ttkBlockAggregator::RequestData(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  // Print status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl
        << "[ttkBlockAggregator] RequestData" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Get Input
  vtkInformation *inInfoObj = inputVector[0]->GetInformationObject(0);
  auto firstInput = inInfoObj->Get(vtkDataObject::DATA_OBJECT());

  // Get iteration information
  auto iterationInformation = vtkDoubleArray::SafeDownCast(
    firstInput->GetFieldData()->GetAbstractArray("_ttk_IterationInfo"));

  bool useStreamingOverTime = iterationInformation != nullptr;

  double iteration = 0;
  double nIterations = 0;
  if(useStreamingOverTime) {
    iteration = iterationInformation->GetValue(0);
    nIterations = iterationInformation->GetValue(1);
  }

  // First timestep
  if(!useStreamingOverTime || this->GetForceReset() || iteration == 0
     || this->AggregatedMultiBlockDataSet == nullptr)
    this->AggregatedMultiBlockDataSet
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  bool useShallowCopy = this->GetForceReset() && nIterations < 1;

  for(size_t i = 0; i < 5; i++) {
    vtkInformationVector *inVector = inputVector[i];
    if(inVector->GetNumberOfInformationObjects() != 1)
      continue;

    vtkInformation *inInfo = inVector->GetInformationObject(0);
    auto input = inInfo->Get(vtkDataObject::DATA_OBJECT());
    if(input) {

      auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast(input);

      if(this->GetFlattenInput() && inputAsMB) {
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
