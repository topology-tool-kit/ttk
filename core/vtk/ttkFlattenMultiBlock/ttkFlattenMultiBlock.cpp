#include <ttkFlattenMultiBlock.h>

#include <vtkFieldData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkUnsignedIntArray.h>

vtkStandardNewMacro(ttkFlattenMultiBlock);

ttkFlattenMultiBlock::ttkFlattenMultiBlock() {
  this->setDebugMsgPrefix("FlattenMultiBlock");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkFlattenMultiBlock::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkFlattenMultiBlock::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkFlattenMultiBlock::RequestData(vtkInformation *ttkNotUsed(request),
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  // Get input data (pointer to vtkDataObject, parent block id)
  std::vector<std::pair<vtkDataObject *, size_t>> blocks{};

  auto input = vtkMultiBlockDataSet::GetData(inputVector[0], 0);

  // Sanity check
  if(input == nullptr || input->GetNumberOfBlocks() == 0) {
    this->printErr("No input block");
    return 0;
  }

  for(size_t i = 0; i < input->GetNumberOfBlocks(); ++i) {
    const auto topLvlBlock
      = vtkMultiBlockDataSet::SafeDownCast(input->GetBlock(i));
    if(topLvlBlock == nullptr) {
      continue;
    }
    for(size_t j = 0; j < topLvlBlock->GetNumberOfBlocks(); ++j) {
      blocks.emplace_back(topLvlBlock->GetBlock(j), i);
    }
  }

  auto outputFlatBlocks = vtkMultiBlockDataSet::GetData(outputVector);

  // Sanity check
  if(blocks.empty()) {
    this->printWrn("Not a hierarchy of vtkMultiBlockDataSet");
    outputFlatBlocks->ShallowCopy(input);
    return 1;
  }

  // Set output
  outputFlatBlocks->SetNumberOfBlocks(blocks.size());

  for(size_t i = 0; i < blocks.size(); ++i) {
    // Set flattened blocks
    outputFlatBlocks->SetBlock(i, blocks[i].first);

    // Set FieldData
    vtkNew<vtkUnsignedIntArray> blockId{};
    blockId->SetName("ParentBlockId");
    blockId->SetNumberOfTuples(1);
    blockId->SetTuple1(0, blocks[i].second);

    blocks[i].first->GetFieldData()->AddArray(blockId);
  }

  return 1;
}
