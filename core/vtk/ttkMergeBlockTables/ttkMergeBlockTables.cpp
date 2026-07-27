#include <ttkMergeBlockTables.h>

#include <vtkFieldData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkTable.h>

vtkStandardNewMacro(ttkMergeBlockTables);

ttkMergeBlockTables::ttkMergeBlockTables() {
  this->setDebugMsgPrefix("MergeBlockTables");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkMergeBlockTables::FillInputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkMergeBlockTables::FillOutputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkMergeBlockTables::RequestData(vtkInformation *ttkNotUsed(request),
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {
  // Get input data
  std::vector<vtkTable *> inputTables;

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);

  // Number of input tables
  size_t nInputs{};

  if(blocks != nullptr) {
    nInputs = blocks->GetNumberOfBlocks();
    for(size_t i = 0; i < nInputs; ++i) {
      inputTables.emplace_back(vtkTable::SafeDownCast(blocks->GetBlock(i)));
    }
  }

  // Sanity check
  for(const auto table : inputTables) {
    if(table == nullptr) {
      this->printErr("Input tables are not all vtkTables");
      return 0;
    }
  }

  // Sanity check
  if(nInputs == 0) {
    this->printErr("No input table");
    return 0;
  }

  // Set output
  auto outputTable = vtkTable::GetData(outputVector);

  // Copy first table to output
  // (deep copy leaves the input untouched)
  outputTable->DeepCopy(inputTables[0]);

  // Clear out the FieldData that may have been copied
  const auto fd = outputTable->GetFieldData();
  if(fd != nullptr) {
    fd->Reset();
  }

  // Insert row of every other table into output
  for(size_t i = 1; i < nInputs; ++i) {
    const auto currTable{inputTables[i]};
    for(vtkIdType j = 0; j < currTable->GetNumberOfRows(); ++j) {
      outputTable->InsertNextRow(currTable->GetRow(j));
    }
  }

  return 1;
}
