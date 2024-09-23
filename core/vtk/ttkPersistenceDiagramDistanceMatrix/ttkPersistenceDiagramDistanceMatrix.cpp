#include <ttkPersistenceDiagramDistanceMatrix.h>
#include <ttkPersistenceDiagramUtils.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkStringArray.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkPersistenceDiagramDistanceMatrix);

ttkPersistenceDiagramDistanceMatrix::ttkPersistenceDiagramDistanceMatrix() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkPersistenceDiagramDistanceMatrix::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagramDistanceMatrix::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramDistanceMatrix::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  ttk::Memory const m;

  // Get input data
  std::vector<vtkUnstructuredGrid *> inputDiagrams;

  auto nBlocks = inputVector[0]->GetNumberOfInformationObjects();
  std::vector<vtkMultiBlockDataSet *> blocks(nBlocks);

  if(nBlocks > 2) {
    this->printWrn("Only dealing with the first two MultiBlockDataSets");
    nBlocks = 2;
  }

  // number of diagrams per input block
  std::array<size_t, 2> nInputs{0, 0};

  for(int i = 0; i < nBlocks; ++i) {
    blocks[i] = vtkMultiBlockDataSet::GetData(inputVector[0], i);
    if(blocks[i] != nullptr) {
      nInputs[i] = blocks[i]->GetNumberOfBlocks();
      for(size_t j = 0; j < nInputs[i]; ++j) {
        inputDiagrams.emplace_back(
          vtkUnstructuredGrid::SafeDownCast(blocks[i]->GetBlock(j)));
      }
    }
  }

  if(nInputs[0] + nInputs[1] == 0) {
    this->printErr("No input detected");
    return 0;
  }

  // total number of diagrams
  const int nDiags = inputDiagrams.size();

  // Sanity check
  for(const auto vtu : inputDiagrams) {
    if(vtu == nullptr) {
      this->printErr("Input diagrams are not all vtkUnstructuredGrid");
      return 0;
    }
  }

  // Set output
  auto diagramsDistTable = vtkTable::GetData(outputVector);

  std::vector<ttk::DiagramType> intermediateDiagrams(nDiags);

  for(int i = 0; i < nDiags; i++) {
    auto &diag{intermediateDiagrams[i]};
    const auto ret = VTUToDiagram(diag, inputDiagrams[i], *this);
    if(ret < 0) {
      this->printErr("Could not read Persistence Diagram");
      return 0;
    }
  }

  const auto diagramsDistMat = this->execute(intermediateDiagrams, nInputs);

  // zero-padd column name to keep Row Data columns ordered
  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string const max{std::to_string(numberCols - 1)};
        std::string const cur{std::to_string(colIdx)};
        std::string const zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  const auto nTuples = nInputs[1] == 0 ? nInputs[0] : nInputs[1];

  // copy diagrams distance matrix to output
  for(size_t i = 0; i < diagramsDistMat.size(); ++i) {
    std::string name{"Diagram"};
    zeroPad(name, diagramsDistMat.size(), i);

    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(nTuples);
    col->SetName(name.c_str());
    for(size_t j = 0; j < diagramsDistMat[i].size(); ++j) {
      col->SetTuple1(j, diagramsDistMat[i][j]);
    }
    diagramsDistTable->AddColumn(col);
  }

  // aggregate the field data arrays from all input diagrams
  vtkNew<vtkFieldData> inputFdA{};

  for(const auto diag : inputDiagrams) {
    const auto fd{diag->GetFieldData()};
    for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
      const auto array{fd->GetAbstractArray(i)};
      if(array->IsA("vtkDataArray") || array->IsA("vtkStringArray")) {
        inputFdA->AddArray(array);
      }
    }
  }

  // avoid modifying input field data
  vtkNew<vtkFieldData> outputFda{};
  outputFda->DeepCopy(inputFdA);

  for(int i = 0; i < outputFda->GetNumberOfArrays(); ++i) {
    const auto array{outputFda->GetAbstractArray(i)};
    array->SetNumberOfTuples(inputDiagrams.size());
    const auto name{array->GetName()};
    for(size_t j = 0; j < inputDiagrams.size(); ++j) {
      const auto fd{inputDiagrams[j]->GetFieldData()};
      const auto inputArray{fd->GetAbstractArray(name)};
      if(inputArray != nullptr) {
        array->SetTuple(j, 0, inputArray);
      } else {
        if(array->IsA("vtkDataArray")) {
          vtkDataArray::SafeDownCast(array)->SetTuple1(j, NAN);
        } else if(array->IsA("vtkStringArray")) {
          vtkStringArray::SafeDownCast(array)->SetValue(j, "");
        }
      }
    }
    // copy "extended" input field data array to output row data
    diagramsDistTable->AddColumn(array);
  }

  return 1;
}
