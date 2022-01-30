#include <ttkLDistanceMatrix.h>
#include <ttkUtils.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkTable.h>

#include <set>

vtkStandardNewMacro(ttkLDistanceMatrix);

ttkLDistanceMatrix::ttkLDistanceMatrix() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkLDistanceMatrix::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkLDistanceMatrix::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

template <typename T>
int ttkLDistanceMatrix::dispatch(
  std::vector<std::vector<double>> &distanceMatrix,
  const std::vector<vtkDataSet *> &inputData,
  const size_t nPoints) {

  std::vector<const T *> inputPtrs(inputData.size());
  for(size_t i = 0; i < inputData.size(); ++i) {
    inputPtrs[i]
      = ttkUtils::GetPointer<T>(this->GetInputArrayToProcess(0, inputData[i]));
  }
  return this->execute(distanceMatrix, inputPtrs, nPoints);
}

int ttkLDistanceMatrix::RequestData(vtkInformation * /*request*/,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  ttk::Timer tm{};

  const auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);

  if(blocks == nullptr) {
    return 0;
  }

  const size_t nInputs{blocks->GetNumberOfBlocks()};

  // Get input data
  std::vector<vtkDataSet *> inputData(nInputs);

  for(size_t i = 0; i < nInputs; ++i) {
    inputData[i] = vtkDataSet::SafeDownCast(blocks->GetBlock(i));
  }

  // sanity check: data non null and same data size
  {
    std::set<vtkIdType> sizes{};
    for(size_t i = 0; i < nInputs; ++i) {
      if(inputData[i] == nullptr) {
        this->printErr("One input block is not a vtkDataSet");
        return 0;
      }
      sizes.emplace(inputData[i]->GetNumberOfPoints());
    }
    if(sizes.size() > 1) {
      this->printErr("Input blocks do not have the same number of points");
      return 0;
    }
  }

  // Get output
  auto DistTable = vtkTable::GetData(outputVector);

  std::vector<std::vector<double>> distMatrix{};

  const auto firstField = this->GetInputArrayToProcess(0, inputData[0]);
  const auto dataType = firstField->GetDataType();
  const size_t nPoints = firstField->GetNumberOfTuples();

  switch(dataType) {
    vtkTemplateMacro(this->dispatch<VTK_TT>(distMatrix, inputData, nPoints));
  }

  // zero-padd column name to keep Row Data columns ordered
  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  // copy distance matrix to output
  for(size_t i = 0; i < nInputs; ++i) {
    std::string name{"Dataset"};
    zeroPad(name, distMatrix.size(), i);

    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(nInputs);
    col->SetName(name.c_str());
    for(size_t j = 0; j < nInputs; ++j) {
      col->SetTuple1(j, distMatrix[i][j]);
    }
    DistTable->AddColumn(col);
  }

  // aggregate input field data
  vtkNew<vtkFieldData> fd{};
  fd->CopyStructure(inputData[0]->GetFieldData());
  fd->SetNumberOfTuples(inputData.size());
  for(size_t i = 0; i < inputData.size(); ++i) {
    fd->SetTuple(i, 0, inputData[i]->GetFieldData());
  }

  // copy input field data to output row data
  for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
    DistTable->AddColumn(fd->GetAbstractArray(i));
  }

  this->printMsg("Complete (#datasets: " + std::to_string(nInputs)
                   + ", #points: " + std::to_string(nPoints) + ")",
                 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 1;
}
