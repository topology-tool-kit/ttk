#include <ttkLDistanceMatrix.h>
#include <LDistance.h>

using namespace std;
using namespace ttk;

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

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkLDistanceMatrix::RequestData(vtkInformation * /*request*/,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  Memory m;

  const auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);

  if(blocks == nullptr) {
    return 0;
  }

  const size_t nInputs{blocks->GetNumberOfBlocks()};

  // Get input data
  std::vector<vtkImageData *> inputData(nInputs);

  for(size_t i = 0; i < nInputs; ++i) {
    inputData[i] = vtkImageData::SafeDownCast(blocks->GetBlock(i));
  }

  // Get output
  auto DistTable = vtkTable::GetData(outputVector);

  std::vector<std::vector<double>> distMatrix(nInputs);

  LDistance worker{};
  worker.setWrapper(this);
  worker.setNumberOfPoints(inputData[0]->GetNumberOfPoints());

  // compute matrix upper triangle
  for(size_t i = 0; i < nInputs; ++i) {
    auto &distCol = distMatrix[i];
    distCol.resize(nInputs);
    // get pointer to scalar field of input i
    const auto inputScalarFieldi
      = inputData[i]->GetPointData()->GetArray(this->ScalarField.data());
    worker.setInputDataPointer1(inputScalarFieldi->GetVoidPointer(0));
    for(size_t j = i + 1; j < nInputs; ++j) {
      // get pointer to scalar field of input jc
      const auto inputScalarFieldj
        = inputData[j]->GetPointData()->GetArray(this->ScalarField.data());
      worker.setInputDataPointer2(inputScalarFieldj->GetVoidPointer(0));
      // call execute
      switch(inputScalarFieldj->GetDataType()) {
        vtkTemplateMacro(worker.execute<VTK_TT>(this->DistanceType));
      }
      // store result
      distMatrix[i][j] = worker.getResult();
    }
  }

  // distance matrix is symmetric
  for(size_t i = 0; i < nInputs; ++i) {
    for(size_t j = i + 1; j < nInputs; ++j) {
      distMatrix[j][i] = distMatrix[i][j];
    }
  }

  // zero-padd column name to keep Row Data columns ordered
  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  // copy diagrams distance matrix to output
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

  {
    stringstream msg;
    msg << "[ttkLDistanceMatrix] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
