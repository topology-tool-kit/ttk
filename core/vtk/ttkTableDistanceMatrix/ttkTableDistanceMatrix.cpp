#include <ttkTableDistanceMatrix.h>
#include <ttkUtils.h>

#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkTable.h>

#include <regex>

vtkStandardNewMacro(ttkTableDistanceMatrix);

ttkTableDistanceMatrix::ttkTableDistanceMatrix() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkTableDistanceMatrix::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkTableDistanceMatrix::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkTableDistanceMatrix::RequestData(vtkInformation * /*request*/,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  ttk::Timer tm{};

  auto *input = vtkTable::GetData(inputVector[0]);
  auto *output = vtkTable::GetData(outputVector);
  output->ShallowCopy(input);

  if(SelectFieldsWithRegexp) {
    // select all input columns whose name is matching the regexp
    ScalarFields.clear();
    const auto n = input->GetNumberOfColumns();
    for(int i = 0; i < n; ++i) {
      const auto &name = input->GetColumnName(i);
      if(std::regex_match(name, std::regex(RegexpString))) {
        ScalarFields.emplace_back(name);
      }
    }
  }

  const auto numberOfRows = input->GetNumberOfRows();
  const auto numberOfColumns = ScalarFields.size();

  if(numberOfRows <= 0 || numberOfColumns <= 0) {
    this->printErr("Input matrix has invalid dimensions (rows: "
                   + std::to_string(numberOfRows)
                   + ", columns: " + std::to_string(numberOfColumns) + ")");
    return 0;
  }

  std::vector<vtkAbstractArray *> arrays{};
  for(const auto &s : ScalarFields) {
    arrays.push_back(input->GetColumnByName(s.data()));
  }

  std::vector<std::vector<double>> inputMatrix(numberOfRows);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < numberOfRows; ++i) {
    for(size_t j = 0; j < numberOfColumns; ++j) {
      inputMatrix[i].emplace_back(arrays[j]->GetVariantValue(i).ToDouble());
    }
  }

  std::vector<const double *> inputPtrs(inputMatrix.size());
  for(size_t i = 0; i < inputMatrix.size(); ++i) {
    const auto &vec{inputMatrix[i]};
    inputPtrs[i] = vec.data();
  }

  std::vector<std::vector<double>> distanceMatrix{};
  this->execute(distanceMatrix, inputPtrs, inputMatrix[0].size());

  // zero-padd column name to keep Row Data columns ordered
  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  // copy distance matrix to output
  for(int i = 0; i < numberOfRows; ++i) {
    std::string name{"Distance"};
    zeroPad(name, distanceMatrix.size(), i);

    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(numberOfRows);
    col->SetName(name.c_str());
    for(int j = 0; j < numberOfRows; ++j) {
      col->SetTuple1(j, distanceMatrix[i][j]);
    }
    output->AddColumn(col);
  }

  this->printMsg("Complete (#dimensions: " + std::to_string(numberOfColumns)
                   + ", #points: " + std::to_string(numberOfRows) + ")",
                 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 1;
}
