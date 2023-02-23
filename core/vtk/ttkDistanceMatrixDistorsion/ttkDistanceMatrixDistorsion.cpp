#include <ttkDistanceMatrixDistorsion.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkVariantArray.h>

#include <regex>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkDistanceMatrixDistorsion);

ttkDistanceMatrixDistorsion::ttkDistanceMatrixDistorsion() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int ttkDistanceMatrixDistorsion::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDistanceMatrixDistorsion::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

inline void fillWithInputColumns(vtkTable *input,
                                 const std::string &regExpStr,
                                 std::vector<std::string> &vectToFill) {
  // Select all input columns whose name is matching the regexp.
  vectToFill.clear();
  const size_t n = input->GetNumberOfColumns();
  for(size_t i = 0; i < n; ++i) {
    const auto &name = input->GetColumnName(i);
    if(std::regex_match(name, std::regex(regExpStr))) {
      vectToFill.emplace_back(name);
    }
  }
}

int ttkDistanceMatrixDistorsion::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  // Get input objects from input vector
  vtkTable *inputHigh = vtkTable::GetData(inputVector[0]);
  vtkTable *inputLow = vtkTable::GetData(inputVector[1]);

  vtkTable *output = vtkTable::GetData(outputVector);

  if(!inputLow || !inputHigh || !output)
    return 0;

  // Selecting the columns names if a regexp was given.
  if(SelectFieldsWithRegexpHigh)
    fillWithInputColumns(inputHigh, RegexpStringHigh, ScalarFieldsHigh);
  if(SelectFieldsWithRegexpLow)
    fillWithInputColumns(inputLow, RegexpStringLow, ScalarFieldsLow);

  const size_t nRowsHigh = inputHigh->GetNumberOfRows();
  const size_t nColsHigh = ScalarFieldsHigh.size();
  const size_t nRowsLow = inputLow->GetNumberOfRows();
  const size_t nColsLow = ScalarFieldsLow.size();
  if(nRowsHigh == 0 || nRowsHigh != nColsHigh) {
    this->printErr("High input matrix is not a valid square matrix (rows: "
                   + std::to_string(nRowsHigh)
                   + ", columns: " + std::to_string(nColsHigh) + ")");
    return 0;
  }

  if(nRowsLow == 0 || nRowsLow != nColsLow) {
    this->printErr("Low input matrix is not a valid square matrix (rows: "
                   + std::to_string(nRowsLow)
                   + ", columns: " + std::to_string(nColsLow) + ")");
    return 0;
  }
  if(nRowsHigh != nRowsLow) {
    this->printErr(
      "High and low input matrices must have same size (rows(high): "
      + std::to_string(nRowsHigh) + ", rows(low): " + std::to_string(nRowsLow)
      + ")");
    return 0;
  }

  int n = nRowsHigh;
  std::vector<double *> vectMatHigh(n),
    vectMatLow(n); // No 2D vectors to avoid copy of data from the VTK layer.

  // Getting the actual input columns.
  std::vector<vtkDataArray *> arraysHigh{}, arraysLow{};
  for(const auto &s : ScalarFieldsHigh)
    arraysHigh.push_back(
      vtkDoubleArray::SafeDownCast(inputHigh->GetColumnByName(s.data())));
  for(const auto &s : ScalarFieldsLow)
    arraysLow.push_back(
      vtkDoubleArray::SafeDownCast(inputLow->GetColumnByName(s.data())));

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < n; i++) {
    vectMatHigh[i] = ttkUtils::GetPointer<double>(arraysHigh[i]);
    vectMatLow[i] = ttkUtils::GetPointer<double>(arraysLow[i]);
  }

  double distorsionValue = 0;
  vtkNew<vtkDoubleArray> tmpCol{}, distorsionValArray{};
  tmpCol->SetNumberOfTuples(n);
  this->printMsg("Starting computation of sim distorsion value...");
  this->execute(vectMatHigh, vectMatLow, distorsionValue,
                ttkUtils::GetPointer<double>(tmpCol));

  tmpCol->SetName("SimValue");
  // No deep copy, makes output->RowData points to the data of tmpCol.
  output->AddColumn(tmpCol);
  distorsionValArray->SetName("DistorsionValue");
  distorsionValArray->SetNumberOfTuples(1);
  distorsionValArray->SetTuple1(0, distorsionValue);
  output->GetFieldData()->AddArray(distorsionValArray);

  return 1;
}
