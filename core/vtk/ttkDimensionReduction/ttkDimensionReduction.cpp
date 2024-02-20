#include <ttkDimensionReduction.h>
#include <ttkUtils.h>

#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkTable.h>

#include <regex>

vtkStandardNewMacro(ttkDimensionReduction);

ttkDimensionReduction::ttkDimensionReduction() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkDimensionReduction::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDimensionReduction::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDimensionReduction::RequestData(vtkInformation *ttkNotUsed(request),
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  using ttk::SimplexId;

  vtkTable *input = vtkTable::GetData(inputVector[0]);
  vtkTable *output = vtkTable::GetData(outputVector);

  // initialize method
  this->setInputMethod(Method);

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

  const SimplexId numberOfRows = input->GetNumberOfRows();
  const SimplexId numberOfColumns = ScalarFields.size();

  if(numberOfRows <= 0 || numberOfColumns <= 0) {
    this->printErr("Input matrix has invalid dimensions (rows: "
                   + std::to_string(numberOfRows)
                   + ", columns: " + std::to_string(numberOfColumns) + ")");
    return 0;
  }

  std::vector<double> inputData;
  std::vector<vtkAbstractArray *> arrays;
  arrays.reserve(ScalarFields.size());
  for(const auto &s : ScalarFields)
    arrays.push_back(input->GetColumnByName(s.data()));
  for(SimplexId i = 0; i < numberOfRows; ++i) {
    for(auto arr : arrays)
      inputData.push_back(arr->GetVariantValue(i).ToDouble());
  }

  outputData_.clear();

  vtkNew<vtkIntArray> insertionTimeCol{};
  int *insertionPtr = nullptr;
  if(this->Method == DimensionReduction::METHOD::TOPOMAP) {
    insertionTimeCol->SetNumberOfTuples(numberOfRows);
    insertionTimeCol->SetName("InsertionTime");
    insertionPtr = ttkUtils::GetPointer<int>(insertionTimeCol);
  }

  const int errorCode = this->execute(
    this->outputData_, inputData, numberOfRows, numberOfColumns, insertionPtr);

  if(!errorCode) {
    if(KeepAllDataArrays)
      output->ShallowCopy(input);

    for(int i = 0; i < NumberOfComponents; ++i) {
      std::string const s = "Component_" + std::to_string(i);
      vtkNew<vtkDoubleArray> arr{};
      ttkUtils::SetVoidArray(arr, outputData_[i].data(), numberOfRows, 1);
      arr->SetName(s.data());
      output->AddColumn(arr);
    }
    if(this->Method == DimensionReduction::METHOD::TOPOMAP) {
      output->AddColumn(insertionTimeCol);
    }
  }

  return 1;
}
