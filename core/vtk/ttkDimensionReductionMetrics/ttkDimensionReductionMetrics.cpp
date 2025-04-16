#include <ttkDimensionReductionMetrics.h>

#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkStringArray.h>
#include <vtkTable.h>

#include <regex>

vtkStandardNewMacro(ttkDimensionReductionMetrics);

ttkDimensionReductionMetrics::ttkDimensionReductionMetrics() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int ttkDimensionReductionMetrics::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDimensionReductionMetrics::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDimensionReductionMetrics::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkTable *input = vtkTable::GetData(inputVector[0]);
  vtkTable *representation = vtkTable::GetData(inputVector[1]);
  if(!input || !representation)
    return 0;

  vtkTable *output = vtkTable::GetData(outputVector, 0);

  /*** input points ***/
  if(SelectInputFieldsWithRegexp) {
    // select all input columns whose name is matching the regexp
    InputScalarFields.clear();
    const auto n = input->GetNumberOfColumns();
    for(int i = 0; i < n; ++i) {
      const auto &name = input->GetColumnName(i);
      if(std::regex_match(name, std::regex(InputRegexpString))) {
        InputScalarFields.emplace_back(name);
      }
    }
  }
  if(input->GetNumberOfRows() <= 0 || InputScalarFields.size() <= 0) {
    this->printErr("Input matrix has invalid dimensions (rows: "
                   + std::to_string(input->GetNumberOfRows()) + ", columns: "
                   + std::to_string(InputScalarFields.size()) + ")");
    return 0;
  }
  std::vector<vtkAbstractArray *> inputArrays;
  inputArrays.reserve(InputScalarFields.size());
  for(const auto &s : InputScalarFields)
    inputArrays.push_back(input->GetColumnByName(s.data()));

  const int numberOfPoints = input->GetNumberOfRows();
  const int dimensionHigh = InputScalarFields.size();
  std::vector<std::vector<double>> inputPoints(numberOfPoints);
  for(int i = 0; i < numberOfPoints; ++i) {
    for(int j = 0; j < dimensionHigh; ++j)
      inputPoints[i].push_back(inputArrays[j]->GetVariantValue(i).ToDouble());
  }

  /*** representation points ***/
  if(SelectRepresentationFieldsWithRegexp) {
    // select all input columns whose name is matching the regexp
    RepresentationScalarFields.clear();
    const auto n = representation->GetNumberOfColumns();
    for(int i = 0; i < n; ++i) {
      const auto &name = representation->GetColumnName(i);
      if(std::regex_match(name, std::regex(RepresentationRegexpString))) {
        RepresentationScalarFields.emplace_back(name);
      }
    }
  }
  if(representation->GetNumberOfRows() != numberOfPoints
     || RepresentationScalarFields.size() <= 0) {
    this->printErr("Representation matrix has invalid dimensions (rows: "
                   + std::to_string(representation->GetNumberOfRows())
                   + ", columns: "
                   + std::to_string(RepresentationScalarFields.size()) + ")");
    return 0;
  }
  std::vector<vtkAbstractArray *> representationArrays;
  representationArrays.reserve(RepresentationScalarFields.size());
  for(const auto &s : RepresentationScalarFields)
    representationArrays.push_back(representation->GetColumnByName(s.data()));

  const int dimensionLow = RepresentationScalarFields.size();
  std::vector<std::vector<double>> representationPoints(numberOfPoints);
  for(int i = 0; i < numberOfPoints; ++i) {
    for(int j = 0; j < dimensionLow; ++j)
      representationPoints[i].push_back(
        representationArrays[j]->GetVariantValue(i).ToDouble());
  }

  execute(inputPoints, representationPoints);

  auto metricNames = vtkSmartPointer<vtkStringArray>::New();
  metricNames->SetName("Metric");
  metricNames->InsertNextValue("0-Wasserstein");
  metricNames->InsertNextValue("1-Wasserstein");
  metricNames->InsertNextValue("Triplet accuracy");
  metricNames->InsertNextValue("Linear correlation");
  metricNames->InsertNextValue("Trustworthiness");
  metricNames->InsertNextValue("Continuity");
  metricNames->InsertNextValue("LCMC");
  metricNames->InsertNextValue("MRRE (input)");
  metricNames->InsertNextValue("MRRE (latent)");
  metricNames->InsertNextValue("RMSE");
  output->AddColumn(metricNames);

  auto metricValues = vtkSmartPointer<vtkDoubleArray>::New();
  metricValues->SetName("Value");
  metricValues->InsertNextValue(w0_);
  metricValues->InsertNextValue(w1_);
  metricValues->InsertNextValue(ta_);
  metricValues->InsertNextValue(lc_);
  metricValues->InsertNextValue(trust_);
  metricValues->InsertNextValue(cont_);
  metricValues->InsertNextValue(lcmc_);
  metricValues->InsertNextValue(mrreh_);
  metricValues->InsertNextValue(mrrel_);
  metricValues->InsertNextValue(rmse_);
  output->AddColumn(metricValues);

  auto metricRanges = vtkSmartPointer<vtkStringArray>::New();
  metricRanges->SetName("Range");
  metricRanges->InsertNextValue("[0,inf)");
  metricRanges->InsertNextValue("[0,inf)");
  metricRanges->InsertNextValue("[0,1]");
  metricRanges->InsertNextValue("[-1,1]");
  metricRanges->InsertNextValue("[0,1]");
  metricRanges->InsertNextValue("[0,1]");
  metricRanges->InsertNextValue("[0,1]");
  metricRanges->InsertNextValue("[0,1]");
  metricRanges->InsertNextValue("[0,1]");
  metricRanges->InsertNextValue("[0,inf)");
  output->AddColumn(metricRanges);

  // return success
  return 1;
}
