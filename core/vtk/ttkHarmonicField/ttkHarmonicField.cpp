#include <ttkHarmonicField.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

// VTK includes
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>

vtkStandardNewMacro(ttkHarmonicField);

ttkHarmonicField::ttkHarmonicField() {
  SetNumberOfInputPorts(2);
  SetNumberOfOutputPorts(1);
}

int ttkHarmonicField::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkHarmonicField::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkHarmonicField::RequestData(vtkInformation *ttkNotUsed(request),
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {

  const auto domain = vtkDataSet::GetData(inputVector[0]);
  const auto identifiers = vtkPointSet::GetData(inputVector[1]);
  auto output = vtkDataSet::GetData(outputVector);
  auto triangulation = ttkAlgorithm::GetTriangulation(domain);

  if(triangulation == nullptr) {
    this->printErr("No triangulation");
    return 0;
  }
  this->preconditionTriangulation(*triangulation, UseCotanWeights);

  vtkDataArray *inputField = this->GetInputArrayToProcess(0, identifiers);
  std::vector<ttk::SimplexId> idSpareStorage{};
  const auto *vertsid = this->GetIdentifierArrayPtr(
    ForceConstraintIdentifiers, 1, ttk::VertexScalarFieldName, identifiers,
    idSpareStorage);

  if(vertsid == nullptr || inputField == nullptr) {
    this->printErr("Input fields are NULL");
    return 0;
  }

  if(inputField->IsA("vtkDoubleArray")) {
    OutputScalarFieldType = FieldType::DOUBLE;
  } else if(inputField->IsA("vtkFloatArray")) {
    OutputScalarFieldType = FieldType::FLOAT;
  } else {
    this->printErr("Filter only supports floating point scalar fields");
    this->printErr(
      "Please select a floating point input scalar field or convert an "
      "existing one with TTKPointDataConverter or TTKArrayEditor");
    return -2;
  }

  const auto nVerts = domain->GetNumberOfPoints();
  const auto nSources = identifiers->GetNumberOfPoints();

  vtkSmartPointer<vtkDataArray> outputField{};

  switch(OutputScalarFieldType) {
    case FieldType::FLOAT:
      outputField = vtkSmartPointer<vtkFloatArray>::New();
      break;
    case FieldType::DOUBLE:
      outputField = vtkSmartPointer<vtkDoubleArray>::New();
      break;
    default:
      this->printErr("Unknown scalar field type");
      return -7;
      break;
  }

  if(outputField == nullptr) {
    this->printErr("vtkArray allocation problem");
    return 0;
  }

  outputField->SetNumberOfComponents(1);
  outputField->SetNumberOfTuples(nVerts);
  outputField->SetName(OutputScalarFieldName.data());
  int res{};

  switch(OutputScalarFieldType) {
    case FieldType::FLOAT:
      res = this->execute<float>(
        *triangulation, nSources, vertsid,
        static_cast<float *>(ttkUtils::GetVoidPointer(inputField)),
        static_cast<float *>(ttkUtils::GetVoidPointer(outputField)),
        UseCotanWeights, SolvingMethod, LogAlpha);
      break;
    case FieldType::DOUBLE:
      res = this->execute<double>(
        *triangulation, nSources, vertsid,
        static_cast<double *>(ttkUtils::GetVoidPointer(inputField)),
        static_cast<double *>(ttkUtils::GetVoidPointer(outputField)),
        UseCotanWeights, SolvingMethod, LogAlpha);
      break;
    default:
      break;
  }

  if(res != 0) {
    this->printErr("HarmonicField execute() error code: "
                   + std::to_string(res));
    return 0;
  }

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(outputField);

  return 1;
}
