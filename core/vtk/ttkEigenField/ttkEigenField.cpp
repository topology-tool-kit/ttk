#include <ttkEigenField.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

// VTK includes
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

vtkStandardNewMacro(ttkEigenField);

ttkEigenField::ttkEigenField() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkEigenField::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkEigenField::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkEigenField::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  const auto domain = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);
  auto triangulation = ttkAlgorithm::GetTriangulation(domain);

  if(triangulation == nullptr) {
    this->printErr("Triangulation is NULL");
    return 0;
  }

  this->preconditionTriangulation(*triangulation);

  int res = 0;

  // array of eigenfunctions
  vtkSmartPointer<vtkDataArray> eigenFunctions{};
  // statistics
  vtkSmartPointer<vtkDataArray> stats{};

  switch(OutputFieldType) {
    case FieldType::FLOAT:
      eigenFunctions = vtkSmartPointer<vtkFloatArray>::New();
      stats = vtkSmartPointer<vtkFloatArray>::New();
      break;
    case FieldType::DOUBLE:
      eigenFunctions = vtkSmartPointer<vtkDoubleArray>::New();
      stats = vtkSmartPointer<vtkDoubleArray>::New();
      break;
    default:
      this->printErr("Unknown field type");
      return 0;
  }

  if(eigenFunctions == nullptr) {
    this->printErr("vtkDataArray allocation problem");
    return 0;
  }

  const auto vertexNumber = triangulation->getNumberOfVertices();

  eigenFunctions->SetNumberOfComponents(EigenNumber);
  eigenFunctions->SetNumberOfTuples(vertexNumber);
  eigenFunctions->SetName(OutputFieldName.data());

  if(ComputeStatistics) {
    stats->SetName("Statistics");
    const int statsComp = 4;
    stats->SetNumberOfComponents(statsComp);
    stats->SetNumberOfTuples(vertexNumber);
    stats->SetComponentName(0, "Min");
    stats->SetComponentName(1, "Max");
    stats->SetComponentName(2, "Sum");
    stats->SetComponentName(3, "EigenMagnitude");
  }

  switch(OutputFieldType) {
    case FieldType::FLOAT:
      res += this->execute<float>(
        *triangulation,
        static_cast<float *>(ttkUtils::GetVoidPointer(eigenFunctions)),
        EigenNumber, ComputeStatistics,
        static_cast<float *>(ttkUtils::GetVoidPointer(stats)));
      break;
    case FieldType::DOUBLE:
      res += this->execute<double>(
        *triangulation,
        static_cast<double *>(ttkUtils::GetVoidPointer(eigenFunctions)),
        EigenNumber, ComputeStatistics,
        static_cast<double *>(ttkUtils::GetVoidPointer(stats)));
      break;
    default:
      break;
  }

  if(res != 0) {
    this->printErr("EigenField execute error code " + std::to_string(res));
    return 0;
  }

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(eigenFunctions);

  if(ComputeStatistics) {
    output->GetPointData()->AddArray(stats);
  }

  return 1;
}
