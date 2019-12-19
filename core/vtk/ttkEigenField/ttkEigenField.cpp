#include <ttkEigenField.h>

// VTK includes
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#define MODULE_S "[ttkEigenField] "
#define MODULE_ERROR_S MODULE_S "Error: "
#ifndef TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)    \
  if(COND) {                            \
    cerr << MODULE_ERROR_S MSG << endl; \
    return RET;                         \
  }
#else // TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)
#endif // TTK_ENABLE_KAMIKAZE

vtkStandardNewMacro(ttkEigenField);

ttkEigenField::ttkEigenField() {
  SetNumberOfInputPorts(1);
}

int ttkEigenField::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  }
  return 1;
}

int ttkEigenField::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  }
  return 1;
}

int ttkEigenField::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  const auto domain = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);
  auto triangulation = ttkAlgorithm::GetTriangulation(domain);

  this->preconditionTriangulation(triangulation);

  int res = 0;

  // array of eigenfunctions
  vtkSmartPointer<vtkDataArray> eigenFunctions{};
  // statistics
  vtkSmartPointer<vtkDataArray> stats{};

  switch(OutputFieldType) {
    case EigenFieldType::Float:
      eigenFunctions = vtkSmartPointer<vtkFloatArray>::New();
      stats = vtkSmartPointer<vtkFloatArray>::New();
      break;
    case EigenFieldType::Double:
      eigenFunctions = vtkSmartPointer<vtkDoubleArray>::New();
      stats = vtkSmartPointer<vtkDoubleArray>::New();
      break;
    default:
      TTK_ABORT_KK(true, "unknown field type", -7);
      break;
  }

  TTK_ABORT_KK(eigenFunctions == nullptr, "vtkArray allocation problem", -3);

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
    case EigenFieldType::Float:
      res += this->execute<float>(
        triangulation, static_cast<float *>(eigenFunctions->GetVoidPointer(0)),
        EigenNumber, ComputeStatistics,
        static_cast<float *>(stats->GetVoidPointer(0)));
      break;
    case EigenFieldType::Double:
      res += this->execute<double>(
        triangulation, static_cast<double *>(eigenFunctions->GetVoidPointer(0)),
        EigenNumber, ComputeStatistics,
        static_cast<double *>(stats->GetVoidPointer(0)));
      break;
    default:
      break;
  }

  TTK_ABORT_KK(res != 0, "EigenField execute error code " << res, -4);

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(eigenFunctions);

  if(ComputeStatistics) {
    output->GetPointData()->AddArray(stats);
  }

  return 1;
}
