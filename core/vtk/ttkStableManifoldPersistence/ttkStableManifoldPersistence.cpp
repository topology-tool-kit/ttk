#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkStableManifoldPersistence.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkStableManifoldPersistence);

ttkStableManifoldPersistence::ttkStableManifoldPersistence() {
}

int ttkStableManifoldPersistence::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkStableManifoldPersistence::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkStableManifoldPersistence::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  using ttk::SimplexId;

  const auto domain = vtkDataSet::GetData(inputVector[0]);
  const auto constraints = vtkPointSet::GetData(inputVector[1]);
  if(!domain || !constraints)
    return !this->printErr("Unable to retrieve required input data objects.");

  auto output = vtkDataSet::GetData(outputVector);

  // domain scalar field
  const auto inputScalars = this->GetInputArrayToProcess(0, domain);
  if(!inputScalars) {
    this->printErr("Input scalar field pointer is null.");
    return -3;
  }

  // create output arrays
  auto outputScalars
    = vtkSmartPointer<vtkDataArray>::Take(inputScalars->NewInstance());
  outputScalars->DeepCopy(inputScalars);

  // constraint identifier field
  int numberOfConstraints = constraints->GetNumberOfPoints();

  int ret = 0;

  // something wrong in baseCode
  if(ret) {
    this->printErr("StableManifoldPersistence.execute() error code: "
                   + std::to_string(ret));
    return -12;
  }

  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(outputScalars);

  return 1;
}
