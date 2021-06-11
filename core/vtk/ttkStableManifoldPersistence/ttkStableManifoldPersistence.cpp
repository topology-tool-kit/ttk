#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <ttkMacros.h>
#include <ttkStableManifoldPersistence.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkStableManifoldPersistence);

ttkStableManifoldPersistence::ttkStableManifoldPersistence() {

  this->setDebugMsgPrefix("StableManifoldPersistence");
  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(1);
}

int ttkStableManifoldPersistence::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }else if(port == 2){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
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

  ttk::Timer t;

  this->printMsg("Go for it!");

  const auto stableManifold = vtkDataSet::GetData(inputVector[0]);
  const auto criticalPoints = vtkPolyData::GetData(inputVector[1]);
  const auto persistenceDiagram = vtkUnstructuredGrid::GetData(inputVector[2]);

//
//
  auto output = vtkDataSet::GetData(outputVector);
//
//   // preparing the outpout data structure
  output->ShallowCopy(stableManifold);
//
// //   vtkNew<vtkDoubleArray> outputPersistence;
//
// //   output->GetCellData()->AddArray(outputPersistence);

  return 1;
}
