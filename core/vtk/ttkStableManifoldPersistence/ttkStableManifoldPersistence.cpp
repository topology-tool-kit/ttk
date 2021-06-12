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

int ttkStableManifoldPersistence::BuildSimplex2PersistenceMap(
  vtkPolyData *criticalPoints,
  vtkUnstructuredGrid *persistenceDiagram,
  std::vector<double> &simplex2persistence) const {

  vtkDataArray *vertexIds
    = persistenceDiagram->GetPointData()->GetArray(ttk::VertexScalarFieldName);

  if(!vertexIds) {
    printErr("The input #2 is not a valid persistence diagram.");
    return -1;
  }

  int maximumVertexId = vertexIds->GetMaxNorm();

  vtkDataArray *criticalPointVertexIds
    = criticalPoints->GetPointData()->GetArray(ttk::VertexScalarFieldName);
  if(!criticalPointVertexIds) {
    printErr("The input #1 is not a valid set of critical points.");
    return -2;
  }
  int maximumCriticalVertexId = criticalPointVertexIds->GetMaxNorm();

  if(maximumCriticalVertexId > maximumVertexId)
    maximumVertexId = maximumCriticalVertexId;

  std::vector<double> vertex2persistence(maximumVertexId, -1);

  return 0;
}

int ttkStableManifoldPersistence::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  ttk::Timer t;

  const auto stableManifold = vtkDataSet::GetData(inputVector[0]);
  const auto criticalPoints = vtkPolyData::GetData(inputVector[1]);
  const auto persistenceDiagram = vtkUnstructuredGrid::GetData(inputVector[2]);

  std::vector<double> simplex2persistence;
  int ret = this->BuildSimplex2PersistenceMap(
    criticalPoints, persistenceDiagram, simplex2persistence);

  if(ret)
    return ret;
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

  // TODO: parallelism
  // TODO: documentation

  return 1;
}
