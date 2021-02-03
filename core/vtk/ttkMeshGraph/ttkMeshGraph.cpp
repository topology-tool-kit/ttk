#include <ttkMeshGraph.h>

#include <vtkAbstractArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataSetTriangleFilter.h>

#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include <vtkArrayDispatch.h>
#include <vtkDataArrayAccessor.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#define ttkTypeMacroErrorCase(idx, type)                          \
  default: {                                                      \
    this->printErr("Unsupported " #idx "-th Template Data Type: " \
                   + std::to_string(type));                       \
  } break;

#define ttkTypeMacroCase(enum, type, number, call) \
  case enum: {                                     \
    typedef type T##number;                        \
    call;                                          \
  } break;

#define ttkTypeMacroR(target, call)                                 \
  switch(target) {                                                  \
    ttkTypeMacroCase(VTK_FLOAT, float, 0, call) ttkTypeMacroCase(   \
      VTK_DOUBLE, double, 0, call) ttkTypeMacroErrorCase(0, target) \
  }

#define ttkTypeMacroA(target, call)                                        \
  switch(target) {                                                         \
    ttkTypeMacroCase(VTK_FLOAT, float, 0, call);                           \
    ttkTypeMacroCase(VTK_DOUBLE, double, 0, call);                         \
    ttkTypeMacroCase(VTK_INT, int, 0, call);                               \
    ttkTypeMacroCase(VTK_UNSIGNED_INT, unsigned int, 0, call);             \
    ttkTypeMacroCase(VTK_CHAR, char, 0, call);                             \
    ttkTypeMacroCase(VTK_SIGNED_CHAR, signed char, 0, call);               \
    ttkTypeMacroCase(VTK_UNSIGNED_CHAR, unsigned char, 0, call);           \
    ttkTypeMacroCase(VTK_LONG, long, 0, call);                             \
    ttkTypeMacroCase(VTK_LONG_LONG, long long, 0, call);                   \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG, unsigned long, 0, call);           \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG_LONG, unsigned long long, 0, call); \
    ttkTypeMacroCase(VTK_ID_TYPE, vtkIdType, 0, call);                     \
    ttkTypeMacroErrorCase(0, target);                                      \
  }

#define ttkTypeMacroRR(target0, target1, call)                             \
  switch(target1) {                                                        \
    ttkTypeMacroCase(VTK_FLOAT, float, 1, ttkTypeMacroR(target0, call));   \
    ttkTypeMacroCase(VTK_DOUBLE, double, 1, ttkTypeMacroR(target0, call)); \
    ttkTypeMacroErrorCase(1, target1);                                     \
  }

#define ttkTypeMacroRRR(target0, target1, target2, call)              \
  switch(target2) {                                                   \
    ttkTypeMacroCase(                                                 \
      VTK_FLOAT, float, 2, ttkTypeMacroRR(target0, target1, call));   \
    ttkTypeMacroCase(                                                 \
      VTK_DOUBLE, double, 2, ttkTypeMacroRR(target0, target1, call)); \
    ttkTypeMacroErrorCase(2, target2);                                \
  }

#define ttkTypeMacroRRI(target0, target1, target2, call)                       \
  switch(target2) {                                                            \
    ttkTypeMacroCase(VTK_INT, int, 2, ttkTypeMacroRR(target0, target1, call)); \
    ttkTypeMacroCase(                                                          \
      VTK_LONG_LONG, long long, 2, ttkTypeMacroRR(target0, target1, call));    \
    ttkTypeMacroCase(                                                          \
      VTK_ID_TYPE, vtkIdType, 2, ttkTypeMacroRR(target0, target1, call));      \
    ttkTypeMacroErrorCase(2, target2);                                         \
  }

#define ttkTypeMacroAI(target0, target1, call)                                 \
  switch(target1) {                                                            \
    ttkTypeMacroCase(VTK_INT, int, 1, ttkTypeMacroA(target0, call));           \
    ttkTypeMacroCase(                                                          \
      VTK_LONG_LONG, long long, 1, ttkTypeMacroA(target0, call));              \
    ttkTypeMacroCase(VTK_ID_TYPE, vtkIdType, 1, ttkTypeMacroA(target0, call)); \
    ttkTypeMacroErrorCase(1, target1);                                         \
  }

vtkStandardNewMacro(ttkMeshGraph);

ttkMeshGraph::ttkMeshGraph() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkMeshGraph::~ttkMeshGraph() {
}

int ttkMeshGraph::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkMeshGraph::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkMeshGraph::RequestData(vtkInformation *request,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {
  ttk::Timer t;

  // ---------------------------------------------------------------------------
  // Get Input
  // ---------------------------------------------------------------------------
  auto input = vtkUnstructuredGrid::GetData(inputVector[0]);
  if(input == nullptr) {
    return -1;
  }

  size_t nInputPoints = input->GetNumberOfPoints();
  size_t nInputCells = input->GetNumberOfCells();

  auto inputPointSizes = this->GetInputArrayToProcess(0, inputVector);
  if(!inputPointSizes) {
    this->printErr("Unable to retrieve point size array.");
    return 0;
  }

  // ---------------------------------------------------------------------------
  // Init Output
  // ---------------------------------------------------------------------------

  // Output Points
  auto inputPoints = input->GetPoints();
  auto outputPoints = vtkSmartPointer<vtkPoints>::New();

  outputPoints->SetDataType(inputPoints->GetDataType());
  auto nOutputPoints = this->computeNumberOfOutputPoints(
    nInputPoints, nInputCells, this->GetUseQuadraticCells(),
    this->GetSubdivisions());
  outputPoints->SetNumberOfPoints(nOutputPoints);

  // Output Topology
  auto inputConnectivityArray = input->GetCells()->GetConnectivityArray();

  auto nOutputCells = this->computeNumberOfOutputCells(
    nInputCells, this->GetUseQuadraticCells());
  auto outputTopologySize = this->computeOutputConnectivityArraySize(
    nInputCells, this->GetUseQuadraticCells(), this->GetSubdivisions());

  auto outputConnectivityArray = vtkSmartPointer<vtkDataArray>::Take(
    inputConnectivityArray->NewInstance());
  outputConnectivityArray->SetNumberOfValues(outputTopologySize);

  auto outputOffsetArray = vtkSmartPointer<vtkDataArray>::Take(
    inputConnectivityArray->NewInstance());
  outputOffsetArray->SetNumberOfValues(nOutputCells + 1);

  // ---------------------------------------------------------------------------
  // Compute cells with base code
  // ---------------------------------------------------------------------------

  int status = 0;
  if(this->GetUseQuadraticCells()) {
    ttkTypeMacroRRI(
      inputPoints->GetDataType(), inputPointSizes->GetDataType(),
      outputConnectivityArray->GetDataType(),
      (status = this->execute<T2, T0, T1>(
         // Output
         static_cast<T0 *>(ttkUtils::GetVoidPointer(outputPoints)),
         static_cast<T2 *>(ttkUtils::GetVoidPointer(outputConnectivityArray)),
         static_cast<T2 *>(ttkUtils::GetVoidPointer(outputOffsetArray)),

         // Input
         static_cast<T0 *>(ttkUtils::GetVoidPointer(inputPoints)),
         static_cast<T2 *>(
           ttkUtils::GetVoidPointer(input->GetCells()->GetConnectivityArray())),
         nInputPoints, nInputCells,
         static_cast<T1 *>(ttkUtils::GetVoidPointer(inputPointSizes)),
         this->GetSizeScale(), this->GetSizeAxis())));
  } else {
    ttkTypeMacroRRI(
      inputPoints->GetDataType(), inputPointSizes->GetDataType(),
      outputConnectivityArray->GetDataType(),
      (status = this->execute2<T2, T0, T1>(
         // Output
         static_cast<T0 *>(ttkUtils::GetVoidPointer(outputPoints)),
         static_cast<T2 *>(ttkUtils::GetVoidPointer(outputConnectivityArray)),
         static_cast<T2 *>(ttkUtils::GetVoidPointer(outputOffsetArray)),

         // Input
         static_cast<T0 *>(ttkUtils::GetVoidPointer(inputPoints)),
         static_cast<T2 *>(
           ttkUtils::GetVoidPointer(input->GetCells()->GetConnectivityArray())),
         nInputPoints, nInputCells, this->GetSubdivisions(),
         static_cast<T1 *>(ttkUtils::GetVoidPointer(inputPointSizes)),
         this->GetSizeScale(), this->GetSizeAxis())));
  }
  if(!status)
    return 0;

  // ---------------------------------------------------------------------------
  // Generate meshed graph as vtkUnstructuredGrid
  // ---------------------------------------------------------------------------

  // Create new vtkUnstructuredGrid for meshed graph
  auto meshedGraph = vtkUnstructuredGrid::GetData(outputVector);

  meshedGraph->SetPoints(outputPoints);

  auto outputCellArray = vtkSmartPointer<vtkCellArray>::New();
  outputCellArray->SetData(outputOffsetArray, outputConnectivityArray);
  meshedGraph->SetCells(
    this->GetUseQuadraticCells() ? VTK_QUADRATIC_QUAD : VTK_POLYGON,
    outputCellArray);

  // Copy input point data to output point data
  {
    auto iPointData = input->GetPointData();
    auto oPointData = meshedGraph->GetPointData();

    for(int i = 0; i < iPointData->GetNumberOfArrays(); i++) {
      auto iArray = iPointData->GetArray(i);
      if(iArray->GetNumberOfComponents() > 1)
        continue;

      auto oArray = vtkSmartPointer<vtkDataArray>::Take(
        vtkDataArray::CreateDataArray(iArray->GetDataType()));
      oArray->SetName(iArray->GetName());
      oArray->SetNumberOfTuples(nOutputPoints);
      oArray->SetNumberOfComponents(1);
      oPointData->AddArray(oArray);

      ttkTypeMacroAI(
        iArray->GetDataType(), inputConnectivityArray->GetDataType(),
        (status = this->mapInputPointDataToOutputPointData<T0, T1>(
           static_cast<T0 *>(ttkUtils::GetVoidPointer(oArray)),

           nInputPoints, nInputCells,
           static_cast<T1 *>(ttkUtils::GetVoidPointer(inputConnectivityArray)),
           static_cast<T0 *>(ttkUtils::GetVoidPointer(iArray)),
           this->GetUseQuadraticCells(), this->GetSubdivisions())));

      if(!status)
        return 0;
    }
  }

  // Copy input cell data to output cell data
  {
    auto iCellData = input->GetCellData();
    auto oCellData = meshedGraph->GetCellData();

    for(int i = 0; i < iCellData->GetNumberOfArrays(); i++) {
      auto iArray = iCellData->GetArray(i);
      if(iArray->GetNumberOfComponents() > 1)
        continue;

      auto oArray = vtkSmartPointer<vtkDataArray>::Take(
        vtkDataArray::CreateDataArray(iArray->GetDataType()));
      oArray->SetName(iArray->GetName());
      oArray->SetNumberOfTuples(nOutputCells);
      oArray->SetNumberOfComponents(1);
      oCellData->AddArray(oArray);

      ttkTypeMacroA(
        iArray->GetDataType(),
        (status = this->mapInputCellDataToOutputCellData<T0>(
           static_cast<T0 *>(ttkUtils::GetVoidPointer(oArray)),

           nInputCells, static_cast<T0 *>(ttkUtils::GetVoidPointer(iArray)),
           this->GetUseQuadraticCells(), this->GetSubdivisions())));
      if(!status)
        return 0;
    }
  }

  // -------------------------------------------------------------------------
  // Print status
  // -------------------------------------------------------------------------
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete", 1, t.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
