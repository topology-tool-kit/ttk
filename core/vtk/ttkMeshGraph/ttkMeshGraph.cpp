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

int ttkMeshGraph::RequestData(vtkInformation *ttkNotUsed(request),
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
         ttkUtils::GetPointer<T0>(outputPoints->GetData()),
         ttkUtils::GetPointer<T2>(outputConnectivityArray),
         ttkUtils::GetPointer<T2>(outputOffsetArray),

         // Input
         ttkUtils::GetPointer<T0>(inputPoints->GetData()),
         ttkUtils::GetPointer<T2>(input->GetCells()->GetConnectivityArray()),
         nInputPoints, nInputCells, ttkUtils::GetPointer<T1>(inputPointSizes),
         this->GetSizeScale(), this->GetSizeAxis())));
  } else {
    ttkTypeMacroRRI(
      inputPoints->GetDataType(), inputPointSizes->GetDataType(),
      outputConnectivityArray->GetDataType(),
      (status = this->execute2<T2, T0, T1>(
         // Output
         ttkUtils::GetPointer<T0>(outputPoints->GetData()),
         ttkUtils::GetPointer<T2>(outputConnectivityArray),
         ttkUtils::GetPointer<T2>(outputOffsetArray),

         // Input
         ttkUtils::GetPointer<T0>(inputPoints->GetData()),
         ttkUtils::GetPointer<T2>(input->GetCells()->GetConnectivityArray()),
         nInputPoints, nInputCells, this->GetSubdivisions(),
         ttkUtils::GetPointer<T1>(inputPointSizes), this->GetSizeScale(),
         this->GetSizeAxis())));
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
           ttkUtils::GetPointer<T0>(oArray),

           nInputPoints, nInputCells,
           ttkUtils::GetPointer<T1>(inputConnectivityArray),
           ttkUtils::GetPointer<T0>(iArray), this->GetUseQuadraticCells(),
           this->GetSubdivisions())));

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
           ttkUtils::GetPointer<T0>(oArray), nInputCells,
           ttkUtils::GetPointer<T0>(iArray), this->GetUseQuadraticCells())));
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
