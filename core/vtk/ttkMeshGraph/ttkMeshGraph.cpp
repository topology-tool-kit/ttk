#include <ttkMeshGraph.h>

#include <vtkAbstractArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

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

  size_t nInputPoints = input->GetNumberOfPoints();
  size_t nInputCells = input->GetNumberOfCells();
  auto inputCells = input->GetCells();

  auto inputPointSizes = this->GetInputArrayToProcess(0, inputVector);
  if(this->GetUseVariableSize() && !inputPointSizes) {
    this->printErr("Unable to retrieve point size array.");
    return 0;
  }
  int sizeType
    = this->GetUseVariableSize() ? inputPointSizes->GetDataType() : VTK_CHAR;

  // ---------------------------------------------------------------------------
  // Init Output
  // ---------------------------------------------------------------------------

  // Output Points
  auto nOutputPoints = this->computeNumberOfOutputPoints(
    nInputPoints, nInputCells, this->GetUseQuadraticCells(),
    this->GetSubdivisions());
  auto outputPoints = vtkSmartPointer<vtkPoints>::New();
  outputPoints->SetNumberOfPoints(nOutputPoints);
  auto outputVertices = (float *)outputPoints->GetVoidPointer(0);

  // Output Topology
  auto nOutputCells = this->computeNumberOfOutputCells(
    nInputCells, this->GetUseQuadraticCells());
  auto outputTopologySize = this->computeOutputConnectivityListSize(
    nInputCells, this->GetUseQuadraticCells(), this->GetSubdivisions());
  auto outputCells = vtkSmartPointer<vtkIdTypeArray>::New();
  outputCells->SetNumberOfValues(outputTopologySize);

  // ---------------------------------------------------------------------------
  // Compute cells with base code
  // ---------------------------------------------------------------------------
  int status = 0;
  if(this->GetUseQuadraticCells()) {
    // Quadratic cells
    switch(sizeType) {
      vtkTemplateMacro(
        (status = this->execute<vtkIdType, VTK_TT>(
           // Output
           outputVertices, (vtkIdType *)outputCells->GetVoidPointer(0),

           // Input
           (float *)input->GetPoints()->GetVoidPointer(0),
           inputCells->GetData()->GetPointer(0), nInputPoints, nInputCells,
           this->GetUseVariableSize()
             ? (VTK_TT *)inputPointSizes->GetVoidPointer(0)
             : nullptr,
           this->GetSizeScale(), this->GetSizeAxis())));
    }
  } else {
    // Linear Polygons
    switch(sizeType) {
      vtkTemplateMacro(
        (status = this->execute2<vtkIdType, VTK_TT>(
           // Output
           outputVertices, (vtkIdType *)outputCells->GetVoidPointer(0),

           // Input
           (float *)input->GetPoints()->GetVoidPointer(0),
           inputCells->GetData()->GetPointer(0), nInputPoints, nInputCells,
           this->GetSubdivisions(),
           this->GetUseVariableSize()
             ? (VTK_TT *)inputPointSizes->GetVoidPointer(0)
             : nullptr,
           this->GetSizeScale(), this->GetSizeAxis())));
    }
  }
  if(status != 1)
    return 0;

  // ---------------------------------------------------------------------------
  // Generate meshed graph as vtkUnstructuredGrid
  // ---------------------------------------------------------------------------

  // Create new vtkUnstructuredGrid for meshed graph
  auto meshedGraph = vtkSmartPointer<vtkUnstructuredGrid>::New();
  meshedGraph->SetPoints(outputPoints);
  auto outputCellArray = vtkSmartPointer<vtkCellArray>::New();
  outputCellArray->SetCells(nOutputCells, outputCells);
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

      auto oArray = vtkDataArray::CreateDataArray(iArray->GetDataType());
      oArray->SetName(iArray->GetName());
      oArray->SetNumberOfTuples(nOutputPoints);
      oArray->SetNumberOfComponents(1);
      oPointData->AddArray(oArray);

      switch(iArray->GetDataType()) {
        vtkTemplateMacro(
          (status = this->mapInputPointDataToOutputPointData<vtkIdType, VTK_TT>(
             inputCells->GetData()->GetPointer(0), nInputPoints, nInputCells,

             (VTK_TT *)iArray->GetVoidPointer(0),
             (VTK_TT *)oArray->GetVoidPointer(0),

             this->GetUseQuadraticCells(), this->GetSubdivisions())));
      }
      if(status != 1)
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

      switch(iArray->GetDataType()) {
        vtkTemplateMacro(
          (status = this->mapInputCellDataToOutputCellData<vtkIdType, VTK_TT>(
             nInputCells,

             (VTK_TT *)iArray->GetVoidPointer(0),
             (VTK_TT *)oArray->GetVoidPointer(0),

             this->GetUseQuadraticCells(), this->GetSubdivisions())));
      }
      if(status != 1)
        return 0;
    }
  }

  // ---------------------------------------------------------------------------
  // Finalize Output
  // ---------------------------------------------------------------------------
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  if(!this->GetUseQuadraticCells() && this->GetSubdivisions() > 1
     && this->GetTetrahedralize()) {
    auto dataSetTriangleFilter
      = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
    dataSetTriangleFilter->SetInputData(meshedGraph);
    dataSetTriangleFilter->Update();

    output->ShallowCopy(dataSetTriangleFilter->GetOutput());
  } else {
    output->ShallowCopy(meshedGraph);
  }

  // -------------------------------------------------------------------------
  // Print status
  // -------------------------------------------------------------------------
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete", 1, t.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
