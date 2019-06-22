#include <ttkMeshGraph.h>

#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkMeshGraph)

  int ttkMeshGraph::RequestData(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {
  Timer t;
  Memory m;

  // Print status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl;
    msg << "[ttkMeshGraph] RequestData" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  // Set Wrapper
  meshGraph.setWrapper(this);

  // -------------------------------------------------------------------------
  // Get Input
  // -------------------------------------------------------------------------
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  auto input = vtkUnstructuredGrid::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  size_t nInputPoints = input->GetNumberOfPoints();
  size_t nInputCells = input->GetNumberOfCells();
  auto inputCells = input->GetCells();

  auto inputPointSizes
    = input->GetPointData()->GetArray(this->GetSizeFieldName().data());
  if(this->GetUseVariableSize() && !inputPointSizes) {
    dMsg(cout,
         "[ttkMeshGraph] ERROR: Input point data does not have array '"
           + this->GetSizeFieldName() + "'.\n",
         fatalMsg);
    return 0;
  }
  int sizeType
    = this->GetUseVariableSize() ? inputPointSizes->GetDataType() : VTK_CHAR;

  // -------------------------------------------------------------------------
  // Init Output
  // -------------------------------------------------------------------------

  // Output Points
  auto nOutputPoints = meshGraph.computeNumberOfOutputPoints(
    nInputPoints, nInputCells, this->GetUseQuadraticCells(),
    this->GetSubdivisions());
  auto outputPoints = vtkSmartPointer<vtkPoints>::New();
  outputPoints->SetNumberOfPoints(nOutputPoints);
  auto outputVertices = (float *)outputPoints->GetVoidPointer(0);

  // Output Topology
  auto nOutputCells = meshGraph.computeNumberOfOutputCells(
    nInputCells, this->GetUseQuadraticCells());
  auto outputTopologySize = meshGraph.computeOutputTopologySize(
    nInputCells, this->GetUseQuadraticCells(), this->GetSubdivisions());
  auto outputCells = vtkSmartPointer<vtkIdTypeArray>::New();
  outputCells->SetNumberOfValues(outputTopologySize);

  // -------------------------------------------------------------------------
  // Compute cells with base code
  // -------------------------------------------------------------------------
  int status = 0;
  if(this->GetUseQuadraticCells()) {
    // Quadratic cells
    switch(sizeType) {
      ttkTemplateMacro({
        status = meshGraph.execute<vtkIdType TTK_COMMA VTK_TT>(
          // Input
          (float *)input->GetPoints()->GetVoidPointer(0),
          inputCells->GetPointer(), nInputPoints, nInputCells,

          this->GetUseVariableSize()
            ? (VTK_TT *)inputPointSizes->GetVoidPointer(0)
            : nullptr,
          this->GetSizeScale(), this->GetSizeAxis(),

          // Output
          outputVertices, (vtkIdType *)outputCells->GetVoidPointer(0));
      });
    }
  } else {
    // Linear Polygons
    switch(sizeType) {
      ttkTemplateMacro({
        status = meshGraph.execute2<vtkIdType TTK_COMMA VTK_TT>(
          // Input
          (float *)input->GetPoints()->GetVoidPointer(0),
          inputCells->GetPointer(), nInputPoints, nInputCells,
          this->GetSubdivisions(),

          this->GetUseVariableSize()
            ? (VTK_TT *)inputPointSizes->GetVoidPointer(0)
            : nullptr,
          this->GetSizeScale(), this->GetSizeAxis(),

          // Output
          outputVertices, (vtkIdType *)outputCells->GetVoidPointer(0));
      });
    }
  }
  if(status != 1)
    return 0;

  // -------------------------------------------------------------------------
  // Generate meshed graph as vtkUnstructuredGrid
  // -------------------------------------------------------------------------

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
        ttkTemplateMacro({
          status
            = meshGraph
                .mapInputPointDataToOutputPointData<vtkIdType TTK_COMMA VTK_TT>(
                  inputCells->GetPointer(), nInputPoints, nInputCells,

                  (VTK_TT *)iArray->GetVoidPointer(0),
                  (VTK_TT *)oArray->GetVoidPointer(0),

                  this->GetUseQuadraticCells(), this->GetSubdivisions());
          if(status != 1)
            return 0;
        });
      }
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
        ttkTemplateMacro({
          status
            = meshGraph
                .mapInputCellDataToOutputCellData<vtkIdType TTK_COMMA VTK_TT>(
                  nInputCells,

                  (VTK_TT *)iArray->GetVoidPointer(0),
                  (VTK_TT *)oArray->GetVoidPointer(0),

                  this->GetUseQuadraticCells(), this->GetSubdivisions());
          if(status != 1)
            return 0;
        });
      }
    }
  }

  // -------------------------------------------------------------------------
  // Finalize Output
  // -------------------------------------------------------------------------
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
  {
    stringstream msg;
    msg << "[ttkMeshGraph] "
           "-----------------------------------------------------------------"
        << endl
        << "[ttkMeshGraph]   Time: " << t.getElapsedTime() << " s" << endl
        << "[ttkMeshGraph] Memory: " << m.getElapsedUsage() << " MB" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}
