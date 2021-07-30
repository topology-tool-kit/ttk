#include <ttkRangePolygon.h>

#include <vtkCleanPolyData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkRangePolygon);

ttkRangePolygon::ttkRangePolygon() {

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  setDebugMsgPrefix("RangePolygon");

  vtkWarningMacro("`TTK RangePolygon' is now deprecated. Please use instead "
                  "`Poly Line Source' followed by `Resample With Dataset'.");
}

ttkRangePolygon::~ttkRangePolygon() {
}

int ttkRangePolygon::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  else
    return 0;

  return 1;
}

int ttkRangePolygon::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  else
    return 0;

  return 1;
}

int ttkRangePolygon::RequestData(vtkInformation *ttkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {

  vtkUnstructuredGrid *input = vtkUnstructuredGrid::GetData(inputVector[0]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector);

  if(input->GetNumberOfCells()) {

    if(input->GetCell(0)->GetCellDimension() == 0) {
      processPoints(input, output);
    } else {
      processTriangles(input, output);
    }
  } else {
    processPoints(input, output);
  }

  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}

int ttkRangePolygon::processPoints(vtkUnstructuredGrid *input,
                                   vtkUnstructuredGrid *output) {

  Timer t;

  vtkSmartPointer<vtkPoints> pointSet = vtkSmartPointer<vtkPoints>::New();
  output->SetPoints(pointSet);

  output->GetPoints()->ShallowCopy(input->GetPoints());
  output->GetPointData()->ShallowCopy(input->GetPointData());

  vtkSmartPointer<vtkCellArray> edgeArray
    = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
  idList->SetNumberOfIds(2);

  for(SimplexId i = 0; i < input->GetNumberOfPoints(); i++) {
    if(i) {
      idList->SetId(0, i - 1);
      idList->SetId(1, i);

      edgeArray->InsertNextCell(idList);
    }
  }
  if(ClosedLoop) {
    idList->SetId(0, input->GetNumberOfPoints() - 1);
    idList->SetId(1, 0);
    edgeArray->InsertNextCell(idList);
  }

  output->SetCells(VTK_LINE, edgeArray);

  printMsg(std::to_string(output->GetNumberOfCells()) + " edges extracted", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}

int ttkRangePolygon::processTriangles(vtkUnstructuredGrid *input,
                                      vtkUnstructuredGrid *output) {

  Timer t;

  vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceMaker
    = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();

  surfaceMaker->SetInputData(input);

  vtkSmartPointer<vtkCleanPolyData> surfaceCleaner
    = vtkSmartPointer<vtkCleanPolyData>::New();

  surfaceCleaner->SetInputConnection(surfaceMaker->GetOutputPort());

  vtkSmartPointer<vtkFeatureEdges> featureEdges
    = vtkSmartPointer<vtkFeatureEdges>::New();

  featureEdges->SetBoundaryEdges(true);
  featureEdges->SetManifoldEdges(false);
  featureEdges->SetFeatureEdges(false);
  featureEdges->SetNonManifoldEdges(false);
  featureEdges->SetColoring(false);
  featureEdges->SetInputConnection(surfaceCleaner->GetOutputPort());

  vtkSmartPointer<vtkDataSetTriangleFilter> triangleMaker
    = vtkSmartPointer<vtkDataSetTriangleFilter>::New();

  triangleMaker->SetInputConnection(featureEdges->GetOutputPort());
  triangleMaker->Update();

  output->ShallowCopy(triangleMaker->GetOutput());

  if(NumberOfIterations > 0) {

    // set up the triangulation
    Triangulation *triangulation = ttkAlgorithm::GetTriangulation(output);

    ScalarFieldSmoother smoother;
    smoother.setDimensionNumber(3);
    smoother.setInputDataPointer(output->GetPoints()->GetVoidPointer(0));
    smoother.setOutputDataPointer(output->GetPoints()->GetVoidPointer(0));
    smoother.preconditionTriangulation(triangulation);

    switch(output->GetPoints()->GetDataType()) {
      vtkTemplateMacro(
        smoother.smooth<VTK_TT>(triangulation, NumberOfIterations));
    }

    for(int i = 0; i < output->GetPointData()->GetNumberOfArrays(); i++) {
      vtkDataArray *field = output->GetPointData()->GetArray(i);

      smoother.setDebugLevel(0);

      smoother.setDimensionNumber(field->GetNumberOfComponents());
      smoother.setInputDataPointer(field->GetVoidPointer(0));

      smoother.setOutputDataPointer(field->GetVoidPointer(0));
      smoother.preconditionTriangulation(triangulation);

      switch(field->GetDataType()) {
        vtkTemplateMacro(
          smoother.smooth<VTK_TT>(triangulation, NumberOfIterations));
      }
    }
  }

  printMsg(std::to_string(output->GetNumberOfCells()) + " edges extracted", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}
