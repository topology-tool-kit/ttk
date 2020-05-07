#include <ttkRangePolygon.h>

#include <vtkPointData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkFeatureEdges.h>
#include <vtkDataSetTriangleFilter.h>

vtkStandardNewMacro(ttkRangePolygon);

ttkRangePolygon::ttkRangePolygon() {
  this->setDebugMsgPrefix("RangePolygon");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkRangePolygon::~ttkRangePolygon() {
}

int ttkRangePolygon::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkRangePolygon::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;

}
int ttkRangePolygon::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {
  auto input = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  if(input->GetNumberOfCells()) {
    if(input->GetCell(0)->GetCellDimension() == 0) {
      processPoints(input, output);
    } else {
      processTriangles(input, output);
    }
  } else {
    processPoints(input, output);
  }

  return 1;

}
int ttkRangePolygon::processPoints(vtkUnstructuredGrid *input,
                                   vtkUnstructuredGrid *output) {

  ttk::Timer t;

  vtkSmartPointer<vtkPoints> pointSet = vtkSmartPointer<vtkPoints>::New();
  output->SetPoints(pointSet);

  output->GetPoints()->ShallowCopy(input->GetPoints());
  output->GetPointData()->ShallowCopy(input->GetPointData());

  vtkSmartPointer<vtkCellArray> edgeArray
    = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
  idList->SetNumberOfIds(2);

  for(ttk::SimplexId i = 0; i < input->GetNumberOfPoints(); i++) {
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

  {
    std::stringstream msg;
    msg << "[ttkRangePolygon] Range polygon extracted in " << t.getElapsedTime()
        << " s. (" << output->GetNumberOfCells() << " edge(s))" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

int ttkRangePolygon::processTriangles(vtkUnstructuredGrid *input,
                                      vtkUnstructuredGrid *output) {

  ttk::Memory m;
  ttk::Timer t;

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
    auto triangulation = ttkAlgorithm::GetTriangulation(output);

    this->setDimensionNumber(3);
    this->setInputDataPointer(output->GetPoints()->GetVoidPointer(0));
    this->setOutputDataPointer(output->GetPoints()->GetVoidPointer(0));
    this->setupTriangulation(triangulation);

    switch(output->GetPoints()->GetDataType()) {
      vtkTemplateMacro(
        this->smooth<VTK_TT>(triangulation, NumberOfIterations));
    }

    for(int i = 0; i < output->GetPointData()->GetNumberOfArrays(); i++) {
      vtkDataArray *field = output->GetPointData()->GetArray(i);

      this->setDimensionNumber(field->GetNumberOfComponents());
      this->setInputDataPointer(field->GetVoidPointer(0));

      this->setOutputDataPointer(field->GetVoidPointer(0));
      this->setupTriangulation(triangulation);

      switch(field->GetDataType()) {
        vtkTemplateMacro(
          this->smooth<VTK_TT>(triangulation, NumberOfIterations));
      }
    }
  }

  {
    std::stringstream msg;
    msg << "[ttkRangePolygon] Range polygon extracted in " << t.getElapsedTime()
        << " s. (" << output->GetNumberOfCells() << " edge(s))" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  {
    std::stringstream msg;
    msg << "[ttkRangePolygon] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
