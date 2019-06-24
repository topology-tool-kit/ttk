#include <ttkRangePolygon.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkRangePolygon)

  ttkRangePolygon::ttkRangePolygon() {

  ClosedLoop = false;
  NumberOfIterations = 0;
  UseAllCores = true;
}

ttkRangePolygon::~ttkRangePolygon() {
}
//

int ttkRangePolygon::doIt(vector<vtkDataSet *> &inputs,
                          vector<vtkDataSet *> &outputs) {

  vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  if(input->GetNumberOfCells()) {

    if(input->GetCell(0)->GetCellDimension() == 0) {
      processPoints(input, output);
    } else {
      processTriangles(input, output);
    }
  } else {
    processPoints(input, output);
  }

  return 0;
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

  {
    stringstream msg;
    msg << "[ttkRangePolygon] Range polygon extracted in " << t.getElapsedTime()
        << " s. (" << output->GetNumberOfCells() << " edge(s))" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

int ttkRangePolygon::processTriangles(vtkUnstructuredGrid *input,
                                      vtkUnstructuredGrid *output) {

  Memory m;
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
    Triangulation *triangulation = ttkTriangulation::getTriangulation(output);
    triangulation->setWrapper(this);

    ScalarFieldSmoother smoother;
    smoother.setWrapper(this);
    smoother.setDimensionNumber(3);
    smoother.setInputDataPointer(output->GetPoints()->GetVoidPointer(0));
    smoother.setOutputDataPointer(output->GetPoints()->GetVoidPointer(0));
    smoother.setupTriangulation(triangulation);

    switch(output->GetPoints()->GetDataType()) {
      vtkTemplateMacro({ smoother.smooth<VTK_TT>(NumberOfIterations); });
    }

    for(int i = 0; i < output->GetPointData()->GetNumberOfArrays(); i++) {
      vtkDataArray *field = output->GetPointData()->GetArray(i);

      switch(field->GetDataType()) {

        vtkTemplateMacro({
          smoother.setWrapper(this);
          smoother.setDebugLevel(0);

          smoother.setDimensionNumber(field->GetNumberOfComponents());
          smoother.setInputDataPointer(field->GetVoidPointer(0));

          smoother.setOutputDataPointer(field->GetVoidPointer(0));
          smoother.setupTriangulation(triangulation);
          smoother.smooth<VTK_TT>(NumberOfIterations);
        });
      }
    }
  }

  {
    stringstream msg;
    msg << "[ttkRangePolygon] Range polygon extracted in " << t.getElapsedTime()
        << " s. (" << output->GetNumberOfCells() << " edge(s))" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  {
    stringstream msg;
    msg << "[ttkRangePolygon] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
