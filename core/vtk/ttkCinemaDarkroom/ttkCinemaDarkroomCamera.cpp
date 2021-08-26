#include <ttkCinemaDarkroomCamera.h>

#include <vtkCamera.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

// #include <pqActiveObjects.h>
// #include <vtkSMRenderViewProxy.h>

#include <vtkPythonInterpreter.h>

#include <ttkUtils.h>

vtkStandardNewMacro(ttkCinemaDarkroomCamera);

ttkCinemaDarkroomCamera::ttkCinemaDarkroomCamera() {
  this->setDebugMsgPrefix("CinemaDarkroomCamera");

  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaDarkroomCamera::~ttkCinemaDarkroomCamera() {
}

int ttkCinemaDarkroomCamera::FillInputPortInformation(int, vtkInformation *) {
  return 0;
}

int ttkCinemaDarkroomCamera::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomCamera::SyncWithParaViewCamera() {
  ttk::Timer timer;
  this->printMsg(
    "Updating Camera Parameters", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  std::string code(R"(
from paraview.simple import GetActiveView
from paraview.simple import FindSource

self = FindSource("TTKDarkroomCamera1")

view = GetActiveView()
if view and self:
  self.Position = view.CameraPosition
  self.Up = view.CameraViewUp
  self.FocalPoint = view.CameraFocalPoint
)");

  vtkPythonInterpreter::RunSimpleString(code.data());

  //   auto& activeObjects = pqActiveObjects::instance();

  //   auto view = activeObjects.activeView();
  //   if(!view){
  //     this->printErr("Unable to retrieve active view.");
  //     return 0;
  //   }

  //   auto proxy = dynamic_cast<vtkSMRenderViewProxy*>(view->getViewProxy());
  //   if(!proxy) {
  //     this->printErr("Unable to cast to vtkSMRenderViewProxy.");
  //     return 0;
  //   }

  //   auto camera = proxy->GetActiveCamera();

  //   this->SetUp( camera->GetViewUp() );
  //   this->SetFocalPoint( camera->GetFocalPoint() );
  //   this->SetPosition( camera->GetPosition() );

  this->printMsg("Updating Camera Parameters", 1, timer.getElapsedTime(), 1);

  this->Modified();

  return 1;
}

int ttkCinemaDarkroomCamera::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **ttkNotUsed(inputVector),
  vtkInformationVector *outputVector) {

  ttk::Timer timer;
  this->printMsg("Generating Camera", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  auto output = vtkUnstructuredGrid::GetData(outputVector);

  // Points
  {
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(this->Position);
    output->SetPoints(points);
  }

  // Cells
  {
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
    {
      offsetArray->SetNumberOfTuples(2);
      auto offsetArrayData
        = static_cast<int *>(ttkUtils::GetVoidPointer(offsetArray));
      offsetArrayData[0] = 0;
      offsetArrayData[1] = 1;
    }

    auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
    {
      connectivityArray->SetNumberOfTuples(1);
      auto connectivityArrayData
        = static_cast<int *>(ttkUtils::GetVoidPointer(connectivityArray));
      connectivityArrayData[0] = 0;
    }

    cells->SetData(offsetArray, connectivityArray);
    output->SetCells(VTK_VERTEX, cells);
  }

  // Point Data
  auto generateArray
    = [](vtkPointData *pd, const std::string &name, const double *data) {
        auto array = vtkSmartPointer<vtkDoubleArray>::New();
        array->SetName(name.data());
        array->SetNumberOfComponents(3);
        array->SetNumberOfTuples(1);
        auto arrayData = static_cast<double *>(ttkUtils::GetVoidPointer(array));
        arrayData[0] = data[0];
        arrayData[1] = data[1];
        arrayData[2] = data[2];
        pd->AddArray(array);
      };
  auto pd = output->GetPointData();
  generateArray(pd, "CamUp", this->Up);
  generateArray(pd, "CamFocalPoint", this->FocalPoint);

  this->printMsg("Generating Camera", 1, timer.getElapsedTime(), 1);

  return 1;
}
