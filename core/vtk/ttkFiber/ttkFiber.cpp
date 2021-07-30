#include <ttkFiber.h>

#include <vtkContourFilter.h>
#include <vtkInformation.h>

vtkStandardNewMacro(ttkFiber);

ttkFiber::ttkFiber() {
  this->setDebugMsgPrefix("Fiber");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkFiber::~ttkFiber() {
}

int ttkFiber::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }

  return 0;
}

int ttkFiber::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }

  return 0;
}

int ttkFiber::RequestData(vtkInformation *ttkNotUsed(request),
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector) {
  ttk::Timer t;

  this->printMsg("Computing Fiber", 0, 0, ttk::debug::LineMode::REPLACE);

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkPolyData::GetData(outputVector);

  auto uArray = this->GetInputArrayToProcess(0, inputVector);
  if(!uArray) {
    this->printErr("Unable to retrieve input array U.");
    return 0;
  }
  std::string uArrayName(uArray->GetName());

  auto vArray = this->GetInputArrayToProcess(1, inputVector);
  if(!vArray) {
    this->printErr("Unable to retrieve input array V.");
    return 0;
  }
  std::string vArrayName(vArray->GetName());

  this->printMsg("Computing Fiber (" + uArrayName + ": "
                   + std::to_string(UValue) + ", " + vArrayName + ": "
                   + std::to_string(VValue) + ")",
                 0.1, t.getElapsedTime(), ttk::debug::LineMode::REPLACE);

  auto isoSurface = vtkSmartPointer<vtkContourFilter>::New();
  isoSurface->SetInputData(input);
  isoSurface->SetComputeScalars(true);
  isoSurface->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, uArray->GetName());
  isoSurface->SetGenerateTriangles(true);
  isoSurface->SetNumberOfContours(1);
  isoSurface->SetValue(0, UValue);
  isoSurface->Update();

  auto isoLine = vtkSmartPointer<vtkContourFilter>::New();
  isoLine->SetInputData(isoSurface->GetOutput());
  isoLine->SetComputeScalars(true);
  isoLine->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, vArray->GetName());
  isoLine->SetNumberOfContours(1);
  isoLine->SetValue(0, VValue);
  isoLine->Update();

  output->ShallowCopy(isoLine->GetOutput());

  this->printMsg("Computing Fiber (" + uArrayName + ": "
                   + std::to_string(UValue) + ", " + vArrayName + ": "
                   + std::to_string(VValue) + ")",
                 1, t.getElapsedTime());

  return 1;
}