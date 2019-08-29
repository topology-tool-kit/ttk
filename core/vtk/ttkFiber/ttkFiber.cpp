#include <ttkFiber.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkFiber)

  ttkFiber::ttkFiber() {
  Uvalue = 0;
  Vvalue = 0;

  UseAllCores = true;
}

ttkFiber::~ttkFiber() {
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkFiber::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkFiber::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkFiber] " << progress * 100 << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkFiber::doIt(vtkDataSet *input, vtkPolyData *output) {

  Timer t;

  vtkDataArray *uField = input->GetPointData()->GetArray(Ucomponent.data());

  if(!uField) {
    stringstream msg;
    msg << "[ttkFiber] Error: cannot find field '" << Ucomponent << "'!"
        << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
    return -1;
  }

  vtkDataArray *vField = input->GetPointData()->GetArray(Vcomponent.data());

  if(!vField) {
    stringstream msg;
    msg << "[ttkFiber] Error: cannot find field '" << Vcomponent << "'!"
        << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
    return -2;
  }

  vtkSmartPointer<vtkContourFilter> isoSurface
    = vtkSmartPointer<vtkContourFilter>::New();

  isoSurface->SetInputData(input);
  isoSurface->SetComputeScalars(true);
  isoSurface->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, Ucomponent.data());
  isoSurface->SetGenerateTriangles(true);
  isoSurface->SetNumberOfContours(1);
  isoSurface->SetValue(0, Uvalue);
  isoSurface->Update();

  vtkDataSet *surface = isoSurface->GetOutput();

  vField = surface->GetPointData()->GetArray(Vcomponent.data());

  if(!vField) {
    stringstream msg;
    msg << "[ttkFiber] Error: cannot find field '" << Vcomponent << "'!"
        << endl;
    msg << "[ttkFiber] U-coordinate may be out of bounds..." << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
    return -3;
  }

  vtkSmartPointer<vtkContourFilter> isoLine
    = vtkSmartPointer<vtkContourFilter>::New();

  isoLine->SetInputData(surface);
  isoLine->SetComputeScalars(true);
  isoLine->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, Vcomponent.data());
  isoLine->SetNumberOfContours(1);
  isoLine->SetValue(0, Vvalue);
  isoLine->Update();

  output->ShallowCopy(isoLine->GetOutput());

  {
    stringstream msg;
    msg << "[ttkFiber] Fiber extracted in " << t.getElapsedTime() << " s. ("
        << output->GetNumberOfCells() << " edge(s))" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkFiber::RequestData(vtkInformation *request,
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector) {

  Memory m;

  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkFiber] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
