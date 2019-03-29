#include <ttkSurfaceQuadrangulation.h>

#define MODULE_S "[ttkSurfaceQuadrangulation] "
#define MODULE_ERROR_S MODULE_S "Error: "

vtkStandardNewMacro(ttkSurfaceQuadrangulation);

ttkSurfaceQuadrangulation::ttkSurfaceQuadrangulation()
  : UseAllCores{true}, ThreadNumber{} {

  // critical points + 1-separatrices
  SetNumberOfInputPorts(2);
  // quad mesh (containing ttkVertexIdentifiers of critical points)
  SetNumberOfOutputPorts(1);
}

int ttkSurfaceQuadrangulation::FillInputPortInformation(int port,
                                                        vtkInformation *info) {

  if(port == 0 || port == 1) {
    // Morse-Smale Complex output
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  }
  return 0;
}

int ttkSurfaceQuadrangulation::getCriticalPoints(vtkUnstructuredGrid *input) {

  // store handler to points
  auto cp = input->GetPoints();

  // store handle to points identifiers
  auto pointData = input->GetPointData();
  auto cpci = pointData->GetArray("CellId");

#ifndef TTK_ENABLE_KAMIKAZE
  if(cp == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale critical points" << endl;
    return -1;
  }
  if(cpci == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale critical points cell identifiers"
         << endl;
    return -3;
  }
#endif // TTK_ENABLE_KAMIKAZE

  surfaceQuadrangulation_.setCriticalPointsNumber(cp->GetNumberOfPoints());
  surfaceQuadrangulation_.setCriticalPointsCellIds(cpci->GetVoidPointer(0));

  return 0;
}

int ttkSurfaceQuadrangulation::getSeparatrices(vtkUnstructuredGrid *input) {

  // get separatrices points
  auto separatrices = input->GetPoints();

  auto cellData = input->GetCellData();
  auto sepSourceId = cellData->GetArray("SourceId");
  auto sepDestId = cellData->GetArray("DestinationId");

#ifndef TTK_ENABLE_KAMIKAZE
  if(separatrices == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale separatrices points." << endl;
    return -1;
  }
  if(sepSourceId == nullptr) {
    cerr << MODULE_ERROR_S "wrong separatrices source id." << endl;
    return -1;
  }
  if(sepDestId == nullptr) {
    cerr << MODULE_ERROR_S "wrong separatrices destination id." << endl;
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  surfaceQuadrangulation_.setSeparatrices(sepSourceId->GetNumberOfValues(),
                                          sepSourceId->GetVoidPointer(0),
                                          sepDestId->GetVoidPointer(0));
  return 0;
}

int ttkSurfaceQuadrangulation::doIt(std::vector<vtkDataSet *> &inputs,
                                    std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  auto cp = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  auto spr = vtkUnstructuredGrid::SafeDownCast(inputs[1]);
  auto output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  int res = 0;

  res += getCriticalPoints(cp);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << MODULE_ERROR_S "wrong points." << endl;
    return -1;
  }
#endif

  res += getSeparatrices(spr);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << MODULE_ERROR_S "wrong separatrices." << endl;
    return -1;
  }
#endif

  surfaceQuadrangulation_.setOutputCells(&outQuadrangles_);

  res += surfaceQuadrangulation_.execute();

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in base/SurfaceQuadrangulation
  if(res != 0) {
    cerr << MODULE_S "SurfaceQuadrangulation.execute() error code: " << res
         << endl;
    return -9;
  }
#endif

  // update result: get critical points from input
  output->ShallowCopy(cp);

  // vtkCellArray of quadrangle values containing outArray
  auto cells = vtkSmartPointer<vtkCellArray>::New();
  for(size_t i = 0; i < outQuadrangles_.size() / 5; i++) {
    cells->InsertNextCell(4, &outQuadrangles_[5 * i + 1]);
  }

  // update output: get quadrangle values
  output->SetCells(VTK_QUAD, cells);

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
