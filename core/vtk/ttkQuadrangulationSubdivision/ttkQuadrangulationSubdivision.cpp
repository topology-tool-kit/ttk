#include <ttkQuadrangulationSubdivision.h>

#define MODULE_S "[ttkQuadrangulationSubdivision] "
#define MODULE_ERROR_S MODULE_S "Error: "

vtkStandardNewMacro(ttkQuadrangulationSubdivision);

ttkQuadrangulationSubdivision::ttkQuadrangulationSubdivision()
  : UseAllCores{true}, ThreadNumber{}, ForceInputIdentifiersField{false},
    ForceInputOffsetIdentifiersField{false}, SubdivisionLevel{5},
    RelaxationIterations{100} {

  InputIdentifiersFieldName
    = std::string(static_cast<const char *>(ttk::VertexScalarFieldName));

  // critical points + 1-separatrices
  SetNumberOfInputPorts(2);
  // quad mesh (containing ttkVertexIdentifiers of critical points)
  SetNumberOfOutputPorts(1);
}

int ttkQuadrangulationSubdivision::FillInputPortInformation(int port,
                                                        vtkInformation *info) {

  if(port == 0 || port == 1) {
    // Morse-Smale Complex output
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  }
  return 0;
}

int ttkQuadrangulationSubdivision::getCriticalPoints(vtkUnstructuredGrid *input) {

  // store handler to points
  auto cp = input->GetPoints();

  // store handle to points identifiers
  auto vsfn = static_cast<const char *>(ttk::VertexScalarFieldName);
  auto pointData = input->GetPointData();
  auto cpi = pointData->GetArray(vsfn);
  auto cpci = pointData->GetArray("CellId");
  auto cpType = pointData->GetArray("CellDimension");

#ifndef TTK_ENABLE_KAMIKAZE
  if(cp == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale critical points" << endl;
    return -1;
  }
  if(cpi == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale critical points identifiers"
         << endl;
    return -2;
  }
  if(cpci == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale critical points cell identifiers"
         << endl;
    return -3;
  }
  if(cpType == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale critical points type" << endl;
    return -4;
  }
#endif // TTK_ENABLE_KAMIKAZE

  quadrangulationSubdivision_.setCriticalPointsNumber(cp->GetNumberOfPoints());
  quadrangulationSubdivision_.setCriticalPoints(cp->GetVoidPointer(0));
  quadrangulationSubdivision_.setCriticalPointsIdentifiers(cpi->GetVoidPointer(0));
  quadrangulationSubdivision_.setCriticalPointsCellIds(cpci->GetVoidPointer(0));
  quadrangulationSubdivision_.setCriticalPointsType(cpType->GetVoidPointer(0));

  return 0;
}

int ttkQuadrangulationSubdivision::getSeparatrices(vtkUnstructuredGrid *input) {

  // get separatrices points
  auto separatrices = input->GetPoints();

  auto cellData = input->GetCellData();
  auto sepId = cellData->GetArray("SeparatrixId");
  auto sepSourceId = cellData->GetArray("SourceId");
  auto sepDestId = cellData->GetArray("DestinationId");

#ifndef TTK_ENABLE_KAMIKAZE
  if(separatrices == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale separatrices points." << endl;
    return -1;
  }
  if(sepId == nullptr) {
    cerr << MODULE_ERROR_S "wrong separatrices id." << endl;
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

  quadrangulationSubdivision_.setSeparatriceNumber(sepId->GetNumberOfValues());
  quadrangulationSubdivision_.setSepId(sepId->GetVoidPointer(0));
  quadrangulationSubdivision_.setSepSourceId(sepSourceId->GetVoidPointer(0));
  quadrangulationSubdivision_.setSepDestId(sepDestId->GetVoidPointer(0));

  return 0;
}

int ttkQuadrangulationSubdivision::doIt(std::vector<vtkDataSet *> &inputs,
                                    std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  quadrangulationSubdivision_.setSubdivisionLevel(SubdivisionLevel);
  quadrangulationSubdivision_.setRelaxationIterations(RelaxationIterations);

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

  quadrangulationSubdivision_.setOutputCells(&outQuadrangles_);

  res += quadrangulationSubdivision_.execute();

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in base/QuadrangulationSubdivision
  if(res != 0) {
    cerr << MODULE_S "QuadrangulationSubdivision.execute() error code: " << res
         << endl;
    return -9;
  }
#endif

  // update result: get critical points from input
  output->SetPoints(cp->GetPoints());
  // number of data in outTempVector
  size_t nOutData = outQuadrangles_.size();

  // vtkDataArray containing quadrangles values
  auto outArray = vtkSmartPointer<vtkIdTypeArray>::New();

  // outArray points to outQuadrangles_ data, but does not deallocate it
  outArray->SetArray(outQuadrangles_.data(), nOutData, 1);

  // vtkCellArray of quadrangle values containing outArray
  auto cells = vtkSmartPointer<vtkCellArray>::New();
  cells->SetCells(nOutData, outArray);

  // update output: get quadrangle values
  output->SetCells(VTK_QUAD, cells);

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
