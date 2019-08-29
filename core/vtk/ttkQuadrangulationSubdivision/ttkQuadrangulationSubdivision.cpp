#include <ttkQuadrangulationSubdivision.h>

#define MODULE_S "[ttkQuadrangulationSubdivision] "
#define MODULE_ERROR_S MODULE_S "Error: "
#ifndef TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)    \
  if(COND) {                            \
    cerr << MODULE_ERROR_S MSG << endl; \
    return RET;                         \
  }
#else // TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)
#endif // TTK_ENABLE_KAMIKAZE

vtkStandardNewMacro(ttkQuadrangulationSubdivision);

ttkQuadrangulationSubdivision::ttkQuadrangulationSubdivision()
  : UseAllCores{true}, ThreadNumber{} {

  // MSC quadrangulation + initial 2D mesh
  SetNumberOfInputPorts(2);
  // quad mesh (containing ttkVertexIdentifiers of critical points)
  SetNumberOfOutputPorts(1);
}

int ttkQuadrangulationSubdivision::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  }
  return 1;
}

int ttkQuadrangulationSubdivision::getTriangulation(
  vtkUnstructuredGrid *const input) {

  auto trg = ttkTriangulation::getTriangulation(input);

  TTK_ABORT_KK(trg == nullptr, "input triangulation pointer is NULL.", -1);

  trg->setWrapper(this);
  baseWorker_.setWrapper(this);
  baseWorker_.setupTriangulation(trg);
  Modified();

  TTK_ABORT_KK(trg->isEmpty(), "ttkTriangulation allocation problem.", -2);

  return 0;
}

int ttkQuadrangulationSubdivision::getQuadVertices(
  vtkUnstructuredGrid *const input) {
  auto cells = input->GetCells();

  TTK_ABORT_KK(cells == nullptr, "invalid input quadrangle cells", -3);
  TTK_ABORT_KK(
    cells->GetData() == nullptr, "invalid input quadrangle cell data", -4);

  auto points = input->GetPoints();

  TTK_ABORT_KK(points == nullptr, "invalid input critical points", -5);
  TTK_ABORT_KK(
    points->GetData() == nullptr, "invalid input quadrangle cell data", -6);

  auto pointData = input->GetPointData();
  auto identifiers = pointData->GetArray(
    static_cast<const char *>(ttk::VertexScalarFieldName));

  TTK_ABORT_KK(pointData == nullptr, "invalid input quadrangle point data", -7);
  TTK_ABORT_KK(identifiers == nullptr,
               "invalid input quadrangle vertices identifiers", -8);

  baseWorker_.setInputQuads(
    cells->GetData()->GetVoidPointer(0), cells->GetNumberOfCells());
  baseWorker_.setInputVertices(
    points->GetData()->GetVoidPointer(0), points->GetNumberOfPoints());
  baseWorker_.setInputVertexIdentifiers(
    identifiers->GetVoidPointer(0), identifiers->GetNumberOfTuples());

  return 0;
}

int ttkQuadrangulationSubdivision::doIt(std::vector<vtkDataSet *> &inputs,
                                        std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  // ensure consistency of dependent options before passing them to
  // base layer
  if(LockAllInputVertices) {
    LockInputExtrema = true;
  }

  baseWorker_.setSubdivisionLevel(SubdivisionLevel);
  baseWorker_.setRelaxationIterations(RelaxationIterations);
  baseWorker_.setLockInputExtrema(LockInputExtrema);
  baseWorker_.setLockAllInputVertices(LockAllInputVertices);
  baseWorker_.setReverseProjection(ReverseProjection);
  baseWorker_.setShowResError(ShowResError);
  baseWorker_.setHausdorffLevel(HausdorffLevel);

  auto quads = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  auto mesh = vtkUnstructuredGrid::SafeDownCast(inputs[1]);
  auto output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  int res = 0;

  res += getTriangulation(mesh);

  TTK_ABORT_KK(res != 0, "Cannot get mesh triangulation", -1);

  res += getQuadVertices(quads);

  TTK_ABORT_KK(res != 0, "Cannot get quad vertices", -2);

  res += baseWorker_.execute();

  if(res != 0) {
    vtkWarningMacro(MODULE_ERROR_S
                    "Please increase the number of relaxation iterations, of "
                    "subdivision levels or consider another function (higher "
                    "eigenfunctions).");
    if(!ShowResError) {
      return res;
    }
  }

  auto &outVertexValences = baseWorker_.outputValences_;
  auto &outVertexType = baseWorker_.outputVertType_;
  auto &outSubdvisionLevel = baseWorker_.outputSubdivision_;

  auto cells = vtkSmartPointer<vtkCellArray>::New();

  for(size_t i = 0; i < baseWorker_.getQuadNumber(); i++) {
    cells->InsertNextCell(4, &baseWorker_.getQuadBuf()[5 * i + 1]);
  }

  // update output: get quadrangle values
  output->SetCells(VTK_QUAD, cells);

  auto points = vtkSmartPointer<vtkPoints>::New();
  for(size_t i = 0; i < baseWorker_.getPointsNumber(); ++i) {
    points->InsertNextPoint(&baseWorker_.getPointsBuf()[3 * i]);
  }

  // update output: get quadrangle vertices
  output->SetPoints(points);

  // add data array of points valences
  auto valences = vtkSmartPointer<vtkIntArray>::New();
  valences->SetName("Valence");
  valences->SetVoidArray(outVertexValences.data(), outVertexValences.size(), 1);
  output->GetPointData()->AddArray(valences);

  // add data array of points infos
  auto infos = vtkSmartPointer<vtkIntArray>::New();
  infos->SetName("Type");
  infos->SetVoidArray(outVertexType.data(), outVertexType.size(), 1);
  output->GetPointData()->AddArray(infos);

  auto subd = vtkSmartPointer<vtkIntArray>::New();
  subd->SetName("Subdivision");
  subd->SetVoidArray(outSubdvisionLevel.data(), outSubdvisionLevel.size(), 1);
  output->GetPointData()->AddArray(subd);

  if(RelaxationIterations > 0) {
    auto &trianglesChecked = baseWorker_.trianglesChecked_;
    auto &projSucceeded = baseWorker_.projSucceeded_;

    // add data array of number of triangles checked
    auto trChecked = vtkSmartPointer<vtkIntArray>::New();
    trChecked->SetName("Triangles checked");
    trChecked->SetVoidArray(
      trianglesChecked.data(), trianglesChecked.size(), 1);
    output->GetPointData()->AddArray(trChecked);

    // add data array of projection success
    auto projSucc = vtkSmartPointer<vtkIntArray>::New();
    projSucc->SetName("Projection");
    projSucc->SetVoidArray(projSucceeded.data(), projSucceeded.size(), 1);
    output->GetPointData()->AddArray(projSucc);
  }

  if(QuadStatistics) {
    auto quadArea = vtkSmartPointer<vtkFloatArray>::New();
    quadArea->SetName("Quad Area");
    quadArea->SetVoidArray(
      baseWorker_.quadArea_.data(), baseWorker_.quadArea_.size(), 1);
    output->GetCellData()->AddArray(quadArea);

    auto diagsRatio = vtkSmartPointer<vtkFloatArray>::New();
    diagsRatio->SetName("Diagonals Ratio");
    diagsRatio->SetVoidArray(baseWorker_.quadDiagsRatio_.data(),
                             baseWorker_.quadDiagsRatio_.size(), 1);
    output->GetCellData()->AddArray(diagsRatio);

    auto edgesRatio = vtkSmartPointer<vtkFloatArray>::New();
    edgesRatio->SetName("Edges Ratio");
    edgesRatio->SetVoidArray(baseWorker_.quadEdgesRatio_.data(),
                             baseWorker_.quadEdgesRatio_.size(), 1);
    output->GetCellData()->AddArray(edgesRatio);

    auto anglesRatio = vtkSmartPointer<vtkFloatArray>::New();
    anglesRatio->SetName("Angles Ratio");
    anglesRatio->SetVoidArray(baseWorker_.quadAnglesRatio_.data(),
                              baseWorker_.quadAnglesRatio_.size(), 1);
    output->GetCellData()->AddArray(anglesRatio);

    auto hausDist = vtkSmartPointer<vtkFloatArray>::New();
    hausDist->SetName("Hausdorff");
    hausDist->SetVoidArray(
      baseWorker_.hausdorff_.data(), baseWorker_.hausdorff_.size(), 1);
    output->GetPointData()->AddArray(hausDist);
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
