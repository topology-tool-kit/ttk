#include <ttkMacros.h>
#include <ttkQuadrangulationSubdivision.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkQuadrangulationSubdivision);

ttkQuadrangulationSubdivision::ttkQuadrangulationSubdivision() {
  // MSC quadrangulation + initial 2D mesh
  SetNumberOfInputPorts(2);
  // quad mesh (containing ttkVertexIdentifiers of critical points)
  SetNumberOfOutputPorts(1);
}

int ttkQuadrangulationSubdivision::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) { // input quadrangulation
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 1) { // triangulated domain
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkQuadrangulationSubdivision::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkQuadrangulationSubdivision::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  auto quads = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto mesh = vtkDataSet::GetData(inputVector[1]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  auto triangulation = ttkAlgorithm::GetTriangulation(mesh);
  if(triangulation == nullptr) {
    return 0;
  }
  this->preconditionTriangulation(triangulation);

  auto inputCells = quads->GetCells();
  if(inputCells == nullptr || inputCells->GetData() == nullptr) {
    this->printErr("Invalid input quadrangle cells");
    return 0;
  }

  auto inputPoints = quads->GetPoints();
  auto pointData = quads->GetPointData();
  if(inputPoints == nullptr || inputPoints->GetData() == nullptr
     || pointData == nullptr) {
    this->printErr("Invalid input quadrangle points");
    return 0;
  }

  auto identifiers = pointData->GetArray(
    static_cast<const char *>(ttk::VertexScalarFieldName));
  if(identifiers == nullptr) {
    this->printErr("Missing point data array named "
                   + std::string(ttk::VertexScalarFieldName));
  }

  this->setInputQuads(ttkUtils::GetVoidPointer(inputCells->GetData()),
                      inputCells->GetNumberOfCells());
  this->setInputVertices(ttkUtils::GetVoidPointer(inputPoints->GetData()),
                         inputPoints->GetNumberOfPoints());
  this->setInputVertexIdentifiers(
    ttkUtils::GetVoidPointer(identifiers), identifiers->GetNumberOfTuples());

#define QUADSUB_EXPLICIT_CALLS(TRIANGL_CASE, TRIANGL_TYPE)                  \
  case TRIANGL_CASE: {                                                      \
    const auto tri = static_cast<TRIANGL_TYPE *>(triangulation->getData()); \
    if(tri != nullptr) {                                                    \
      res = this->execute<TRIANGL_TYPE>(*tri);                              \
    }                                                                       \
    break;                                                                  \
  }

  int res{-1};
  switch(triangulation->getType()) {
    QUADSUB_EXPLICIT_CALLS(
      ttk::Triangulation::Type::EXPLICIT, ttk::ExplicitTriangulation);
    QUADSUB_EXPLICIT_CALLS(
      ttk::Triangulation::Type::IMPLICIT, ttk::ImplicitTriangulation);
    QUADSUB_EXPLICIT_CALLS(
      ttk::Triangulation::Type::PERIODIC, ttk::PeriodicImplicitTriangulation);
  }

  if(res != 0) {
    this->printWrn("Please increase the number of relaxation iterations, of "
                   "subdivision levels or consider another function (higher "
                   "eigenfunctions).");
    if(!ShowResError) {
      return 0;
    }
  }

  auto cells = vtkSmartPointer<vtkCellArray>::New();

  for(size_t i = 0; i < getQuadNumber(); i++) {
    cells->InsertNextCell(4, &this->getQuadBuf()[5 * i + 1]);
  }

  // update output: get quadrangle values
  output->SetCells(VTK_QUAD, cells);

  auto points = vtkSmartPointer<vtkPoints>::New();
  for(size_t i = 0; i < getPointsNumber(); ++i) {
    points->InsertNextPoint(&this->getPointsBuf()[3 * i]);
  }

  // update output: get quadrangle vertices
  output->SetPoints(points);

  // add data array of points valences
  auto valences = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  valences->SetName("Valence");
  ttkUtils::SetVoidArray(
    valences, outputValences_.data(), outputValences_.size(), 1);
  output->GetPointData()->AddArray(valences);

  // add data array of points infos
  auto infos = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  infos->SetName("Type");
  ttkUtils::SetVoidArray(
    infos, outputVertType_.data(), outputVertType_.size(), 1);
  output->GetPointData()->AddArray(infos);

  auto subd = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  subd->SetName("Subdivision");
  ttkUtils::SetVoidArray(
    subd, outputSubdivision_.data(), outputSubdivision_.size(), 1);
  output->GetPointData()->AddArray(subd);

  if(RelaxationIterations > 0) {

    // add data array of number of triangles checked
    auto trChecked = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
    trChecked->SetName("Triangles checked");
    ttkUtils::SetVoidArray(
      trChecked, trianglesChecked_.data(), trianglesChecked_.size(), 1);
    output->GetPointData()->AddArray(trChecked);

    // add data array of projection success
    auto projSucc = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
    projSucc->SetName("Projection");
    ttkUtils::SetVoidArray(
      projSucc, projSucceeded_.data(), projSucceeded_.size(), 1);
    output->GetPointData()->AddArray(projSucc);
  }

  if(QuadStatistics) {
    auto quadArea = vtkSmartPointer<vtkFloatArray>::New();
    quadArea->SetName("Quad Area");
    ttkUtils::SetVoidArray(quadArea, quadArea_.data(), quadArea_.size(), 1);
    output->GetCellData()->AddArray(quadArea);

    auto diagsRatio = vtkSmartPointer<vtkFloatArray>::New();
    diagsRatio->SetName("Diagonals Ratio");
    ttkUtils::SetVoidArray(
      diagsRatio, quadDiagsRatio_.data(), quadDiagsRatio_.size(), 1);
    output->GetCellData()->AddArray(diagsRatio);

    auto edgesRatio = vtkSmartPointer<vtkFloatArray>::New();
    edgesRatio->SetName("Edges Ratio");
    ttkUtils::SetVoidArray(
      edgesRatio, quadEdgesRatio_.data(), quadEdgesRatio_.size(), 1);
    output->GetCellData()->AddArray(edgesRatio);

    auto anglesRatio = vtkSmartPointer<vtkFloatArray>::New();
    anglesRatio->SetName("Angles Ratio");
    ttkUtils::SetVoidArray(
      anglesRatio, quadAnglesRatio_.data(), quadAnglesRatio_.size(), 1);
    output->GetCellData()->AddArray(anglesRatio);

    auto hausDist = vtkSmartPointer<vtkFloatArray>::New();
    hausDist->SetName("Hausdorff");
    ttkUtils::SetVoidArray(hausDist, hausdorff_.data(), hausdorff_.size(), 1);
    output->GetPointData()->AddArray(hausDist);
  }

  // shallow copy input field data
  output->GetFieldData()->ShallowCopy(mesh->GetFieldData());

  return 1;
}
