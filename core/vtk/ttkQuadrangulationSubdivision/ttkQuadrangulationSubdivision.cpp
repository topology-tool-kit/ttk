#include <ttkMacros.h>
#include <ttkQuadrangulationSubdivision.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

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
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
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
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkQuadrangulationSubdivision::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  auto quads = vtkPolyData::GetData(inputVector[0]);
  auto mesh = vtkDataSet::GetData(inputVector[1]);
  auto output = vtkPolyData::GetData(outputVector);

  auto triangulation = ttkAlgorithm::GetTriangulation(mesh);
  if(triangulation == nullptr) {
    return 0;
  }
  this->preconditionTriangulation(triangulation);

  auto inputCells = quads->GetPolys();
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
    return 0;
  }

  this->setInputQuads(
    // get quads from PolyData's connectivity array
    ttkUtils::GetVoidPointer(inputCells->GetConnectivityArray()),
    inputCells->GetNumberOfCells());
  this->setInputVertices(ttkUtils::GetVoidPointer(inputPoints->GetData()),
                         inputPoints->GetNumberOfPoints());
  this->setInputVertexIdentifiers(
    ttkUtils::GetVoidPointer(identifiers), identifiers->GetNumberOfTuples());

  int res{-1};
  ttkTemplateMacro(
    triangulation->getType(),
    res = this->execute(*static_cast<TTK_TT *>(triangulation->getData())));

  if(res != 0) {
    this->printWrn("Please increase the number of relaxation iterations, of "
                   "subdivision levels or consider another function (higher "
                   "eigenfunctions).");
    if(!ShowResError) {
      return 0;
    }
  }

  vtkNew<vtkCellArray> cells{};

  for(size_t i = 0; i < outputQuads_.size(); i++) {
    cells->InsertNextCell(4, this->outputQuads_[i].data());
  }

  // update output: get quadrangle values
  output->SetPolys(cells);

  vtkNew<vtkPoints> points{};
  for(size_t i = 0; i < outputPoints_.size(); ++i) {
    points->InsertNextPoint(&outputPoints_[i].x);
  }

  // update output: get quadrangle vertices
  output->SetPoints(points);

  // add data array of points valences
  vtkNew<ttkSimplexIdTypeArray> valences{};
  valences->SetName("Valence");
  ttkUtils::SetVoidArray(
    valences, outputValences_.data(), outputValences_.size(), 1);
  output->GetPointData()->AddArray(valences);

  // add data array of points infos
  vtkNew<ttkSimplexIdTypeArray> infos{};
  infos->SetName("Type");
  ttkUtils::SetVoidArray(
    infos, outputVertType_.data(), outputVertType_.size(), 1);
  output->GetPointData()->AddArray(infos);

  vtkNew<ttkSimplexIdTypeArray> subd{};
  subd->SetName("Subdivision");
  ttkUtils::SetVoidArray(
    subd, outputSubdivision_.data(), outputSubdivision_.size(), 1);
  output->GetPointData()->AddArray(subd);

  if(RelaxationIterations > 0) {

    // add data array of number of triangles checked
    vtkNew<ttkSimplexIdTypeArray> trChecked{};
    trChecked->SetName("Triangles checked");
    ttkUtils::SetVoidArray(
      trChecked, trianglesChecked_.data(), trianglesChecked_.size(), 1);
    output->GetPointData()->AddArray(trChecked);

    // add data array of projection success
    vtkNew<ttkSimplexIdTypeArray> projSucc{};
    projSucc->SetName("Projection");
    ttkUtils::SetVoidArray(
      projSucc, projSucceeded_.data(), projSucceeded_.size(), 1);
    output->GetPointData()->AddArray(projSucc);
  }

  if(QuadStatistics) {
    vtkNew<vtkFloatArray> quadArea{};
    quadArea->SetName("Quad Area");
    ttkUtils::SetVoidArray(quadArea, quadArea_.data(), quadArea_.size(), 1);
    output->GetCellData()->AddArray(quadArea);

    vtkNew<vtkFloatArray> diagsRatio{};
    diagsRatio->SetName("Diagonals Ratio");
    ttkUtils::SetVoidArray(
      diagsRatio, quadDiagsRatio_.data(), quadDiagsRatio_.size(), 1);
    output->GetCellData()->AddArray(diagsRatio);

    vtkNew<vtkFloatArray> edgesRatio{};
    edgesRatio->SetName("Edges Ratio");
    ttkUtils::SetVoidArray(
      edgesRatio, quadEdgesRatio_.data(), quadEdgesRatio_.size(), 1);
    output->GetCellData()->AddArray(edgesRatio);

    vtkNew<vtkFloatArray> anglesRatio{};
    anglesRatio->SetName("Angles Ratio");
    ttkUtils::SetVoidArray(
      anglesRatio, quadAnglesRatio_.data(), quadAnglesRatio_.size(), 1);
    output->GetCellData()->AddArray(anglesRatio);

    vtkNew<vtkFloatArray> hausDist{};
    hausDist->SetName("Hausdorff");
    ttkUtils::SetVoidArray(hausDist, hausdorff_.data(), hausdorff_.size(), 1);
    output->GetPointData()->AddArray(hausDist);
  }

  // shallow copy input field data
  output->GetFieldData()->ShallowCopy(mesh->GetFieldData());

  return 1;
}
