#include <ttkMacros.h>
#include <ttkMorseSmaleQuadrangulation.h>
#include <ttkUtils.h>

#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

vtkStandardNewMacro(ttkMorseSmaleQuadrangulation);

ttkMorseSmaleQuadrangulation::ttkMorseSmaleQuadrangulation() {
  // critical points + 1-separatrices + segmentation
  SetNumberOfInputPorts(3);
  // quad mesh (containing ttkVertexIdentifiers of critical points)
  SetNumberOfOutputPorts(1);
}

int ttkMorseSmaleQuadrangulation::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) { // critical points
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  } else if(port == 1) { // MSC separatrices
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  } else if(port == 2) { // triangulated domain
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkMorseSmaleQuadrangulation::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkMorseSmaleQuadrangulation::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  auto critpoints = vtkPointSet::GetData(inputVector[0]);
  auto seprs = vtkPolyData::GetData(inputVector[1]);
  auto domain = vtkDataSet::GetData(inputVector[2]);
  auto output = vtkPolyData::GetData(outputVector);

  auto triangulation = ttkAlgorithm::GetTriangulation(domain);
  if(triangulation == nullptr) {
    return 0;
  }
  this->preconditionTriangulation(triangulation);

  auto cpPoints = critpoints->GetPoints();
  auto cpData = critpoints->GetPointData();
  auto seprsPoints = seprs->GetPoints();
  auto seprsData = seprs->GetPointData();

  if(seprsPoints == nullptr || seprsData == nullptr || cpPoints == nullptr
     || cpData == nullptr) {
    this->printErr("Invalid input");
    return 0;
  }

  auto cpci = cpData->GetArray("CellId");
  auto cpcd = cpData->GetArray("CellDimension");
  auto cpid = cpData->GetArray(ttk::VertexScalarFieldName);
  auto sepid = seprsData->GetArray("CellId");
  auto sepdim = seprsData->GetArray("CellDimension");
  auto sepmask = seprsData->GetArray(ttk::MaskScalarFieldName);

  if(cpci == nullptr || cpcd == nullptr || cpid == nullptr || sepid == nullptr
     || sepdim == nullptr || sepmask == nullptr) {
    this->printErr("Missing data arrays");
    return 0;
  }

  this->setCriticalPoints(
    cpPoints->GetNumberOfPoints(), ttkUtils::GetVoidPointer(cpPoints),
    ttkUtils::GetVoidPointer(cpid), ttkUtils::GetVoidPointer(cpci),
    ttkUtils::GetVoidPointer(cpcd));

  this->setSeparatrices(
    sepid->GetNumberOfTuples(), ttkUtils::GetVoidPointer(sepid),
    ttkUtils::GetVoidPointer(sepdim), ttkUtils::GetVoidPointer(sepmask),
    ttkUtils::GetVoidPointer(seprsPoints));

  int res{-1};
  ttkTemplateMacro(
    triangulation->getType(),
    res = this->execute(*static_cast<TTK_TT *>(triangulation->getData())););

  if(res != 0) {
    this->printWrn("Consider another (eigen) function, persistence threshold "
                   "or refine your input triangulation");
    if(!ShowResError) {
      return 0;
    }
  }

  // output points: critical points + generated separatrices middles
  vtkNew<vtkPoints> outQuadPoints{};
  for(size_t i = 0; i < outputPoints_.size() / 3; i++) {
    outQuadPoints->InsertNextPoint(&outputPoints_[3 * i]);
  }
  output->SetPoints(outQuadPoints);

  // quad vertices identifiers
  vtkNew<ttkSimplexIdTypeArray> identifiers{};
  identifiers->SetName(ttk::VertexScalarFieldName);
  ttkUtils::SetVoidArray(
    identifiers, outputPointsIds_.data(), outputPointsIds_.size(), 1);
  output->GetPointData()->AddArray(identifiers);

  // quad vertices type
  vtkNew<ttkSimplexIdTypeArray> type{};
  type->SetName("QuadVertType");
  ttkUtils::SetVoidArray(
    type, outputPointsTypes_.data(), outputPointsTypes_.size(), 1);
  output->GetPointData()->AddArray(type);

  // quad vertices cells
  vtkNew<ttkSimplexIdTypeArray> cellid{};
  cellid->SetName("QuadCellId");
  ttkUtils::SetVoidArray(
    cellid, outputPointsCells_.data(), outputPointsCells_.size(), 1);
  output->GetPointData()->AddArray(cellid);

  // vtkCellArray of quadrangle values containing outArray
  vtkNew<vtkCellArray> cells{};
  for(size_t i = 0; i < outputCells_.size(); i++) {
    cells->InsertNextCell(4, outputCells_[i].data());
  }

  // update output: get quadrangle values
  output->SetPolys(cells);

  // shallow copy input field data
  output->GetFieldData()->ShallowCopy(domain->GetFieldData());

  return 1;
}
