#include <ttkBarycentricSubdivision.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkBarycentricSubdivision);

ttkBarycentricSubdivision::ttkBarycentricSubdivision()
  : BarycentricSubdivision(points_, cells_, pointId_, pointDim_) {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkBarycentricSubdivision::~ttkBarycentricSubdivision() {
}

int ttkBarycentricSubdivision::FillInputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  else
    return 0;

  return 1;
}

int ttkBarycentricSubdivision::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;

  return 1;
}

vtkSmartPointer<vtkDataArray> ttkBarycentricSubdivision::AllocateScalarField(
  vtkDataArray *const inputScalarField, int ntuples) const {

  vtkSmartPointer<vtkDataArray> res{};

  // allocate the memory for the output scalar field
  switch(inputScalarField->GetDataType()) {
    case VTK_CHAR:
    case VTK_DOUBLE:
    case VTK_FLOAT:
    case VTK_INT:
    case VTK_ID_TYPE:
    case VTK_LONG:
      res = inputScalarField->NewInstance();
      break;
    default:
      this->printErr("Unsupported data array type");
      break;
  }
  res->SetNumberOfComponents(1);
  res->SetNumberOfTuples(ntuples);
  res->SetName(inputScalarField->GetName());
  return res;
}

int ttkBarycentricSubdivision::InterpolateScalarFields(
  vtkUnstructuredGrid *const input, vtkUnstructuredGrid *const output) const {

  const size_t npointdata = input->GetPointData()->GetNumberOfArrays();
  const size_t ncelldata = input->GetCellData()->GetNumberOfArrays();

  const auto outPointsNumber = this->getNumberOfVertices();

  for(size_t i = 0; i < npointdata; ++i) {
    auto inputScalarField = input->GetPointData()->GetArray(i);
    if(inputScalarField == nullptr) {
      return -2;
    }

#define DISPATCH_INTERPOLATE_DIS(CASE, TYPE)                             \
  case CASE:                                                             \
    this->interpolateDiscreteScalarField<TYPE>(                          \
      static_cast<TYPE *>(ttkUtils::GetVoidPointer(inputScalarField)),   \
      static_cast<TYPE *>(ttkUtils::GetVoidPointer(outputScalarField))); \
    break
#define DISPATCH_INTERPOLATE_CONT(CASE, TYPE)                            \
  case CASE:                                                             \
    this->interpolateContinuousScalarField<TYPE>(                        \
      static_cast<TYPE *>(ttkUtils::GetVoidPointer(inputScalarField)),   \
      static_cast<TYPE *>(ttkUtils::GetVoidPointer(outputScalarField))); \
    break

    auto outputScalarField
      = AllocateScalarField(inputScalarField, outPointsNumber);
    if(outputScalarField == nullptr) {
      return -3;
    }

    // only for scalar fields
    switch(inputScalarField->GetDataType()) {
      DISPATCH_INTERPOLATE_DIS(VTK_CHAR, char);
      DISPATCH_INTERPOLATE_DIS(VTK_INT, int);
      DISPATCH_INTERPOLATE_DIS(VTK_LONG, long);
      DISPATCH_INTERPOLATE_DIS(VTK_ID_TYPE, vtkIdType);
      DISPATCH_INTERPOLATE_CONT(VTK_FLOAT, float);
      DISPATCH_INTERPOLATE_CONT(VTK_DOUBLE, double);
    }
    output->GetPointData()->AddArray(outputScalarField);
  }

  const auto outCellsNumber = this->getNumberOfTriangles();

  for(size_t i = 0; i < ncelldata; ++i) {
    auto inputScalarField = input->GetCellData()->GetArray(i);
    if(inputScalarField == nullptr) {
      return -2;
    }

    auto outputScalarField
      = AllocateScalarField(inputScalarField, outCellsNumber);
    if(outputScalarField == nullptr) {
      return -3;
    }

    // only for scalar fields
    switch(inputScalarField->GetDataType()) {
      vtkTemplateMacro(this->interpolateCellDataField<VTK_TT>(
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalarField)),
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalarField))));
    }
    output->GetCellData()->AddArray(outputScalarField);
  }

  return 0;
}

int ttkBarycentricSubdivision::RequestData(vtkInformation *request,
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector) {

  ttk::Memory m;
  ttk::Timer t;

  auto input = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  // early return: copy input if no subdivision
  if(SubdivisionLevel == 0) {
    output->ShallowCopy(input);
    return 0;
  }

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    return -1;
  }
  this->setupTriangulation(triangulation);

  ttk::Triangulation triangulationSubdivision;
  this->setOutputTriangulation(&triangulationSubdivision);

  this->setInputPoints(ttkUtils::GetVoidPointer(input->GetPoints()));

  // first iteration: generate the new triangulation
  this->execute();

  // first iteration: interpolate input scalar fields
  int ret = InterpolateScalarFields(input, output);
  if(ret < 0) {
    this->printErr("Error interpolating input data array(s)");
    return -1;
  }

  for(unsigned int i = 1; i < SubdivisionLevel; ++i) {
    // move previous points to temp vector
    decltype(points_) tmpPoints{};
    std::swap(points_, tmpPoints);
    this->setInputPoints(tmpPoints.data());

    // move previous triangulation cells to temp vector
    decltype(cells_) tmpCells{};
    std::swap(cells_, tmpCells);

    // move previous triangulation to temp triangulation
    decltype(triangulationSubdivision) tmpTr{};
    std::swap(triangulationSubdivision, tmpTr);

    tmpTr.setInputCells(tmpCells.size() / 4, tmpCells.data());
    tmpTr.setInputPoints(tmpPoints.size() / 3, tmpPoints.data());
    this->setupTriangulation(&tmpTr);
    this->setOutputTriangulation(&triangulationSubdivision);

    // generate the new triangulation
    this->execute();

    // temporary vtkUnstructuredGrid moved from output
    vtkSmartPointer<vtkUnstructuredGrid> tmp(std::move(output));
    // interpolate from tmp to output
    InterpolateScalarFields(tmp, output);
  }

  // generated 3D coordinates
  auto points = vtkSmartPointer<vtkPoints>::New();
  for(size_t i = 0; i < points_.size() / 3; i++) {
    points->InsertNextPoint(&points_[3 * i]);
  }
  output->SetPoints(points);

  // generated triangles
  const size_t dataPerCell = 4;
  auto cells = vtkSmartPointer<vtkCellArray>::New();
  for(size_t i = 0; i < cells_.size() / dataPerCell; i++) {
    cells->InsertNextCell(3, &cells_[dataPerCell * i + 1]);
  }
  output->SetCells(VTK_TRIANGLE, cells);

  // cell id
  auto cellId = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  cellId->SetName("CellId");
  ttkUtils::SetVoidArray(cellId, pointId_.data(), pointId_.size(), 1);
  output->GetPointData()->AddArray(cellId);

  // cell dimension
  auto cellDim = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  cellDim->SetName("CellDimension");
  ttkUtils::SetVoidArray(cellDim, pointDim_.data(), pointDim_.size(), 1);
  output->GetPointData()->AddArray(cellDim);

  // {
  //   std::stringstream msg;
  //   msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." <<
  //   endl; dMsg(std::cout, msg.str(), memoryMsg);
  // }

  return 1;
}
