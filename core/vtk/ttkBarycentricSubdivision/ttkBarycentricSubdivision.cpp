#include <ttkBarycentricSubdivision.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkBarycentricSubdivision);

ttkBarycentricSubdivision::ttkBarycentricSubdivision() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkBarycentricSubdivision::FillInputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkBarycentricSubdivision::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

vtkSmartPointer<vtkDataArray> ttkBarycentricSubdivision::AllocateScalarField(
  vtkDataArray *const inputScalarField, int ntuples) const {

  vtkSmartPointer<vtkDataArray> res;

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
  vtkDataSet *const input,
  vtkUnstructuredGrid *const output,
  ttk::Triangulation &inputTriangulation) const {

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
#define DISPATCH_INTERPOLATE_CONT(CASE, TYPE)                                 \
  case CASE:                                                                  \
    switch(inputTriangulation.getType()) {                                    \
      BARYSUBD_TRIANGL_CALLS(                                                 \
        TYPE, ttk::Triangulation::Type::EXPLICIT, ttk::ExplicitTriangulation) \
      BARYSUBD_TRIANGL_CALLS(                                                 \
        TYPE, ttk::Triangulation::Type::IMPLICIT, ttk::ImplicitTriangulation) \
      BARYSUBD_TRIANGL_CALLS(                                                 \
        TYPE, ttk::Triangulation::Type::COMPACT, ttk::CompactTriangulation)   \
      BARYSUBD_TRIANGL_CALLS(TYPE, ttk::Triangulation::Type::PERIODIC,        \
                             ttk::PeriodicImplicitTriangulation)              \
    }                                                                         \
    break;
#define BARYSUBD_TRIANGL_CALLS(DATATYPE, TRIANGL_CASE, TRIANGL_TYPE)          \
  case TRIANGL_CASE: {                                                        \
    const auto inpTri                                                         \
      = static_cast<TRIANGL_TYPE *>(inputTriangulation.getData());            \
    if(inpTri != nullptr) {                                                   \
      this->interpolateContinuousScalarField<DATATYPE, TRIANGL_TYPE>(         \
        static_cast<DATATYPE *>(ttkUtils::GetVoidPointer(inputScalarField)),  \
        static_cast<DATATYPE *>(ttkUtils::GetVoidPointer(outputScalarField)), \
        *inpTri);                                                             \
    }                                                                         \
    break;                                                                    \
  }

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

int ttkBarycentricSubdivision::RequestData(vtkInformation *ttkNotUsed(request),
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector) {

  ttk::Timer tm;

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  ttk::ExplicitTriangulation triangulationSubdivision{};

  if(triangulation == nullptr) {
    printMsg("Error, internal triangulation is empty.");
    return 0;
  }

  // early return: copy input if no subdivision
  if(SubdivisionLevel == 0) {
    output->ShallowCopy(input);
    return 0;
  }

  this->preconditionTriangulation(triangulation);

  // first iteration: generate the new triangulation
  int ret = this->execute(*triangulation, triangulationSubdivision);
  if(ret != 0) {
    this->printErr("Could not subdivide input mesh");
    return 0;
  }

  // first iteration: interpolate input scalar fields
  ret = InterpolateScalarFields(input, output, *triangulation);
  if(ret != 0) {
    this->printErr("Error interpolating input data array(s)");
    return 0;
  }

  for(unsigned int i = 1; i < SubdivisionLevel; ++i) {
    // move previous points to temp vector
    decltype(points_) tmpPoints{};
    std::swap(points_, tmpPoints);

    // move previous triangulation cells to temp vector
    decltype(cells_connectivity_) tmpCellsCo{};
    std::swap(cells_connectivity_, tmpCellsCo);
    decltype(cells_offsets_) tmpCellsOff{};
    std::swap(cells_offsets_, tmpCellsOff);

    // move previous triangulation to temp triangulation
    decltype(triangulationSubdivision) tmpTr{};
    std::swap(triangulationSubdivision, tmpTr);

#ifdef TTK_CELL_ARRAY_NEW
    tmpTr.setInputCells(
      tmpCellsOff.size() - 1, tmpCellsCo.data(), tmpCellsOff.data());
#else
    ttk::LongSimplexId *tmpCells = nullptr;
    ttk::CellArray::TranslateToFlatLayout(tmpCellsCo, tmpCellsOff, tmpCells);
    tmpTr.setInputCells(tmpCellsCo.size() - 1, tmpCells);
#endif
    tmpTr.setInputPoints(tmpPoints.size() / 3, tmpPoints.data());
    this->preconditionTriangulation(&tmpTr);

    // generate the new triangulation
    this->execute(*triangulation, triangulationSubdivision);

    // interpolate output scalar fields
    InterpolateScalarFields(output, output, *triangulation);
  }

  // generated 3D coordinates
  auto points = vtkSmartPointer<vtkPoints>::New();
  for(size_t i = 0; i < points_.size() / 3; i++) {
    points->InsertNextPoint(&points_[3 * i]);
  }
  output->SetPoints(points);

  // generated triangles
  auto cells = vtkSmartPointer<vtkCellArray>::New();
  for(size_t i = 0; i < cells_offsets_.size() - 1; i++) {
    cells->InsertNextCell(3, &cells_connectivity_[cells_offsets_[i]]);
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

  // shallow copy input field data
  output->GetFieldData()->ShallowCopy(input->GetFieldData());

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 1;
}
