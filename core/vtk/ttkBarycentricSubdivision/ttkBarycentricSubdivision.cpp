#include <ttkBarycentricSubdivision.h>

#define MODULE_S "[ttkBarycentricSubdivision] "
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

vtkStandardNewMacro(ttkBarycentricSubdivision);

int ttkBarycentricSubdivision::doIt(std::vector<vtkDataSet *> &inputs,
                                    std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  auto input = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  auto output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  auto triangulation = ttkTriangulation::getTriangulation(input);
  ttk::Triangulation triangulationSubdivision;

  if(triangulation == nullptr) {
    return -1;
  }

  triangulation->setWrapper(this);
  baseWorker_.setupTriangulation(triangulation);
  baseWorker_.setWrapper(this);
  baseWorker_.setOutputTriangulation(&triangulationSubdivision);
  baseWorker_.setInputPoints(input->GetPoints()->GetVoidPointer(0));

  // generate the new triangulation
  baseWorker_.execute();

  size_t npointdata = input->GetPointData()->GetNumberOfArrays();
  size_t ncelldata = input->GetCellData()->GetNumberOfArrays();

  auto allocateScalarField = [&](vtkDataArray *const inputScalarField) {
    vtkSmartPointer<vtkDataArray> outputScalarField{};

    // allocate the memory for the output scalar field
    switch(inputScalarField->GetDataType()) {
      case VTK_CHAR:
        outputScalarField = vtkSmartPointer<vtkCharArray>::New();
        break;
      case VTK_DOUBLE:
        outputScalarField = vtkSmartPointer<vtkDoubleArray>::New();
        break;
      case VTK_FLOAT:
        outputScalarField = vtkSmartPointer<vtkFloatArray>::New();
        break;
      case VTK_INT:
        outputScalarField = vtkSmartPointer<vtkIntArray>::New();
        break;
      case VTK_ID_TYPE:
        outputScalarField = vtkSmartPointer<vtkIdTypeArray>::New();
        break;
      default:
        break;
    }
    return outputScalarField;
  };

  const auto outPointsNumber = baseWorker_.getNumberOfVertices();

  for(size_t i = 0; i < npointdata; ++i) {
    auto inputScalarField = input->GetPointData()->GetArray(i);
    if(inputScalarField == nullptr) {
      return -2;
    }

#define DISPATCH_INTERPOLATE_DIS(CASE, TYPE)                      \
  case CASE:                                                      \
    baseWorker_.interpolateDiscreteScalarField<TYPE>(             \
      static_cast<TYPE *>(inputScalarField->GetVoidPointer(0)),   \
      static_cast<TYPE *>(outputScalarField->GetVoidPointer(0))); \
    break
#define DISPATCH_INTERPOLATE_CONT(CASE, TYPE)                     \
  case CASE:                                                      \
    baseWorker_.interpolateContinuousScalarField<TYPE>(           \
      static_cast<TYPE *>(inputScalarField->GetVoidPointer(0)),   \
      static_cast<TYPE *>(outputScalarField->GetVoidPointer(0))); \
    break

    auto outputScalarField = allocateScalarField(inputScalarField);
    // only for scalar fields
    outputScalarField->SetNumberOfComponents(1);
    outputScalarField->SetNumberOfTuples(outPointsNumber);
    outputScalarField->SetName(inputScalarField->GetName());
    output->GetPointData()->AddArray(outputScalarField);
    switch(inputScalarField->GetDataType()) {
      DISPATCH_INTERPOLATE_DIS(VTK_CHAR, char);
      DISPATCH_INTERPOLATE_DIS(VTK_INT, int);
      DISPATCH_INTERPOLATE_DIS(VTK_ID_TYPE, vtkIdType);
      DISPATCH_INTERPOLATE_CONT(VTK_FLOAT, float);
      DISPATCH_INTERPOLATE_CONT(VTK_DOUBLE, double);
    }
    output->GetPointData()->AddArray(outputScalarField);
  }

  for(size_t i = 0; i < ncelldata; ++i) {
    auto inputScalarField = input->GetCellData()->GetArray(i);
    if(inputScalarField == nullptr) {
      return -2;
    }

    auto outputScalarField = allocateScalarField(inputScalarField);
    // only for scalar fields
    outputScalarField->SetNumberOfComponents(1);
    outputScalarField->SetNumberOfTuples(outPointsNumber);
    outputScalarField->SetName(inputScalarField->GetName());
    output->GetPointData()->AddArray(outputScalarField);
    switch(inputScalarField->GetDataType()) {
      ttkTemplateMacro(baseWorker_.interpolateCellDataField<VTK_TT>(
        static_cast<VTK_TT *>(inputScalarField->GetVoidPointer(0)),
        static_cast<VTK_TT *>(outputScalarField->GetVoidPointer(0))));
    }
    output->GetCellData()->AddArray(outputScalarField);
  }

  // output variables
  auto &outPoints = baseWorker_.points_;
  auto &outCells = baseWorker_.cells_;
  auto &outPointId = baseWorker_.pointId_;
  auto &outPointDim = baseWorker_.pointDim_;

  // generated 3D coordinates
  auto points = vtkSmartPointer<vtkPoints>::New();
  for(size_t i = 0; i < outPoints.size() / 3; i++) {
    points->InsertNextPoint(&outPoints[3 * i]);
  }
  output->SetPoints(points);

  // generated triangles
  const size_t dataPerCell = 4;
  auto cells = vtkSmartPointer<vtkCellArray>::New();
  for(size_t i = 0; i < outCells.size() / dataPerCell; i++) {
    cells->InsertNextCell(3, &outCells[dataPerCell * i + 1]);
  }
  output->SetCells(VTK_TRIANGLE, cells);

  // cell id
  auto cellId = vtkSmartPointer<vtkIntArray>::New();
  cellId->SetName("CellId");
  cellId->SetVoidArray(outPointId.data(), outPointId.size(), 1);
  output->GetPointData()->AddArray(cellId);

  // cell dimension
  auto cellDim = vtkSmartPointer<vtkIntArray>::New();
  cellDim->SetName("CellDimension");
  cellDim->SetVoidArray(outPointDim.data(), outPointDim.size(), 1);
  output->GetPointData()->AddArray(cellDim);

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(std::cout, msg.str(), memoryMsg);
  }

  return 0;
}
