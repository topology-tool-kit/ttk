#include <ttkProjectionFromField.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkUnstructuredGrid.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkProjectionFromField);

ttkProjectionFromField::ttkProjectionFromField() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  this->setDebugMsgPrefix("ProjectionFromField");
}

int ttkProjectionFromField::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkProjectionFromField::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkProjectionFromField::projectPersistenceDiagram(
  vtkUnstructuredGrid *const inputDiagram,
  vtkUnstructuredGrid *const outputDiagram) {

  ttk::Timer tm{};

  auto pointData = inputDiagram->GetPointData();

  const auto vertexIdentifierScalars = vtkIntArray::SafeDownCast(
    pointData->GetAbstractArray(ttk::VertexScalarFieldName));
  const auto nodeTypeScalars
    = vtkIntArray::SafeDownCast(pointData->GetAbstractArray("CriticalType"));
  const auto critCoordinates
    = vtkFloatArray::SafeDownCast(pointData->GetAbstractArray("Coordinates"));

  auto cellData = inputDiagram->GetCellData();

  const auto pairIdentifierScalars
    = vtkIntArray::SafeDownCast(cellData->GetAbstractArray("PairIdentifier"));
  const auto extremumIndexScalars
    = vtkIntArray::SafeDownCast(cellData->GetAbstractArray("PairType"));
  const auto persistenceScalars
    = vtkDoubleArray::SafeDownCast(cellData->GetAbstractArray("Persistence"));

  // ensure we have a Coordinates array
  if(critCoordinates == nullptr) {
    this->printErr("Missing `Coordinates' vtkPointData array");
    return 0;
  }
  if(critCoordinates->GetNumberOfComponents() != 3) {
    this->printErr("`Coordinates' array should have 3 components");
    return 0;
  }
  if(vertexIdentifierScalars == nullptr || nodeTypeScalars == nullptr
     || critCoordinates == nullptr || pairIdentifierScalars == nullptr
     || extremumIndexScalars == nullptr || persistenceScalars == nullptr) {
    this->printErr("Missing at least one data array");
    return 0;
  }

  // set new points from Coordinates array
  vtkNew<vtkFloatArray> coords{};
  coords->DeepCopy(critCoordinates);
  coords->SetName("Points");
  vtkNew<vtkPoints> points{};
  points->SetData(coords);

  const auto inputPoints = inputDiagram->GetPoints();

  // generate a birth and death arrays from diagram points coordinates
  vtkNew<vtkFloatArray> births{}, deaths{};
  births->SetNumberOfComponents(1);
  births->SetName("Birth");
  births->SetNumberOfTuples(inputPoints->GetNumberOfPoints());
  deaths->SetNumberOfComponents(1);
  deaths->SetName("Death");
  deaths->SetNumberOfTuples(inputPoints->GetNumberOfPoints());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < inputPoints->GetNumberOfPoints(); ++i) {
    births->SetTuple1(i, inputPoints->GetPoint(i)[0]);
    deaths->SetTuple1(i, inputPoints->GetPoint(i)[1]);
  }

  // remove diagonal data from cell arrays
  const auto inputCells = inputDiagram->GetCells();
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->DeepCopy(inputCells->GetOffsetsArray());
  connectivity->DeepCopy(inputCells->GetConnectivityArray());
  // remove last entry (diagonal data)
  offsets->SetNumberOfTuples(offsets->GetNumberOfTuples() - 1);
  connectivity->SetNumberOfTuples(connectivity->GetNumberOfTuples() - 1);
  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);

  // copy cell data arrays, removing diagonal data
  vtkNew<vtkIntArray> pairIds{}, pairTypes{};
  pairIds->DeepCopy(pairIdentifierScalars);
  pairTypes->DeepCopy(extremumIndexScalars);
  pairIds->SetNumberOfTuples(cells->GetNumberOfCells());
  pairTypes->SetNumberOfTuples(cells->GetNumberOfCells());
  vtkNew<vtkDoubleArray> persistence{};
  persistence->DeepCopy(persistenceScalars);
  persistence->SetNumberOfTuples(cells->GetNumberOfCells());

  // create a new vtkUnstructuredGrid
  vtkNew<vtkUnstructuredGrid> persistenceDiagram{};
  persistenceDiagram->SetPoints(points);
  persistenceDiagram->SetCells(VTK_LINE, cells);

  // add  data arrays
  persistenceDiagram->GetPointData()->AddArray(vertexIdentifierScalars);
  persistenceDiagram->GetPointData()->AddArray(nodeTypeScalars);
  persistenceDiagram->GetPointData()->AddArray(births);
  persistenceDiagram->GetPointData()->AddArray(deaths);

  persistenceDiagram->GetCellData()->AddArray(pairIds);
  persistenceDiagram->GetCellData()->AddArray(pairTypes);
  persistenceDiagram->GetCellData()->AddArray(persistence);

  outputDiagram->ShallowCopy(persistenceDiagram);

  // don't forget to forward the Field Data
  outputDiagram->GetFieldData()->ShallowCopy(inputDiagram->GetFieldData());

  this->printMsg("Projected Persistence Diagram inside domain", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);

  return 1;
}

int ttkProjectionFromField::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  ttk::Timer t;

  vtkPointSet *input = vtkPointSet::GetData(inputVector[0]);
  vtkPointSet *output = vtkPointSet::GetData(outputVector, 0);

  if(this->ProjectPersistenceDiagram) {
    auto inputGrid = vtkUnstructuredGrid::SafeDownCast(input);
    auto outputGrid = vtkUnstructuredGrid::SafeDownCast(output);
    if(inputGrid != nullptr && outputGrid != nullptr) {
      return projectPersistenceDiagram(inputGrid, outputGrid);
    }
    this->printErr("Input should be a vtkUnstructuredGrid");
    return 0;
  }

  output->ShallowCopy(input);

  vtkDataArray *inputScalarFieldU = nullptr;
  vtkDataArray *inputScalarFieldV = nullptr;
  vtkDataArray *textureCoordinates = nullptr;

  if(UseTextureCoordinates) {
    textureCoordinates = input->GetPointData()->GetTCoords();

    if(!textureCoordinates)
      return -3;

    printMsg("Starting computation with texture coordinates...");
  } else {

    inputScalarFieldU = this->GetInputArrayToProcess(0, inputVector);
    inputScalarFieldV = this->GetInputArrayToProcess(1, inputVector);

    if(!inputScalarFieldU)
      return -1;

    if(!inputScalarFieldV)
      return -2;

    printMsg("Starting computation...");
    printMsg({{"  U-component", inputScalarFieldU->GetName()},
              {"  V-component", inputScalarFieldV->GetName()}});
  }

  vtkNew<vtkPoints> pointSet{};

  if(pointSet->GetNumberOfPoints() != input->GetNumberOfPoints()) {
    pointSet->SetNumberOfPoints(input->GetNumberOfPoints());
  }

  std::vector<std::vector<double>> points(threadNumber_);

  for(int i = 0; i < threadNumber_; i++) {
    points[i].resize(3);
    points[i][2] = 0;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < input->GetNumberOfPoints(); i++) {

#ifdef TTK_ENABLE_OPENMP
    const int threadId = omp_get_thread_num();
#else
    const int threadId = 0;
#endif

    if(UseTextureCoordinates) {
      textureCoordinates->GetTuple(i, points[threadId].data());
    } else {
      points[threadId][0] = inputScalarFieldU->GetComponent(i, 0);
      points[threadId][1] = inputScalarFieldV->GetComponent(i, 0);
    }

    pointSet->SetPoint(
      i, points[threadId][0], points[threadId][1], points[threadId][2]);
  }

  output->SetPoints(pointSet);

  printMsg(std::to_string(input->GetNumberOfPoints()) + " points projected", 1,
           t.getElapsedTime(), threadNumber_);

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
