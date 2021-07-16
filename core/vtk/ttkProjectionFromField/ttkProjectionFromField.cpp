#include <ttkProjectionFromField.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkThreshold.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <array>

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

int ttkProjectionFromField::projectDiagramInsideDomain(
  vtkUnstructuredGrid *const inputDiagram,
  vtkUnstructuredGrid *const outputDiagram) {

  ttk::Timer tm{};

  // use vtkThreshold to remove diagonal (PairIdentifier == -1)
  vtkNew<vtkThreshold> threshold{};
  threshold->SetInputDataObject(0, inputDiagram);
  threshold->SetInputArrayToProcess(0, 0, 0,
                                    vtkDataObject::FIELD_ASSOCIATION_CELLS,
                                    ttk::PersistencePairIdentifierName);
  threshold->ThresholdByUpper(0);
  threshold->Update();

  auto diagonalLess = threshold->GetOutput();
  auto diagonalLessData = diagonalLess->GetPointData();

  const auto critCoordinates = vtkFloatArray::SafeDownCast(
    diagonalLessData->GetAbstractArray(ttk::PersistenceCoordinatesName));

  // set new points from Coordinates array
  vtkNew<vtkFloatArray> coords{};
  coords->DeepCopy(critCoordinates);
  coords->SetName("Points");
  diagonalLess->GetPoints()->SetData(coords);
  diagonalLessData->RemoveArray(ttk::PersistenceCoordinatesName);

  outputDiagram->ShallowCopy(diagonalLess);

  // don't forget to forward the Field Data
  outputDiagram->GetFieldData()->ShallowCopy(inputDiagram->GetFieldData());

  this->printMsg("Projected Persistence Diagram inside domain", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);

  return 1;
}

template <typename VTK_TT>
int ttkProjectionFromField::projectDiagramIn2D(
  vtkUnstructuredGrid *const inputDiagram,
  vtkUnstructuredGrid *const outputDiagram,
  const VTK_TT *const births,
  const VTK_TT *const deaths) {

  ttk::Timer tm{};

  outputDiagram->ShallowCopy(inputDiagram);

  auto pointData = outputDiagram->GetPointData();

  vtkNew<vtkFloatArray> coords{};
  coords->DeepCopy(inputDiagram->GetPoints()->GetData());
  coords->SetName(ttk::PersistenceCoordinatesName);
  pointData->AddArray(coords);

  vtkNew<vtkPoints> points{};
  const auto nPoints = inputDiagram->GetNumberOfPoints();
  points->SetNumberOfPoints(nPoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < nPoints / 2; ++i) {
    std::array<float, 3> pt0{
      static_cast<float>(births[i]),
      static_cast<float>(births[i]),
      0,
    };
    std::array<float, 3> pt1{
      static_cast<float>(births[i]),
      static_cast<float>(deaths[i]),
      0,
    };
    points->SetPoint(2 * i + 0, pt0.data());
    points->SetPoint(2 * i + 1, pt1.data());
  }

  outputDiagram->SetPoints(points);

  // add diagonal(first point -> last birth/penultimate point)
  std::array<vtkIdType, 2> diag{0, 2 * (outputDiagram->GetNumberOfCells() - 1)};
  outputDiagram->InsertNextCell(VTK_LINE, 2, diag.data());

  // add diagonal data
  auto cellData = outputDiagram->GetCellData();
  auto pairIdentifierScalars = vtkIntArray::SafeDownCast(
    cellData->GetArray(ttk::PersistencePairIdentifierName));
  auto extremumIndexScalars = vtkIntArray::SafeDownCast(
    cellData->GetArray(ttk::PersistencePairTypeName));
  auto persistenceScalars = cellData->GetArray(ttk::PersistenceName);
  auto birthScalars = cellData->GetArray(ttk::PersistenceBirthName);
  auto deathScalars = cellData->GetArray(ttk::PersistenceDeathName);
  auto isFinite = cellData->GetArray(ttk::PersistenceIsFinite);

  pairIdentifierScalars->InsertNextTuple1(-1);
  extremumIndexScalars->InsertNextTuple1(-1);
  isFinite->InsertNextTuple1(0);
  // 2 * persistence of min-max pair
  persistenceScalars->InsertNextTuple1(2 * persistenceScalars->GetTuple1(0));
  // birth == death == 0
  birthScalars->InsertNextTuple1(0);
  deathScalars->InsertNextTuple1(0);

  // don't forget to forward the Field Data
  outputDiagram->GetFieldData()->ShallowCopy(inputDiagram->GetFieldData());

  this->printMsg("Projected Persistence Diagram back to 2D", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);

  return 1;
}

int ttkProjectionFromField::projectPersistenceDiagram(
  vtkUnstructuredGrid *const inputDiagram,
  vtkUnstructuredGrid *const outputDiagram) {

  auto pointData = inputDiagram->GetPointData();
  auto cellData = inputDiagram->GetCellData();

  // ensure we have the right arrays
  const auto critCoordinates = vtkFloatArray::SafeDownCast(
    pointData->GetAbstractArray(ttk::PersistenceCoordinatesName));
  const auto inputBirths = vtkDataArray::SafeDownCast(
    cellData->GetAbstractArray(ttk::PersistenceBirthName));
  const auto inputDeaths = vtkDataArray::SafeDownCast(
    cellData->GetAbstractArray(ttk::PersistenceDeathName));

  bool embed = critCoordinates == nullptr;

  if(embed) {
    switch(inputBirths->GetDataType()) {
      vtkTemplateMacro(return this->projectDiagramIn2D(
        inputDiagram, outputDiagram,
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputBirths)),
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputDeaths))));
    }
  } else {
    if(critCoordinates == nullptr) {
      this->printErr("Missing `Coordinates' vtkPointData array");
      return 0;
    }
    if(critCoordinates->GetNumberOfComponents() != 3) {
      this->printErr("`Coordinates' array should have 3 components");
      return 0;
    }
    return this->projectDiagramInsideDomain(inputDiagram, outputDiagram);
  }

  return 1;
}

int ttkProjectionFromField::RequestData(vtkInformation *ttkNotUsed(request),
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

  vtkNew<vtkPoints> pointSet{};
  pointSet->SetNumberOfPoints(input->GetNumberOfPoints());

  if(UseTextureCoordinates) {

    const auto textureCoordinates = input->GetPointData()->GetTCoords();
    if(textureCoordinates == nullptr) {
      return 0;
    }
    printMsg("Starting computation with texture coordinates...");

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < input->GetNumberOfPoints(); i++) {
      std::array<double, 3> pt{};
      textureCoordinates->GetTuple(i, pt.data());
      pointSet->SetPoint(i, pt[0], pt[1], pt[2]);
    }

  } else if(this->Use3DCoordinatesArray) {

    const auto inputCoordsArray = this->GetInputArrayToProcess(2, inputVector);
    if(inputCoordsArray == nullptr) {
      return 0;
    }
    printMsg("Starting computation...");
    printMsg(std::vector<std::vector<std::string>>{
      {"  Coordinates Array", inputCoordsArray->GetName()}});

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < input->GetNumberOfPoints(); i++) {
      std::array<double, 3> pt{};
      inputCoordsArray->GetTuple(i, pt.data());
      pointSet->SetPoint(i, pt[0], pt[1], pt[2]);
    }

  } else {

    const auto inputScalarFieldU = this->GetInputArrayToProcess(0, inputVector);
    const auto inputScalarFieldV = this->GetInputArrayToProcess(1, inputVector);

    if(inputScalarFieldU == nullptr || inputScalarFieldV == nullptr) {
      return 0;
    }

    printMsg("Starting computation...");
    printMsg({{"  U-component", inputScalarFieldU->GetName()},
              {"  V-component", inputScalarFieldV->GetName()}});

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < input->GetNumberOfPoints(); i++) {
      pointSet->SetPoint(i, inputScalarFieldU->GetComponent(i, 0),
                         inputScalarFieldV->GetComponent(i, 0), 0);
    }
  }

  output->SetPoints(pointSet);

  printMsg(std::to_string(input->GetNumberOfPoints()) + " points projected", 1,
           t.getElapsedTime(), threadNumber_);

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
