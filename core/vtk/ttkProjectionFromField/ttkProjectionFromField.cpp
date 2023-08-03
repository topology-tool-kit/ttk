#include <ttkPersistenceDiagramUtils.h>
#include <ttkProjectionFromField.h>

#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

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

int ttkProjectionFromField::projectPersistenceDiagram(
  vtkUnstructuredGrid *const inputDiagram,
  vtkUnstructuredGrid *const outputDiagram) {

  auto pointData = inputDiagram->GetPointData();

  // ensure we have the right arrays
  const auto critCoordinates = vtkFloatArray::SafeDownCast(
    pointData->GetAbstractArray(ttk::PersistenceCoordinatesName));
  bool const embed = critCoordinates == nullptr;
  int ret{0};

  if(embed) {
    ret = ProjectDiagramIn2D(inputDiagram, outputDiagram, *this);
  } else {
    if(critCoordinates == nullptr) {
      this->printErr("Missing `Coordinates' vtkPointData array");
      return 0;
    }
    if(critCoordinates->GetNumberOfComponents() != 3) {
      this->printErr("`Coordinates' array should have 3 components");
      return 0;
    }
    ret = ProjectDiagramInsideDomain(inputDiagram, outputDiagram, *this);
  }

  return ret == 0 ? 1 : 0;
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
