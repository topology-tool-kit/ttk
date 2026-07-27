#include <ttkSurfaceGeometrySmoother.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>

vtkStandardNewMacro(ttkSurfaceGeometrySmoother);

ttkSurfaceGeometrySmoother::ttkSurfaceGeometrySmoother() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int ttkSurfaceGeometrySmoother::FillInputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkSurfaceGeometrySmoother::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkSurfaceGeometrySmoother::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  auto inputPointSet = vtkPointSet::GetData(inputVector[0]);
  auto inputMesh = vtkPointSet::GetData(inputVector[1]);
  auto outputPointSet = vtkPointSet::GetData(outputVector);

  auto triangulationToSmooth = ttkAlgorithm::GetTriangulation(inputPointSet);
  if(triangulationToSmooth == nullptr) {
    return 0;
  }
  this->preconditionTriangulationToSmooth(triangulationToSmooth);

  auto triangulationSurface = ttkAlgorithm::GetTriangulation(inputMesh);
  if(triangulationSurface == nullptr) {
    return 0;
  }
  this->preconditionTriangulationSurface(triangulationSurface);

  std::vector<ttk::SimplexId> idSpareStorage{};
  const auto *vertsId = this->GetIdentifierArrayPtr(
    ForceIdentifiersField, 0, ttk::VertexScalarFieldName, inputPointSet,
    idSpareStorage, 0, false);
  if(vertsId == nullptr) {
    this->printWrn("No vertex scalar field detected on input");
  }

  vtkDataArray *inputMaskField = ttkAlgorithm::GetOptionalArray(
    ForceInputMaskScalarField, 1, ttk::MaskScalarFieldName, inputPointSet);

  // Copy all input points/cells + scalar fields
  outputPointSet->DeepCopy(inputPointSet);

  // calling the smoothing package
  auto inputPoints = inputPointSet->GetPoints();
  auto outputPoints = outputPointSet->GetPoints();

  const auto hasMask{this->UseMaskScalarField && inputMaskField != nullptr};

  if(triangulationSurface->getType() != triangulationToSmooth->getType()) {
    this->printErr("Triangulations should have the same type");
    return 0;
  }

  // here we assume that the two triangulation objects have the same
  // underlying type
  ttkTemplateMacro(
    triangulationToSmooth->getType(),
    this->execute(
      ttkUtils::GetPointer<float>(outputPoints->GetData()),
      ttkUtils::GetPointer<float>(inputPoints->GetData()),
      hasMask ? ttkUtils::GetPointer<char>(inputMaskField) : nullptr, vertsId,
      this->NumberOfIterations,
      *static_cast<TTK_TT *>(triangulationToSmooth->getData()),
      *static_cast<TTK_TT *>(triangulationSurface->getData())));

  return 1;
}
