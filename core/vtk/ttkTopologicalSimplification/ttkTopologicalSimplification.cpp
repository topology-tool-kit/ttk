#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkTopologicalSimplification.h>
#include <ttkUtils.h>

#include <LocalizedTopologicalSimplification.h>

vtkStandardNewMacro(ttkTopologicalSimplification);

ttkTopologicalSimplification::ttkTopologicalSimplification() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int ttkTopologicalSimplification::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkTopologicalSimplification::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkTopologicalSimplification::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  using ttk::SimplexId;

  const auto domain = vtkDataSet::GetData(inputVector[0]);
  const auto constraints = vtkPointSet::GetData(inputVector[1]);
  if(!domain || !constraints)
    return !this->printErr("Unable to retrieve required input data objects.");

  auto output = vtkDataSet::GetData(outputVector);

  // triangulation
  auto triangulation = ttkAlgorithm::GetTriangulation(domain);

  if(!triangulation) {
    this->printErr("Input triangulation pointer is NULL.");
    return -1;
  }

  this->preconditionTriangulation(triangulation);

  if(triangulation->isEmpty()) {
    this->printErr("Triangulation allocation problem.");
    return -1;
  }

  if(!domain) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }

  const auto numberOfVertices = domain->GetNumberOfPoints();
  if(numberOfVertices <= 0) {
    this->printErr("Domain has no points.");
    return -5;
  }

  // domain scalar field
  const auto inputScalars = this->GetInputArrayToProcess(0, domain);
  if(!inputScalars) {
    this->printErr("Input scalar field pointer is null.");
    return -3;
  }

  // domain offset field
  const auto inputOrder
    = this->GetOrderArray(domain, 0, 2, ForceInputOffsetScalarField);
  if(!inputOrder) {
    this->printErr("Wrong input offset scalar field.");
    return -1;
  }

  // create output arrays
  auto outputScalars
    = vtkSmartPointer<vtkDataArray>::Take(inputScalars->NewInstance());
  outputScalars->DeepCopy(inputScalars);
  auto outputOrder
    = vtkSmartPointer<vtkDataArray>::Take(inputOrder->NewInstance());
  outputOrder->DeepCopy(inputOrder);

  // constraint identifier field
  int numberOfConstraints = constraints->GetNumberOfPoints();

  std::vector<ttk::SimplexId> idSpareStorage{};
  auto identifiers = this->GetIdentifierArrayPtr(ForceInputVertexScalarField, 1,
                                                 ttk::VertexScalarFieldName,
                                                 constraints, idSpareStorage);

  // NOTE it'd be better if the two backends were inheriting from the same API
  // (the switch would then happen in the base code)
  int ret{};
  if(this->UseLTS) {
    ttk::lts::LocalizedTopologicalSimplification lts{};
    lts.setDebugLevel(this->debugLevel_);
    lts.setThreadNumber(this->threadNumber_);

    lts.preconditionTriangulation(triangulation);

    ttkVtkTemplateMacro(
      inputScalars->GetDataType(), triangulation->getType(),
      (ret = lts.removeUnauthorizedExtrema<VTK_TT, ttk::SimplexId, TTK_TT>(
         ttkUtils::GetPointer<VTK_TT>(outputScalars),
         ttkUtils::GetPointer<SimplexId>(outputOrder),

         static_cast<TTK_TT *>(triangulation->getData()), identifiers,
         numberOfConstraints, this->AddPerturbation)));

    // TODO: fix convention in original ttk module
    ret = !ret;
  } else {
    switch(inputScalars->GetDataType()) {
      vtkTemplateMacro(
        ret = this->execute(ttkUtils::GetPointer<VTK_TT>(inputScalars),
                            ttkUtils::GetPointer<VTK_TT>(outputScalars),
                            identifiers,
                            ttkUtils::GetPointer<SimplexId>(inputOrder),
                            ttkUtils::GetPointer<SimplexId>(outputOrder),
                            numberOfConstraints, *triangulation->getData()));
    }
  }

  // something wrong in baseCode
  if(ret) {
    this->printErr("TopologicalSimplification.execute() error code: "
                   + std::to_string(ret));
    return -12;
  }

  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(outputOrder);
  output->GetPointData()->AddArray(outputScalars);

  return 1;
}
