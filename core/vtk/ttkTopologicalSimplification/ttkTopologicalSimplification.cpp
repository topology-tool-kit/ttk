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
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  using ttk::SimplexId;

  const auto domain = vtkDataSet::GetData(inputVector[0]);
  const auto constraints = vtkPointSet::GetData(inputVector[1]);
  auto output = vtkDataSet::GetData(outputVector);

  // triangulation
  auto triangulation = ttkAlgorithm::GetTriangulation(domain);

  if(!triangulation) {
    this->printErr("Input triangulation pointer is NULL.");
    return -1;
  }

  this->preconditionTriangulation(triangulation);
  Modified();

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

  // constraint identifier field
  std::vector<SimplexId> idSpareStorage{};
  const auto identifiers = this->GetIdentifierArrayPtr(
    ForceInputVertexScalarField, 1, ttk::VertexScalarFieldName, constraints,
    idSpareStorage);

  if(!identifiers) {
    this->printErr("Wrong vertex identifier scalar field.");
    return -1;
  }

  // domain offset field
  const auto offsets
    = this->GetOrderArray(domain, 0, 2, ForceInputOffsetScalarField);
  if(!offsets) {
    this->printErr("Wrong input offset scalar field.");
    return -1;
  }

  if(offsets->GetDataType() != VTK_INT
     and offsets->GetDataType() != VTK_ID_TYPE) {
    this->printErr("Input offset field type not supported.");
    return -1;
  }

  vtkNew<ttkSimplexIdTypeArray> outputOffsets{};
  if(outputOffsets) {
    outputOffsets->SetNumberOfComponents(1);
    outputOffsets->SetNumberOfTuples(numberOfVertices);
    outputOffsets->SetName(offsets->GetName());
  } else {
    this->printErr("ttkSimplexIdTypeArray allocation problem.");
    return -7;
  }

  vtkSmartPointer<vtkDataArray> outputScalars{inputScalars->NewInstance()};
  if(outputScalars) {
    outputScalars->SetNumberOfTuples(numberOfVertices);
    outputScalars->SetName(inputScalars->GetName());
  } else {
    this->printErr("vtkDataArray allocation problem.");
    return -9;
  }

  const auto numberOfConstraints = constraints->GetNumberOfPoints();
  if(numberOfConstraints <= 0) {
    this->printErr("Input has no constraints.");
    return -10;
  }

  int ret{};
  if(this->UseLTS) {
    ttk::LocalizedTopologicalSimplification lts{};
    lts.setDebugLevel(this->debugLevel_);
    lts.setThreadNumber(this->threadNumber_);

    lts.preconditionTriangulation(triangulation);

    //   switch(inputScalars->GetDataType()) {
    ttkVtkTemplateMacro(
      inputScalars->GetDataType(), triangulation->getType(),
      (ret = lts.removeUnauthorizedExtrema<VTK_TT, ttk::SimplexId, TTK_TT>(
         static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
         static_cast<SimplexId *>(ttkUtils::GetVoidPointer(outputOffsets)),

         static_cast<TTK_TT *>(triangulation->getData()),
         static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars)),
         static_cast<SimplexId *>(ttkUtils::GetVoidPointer(offsets)),
         identifiers, numberOfConstraints, this->AddPerturbation)));

    // TODO: fix convention in original ttk module
    ret = !ret;
  } else {
    switch(inputScalars->GetDataType()) {
      vtkTemplateMacro(
        ret = this->execute(
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
          identifiers,
          static_cast<SimplexId *>(ttkUtils::GetVoidPointer(offsets)),
          static_cast<SimplexId *>(ttkUtils::GetVoidPointer(outputOffsets)),
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
  output->GetPointData()->AddArray(outputOffsets);
  output->GetPointData()->AddArray(outputScalars);

  return 1;
}
