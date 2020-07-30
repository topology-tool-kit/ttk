#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>

#include <ttkMacros.h>
#include <ttkTopologicalSimplification.h>
#include <ttkUtils.h>

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

  // use the GetOptionalArray variant here to fix a segfault occuring
  // when changing a threshold bound on the Threshold filter usually
  // connected to the second input port
  const auto identifiers = this->GetOptionalArray(
    ForceInputVertexScalarField, 1, ttk::VertexScalarFieldName, constraints);

  if(!identifiers) {
    this->printErr("Wrong vertex identifier scalar field.");
    return -1;
  }

  // domain offset field
  const auto inputOffsets = this->GetOptionalArray(
    ForceInputOffsetScalarField, 2, ttk::OffsetScalarFieldName, inputVector);
  vtkNew<ttkSimplexIdTypeArray> offsets{};
  if(inputOffsets == nullptr) {
    // fill in with default offset
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(numberOfVertices);
    offsets->SetName(ttk::OffsetScalarFieldName);
    for(int i = 0; i < numberOfVertices; ++i) {
      offsets->SetTuple1(i, i);
    }
  } else {
    offsets->ShallowCopy(inputOffsets);
  }

  if(!offsets) {
    this->printErr("Wrong input offset scalar field.");
    return -1;
  }

  if(offsets->GetDataType() != VTK_INT
     and offsets->GetDataType() != VTK_ID_TYPE) {
    this->printErr("Input offset field type not supported.");
    return -1;
  }

  if(OutputOffsetScalarFieldName.length() <= 0)
    OutputOffsetScalarFieldName = ttk::OffsetScalarFieldName;

  vtkNew<ttkSimplexIdTypeArray> outputOffsets{};
  if(outputOffsets) {
    outputOffsets->SetNumberOfComponents(1);
    outputOffsets->SetNumberOfTuples(numberOfVertices);
    outputOffsets->SetName(OutputOffsetScalarFieldName.data());
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

  if(identifiers->GetDataType() != offsets->GetDataType()) {
    this->printErr("Type of identifiers and offsets are different.");
    return -11;
  }

  int ret{};
  if(offsets->GetDataType() == VTK_INT) {
    switch(inputScalars->GetDataType()) {
      vtkTemplateMacro(
        ret = this->execute(
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
          static_cast<int *>(ttkUtils::GetVoidPointer(identifiers)),
          static_cast<int *>(ttkUtils::GetVoidPointer(offsets)),
          static_cast<SimplexId *>(ttkUtils::GetVoidPointer(outputOffsets)),
          numberOfConstraints, *triangulation->getData()));
    }
  } else if(offsets->GetDataType() == VTK_ID_TYPE) {
    switch(inputScalars->GetDataType()) {
      vtkTemplateMacro(
        ret = this->execute(
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
          static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(identifiers)),
          static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(offsets)),
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
