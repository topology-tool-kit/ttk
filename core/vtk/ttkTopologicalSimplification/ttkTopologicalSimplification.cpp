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

  int ret{};

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

  // domain scalar field

  if(!domain) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }

  if(!domain->GetNumberOfPoints()) {
    this->printErr("Input has no point.");
    return -1;
  }

  const auto pointData = domain->GetPointData();

  if(!pointData) {
    this->printErr("Input has no point data.");
    return -1;
  }

  vtkDataArray *inputScalars{};
  if(ScalarField.length()) {
    inputScalars = pointData->GetArray(ScalarField.data());
  } else {
    inputScalars = pointData->GetArray(ScalarFieldId);
    if(inputScalars)
      ScalarField = inputScalars->GetName();
  }

  if(!inputScalars) {
    this->printErr("Input scalar field pointer is null.");
    return -3;
  }

  // constraint identifier field

  vtkDataArray *identifiers{};
  if(ForceInputVertexScalarField && InputVertexScalarFieldName.length())
    identifiers = constraints->GetPointData()->GetArray(
      InputVertexScalarFieldName.data());
  else if(constraints->GetPointData()->GetArray(ttk::VertexScalarFieldName))
    identifiers
      = constraints->GetPointData()->GetArray(ttk::VertexScalarFieldName);

  if(!identifiers) {
    this->printErr("Wrong vertex identifier scalar field.");
    return -1;
  }

  // domain offset field

  vtkDataArray *inputOffsets{};
  if(ForceInputOffsetScalarField && InputOffsetScalarFieldName.length()) {
    inputOffsets
      = domain->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  } else if(OffsetFieldId != -1
            && domain->GetPointData()->GetArray(OffsetFieldId)) {
    inputOffsets = domain->GetPointData()->GetArray(OffsetFieldId);
  } else if(domain->GetPointData()->GetArray(ttk::OffsetScalarFieldName)) {
    inputOffsets = domain->GetPointData()->GetArray(ttk::OffsetScalarFieldName);
  } else {
    if(offsets_ != nullptr) {
      offsets_->Delete();
      offsets_ = nullptr;
    }

    if(!offsets_) {
      const auto numberOfVertices = domain->GetNumberOfPoints();

      offsets_ = ttkSimplexIdTypeArray::New();
      offsets_->SetNumberOfComponents(1);
      offsets_->SetNumberOfTuples(numberOfVertices);
      offsets_->SetName(ttk::OffsetScalarFieldName);
      for(int i = 0; i < numberOfVertices; ++i)
        offsets_->SetTuple1(i, i);
    }

    inputOffsets = offsets_;
  }

  if(!inputOffsets) {
    this->printErr("Wrong input offset scalar field.");
    return -1;
  }

  if(inputOffsets->GetDataType() != VTK_INT
     and inputOffsets->GetDataType() != VTK_ID_TYPE) {
    this->printErr("Input offset field type not supported.");
    return -1;
  }

  const auto numberOfVertices = domain->GetNumberOfPoints();
  if(numberOfVertices <= 0) {
    this->printErr("Domain has no points.");
    return -5;
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

  if(identifiers->GetDataType() != inputOffsets->GetDataType()) {
    this->printErr("Type of identifiers and offsets are different.");
    return -11;
  }

  if(inputOffsets->GetDataType() == VTK_INT) {
    switch(inputScalars->GetDataType()) {
      vtkTemplateMacro(
        ret = this->execute(
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
          static_cast<int *>(ttkUtils::GetVoidPointer(identifiers)),
          static_cast<int *>(ttkUtils::GetVoidPointer(inputOffsets)),
          static_cast<SimplexId *>(ttkUtils::GetVoidPointer(outputOffsets)),
          numberOfConstraints, *triangulation->getData()));
    }
  } else if(inputOffsets->GetDataType() == VTK_ID_TYPE) {
    switch(inputScalars->GetDataType()) {
      vtkTemplateMacro(
        ret = this->execute(
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
          static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(identifiers)),
          static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(inputOffsets)),
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
