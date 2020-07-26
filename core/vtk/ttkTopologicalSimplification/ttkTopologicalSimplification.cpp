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

  triangulation_ = ttkAlgorithm::GetTriangulation(domain);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_) {
    this->printErr("Input triangulation pointer is NULL.");
    return -1;
  }
#endif

  this->preconditionTriangulation(triangulation_);
  Modified();
  hasUpdatedMesh_ = true;

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_->isEmpty()) {
    this->printErr("Triangulation allocation problem.");
    return -1;
  }
#endif

#ifndef TTK_ENABLE_KAMIKAZE
  if(!domain) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }

  if(!domain->GetNumberOfPoints()) {
    this->printErr("Input has no point.");
    return -1;
  }
#endif

  const auto pointData = domain->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    this->printErr("Input has no point data.");
    return -1;
  }
#endif

  if(ScalarField.length()) {
    inputScalars_ = pointData->GetArray(ScalarField.data());
  } else {
    inputScalars_ = pointData->GetArray(ScalarFieldId);
    if(inputScalars_)
      ScalarField = inputScalars_->GetName();
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars_) {
    this->printErr("Input scalar field pointer is null.");
    return -3;
  }
#endif

  if(ForceInputVertexScalarField && InputVertexScalarFieldName.length())
    identifiers_ = constraints->GetPointData()->GetArray(
      InputVertexScalarFieldName.data());
  else if(constraints->GetPointData()->GetArray(ttk::VertexScalarFieldName))
    identifiers_
      = constraints->GetPointData()->GetArray(ttk::VertexScalarFieldName);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!identifiers_) {
    this->printErr("Wrong vertex identifier scalar field.");
    return -1;
  }
#endif

  if(ForceInputOffsetScalarField && InputOffsetScalarFieldName.length()) {
    inputOffsets_
      = domain->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  } else if(OffsetFieldId != -1
            && domain->GetPointData()->GetArray(OffsetFieldId)) {
    inputOffsets_ = domain->GetPointData()->GetArray(OffsetFieldId);
  } else if(domain->GetPointData()->GetArray(ttk::OffsetScalarFieldName)) {
    inputOffsets_
      = domain->GetPointData()->GetArray(ttk::OffsetScalarFieldName);
  } else {
    if(hasUpdatedMesh_ && offsets_) {
      offsets_->Delete();
      offsets_ = nullptr;
      hasUpdatedMesh_ = false;
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

    inputOffsets_ = offsets_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOffsets_) {
    this->printErr("Wrong input offset scalar field.");
    return -1;
  }
#endif

#ifndef TTK_ENABLE_KAMIKAZE
  if(inputOffsets_->GetDataType() != VTK_INT
     and inputOffsets_->GetDataType() != VTK_ID_TYPE) {
    this->printErr("Input offset field type not supported.");
    return -1;
  }
#endif

  const auto numberOfVertices = domain->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfVertices <= 0) {
    this->printErr("Domain has no points.");
    return -5;
  }
#endif

  if(OutputOffsetScalarFieldName.length() <= 0)
    OutputOffsetScalarFieldName = ttk::OffsetScalarFieldName;

  vtkNew<ttkSimplexIdTypeArray> outputOffsets{};
  if(outputOffsets) {
    outputOffsets->SetNumberOfComponents(1);
    outputOffsets->SetNumberOfTuples(numberOfVertices);
    outputOffsets->SetName(OutputOffsetScalarFieldName.data());
  }
#ifndef TTK_ENABLE_KAMIKAZE
  else {
    this->printErr("ttkSimplexIdTypeArray allocation problem.");
    return -7;
  }
#endif

  vtkSmartPointer<vtkDataArray> outputScalars{inputScalars_->NewInstance()};
  if(outputScalars) {
    outputScalars->SetNumberOfTuples(numberOfVertices);
    outputScalars->SetName(inputScalars_->GetName());
  }
#ifndef TTK_ENABLE_KAMIKAZE
  else {
    this->printErr("vtkDataArray allocation problem.");
    return -9;
  }
#endif

  const auto numberOfConstraints = constraints->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfConstraints <= 0) {
    this->printErr("Input has no constraints.");
    return -10;
  }
#endif

#ifndef TTK_ENABLE_KAMIKAZE
  if(identifiers_->GetDataType() != inputOffsets_->GetDataType()) {
    this->printErr("Type of identifiers and offsets are different.");
    return -11;
  }
#endif

  if(inputOffsets_->GetDataType() == VTK_INT) {
    switch(inputScalars_->GetDataType()) {
      vtkTemplateMacro(
        ret = this->execute(
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars_)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
          static_cast<int *>(ttkUtils::GetVoidPointer(identifiers_)),
          static_cast<int *>(ttkUtils::GetVoidPointer(inputOffsets_)),
          static_cast<SimplexId *>(ttkUtils::GetVoidPointer(outputOffsets)),
          numberOfConstraints, *triangulation_->getData()));
    }
  } else if(inputOffsets_->GetDataType() == VTK_ID_TYPE) {
    switch(inputScalars_->GetDataType()) {
      vtkTemplateMacro(
        ret = this->execute(
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars_)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
          static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(identifiers_)),
          static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(inputOffsets_)),
          static_cast<SimplexId *>(ttkUtils::GetVoidPointer(outputOffsets)),
          numberOfConstraints, *triangulation_->getData()));
    }
  }
#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(ret) {
    this->printErr("TopologicalSimplification.execute() error code: "
                   + std::to_string(ret));
    return -12;
  }
#endif

  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(outputOffsets);
  output->GetPointData()->AddArray(outputScalars);

  return 1;
}
