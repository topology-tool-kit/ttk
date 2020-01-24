#include <ttkHarmonicField.h>

// VTK includes
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#ifndef TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET) \
  if(COND) {                         \
    this->printErr(MSG);             \
    return RET;                      \
  }
#else // TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)
#endif // TTK_ENABLE_KAMIKAZE

vtkStandardNewMacro(ttkHarmonicField);

ttkHarmonicField::ttkHarmonicField() {
  SetNumberOfInputPorts(2);
  SetNumberOfOutputPorts(1);
}

int ttkHarmonicField::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  } else if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
  }
  return 1;
}

int ttkHarmonicField::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkHarmonicField::getIdentifiers(vtkPointSet *input) {
  auto vsfn = static_cast<const char *>(ttk::VertexScalarFieldName);

  if(ForceConstraintIdentifiers && InputIdentifiersFieldName.length() != 0) {
    identifiers_
      = input->GetPointData()->GetArray(InputIdentifiersFieldName.data());
  } else if(input->GetPointData()->GetArray(vsfn) != nullptr) {
    identifiers_ = input->GetPointData()->GetArray(vsfn);
  }

  TTK_ABORT_KK(
    identifiers_ == nullptr, "wrong vertex identifier scalar field", -1);
  return 0;
}

int ttkHarmonicField::getConstraints(vtkPointSet *input) {
  auto osfn = static_cast<const char *>(ttk::OffsetScalarFieldName);
  auto pointData = input->GetPointData();

  if(InputScalarFieldName.length() != 0) {
    constraints_ = pointData->GetArray(InputScalarFieldName.data());
  } else if(pointData->GetArray(osfn) != nullptr) {
    constraints_ = pointData->GetArray(osfn);
  }

  TTK_ABORT_KK(constraints_ == nullptr, "wrong constraints scalar field", -1);

  if(constraints_->IsA("vtkDoubleArray")) {
    OutputScalarFieldType = FieldType::DOUBLE;
  } else if(constraints_->IsA("vtkFloatArray")) {
    OutputScalarFieldType = FieldType::FLOAT;
  } else {
    return -2;
  }

  return 0;
}

int ttkHarmonicField::RequestData(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {

  const auto domain = vtkDataSet::GetData(inputVector[0]);
  const auto identifiers = vtkPointSet::GetData(inputVector[1]);
  auto output = vtkDataSet::GetData(outputVector);
  auto triangulation = ttkAlgorithm::GetTriangulation(domain);

  TTK_ABORT_KK(triangulation == nullptr, "wrong triangulation", -1);

  this->preconditionTriangulation(*triangulation, UseCotanWeights);

  int res = this->getIdentifiers(identifiers);

  TTK_ABORT_KK(res != 0, "wrong identifiers", -2);

  res += this->getConstraints(identifiers);

  TTK_ABORT_KK(res != 0, "wrong constraints", -2);

  const auto numberOfPointsInDomain = domain->GetNumberOfPoints();

  TTK_ABORT_KK(numberOfPointsInDomain == 0, "domain has no points", -3);

  const auto numberOfPointsInSources = identifiers->GetNumberOfPoints();

  TTK_ABORT_KK(numberOfPointsInSources == 0, "sources has no points", -3);

  vtkSmartPointer<vtkDataArray> harmonicScalarField{};

  switch(OutputScalarFieldType) {
    case FieldType::FLOAT:
      harmonicScalarField = vtkSmartPointer<vtkFloatArray>::New();
      break;
    case FieldType::DOUBLE:
      harmonicScalarField = vtkSmartPointer<vtkDoubleArray>::New();
      break;
    default:
      this->printErr("Unknown scalar field type");
      return -7;
      break;
  }

  TTK_ABORT_KK(
    harmonicScalarField == nullptr, "vtkArray allocation problem", -8);

  harmonicScalarField->SetNumberOfComponents(1);
  harmonicScalarField->SetNumberOfTuples(numberOfPointsInDomain);
  harmonicScalarField->SetName(OutputScalarFieldName.data());

  switch(OutputScalarFieldType) {
    case FieldType::FLOAT:
      res += this->execute<float>(
        *triangulation, numberOfPointsInSources,
        static_cast<ttk::SimplexId *>(identifiers_->GetVoidPointer(0)),
        static_cast<float *>(constraints_->GetVoidPointer(0)),
        static_cast<float *>(harmonicScalarField->GetVoidPointer(0)),
        UseCotanWeights, SolvingMethod, LogAlpha);
      break;
    case FieldType::DOUBLE:
      res += this->execute<double>(
        *triangulation, numberOfPointsInSources,
        static_cast<ttk::SimplexId *>(identifiers_->GetVoidPointer(0)),
        static_cast<double *>(constraints_->GetVoidPointer(0)),
        static_cast<double *>(harmonicScalarField->GetVoidPointer(0)),
        UseCotanWeights, SolvingMethod, LogAlpha);
      break;
    default:
      break;
  }

  TTK_ABORT_KK(
    res != 0, "HarmonicField execute() error code: " + std::to_string(res), -9);

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(harmonicScalarField);

  return 1;
}
