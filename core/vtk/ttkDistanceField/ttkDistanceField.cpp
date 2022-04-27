#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>

#include <ttkDistanceField.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkDistanceField);

ttkDistanceField::ttkDistanceField() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int ttkDistanceField::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkDistanceField::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkDistanceField::RequestData(vtkInformation *ttkNotUsed(request),
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  ttk::Timer globalTimer;

  vtkDataSet *domain = vtkDataSet::GetData(inputVector[0]);
  vtkPointSet *sources = vtkPointSet::GetData(inputVector[1]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  auto triangulation = ttkAlgorithm::GetTriangulation(domain);
  if(triangulation == nullptr) {
    this->printErr("Wrong triangulation.");
    return 0;
  }

  this->preconditionTriangulation(triangulation);

  std::vector<ttk::SimplexId> idSpareStorage{};
  auto *identifiers = this->GetIdentifierArrayPtr(ForceInputVertexScalarField,
                                                  0, ttk::VertexScalarFieldName,
                                                  sources, idSpareStorage);
  if(identifiers == nullptr) {
    printErr("Wrong identifiers.");
    return 0;
  }

  const int numberOfPointsInDomain = domain->GetNumberOfPoints();
  if(numberOfPointsInDomain == 0) {
    printErr("Domain has no points.");
    return 0;
  }

  const int numberOfPointsInSources = sources->GetNumberOfPoints();
  if(numberOfPointsInSources == 0) {
    printErr("Sources have no points.");
    return 0;
  }

  vtkNew<ttkSimplexIdTypeArray> origin{};
  origin->SetNumberOfComponents(1);
  origin->SetNumberOfTuples(numberOfPointsInDomain);
  origin->SetName(ttk::VertexScalarFieldName);

  vtkNew<ttkSimplexIdTypeArray> seg{};
  seg->SetNumberOfComponents(1);
  seg->SetNumberOfTuples(numberOfPointsInDomain);
  seg->SetName("SeedIdentifier");

  this->setVertexNumber(numberOfPointsInDomain);
  this->setSourceNumber(numberOfPointsInSources);
  this->setVertexIdentifierScalarFieldPointer(identifiers);
  this->setOutputIdentifiers(ttkUtils::GetVoidPointer(origin));
  this->setOutputSegmentation(ttkUtils::GetVoidPointer(seg));

  vtkSmartPointer<vtkDataArray> distanceScalars{};
  int ret{};
  switch(static_cast<DistanceType>(OutputScalarFieldType)) {
    case DistanceType::Float:
      distanceScalars = vtkFloatArray::New();
      distanceScalars->SetNumberOfComponents(1);
      distanceScalars->SetNumberOfTuples(numberOfPointsInDomain);
      distanceScalars->SetName(OutputScalarFieldName.data());

      this->setOutputScalarFieldPointer(
        ttkUtils::GetVoidPointer(distanceScalars));
      ttkTemplateMacro(
        triangulation->getType(), (ret = this->execute<float, TTK_TT>(
                                     (TTK_TT *)triangulation->getData())));
      break;

    case DistanceType::Double:
      distanceScalars = vtkDoubleArray::New();
      distanceScalars->SetNumberOfComponents(1);
      distanceScalars->SetNumberOfTuples(numberOfPointsInDomain);
      distanceScalars->SetName(OutputScalarFieldName.data());

      this->setOutputScalarFieldPointer(
        ttkUtils::GetVoidPointer(distanceScalars));

      ttkTemplateMacro(
        triangulation->getType(), (ret = this->execute<double, TTK_TT>(
                                     (TTK_TT *)triangulation->getData())));
      break;

    default:
      printErr("Invalid scalar field type.");
      return 0;
      break;
  }

  // something wrong in baseCode
  if(ret != 0) {
    printErr("DistanceField.execute() error code : " + std::to_string(ret));
    return 0;
  }

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(distanceScalars);
  output->GetPointData()->AddArray(origin);
  output->GetPointData()->AddArray(seg);
  distanceScalars->Delete();

  return 1;
}
