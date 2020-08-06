#include <string>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <ttkDistanceField.h>
#include <vtkPointSet.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkDistanceField);

ttkDistanceField::ttkDistanceField() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkDistanceField::~ttkDistanceField() {
}

int ttkDistanceField::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  if(port == 1)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

int ttkDistanceField::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkDistanceField::RequestData(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  ttk::Timer globalTimer;

  vtkDataSet *domain = vtkDataSet::GetData(inputVector[0]);
  vtkPointSet *sources = vtkPointSet::GetData(inputVector[1]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(domain);
  int ret = this->preconditionTriangulation(triangulation);
  Modified();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    this->printErr("wrong triangulation.");
    return -1;
  }
#endif

  vtkDataArray *identifiers = this->GetOptionalArray(
    ForceInputVertexScalarField, 2, ttk::VertexScalarFieldName, inputVector, 1);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!identifiers) {
    printErr("wrong identifiers.");
    return -2;
  }
#endif

  const SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfPointsInDomain) {
    printErr("domain has no points.");
    return -3;
  }
#endif

  const SimplexId numberOfPointsInSources = sources->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfPointsInSources) {
    printErr("sources have no points.");
    return -4;
  }
#endif

  vtkSmartPointer<ttkSimplexIdTypeArray> origin
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();

  if(origin) {
    origin->SetNumberOfComponents(1);
    origin->SetNumberOfTuples(numberOfPointsInDomain);
    origin->SetName(ttk::VertexScalarFieldName);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  else {
    printErr("ttkSimplexIdTypeArray allocation problem.");
    return -5;
  }
#endif

  vtkSmartPointer<ttkSimplexIdTypeArray> seg
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  if(seg) {
    seg->SetNumberOfComponents(1);
    seg->SetNumberOfTuples(numberOfPointsInDomain);
    seg->SetName("SeedIdentifier");
  }
#ifndef TTK_ENABLE_KAMIKAZE
  else {
    printErr("ttkSimplexIdTypeArray allocation problem.");
    return -6;
  }
#endif

  this->setVertexNumber(numberOfPointsInDomain);
  this->setSourceNumber(numberOfPointsInSources);
  this->setVertexIdentifierScalarFieldPointer(
    ttkUtils::GetVoidPointer(identifiers));
  this->setOutputIdentifiers(ttkUtils::GetVoidPointer(origin));
  this->setOutputSegmentation(ttkUtils::GetVoidPointer(seg));

  vtkDataArray *distanceScalars{};
  switch(OutputScalarFieldType) {
    case DistanceType::Float:
      distanceScalars = vtkFloatArray::New();
      if(distanceScalars) {
        distanceScalars->SetNumberOfComponents(1);
        distanceScalars->SetNumberOfTuples(numberOfPointsInDomain);
        distanceScalars->SetName(OutputScalarFieldName.data());
      }
#ifndef TTK_ENABLE_KAMIKAZE
      else {
        printErr("vtkFloatArray allocation problem.");
        return -7;
      }
#endif
      this->setOutputScalarFieldPointer(
        ttkUtils::GetVoidPointer(distanceScalars));
      ttkTemplateMacro(
        triangulation->getType(), (ret = this->execute<float, TTK_TT>(
                                     (TTK_TT *)triangulation->getData())));
      break;

    case DistanceType::Double:
      distanceScalars = vtkDoubleArray::New();
      if(distanceScalars) {
        distanceScalars->SetNumberOfComponents(1);
        distanceScalars->SetNumberOfTuples(numberOfPointsInDomain);
        distanceScalars->SetName(OutputScalarFieldName.data());
      }
#ifndef TTK_ENABLE_KAMIKAZE
      else {
        printErr("vtkDoubleArray allocation problem.");
        return -8;
      }
#endif
      this->setOutputScalarFieldPointer(
        ttkUtils::GetVoidPointer(distanceScalars));

      ttkTemplateMacro(
        triangulation->getType(), (ret = this->execute<double, TTK_TT>(
                                     (TTK_TT *)triangulation->getData())));
      break;

    default:
#ifndef TTK_ENABLE_KAMIKAZE
      printErr("Scalar field type problem.");
      return -9;
#endif
      break;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(ret) {
    printErr("DistanceField.execute() error code : " + std::to_string(ret));
    return -10;
  }
#endif

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(distanceScalars);
  output->GetPointData()->AddArray(origin);
  output->GetPointData()->AddArray(seg);
  distanceScalars->Delete();

  {
    this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
    this->printMsg(
      "Complete", 1, globalTimer.getElapsedTime()); // global progress, time);
    this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
  }

  return 1;
}
