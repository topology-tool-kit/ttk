#include <string>
#include <ttkUtils.h>
#include <ttkMacros.h>

#include <ttkDistanceField.h>
#include <vtkPointSet.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkDistanceField)

  ttkDistanceField::ttkDistanceField()
  : identifiers_{} {
  //OutputScalarFieldType = 0;
  //OutputScalarFieldName = "OutputDistanceField";
  //ForceInputVertexScalarField = false;
  //InputVertexScalarFieldName = ttk::VertexScalarFieldName;
  //UseAllCores = true;
  SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);

  triangulation_ = NULL;
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

int ttkDistanceField::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkDistanceField::getTriangulation(vtkDataSet *input) {

  //triangulation_ = ttkTriangulation::getTriangulation(input);
  triangulation_ = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation_) return -1;

  //triangulation_->setWrapper(this);
  //this->setWrapper(this);

  this->setupTriangulation(triangulation_);
  Modified();

//#ifndef TTK_ENABLE_KAMIKAZE
  //// allocation problem
  //if(triangulation_->isEmpty()) {
    //cerr << "[ttkDistanceField] Error : ttkTriangulation allocation problem."
         //<< endl;
    //return -1;
  //}
//#endif

  return 0;
}

int ttkDistanceField::getIdentifiers(vtkDataSet *input) {
  if(ForceInputVertexScalarField and InputVertexScalarFieldName.length())
    identifiers_
      = input->GetPointData()->GetArray(InputVertexScalarFieldName.data());
  else if(input->GetPointData()->GetArray(ttk::VertexScalarFieldName))
    identifiers_ = input->GetPointData()->GetArray(ttk::VertexScalarFieldName);

#ifndef TTK_ENABLE_KAMIKAZE
  // allocation problem
  if(!identifiers_) {
    printErr("[ttkDistanceField] Error : wrong vertex identifiers.");
    return -1;
  }
#endif

  return 0;
}

int ttkDistanceField::RequestData(vtkInformation *request,vtkInformationVector **inputVector,vtkInformationVector *outputVector) {

  vtkDataSet *domain = vtkDataSet::GetData(inputVector[0]);
  vtkPointSet *sources = vtkPointSet::GetData(inputVector[1]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  int ret{};

  ret = getTriangulation(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    printErr("[ttkDistanceField] Error : wrong triangulation.");
    return -1;
  }
#endif

  ret = getIdentifiers(sources);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    printErr("[ttkDistanceField] Error : wrong identifiers.");
    return -2;
  }
#endif

  const SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfPointsInDomain) {
    printErr("[ttkDistanceField] Error : domain has no points.");
    return -3;
  }
#endif

  const SimplexId numberOfPointsInSources = sources->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfPointsInSources) {
    printErr("[ttkDistanceField] Error : sources have no points.");
    return -4;
  }
#endif

  vtkSmartPointer<ttkSimplexIdTypeArray> origin = vtkSmartPointer<ttkSimplexIdTypeArray>::New();

  if(origin) {
    origin->SetNumberOfComponents(1);
    origin->SetNumberOfTuples(numberOfPointsInDomain);
    origin->SetName(ttk::VertexScalarFieldName);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  else {
    printErr("[ttkDistanceField] Error : ttkSimplexIdTypeArray allocation problem.");
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
    printErr("[ttkDistanceField] Error : ttkSimplexIdTypeArray allocation problem.");
    return -6;
  }
#endif

  this->setVertexNumber(numberOfPointsInDomain);
  this->setSourceNumber(numberOfPointsInSources);

  this->setVertexIdentifierScalarFieldPointer(ttkUtils::GetVoidPointer(identifiers_));
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
        printErr("[ttkDistanceField] Error : vtkFloatArray allocation problem.");
        return -7;
      }
#endif

      this->setOutputScalarFieldPointer(ttkUtils::GetVoidPointer(distanceScalars));
      // @PETER I don't think there's any need to use ttkTemplateMacro here because the triangulation is not part of the function call
      ret = this->execute<float>();
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
        printErr("[ttkDistanceField] Error : vtkDoubleArray allocation problem.");
        return -8;
      }
#endif
      this->setOutputScalarFieldPointer(ttkUtils::GetVoidPointer(distanceScalars));
      // @PETER I don't think there's any need to use ttkTemplateMacro here because the triangulation is not part of the function call
      ret = this->execute<double>();
      break;

    default:
#ifndef TTK_ENABLE_KAMIKAZE
      printErr("[ttkDistanceField] Error : Scalar field type problem.");
      return -9;
#endif
      break;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(ret) {
    printErr("[ttkDistanceField] DistanceField.execute() error code : " + std::to_string(ret));
    return -10;
  }
#endif

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(distanceScalars);
  output->GetPointData()->AddArray(origin);
  output->GetPointData()->AddArray(seg);
  distanceScalars->Delete();

  printMsg("The ret is " + std::to_string(ret) , ttk::Debug::debugPriority::memoryMsg);

  return 1;
}
