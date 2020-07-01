#include <ttkUtils.h>
#include <ttkMacros.h>

#include <ttkDistanceField.h>
#include <vtkPointSet.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkDistanceField)

ttkDistanceField::ttkDistanceField() : identifiers_{} {
  //OutputScalarFieldType = 0;
  //OutputScalarFieldName = "OutputDistanceField";
  //ForceInputVertexScalarField = false;
  //InputVertexScalarFieldName = ttk::VertexScalarFieldName;
  //UseAllCores = true;
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);

  //triangulation_ = NULL;
}

ttkDistanceField::~ttkDistanceField() {
}

int ttkDistanceField::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");

  return 1;
}

int ttkDistanceField::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
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
    cerr << "[ttkDistanceField] Error : wrong vertex identifiers." << endl;
    return -1;
  }
#endif

  return 0;
}

//int ttkDistanceField::doIt(vector<vtkDataSet *> &inputs,
                           //vector<vtkDataSet *> &outputs) {
int ttkDistanceField::RequestData(vtkInformation *request,vtkInformationVector **inputVector,vtkInformationVector *outputVector) {

  Memory m;

  //vtkDataSet *domain = vtkDataSet::SafeDownCast(inputs[0]);
  //vtkPointSet *sources = vtkPointSet::SafeDownCast(inputs[1]);
  //vtkDataSet *output = vtkDataSet::SafeDownCast(outputs[0]);

  vtkDataSet *domain = vtkDataSet::GetData(inputVector[0]);
  vtkPointSet *sources = vtkPointSet::GetData(inputVector[1]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  int ret{};

  ret = getTriangulation(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkDistanceField] Error : wrong triangulation." << endl;
    return -1;
  }
#endif

  ret = getIdentifiers(sources);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkDistanceField] Error : wrong identifiers." << endl;
    return -2;
  }
#endif

  const SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfPointsInDomain) {
    cerr << "[ttkDistanceField] Error : domain has no points." << endl;
    return -3;
  }
#endif

  const SimplexId numberOfPointsInSources = sources->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfPointsInSources) {
    cerr << "[ttkDistanceField] Error : sources have no points." << endl;
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
    cerr
      << "[ttkDistanceField] Error : ttkSimplexIdTypeArray allocation problem."
      << endl;
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
    cerr
      << "[ttkDistanceField] Error : ttkSimplexIdTypeArray allocation problem."
      << endl;
    return -6;
  }
#endif

  this->setVertexNumber(numberOfPointsInDomain);
  this->setSourceNumber(numberOfPointsInSources);

  this->setVertexIdentifierScalarFieldPointer(ttkUtils::GetVoidPointer(identifiers_));
  this->setOutputIdentifiers(ttkUtils::GetVoidPointer(origin));
  this->setOutputSegmentation(ttkUtils::GetVoidPointer(seg));

  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);

  vtkSmartPointer<vtkDataArray> distanceScalars = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  distanceScalars->SetNumberOfComponents(1); // only one component per tuple
  distanceScalars->SetNumberOfTuples(numberOfPointsInDomain);
  distanceScalars->SetName(OutputScalarFieldName.data());
  this->setOutputScalarFieldPointer(ttkUtils::GetVoidPointer(distanceScalars));

  int status = 0;
  ttkVtkTemplateMacro(this->triangulation_->getType(), inputArray->GetDataType(), (status = this->execute<VTK_TT>()));

  // On error cancel filter execution
  if(status == 0)
    return 0;

  //vtkDataArray *distanceScalars{};

  //vtkDataArray *distanceScalars{};

  //distanceScalars = vtkFloatArray::New();
  //if(distanceScalars) {
      //distanceScalars->SetNumberOfComponents(1);
      //distanceScalars->SetNumberOfTuples(numberOfPointsInDomain);
      //distanceScalars->SetName(OutputScalarFieldName.data());

  //this->setOutputScalarFieldPointer(ttkUtils::GetVoidPointer(distanceScalars));

  //switch(OutputScalarFieldType) {

    //case DistanceType::Float:
      //distanceScalars = vtkFloatArray::New();
      //if(distanceScalars) {
        //distanceScalars->SetNumberOfComponents(1);
        //distanceScalars->SetNumberOfTuples(numberOfPointsInDomain);
        //distanceScalars->SetName(OutputScalarFieldName.data());
      //}
//#ifndef TTK_ENABLE_KAMIKAZE
      //else {
        //cerr << "[ttkDistanceField] Error : vtkFloatArray allocation problem."
             //<< endl;
        //return -7;
      //}
//#endif

      //this->setOutputScalarFieldPointer(ttkUtils::GetVoidPointer(distanceScalars));
      //ret = this->execute<float>();
      //break;

    //case DistanceType::Double:
      //distanceScalars = vtkDoubleArray::New();
      //if(distanceScalars) {
        //distanceScalars->SetNumberOfComponents(1);
        //distanceScalars->SetNumberOfTuples(numberOfPointsInDomain);
        //distanceScalars->SetName(OutputScalarFieldName.data());
      //}
//#ifndef TTK_ENABLE_KAMIKAZE
      //else {
        //cerr << "[ttkDistanceField] Error : vtkDoubleArray allocation problem."
             //<< endl;
        //return -8;
      //}
//#endif

      //this->setOutputScalarFieldPointer(ttkUtils::GetVoidPointer(distanceScalars));
      //ret = this->execute<double>();
      //break;

    //default:
//#ifndef TTK_ENABLE_KAMIKAZE
      //cerr << "[ttkDistanceField] Error : Scalar field type problem." << endl;
      //return -9;
//#endif
      //break;
  //}

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(ret) {
    cerr << "[ttkDistanceField] DistanceField.execute() error code : " << ret
         << endl;
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
    stringstream msg;
    msg << "[ttkDistanceField] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    //dMsg(cout, msg.str(), memoryMsg);
    printMsg(msg.str(), ttk::Debug::debugPriority::memoryMsg);
  }

  return ret;
}
