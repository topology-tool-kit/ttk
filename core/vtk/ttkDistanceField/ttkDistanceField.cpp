#include <ttkDistanceField.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkDistanceField)

  ttkDistanceField::ttkDistanceField()
  : identifiers_{} {
  OutputScalarFieldType = 0;
  OutputScalarFieldName = "OutputDistanceField";
  ForceInputVertexScalarField = false;
  InputVertexScalarFieldName = ttk::VertexScalarFieldName;
  UseAllCores = true;
  SetNumberOfInputPorts(2);

  triangulation_ = NULL;
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

int ttkDistanceField::getTriangulation(vtkDataSet *input) {

  triangulation_ = ttkTriangulation::getTriangulation(input);
  triangulation_->setWrapper(this);
  distanceField_.setWrapper(this);

  distanceField_.setupTriangulation(triangulation_);
  Modified();

#ifndef TTK_ENABLE_KAMIKAZE
  // allocation problem
  if(triangulation_->isEmpty()) {
    cerr << "[ttkDistanceField] Error : ttkTriangulation allocation problem."
         << endl;
    return -1;
  }
#endif

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

int ttkDistanceField::doIt(vector<vtkDataSet *> &inputs,
                           vector<vtkDataSet *> &outputs) {

  Memory m;

  vtkDataSet *domain = vtkDataSet::SafeDownCast(inputs[0]);
  vtkPointSet *sources = vtkPointSet::SafeDownCast(inputs[1]);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outputs[0]);

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

  vtkSmartPointer<ttkSimplexIdTypeArray> origin
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
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

  distanceField_.setVertexNumber(numberOfPointsInDomain);
  distanceField_.setSourceNumber(numberOfPointsInSources);

  distanceField_.setVertexIdentifierScalarFieldPointer(
    identifiers_->GetVoidPointer(0));
  distanceField_.setOutputIdentifiers(origin->GetVoidPointer(0));
  distanceField_.setOutputSegmentation(seg->GetVoidPointer(0));

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
        cerr << "[ttkDistanceField] Error : vtkFloatArray allocation problem."
             << endl;
        return -7;
      }
#endif

      distanceField_.setOutputScalarFieldPointer(
        distanceScalars->GetVoidPointer(0));
      ret = distanceField_.execute<float>();
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
        cerr << "[ttkDistanceField] Error : vtkDoubleArray allocation problem."
             << endl;
        return -8;
      }
#endif

      distanceField_.setOutputScalarFieldPointer(
        distanceScalars->GetVoidPointer(0));
      ret = distanceField_.execute<double>();
      break;

    default:
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkDistanceField] Error : Scalar field type problem." << endl;
      return -9;
#endif
      break;
  }

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
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
