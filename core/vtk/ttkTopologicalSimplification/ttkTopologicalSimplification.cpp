#include <ttkTopologicalSimplification.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkTopologicalSimplification)

  ttkTopologicalSimplification::ttkTopologicalSimplification()
  : hasUpdatedMesh_{false}, identifiers_{}, inputScalars_{}, offsets_{},
    inputOffsets_{} {
  SetNumberOfInputPorts(2);
  triangulation_ = NULL;

  ScalarFieldId = 0;
  OffsetFieldId = -1;
  ForceInputOffsetScalarField = false;
  AddPerturbation = false;
  OutputOffsetScalarFieldName = ttk::OffsetScalarFieldName;
  ForceInputVertexScalarField = false;
  InputVertexScalarFieldName = ttk::VertexScalarFieldName;
  ConsiderIdentifierAsBlackList = false;
  InputOffsetScalarFieldName = ttk::OffsetScalarFieldName;

  UseAllCores = true;
}

ttkTopologicalSimplification::~ttkTopologicalSimplification() {
  if(offsets_)
    offsets_->Delete();
}

int ttkTopologicalSimplification::FillInputPortInformation(
  int port, vtkInformation *info) {

  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");

  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");

  return 1;
}

int ttkTopologicalSimplification::getTriangulation(vtkDataSet *input) {

  triangulation_ = ttkTriangulation::getTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_) {
    cerr << "[ttkTopologicalSimplification] Error : input triangulation "
            "pointer is NULL."
         << endl;
    return -1;
  }
#endif

  triangulation_->setWrapper(this);
  topologicalSimplification_.setWrapper(this);
  topologicalSimplification_.setupTriangulation(triangulation_);
  Modified();
  hasUpdatedMesh_ = true;

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_->isEmpty()) {
    cerr << "[ttkTopologicalSimplification] Error : ttkTriangulation "
            "allocation problem."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkTopologicalSimplification::getScalars(vtkDataSet *input) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    cerr << "[ttkTopologicalSimplification] Error : input pointer is NULL."
         << endl;
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    cerr << "[ttkTopologicalSimplification] Error : input has no point."
         << endl;
    return -1;
  }
#endif

  vtkPointData *pointData = input->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    cerr << "[ttkTopologicalSimplification] Error : input has no point data."
         << endl;
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
    cerr << "[ttkTopologicalSimplification] Error : input scalar field pointer "
            "is null."
         << endl;
    return -3;
  }
#endif

  return 0;
}

int ttkTopologicalSimplification::getIdentifiers(vtkPointSet *input) {
  if(ForceInputVertexScalarField and InputVertexScalarFieldName.length())
    identifiers_
      = input->GetPointData()->GetArray(InputVertexScalarFieldName.data());
  else if(input->GetPointData()->GetArray(ttk::VertexScalarFieldName))
    identifiers_ = input->GetPointData()->GetArray(ttk::VertexScalarFieldName);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!identifiers_) {
    cerr << "[ttkTopologicalSimplification] Error : wrong vertex identifier "
            "scalar field."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkTopologicalSimplification::getOffsets(vtkDataSet *input) {
  if(ForceInputOffsetScalarField and InputOffsetScalarFieldName.length()) {
    inputOffsets_
      = input->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  } else if(OffsetFieldId != -1
            and input->GetPointData()->GetArray(OffsetFieldId)) {
    inputOffsets_ = input->GetPointData()->GetArray(OffsetFieldId);
  } else if(input->GetPointData()->GetArray(ttk::OffsetScalarFieldName)) {
    inputOffsets_ = input->GetPointData()->GetArray(ttk::OffsetScalarFieldName);
  } else {
    if(hasUpdatedMesh_ and offsets_) {
      offsets_->Delete();
      offsets_ = nullptr;
      hasUpdatedMesh_ = false;
    }

    if(!offsets_) {
      const SimplexId numberOfVertices = input->GetNumberOfPoints();

      offsets_ = ttkSimplexIdTypeArray::New();
      offsets_->SetNumberOfComponents(1);
      offsets_->SetNumberOfTuples(numberOfVertices);
      offsets_->SetName(ttk::OffsetScalarFieldName);
      for(SimplexId i = 0; i < numberOfVertices; ++i)
        offsets_->SetTuple1(i, i);
    }

    inputOffsets_ = offsets_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOffsets_) {
    cerr << "[ttkTopologicalSimplification] Error : wrong input offset scalar "
            "field."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkTopologicalSimplification::doIt(vector<vtkDataSet *> &inputs,
                                       vector<vtkDataSet *> &outputs) {

  Memory m;

  vtkDataSet *domain = inputs[0];
  vtkPointSet *constraints = vtkPointSet::SafeDownCast(inputs[1]);
  vtkDataSet *output = outputs[0];

  int ret{};

  ret = getTriangulation(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkTopologicalSimplification] Error : wrong triangulation."
         << endl;
    return -1;
  }
#endif

  ret = getScalars(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkTopologicalSimplification] Error : wrong scalars." << endl;
    return -2;
  }
#endif

  ret = getIdentifiers(constraints);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkTopologicalSimplification] Error : wrong identifiers." << endl;
    return -3;
  }
#endif

  ret = getOffsets(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkTopologicalSimplification] Error : wrong offsets." << endl;
    return -4;
  }
#endif
#ifndef TTK_ENABLE_KAMIKAZE
  if(inputOffsets_->GetDataType() != VTK_INT
     and inputOffsets_->GetDataType() != VTK_ID_TYPE) {
    cerr << "[ttkTopologicalSimplification] Error : input offset field type "
            "not supported."
         << endl;
    return -1;
  }
#endif

  const SimplexId numberOfVertices = domain->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfVertices <= 0) {
    cerr << "[ttkTopologicalSimplification] Error : domain has no points."
         << endl;
    return -5;
  }
#endif

  if(OutputOffsetScalarFieldName.length() <= 0)
    OutputOffsetScalarFieldName = ttk::OffsetScalarFieldName;

  vtkSmartPointer<ttkSimplexIdTypeArray> outputOffsets
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  if(outputOffsets) {
    outputOffsets->SetNumberOfComponents(1);
    outputOffsets->SetNumberOfTuples(numberOfVertices);
    outputOffsets->SetName(OutputOffsetScalarFieldName.data());
  }
#ifndef TTK_ENABLE_KAMIKAZE
  else {
    cerr << "[ttkTopologicalSimplification] Error : ttkSimplexIdTypeArray "
            "allocation problem."
         << endl;
    return -7;
  }
#endif

  vtkDataArray *outputScalars{};
  switch(inputScalars_->GetDataType()) {
    case VTK_DOUBLE:
      outputScalars = vtkDoubleArray::New();
      break;

    case VTK_FLOAT:
      outputScalars = vtkFloatArray::New();
      break;

    case VTK_INT:
      outputScalars = vtkIntArray::New();
      break;

    case VTK_ID_TYPE:
      outputScalars = vtkIdTypeArray::New();
      break;

    case VTK_SHORT:
      outputScalars = vtkShortArray::New();
      break;

    case VTK_UNSIGNED_SHORT:
      outputScalars = vtkUnsignedShortArray::New();
      break;

    case VTK_CHAR:
      outputScalars = vtkCharArray::New();
      break;

    case VTK_UNSIGNED_CHAR:
      outputScalars = vtkUnsignedCharArray::New();
      break;

#ifndef TTK_ENABLE_KAMIKAZE
    default:
      cerr << "[ttkTopologicalSimplification] Error : Unsupported data type."
           << endl;
      return -8;
#endif
  }
  if(outputScalars) {
    outputScalars->SetNumberOfTuples(numberOfVertices);
    outputScalars->SetName(inputScalars_->GetName());
  }
#ifndef TTK_ENABLE_KAMIKAZE
  else {
    cerr << "[ttkTopologicalSimplification] Error : vtkDataArray allocation "
            "problem."
         << endl;
    return -9;
  }
#endif

  const SimplexId numberOfConstraints = constraints->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfConstraints <= 0) {
    cerr << "[ttkTopologicalSimplification] Error : input has no constraints."
         << endl;
    return -10;
  }
#endif

  topologicalSimplification_.setVertexNumber(numberOfVertices);
  topologicalSimplification_.setConstraintNumber(numberOfConstraints);
  topologicalSimplification_.setInputScalarFieldPointer(
    inputScalars_->GetVoidPointer(0));
  topologicalSimplification_.setVertexIdentifierScalarFieldPointer(
    identifiers_->GetVoidPointer(0));
  topologicalSimplification_.setConsiderIdentifierAsBlackList(
    ConsiderIdentifierAsBlackList);
  topologicalSimplification_.setAddPerturbation(AddPerturbation);

  topologicalSimplification_.setInputOffsetScalarFieldPointer(
    inputOffsets_->GetVoidPointer(0));

  topologicalSimplification_.setOutputScalarFieldPointer(
    outputScalars->GetVoidPointer(0));

  topologicalSimplification_.setOutputOffsetScalarFieldPointer(
    outputOffsets->GetVoidPointer(0));

#ifndef TTK_ENABLE_KAMIKAZE
  if(identifiers_->GetDataType() != inputOffsets_->GetDataType()) {
    cerr << "[ttkTopologicalSimplification] Error : type of identifiers and "
            "offsets are different."
         << endl;
    return -11;
  }
#endif

  switch(inputScalars_->GetDataType()) {
    ttkTemplateMacro({
      if(inputOffsets_->GetDataType() == VTK_INT)
        ret = topologicalSimplification_.execute<VTK_TT TTK_COMMA int>();
      if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
        ret = topologicalSimplification_.execute<VTK_TT TTK_COMMA vtkIdType>();
    });
  }
#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(ret) {
    cerr << "[ttkTopologicalSimplification] "
            "TopologicalSimplification.execute() error code : "
         << ret << endl;
    return -12;
  }
#endif

  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(outputOffsets);
  output->GetPointData()->AddArray(outputScalars);
  outputScalars->Delete();

  {
    stringstream msg;
    msg << "[ttkTopologicalSimplification] Memory usage: "
        << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
