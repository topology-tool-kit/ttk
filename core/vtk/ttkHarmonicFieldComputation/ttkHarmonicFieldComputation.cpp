#include <ttkHarmonicFieldComputation.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkHarmonicFieldComputation)

    ttkHarmonicFieldComputation::ttkHarmonicFieldComputation()
    : hasUpdatedMesh_{false}, identifiers_{}, inputScalars_{}, offsets_{},
      inputOffsets_{} {
  SetNumberOfInputPorts(2);
  triangulation_ = NULL;

  ScalarFieldId = 0;
  OffsetFieldId = -1;
  OutputOffsetScalarFieldName = ttk::OffsetScalarFieldName;
  InputVertexScalarFieldName = ttk::VertexScalarFieldName;
  InputOffsetScalarFieldName = ttk::OffsetScalarFieldName;

  UseAllCores = true;
}

ttkHarmonicFieldComputation::~ttkHarmonicFieldComputation() {
  if (offsets_)
    offsets_->Delete();
}

int ttkHarmonicFieldComputation::FillInputPortInformation(
    int port, vtkInformation *info) {

  if (port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");

  if (port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");

  return 1;
}

int ttkHarmonicFieldComputation::getTriangulation(vtkDataSet *input) {

  triangulation_ = ttkTriangulation::getTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE
  if (!triangulation_) {
    cerr << "[ttkHarmonicFieldComputation] Error : input triangulation pointer "
            "is NULL."
         << endl;
    return -1;
  }
#endif

  triangulation_->setWrapper(this);
  harmonicFieldComputation_.setWrapper(this);
  harmonicFieldComputation_.setupTriangulation(triangulation_);
  Modified();
  hasUpdatedMesh_ = true;

#ifndef TTK_ENABLE_KAMIKAZE
  if (triangulation_->isEmpty()) {
    cerr << "[ttkHarmonicFieldComputation] Error : ttkTriangulation allocation "
            "problem."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkHarmonicFieldComputation::getScalars(vtkDataSet *input) {
#ifndef TTK_ENABLE_KAMIKAZE
  if (!input) {
    cerr << "[ttkHarmonicFieldComputation] Error : input pointer is NULL."
         << endl;
    return -1;
  }

  if (!input->GetNumberOfPoints()) {
    cerr << "[ttkHarmonicFieldComputation] Error : input has no point." << endl;
    return -1;
  }
#endif

  vtkPointData *pointData = input->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if (!pointData) {
    cerr << "[ttkHarmonicFieldComputation] Error : input has no point data."
         << endl;
    return -1;
  }
#endif

  if (ScalarField.length()) {
    inputScalars_ = pointData->GetArray(ScalarField.data());
  } else {
    inputScalars_ = pointData->GetArray(ScalarFieldId);
    if (inputScalars_)
      ScalarField = inputScalars_->GetName();
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if (!inputScalars_) {
    cerr << "[ttkHarmonicFieldComputation] Error : input scalar field pointer "
            "is null."
         << endl;
    return -3;
  }
#endif

  return 0;
}

int ttkHarmonicFieldComputation::getIdentifiers(vtkPointSet *input) {
  if (InputVertexScalarFieldName.length())
    identifiers_ =
        input->GetPointData()->GetArray(InputVertexScalarFieldName.data());
  else if (input->GetPointData()->GetArray(ttk::VertexScalarFieldName))
    identifiers_ = input->GetPointData()->GetArray(ttk::VertexScalarFieldName);

#ifndef TTK_ENABLE_KAMIKAZE
  if (!identifiers_) {
    cerr << "[ttkHarmonicFieldComputation] Error : wrong vertex identifier "
            "scalar field."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkHarmonicFieldComputation::getOffsets(vtkDataSet *input) {
  if (InputOffsetScalarFieldName.length()) {
    inputOffsets_ =
        input->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  } else if (OffsetFieldId != -1 and
             input->GetPointData()->GetArray(OffsetFieldId)) {
    inputOffsets_ = input->GetPointData()->GetArray(OffsetFieldId);
  } else if (input->GetPointData()->GetArray(ttk::OffsetScalarFieldName)) {
    inputOffsets_ = input->GetPointData()->GetArray(ttk::OffsetScalarFieldName);
  } else {
    if (hasUpdatedMesh_ and offsets_) {
      offsets_->Delete();
      offsets_ = nullptr;
      hasUpdatedMesh_ = false;
    }

    if (!offsets_) {
      const SimplexId numberOfVertices = input->GetNumberOfPoints();

      offsets_ = ttkSimplexIdTypeArray::New();
      offsets_->SetNumberOfComponents(1);
      offsets_->SetNumberOfTuples(numberOfVertices);
      offsets_->SetName(ttk::OffsetScalarFieldName);
      for (SimplexId i = 0; i < numberOfVertices; ++i)
        offsets_->SetTuple1(i, i);
    }

    inputOffsets_ = offsets_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if (!inputOffsets_) {
    cerr << "[ttkHarmonicFieldComputation] Error : wrong input offset scalar "
            "field."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkHarmonicFieldComputation::doIt(vector<vtkDataSet *> &inputs,
                                      vector<vtkDataSet *> &outputs) {

  Memory m;

  vtkDataSet *domain = inputs[0];
  vtkPointSet *constraints = vtkPointSet::SafeDownCast(inputs[1]);
  vtkDataSet *output = outputs[0];

  int ret{};

  ret = getTriangulation(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  if (ret) {
    cerr << "[ttkHarmonicFieldComputation] Error : wrong triangulation."
         << endl;
    return -1;
  }
#endif

  ret = getScalars(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  if (ret) {
    cerr << "[ttkHarmonicFieldComputation] Error : wrong scalars." << endl;
    return -2;
  }
#endif

  ret = getIdentifiers(constraints);
#ifndef TTK_ENABLE_KAMIKAZE
  if (ret) {
    cerr << "[ttkHarmonicFieldComputation] Error : wrong identifiers." << endl;
    return -3;
  }
#endif

  ret = getOffsets(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  if (ret) {
    cerr << "[ttkHarmonicFieldComputation] Error : wrong offsets." << endl;
    return -4;
  }
#endif
#ifndef TTK_ENABLE_KAMIKAZE
  if (inputOffsets_->GetDataType() != VTK_INT and
      inputOffsets_->GetDataType() != VTK_ID_TYPE) {
    cerr << "[ttkHarmonicFieldComputation] Error : input offset field type not "
            "supported."
         << endl;
    return -1;
  }
#endif

  const SimplexId numberOfVertices = domain->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if (numberOfVertices <= 0) {
    cerr << "[ttkHarmonicFieldComputation] Error : domain has no points."
         << endl;
    return -5;
  }
#endif

  if (OutputOffsetScalarFieldName.length() <= 0)
    OutputOffsetScalarFieldName = ttk::OffsetScalarFieldName;

  vtkSmartPointer<ttkSimplexIdTypeArray> outputOffsets =
      vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  if (outputOffsets) {
    outputOffsets->SetNumberOfComponents(1);
    outputOffsets->SetNumberOfTuples(numberOfVertices);
    outputOffsets->SetName(OutputOffsetScalarFieldName.data());
  }
#ifndef TTK_ENABLE_KAMIKAZE
  else {
    cerr << "[ttkHarmonicFieldComputation] Error : ttkSimplexIdTypeArray "
            "allocation problem."
         << endl;
    return -7;
  }
#endif

  vtkDataArray *outputScalars{};
  switch (inputScalars_->GetDataType()) {
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
    cerr << "[ttkHarmonicFieldComputation] Error : Unsupported data type."
         << endl;
    return -8;
#endif
  }
  if (outputScalars) {
    outputScalars->SetNumberOfTuples(numberOfVertices);
    outputScalars->SetName(inputScalars_->GetName());
  }
#ifndef TTK_ENABLE_KAMIKAZE
  else {
    cerr << "[ttkHarmonicFieldComputation] Error : vtkDataArray allocation "
            "problem."
         << endl;
    return -9;
  }
#endif

  const SimplexId numberOfConstraints = constraints->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if (numberOfConstraints <= 0) {
    cerr << "[ttkHarmonicFieldComputation] Error : input has no constraints."
         << endl;
    return -10;
  }
#endif

  harmonicFieldComputation_.setVertexNumber(numberOfVertices);
  harmonicFieldComputation_.setConstraintNumber(numberOfConstraints);
  harmonicFieldComputation_.setInputScalarFieldPointer(
      inputScalars_->GetVoidPointer(0));
  harmonicFieldComputation_.setVertexIdentifierScalarFieldPointer(
      identifiers_->GetVoidPointer(0));

  harmonicFieldComputation_.setOutputScalarFieldPointer(
      outputScalars->GetVoidPointer(0));

#ifndef TTK_ENABLE_KAMIKAZE
  if (identifiers_->GetDataType() != inputOffsets_->GetDataType()) {
    cerr << "[ttkHarmonicFieldComputation] Error : type of identifiers and "
            "offsets are different."
         << endl;
    return -11;
  }
#endif

  switch (inputScalars_->GetDataType()) {
    ttkTemplateMacro({
      if (inputOffsets_->GetDataType() == VTK_INT)
        ret = harmonicFieldComputation_.execute<VTK_TT>();
      if (inputOffsets_->GetDataType() == VTK_ID_TYPE)
        ret = harmonicFieldComputation_.execute<VTK_TT>();
    });
  }
#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if (ret) {
    cerr << "[ttkHarmonicFieldComputation] HarmonicFieldComputation.execute() "
            "error code : "
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
    msg << "[ttkHarmonicFieldComputation] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
