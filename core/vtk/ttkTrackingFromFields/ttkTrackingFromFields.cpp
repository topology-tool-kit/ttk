#include                  <ttkTrackingFromFields.h>

vtkStandardNewMacro(ttkTrackingFromFields)

constexpr unsigned int str2int(const char* str, int h = 0)
{
  return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

int ttkTrackingFromFields::FillOutputPortInformation(int port, vtkInformation* info)
{
  switch (port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    // default:
    //   break;
  }

  return 1;
}

int ttkTrackingFromFields::doIt(
  std::vector<vtkDataSet *> &inputs,
  std::vector<vtkDataSet *> &outputs)
{
  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  // triangulation_.setWrapper(this);
  // triangulation_.setInputData(input);
  internalTriangulation_ = ttkTriangulation::getTriangulation(input);
  internalTriangulation_->setWrapper(this);

  // Test validity of datasets (must present the same number of points).
  if (!input)
    return -1;

  // Get number and list of inputs.
  std::vector<vtkDataArray*> inputScalarFieldsRaw;
  std::vector<vtkDataArray*> inputScalarFields;
  int numberOfInputFields = input->GetPointData()->GetNumberOfArrays();
  if (numberOfInputFields < 3) {
    std::stringstream msg;
    msg << "[ttkTrackingFromField] not enough input fields to perform tracking." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  vtkDataArray* firstScalarField = input->GetPointData()->GetArray(0);
  // inputScalarFields.push_back(firstScalarField);

  for (int i = 0; i < numberOfInputFields; ++i) {
    vtkDataArray* currentScalarField = input->GetPointData()->GetArray(i);
    if (!currentScalarField ||
        firstScalarField->GetDataType() != currentScalarField->GetDataType())
    {
      std::stringstream msg;
      msg << "[ttkTrackingFromField] inconsistent field data type or size (" << i << ")." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      return -1;
    }
    inputScalarFieldsRaw.push_back(currentScalarField);
  }

  std::sort(
      inputScalarFieldsRaw.begin(),
      inputScalarFieldsRaw.end(),
      [] (vtkDataArray* a,
          vtkDataArray* b)
      {
          std::string s1 = a->GetName();
          std::string s2 = b->GetName();
          return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end());
      });

  int end = EndTimestep <= 0 ? numberOfInputFields : std::min(numberOfInputFields, EndTimestep);
  for (int i = StartTimestep; i < end; i += Sampling) {
    vtkDataArray* currentScalarField = inputScalarFieldsRaw[i];
    std::cout << currentScalarField->GetName() << std::endl;
    inputScalarFields.push_back(currentScalarField);
  }

  // Input -> persistence filter.
  std::string algorithm = DistanceAlgorithm;
  int pvalg = PVAlgorithm;
  bool useTTKMethod = false;

  if (pvalg >= 0) {
    switch (pvalg) {
      case 0:
      case 1:
      case 2:
      case 3:
        useTTKMethod = true;
        break;
      case 4:
        break;
      default:
        std::stringstream msg;
        msg << "[ttkTrackingFromFieldsFromField] Unrecognized tracking method." << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
        break;
    }
  } else {
    switch (str2int(algorithm.c_str())) {
      case str2int("0"):
      case str2int("ttk"):
      case str2int("1"):
      case str2int("legacy"):
      case str2int("2"):
      case str2int("geometric"):
      case str2int("3"):
      case str2int("parallel"):
        useTTKMethod = true;
        break;
      case str2int("4"):
      case str2int("greedy"):
        break;
      default:
        std::stringstream msg;
        msg << "[ttkTrackingFromField] Unrecognized tracking method." << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
        break;
    }
  }

  int res = 0;
  if (useTTKMethod)
    res = trackWithPersistenceMatching(
        input, output, inputScalarFields);
  else
  {
    std::stringstream msg;
    msg << "[ttkTrackingFromField] The specified matching method does not." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return res;
}

// (*) Persistence-driven approach
int ttkTrackingFromFields::trackWithPersistenceMatching(
    vtkDataSet *input,
    vtkUnstructuredGrid *output,
    std::vector<vtkDataArray*> inputScalarFields)
{
  using dataType = double;
  unsigned long fieldNumber = inputScalarFields.size();

  // 0. get data
  trackingF_.setTriangulation(internalTriangulation_);
  std::vector<void *> inputFields(fieldNumber);
  for (int i = 0; i < (int) fieldNumber; ++i)
    inputFields[i] = inputScalarFields[i]->GetVoidPointer(0);
  trackingF_.setInputScalars(inputFields);

  // 0'. get offsets
  auto numberOfVertices = (int) input->GetNumberOfPoints();
  vtkIdTypeArray* offsets_ = vtkIdTypeArray::New();
  offsets_->SetNumberOfComponents(1);
  offsets_->SetNumberOfTuples(numberOfVertices);
  offsets_->SetName("OffsetScalarField");
  for(int i = 0; i < numberOfVertices; ++i)
    offsets_->SetTuple1(i,i);
  trackingF_.setInputOffsets(offsets_->GetVoidPointer(0));

  // 1. get persistence diagrams.
  auto persistenceDiagrams = new std::vector<std::vector<diagramTuple>*>(fieldNumber);
//  std::vector<vtkUnstructuredGrid*> persistenceDiagrams(fieldNumber);

  trackingF_.performDiagramComputation((int) fieldNumber, persistenceDiagrams, this);

  // 2. call feature tracking with threshold.
  auto outputMatchings =
    new std::vector<std::vector<matchingTuple>*>(fieldNumber - 1);
  for (unsigned long i = 0; i < fieldNumber - 1; ++i) {
      outputMatchings->at(i) = new std::vector<matchingTuple>();
  }

  double spacing = Spacing;
  std::string algorithm = DistanceAlgorithm;
  double alpha = Alpha;
  double tolerance = Tolerance;
  bool is3D = Is3D;
  std::string wasserstein = WassersteinMetric;

  tracking_.performMatchings(
      (int) fieldNumber,
      persistenceDiagrams,
      outputMatchings,
      algorithm, // Not from paraview, from enclosing tracking plugin
      wasserstein,
      tolerance,
      is3D,
      alpha, // Blending
      PX, PY, PZ, PS, PE, // Coefficients
      this // Wrapper for accessing threadNumber
  );

  outputMesh_ = vtkUnstructuredGrid::New();
  vtkUnstructuredGrid *outputMesh = vtkUnstructuredGrid::SafeDownCast(output);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkDoubleArray> persistenceScalars = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> valueScalars = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> matchingIdScalars = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> timeScalars = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> componentIds = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> pointTypeScalars = vtkSmartPointer<vtkIntArray>::New();
  persistenceScalars->SetName("Cost");
  valueScalars->SetName("Scalar");
  matchingIdScalars->SetName("MatchingIdentifier");
  timeScalars->SetName("TimeStep");
  componentIds->SetName("ConnectedComponentId");
  pointTypeScalars->SetName("NodeType");

  // (+ vertex id)
  auto trackingsBase = new std::vector<trackingTuple>(); // trackings;
  tracking_.performTracking(
      persistenceDiagrams, outputMatchings,
      trackingsBase);

  auto trackingTupleToMerged = new std::vector<std::set<int>>(trackingsBase->size(), std::set<int>());
  if (DoPostProc)
    tracking_.performPostProcess(
      persistenceDiagrams,
      trackingsBase,
      trackingTupleToMerged,
      PostProcThresh);

  bool useGeometricSpacing = UseGeometricSpacing;
  std::vector<trackingTuple> trackings = *trackingsBase;

  // Build mesh.
  ttkTrackingFromPersistenceDiagrams::buildMesh(
    trackings, outputMatchings, persistenceDiagrams,
    useGeometricSpacing, spacing, DoPostProc, trackingTupleToMerged,
    points, persistenceDiagram,
    persistenceScalars, valueScalars, matchingIdScalars, timeScalars,
    componentIds, pointTypeScalars);

  outputMesh_->ShallowCopy(persistenceDiagram);
  outputMesh->ShallowCopy(outputMesh_);

  return 0;
}
