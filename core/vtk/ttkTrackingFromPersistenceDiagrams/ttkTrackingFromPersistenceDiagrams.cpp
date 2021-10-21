#include <ttkTrackingFromPersistenceDiagrams.h>

vtkStandardNewMacro(ttkTrackingFromPersistenceDiagrams)

  using dataType = double;

ttkTrackingFromPersistenceDiagrams::ttkTrackingFromPersistenceDiagrams() {

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

ttkTrackingFromPersistenceDiagrams::~ttkTrackingFromPersistenceDiagrams() {
  if(outputMesh_)
    outputMesh_->Delete();
}

int ttkTrackingFromPersistenceDiagrams::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

int ttkTrackingFromPersistenceDiagrams::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkTrackingFromPersistenceDiagrams::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  // Output pointers informations
  vtkInformation *outInfo;

  // Unified bound fields
  outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid *mesh = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkUnstructuredGrid::DATA_OBJECT()));

  // Number of input files
  int numInputs = inputVector[0]->GetNumberOfInformationObjects();
  this->printMsg("Number of inputs: " + std::to_string(numInputs));

  // Get input data
  std::vector<vtkDataSet *> input(numInputs);
  for(int i = 0; i < numInputs; i++) {
    input[i] = vtkDataSet::GetData(inputVector[0], i);
  }

  std::vector<std::vector<diagramTuple>> inputPersistenceDiagrams(
    (unsigned long)numInputs, std::vector<diagramTuple>());

  std::vector<vtkSmartPointer<vtkUnstructuredGrid>> outputPersistenceDiagrams(
    (unsigned long)2 * numInputs - 2,
    vtkSmartPointer<vtkUnstructuredGrid>::New());

  std::vector<std::vector<matchingTuple>> outputMatchings(
    (unsigned long)numInputs - 1, std::vector<matchingTuple>());

  // Input parameters.
  double spacing = Spacing;
  std::string algorithm = DistanceAlgorithm;
  double alpha = Alpha;
  double tolerance = Tolerance;
  bool is3D = Is3D;
  std::string wasserstein = WassersteinMetric;

  // Transform inputs into the right structure.
  for(int i = 0; i < numInputs; ++i) {
    vtkUnstructuredGrid *grid1 = vtkUnstructuredGrid::New();
    grid1->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(input[i]));
    this->getPersistenceDiagram(inputPersistenceDiagrams[i], grid1, spacing, 0);
  }

  this->performMatchings<dataType>(
    numInputs, inputPersistenceDiagrams, outputMatchings,
    algorithm, // Not from paraview, from enclosing tracking plugin
    wasserstein, tolerance, is3D,
    alpha, // Blending
    PX, PY, PZ, PS, PE // Coefficients
  );

  // Get back meshes.
  //  #pragma omp parallel for num_threads(ThreadNumber)
  for(int i = 0; i < numInputs - 1; ++i) {
    vtkUnstructuredGrid *grid1 = vtkUnstructuredGrid::New();
    grid1->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(input[i]));
    vtkUnstructuredGrid *grid2 = vtkUnstructuredGrid::New();
    grid2->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(input[i + 1]));

    vtkUnstructuredGrid *CTPersistenceDiagram1_
      = vtkUnstructuredGrid::SafeDownCast(grid1);
    vtkUnstructuredGrid *CTPersistenceDiagram2_
      = vtkUnstructuredGrid::SafeDownCast(grid2);

    int status = this->augmentPersistenceDiagrams<dataType>(
      inputPersistenceDiagrams[i], inputPersistenceDiagrams[i + 1],
      outputMatchings[i], CTPersistenceDiagram1_, CTPersistenceDiagram2_);
    if(status < 0)
      return -2;

    outputPersistenceDiagrams[2 * i]->ShallowCopy(CTPersistenceDiagram1_);
    outputPersistenceDiagrams[2 * i + 1]->ShallowCopy(CTPersistenceDiagram2_);
  }

  auto numPersistenceDiagramsInput
    = (int)outputPersistenceDiagrams.size(); // numInputs;

  for(int i = 0; i < numPersistenceDiagramsInput; ++i) {
    vtkSmartPointer<vtkUnstructuredGrid> grid = outputPersistenceDiagrams[i];
    if(!grid || !grid->GetCellData()
       || !grid->GetCellData()->GetArray("Persistence")) {
      this->printErr("Inputs are not persistence diagrams");
      return 0;
    }

    // Check if inputs have the same data type and the same number of points
    if(grid->GetCellData()->GetArray("Persistence")->GetDataType()
       != outputPersistenceDiagrams[0]
            ->GetCellData()
            ->GetArray("Persistence")
            ->GetDataType()) {
      this->printErr("Inputs of different data types");
      return 0;
    }
  }

  for(int i = 0; i < numPersistenceDiagramsInput - 1; ++i) {
    vtkSmartPointer<vtkUnstructuredGrid> grid1 = outputPersistenceDiagrams[i];
    vtkSmartPointer<vtkUnstructuredGrid> grid2
      = outputPersistenceDiagrams[i + 1];
    if(i % 2 == 1 && i < numPersistenceDiagramsInput - 1
       && grid1->GetCellData()->GetNumberOfTuples()
            != grid2->GetCellData()->GetNumberOfTuples()) {
      this->printErr("Inconsistent length or order of input diagrams.");
      return 0;
    }
  }

  outputMesh_ = vtkUnstructuredGrid::New();
  vtkUnstructuredGrid *outputMesh = vtkUnstructuredGrid::SafeDownCast(mesh);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> valueScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> matchingIdScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> lengthScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> timeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> componentIds
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> pointTypeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  persistenceScalars->SetName("Cost");
  valueScalars->SetName("Scalar");
  matchingIdScalars->SetName("MatchingIdentifier");
  lengthScalars->SetName("ComponentLength");
  timeScalars->SetName("TimeStep");
  componentIds->SetName("ConnectedComponentId");
  pointTypeScalars->SetName("CriticalType");

  // (+ vertex id)
  std::vector<trackingTuple>
    trackingsBase; // structure containing all trajectories
  this->performTracking<dataType>(
    inputPersistenceDiagrams, outputMatchings, trackingsBase);

  std::vector<std::set<int>> trackingTupleToMerged(
    trackingsBase.size(), std::set<int>());
  if(DoPostProc)
    this->performPostProcess<dataType>(inputPersistenceDiagrams, trackingsBase,
                                       trackingTupleToMerged, PostProcThresh);

  // bool Is3D = true;
  bool useGeometricSpacing = UseGeometricSpacing;
  // auto spacing = (float) Spacing;

  // std::vector<trackingTuple> trackings = *trackingsBase;
  // Row = iteration number
  // Col = (id in pd 2, arrival point id)

  // Build mesh.
  buildMesh<dataType>(
    trackingsBase, outputMatchings, inputPersistenceDiagrams,
    useGeometricSpacing, spacing, DoPostProc, trackingTupleToMerged, points,
    persistenceDiagram, persistenceScalars, valueScalars, matchingIdScalars,
    lengthScalars, timeScalars, componentIds, pointTypeScalars);

  outputMesh_->ShallowCopy(persistenceDiagram);
  outputMesh->ShallowCopy(outputMesh_);

  return 1;
}
