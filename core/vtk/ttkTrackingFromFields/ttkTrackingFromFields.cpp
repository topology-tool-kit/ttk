#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkTrackingFromFields.h>
#include <ttkTrackingFromPersistenceDiagrams.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkTrackingFromFields);

ttkTrackingFromFields::ttkTrackingFromFields() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkTrackingFromFields::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}
int ttkTrackingFromFields::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

// (*) Persistence-driven approach
template <class dataType, class triangulationType>
int ttkTrackingFromFields::trackWithPersistenceMatching(
  vtkUnstructuredGrid *output,
  unsigned long fieldNumber,
  const triangulationType *triangulation) {

  std::vector<ttk::DiagramType> persistenceDiagrams(fieldNumber);

  this->performDiagramComputation<dataType, triangulationType>(
    (int)fieldNumber, persistenceDiagrams, triangulation);

  std::vector<std::vector<ttk::MatchingType>> outputMatchings(fieldNumber - 1);

  double const spacing = Spacing;
  std::string const algorithm = DistanceAlgorithm;
  double const tolerance = Tolerance;
  std::string const wasserstein = WassersteinMetric;

  ttk::TrackingFromPersistenceDiagrams tfp{};
  tfp.setThreadNumber(this->threadNumber_);
  tfp.setDebugLevel(this->debugLevel_);
  tfp.performMatchings(
    (int)fieldNumber, persistenceDiagrams, outputMatchings,
    algorithm, // Not from paraview, from enclosing tracking plugin
    wasserstein, tolerance, PX, PY, PZ, PS, PE // Coefficients
  );

  vtkNew<vtkPoints> const points{};
  vtkNew<vtkUnstructuredGrid> const persistenceDiagram{};

  vtkNew<vtkDoubleArray> persistenceScalars{};
  vtkNew<vtkDoubleArray> valueScalars{};
  vtkNew<vtkIntArray> matchingIdScalars{};
  vtkNew<vtkIntArray> lengthScalars{};
  vtkNew<vtkIntArray> timeScalars{};
  vtkNew<vtkIntArray> componentIds{};
  vtkNew<vtkIntArray> pointTypeScalars{};

  persistenceScalars->SetName("Cost");
  valueScalars->SetName("Scalar");
  matchingIdScalars->SetName("MatchingIdentifier");
  lengthScalars->SetName("ComponentLength");
  timeScalars->SetName("TimeStep");
  componentIds->SetName("ConnectedComponentId");
  pointTypeScalars->SetName("CriticalType");

  std::vector<ttk::trackingTuple> trackingsBase;
  tfp.performTracking(persistenceDiagrams, outputMatchings, trackingsBase);

  std::vector<std::set<int>> trackingTupleToMerged(trackingsBase.size());

  if(DoPostProc) {
    tfp.performPostProcess(persistenceDiagrams, trackingsBase,
                           trackingTupleToMerged, PostProcThresh);
  }

  bool const useGeometricSpacing = UseGeometricSpacing;

  // Build mesh.
  ttkTrackingFromPersistenceDiagrams::buildMesh(
    trackingsBase, outputMatchings, persistenceDiagrams, useGeometricSpacing,
    spacing, DoPostProc, trackingTupleToMerged, points, persistenceDiagram,
    persistenceScalars, valueScalars, matchingIdScalars, lengthScalars,
    timeScalars, componentIds, pointTypeScalars, *this);

  output->ShallowCopy(persistenceDiagram);

  return 1;
}

template <class dataType, class triangulationType>
int ttkTrackingFromFields::trackWithCriticalPointMatching(
  vtkUnstructuredGrid *output,
  unsigned long fieldNumber,
  const triangulationType *triangulation) {

  float x = 0, y = 0, z = 0;
  float maxX = 0, minX = 0, maxY = 0, minY = 0, maxZ = 0, minZ = 0;
  triangulation->getVertexPoint(0, minX, minY, minZ);
  triangulation->getVertexPoint(0, maxX, maxY, maxZ);

  for(int i = 0; i < triangulation->getNumberOfVertices(); i++) {
    triangulation->getVertexPoint(i, x, y, z);
    maxX = std::max(x, maxX);
    maxX = std::min(x, minX);
    maxY = std::max(y, maxX);
    minY = std::min(y, minY);
    maxZ = std::max(z, maxZ);
    minZ = std::min(z, minZ);
  }

  double const relativeDestructionCost = RelativeDestructionCost;
  double const tolerance = (double)Tolerance;
  float meshDiameter
    = std::sqrt(std::pow(maxX - minX, 2) + std::pow(maxY - minY, 2)
                + std::pow(maxZ - minZ, 2));
  int assignmentMethod = AssignmentMethod;
  bool adaptDeathBirthCost = AdaptDeathBirthCost;
  double epsilonAdapt = EpsilonAdapt;

  ttk::TrackingFromCriticalPoints tracker;
  tracker.setMeshDiameter(meshDiameter);
  tracker.setTolerance(tolerance);
  tracker.setEpsilon(relativeDestructionCost);
  tracker.setAdaptDeathBirthCost(adaptDeathBirthCost);
  tracker.setAssignmentMethod(assignmentMethod);
  tracker.setEpsilonAdapt(epsilonAdapt);
  tracker.setWeights(PX, PY, PZ, PF);

  tracker.setThreadNumber(this->threadNumber_);
  tracker.setDebugLevel(this->debugLevel_);

  std::vector<ttk::DiagramType> persistenceDiagrams(fieldNumber);
  this->performDiagramComputation<dataType, triangulationType>(
    (int)fieldNumber, persistenceDiagrams, triangulation);

  std::vector<std::vector<ttk::MatchingType>> maximaMatchings(fieldNumber - 1);
  std::vector<std::vector<ttk::MatchingType>> sad_1_Matchings(fieldNumber - 1);
  std::vector<std::vector<ttk::MatchingType>> sad_2_Matchings(fieldNumber - 1);
  std::vector<std::vector<ttk::MatchingType>> minimaMatchings(fieldNumber - 1);

  std::vector<std::vector<ttk::SimplexId>> maxMap(fieldNumber);
  std::vector<std::vector<ttk::SimplexId>> sad_1Map(fieldNumber);
  std::vector<std::vector<ttk::SimplexId>> sad_2Map(fieldNumber);
  std::vector<std::vector<ttk::SimplexId>> minMap(fieldNumber);

  tracker.performMatchings(persistenceDiagrams, maximaMatchings,
                           sad_1_Matchings, sad_2_Matchings, minimaMatchings,
                           maxMap, sad_1Map, sad_2Map, minMap);

  vtkNew<vtkPoints> const points{};
  vtkNew<vtkUnstructuredGrid> const outputMesh{};

  vtkNew<vtkDoubleArray> costs{};
  vtkNew<vtkDoubleArray> averagePersistences{};
  vtkNew<vtkDoubleArray> valueScalars{};
  vtkNew<vtkIntArray> globalVertexIds{};
  vtkNew<vtkIntArray> lengthScalars{};
  vtkNew<vtkIntArray> timeScalars{};
  vtkNew<vtkIntArray> connectedComponentIds{};
  vtkNew<vtkIntArray> pointsCriticalType{};

  costs->SetName("Costs");
  averagePersistences->SetName("AveragePersistence");
  valueScalars->SetName("Scalar");
  globalVertexIds->SetName("VertexGlobalId");
  lengthScalars->SetName("ComponentLength");
  timeScalars->SetName("TimeStep");
  connectedComponentIds->SetName("ConnectedComponentId");
  pointsCriticalType->SetName("CriticalType");

  std::vector<ttk::trackingTuple> allTrackings;
  std::vector<double> allTrackingsMeanPersistence;
  std::vector<std::vector<double>> allTrackingsCosts;

  unsigned int typesArrayLimits[3] = {};

  tracker.performTrackings(
    persistenceDiagrams, maximaMatchings, sad_1_Matchings, sad_2_Matchings,
    minimaMatchings, maxMap, sad_1Map, sad_2Map, minMap, allTrackings,
    allTrackingsCosts, allTrackingsMeanPersistence, typesArrayLimits);

  double const spacing = Spacing;
  bool const useGeometricSpacing = UseGeometricSpacing;

  ttkTrackingFromPersistenceDiagrams::buildMesh(
    triangulation, allTrackings, allTrackingsCosts, allTrackingsMeanPersistence,
    useGeometricSpacing, spacing, points, outputMesh, pointsCriticalType,
    timeScalars, lengthScalars, globalVertexIds, connectedComponentIds, costs,
    averagePersistences, typesArrayLimits);

  output->ShallowCopy(outputMesh);

  return 1;
}

int ttkTrackingFromFields::RequestData(vtkInformation *ttkNotUsed(request),
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);
  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);

  // Test validity of datasets
  if(input == nullptr || output == nullptr) {
    return -1;
  }

  // Get number and list of inputs.
  std::vector<vtkDataArray *> inputScalarFieldsRaw;
  std::vector<vtkDataArray *> inputScalarFields;
  const auto pointData = input->GetPointData();
  int numberOfInputFields = pointData->GetNumberOfArrays();
  if(numberOfInputFields < 3) {
    this->printErr("Not enough input fields to perform tracking.");
  }

  vtkDataArray *firstScalarField = pointData->GetArray(0);

  for(int i = 0; i < numberOfInputFields; ++i) {
    vtkDataArray *currentScalarField = pointData->GetArray(i);
    if(currentScalarField == nullptr
       || currentScalarField->GetName() == nullptr) {
      continue;
    }
    std::string const sfname{currentScalarField->GetName()};
    if(sfname.rfind("_Order") == (sfname.size() - 6)) {
      continue;
    }
    if(firstScalarField->GetDataType() != currentScalarField->GetDataType()) {
      this->printErr("Inconsistent field data type or size between fields `"
                     + std::string{firstScalarField->GetName()} + "' and `"
                     + sfname + "'");
      return -1;
    }
    inputScalarFieldsRaw.push_back(currentScalarField);
  }

  std::sort(inputScalarFieldsRaw.begin(), inputScalarFieldsRaw.end(),
            [](vtkDataArray *a, vtkDataArray *b) {
              std::string s1 = a->GetName();
              std::string s2 = b->GetName();
              return std::lexicographical_compare(
                s1.begin(), s1.end(), s2.begin(), s2.end());
            });

  numberOfInputFields = inputScalarFieldsRaw.size();
  int const end = EndTimestep <= 0 ? numberOfInputFields
                                   : std::min(numberOfInputFields, EndTimestep);
  for(int i = StartTimestep; i < end; i += Sampling) {
    vtkDataArray *currentScalarField = inputScalarFieldsRaw[i];
    // Print scalar field names:
    // std::cout << currentScalarField->GetName() << std::endl;
    inputScalarFields.push_back(currentScalarField);
  }

  // Input -> persistence filter.
  std::string const algorithm = DistanceAlgorithm;
  int const pvalg = PVAlgorithm;
  bool useTTKMethod = false;
  bool trackWithCriticalPoints = (pvalg == 2);

  if(pvalg >= 0) {
    switch(pvalg) {
      case 0:
      case 1:
      case 2:
      case 3:
        useTTKMethod = true;
        break;
      case 4:
        break;
      default:
        this->printMsg("Unrecognized tracking method.");
        break;
    }
  } else {
    using ttk::str2int;
    switch(str2int(algorithm.c_str())) {
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
        this->printMsg("Unrecognized tracking method.");
        break;
    }
  }

  // 0. get data
  int const fieldNumber = inputScalarFields.size();
  std::vector<void *> inputFields(fieldNumber);
  for(int i = 0; i < fieldNumber; i++) {
    inputFields[i] = ttkUtils::GetVoidPointer(inputScalarFields[i]);
  }
  this->setInputScalars(inputFields);

  // 0'. get offsets
  std::vector<ttk::SimplexId *> inputOrders(fieldNumber);
  for(int i = 0; i < fieldNumber; ++i) {
    this->SetInputArrayToProcess(0, 0, 0, 0, inputScalarFields[i]->GetName());
    auto orderArray
      = this->GetOrderArray(input, 0, triangulation, false, 0, false);
    inputOrders[i]
      = static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(orderArray));
  }
  this->setInputOffsets(inputOrders);

  int status = 0;
  if(useTTKMethod && !trackWithCriticalPoints) {
    ttkVtkTemplateMacro(
      inputScalarFields[0]->GetDataType(), triangulation->getType(),
      (status = this->trackWithPersistenceMatching<VTK_TT, TTK_TT>(
         output, fieldNumber, (TTK_TT *)triangulation->getData())));
  } else if(useTTKMethod && trackWithCriticalPoints) {
    ttkVtkTemplateMacro(
      inputScalarFields[0]->GetDataType(), triangulation->getType(),
      (status = this->trackWithCriticalPointMatching<VTK_TT, TTK_TT>(
         output, fieldNumber, (TTK_TT *)triangulation->getData())));
  } else {
    this->printMsg("The specified matching method is not supported.");
  }

  return status;
}
