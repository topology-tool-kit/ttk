#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>
#include <vtkLine.h>

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

  ttk::Timer t{};

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

  ttk::TrackingFromCriticalPoints tracker;
  tracker.setMeshDiameter(meshDiameter);
  tracker.setTolerance(tolerance);
  tracker.setEpsilon(relativeDestructionCost);
  tracker.setAssignmentMethod(assignmentMethod);
  tracker.setWeights(PX, PY, PZ, PF);
  tracker.setUsePersistenceForDistance(UsePersistenceForDistance);
  tracker.setAssignmentPrecision(AssignmentPrecision);

  tracker.setThreadNumber(this->threadNumber_);
  tracker.setDebugLevel(this->debugLevel_);

  std::vector<ttk::DiagramType> persistenceDiagrams(fieldNumber);
  this->performDiagramComputation<dataType, triangulationType>(
    (int)fieldNumber, persistenceDiagrams, triangulation);

  this->printMsg("Diagram computed", 1, t.getElapsedTime(), threadNumber_);
  double previousStepTime = t.getElapsedTime();

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

  this->printMsg("Matchings computed", 1, t.getElapsedTime() - previousStepTime,
                 threadNumber_);
  previousStepTime = t.getElapsedTime();

  vtkNew<vtkPoints> const points{};
  vtkNew<vtkUnstructuredGrid> const outputMesh{};

  vtkNew<vtkDoubleArray> costs{};
  vtkNew<vtkDoubleArray> averagePersistences{};
  vtkNew<vtkDoubleArray> integratedPersistences{};
  vtkNew<vtkDoubleArray> maximalPersistences{};
  vtkNew<vtkDoubleArray> minimalPersistences{};
  vtkNew<vtkDoubleArray> instantPersistences{};
  vtkNew<vtkDoubleArray> valueScalars{};
  vtkNew<vtkIntArray> globalVertexIds{};
  vtkNew<vtkIntArray> lengthScalars{};
  vtkNew<vtkIntArray> timeScalars{};
  vtkNew<vtkIntArray> connectedComponentIds{};
  vtkNew<vtkIntArray> pointsCriticalType{};

  costs->SetName("Costs");
  averagePersistences->SetName("AveragePersistence");
  integratedPersistences->SetName("IntegratedPersistence");
  maximalPersistences->SetName("MaximalPersistence");
  minimalPersistences->SetName("MinimalPersistence");
  instantPersistences->SetName("InstantPersistence");
  valueScalars->SetName("Scalar");
  globalVertexIds->SetName("VertexGlobalId");
  lengthScalars->SetName("ComponentLength");
  timeScalars->SetName("TimeStep");
  connectedComponentIds->SetName("ConnectedComponentId");
  pointsCriticalType->SetName("CriticalType");

  std::vector<ttk::trackingTuple> allTrackings;
  std::vector<std::vector<double>> allTrackingsCosts;
  std::vector<std::vector<double>> allTrackingsInstantPersistence;

  unsigned int typesArrayLimits[3] = {};

  tracker.performTrackings(
    persistenceDiagrams, maximaMatchings, sad_1_Matchings, sad_2_Matchings,
    minimaMatchings, maxMap, sad_1Map, sad_2Map, minMap, allTrackings,
    allTrackingsCosts, allTrackingsInstantPersistence, typesArrayLimits);

  this->printMsg("Trackings computed", 1, t.getElapsedTime() - previousStepTime,
                 threadNumber_);
  previousStepTime = t.getElapsedTime();

  double const spacing = Spacing;
  bool const useGeometricSpacing = UseGeometricSpacing;

  ttkTrackingFromPersistenceDiagrams::buildMeshAlt(
    triangulation, allTrackings, allTrackingsCosts,
    allTrackingsInstantPersistence, useGeometricSpacing, spacing, points,
    outputMesh, pointsCriticalType, timeScalars, lengthScalars, globalVertexIds,
    connectedComponentIds, costs, averagePersistences, integratedPersistences,
    maximalPersistences, minimalPersistences, instantPersistences,
    typesArrayLimits);

  this->printMsg(
    "Mesh built", 1, t.getElapsedTime() - previousStepTime, threadNumber_);
  this->printMsg("Total run time ", 1, t.getElapsedTime(), this->threadNumber_);

  output->ShallowCopy(outputMesh);

  return 1;
}

template <class dataType, class triangulationType>
int ttkTrackingFromFields::applyPostProcessing(
  vtkUnstructuredGrid *output,
  vtkDataSet *input,
  const std::vector<vtkDataArray *> &inputScalarFields,
  const triangulationType *triangulation) {

  ttk::Timer timer;
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Post-processing (linearize="
                 + std::to_string(DoLinearize)
                 + ", fuse=" + std::to_string(DoFusion)
                 + ", mergeTree=" + std::to_string(DoMergeTree) + ")");

  vtkIntArray *compIdArray = vtkIntArray::SafeDownCast(
    output->GetCellData()->GetArray("ConnectedComponentId"));
  vtkIntArray *timeArray = vtkIntArray::SafeDownCast(
    output->GetPointData()->GetArray("TimeStep"));
  vtkIntArray *vertexGlobalIdArray = vtkIntArray::SafeDownCast(
    output->GetPointData()->GetArray("VertexGlobalId"));
  vtkIntArray *criticalTypeArray = vtkIntArray::SafeDownCast(
    output->GetPointData()->GetArray("CriticalType"));

  if(!compIdArray || !timeArray || !vertexGlobalIdArray) {
    this->printErr("Tracking mesh is missing "
                   "ConnectedComponentId/TimeStep/VertexGlobalId; "
                   "skipping post-processing.");
    return 0;
  }

  const vtkIdType numCells = output->GetNumberOfCells();
  std::map<int, std::vector<vtkIdType>> cellsByTraj;
  for(vtkIdType cellId = 0; cellId < numCells; ++cellId)
    cellsByTraj[compIdArray->GetValue(cellId)].push_back(cellId);

  const int numTraj = static_cast<int>(cellsByTraj.size());
  std::vector<std::vector<int>> trajTime(numTraj);
  std::vector<std::vector<int>> trajVertexId(numTraj);
  std::vector<std::vector<double>> trajX(numTraj), trajY(numTraj);
  std::vector<int> trajCriticalType(numTraj, -1);

  vtkNew<vtkIdList> cellPointIds;
  auto collectUniqueSortedPointIds
    = [&](const std::vector<vtkIdType> &cellIds,
          std::vector<vtkIdType> &pointIds) {
        pointIds.clear();
        pointIds.reserve(cellIds.size() * 2);
        for(const vtkIdType cId : cellIds) {
          cellPointIds->Reset();
          output->GetCellPoints(cId, cellPointIds);
          const vtkIdType n = cellPointIds->GetNumberOfIds();
          for(vtkIdType k = 0; k < n; ++k)
            pointIds.push_back(cellPointIds->GetId(k));
        }
        std::sort(pointIds.begin(), pointIds.end());
        pointIds.erase(std::unique(pointIds.begin(), pointIds.end()),
                       pointIds.end());
        std::sort(pointIds.begin(), pointIds.end(),
                  [&](vtkIdType a, vtkIdType b) {
                    return timeArray->GetValue(a) < timeArray->GetValue(b);
                  });
      };

  size_t tIdx = 0;
  for(const auto &kv : cellsByTraj) {
    std::vector<vtkIdType> pointIds;
    collectUniqueSortedPointIds(kv.second, pointIds);

    auto &ts = trajTime[tIdx];
    auto &vid = trajVertexId[tIdx];
    auto &cx = trajX[tIdx];
    auto &cy = trajY[tIdx];
    ts.reserve(pointIds.size());
    vid.reserve(pointIds.size());
    cx.reserve(pointIds.size());
    cy.reserve(pointIds.size());

    double xyz[3];
    for(const vtkIdType pId : pointIds) {
      ts.push_back(timeArray->GetValue(pId));
      vid.push_back(vertexGlobalIdArray->GetValue(pId));
      output->GetPoint(pId, xyz);
      cx.push_back(xyz[0]);
      cy.push_back(xyz[1]);
    }
    if(criticalTypeArray && !pointIds.empty())
      trajCriticalType[tIdx] = criticalTypeArray->GetValue(pointIds.front());
    ++tIdx;
  }

  ttk::PostProcessingTracking ppt;
  ppt.setThreadNumber(this->threadNumber_);
  ppt.setDebugLevel(this->debugLevel_);

  ppt.setDoLinearize(DoLinearize);
  ppt.setDoFusion(DoFusion);
  ppt.setDoMergeTree(DoMergeTree);

  const double pi = M_PI;
  ppt.setCosCol(std::cos(CosColDegrees * pi / 180.0));
  ppt.setMaxRadius(MaxLinkRadius);
  ppt.setMaxFrameDist(MaxFrameDist);
  ppt.setPersistenceThreshold(Tolerance);
  ppt.setMaxSurfSize(MaxSurfSize);
  ppt.setUseOtsuSimplification(UseOtsuSimplification);
  ppt.setOtsuBins(OtsuBins);

  double *bounds = input->GetBounds();
  ppt.setBoundaryXMin(bounds[0]); ppt.setBoundaryXMax(bounds[1]); ppt.setBoundaryYMin(bounds[2]); ppt.setBoundaryYMax(bounds[3]);

  ppt.preconditionTriangulation(
    const_cast<triangulationType *>(triangulation));

  if(DoMergeTree) {
    std::vector<void *> inputFields;
    inputFields.reserve(inputScalarFields.size());
    for(vtkDataArray *a : inputScalarFields)
      inputFields.push_back(ttkUtils::GetVoidPointer(a));
    ppt.setInputScalars(inputFields);
  }

  std::vector<ttk::PostProcessingTracking::LinearTrajectory> linearTraj;
  std::vector<ttk::PostProcessingTracking::LinearTrajectory> finalTraj;
  std::vector<ttk::PostProcessingTracking::FuseRecord> fuseRecords;
  std::vector<double> surfMin, surfMax, surfMean;

  const int status = ppt.execute<dataType, triangulationType>(
    trajTime, trajVertexId, trajX, trajY, trajCriticalType,
    linearTraj, finalTraj, fuseRecords,
    surfMin, surfMax, surfMean, triangulation);
  if(status != 1) {
    this->printWrn("Post-processing returned non-success status; "
                   "keeping the raw tracking mesh.");
    return 0;
  }

  const vtkIdType nOut = static_cast<vtkIdType>(finalTraj.size());

  std::vector<int> finalCriticalType(nOut, -1);
  for(size_t i = 0; i < linearTraj.size() && i < trajCriticalType.size(); ++i) {
    const int cid = linearTraj[i].finalChainId;
    if(cid >= 0 && cid < (int)nOut && finalCriticalType[cid] < 0)
      finalCriticalType[cid] = trajCriticalType[i];
  }

  vtkNew<vtkUnstructuredGrid> newGrid{};
  vtkNew<vtkPoints> newPoints{};
  vtkNew<vtkCellArray> newLines{};
  newPoints->SetNumberOfPoints(2 * nOut);

  auto makeIntArr = [](const char *name, vtkIdType n) {
    auto a = vtkSmartPointer<vtkIntArray>::New();
    a->SetName(name);
    a->SetNumberOfTuples(n);
    return a;
  };
  auto makeDblArr = [](const char *name, vtkIdType n) {
    auto a = vtkSmartPointer<vtkDoubleArray>::New();
    a->SetName(name);
    a->SetNumberOfTuples(n);
    return a;
  };

  auto trajIdArr = makeIntArr("TrajId", nOut);
  auto startFrameArr = makeIntArr("StartFrame", nOut);
  auto endFrameArr = makeIntArr("EndFrame", nOut);
  auto durationArr = makeIntArr("Duration", nOut);
  auto lengthArr = makeIntArr("ComponentLength", nOut);
  auto criticalTypeOut = makeIntArr("CriticalType", nOut);
  auto axArr = makeDblArr("ax", nOut);
  auto bxArr = makeDblArr("bx", nOut);
  auto ayArr = makeDblArr("ay", nOut);
  auto byArr = makeDblArr("by", nOut);
  auto surfMinArr = makeDblArr("SurfaceMin", nOut);
  auto surfMaxArr = makeDblArr("SurfaceMax", nOut);
  auto surfMeanArr = makeDblArr("SurfaceMean", nOut);
  auto compIdOut = makeIntArr("ConnectedComponentId", nOut);

  for(vtkIdType i = 0; i < nOut; ++i) {
    const auto &c = finalTraj[i];

    double x0, y0, x1, y1;
    const int sF = c.startFrame;
    const int eF = c.endFrame;
    if(DoLinearize) {
      x0 = c.evalX(sF); y0 = c.evalY(sF); x1 = c.evalX(eF); y1 = c.evalY(eF);
    } else if(!c.criticalPoints.empty()) {
      x0 = c.evalX(sF); y0 = c.evalY(sF); x1 = c.evalX(eF); y1 = c.evalY(eF);
      const ttk::SimplexId v0 = c.criticalPoints.front().second;
      const ttk::SimplexId v1 = c.criticalPoints.back().second;
      if(v0 >= 0 && v0 < triangulation->getNumberOfVertices()) {
        float a, b, cZ;
        triangulation->getVertexPoint(v0, a, b, cZ);
        x0 = a; y0 = b;
      }
      if(v1 >= 0 && v1 < triangulation->getNumberOfVertices()) {
        float a, b, cZ;
        triangulation->getVertexPoint(v1, a, b, cZ);
        x1 = a; y1 = b;
      }
    } else {
      x0 = y0 = x1 = y1 = 0.0;
    }

    const vtkIdType p0 = 2 * i + 0;
    const vtkIdType p1 = 2 * i + 1;
    double const spacing = Spacing;
    newPoints->SetPoint(p0, x0, y0, static_cast<double>(sF*spacing));
    newPoints->SetPoint(p1, x1, y1, static_cast<double>(eF*spacing));

    vtkNew<vtkLine> line{};
    line->GetPointIds()->SetId(0, p0);
    line->GetPointIds()->SetId(1, p1);
    newLines->InsertNextCell(line);

    trajIdArr->SetValue(i, static_cast<int>(i));
    startFrameArr->SetValue(i, sF);
    endFrameArr->SetValue(i, eF);
    durationArr->SetValue(i, eF - sF);
    lengthArr->SetValue(i, static_cast<int>(c.criticalPoints.size()));
    criticalTypeOut->SetValue(i, finalCriticalType[i]);
    axArr->SetValue(i, c.ax);
    bxArr->SetValue(i, c.bx);
    ayArr->SetValue(i, c.ay);
    byArr->SetValue(i, c.by);
    surfMinArr->SetValue(i, surfMin[i]);
    surfMaxArr->SetValue(i, surfMax[i]);
    surfMeanArr->SetValue(i, surfMean[i]);
    compIdOut->SetValue(i, static_cast<int>(i));
  }

  newGrid->SetPoints(newPoints);
  newGrid->SetCells(VTK_LINE, newLines);
  newGrid->GetCellData()->AddArray(trajIdArr);
  newGrid->GetCellData()->AddArray(compIdOut);
  newGrid->GetCellData()->AddArray(startFrameArr);
  newGrid->GetCellData()->AddArray(endFrameArr);
  newGrid->GetCellData()->AddArray(durationArr);
  newGrid->GetCellData()->AddArray(lengthArr);
  newGrid->GetCellData()->AddArray(criticalTypeOut);
  newGrid->GetCellData()->AddArray(axArr);
  newGrid->GetCellData()->AddArray(bxArr);
  newGrid->GetCellData()->AddArray(ayArr);
  newGrid->GetCellData()->AddArray(byArr);
  if(DoMergeTree) {
    newGrid->GetCellData()->AddArray(surfMinArr);
    newGrid->GetCellData()->AddArray(surfMaxArr);
    newGrid->GetCellData()->AddArray(surfMeanArr);
  }

  output->ShallowCopy(newGrid);

  this->printMsg("Post-processing ("
                   + std::to_string(nOut) + " output trajectories)",
                 1.0, timer.getElapsedTime(), this->threadNumber_);
  this->printMsg(ttk::debug::Separator::L2);
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
  this->printMsg("Tracking trajectories over " + std::to_string(fieldNumber)
                 + " timesteps");
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

  if(status == 1 && EnablePostProc) {
    ttkVtkTemplateMacro(
      inputScalarFields[0]->GetDataType(), triangulation->getType(),
      (this->applyPostProcessing<VTK_TT, TTK_TT>(
        output, input, inputScalarFields,
        (TTK_TT *)triangulation->getData())));
  }
  return status;
}
