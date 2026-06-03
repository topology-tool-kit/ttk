#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkLine.h>
#include <vtkPointData.h>

#include <cstdio>
// #include <vtkNew.h>
// #include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkTrackingFromFields.h>
#include <ttkTrackingFromPersistenceDiagrams.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkTrackingFromFields);

ttkTrackingFromFields::ttkTrackingFromFields() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
}

int ttkTrackingFromFields::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
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

int ttkTrackingFromFields::RequestDataObject(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  if(outInfo
     && !vtkUnstructuredGrid::SafeDownCast(
       outInfo->Get(vtkDataObject::DATA_OBJECT()))) {
    vtkNew<vtkUnstructuredGrid> ug;
    outInfo->Set(vtkDataObject::DATA_OBJECT(), ug);
  }

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
  if(!inInfo || !outInfo1)
    return 0;

  vtkDataObject *inputDO = inInfo->Get(vtkDataObject::DATA_OBJECT());
  vtkDataObject *currentDO = outInfo1->Get(vtkDataObject::DATA_OBJECT());

  if(inputDO == nullptr)
    return 0;

  if(currentDO == nullptr || !currentDO->IsA(inputDO->GetClassName())) {
    vtkSmartPointer<vtkDataObject> newDO
      = vtkSmartPointer<vtkDataObject>::Take(inputDO->NewInstance());
    outInfo1->Set(vtkDataObject::DATA_OBJECT(), newDO);
  }

  return 1;
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
  vtkDataSet *segOutput,
  vtkDataSet *input,
  const std::vector<vtkDataArray *> &inputScalarFields,
  const triangulationType *triangulation) {

  ttk::Timer timer;
  this->printMsg(ttk::debug::Separator::L2);

  const bool rebuildMesh = (DoLinearize || DoFusion);
  if(!rebuildMesh && !DoMergeTree) {
    this->printMsg("Nothing to do: all post-processing flags disabled.");
    return 1;
  }

  vtkIntArray *compIdArray = vtkIntArray::SafeDownCast(
    output->GetCellData()->GetArray("ConnectedComponentId"));
  vtkIntArray *timeArray
    = vtkIntArray::SafeDownCast(output->GetPointData()->GetArray("TimeStep"));
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
  std::vector<int> originalCCIds(numTraj, -1);
  std::vector<std::vector<int>> trajTime(numTraj);
  std::vector<std::vector<int>> trajVertexId(numTraj);
  std::vector<std::vector<double>> trajX, trajY;
  std::vector<int> trajCriticalType(numTraj, -1);
  if(rebuildMesh) {
    trajX.assign(numTraj, {});
    trajY.assign(numTraj, {});
  }

  vtkNew<vtkIdList> cellPointIds;
  auto collectUniqueSortedPointIds = [&](const std::vector<vtkIdType> &cellIds,
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
    pointIds.erase(
      std::unique(pointIds.begin(), pointIds.end()), pointIds.end());
    std::sort(pointIds.begin(), pointIds.end(), [&](vtkIdType a, vtkIdType b) {
      return timeArray->GetValue(a) < timeArray->GetValue(b);
    });
  };

  size_t tIdx = 0;
  for(const auto &kv : cellsByTraj) {
    originalCCIds[tIdx] = kv.first;
    std::vector<vtkIdType> pointIds;
    collectUniqueSortedPointIds(kv.second, pointIds);

    auto &ts = trajTime[tIdx];
    auto &vid = trajVertexId[tIdx];
    ts.reserve(pointIds.size());
    vid.reserve(pointIds.size());
    if(rebuildMesh) {
      trajX[tIdx].reserve(pointIds.size());
      trajY[tIdx].reserve(pointIds.size());
    }

    double xyz[3];
    for(const vtkIdType pId : pointIds) {
      ts.push_back(timeArray->GetValue(pId));
      vid.push_back(vertexGlobalIdArray->GetValue(pId));
      if(rebuildMesh) {
        output->GetPoint(pId, xyz);
        trajX[tIdx].push_back(xyz[0]);
        trajY[tIdx].push_back(xyz[1]);
      }
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
  ppt.setDoLinearizeFuse(LinearizeFuse);
  ppt.setDoMergeTree(DoMergeTree);
  ppt.setUseSplitTree(UseSplitTree);

  ppt.setCosCol(std::cos(CosColDegrees * M_PI / 180.0));
  ppt.setMaxRadius(MaxLinkRadius);
  ppt.setMaxFrameDist(MaxFrameDist);
  ppt.setPersistenceThreshold(Tolerance);
  ppt.setMaxSurfSize(MaxSurfSize);
  ppt.setUseOtsuSimplification(UseOtsuSimplification);
  ppt.setOtsuBins(OtsuBins);

  double *bounds = input->GetBounds();
  ppt.setBoundaryXMin(bounds[0]);
  ppt.setBoundaryXMax(bounds[1]);
  ppt.setBoundaryYMin(bounds[2]);
  ppt.setBoundaryYMax(bounds[3]);

  ppt.preconditionTriangulation(const_cast<triangulationType *>(triangulation));

  if(DoMergeTree) {
    std::vector<void *> inputFields;
    inputFields.reserve(inputScalarFields.size());
    for(vtkDataArray *a : inputScalarFields)
      inputFields.push_back(ttkUtils::GetVoidPointer(a));
    ppt.setInputScalars(inputFields);
  }

  std::vector<std::vector<int>> vertexTrajPerFrame;

  // merge-tree segmentation only
  if(!rebuildMesh) {
    std::vector<ttk::PostProcessingTracking::LinearTrajectory> rawTraj;
    rawTraj.reserve(numTraj);
    for(int i = 0; i < numTraj; ++i) {
      if(trajTime[i].empty())
        continue;
      ttk::PostProcessingTracking::LinearTrajectory lt{};
      lt.isLinearized = false;
      lt.startFrame = trajTime[i].front();
      lt.endFrame = trajTime[i].back();
      lt.finalChainId = originalCCIds[i];
      lt.originalTrajId = originalCCIds[i];
      lt.criticalPoints.reserve(trajTime[i].size());
      for(size_t k = 0; k < trajTime[i].size(); ++k)
        lt.criticalPoints.emplace_back(
          trajTime[i][k], static_cast<ttk::SimplexId>(trajVertexId[i][k]));
      rawTraj.push_back(std::move(lt));
    }

    std::vector<double> surfMin, surfMax, surfMean;
    const int mtStatus = ppt.computeMergeTree<dataType, triangulationType>(
      triangulation, rawTraj, surfMin, surfMax, surfMean, vertexTrajPerFrame);
    if(mtStatus < 0) {
      this->printWrn("Merge-tree segmentation failed; "
                     "keeping the raw tracking mesh.");
      return 0;
    }

    writeSegmentationArrays(segOutput, vertexTrajPerFrame);

    this->printMsg("Post-processing (merge-tree only, "
                     + std::to_string(rawTraj.size()) + " trajectories)",
                   1.0, timer.getElapsedTime(), this->threadNumber_);
    this->printMsg(ttk::debug::Separator::L2);
    return 1;
  }

  // full postprocess pipeline (+ optional merge-tree)
  std::vector<ttk::PostProcessingTracking::LinearTrajectory> linearTraj;
  std::vector<ttk::PostProcessingTracking::LinearTrajectory> finalTraj;
  std::vector<ttk::PostProcessingTracking::FuseRecord> fuseRecords;
  std::vector<double> surfMin, surfMax, surfMean;

  const int status = ppt.execute<dataType, triangulationType>(
    trajTime, trajVertexId, trajX, trajY, trajCriticalType, linearTraj,
    finalTraj, fuseRecords, surfMin, surfMax, surfMean, vertexTrajPerFrame,
    triangulation);
  if(status != 1) {
    this->printWrn("Post-processing returned non-success status; "
                   "keeping the raw tracking mesh.");
    return 0;
  }

  const vtkIdType nOut = static_cast<vtkIdType>(finalTraj.size());

  int maxChainId = -1;
  for(const auto &c : finalTraj) {
    if(c.finalChainId > maxChainId)
      maxChainId = c.finalChainId;
  }
  const int nChains = maxChainId + 1;
  std::vector<int> chainCriticalType(std::max(nChains, 0), -1);
  for(size_t i = 0; i < linearTraj.size() && i < trajCriticalType.size(); ++i) {
    const int cid = linearTraj[i].finalChainId;
    if(cid >= 0 && cid < nChains && chainCriticalType[cid] < 0)
      chainCriticalType[cid] = trajCriticalType[i];
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
  auto criticalTypeOut = makeIntArr("CriticalType", nOut);
  auto axArr = makeDblArr("ax", nOut);
  auto bxArr = makeDblArr("bx", nOut);
  auto ayArr = makeDblArr("ay", nOut);
  auto byArr = makeDblArr("by", nOut);
  auto surfMinArr = makeDblArr("SegmentationMin", nOut);
  auto surfMaxArr = makeDblArr("SegmentationMax", nOut);
  auto surfMeanArr = makeDblArr("SegmentationMean", nOut);
  vtkSmartPointer<vtkIntArray> compIdOut;
  if(!LinearizeFuse)
    compIdOut = makeIntArr("ConnectedComponentId", nOut);

  for(vtkIdType i = 0; i < nOut; ++i) {
    const auto &c = finalTraj[i];

    double x0, y0, x1, y1;
    int sF;
    if(DoStartFrame && (DoLinearize || LinearizeFuse)) {
      sF = StartFrame;
    } else {
      sF = c.startFrame;
    }
    const int eF = c.endFrame;
    if(DoLinearize) {
      x0 = c.evalX(sF);
      y0 = c.evalY(sF);
      x1 = c.evalX(eF);
      y1 = c.evalY(eF);
    } else if(!c.criticalPoints.empty()) {
      x0 = c.evalX(sF);
      y0 = c.evalY(sF);
      x1 = c.evalX(eF);
      y1 = c.evalY(eF);
      const ttk::SimplexId v0 = c.criticalPoints.front().second;
      const ttk::SimplexId v1 = c.criticalPoints.back().second;
      if(v0 >= 0 && v0 < triangulation->getNumberOfVertices()) {
        float a, b, cZ;
        triangulation->getVertexPoint(v0, a, b, cZ);
        x0 = a;
        y0 = b;
      }
      if(v1 >= 0 && v1 < triangulation->getNumberOfVertices()) {
        float a, b, cZ;
        triangulation->getVertexPoint(v1, a, b, cZ);
        x1 = a;
        y1 = b;
      }
    } else {
      x0 = c.evalX(sF);
      y0 = c.evalY(sF);
      x1 = c.evalX(eF);
      y1 = c.evalY(eF);
    }

    const vtkIdType p0 = 2 * i + 0;
    const vtkIdType p1 = 2 * i + 1;
    const double spacing = Spacing;

    newPoints->SetPoint(p0, x0, y0, static_cast<double>(sF * spacing));
    newPoints->SetPoint(p1, x1, y1, static_cast<double>(eF * spacing));

    vtkNew<vtkLine> line{};
    line->GetPointIds()->SetId(0, p0);
    line->GetPointIds()->SetId(1, p1);
    newLines->InsertNextCell(line);

    trajIdArr->SetValue(i, c.finalChainId);
    startFrameArr->SetValue(i, c.startFrame);
    endFrameArr->SetValue(i, eF);
    durationArr->SetValue(i, (eF - c.startFrame) + 1);
    {
      const int cid = c.finalChainId;
      const int t = (cid >= 0 && cid < nChains) ? chainCriticalType[cid] : -1;
      criticalTypeOut->SetValue(i, t);
    }
    axArr->SetValue(i, c.ax);
    bxArr->SetValue(i, c.bx);
    ayArr->SetValue(i, c.ay);
    byArr->SetValue(i, c.by);
    surfMinArr->SetValue(i, surfMin[i]);
    surfMaxArr->SetValue(i, surfMax[i]);
    surfMeanArr->SetValue(i, surfMean[i]);
    if(compIdOut)
      compIdOut->SetValue(i, c.originalTrajId);
  }

  newGrid->SetPoints(newPoints);
  newGrid->SetCells(VTK_LINE, newLines);
  newGrid->GetCellData()->AddArray(trajIdArr);
  if(compIdOut)
    newGrid->GetCellData()->AddArray(compIdOut);
  newGrid->GetCellData()->AddArray(startFrameArr);
  newGrid->GetCellData()->AddArray(endFrameArr);
  newGrid->GetCellData()->AddArray(durationArr);
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

  if(DoMergeTree)
    writeSegmentationArrays(segOutput, vertexTrajPerFrame);

  this->printMsg(ttk::debug::Separator::L2);
  return 1;
}

void ttkTrackingFromFields::writeSegmentationArrays(
  vtkDataSet *segOutput,
  const std::vector<std::vector<int>> &vertexTrajPerFrame) {

  if(vertexTrajPerFrame.empty())
    return;

  const vtkIdType nPts = segOutput->GetNumberOfPoints();
  const int nFrames = static_cast<int>(vertexTrajPerFrame.size());

  for(int frame = 0; frame < nFrames; ++frame) {
    const auto &labels = vertexTrajPerFrame[frame];
    if(static_cast<vtkIdType>(labels.size()) != nPts) {
      this->printWrn("Segmentation output size mismatch on frame "
                     + std::to_string(frame));
      continue;
    }
    char segName[20];
    std::snprintf(segName, sizeof(segName), "Seg_%04d", frame);
    vtkNew<vtkIntArray> segArr;
    segArr->SetName(segName);
    segArr->SetNumberOfComponents(1);
    segArr->SetNumberOfTuples(nPts);
    for(vtkIdType v = 0; v < nPts; ++v)
      segArr->SetValue(v, labels[v]);
    segOutput->GetPointData()->AddArray(segArr);
  }
}

int ttkTrackingFromFields::RequestData(vtkInformation *ttkNotUsed(request),
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto segOutput = vtkDataSet::GetData(outputVector, 1);
  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);

  if(input == nullptr || output == nullptr || segOutput == nullptr) {
    return -1;
  }
  segOutput->ShallowCopy(input);
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
    ttkVtkTemplateMacro(inputScalarFields[0]->GetDataType(),
                        triangulation->getType(),
                        (this->applyPostProcessing<VTK_TT, TTK_TT>(
                          output, segOutput, input, inputScalarFields,
                          (TTK_TT *)triangulation->getData())));
  }
  return status;
}
