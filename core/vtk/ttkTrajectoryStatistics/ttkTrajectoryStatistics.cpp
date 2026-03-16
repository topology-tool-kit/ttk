
#include <ttkTrajectoryStatistics.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCharArray.h>
#include <vtkTable.h> 
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkGradientFilter.h>
#include <vtkLine.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <algorithm> // std::sort, std::unique
#include <vector> 
#include <map>

// A VTK macro that enables the instantiation of this class via ::New()
vtkStandardNewMacro(ttkTrajectoryStatistics);


ttkTrajectoryStatistics::ttkTrajectoryStatistics() {
  this->setDebugMsgPrefix("TrajectoryStatistics");
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(4);
}


int ttkTrajectoryStatistics::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTrajectoryStatistics::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  if (port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 1);
    return 1;
  }
  if (port == 3) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  /*
  if (port == 5) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  */
  return 0;
}


int ttkTrajectoryStatistics::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  vtkUnstructuredGrid *inputGrid = vtkUnstructuredGrid::GetData(inputVector[0]);
  if (!inputGrid) {
    this->printErr("Invalid input : vtkUnstructuredGrid missing");
    return 0;
  }
  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[1]);
  if (!inputDataSet) {
    this->printErr("Invalid input : vtkDataSet missing");
    return 0;
  }

    
  //trackingFromFields data
  
  vtkIntArray *compIdArray = vtkIntArray::SafeDownCast(
            inputGrid->GetCellData()->GetArray("ConnectedComponentId"));
  vtkIntArray *timeArray = vtkIntArray::SafeDownCast(
            inputGrid->GetPointData()->GetArray("TimeStep"));
  vtkIntArray *vertexGlobalIdArray = vtkIntArray::SafeDownCast(
            inputGrid->GetPointData()->GetArray("VertexGlobalId"));
  vtkIntArray *compLength = vtkIntArray::SafeDownCast(
            inputGrid->GetCellData()->GetArray("ComponentLength"));
  vtkDoubleArray *persistanceArray = vtkDoubleArray::SafeDownCast(
            inputGrid->GetPointData()->GetArray("InstantPersistence"));

  if (!compIdArray || !timeArray || !vertexGlobalIdArray || !compLength){
    this->printErr("Missing data in input vtu");
    return 0;
  }

  const vtkIdType numCells = inputGrid->GetNumberOfCells();
  std::map<int, std::vector<vtkIdType>> cellsByTraj;

  for(vtkIdType cellId = 0; cellId < numCells; ++cellId) {
    const int trajId = compIdArray->GetValue(cellId);
    cellsByTraj[trajId].push_back(cellId);
  }
  
  const int numTraj = static_cast<int>(cellsByTraj.size());
  
  // Outputs per trajectory
  std::vector<std::vector<int>>           trajTime(numTraj);
  std::vector<std::vector<vtkIdType>>     cellsPerTraj(numTraj);
  std::vector<std::vector<double>>        trajX(numTraj), trajY(numTraj), instantPersistance(numTraj);
  std::vector<std::vector<int>>           trajVertexId(numTraj); // VertexGlobalId
  
  vtkNew<vtkIdList> cellPointIds;
  
  auto collectUniqueSortedPointIds = [&](const std::vector<vtkIdType> &cellIds,
                                         std::vector<vtkIdType> &pointIds) {
    pointIds.clear();
    pointIds.reserve(cellIds.size() * 2);
    for(const vtkIdType cId : cellIds) {
      cellPointIds->Reset();
      inputGrid->GetCellPoints(cId, cellPointIds);
      const vtkIdType n = cellPointIds->GetNumberOfIds();
      for(vtkIdType k = 0; k < n; ++k)
        pointIds.push_back(cellPointIds->GetId(k));
    }
    std::sort(pointIds.begin(), pointIds.end());
    pointIds.erase(std::unique(pointIds.begin(), pointIds.end()), pointIds.end());
  
    std::sort(pointIds.begin(), pointIds.end(),
              [&](vtkIdType a, vtkIdType b) { return timeArray->GetValue(a) < timeArray->GetValue(b); });
  };
  
  size_t tIdx = 0;
  for(const auto &kv : cellsByTraj) {
    const auto &cellIds = kv.second;
    cellsPerTraj[tIdx] = cellIds;
  
    std::vector<vtkIdType> pointIds;
    collectUniqueSortedPointIds(cellIds, pointIds);
  
    auto &TT   = trajTime[tIdx];
    auto &VID  = trajVertexId[tIdx];
    auto &IP   = instantPersistance[tIdx];
    auto &XX   = trajX[tIdx];
    auto &YY   = trajY[tIdx];
  
    TT.reserve(pointIds.size());
    VID.reserve(pointIds.size());
    IP.reserve(pointIds.size());
    XX.reserve(pointIds.size());
    YY.reserve(pointIds.size());
  
    double xyz[3];
    for(const vtkIdType pId : pointIds) {
      const int    t   = timeArray->GetValue(pId);
      const double per = persistanceArray->GetValue(pId);
      const int    gid = vertexGlobalIdArray->GetValue(pId);
  
      TT.push_back(t);
      VID.push_back(gid);
      IP.push_back(per);
  
      inputGrid->GetPoint(pId, xyz);
      XX.push_back(xyz[0]);
      YY.push_back(xyz[1]);
    }
  
    ++tIdx;
  }
  
  //scalar dataset
    
  vtkPointData *pd = inputDataSet->GetPointData();
  if(!pd) {
    this->printErr("scalarArray missing");
    return 0;
  }
  
  const int nArrays = pd->GetNumberOfArrays();
  vtkDataArray *ref = pd->GetArray(0);  
  std::vector<vtkDataArray*> fields;
  fields.reserve(nArrays);
  
  for(int i = 0; i < nArrays; ++i) {
    vtkDataArray *a = pd->GetArray(i);
    if(!a || !a->GetName()) continue;
  
    const std::string name{a->GetName()};
    // skip arrays ending with "_Order"
    if(name.size() >= 6 && name.compare(name.size() - 6, 6, "_Order") == 0)
      continue;
  
    if(ref->GetDataType() != a->GetDataType()) {
      this->printErr("Inconsistent field data type or size between fields `"
                     + std::string{ref->GetName()} + "' and `" + name + "'");
      return -1;
    }
  
    fields.push_back(a);
  }
  
  std::sort(fields.begin(), fields.end(),
            [](vtkDataArray *a, vtkDataArray *b) {
              const std::string s1 = a->GetName();
              const std::string s2 = b->GetName();
              return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end());
            });
  
  std::vector<void*> inputFields(fields.size());
  for(size_t i = 0; i < fields.size(); ++i)
    inputFields[i] = ttkUtils::GetVoidPointer(fields[i]);
  
  
  this->setInputScalars(inputFields);
  this->setInstantPersistence(instantPersistance);
  this->setFiltreY(filtreY);
  this->setCosCol(cosCol);
  this->setMaxRadius(maxRadus);
  this->setMaxFrameDist(maxFrameDist);
  this->setSpatialScale(1/(spatialScale*1000));
  this->setInterFrame(interFrame*std::pow(10.0,-6.0));
  this->setConvertDur(convertDur);
  this->setMinVx(minVx);
  this->setMaxVx(maxVx);
  this->setEnableFilteringMinVx(enableFilteringMinVx);
  this->setEnableFilteringTimeOrigin(enableFilteringTimeOrigin);
  this->setEnableFilteringCosY(enableFilteringCosY);
  this->SetEnableFilteringDuration(enableFilteringDuration);
  this->setDuraMin(duraMin);
  this->setXOrigin(xOrigin);
  this->setMinTimeOrigin(minTimeOrigin);
  this->setMinYTimeOrigin(minYTimeOrigin);
  this->setMaxYTimeOrigin(maxYTimeOrigin);
  this->setMaxX(maxX);
  this->setMaxY(maxY);
  this->setMinY(minY);
  this->setMinX(minX);
  this->setSurfaceMethod(surfaceMethod);
  this->setPersisThresh(persisThresh);
  std::vector<ttk::SimplexId> minSeg(numTraj);
  std::vector<ttk::SimplexId> saddleSeg(numTraj);
  this->setMinSeg(minSeg);
  this->setSaddleSeg(saddleSeg);
  this->setErrSurf(errSurf);
  this->setOnlyFrameSurface(onlyFrameSurface);
  this->setMaxSurfSize(maxSurfSize);
  double *bounds = inputDataSet->GetBounds();
  this->setBoundaryX(bounds[1]);
  this->setBoundaryXMin_(bounds[0]);
  this->setBoundaryYMin_(bounds[2]);
  this->setBoundaryY(bounds[3]);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation) return 0;
  this->preconditionTriangulation(triangulation);

  std::vector<ttk::TrajectoryStatistics::LinearTrajectory> linearTraj(numTraj);
  std::vector<ttk::TrajectoryStatistics::LinearTrajectory> finalTraj;
  std::vector<ttk::TrajectoryStatistics::FuseRecord> fuseRecords;

  this->correctTrajectory(trajTime, trajVertexId, trajX, trajY, finalTraj,linearTraj, fuseRecords); 
  std::vector<int>  durations(numTraj);
  std::vector<double> VX(numTraj), VY(numTraj), 
                      surfMin(numTraj), surfMax(numTraj), surfMean(numTraj);
  std::vector<std::vector<ttk::SimplexId>> allVertexDebris(inputFields.size());
  for (size_t frame = 0; frame < allVertexDebris.size(); frame ++){
  	allVertexDebris[frame].assign(triangulation->getNumberOfVertices(), -1);
  }
  this->printMsg("taille = "+ std::to_string(allVertexDebris.size()) + " " + std::to_string(allVertexDebris[0].size()));


  std::vector<std::vector<double>> gradientNorms;

  if(!computeAllGradientMagnitudes(inputDataSet,
                                   fields,
                                   gradientNorms)) {
    this->printErr("Gradient Magnitudes fails");
    return 0;
  }


  int status = 0;
  ttkVtkTemplateMacro(fields[0]->GetDataType(), triangulation->getType(),
      (status = this->execute<VTK_TT, TTK_TT>(
                        trajTime, 
                        trajVertexId,
                        durations,
                        VX,
                        VY,
                        surfMin, 
                        surfMax,
                        surfMean,
                        allVertexDebris,
                        frameSurface,
                        gradientNorms,
                        finalTraj,
                        (TTK_TT *)triangulation->getData()
                        )));
  if (status != 1) return 0;

  // OUTPUT 
  
  vtkTable *outputTable = vtkTable::GetData(outputVector, 0);
  if(!outputTable) { this->printErr("output vtkTable"); return 0; }
  
  vtkUnstructuredGrid *outputTraj = vtkUnstructuredGrid::GetData(outputVector, 1);
  if(!outputTraj) { this->printErr("Null output Trajectories."); return 0; }

  vtkDataSet *outputSurface = vtkDataSet::GetData(outputVector, 2);
  if(!outputSurface) { this->printErr("Null output Surfaces."); return 0; }
  
  vtkUnstructuredGrid *outputLinear = vtkUnstructuredGrid::GetData(outputVector, 3);
  if(!outputLinear) { this->printErr("Null output Linear."); return 0; }
  
/*
  vtkUnstructuredGrid *outputInitial = vtkUnstructuredGrid::GetData(outputVector, 5);
  if(!outputSegmentsLabeled) { this->printErr("Null output "); return 0; }
*/
  auto makeIntCol = [&](const char *name, vtkIdType n) {
    auto arr = vtkSmartPointer<vtkIntArray>::New();
    arr->SetName(name);
    arr->SetNumberOfTuples(n);
    return arr;
  };

  auto makeDblCol = [&](const char *name, vtkIdType n) {
    auto arr = vtkSmartPointer<vtkDoubleArray>::New();
    arr->SetName(name);
    arr->SetNumberOfTuples(n);
    return arr;
  };

  auto evalX = [](const std::vector<double> &c, int t) { return c[0] * t + c[2]; };
  auto evalY = [](const std::vector<double> &c, int t) { return c[1] * t + c[3]; };

  auto addSegment = [&](vtkPoints *pts, vtkCellArray *lines, vtkIdType segIdx,
                        double x0, double y0, double t0,
                        double x1, double y1, double t1) {
    const vtkIdType p0 = 2 * segIdx + 0;
    const vtkIdType p1 = 2 * segIdx + 1;

    pts->SetPoint(p0, x0, y0, t0);
    pts->SetPoint(p1, x1, y1, t1);

    auto line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, p0);
    line->GetPointIds()->SetId(1, p1);
    lines->InsertNextCell(line);
  };
  // ---------------------- VTK TABLE -------------------------------

  const int numMerge = static_cast<int>(finalTraj.size());
  
  auto colStartF = makeIntCol("StartFrame", numMerge);
  auto colEndF   = makeIntCol("EndFrame",   numMerge);
  auto colDur    = makeIntCol("Duration",   numMerge);
  auto colVX     = makeDblCol("VX",         numMerge);
  auto colVY     = makeDblCol("VY",         numMerge);
  auto colSurfMin  = makeDblCol("SurfaceMin",  numMerge);
  auto colSurfMax  = makeDblCol("SurfaceMax",  numMerge);
  auto colSurfMean = makeDblCol("SurfaceMean", numMerge);
  auto colVolMean = makeDblCol("VolumeMean", numMerge);
  auto colTrajId = makeIntCol("TrajId", numMerge);
  double scale_pixel_to_meter = 1/(spatialScale*1000);

  for(int i = 0; i < numMerge; ++i) {
    colStartF->SetValue(i, finalTraj[i].startFrame);
    colEndF  ->SetValue(i, finalTraj[i].endFrame);
    colDur   ->SetValue(i, durations[i]); 
    colVX    ->SetValue(i, VX[i]);
    colVY    ->SetValue(i, VY[i]);
	colSurfMin -> SetValue(i, surfMin[i]);
	colSurfMax -> SetValue(i, surfMax[i]);
	colSurfMean -> SetValue(i, surfMean[i]);
	colTrajId -> SetValue(i, i);

	double vol = 0.0;
	if (surfMean[i] > 0.0 && surfMean[i]<100) {
		constexpr double pi = 3.14159265358979323846;
		vol = std::pow(surfMean[i]*std::pow(scale_pixel_to_meter, 2), 1.5) / (6.0 * std::sqrt(pi));
	}
	colVolMean->SetValue(i,vol);
  }

  outputTable->AddColumn(colTrajId);
  outputTable->AddColumn(colStartF);
  outputTable->AddColumn(colEndF);
  outputTable->AddColumn(colDur);
  outputTable->AddColumn(colVX);
  outputTable->AddColumn(colVY);
  outputTable->AddColumn(colSurfMin);
  outputTable->AddColumn(colSurfMax);
  outputTable->AddColumn(colSurfMean);
  outputTable->AddColumn(colVolMean);
  // ------------------------- SURFACE POINT DATA -----------------------------
  outputSurface->CopyStructure(inputDataSet);

  const vtkIdType nPts = outputSurface->GetNumberOfPoints();
  vtkIdType nCells = outputSurface->GetNumberOfCells();
  if (nPts == 0) {
    this->printErr("Input dataset has no points.");
    return 0;
  }
  auto pts = vtkSmartPointer<vtkIdList>::New();
 
  for(size_t frame = 0; frame < allVertexDebris.size(); ++frame) {

    if(static_cast<vtkIdType>(allVertexDebris[frame].size()) != nPts){
      return 0;
	}
  
    char name[32];
    std::snprintf(name, sizeof(name), "%04zu", frame);
  
    auto col = makeIntCol(name, allVertexDebris[frame].size());
    col->FillValue(-1);
  
    auto base = vtkSmartPointer<vtkIntArray>::New();
    base->SetNumberOfTuples(nPts);
    base->FillValue(-1);
  
    for(vtkIdType pid = 0; pid < nPts; ++pid) {
      const int lab = allVertexDebris[frame][static_cast<size_t>(pid)];
      base->SetValue(pid, lab);
      col->SetValue(pid, lab);
    }
  
    for(vtkIdType cid = 0; cid < nCells; ++cid) {
      outputSurface->GetCellPoints(cid, pts);
      const vtkIdType m = pts->GetNumberOfIds();
  
      int chosen = -1;
      bool conflict = false;
      bool hasLabeled = false;
      bool hasUnlabeled = false;
  
      for(vtkIdType k = 0; k < m; ++k) {
        const vtkIdType pid = pts->GetId(k);
        const int lab = base->GetValue(pid);
  
        if(lab < 0) {
          hasUnlabeled = true;
          continue;
        }
  
        hasLabeled = true;
        if(chosen < 0) {
          chosen = lab;
        } else if(lab != chosen) {
          conflict = true;
          break;
        }
      }
  
      if(!(hasLabeled && hasUnlabeled)) {
        continue;
      }
  
      const int outLab = conflict ? -2 : chosen;
      if(outLab < 0) {
        continue; 
      }
  
      for(vtkIdType k = 0; k < m; ++k) {
        const vtkIdType pid = pts->GetId(k);
        if(base->GetValue(pid) < 0) {
          col->SetValue(pid, outLab);
        }
      }
    }
    outputSurface->GetPointData()->AddArray(col);
  }

  // --------------------------- LINEAR REG && ADDED --------------------------
  
  const vtkIdType nInit  = static_cast<vtkIdType>(linearTraj.size());
  const vtkIdType nLinks = static_cast<vtkIdType>(fuseRecords.size());

  nCells = nInit + nLinks;
  
  auto ppts   = vtkSmartPointer<vtkPoints>::New();
  auto lines = vtkSmartPointer<vtkCellArray>::New();
  ppts->SetNumberOfPoints(2 * nCells);
  
  auto finalChainId = vtkSmartPointer<vtkIntArray>::New(); // chain id 
  finalChainId->SetName("FinalChainId");
  finalChainId->SetNumberOfTuples(nCells);
  
  auto inputTrajId = vtkSmartPointer<vtkIntArray>::New();  // input traj id for initial segments; -1 for links
  inputTrajId->SetName("InputTrajId");
  inputTrajId->SetNumberOfTuples(nCells);
  
  auto segmentKind = vtkSmartPointer<vtkIntArray>::New();  // 0 = initial segment, 1 = fusion link
  segmentKind->SetName("SegmentKind");
  segmentKind->SetNumberOfTuples(nCells);
  
  for(vtkIdType i = 0; i < nInit; ++i) {
	const auto &traj = linearTraj[static_cast<size_t>(i)];
	const int startF = trajTime[static_cast<size_t>(i)].front();
    const int endF   = trajTime[static_cast<size_t>(i)].back();
    addSegment(ppts, lines, i,
               traj.evalX(startF), traj.evalY(startF), startF,
               traj.evalX(endF),   traj.evalY(endF),   endF);
  
    finalChainId->SetValue(i, traj.finalChainId);
    inputTrajId->SetValue(i, static_cast<int>(i));
    segmentKind->SetValue(i, 0);
  }
  
  for(vtkIdType k = 0; k < nLinks; ++k) {
    const auto &f = fuseRecords[static_cast<size_t>(k)];
    const LinearTrajectory &ci = linearTraj[f.i];
    const LinearTrajectory &cj = linearTraj[f.j];
    
	const vtkIdType segIdx = nInit + k;
    addSegment(ppts, lines, segIdx,
               ci.evalX(f.endFrame),   ci.evalY(f.endFrame),   f.endFrame,
               cj.evalX(f.startFrame), cj.evalY(f.startFrame), f.startFrame);

    finalChainId->SetValue(segIdx, ci.finalChainId);
    inputTrajId->SetValue(segIdx, -1);
    segmentKind->SetValue(segIdx, 1);
  }
  
  
  outputLinear->SetPoints(ppts);
  outputLinear->SetCells(VTK_LINE, lines);
  outputLinear->GetCellData()->AddArray(finalChainId);
  outputLinear->GetCellData()->AddArray(inputTrajId);
  outputLinear->GetCellData()->AddArray(segmentKind);
  //--------------------------- TRAJECTORIES -----------------------
  
  const vtkIdType n = static_cast<vtkIdType>(finalTraj.size());
  
  auto mergePoints = vtkSmartPointer<vtkPoints>::New();
  auto mergeLines  = vtkSmartPointer<vtkCellArray>::New();
  mergePoints->SetNumberOfPoints(2 * n);
  
  auto mergeIdArr = makeIntCol("TrajId",   n);
  auto durArr     = makeIntCol("Duration",  n);
  auto ejecArr    = makeDblCol("AngleEjection", n);
  for(vtkIdType i = 0; i < n; ++i) {
    const ttk::TrajectoryStatistics::LinearTrajectory traj = finalTraj[i];
    double startF = extendTraj ? 0.0 : static_cast<double>(traj.startFrame);
    double endF   = static_cast<double>(traj.endFrame);
    constexpr double pi = 3.14159265358979323846;
    double ejection = atan(traj.ay/traj.ax) *180/pi;

    // Avance startF tant que le point est hors boundary
    while(startF < endF) {
      const double x = traj.evalX(startF);
      const double y = traj.evalY(startF);
      if(x >= bounds[0] && x <= bounds[1] && y >= bounds[2] && y <= bounds[3])
        break;
      startF += 1.0;
    }

    // Recule endF tant que le point est hors boundary
    while(endF > startF) {
      const double x = traj.evalX(endF);
      const double y = traj.evalY(endF);
      if(x >= bounds[0] && x <= bounds[1] && y >= bounds[2] && y <= bounds[3])
        break;
      endF -= 1.0;
    }

    addSegment(mergePoints, mergeLines, i,
               traj.evalX(startF), traj.evalY(startF), startF,
               traj.evalX(endF),   traj.evalY(endF),   endF);

    ejecArr->SetValue(i, ejection);
    mergeIdArr->SetValue(i, static_cast<int>(i));
    durArr    ->SetValue(i, static_cast<int>(traj.endFrame - startF));
  } 
  outputTraj->SetPoints(mergePoints);
  outputTraj->SetCells(VTK_LINE, mergeLines);
  outputTraj->GetCellData()->AddArray(mergeIdArr);
  outputTraj->GetCellData()->AddArray(durArr);
  outputTraj->GetCellData()->AddArray(ejecArr);
  
  this->printMsg("End TrajectoryStatistic");

  return 1;
}

int ttkTrajectoryStatistics::computeAllGradientMagnitudes(
  vtkDataSet *inputDataSet,
  const std::vector<vtkDataArray *> &inputScalarFields,
  std::vector<std::vector<double>> &gradientNorms
) {
  if(!inputDataSet) {
    this->printErr("inputDataSet is nullptr.");
    return 0;
  }

  const size_t nFields = inputScalarFields.size();
  if(nFields == 0) {
    this->printErr("pas de scalar fields fournis.");
    return 0;
  }

  vtkIdType nPts = inputDataSet->GetNumberOfPoints();
  if(nPts <= 0) {
    this->printErr("maillage vide (nPts <= 0).");
    return 0;
  }

  gradientNorms.clear();
  gradientNorms.resize(nFields);
  for(size_t f = 0; f < nFields; ++f) {
    gradientNorms[f].assign(static_cast<size_t>(nPts), 0.0);
  }

  for(size_t f = 0; f < nFields; f++) {
    vtkDataArray *currScalar = inputScalarFields[f];
    if(!currScalar || !currScalar->GetName()) {
      this->printErr("scalar array invalide en frame " + std::to_string(f));
      return 0;
    }
    const char *scalarName = currScalar->GetName();

    vtkSmartPointer<vtkGradientFilter> gradFilter = vtkSmartPointer<vtkGradientFilter>::New();
    gradFilter->SetInputData(inputDataSet);
    gradFilter->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS, scalarName);


    gradFilter->Update();

    vtkDataSet *gradOutput = gradFilter->GetOutput();
    if (!gradOutput){
        this->printErr("grad output missing");
    }

    vtkDataArray *gradArray = gradOutput->GetPointData()->GetArray("Gradients");
    if(!gradArray) {
      this->printErr(
            "recup gradArray echec"
      );
      return 0;
    }

    for(vtkIdType pid = 0; pid < nPts; ++pid) {
      double gx = gradArray->GetComponent(pid, 0);
      double gy = gradArray->GetComponent(pid, 1);
      double gz = gradArray->GetComponent(pid, 2);
      gradientNorms[f][static_cast<size_t>(pid)] = std::sqrt(gx * gx + gy * gy + gz * gz);
    }
  }

  return 1;
}

