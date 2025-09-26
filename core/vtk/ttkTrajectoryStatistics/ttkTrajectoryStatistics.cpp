
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
//  this->setCoordCratere(coordCratere);
//  this->setCraterAngle(threshCratereAngle);
  this->setMaxX(maxX);
  this->setMaxY(maxY);
  this->setMinY(minY);
  this->setMinX(minX);
  this->setSurfaceMethod(surfaceMethod);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation) return 0;
  this->preconditionTriangulation(triangulation);

  std::vector<std::vector<double>> newTraj(numTraj);
  std::vector<std::vector<double>> finalTraj;
  std::vector<FuseRecord> fuseRecords;

  this->correctTrajectory(trajTime, trajX, trajY, finalTraj,newTraj, fuseRecords); 

  std::vector<int>  durations(numTraj);
  std::vector<double> VX(numTraj), VY(numTraj), 
                      surfMin(numTraj), surfMax(numTraj), surfMean(numTraj);
  std::vector<std::vector<ttk::SimplexId>> allVertexDebris(numTraj);
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
                        errSurf,
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

  for(int i = 0; i < numMerge; ++i) {
    colStartF->SetValue(i, static_cast<int>(finalTraj[i][4]));
    colEndF  ->SetValue(i, static_cast<int>(finalTraj[i][5]));
    colDur   ->SetValue(i, durations[i]); 
    colVX    ->SetValue(i, VX[i]);
    colVY    ->SetValue(i, VY[i]);
  }
  
  const size_t F = finalTraj.size();
  std::vector<double> aggMin(F, std::numeric_limits<double>::infinity());
  std::vector<double> aggMax(F, -std::numeric_limits<double>::infinity());
  std::vector<double> aggSum(F, 0.0);
  std::vector<size_t> aggCnt(F, 0);
  
  for(size_t ti = 0; ti < newTraj.size(); ++ti) {
    if(newTraj[ti].size() > 4) {
      const int fid = static_cast<int>(newTraj[ti][4]);
      if(fid >= 0 && fid < static_cast<int>(F)) {
        aggMin[fid] = std::min(aggMin[fid], surfMin[ti]);
        aggMax[fid] = std::max(aggMax[fid], surfMax[ti]);
        aggSum[fid] += surfMean[ti];
        aggCnt[fid] += 1;
      }
    }
  }
  
  double sMean = 0.0;
  for(int i = 0; i < numMerge; ++i) {
    double mean = (aggCnt[i] > 0 ? (aggSum[i] / static_cast<double>(aggCnt[i])) : 0.0);
    const double vMin = (aggCnt[i] > 0 ? aggMin[i] : 0.0);
    const double vMax = (aggCnt[i] > 0 ? aggMax[i] : 0.0);
    colSurfMin ->SetValue(i, vMin);
    colSurfMax ->SetValue(i, vMax);
    colSurfMean->SetValue(i, mean);
    sMean += mean;
  }

  outputTable->AddColumn(colStartF);
  outputTable->AddColumn(colEndF);
  outputTable->AddColumn(colDur);
  outputTable->AddColumn(colVX);
  outputTable->AddColumn(colVY);
  outputTable->AddColumn(colSurfMin);
  outputTable->AddColumn(colSurfMax);
  outputTable->AddColumn(colSurfMean);
  
  // ------------------------- SURFACE POINT DATA -----------------------------

  vtkSmartPointer<vtkDataSet> surfOutput =
    vtkSmartPointer<vtkDataSet>::Take(inputDataSet->NewInstance());
  surfOutput->ShallowCopy(inputDataSet);
  outputSurface->ShallowCopy(surfOutput);
  
  //  -1 = no surface , >=0 = finalId, -2 = double traj 
  auto surfaceFinalId = vtkSmartPointer<vtkIntArray>::New();
  surfaceFinalId->SetName("FinalTrajId"); 
  const vtkIdType nPts = inputDataSet->GetNumberOfPoints();
  surfaceFinalId->SetNumberOfTuples(nPts);
  surfaceFinalId->FillComponent(0, -1);
  
  for(size_t i = 0; i < allVertexDebris.size(); ++i) {
    const auto &trajSurface = allVertexDebris[i];
  
    int finalId = -1;
    if(i < newTraj.size() && newTraj[i].size() > 4 && static_cast<int>(newTraj[i][4]) != -1) {
      finalId = static_cast<int>(newTraj[i][4]);
    }
  
    for(const auto v : trajSurface) {
      if(v < 0 || v >= nPts) continue;
  
      const int current = surfaceFinalId->GetValue(v);
      if(current == -1) {
        surfaceFinalId->SetValue(v, finalId);
      } else {
        surfaceFinalId->SetValue(v, -2);
      }
    }
  }

  outputSurface->GetPointData()->AddArray(surfaceFinalId);
 
  // --------------------------- LINEAR REG && ADDED --------------------------
  
  const vtkIdType nInit  = static_cast<vtkIdType>(newTraj.size());
  const vtkIdType nLinks = static_cast<vtkIdType>(fuseRecords.size());
  const vtkIdType nCells = nInit + nLinks;
  
  auto pts   = vtkSmartPointer<vtkPoints>::New();
  auto lines = vtkSmartPointer<vtkCellArray>::New();
  pts->SetNumberOfPoints(2 * nCells);
  
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
    const auto &coef = newTraj[static_cast<size_t>(i)];
    const int startF = trajTime[static_cast<size_t>(i)].front();
    const int endF   = trajTime[static_cast<size_t>(i)].back();
  
    const double x0 = evalX(coef, startF);
    const double y0 = evalY(coef, startF);
    const double x1 = evalX(coef, endF);
    const double y1 = evalY(coef, endF);
  
    addSegment(pts, lines, i, x0, y0, startF, x1, y1, endF);
  
    finalChainId->SetValue(i, static_cast<int>(coef[4]));
    inputTrajId->SetValue(i, static_cast<int>(i));
    segmentKind->SetValue(i, 0); // initial
  }
  
  for(vtkIdType k = 0; k < nLinks; ++k) {
    const auto &f = fuseRecords[static_cast<size_t>(k)];
    const auto &ci = newTraj[static_cast<size_t>(f.i)];
    const auto &cj = newTraj[static_cast<size_t>(f.j)];
  
    const int startF = f.startFrame;
    const int endF   = f.endFrame;
  
    const double xEndI = evalX(ci, endF);
    const double yEndI = evalY(ci, endF);
    const double xStartJ = evalX(cj, startF);
    const double yStartJ = evalY(cj, startF);
  
    const vtkIdType segIdx = nInit + k;
    addSegment(pts, lines, segIdx, xEndI, yEndI, endF, xStartJ, yStartJ, startF);
  
    finalChainId->SetValue(segIdx, static_cast<int>(ci[4]));
    inputTrajId->SetValue(segIdx, -1);  
    segmentKind->SetValue(segIdx, 1);   
  }
  
  
  outputLinear->SetPoints(pts);
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
  
  for(vtkIdType i = 0; i < n; ++i) {
    const auto &coef = finalTraj[static_cast<size_t>(i)];
    const double startF = extendTraj ? 0.0 : static_cast<double>(coef[4]);
    const double endF   = static_cast<double>(coef[5]);
  
    const double x0 = evalX(coef, startF);
    const double y0 = evalY(coef, startF);
    const double x1 = evalX(coef, endF);
    const double y1 = evalY(coef, endF);
  
    addSegment(mergePoints, mergeLines, i, x0, y0, startF, x1, y1, endF);
  
    mergeIdArr->SetValue(i, static_cast<int>(i));
    durArr    ->SetValue(i, static_cast<int>(endF - startF));
  }
  
  outputTraj->SetPoints(mergePoints);
  outputTraj->SetCells(VTK_LINE, mergeLines);
  outputTraj->GetCellData()->AddArray(mergeIdArr);
  outputTraj->GetCellData()->AddArray(durArr);
  
  this->printMsg("End TrajectoryStatistic");

  //  copy of the input grid with per-segment final trajectory id ---
  /*
  vtkSmartPointer<vtkUnstructuredGrid> segCopy = vtkSmartPointer<vtkUnstructuredGrid>::Take(inputGrid->NewInstance());
  segCopy->ShallowCopy(inputGrid);

  vtkIdType inNumCells = inputGrid->GetNumberOfCells();
  vtkSmartPointer<vtkIntArray> cellLinearFinalId = vtkSmartPointer<vtkIntArray>::New();
  cellLinearFinalId->SetName("linearFinalId");
  cellLinearFinalId->SetNumberOfTuples(inNumCells);
  cellLinearFinalId->FillComponent(0, -1);

  for(size_t i = 0; i < cellsPerTraj.size(); ++i) {
    int finalId = -1;
    if(i < newTraj.size() && newTraj[i].size() > 4) {
      finalId = static_cast<int>(newTraj[i][4]);
    }
    const auto &cells = cellsPerTraj[i];
    for(const auto cId : cells) {
      if(cId >= 0 && cId < inNumCells) {
        cellLinearFinalId->SetValue(cId, finalId);
      }
    }
  }

  outputSegmentsLabeled->ShallowCopy(segCopy);
  outputSegmentsLabeled->GetCellData()->AddArray(cellLinearFinalId);
  */
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

