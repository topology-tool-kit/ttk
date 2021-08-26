#include <ttkFTMTree.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

// VTK inclues
#include <vtkConnectivityFilter.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkFTMTree);

ttkFTMTree::ttkFTMTree() {
  this->setDebugMsgPrefix("FTMTree");
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
}

int ttkFTMTree::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkFTMTree::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkFTMTree::RequestData(vtkInformation *ttkNotUsed(request),
                            vtkInformationVector **inputVector,
                            vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputSkeletonNodes = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputSkeletonArcs = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto outputSegmentation = vtkDataSet::GetData(outputVector, 2);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    this->printErr("Error: input pointer is NULL.");
    return 0;
  }

  if(!input->GetNumberOfPoints()) {
    this->printErr("Error: input has no point.");
    return 0;
  }

  if(!outputSkeletonNodes || !outputSkeletonArcs || !outputSegmentation) {
    this->printErr("Error: output pointer is NULL.");
    return 0;
  }
#endif

  // Arrays

  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);
  if(!inputArray)
    return 0;

  // Connected components
  if(input->IsA("vtkUnstructuredGrid")) {
    // This data set may have several connected components,
    // we need to apply the FTM Tree for each one of these components
    // We then reconstruct the global tree using an offest mecanism
    auto inputWithId = vtkSmartPointer<vtkUnstructuredGrid>::New();
    inputWithId->ShallowCopy(input);
    identify(inputWithId);

    vtkNew<vtkConnectivityFilter> connectivity{};
    connectivity->SetInputData(inputWithId);
    connectivity->SetExtractionModeToAllRegions();
    connectivity->ColorRegionsOn();
    connectivity->Update();

    nbCC_ = connectivity->GetOutput()
              ->GetCellData()
              ->GetArray("RegionId")
              ->GetRange()[1]
            + 1;
    connected_components_.resize(nbCC_);

    if(nbCC_ > 1) {
      // Warning, in case of several connected components, the ids seen by
      // the base code will not be consistent with those of the original
      // mesh
      for(int cc = 0; cc < nbCC_; cc++) {
        vtkNew<vtkThreshold> threshold{};
        threshold->SetInputConnection(connectivity->GetOutputPort());
        threshold->SetInputArrayToProcess(
          0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionId");
        threshold->ThresholdBetween(cc, cc);
        threshold->Update();
        connected_components_[cc] = vtkSmartPointer<vtkUnstructuredGrid>::New();
        connected_components_[cc]->ShallowCopy(threshold->GetOutput());
      }
    } else {
      connected_components_[0] = inputWithId;
    }
  } else if(input->IsA("vtkPolyData")) {
    // NOTE: CC check should not be implemented on a per vtk module layer.
    nbCC_ = 1;
    connected_components_.resize(nbCC_);
    connected_components_[0] = vtkSmartPointer<vtkPolyData>::New();
    connected_components_[0]->ShallowCopy(input);
    identify(connected_components_[0]);
  } else {
    nbCC_ = 1;
    connected_components_.resize(nbCC_);
    connected_components_[0] = vtkSmartPointer<vtkImageData>::New();
    connected_components_[0]->ShallowCopy(input);
    identify(connected_components_[0]);
  }

  // now proceed for each triangulation obtained.

  if(preconditionTriangulation() == 0) {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printErr("Error : wrong triangulation.");
    return 0;
#endif
  }

  // Fill the vector of scalar/offset, cut the array in pieces if needed
  if(getScalars() == 0) {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printErr("Error : wrong input scalars.");
    return 0;
#endif
  }
  getOffsets();

  this->printMsg("Launching on field "
                 + std::string{inputScalars_[0]->GetName()});

  ttk::ftm::idNode acc_nbNodes = 0;

  // Build tree
  for(int cc = 0; cc < nbCC_; cc++) {
    ftmTree_[cc].tree.setVertexScalars(
      ttkUtils::GetVoidPointer(inputScalars_[cc]));
    ftmTree_[cc].tree.setVertexSoSoffsets(offsets_[cc].data());
    ftmTree_[cc].tree.setTreeType(GetTreeType());
    ftmTree_[cc].tree.setSegmentation(GetWithSegmentation());
    ftmTree_[cc].tree.setNormalizeIds(GetWithNormalize());

    ttkVtkTemplateMacro(inputArray->GetDataType(),
                        triangulation_[cc]->getType(),
                        (ftmTree_[cc].tree.build<VTK_TT, TTK_TT>(
                          (TTK_TT *)triangulation_[cc]->getData())));

    ftmTree_[cc].offset = acc_nbNodes;
    acc_nbNodes += ftmTree_[cc].tree.getTree(GetTreeType())->getNumberOfNodes();
  }

  UpdateProgress(0.50);

  // Construct output
  if(getSkeletonNodes(outputSkeletonNodes) == 0) {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printErr("Error : wrong properties on skeleton nodes.");
    return 0;
#endif
  }

  if(getSkeletonArcs(outputSkeletonArcs) == 0) {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printErr("Error : wrong properties on skeleton arcs.");
    return 0;
#endif
  }

  if(GetWithSegmentation()) {
    outputSegmentation->ShallowCopy(input);
    if(getSegmentation(outputSegmentation) == 0) {
#ifndef TTK_ENABLE_KAMIKAZE
      this->printErr("Error : wrong properties on segmentation.");
      return 0;
#endif
    }
  }

  UpdateProgress(1);

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
  printCSVStats();
#endif

  return 1;
}

int ttkFTMTree::addCompleteSkeletonArc(const ttk::ftm::idSuperArc arcId,
                                       const int cc,
                                       vtkPoints *points,
                                       vtkUnstructuredGrid *skeletonArcs,
                                       ttk::ftm::ArcData &arcData) {
  auto tree = ftmTree_[cc].tree.getTree(GetTreeType());
  vtkDataArray *idMapper = connected_components_[cc]->GetPointData()->GetArray(
    ttk::VertexScalarFieldName);
  auto arc = tree->getSuperArc(arcId);
  float point[3];
  vtkIdType pointIds[2];

  using ttk::SimplexId;
  const SimplexId downNodeId = tree->getLowerNodeId(arc);
  const SimplexId l_downVertexId = tree->getNode(downNodeId)->getVertexId();
  const SimplexId g_downVertexId = idMapper->GetTuple1(l_downVertexId);
  triangulation_[cc]->getVertexPoint(
    l_downVertexId, point[0], point[1], point[2]);
  const double scalarMin = inputScalars_[cc]->GetTuple1(l_downVertexId);

  // Get or create first point of the arc
  SimplexId nextPointId;
  if(!arcData.hasPoint(g_downVertexId)) {
    nextPointId = points->InsertNextPoint(point);
    arcData.addPoint(g_downVertexId, nextPointId, scalarMin, false);
  } else {
    nextPointId = arcData.point_ids[g_downVertexId];
  }

  pointIds[0] = nextPointId;

  for(const SimplexId vertexId : *arc) {
    triangulation_[cc]->getVertexPoint(vertexId, point[0], point[1], point[2]);
    pointIds[1] = points->InsertNextPoint(point);
    const double scalar = inputScalars_[cc]->GetTuple1(vertexId);
    arcData.setPoint(pointIds[1], scalar, true);

    const SimplexId nextCell
      = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
    arcData.fillArrayCell(
      nextCell, arcId, ftmTree_[cc], triangulation_[cc], params_);

    pointIds[0] = pointIds[1];
  }

  const SimplexId upNodeId = tree->getUpperNodeId(arc);
  const SimplexId l_upVertexId = tree->getNode(upNodeId)->getVertexId();
  const SimplexId g_upVertexId = idMapper->GetTuple1(l_upVertexId);
  triangulation_[cc]->getVertexPoint(
    l_upVertexId, point[0], point[1], point[2]);
  const double scalarMax = inputScalars_[cc]->GetTuple1(l_upVertexId);

  // Get or create last point of the arc
  if(!arcData.hasPoint(g_upVertexId)) {
    nextPointId = points->InsertNextPoint(point);
    arcData.addPoint(g_upVertexId, nextPointId, scalarMax, false);
  } else {
    nextPointId = arcData.point_ids[g_upVertexId];
  }

  pointIds[1] = nextPointId;

  const SimplexId nextCell
    = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
  arcData.fillArrayCell(
    nextCell, arcId, ftmTree_[cc], triangulation_[cc], params_);

  return 1;
}

int ttkFTMTree::addDirectSkeletonArc(const ttk::ftm::idSuperArc arcId,
                                     const int cc,
                                     vtkPoints *points,
                                     vtkUnstructuredGrid *skeletonArcs,
                                     ttk::ftm::ArcData &arcData) {
  auto tree = ftmTree_[cc].tree.getTree(GetTreeType());
  vtkDataArray *idMapper = connected_components_[cc]->GetPointData()->GetArray(
    ttk::VertexScalarFieldName);
  auto arc = tree->getSuperArc(arcId);
  float point[3];
  vtkIdType pointIds[2];

  using ttk::SimplexId;
  const SimplexId downNodeId = tree->getLowerNodeId(arc);
  const SimplexId l_downVertexId = tree->getNode(downNodeId)->getVertexId();
  const SimplexId g_downVertexId = idMapper->GetTuple1(l_downVertexId);
  triangulation_[cc]->getVertexPoint(
    l_downVertexId, point[0], point[1], point[2]);
  const double scalarMin = inputScalars_[cc]->GetTuple1(l_downVertexId);
  // Get or create first point of the arc
  if(!arcData.hasPoint(g_downVertexId)) {
    const SimplexId nextPointId = points->InsertNextPoint(point);
    pointIds[0] = nextPointId;
    arcData.addPoint(g_downVertexId, nextPointId, scalarMin, false);
  } else {
    pointIds[0] = arcData.point_ids[g_downVertexId];
  }

  const SimplexId upNodeId = tree->getUpperNodeId(arc);
  const SimplexId l_upVertexId = tree->getNode(upNodeId)->getVertexId();
  const SimplexId g_upVertexId = idMapper->GetTuple1(l_upVertexId);
  triangulation_[cc]->getVertexPoint(
    l_upVertexId, point[0], point[1], point[2]);
  const double scalarMax = inputScalars_[cc]->GetTuple1(l_upVertexId);
  // Get or create last point of the arc
  if(!arcData.hasPoint(g_upVertexId)) {
    const SimplexId nextPointId = points->InsertNextPoint(point);
    pointIds[1] = nextPointId;
    arcData.addPoint(g_upVertexId, nextPointId, scalarMax, false);
  } else {
    pointIds[1] = arcData.point_ids[g_upVertexId];
  }

  const SimplexId nextCell
    = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
  arcData.fillArrayCell(
    nextCell, arcId, ftmTree_[cc], triangulation_[cc], params_);

  return 1;
}

int ttkFTMTree::addSampledSkeletonArc(const ttk::ftm::idSuperArc arcId,
                                      const int cc,
                                      vtkPoints *points,
                                      vtkUnstructuredGrid *skeletonArcs,
                                      ttk::ftm::ArcData &arcData) {
  auto tree = ftmTree_[cc].tree.getTree(GetTreeType());
  vtkDataArray *idMapper = connected_components_[cc]->GetPointData()->GetArray(
    ttk::VertexScalarFieldName);
  auto arc = tree->getSuperArc(arcId);
  float point[3];
  vtkIdType pointIds[2];

  using ttk::SimplexId;
  const SimplexId downNodeId = tree->getLowerNodeId(arc);
  const SimplexId l_downVertexId = tree->getNode(downNodeId)->getVertexId();
  const SimplexId g_downVertexId = idMapper->GetTuple1(l_downVertexId);
  triangulation_[cc]->getVertexPoint(
    l_downVertexId, point[0], point[1], point[2]);
  const double scalarMin = inputScalars_[cc]->GetTuple1(l_downVertexId);

  // Get or create first point of the arc
  SimplexId nextPointId;
  if(!arcData.hasPoint(g_downVertexId)) {
    nextPointId = points->InsertNextPoint(point);
    arcData.addPoint(g_downVertexId, nextPointId, scalarMin, false);
  } else {
    nextPointId = arcData.point_ids[g_downVertexId];
  }

  pointIds[0] = nextPointId;

  const SimplexId upNodeId = tree->getUpperNodeId(arc);
  const SimplexId l_upVertexId = tree->getNode(upNodeId)->getVertexId();
  const SimplexId g_upVertexId = idMapper->GetTuple1(l_upVertexId);
  triangulation_[cc]->getVertexPoint(
    l_upVertexId, point[0], point[1], point[2]);
  const double scalarMax = inputScalars_[cc]->GetTuple1(l_upVertexId);

  const double delta = (scalarMax - scalarMin) / (params_.samplingLvl + 1);
  double scalarLimit = scalarMin + delta;
  double scalarAvg = 0;

  // Get or create last point of the arc
  if(!arcData.hasPoint(g_upVertexId)) {
    nextPointId = points->InsertNextPoint(point);
    arcData.addPoint(g_upVertexId, nextPointId, scalarMax, false);
  } else {
    nextPointId = arcData.point_ids[g_upVertexId];
  }

  SimplexId c = 0;
  float sum[3]{0, 0, 0};
  for(const SimplexId vertexId : *arc) {
    triangulation_[cc]->getVertexPoint(vertexId, point[0], point[1], point[2]);
    const double scalarVertex = inputScalars_[cc]->GetTuple1(vertexId);

    if(scalarVertex < scalarLimit) {
      sum[0] += point[0];
      sum[1] += point[1];
      sum[2] += point[2];
      scalarAvg += scalarVertex;
      ++c;
    } else {
      if(c) {
        sum[0] /= c;
        sum[1] /= c;
        sum[2] /= c;
        scalarAvg /= c;

        pointIds[1] = points->InsertNextPoint(sum);
        arcData.setPoint(pointIds[1], scalarAvg, true);
        const SimplexId nextCell
          = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
        arcData.fillArrayCell(
          nextCell, arcId, ftmTree_[cc], triangulation_[cc], params_);

        pointIds[0] = pointIds[1];
      }

      scalarLimit += delta;
      sum[0] = 0;
      sum[1] = 0;
      sum[2] = 0;
      scalarAvg = 0;
      c = 0;
    }
  }

  // The up id
  pointIds[1] = nextPointId;

  const SimplexId nextCell
    = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
  arcData.fillArrayCell(
    nextCell, arcId, ftmTree_[cc], triangulation_[cc], params_);

  return 1;
}

int ttkFTMTree::getOffsets() {
  // should be called after getScalars for inputScalars_ needs to be filled

  offsets_.resize(nbCC_);
  for(int cc = 0; cc < nbCC_; cc++) {
    const auto offsets = this->GetOrderArray(
      connected_components_[cc], 0, 1, ForceInputOffsetScalarField);

    offsets_[cc].resize(connected_components_[cc]->GetNumberOfPoints());

    for(size_t i = 0; i < offsets_[cc].size(); i++) {
      offsets_[cc][i] = offsets->GetTuple1(i);
    }

#ifndef TTK_ENABLE_KAMIKAZE
    if(offsets_[cc].empty()) {
      this->printMsg(
        {"Error : wrong input offset scalar field for ", std::to_string(cc)},
        ttk::debug::Priority::ERROR);
      return 0;
    }
#endif
  }

  return 1;
}

int ttkFTMTree::getScalars() {
  inputScalars_.resize(nbCC_);
  for(int cc = 0; cc < nbCC_; cc++) {
    inputScalars_[cc]
      = this->GetInputArrayToProcess(0, connected_components_[cc]);

#ifndef TTK_ENABLE_KAMIKAZE
    if(!inputScalars_[cc]) {
      this->printMsg({"Error : input scalar ", std::to_string(cc),
                      " field pointer is null."},
                     ttk::debug::Priority::ERROR);
      return 0;
    }
#endif
  }

  return 1;
}

int ttkFTMTree::getSegmentation(vtkDataSet *outputSegmentation) {
  ttk::ftm::VertData vertData;
  vertData.init(ftmTree_, params_);

  for(int cc = 0; cc < nbCC_; cc++) {
    auto tree = ftmTree_[cc].tree.getTree(GetTreeType());
    vtkDataArray *idMapper
      = connected_components_[cc]->GetPointData()->GetArray(
        ttk::VertexScalarFieldName);
    const auto numberOfSuperArcs = tree->getNumberOfSuperArcs();
    // #pragma omp for
    for(ttk::ftm::idSuperArc arcId = 0; arcId < numberOfSuperArcs; ++arcId) {
      vertData.fillArrayPoint(
        arcId, ftmTree_[cc], triangulation_[cc], idMapper, params_);
    }
  }

  vtkPointData *pointData = outputSegmentation->GetPointData();
  vertData.addArray(pointData, params_);
  // vertex identifier field should not be propagated (it has caused
  // downstream bugs in distanceField in oceanVortices)
  pointData->RemoveArray(ttk::VertexScalarFieldName);

  return 1;
}

int ttkFTMTree::getSkeletonArcs(vtkUnstructuredGrid *outputSkeletonArcs) {
  vtkNew<vtkUnstructuredGrid> skeletonArcs{};
  vtkNew<vtkPoints> points{};

  ttk::ftm::ArcData arcData;
  arcData.init(ftmTree_, params_);

  const int samplingLevel = params_.samplingLvl;
  for(int cc = 0; cc < nbCC_; cc++) {
    auto tree = ftmTree_[cc].tree.getTree(GetTreeType());

    using ttk::SimplexId;
    const SimplexId numberOfSuperArcs = tree->getNumberOfSuperArcs();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!numberOfSuperArcs) {
      this->printErr("Error : tree has no super arcs.");
      return 0;
    }
#endif

    for(SimplexId arcId = 0; arcId < numberOfSuperArcs; ++arcId) {
      const SimplexId numberOfRegularNodes = tree->getArcSize(arcId);
      if(numberOfRegularNodes > 0 and samplingLevel > 0) {
        addSampledSkeletonArc(arcId, cc, points, skeletonArcs, arcData);
      } else if(samplingLevel == -1) {
        addCompleteSkeletonArc(arcId, cc, points, skeletonArcs, arcData);
      } else {
        addDirectSkeletonArc(arcId, cc, points, skeletonArcs, arcData);
      }
    }
  }

  skeletonArcs->SetPoints(points);
  arcData.addArray(skeletonArcs, params_);
  outputSkeletonArcs->ShallowCopy(skeletonArcs);

  // const SimplexId p_size = points->GetNumberOfPoints();
  // const SimplexId s_size = tree->getNumberOfVertices();
  // cout << "arcs points " << p_size << endl;
  // cout << "scal points " << s_size << endl;
  // cout << "nb arcs     " << tree->getNumberOfSuperArcs()<< endl;
  // if(p_size != s_size){
  //    exit(3);
  // }

  return 1;
}

int ttkFTMTree::getSkeletonNodes(vtkUnstructuredGrid *outputSkeletonNodes) {
  vtkNew<vtkUnstructuredGrid> skeletonNodes{};
  vtkNew<vtkPoints> points{};

  ttk::ftm::NodeData nodeData;
  nodeData.init(ftmTree_, params_);
  nodeData.setScalarType(inputScalars_[0]->GetDataType());

  for(int cc = 0; cc < nbCC_; cc++) {
    auto tree = ftmTree_[cc].tree.getTree(GetTreeType());
    vtkDataArray *idMapper
      = connected_components_[cc]->GetPointData()->GetArray(
        ttk::VertexScalarFieldName);

    const auto numberOfNodes = tree->getNumberOfNodes();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!numberOfNodes) {
      this->printErr("Error : tree has no nodes.");
      return 0;
    }
#endif

    for(ttk::ftm::idNode nodeId = 0; nodeId < numberOfNodes; ++nodeId) {
      const auto node = tree->getNode(nodeId);
#ifndef TTK_ENABLE_KAMIKAZE
      if(!node) {
        this->printMsg({"Error : node ", std::to_string(nodeId), " is null."},
                       ttk::debug::Priority::ERROR);
        return 0;
      }
#endif
      const auto local_vertId = node->getVertexId();
      float point[3];
      triangulation_[cc]->getVertexPoint(
        local_vertId, point[0], point[1], point[2]);
      const auto nextPoint = points->InsertNextPoint(point);
      nodeData.fillArrayPoint(
        nextPoint, nodeId, ftmTree_[cc], idMapper, triangulation_[cc], params_);
    }
  }

  skeletonNodes->SetPoints(points);
  vtkPointData *pointData = skeletonNodes->GetPointData();
  nodeData.addArray(pointData, params_);
  outputSkeletonNodes->ShallowCopy(skeletonNodes);

  return 1;
}

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
void ttkFTMTree::printCSVStats() {
  for(auto &t : ftmTree_) {
    switch(GetTreeType()) {
      case ftm::TreeType::Join:
        cout << "JT" << endl;
        printCSVTree(t.tree.getTree(ftm::TreeType::Join));
        break;
      case ftm::TreeType::Split:
        cout << "ST" << endl;
        printCSVTree(t.tree.getTree(ftm::TreeType::Split));
        break;
      default:
        cout << "JT" << endl;
        printCSVTree(t.tree.getTree(ftm::TreeType::Join));
        cout << "ST" << endl;
        printCSVTree(t.tree.getTree(ftm::TreeType::Split));
        break;
    }
  }
}

void ttkFTMTree::printCSVTree(const ftm::FTMTree_MT *const tree) const {
  const idSuperArc nbArc = tree->getNumberOfLeaves();
  cout << "begin; ";
  for(idSuperArc a = 0; a < nbArc; a++) {
    cout << tree->getActiveTasks(a).begin << "; ";
  }
  cout << endl << "end; ";
  for(idSuperArc a = 0; a < nbArc; a++) {
    cout << tree->getActiveTasks(a).end << "; ";
  }
  cout << endl << "origing; ";
  for(idSuperArc a = 0; a < nbArc; a++) {
    cout << tree->getActiveTasks(a).origin << "; ";
  }
  cout << endl;
}
#endif

int ttkFTMTree::preconditionTriangulation() {
  triangulation_.resize(nbCC_);
  ftmTree_.resize(nbCC_);

  for(int cc = 0; cc < nbCC_; cc++) {
    triangulation_[cc]
      = ttkAlgorithm::GetTriangulation(connected_components_[cc]);
#ifndef TTK_ENABLE_KAMIKAZE
    if(!triangulation_[cc]) {
      this->printErr("Error : ttkTriangulation::getTriangulation() is null.");
      return 0;
    }
#endif

    ftmTree_[cc].tree.setDebugLevel(debugLevel_);
    ftmTree_[cc].tree.setThreadNumber(threadNumber_);
    ftmTree_[cc].tree.preconditionTriangulation(triangulation_[cc]);

#ifndef TTK_ENABLE_KAMIKAZE
    if(triangulation_[cc]->isEmpty()) {
      this->printMsg({"Error : ttkTriangulation on connected component",
                      std::to_string(cc), " allocation problem."},
                     ttk::debug::Priority::ERROR);
      return 0;
    }
#endif
  }
  return 1;
}

// protected

void ttkFTMTree::identify(vtkDataSet *ds) const {
  vtkNew<ttkSimplexIdTypeArray> identifiers{};
  const auto nbPoints = ds->GetNumberOfPoints();
  identifiers->SetName(ttk::VertexScalarFieldName);
  identifiers->SetNumberOfComponents(1);
  identifiers->SetNumberOfTuples(nbPoints);

  for(int i = 0; i < nbPoints; i++) {
    identifiers->SetTuple1(i, i);
  }

  ds->GetPointData()->AddArray(identifiers);
}
