#include <ttkFTMTree.h>

// only used on the cpp
#include <vtkConnectivityFilter.h>
#include <vtkDataObject.h>
#include <vtkThreshold.h>

using namespace std;
using namespace ttk;

using namespace ftm;

vtkStandardNewMacro(ttkFTMTree);

int ttkFTMTree::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}

int ttkFTMTree::FillOutputPortInformation(int port, vtkInformation *info) {
  switch(port) {
    case 0:
    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;

    case 2:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
      break;
  }

  return 1;
}

int ttkFTMTree::addCompleteSkeletonArc(const ftm::idSuperArc arcId,
                                       const int cc,
                                       vtkPoints *points,
                                       vtkUnstructuredGrid *skeletonArcs,
                                       ttk::ftm::ArcData &arcData) {
  FTMTree_MT *tree = ftmTree_[cc].tree.getTree(GetTreeType());
  vtkDataArray *idMapper = connected_components_[cc]->GetPointData()->GetArray(
    ttk::VertexScalarFieldName);
  SuperArc *arc = tree->getSuperArc(arcId);
  float point[3];
  vtkIdType pointIds[2];

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

  return 0;
}

int ttkFTMTree::addDirectSkeletonArc(const idSuperArc arcId,
                                     const int cc,
                                     vtkPoints *points,
                                     vtkUnstructuredGrid *skeletonArcs,
                                     ttk::ftm::ArcData &arcData) {
  FTMTree_MT *tree = ftmTree_[cc].tree.getTree(GetTreeType());
  vtkDataArray *idMapper = connected_components_[cc]->GetPointData()->GetArray(
    ttk::VertexScalarFieldName);
  SuperArc *arc = tree->getSuperArc(arcId);
  float point[3];
  vtkIdType pointIds[2];

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

  return 0;
}

int ttkFTMTree::addSampledSkeletonArc(const idSuperArc arcId,
                                      const int cc,
                                      vtkPoints *points,
                                      vtkUnstructuredGrid *skeletonArcs,
                                      ttk::ftm::ArcData &arcData) {
  FTMTree_MT *tree = ftmTree_[cc].tree.getTree(GetTreeType());
  vtkDataArray *idMapper = connected_components_[cc]->GetPointData()->GetArray(
    ttk::VertexScalarFieldName);
  SuperArc *arc = tree->getSuperArc(arcId);
  float point[3];
  vtkIdType pointIds[2];

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

  return 0;
}

int ttkFTMTree::doIt(vector<vtkDataSet *> &inputs,
                     vector<vtkDataSet *> &outputs) {
  Memory m;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputs.size()) {
    cerr << "[ttkFTMTree] Error: not enough input information." << endl;
    return -1;
  }
#endif

  vtkUnstructuredGrid *outputSkeletonNodes
    = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  vtkUnstructuredGrid *outputSkeletonArcs
    = vtkUnstructuredGrid::SafeDownCast(outputs[1]);
  vtkDataSet *outputSegmentation = outputs[2];

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputs[0]) {
    cerr << "[ttkFTMTree] Error: input pointer is NULL." << endl;
    return -1;
  }

  if(!inputs[0]->GetNumberOfPoints()) {
    cerr << "[ttkFTMTree] Error: input has no point." << endl;
    return -1;
  }

  if(!outputSkeletonNodes || !outputSkeletonArcs || !outputSegmentation) {
    cerr << "[ttkFTMTree] Error: output pointer is NULL." << endl;
    return -1;
  }
#endif

  if(inputs[0]->IsA("vtkUnstructuredGrid")) {
    // This data set may have several connected components,
    // we need to apply the FTM Tree for each one of these components
    // We then reconstruct the global tree using an offest mecanism
    vtkSmartPointer<ttkUnstructuredGrid> input
      = vtkSmartPointer<ttkUnstructuredGrid>::New();
    input->ShallowCopy(inputs[0]);
    identify(input);

    vtkSmartPointer<vtkConnectivityFilter> connectivity
      = vtkSmartPointer<vtkConnectivityFilter>::New();
    connectivity->SetInputData(input);
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
        vtkSmartPointer<vtkThreshold> threshold
          = vtkSmartPointer<vtkThreshold>::New();
        threshold->SetInputConnection(connectivity->GetOutputPort());
        threshold->SetInputArrayToProcess(
          0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionId");
        threshold->ThresholdBetween(cc, cc);
        threshold->Update();
        connected_components_[cc] = vtkSmartPointer<ttkUnstructuredGrid>::New();
        connected_components_[cc]->ShallowCopy(threshold->GetOutput());
      }
    } else {
      connected_components_[0] = vtkSmartPointer<ttkUnstructuredGrid>::New();
      connected_components_[0]->ShallowCopy(input);
    }
  } else if(inputs[0]->IsA("vtkPolyData")) {
    // NOTE: CC check should not be implemented on a per vtk module layer.
    nbCC_ = 1;
    connected_components_.resize(nbCC_);
    connected_components_[0] = vtkSmartPointer<ttkPolyData>::New();
    connected_components_[0]->ShallowCopy(inputs[0]);
    identify(connected_components_[0]);
  } else {
    nbCC_ = 1;
    connected_components_.resize(nbCC_);
    connected_components_[0] = vtkSmartPointer<ttkImageData>::New();
    connected_components_[0]->ShallowCopy(inputs[0]);
    identify(connected_components_[0]);
  }

  // now proceed for each triangulation obtained.

  if(setupTriangulation()) {
#ifndef TTK_ENABLE_KAMIKAZE
    cerr << "[ttkFTMTree] Error : wrong triangulation." << endl;
    return -1;
#endif
  }

  // Fill the vector of scalar/offset, cut the array in pieces if needed
  if(getScalars()) {
#ifndef TTK_ENABLE_KAMIKAZE
    cerr << "[ttkFTMTree] Error : wrong input scalars." << endl;
    return -1;
#endif
  }
  getOffsets();

  {
    stringstream msg;
    msg << "[ttkFTMTree] Launch on field: " << ScalarField << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  ftm::idNode acc_nbNodes = 0;

  // Build tree
  for(int cc = 0; cc < nbCC_; cc++) {
    ftmTree_[cc].tree.setVertexScalars(inputScalars_[cc]->GetVoidPointer(0));
    ftmTree_[cc].tree.setVertexSoSoffsets(offsets_[cc].data());
    ftmTree_[cc].tree.setTreeType(GetTreeType());
    ftmTree_[cc].tree.setSegmentation(GetWithSegmentation());
    ftmTree_[cc].tree.setNormalizeIds(GetWithNormalize());

    switch(inputScalars_[cc]->GetDataType()) {
      vtkTemplateMacro((ftmTree_[cc].tree.build<VTK_TT, SimplexId>()));
    }

    ftmTree_[cc].offset = acc_nbNodes;
    acc_nbNodes += ftmTree_[cc].tree.getTree(GetTreeType())->getNumberOfNodes();
  }

  UpdateProgress(0.50);

  // Construct output
  if(getSkeletonNodes(outputSkeletonNodes)) {
#ifndef TTK_ENABLE_KAMIKAZE
    cerr << "[ttkFTMTree] Error : wrong properties on skeleton nodes." << endl;
    return -7;
#endif
  }

  if(getSkeletonArcs(outputSkeletonArcs)) {
#ifndef TTK_ENABLE_KAMIKAZE
    cerr << "[ttkFTMTree] Error : wrong properties on skeleton arcs." << endl;
    return -8;
#endif
  }

  if(GetWithSegmentation()) {
    outputSegmentation->ShallowCopy(inputs[0]);
    if(getSegmentation(outputSegmentation)) {
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkFTMTree] Error : wrong properties on segmentation." << endl;
      return -9;
#endif
    }
  }

  UpdateProgress(1);

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
  printCSVStats();
#endif

  {
    stringstream msg;
    msg << "[ttkFTMTree] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkFTMTree::getOffsets() {
  offsets_.resize(nbCC_);
  for(int cc = 0; cc < nbCC_; cc++) {
    vtkDataArray *inputOffsets;
    if(OffsetFieldId != -1) {
      inputOffsets
        = connected_components_[cc]->GetPointData()->GetArray(OffsetFieldId);
      if(inputOffsets) {
        InputOffsetScalarFieldName = inputOffsets->GetName();
        ForceInputOffsetScalarField = true;
      }
    }

    const SimplexId numberOfVertices
      = connected_components_[cc]->GetNumberOfPoints();

    if(ForceInputOffsetScalarField and InputOffsetScalarFieldName.length()) {
      inputOffsets = connected_components_[cc]->GetPointData()->GetArray(
        InputOffsetScalarFieldName.data());
      offsets_[cc].resize(numberOfVertices);
      for(SimplexId i = 0; i < numberOfVertices; i++) {
        offsets_[cc][i] = inputOffsets->GetTuple1(i);
      }
    } else if(connected_components_[cc]->GetPointData()->GetArray(
                ttk::OffsetScalarFieldName)) {
      inputOffsets = connected_components_[cc]->GetPointData()->GetArray(
        ttk::OffsetScalarFieldName);
      offsets_[cc].resize(numberOfVertices);
      for(SimplexId i = 0; i < numberOfVertices; i++) {
        offsets_[cc][i] = inputOffsets->GetTuple1(i);
      }
    } else {
      if(hasUpdatedMesh_ and offsets_[cc].size()) {
        // don't keep an out-dated offset array
        offsets_[cc].clear();
      }

      if(offsets_[cc].empty()) {
        offsets_[cc].resize(numberOfVertices);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
        for(SimplexId i = 0; i < numberOfVertices; i++) {
          offsets_[cc][i] = i;
        }
      }
    }

#ifndef TTK_ENABLE_KAMIKAZE
    if(offsets_[cc].empty()) {
      cerr << "[ttkFTMTree] Error : wrong input offset scalar field for " << cc
           << endl;
      return -1;
    }
#endif
  }

  return 0;
}

int ttkFTMTree::getScalars() {
  inputScalars_.resize(nbCC_);
  for(int cc = 0; cc < nbCC_; cc++) {
    vtkPointData *pointData = connected_components_[cc]->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData) {
      cerr << "[ttkFTMTree] Error : input cc " << cc << " has no point data."
           << endl;
      return -1;
    }
#endif

    if(ScalarField.length()) {
      inputScalars_[cc] = pointData->GetArray(ScalarField.data());
    } else {
      inputScalars_[cc] = pointData->GetArray(ScalarFieldId);
      if(inputScalars_[cc])
        ScalarField = inputScalars_[cc]->GetName();
    }

#ifndef TTK_ENABLE_KAMIKAZE
    if(!inputScalars_[cc]) {
      cerr << "[ttkFTMTree] Error : input scalar " << cc
           << " field pointer is null." << endl;
      return -3;
    }
#endif
  }

  return 0;
}

int ttkFTMTree::getSegmentation(vtkDataSet *outputSegmentation) {
  VertData vertData;
  vertData.init(ftmTree_, params_);

  for(int cc = 0; cc < nbCC_; cc++) {
    FTMTree_MT *tree = ftmTree_[cc].tree.getTree(GetTreeType());
    vtkDataArray *idMapper
      = connected_components_[cc]->GetPointData()->GetArray(
        ttk::VertexScalarFieldName);
    const idSuperArc numberOfSuperArcs = tree->getNumberOfSuperArcs();
    // #pragma omp for
    for(idSuperArc arcId = 0; arcId < numberOfSuperArcs; ++arcId) {
      vertData.fillArrayPoint(
        arcId, ftmTree_[cc], triangulation_[cc], idMapper, params_);
    }
  }

  vtkPointData *pointData = outputSegmentation->GetPointData();
  vertData.addArray(pointData, params_);

  return 0;
}

int ttkFTMTree::getSkeletonArcs(vtkUnstructuredGrid *outputSkeletonArcs) {
  vtkSmartPointer<vtkUnstructuredGrid> skeletonArcs
    = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  ttk::ftm::ArcData arcData;
  arcData.init(ftmTree_, params_);

  const int samplingLevel = params_.samplingLvl;
  for(int cc = 0; cc < nbCC_; cc++) {
    FTMTree_MT *tree = ftmTree_[cc].tree.getTree(GetTreeType());

    const SimplexId numberOfSuperArcs = tree->getNumberOfSuperArcs();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!numberOfSuperArcs) {
      cerr << "[ttkFTMTree] Error : tree has no super arcs." << endl;
      return -2;
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

  return 0;
}

int ttkFTMTree::getSkeletonNodes(vtkUnstructuredGrid *outputSkeletonNodes) {
  vtkSmartPointer<vtkUnstructuredGrid> skeletonNodes
    = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  NodeData nodeData;
  nodeData.init(ftmTree_, params_);
  nodeData.setScalarType(inputScalars_[0]->GetDataType());

  for(int cc = 0; cc < nbCC_; cc++) {
    FTMTree_MT *tree = ftmTree_[cc].tree.getTree(GetTreeType());
    vtkDataArray *idMapper
      = connected_components_[cc]->GetPointData()->GetArray(
        ttk::VertexScalarFieldName);

    const idNode numberOfNodes = tree->getNumberOfNodes();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!numberOfNodes) {
      cerr << "[ttkFTMTree] Error : tree has no nodes." << endl;
      return -2;
    }
#endif

    for(idNode nodeId = 0; nodeId < numberOfNodes; ++nodeId) {
      const Node *node = tree->getNode(nodeId);
#ifndef TTK_ENABLE_KAMIKAZE
      if(!node) {
        cerr << "[ttkFTMTree] Error : node " << nodeId << " is null." << endl;
        return -7;
      }
#endif
      const SimplexId local_vertId = node->getVertexId();
      float point[3];
      triangulation_[cc]->getVertexPoint(
        local_vertId, point[0], point[1], point[2]);
      const SimplexId nextPoint = points->InsertNextPoint(point);
      nodeData.fillArrayPoint(
        nextPoint, nodeId, ftmTree_[cc], idMapper, triangulation_[cc], params_);
    }
  }

  skeletonNodes->SetPoints(points);
  vtkPointData *pointData = skeletonNodes->GetPointData();
  nodeData.addArray(pointData, params_);
  outputSkeletonNodes->ShallowCopy(skeletonNodes);

  return 0;
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

int ttkFTMTree::setupTriangulation() {
  triangulation_.resize(nbCC_);
  ftmTree_.resize(nbCC_);

  for(int cc = 0; cc < nbCC_; cc++) {
    triangulation_[cc]
      = ttkTriangulation::getTriangulation(connected_components_[cc]);
#ifndef TTK_ENABLE_KAMIKAZE
    if(!triangulation_[cc]) {
      cerr
        << "[ttkFTMTree] Error : ttkTriangulation::getTriangulation() is null."
        << endl;
      return -1;
    }
#endif

    triangulation_[cc]->setPeriodicBoundaryConditions(
      PeriodicBoundaryConditions);
    triangulation_[cc]->setWrapper(this);
    ftmTree_[cc].tree.setDebugLevel(debugLevel_);
    ftmTree_[cc].tree.setThreadNumber(threadNumber_);
    ftmTree_[cc].tree.setupTriangulation(triangulation_[cc]);

    hasUpdatedMesh_ = ttkTriangulation::hasChangedConnectivity(
      triangulation_[cc], connected_components_[cc], this);

#ifndef TTK_ENABLE_KAMIKAZE
    if(triangulation_[cc]->isEmpty()) {
      cerr << "[ttkFTMTree] Error : ttkTriangulation on connected component"
           << cc << " allocation problem." << endl;
      return -1;
    }
#endif
  }
  return 0;
}

// protected

ttkFTMTree::ttkFTMTree()
  : ScalarField{}, ForceInputOffsetScalarField{false},
    InputOffsetScalarFieldName{ttk::OffsetScalarFieldName}, ScalarFieldId{},
    OffsetFieldId{-1}, PeriodicBoundaryConditions{false}, params_{},
    triangulation_{}, inputScalars_{}, offsets_{}, hasUpdatedMesh_{} {
  SetSuperArcSamplingLevel(0);
  SetWithNormalize(true);
  SetWithAdvStats(true);
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
  UseAllCores = true;
}

ttkFTMTree::~ttkFTMTree() {
}

void ttkFTMTree::identify(vtkDataSet *ds) const {
  vtkSmartPointer<ttkSimplexIdTypeArray> identifiers
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  const SimplexId nbPoints = ds->GetNumberOfPoints();
  identifiers->SetName(ttk::VertexScalarFieldName);
  identifiers->SetNumberOfComponents(1);
  identifiers->SetNumberOfTuples(nbPoints);

  for(SimplexId i = 0; i < nbPoints; i++) {
    identifiers->SetTuple1(i, i);
  }

  ds->GetPointData()->AddArray(identifiers);
}
