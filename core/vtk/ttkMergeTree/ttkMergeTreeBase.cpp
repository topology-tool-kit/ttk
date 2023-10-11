#include <ttkMergeTreeBase.h>
#include <ttkUtils.h>

ttkMergeTreeBase::ttkMergeTreeBase() {
  this->setDebugMsgPrefix("MergeTreeBase");
}

int ttkMergeTreeBase::addCompleteSkeletonArc(const ttk::ftm::idSuperArc arcId,
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

int ttkMergeTreeBase::addDirectSkeletonArc(const ttk::ftm::idSuperArc arcId,
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

int ttkMergeTreeBase::addSampledSkeletonArc(const ttk::ftm::idSuperArc arcId,
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

int ttkMergeTreeBase::getSegmentation(vtkDataSet *outputSegmentation) {
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

int ttkMergeTreeBase::getSkeletonArcs(vtkUnstructuredGrid *outputSkeletonArcs) {
  vtkNew<vtkUnstructuredGrid> skeletonArcs{};
  vtkNew<vtkPoints> const points{};

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

int ttkMergeTreeBase::getSkeletonNodes(
  vtkUnstructuredGrid *outputSkeletonNodes) {
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

  ttkUtils::CellVertexFromPoints(skeletonNodes, points);

  vtkPointData *pointData = skeletonNodes->GetPointData();
  nodeData.addArray(pointData, params_);

  outputSkeletonNodes->ShallowCopy(skeletonNodes);

  return 1;
}

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
void ttkMergeTreeBase::printCSVStats() {
  using namespace ttk;

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

void ttkMergeTreeBase::printCSVTree(
  const ttk::ftm::FTMTree_MT *const tree) const {
  using namespace ttk::ftm;

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

// protected

void ttkMergeTreeBase::identify(vtkDataSet *ds) const {
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
