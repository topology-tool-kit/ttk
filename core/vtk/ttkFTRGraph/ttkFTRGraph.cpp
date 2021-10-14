#include "ttkFTRGraph.h"
#include <ttkMacros.h>
#include <ttkUtils.h>

// only used on the cpp
#include <vtkConnectivityFilter.h>
#include <vtkDataObject.h>
#include <vtkInformation.h>

using namespace ttk::ftr;

vtkStandardNewMacro(ttkFTRGraph);

ttkFTRGraph::ttkFTRGraph() {
  this->setDebugMsgPrefix("FTRGraph");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(3);
}

int ttkFTRGraph::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkFTRGraph::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkFTRGraph::addCompleteSkeletonArc(const Graph &graph,
                                        const idSuperArc arcId,
                                        vtkPoints *points,
                                        vtkUnstructuredGrid *skeletonArcs,
                                        ttk::ftr::ArcData &arcData) {
  const SuperArc &arc = graph.getArc(arcId);
  float pointCoord[3];
  const idNode upNodeId = arc.getUpNodeId();
  const idNode downNodeId = arc.getDownNodeId();

#ifndef TTK_ENABLE_KAMIKAZE
  if(upNodeId == nullNode || downNodeId == nullNode) {
    this->printErr("NULL NODES IN SKELETON " + graph.printArc(arcId));
    return 1;
  }
#endif

  const idVertex upVertId = graph.getNode(upNodeId).getVertexIdentifier();
  const idVertex downVertId = graph.getNode(downNodeId).getVertexIdentifier();

  vtkIdType pointIds[2];
  if(arcData.points.count(downVertId)) {
    pointIds[0] = arcData.points[downVertId];
  } else {
    triangulation_->getVertexPoint(
      downVertId, pointCoord[0], pointCoord[1], pointCoord[2]);
    pointIds[0] = points->InsertNextPoint(pointCoord);
    arcData.points.emplace(downVertId, pointIds[0]);
    arcData.setPointInfo(graph, arcId, pointIds[0]);
  }

  for(const idVertex regV : arc.segmentation()) {

    // for regular mask to be set by the boundary vertices correctly
    if(regV == downVertId || regV == upVertId)
      continue;

    if(arcData.points.count(regV)) {
      pointIds[1] = arcData.points[regV];
    } else {
      triangulation_->getVertexPoint(
        regV, pointCoord[0], pointCoord[1], pointCoord[2]);
      pointIds[1] = points->InsertNextPoint(pointCoord);
      arcData.points.emplace(regV, pointIds[1]);
      arcData.setPointInfo(graph, arcId, pointIds[1], true);
    }
    if(pointIds[0] != pointIds[1]) {
      vtkIdType nextCellId
        = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
      arcData.setArcInfo(graph, arcId, nextCellId);
    }
    pointIds[0] = pointIds[1];
  }

  if(arcData.points.count(upVertId)) {
    pointIds[1] = arcData.points[upVertId];
  } else {
    triangulation_->getVertexPoint(
      upVertId, pointCoord[0], pointCoord[1], pointCoord[2]);
    pointIds[1] = points->InsertNextPoint(pointCoord);
    arcData.points.emplace(upVertId, pointIds[1]);
    arcData.setPointInfo(graph, arcId, pointIds[1]);
  }
  if(pointIds[0] != pointIds[1]) {
    vtkIdType nextCellId = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
    arcData.setArcInfo(graph, arcId, nextCellId);
  }

  return 0;
}

int ttkFTRGraph::addDirectSkeletonArc(const Graph &graph,
                                      const idSuperArc arcId,
                                      vtkPoints *points,
                                      vtkUnstructuredGrid *skeletonArcs,
                                      ttk::ftr::ArcData &arcData) {
  float pointCoord[3];
  const idNode upNodeId = graph.getArc(arcId).getUpNodeId();
  const idNode downNodeId = graph.getArc(arcId).getDownNodeId();

#ifndef TTK_ENABLE_KAMIKAZE
  if(upNodeId == nullNode || downNodeId == nullNode) {
    this->printErr("NULL NODES IN SKELETON " + graph.printArc(arcId));
    return 1;
  }
#endif

  const idVertex upVertId = graph.getNode(upNodeId).getVertexIdentifier();
  const idVertex downVertId = graph.getNode(downNodeId).getVertexIdentifier();

  // Mesh skeleton
  vtkIdType pointIds[2];
  if(arcData.points.count(downVertId)) {
    pointIds[0] = arcData.points[downVertId];
  } else {
    triangulation_->getVertexPoint(
      downVertId, pointCoord[0], pointCoord[1], pointCoord[2]);
    pointIds[0] = points->InsertNextPoint(pointCoord);
    arcData.points.emplace(downVertId, pointIds[0]);
    arcData.setPointInfo(graph, arcId, pointIds[0]);
  }
  if(arcData.points.count(upVertId)) {
    pointIds[1] = arcData.points[upVertId];
  } else {
    triangulation_->getVertexPoint(
      upVertId, pointCoord[0], pointCoord[1], pointCoord[2]);
    pointIds[1] = points->InsertNextPoint(pointCoord);
    arcData.points.emplace(upVertId, pointIds[1]);
    arcData.setPointInfo(graph, arcId, pointIds[1]);
  }
  vtkIdType nextCellId = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);

  // Scalar arrays
  arcData.setArcInfo(graph, arcId, nextCellId);
  return 0;
}

int ttkFTRGraph::addSampledSkeletonArc(const Graph &graph,
                                       const idSuperArc arcId,
                                       vtkPoints *points,
                                       vtkUnstructuredGrid *skeletonArcs,
                                       ttk::ftr::ArcData &arcData) {
  const SuperArc &arc = graph.getArc(arcId);
  float pointCoord[3];
  const idNode upNodeId = arc.getUpNodeId();
  const idNode downNodeId = arc.getDownNodeId();

  const double scalarMin
    = inputScalars_->GetTuple1(graph.getNode(downNodeId).getVertexIdentifier());
  const double scalarMax
    = inputScalars_->GetTuple1(graph.getNode(upNodeId).getVertexIdentifier());
  const double delta = (scalarMax - scalarMin) / (params_.samplingLvl + 1);

#ifndef TTK_ENABLE_KAMIKAZE
  if(upNodeId == nullNode || downNodeId == nullNode) {
    this->printErr("NULL NODES IN SKELETON " + graph.printArc(arcId));
    return 1;
  }
#endif

  const idVertex upVertId = graph.getNode(upNodeId).getVertexIdentifier();
  const idVertex downVertId = graph.getNode(downNodeId).getVertexIdentifier();

  vtkIdType pointIds[2];
  if(arcData.points.count(downVertId)) {
    pointIds[0] = arcData.points[downVertId];
  } else {
    triangulation_->getVertexPoint(
      downVertId, pointCoord[0], pointCoord[1], pointCoord[2]);
    pointIds[0] = points->InsertNextPoint(pointCoord);
    arcData.points.emplace(downVertId, pointIds[0]);
    arcData.setPointInfo(graph, arcId, pointIds[0]);
  }

  float sum[3]{};
  int chunk = 0;
  double scalarLimit = scalarMin + delta;

  for(const idVertex regV : arc.segmentation()) {

    // for regular mask to be set by the boundary vertices correctly
    if(regV == downVertId || regV == upVertId)
      continue;

    triangulation_->getVertexPoint(
      regV, pointCoord[0], pointCoord[1], pointCoord[2]);
    const double scalar = inputScalars_->GetTuple1(regV);
    if(scalar < scalarLimit) {
      sum[0] += pointCoord[0];
      sum[1] += pointCoord[1];
      sum[2] += pointCoord[2];
      ++chunk;
    } else {
      if(chunk) {
        sum[0] /= chunk;
        sum[1] /= chunk;
        sum[2] /= chunk;

        // do not use memorized points as even if a point already exisit it is
        // from another vertex and should not be used
        pointIds[1] = points->InsertNextPoint(sum);
        arcData.points.emplace(regV, pointIds[1]);
        arcData.setPointInfo(graph, arcId, pointIds[1], true);

        if(pointIds[0] != pointIds[1]) {
          vtkIdType nextCellId
            = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
          arcData.setArcInfo(graph, arcId, nextCellId);
        }

        pointIds[0] = pointIds[1];
      }

      // reset for next iteration
      scalarLimit += delta;
      sum[0] = 0;
      sum[1] = 0;
      sum[2] = 0;
      chunk = 0;
    }
  }

  if(arcData.points.count(upVertId)) {
    pointIds[1] = arcData.points[upVertId];
  } else {
    triangulation_->getVertexPoint(
      upVertId, pointCoord[0], pointCoord[1], pointCoord[2]);
    pointIds[1] = points->InsertNextPoint(pointCoord);
    arcData.points.emplace(upVertId, pointIds[1]);
    arcData.setPointInfo(graph, arcId, pointIds[1]);
  }
  if(pointIds[0] != pointIds[1]) {
    vtkIdType nextCellId = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
    arcData.setArcInfo(graph, arcId, nextCellId);
  }

  return 0;
}

template <typename VTK_TT, typename TTK_TT>
int ttkFTRGraph::dispatch(Graph &graph) {
  ttk::ftr::FTRGraph<VTK_TT, TTK_TT> ftrGraph_(
    static_cast<TTK_TT *>(triangulation_->getData()));

  // common parameters
  ftrGraph_.setParams(params_);
  // reeb graph parameters
  ftrGraph_.setScalars(ttkUtils::GetVoidPointer(inputScalars_));
  // TODO: SimplexId -> template to int + long long?
  ftrGraph_.setVertexSoSoffsets(
    static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(offsets_)));
  // build
  this->printMsg("Starting computation on field: "
                 + std::string{inputScalars_->GetName()});
  ftrGraph_.build();
  // get output
  graph = std::move(ftrGraph_.extractOutputGraph());

  return 0;
}

int ttkFTRGraph::RequestData(vtkInformation *ttkNotUsed(request),
                             vtkInformationVector **inputVector,
                             vtkInformationVector *outputVector) {

  // Input
  mesh_ = vtkDataSet::GetData(inputVector[0]);

  // Output
  auto outputSkeletonNodes = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputSkeletonArcs = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto outputSegmentation = vtkDataSet::GetData(outputVector, 2);

  // Skeleton
  Graph graph;

  triangulation_ = ttkAlgorithm::GetTriangulation(mesh_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_) {
    this->printErr("ttkTriangulation::getTriangulation() is null.");
    return -1;
  }
  if(triangulation_->isEmpty()) {
    this->printErr(
      "ttkTriangulation on connected component allocation problem.");
    return -1;
  }
#endif

  // gives parameters to the params_ structure
  params_.debugLevel = debugLevel_;
  params_.threadNumber = threadNumber_;

  // Scalar field related
  inputScalars_ = this->GetInputArrayToProcess(0, inputVector);
  if(inputScalars_ == nullptr) {
    this->printErr("input scalar field pointer is null.");
    return -3;
  }

  offsets_ = this->GetOrderArray(mesh_, 0, 1, ForceInputOffsetScalarField);

  // compute graph
  ttkVtkTemplateMacro(inputScalars_->GetDataType(), triangulation_->getType(),
                      (dispatch<VTK_TT, TTK_TT>(graph)));

  UpdateProgress(0.50);

  // Construct output
  if(getSkeletonNodes(graph, outputSkeletonNodes)) {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printErr("wrong properties on skeleton nodes.");
    return -7;
#endif
  }

  if(getSkeletonArcs(graph, outputSkeletonArcs)) {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printErr("wrong properties on skeleton arcs.");
    return -8;
#endif
  }

  if(GetWithSegmentation() && getSegmentation(graph, outputSegmentation)) {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printErr("wrong properties on segmentation.");
    return -9;
#endif
  }

  UpdateProgress(1);

  return 1;
}

int ttkFTRGraph::getSegmentation(const ttk::ftr::Graph &graph,
                                 vtkDataSet *outputSegmentation) {
  outputSegmentation->ShallowCopy(mesh_);

  const idVertex numberOfVertices = mesh_->GetNumberOfPoints();
  ttk::ftr::VertData vertData(numberOfVertices);

  // TODO parallel
  for(idVertex v = 0; v < numberOfVertices; ++v) {
    vertData.setVertexInfo(graph, v);
  }

  vertData.addArrays(outputSegmentation, params_);

  return 0;
}

int ttkFTRGraph::getSkeletonArcs(const ttk::ftr::Graph &graph,
                                 vtkUnstructuredGrid *outputSkeletonArcs) {
  const idSuperArc nbArcs = graph.getNumberOfArcs();

  idSuperArc nbFinArc = 0;
  switch(params_.samplingLvl) {
    case -1:
      // loops may create more arcs
      nbFinArc = triangulation_->getNumberOfVertices() * 1.5;
      break;
    case 0:
      nbFinArc = graph.getNumberOfVisibleArcs();
      break;
    default:
      nbFinArc = graph.getNumberOfVisibleArcs() * (params_.samplingLvl + 1);
      break;
  }

  ttk::ftr::ArcData arcData(nbFinArc);
  vtkNew<vtkUnstructuredGrid> arcs{};
  vtkNew<vtkPoints> points{};

  for(idSuperArc arcId = 0; arcId < nbArcs; ++arcId) {
    if(!graph.getArc(arcId).isVisible())
      continue;

    switch(params_.samplingLvl) {
      case -1:
        addCompleteSkeletonArc(graph, arcId, points, arcs, arcData);
        break;
      case 0:
        addDirectSkeletonArc(graph, arcId, points, arcs, arcData);
        break;
      default:
        addSampledSkeletonArc(graph, arcId, points, arcs, arcData);
        break;
    }
  }

  arcs->SetPoints(points);
  arcData.addArrays(arcs, params_);
  outputSkeletonArcs->ShallowCopy(arcs);

  return 0;
}

int ttkFTRGraph::getSkeletonNodes(const Graph &graph,
                                  vtkUnstructuredGrid *outputSkeletonNodes) {
  const idNode nbNodes = graph.getNumberOfNodes();

  ttk::ftr::NodeData nodeData(nbNodes);
  vtkNew<vtkUnstructuredGrid> nodes{};
  vtkNew<vtkPoints> points{};

  for(idNode nodeId = 0; nodeId < nbNodes; nodeId++) {
    const ttk::ftr::idVertex vertId
      = graph.getNode(nodeId).getVertexIdentifier();
    float point[3];
    triangulation_->getVertexPoint(vertId, point[0], point[1], point[2]);
    points->InsertNextPoint(point);

    double scalar = inputScalars_->GetTuple1(vertId);
    nodeData.addNode(graph, nodeId, scalar);
  }

  nodes->SetPoints(points);
  nodeData.addArrays(nodes->GetPointData(), params_);

  outputSkeletonNodes->ShallowCopy(nodes);

  return 0;
}

// protected

void ttkFTRGraph::identify(vtkDataSet *ds) const {
  vtkNew<vtkIntArray> identifiers{};
  const vtkIdType nbPoints = ds->GetNumberOfPoints();
  identifiers->SetName("VertexIdentifier");
  identifiers->SetNumberOfComponents(1);
  identifiers->SetNumberOfTuples(nbPoints);

  for(int i = 0; i < nbPoints; i++) {
    identifiers->SetTuple1(i, i);
  }

  ds->GetPointData()->AddArray(identifiers);
}
