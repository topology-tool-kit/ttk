#include "ttkFTRGraph.h"

// only used on the cpp
#include <vtkConnectivityFilter.h>
#include <vtkDataObject.h>
#include <vtkThreshold.h>

using namespace ttk::ftr;

vtkStandardNewMacro(ttkFTRGraph);

int ttkFTRGraph::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}

int ttkFTRGraph::FillOutputPortInformation(int port, vtkInformation *info) {
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
    std::cerr << "NULL NODES IN SKELETON " << graph.printArc(arcId)
              << std::endl;
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
    std::cerr << "NULL NODES IN SKELETON " << graph.printArc(arcId)
              << std::endl;
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
    std::cerr << "NULL NODES IN SKELETON " << graph.printArc(arcId)
              << std::endl;
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

int ttkFTRGraph::doIt(std::vector<vtkDataSet *> &inputs,
                      std::vector<vtkDataSet *> &outputs) {
  ttk::Memory m;

  // Input
  mesh_ = inputs[0];

  // Output
  vtkUnstructuredGrid *outputSkeletonNodes
    = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  vtkUnstructuredGrid *outputSkeletonArcs
    = vtkUnstructuredGrid::SafeDownCast(outputs[1]);
  vtkDataSet *outputSegmentation = outputs[2];
  outputSegmentation->ShallowCopy(inputs[0]);

  // Skeleton
  Graph graph;

  if(setupTriangulation()) {
#ifndef TTK_ENABLE_KAMIKAZE
    cerr << "[ttkFTRGraph] Error : wrong triangulation." << endl;
    return -1;
#endif
  }

  // gives parameters to the params_ structure
  params_.debugLevel = debugLevel_;
  params_.threadNumber = threadNumber_;

  // Scalar field related
  getScalars();
  getOffsets();

  // compute graph
  switch(inputScalars_->GetDataType()) {
    vtkTemplateMacro({
      ttk::ftr::FTRGraph<VTK_TT> ftrGraph_(triangulation_);
      // common parameters
      ftrGraph_.setWrapper(this);
      ftrGraph_.setParams(params_);
      // reeb graph parameters
      ftrGraph_.setScalars(inputScalars_->GetVoidPointer(0));
      ftrGraph_.setVertexSoSoffsets(&offsets_);
      // build
      {
        std::stringstream msg;
        msg << "[ttkFTRGraph] Starting computation on field: " << ScalarField
            << std::endl;
        dMsg(cout, msg.str(), infoMsg);
      }
      ftrGraph_.build();
      // get output
      graph = std::move(ftrGraph_.extractOutputGraph());
    });
  }

  UpdateProgress(0.50);

  // Construct output
  if(getSkeletonNodes(graph, outputSkeletonNodes)) {
#ifndef TTK_ENABLE_KAMIKAZE
    cerr << "[ttkFTRGraph] Error : wrong properties on skeleton nodes." << endl;
    return -7;
#endif
  }

  if(getSkeletonArcs(graph, outputSkeletonArcs)) {
#ifndef TTK_ENABLE_KAMIKAZE
    cerr << "[ttkFTRGraph] Error : wrong properties on skeleton arcs." << endl;
    return -8;
#endif
  }

  if(GetWithSegmentation() && getSegmentation(graph, outputSegmentation)) {
#ifndef TTK_ENABLE_KAMIKAZE
    cerr << "[ttkFTRGraph] Error : wrong properties on segmentation." << endl;
    return -9;
#endif
  }

  UpdateProgress(1);

  {
    std::stringstream msg;
    msg << "[ttkFTRGraph] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkFTRGraph::getOffsets() {

  vtkDataArray *inputOffsets;

  inputOffsets = mesh_->GetPointData()->GetArray(ttk::OffsetScalarFieldName);
  if(inputOffsets) {
    InputOffsetScalarFieldName = inputOffsets->GetName();
    UseInputOffsetScalarField = true;
  } else {
    inputOffsets = mesh_->GetPointData()->GetArray(OffsetFieldId);
    if(inputOffsets) {
      InputOffsetScalarFieldName = inputOffsets->GetName();
      UseInputOffsetScalarField = true;
    }
  }

  const idVertex numberOfVertices = mesh_->GetNumberOfPoints();

  if(UseInputOffsetScalarField and InputOffsetScalarFieldName.length()) {
    inputOffsets
      = mesh_->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
    offsets_.resize(numberOfVertices);
    for(int i = 0; i < numberOfVertices; i++) {
      offsets_[i] = inputOffsets->GetTuple1(i);
    }
  } else {
    if(offsets_.size()) {
      // don't keep an out-dated offset array
      offsets_.clear();
    }

    if(offsets_.empty()) {
      offsets_.resize(numberOfVertices);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  schedule(static, numberOfVertices / threadNumber_)
#endif
      for(int i = 0; i < numberOfVertices; i++) {
        offsets_[i] = i;
      }
    }
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(offsets_.empty()) {
    cerr << "[ttkFTRGraph] Error : wrong input offset scalar field " << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkFTRGraph::getScalars() {
  vtkPointData *pointData = mesh_->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    cerr << "[ttkFTRGraph] Error : input has no point data." << endl;
    return -1;
  }
#endif

  if(ScalarField.length()) {
    inputScalars_ = pointData->GetArray(ScalarField.data());
  } else {
    inputScalars_ = pointData->GetArray(ScalarFieldId);
    if(inputScalars_)
      ScalarField = inputScalars_->GetName();
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars_) {
    cerr << "[ttkFTRGraph] Error : input scalar field pointer is null." << endl;
    return -3;
  }
#endif

  return 0;
}

int ttkFTRGraph::getSegmentation(const ttk::ftr::Graph &graph,
                                 vtkDataSet *outputSegmentation) {
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
  vtkSmartPointer<vtkUnstructuredGrid> arcs
    = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

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
  vtkSmartPointer<vtkUnstructuredGrid> nodes
    = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

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

int ttkFTRGraph::setupTriangulation() {
  triangulation_ = ttkTriangulation::getTriangulation(mesh_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_) {
    cerr
      << "[ttkFTRGraph] Error : ttkTriangulation::getTriangulation() is null."
      << endl;
    return -1;
  }
#endif

  triangulation_->setWrapper(this);
  hasUpdatedMesh_
    = ttkTriangulation::hasChangedConnectivity(triangulation_, mesh_, this);

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_->isEmpty()) {
    cerr << "[ttkFTRGraph] Error : ttkTriangulation on connected component "
            "allocation problem."
         << endl;
    return -1;
  }
#endif
  return 0;
}

// protected

ttkFTRGraph::ttkFTRGraph()
  : ScalarField{}, UseInputOffsetScalarField{}, InputOffsetScalarFieldName{},
    ScalarFieldId{}, OffsetFieldId{-1}, params_{}, mesh_{}, triangulation_{},
    inputScalars_{}, offsets_{}, hasUpdatedMesh_{} {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
}

ttkFTRGraph::~ttkFTRGraph() {
}

void ttkFTRGraph::identify(vtkDataSet *ds) const {
  vtkSmartPointer<vtkIntArray> identifiers
    = vtkSmartPointer<vtkIntArray>::New();
  const vtkIdType nbPoints = ds->GetNumberOfPoints();
  identifiers->SetName("VertexIdentifier");
  identifiers->SetNumberOfComponents(1);
  identifiers->SetNumberOfTuples(nbPoints);

  for(int i = 0; i < nbPoints; i++) {
    identifiers->SetTuple1(i, i);
  }

  ds->GetPointData()->AddArray(identifiers);
}
