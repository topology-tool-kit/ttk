#include <ttkMergeTree.h>

#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkSignedCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertex.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <PathCompression.h>

vtkStandardNewMacro(ttkMergeTree);

ttkMergeTree::ttkMergeTree() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(3);
}

ttkMergeTree::~ttkMergeTree() = default;

int ttkMergeTree::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkMergeTree::FillOutputPortInformation(int port, vtkInformation *info) {
  switch(port) {
    case 0:
    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      return 1;
    case 2:
      info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
      return 1;
    default:
      return 0;
  }
}

// TODO: future optimization: flag critical vertices based on if they appear in
// the merge tree, count the number of critical points and allocate the arrays
// accordingly
template <class triangulationType>
int ttkMergeTree::getMergeTree(vtkUnstructuredGrid *outputSkeletonArcs,
                               std::vector<ExTreeM::Branch> &mergeTree,
                               vtkDataArray *inputScalars,
                               const triangulationType *triangulation) {
  vtkNew<vtkUnstructuredGrid> skeletonArcs{};
  ttk::SimplexId pointIds[2];
  ttk::SimplexId pointOrders[2];
  vtkNew<vtkPoints> points{};
  vtkNew<vtkLongLongArray> data{};
  data->SetNumberOfComponents(1);
  data->SetName("Order");
  vtkNew<vtkIdTypeArray> gIdArray{};
  gIdArray->SetNumberOfComponents(1);
  gIdArray->SetName("GlobalPointIds");
  vtkNew<vtkIdTypeArray> downId{};
  downId->SetNumberOfComponents(1);
  downId->SetName("downNodeId");
  vtkNew<vtkIdTypeArray> upId{};
  upId->SetNumberOfComponents(1);
  upId->SetName("upNodeId");
  float point[3];
  vtkSmartPointer<vtkDataArray> const scalarArray
    = vtkSmartPointer<vtkDataArray>::Take(inputScalars->NewInstance());
  scalarArray->SetNumberOfComponents(1);
  scalarArray->SetName("Scalar");
  std::map<ttk::SimplexId, ttk::SimplexId> addedPoints;
  ttk::SimplexId currentId = 0;
  for(auto const &b : mergeTree) {
    auto &vertices = b.vertices;
    for(size_t p = 0; p < vertices.size() - 1; p++) {
      pointIds[0] = vertices[p].second;
      pointIds[1] = vertices[p + 1].second;
      pointOrders[0] = vertices[p].first;
      pointOrders[1] = vertices[p + 1].first;
      // add each point only once to the vtkPoints
      // addedPoints.insert(x).second inserts x and is true if x was not in
      // addedPoints beforehand
      if(addedPoints.insert({pointIds[0], currentId}).second) {
        triangulation->getVertexPoint(
          pointIds[0], point[0], point[1], point[2]);
        points->InsertNextPoint(point);
        data->InsertNextTuple1(pointOrders[0]);
        gIdArray->InsertNextTuple1(pointIds[0]);
        scalarArray->InsertNextTuple1(inputScalars->GetTuple1(pointIds[0]));
        currentId++;
      }
      if(addedPoints.insert({pointIds[1], currentId}).second) {
        triangulation->getVertexPoint(
          pointIds[1], point[0], point[1], point[2]);
        points->InsertNextPoint(point);
        data->InsertNextTuple1(pointOrders[1]);
        gIdArray->InsertNextTuple1(pointIds[1]);
        scalarArray->InsertNextTuple1(inputScalars->GetTuple1(pointIds[1]));
        currentId++;
      }
      downId->InsertNextTuple1(pointIds[0]);
      upId->InsertNextTuple1(pointIds[1]);
      vtkIdType pointIdsASVTKIDTYPE[2];
      pointIdsASVTKIDTYPE[0] = addedPoints.at(pointIds[0]);
      pointIdsASVTKIDTYPE[1] = addedPoints.at(pointIds[1]);
      skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIdsASVTKIDTYPE);
    }
  }
  skeletonArcs->SetPoints(points);
  outputSkeletonArcs->ShallowCopy(skeletonArcs);
  outputSkeletonArcs->GetPointData()->AddArray(data);
  outputSkeletonArcs->GetPointData()->AddArray(gIdArray);
  outputSkeletonArcs->GetPointData()->AddArray(scalarArray);
  outputSkeletonArcs->GetCellData()->AddArray(upId);
  outputSkeletonArcs->GetCellData()->AddArray(downId);

  return 1;
}

template <class triangulationType>
int ttkMergeTree::getMergeTreePoints(
  vtkUnstructuredGrid *outputSkeletonNodes,
  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &persistencePairs,
  vtkDataArray *inputScalars,
  const triangulationType *triangulation) {
  vtkNew<vtkUnstructuredGrid> skeletonNodes{};
  vtkNew<vtkPoints> points{};
  vtkNew<vtkCellArray> cells{};
  vtkNew<vtkLongLongArray> gIdArray{};
  gIdArray->SetNumberOfComponents(1);
  gIdArray->SetName("VertexId");
  vtkSmartPointer<vtkDataArray> const scalarArray
    = vtkSmartPointer<vtkDataArray>::Take(inputScalars->NewInstance());
  scalarArray->SetNumberOfComponents(1);
  scalarArray->SetName("Scalar");
  float point[3];
  long long pointId = 0;

  for(auto const &pair : persistencePairs) {
    triangulation->getVertexPoint(pair.first, point[0], point[1], point[2]);
    points->InsertNextPoint(point);
    gIdArray->InsertNextTuple1(pair.first);
    scalarArray->InsertNextTuple1(inputScalars->GetTuple1(pair.first));
    triangulation->getVertexPoint(pair.second, point[0], point[1], point[2]);
    points->InsertNextPoint(point);
    gIdArray->InsertNextTuple1(pair.second);
    scalarArray->InsertNextTuple1(inputScalars->GetTuple1(pair.second));
    skeletonNodes->InsertNextCell(VTK_VERTEX, 1, &pointId);
    pointId++;
    skeletonNodes->InsertNextCell(VTK_VERTEX, 1, &pointId);
    pointId++;
  }
  skeletonNodes->SetPoints(points);
  outputSkeletonNodes->ShallowCopy(skeletonNodes);
  outputSkeletonNodes->GetPointData()->AddArray(gIdArray);
  outputSkeletonNodes->GetPointData()->AddArray(scalarArray);

  return 1;
}

int ttkMergeTree::RequestData(vtkInformation *,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {
#ifndef TTK_ENABLE_OPENMP
  printMsg(ttk::debug::Separator::L2);
  this->printErr("ExTreeM requires OpenMP.");
  printMsg(ttk::debug::Separator::L2);
  return 0;
#endif
  // Get the input
  auto input = vtkImageData::GetData(inputVector[0]);
  if(!input) {
    this->printErr("Unable to retrieve input data object.");
    return 0;
  }
  const size_t nVertices = input->GetNumberOfPoints();

  // Get triangulation of the input object
  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  triangulation->preconditionVertexNeighbors();

  // Get input array
  auto scalarArray = this->GetInputArrayToProcess(0, inputVector);
  if(!scalarArray) {
    this->printErr("Unable to retrieve scalar array.");
    return 0;
  }

  // Order Array
  auto orderArray = this->GetOrderArray(input, 0);
  auto orderArrayData = ttkUtils::GetPointer<ttk::SimplexId>(orderArray);

  auto segmentation = vtkDataSet::GetData(outputVector, 2);
  segmentation->ShallowCopy(input);
  auto segmentationPD = segmentation->GetPointData();

  // enforce that ascending and descending manifolds exist
  if(!segmentationPD->HasArray(ttk::MorseSmaleAscendingName)
     || !segmentationPD->HasArray(ttk::MorseSmaleDescendingName)) {
    printMsg(ttk::debug::Separator::L2);
    this->printWrn("TIP: run `ttkPathCompression` first");
    this->printWrn("for improved performances :)");
    printMsg(ttk::debug::Separator::L2);
    bool doAscend = false;
    bool doDescend = false;
    if(!segmentationPD->HasArray(ttk::MorseSmaleAscendingName)) {
      doAscend = true;
      auto ascendingManifold = vtkSmartPointer<vtkIdTypeArray>::New();
      ascendingManifold->SetNumberOfComponents(1);
      ascendingManifold->SetNumberOfTuples(nVertices);
      ascendingManifold->SetName(ttk::MorseSmaleAscendingName);
      segmentationPD->AddArray(ascendingManifold);
    }
    if(!segmentationPD->HasArray(ttk::MorseSmaleDescendingName)) {
      doDescend = true;
      auto descendingManifold = vtkSmartPointer<vtkIdTypeArray>::New();
      descendingManifold->SetNumberOfComponents(1);
      descendingManifold->SetNumberOfTuples(nVertices);
      descendingManifold->SetName(ttk::MorseSmaleDescendingName);
      segmentationPD->AddArray(descendingManifold);
    }
    ttk::PathCompression subModule;
    subModule.setThreadNumber(this->threadNumber_);
    subModule.setDebugLevel(this->debugLevel_);
    // only compute the segmentation which doesn't exist (maybe both)
    subModule.setComputeSegmentation(doAscend, doDescend, false);

    ttk::PathCompression::OutputSegmentation om{
      ttkUtils::GetPointer<ttk::SimplexId>(
        segmentationPD->GetArray(ttk::MorseSmaleAscendingName)),
      ttkUtils::GetPointer<ttk::SimplexId>(
        segmentationPD->GetArray(ttk::MorseSmaleDescendingName)),
      nullptr};

    int status = 0;
    ttkTypeMacroT(triangulation->getType(),
                  (status = subModule.execute<T0>(
                     om, ttkUtils::GetPointer<const ttk::SimplexId>(orderArray),
                     *(T0 *)triangulation->getData())));
    if(status != 0)
      return 0;
  }

  auto ascendingManifold
    = segmentationPD->GetArray(ttk::MorseSmaleAscendingName);
  auto descendingManifold
    = segmentationPD->GetArray(ttk::MorseSmaleDescendingName);

  vtkNew<ttkSimplexIdTypeArray> segmentationId{};
  segmentationId->SetNumberOfComponents(1);
  segmentationId->SetNumberOfTuples(nVertices);
  segmentationId->SetName("SegmentationId");

  vtkNew<vtkUnsignedCharArray> isLeaf{};
  isLeaf->SetNumberOfComponents(1);
  isLeaf->SetNumberOfTuples(nVertices);
  isLeaf->SetName("IsLeaf");

  // compute joinTree
  if(this->Type == 1) {
    std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
      persistencePairsJoin{};
    std::vector<ExTreeM::Branch> mergeTreeJoin{};

    int status = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
    for(size_t i = 0; i < nVertices; i++) {
      orderArrayData[i] = nVertices - orderArrayData[i] - 1;
    }

    ttkTypeMacroT(triangulation->getType(),
                  (status = this->computePairs<T0>(
                     persistencePairsJoin, mergeTreeJoin,
                     ttkUtils::GetPointer<ttk::SimplexId>(segmentationId),
                     ttkUtils::GetPointer<unsigned char>(isLeaf),
                     ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold),
                     ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
                     orderArrayData, (T0 *)triangulation->getData())));

    if(status != 1)
      return 0;

    auto outputPoints = vtkUnstructuredGrid::GetData(outputVector, 0);
    auto outputMergeTreeJoin = vtkUnstructuredGrid::GetData(outputVector, 1);
    ttkTypeMacroT(
      triangulation->getType(),
      getMergeTree<T0>(outputMergeTreeJoin, mergeTreeJoin, scalarArray,
                       (T0 *)triangulation->getData()));
    ttkTypeMacroT(
      triangulation->getType(),
      getMergeTreePoints<T0>(outputPoints, persistencePairsJoin, scalarArray,
                             (T0 *)triangulation->getData()));

  } else {
    std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
      persistencePairsSplit{};
    std::vector<ExTreeM::Branch> mergeTreeSplit{};

    int status = 0;

    ttkTypeMacroT(triangulation->getType(),
                  (status = this->computePairs<T0>(
                     persistencePairsSplit, mergeTreeSplit,
                     ttkUtils::GetPointer<ttk::SimplexId>(segmentationId),
                     ttkUtils::GetPointer<unsigned char>(isLeaf),
                     ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
                     ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold),
                     orderArrayData, (T0 *)triangulation->getData())));

    if(status != 1)
      return 0;

    auto outputPoints = vtkUnstructuredGrid::GetData(outputVector, 0);
    auto outputMergeTreeSplit = vtkUnstructuredGrid::GetData(outputVector, 1);

    ttkTypeMacroT(
      triangulation->getType(),
      getMergeTree<T0>(outputMergeTreeSplit, mergeTreeSplit, scalarArray,
                       (T0 *)triangulation->getData()));
    ttkTypeMacroT(
      triangulation->getType(),
      getMergeTreePoints<T0>(outputPoints, persistencePairsSplit, scalarArray,
                             (T0 *)triangulation->getData()));
  }

  // Create segmentation output
  {
    segmentationPD->AddArray(segmentationId);
    segmentationPD->AddArray(isLeaf);
  }

  return 1;
}
