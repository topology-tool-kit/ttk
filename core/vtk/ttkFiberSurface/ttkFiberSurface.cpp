#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <ttkFiberSurface.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkFiberSurface);

ttkFiberSurface::ttkFiberSurface() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int ttkFiberSurface::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkFiberSurface::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

template <typename VTK_T1, typename VTK_T2>
int ttkFiberSurface::dispatch(ttk::Triangulation *const triangulation) {

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
  if(RangeOctree) {
    ttkTemplateMacro(triangulation->getType(),
                     (this->buildOctree<VTK_T1, VTK_T2>(
                       static_cast<TTK_TT *>(triangulation->getData()))));
  }
#endif // TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE

  ttkTemplateMacro(triangulation->getType(),
                   (this->computeSurface<VTK_T1, VTK_T2>(
                     static_cast<TTK_TT *>(triangulation->getData()))));
  return 0;
}

int ttkFiberSurface::RequestData(vtkInformation *ttkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {

  using ttk::SimplexId;
  ttk::Timer t;

  const auto input = vtkDataSet::GetData(inputVector[0]);
  const auto polygon = vtkUnstructuredGrid::GetData(inputVector[1]);
  auto output = vtkPolyData::GetData(outputVector);

  const auto dataUfield = this->GetInputArrayToProcess(0, input);
  const auto dataVfield = this->GetInputArrayToProcess(1, input);
  const auto polygonUfield = this->GetInputArrayToProcess(2, polygon);
  const auto polygonVfield = this->GetInputArrayToProcess(3, polygon);

  if(dataUfield == nullptr || dataVfield == nullptr || polygonUfield == nullptr
     || polygonVfield == nullptr) {
    this->printErr("Could not find data array");
    return -1;
  }

  if(!(input->GetDataObjectType() == VTK_UNSTRUCTURED_GRID
       || input->GetDataObjectType() == VTK_IMAGE_DATA)) {
    this->printErr("Unsupported VTK data structure");
    return -5;
  }

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    return -1;
  }

  this->preconditionTriangulation(triangulation);
  outputVertexList_.clear();
  this->setGlobalVertexList(&outputVertexList_);
  this->setInputField(
    ttkUtils::GetVoidPointer(dataUfield), ttkUtils::GetVoidPointer(dataVfield));
  this->setPolygonEdgeNumber(polygon->GetNumberOfCells());
  threadedTriangleList_.resize(polygon->GetNumberOfCells());
  threadedVertexList_.resize(polygon->GetNumberOfCells());
  this->setPolygon(&inputPolygon_);

  this->setPointMerging(PointMerge);
  this->setPointMergingThreshold(PointMergeDistanceThreshold);

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
  if((!RangeOctree) || (dataUfield->GetMTime() > GetMTime())
     || (dataVfield->GetMTime() > GetMTime())) {

    this->printMsg("Resetting octree...");

    this->flushOctree();
    Modified();
  }
#endif

  inputPolygon_.clear();

#if !defined(_WIN32) || defined(_WIN32) && defined(VTK_USE_64BIT_IDS)
  const long long int *cellArray
    = polygon->GetCells()->GetData()->GetPointer(0);
#else
  int *pt = polygon->GetCells()->GetPointer();
  long long extra_pt = *pt;
  const long long int *cellArray = &extra_pt;
#endif
  SimplexId cellNumber = polygon->GetNumberOfCells();

  SimplexId vertexId0, vertexId1;
  std::pair<std::pair<double, double>, std::pair<double, double>> rangeEdge;

  for(SimplexId i = 0; i < cellNumber; i++) {

    vertexId0 = cellArray[3 * i + 1];
    vertexId1 = cellArray[3 * i + 2];

    rangeEdge.first.first = polygonUfield->GetTuple1(vertexId0);
    rangeEdge.first.second = polygonVfield->GetTuple1(vertexId0);

    rangeEdge.second.first = polygonUfield->GetTuple1(vertexId1);
    rangeEdge.second.second = polygonVfield->GetTuple1(vertexId1);

    inputPolygon_.push_back(rangeEdge);
  }

  for(size_t i = 0; i < threadedTriangleList_.size(); i++) {
    threadedTriangleList_[i].clear();
    this->setTriangleList(i, &(threadedTriangleList_[i]));
    threadedVertexList_[i].clear();
    this->setVertexList(i, &(threadedVertexList_[i]));
  }

#ifndef TTK_ENABLE_DOUBLE_TEMPLATING
  if(dataUfield->GetDataType() != dataVfield->GetDataType()) {
    this->printErr(
      "Scalar fields should have same input type. Use TTKPointDataConverter or "
      "TTKArrayEditor to convert array types.");
    return 0;
  }
  switch(dataUfield->GetDataType()) {
    vtkTemplateMacro((dispatch<VTK_TT, VTK_TT>(triangulation)));
  }
#else
  switch(vtkTemplate2PackMacro(
    dataUfield->GetDataType(), dataVfield->GetDataType())) {
    vtkTemplate2Macro((dispatch<VTK_T1, VTK_T2>(triangulation)));
  }
#endif // TTK_ENABLE_DOUBLE_TEMPLATING

  // prepare the VTK output
  // NOTE: right now, there is a copy of the output data. this is no good.
  // to fix.

  size_t triangleNumber = 0;

  for(size_t i = 0; i < threadedTriangleList_.size(); i++) {
    triangleNumber += threadedTriangleList_[i].size();
  }

  vtkNew<vtkPoints> outputVertexList{};
  vtkNew<vtkDoubleArray> outputU{};
  vtkNew<vtkDoubleArray> outputV{};
  vtkNew<vtkDoubleArray> outputParameterization{};
  vtkNew<vtkCellArray> outputTriangleList{};
  vtkNew<ttkSimplexIdTypeArray> outputEdgeIds{};
  vtkNew<ttkSimplexIdTypeArray> outputTetIds{};
  vtkNew<ttkSimplexIdTypeArray> outputCaseIds{};

  if(RangeCoordinates) {
    outputU->SetName(dataUfield->GetName());
    outputU->SetNumberOfTuples(outputVertexList_.size());

    outputV->SetName(dataVfield->GetName());
    outputV->SetNumberOfTuples(outputVertexList_.size());
  }

  if(EdgeParameterization) {
    outputParameterization->SetName("EdgeParameterization");
    outputParameterization->SetNumberOfTuples(outputVertexList_.size());
  }

  outputVertexList->SetNumberOfPoints(outputVertexList_.size());
  output->SetPoints(outputVertexList);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < outputVertexList_.size(); i++) {
    outputVertexList->SetPoint(i, outputVertexList_[i].p_[0],
                               outputVertexList_[i].p_[1],
                               outputVertexList_[i].p_[2]);
    if(RangeCoordinates) {
      outputU->SetTuple1(i, outputVertexList_[i].uv_.first);
      outputV->SetTuple1(i, outputVertexList_[i].uv_.second);
    }
    if(EdgeParameterization) {
      outputParameterization->SetTuple1(i, outputVertexList_[i].t_);
    }
  }
  if(RangeCoordinates) {
    output->GetPointData()->AddArray(outputU);
    output->GetPointData()->AddArray(outputV);
  } else {
    output->GetPointData()->RemoveArray(dataUfield->GetName());
    output->GetPointData()->RemoveArray(dataVfield->GetName());
  }
  if(EdgeParameterization) {
    output->GetPointData()->AddArray(outputParameterization);
  } else {
    output->GetPointData()->RemoveArray("EdgeParameterization");
  }

  if(EdgeIds) {
    outputEdgeIds->SetName("EdgeIds");
    outputEdgeIds->SetNumberOfTuples(triangleNumber);
  }

  if(TetIds) {
    outputTetIds->SetName("TetIds");
    outputTetIds->SetNumberOfTuples(triangleNumber);
  }

  if(CaseIds) {
    outputCaseIds->SetName("CaseIds");
    outputCaseIds->SetNumberOfTuples(triangleNumber);
  }

  vtkNew<vtkIdList> idList{};
  idList->SetNumberOfIds(3);

  triangleNumber = 0;
  for(size_t i = 0; i < threadedTriangleList_.size(); i++) {
    for(size_t j = 0; j < threadedTriangleList_[i].size(); j++) {
      for(int k = 0; k < 3; k++) {
        idList->SetId(k, threadedTriangleList_[i][j].vertexIds_[k]);
      }
      outputTriangleList->InsertNextCell(idList);
      if(EdgeIds) {
        outputEdgeIds->SetTuple1(triangleNumber, i);
      }
      if(TetIds) {
        outputTetIds->SetTuple1(
          triangleNumber, threadedTriangleList_[i][j].tetId_);
      }
      if(CaseIds) {
        outputCaseIds->SetTuple1(
          triangleNumber, threadedTriangleList_[i][j].caseId_);
      }
      triangleNumber++;
    }
  }
  output->SetPolys(outputTriangleList);
  if(EdgeIds) {
    output->GetCellData()->AddArray(outputEdgeIds);
  } else {
    output->GetCellData()->RemoveArray("EdgeIds");
  }
  if(TetIds) {
    output->GetCellData()->AddArray(outputTetIds);
  } else {
    output->GetCellData()->RemoveArray("TetIds");
  }
  if(CaseIds) {
    output->GetCellData()->AddArray(outputCaseIds);
  } else {
    output->GetCellData()->RemoveArray("CaseIds");
  }

  return 1;
}
