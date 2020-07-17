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

using namespace std;
using namespace ttk;

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

template <typename VTK_T1, typename VTK_T2, typename triangulationType>
int ttkFiberSurface::dispatch(const triangulationType *const triangulation) {
#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
  if(RangeOctree) {
    this->buildOctree<VTK_T1, VTK_T2>(triangulation);
  }
#endif // TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
  this->computeSurface<VTK_T1, VTK_T2>(triangulation);
  return 0;
}

int ttkFiberSurface::RequestData(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {

  Timer t;

  const auto input = vtkDataSet::GetData(inputVector[0]);
  const auto polygon = vtkUnstructuredGrid::GetData(inputVector[1]);
  auto output = vtkPolyData::GetData(outputVector);

  vtkDataArray *dataUfield = NULL, *dataVfield = NULL, *polygonUfield = NULL,
               *polygonVfield = NULL;

  if(DataUcomponent.length()) {
    dataUfield = input->GetPointData()->GetArray(DataUcomponent.data());
  } else {
    // default
    dataUfield = input->GetPointData()->GetArray(0);
  }
  if(!dataUfield) {
    this->printErr("Could not find data array");
    return -1;
  }

  if(DataVcomponent.length()) {
    dataVfield = input->GetPointData()->GetArray(DataVcomponent.data());
  } else {
    // default
    dataVfield = input->GetPointData()->GetArray(0);
  }
  if(!dataVfield) {
    this->printErr("Could not find data array");
    return -2;
  }

  if(PolygonUcomponent.length()) {
    polygonUfield = polygon->GetPointData()->GetArray(PolygonUcomponent.data());
  } else {
    // default
    polygonUfield = polygon->GetPointData()->GetArray(0);
  }
  if(!polygonUfield) {
    this->printErr("Could not find data array");
    return -3;
  }

  if(PolygonVcomponent.length()) {
    polygonVfield = polygon->GetPointData()->GetArray(PolygonVcomponent.data());
  } else {
    // default
    polygonVfield = polygon->GetPointData()->GetArray(0);
  }
  if(!polygonVfield) {
    this->printErr("Could not find data array");
    return -4;
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
    dataUfield->GetVoidPointer(0), dataVfield->GetVoidPointer(0));
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
  pair<pair<double, double>, pair<double, double>> rangeEdge;

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

  switch(vtkTemplate2PackMacro(
    dataUfield->GetDataType(), dataVfield->GetDataType())) {
    vtkTemplate2Macro((dispatch<VTK_T1, VTK_T2>(triangulation->getData())));
  }

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
    outputU->SetName(DataUcomponent.data());
    outputU->SetNumberOfTuples(outputVertexList_.size());

    outputV->SetName(DataVcomponent.data());
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
    output->GetPointData()->RemoveArray(DataUcomponent.data());
    output->GetPointData()->RemoveArray(DataVcomponent.data());
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
