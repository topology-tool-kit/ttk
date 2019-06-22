#include <ttkFiberSurface.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkFiberSurface)

  ttkFiberSurface::ttkFiberSurface() {
  UseAllCores = true;
  RangeCoordinates = true;
  EdgeParameterization = true;
  EdgeIds = true;
  TetIds = true;
  CaseIds = true;
  PointMerge = false;
  RangeOctree = true;
  PointMergeDistanceThreshold = 0.000001;
  SetNumberOfInputPorts(2);
}

ttkFiberSurface::~ttkFiberSurface() {
}
int ttkFiberSurface::doIt(vector<vtkDataSet *> &inputs,
                          vector<vtkDataSet *> &outputs) {

  Memory m;
  Timer t;

  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *polygon = vtkUnstructuredGrid::SafeDownCast(inputs[1]);
  vtkPolyData *output = vtkPolyData::SafeDownCast(outputs[0]);

  vtkDataArray *dataUfield = NULL, *dataVfield = NULL, *polygonUfield = NULL,
               *polygonVfield = NULL;

  if(DataUcomponent.length()) {
    dataUfield = input->GetPointData()->GetArray(DataUcomponent.data());
  } else {
    // default
    dataUfield = input->GetPointData()->GetArray(0);
  }
  if(!dataUfield) {
    stringstream msg;
    msg << "[ttkFiberSurface] Error1: Could not find data array '"
        << DataUcomponent << "'!" << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -1;
  }

  if(DataVcomponent.length()) {
    dataVfield = input->GetPointData()->GetArray(DataVcomponent.data());
  } else {
    // default
    dataVfield = input->GetPointData()->GetArray(0);
  }
  if(!dataVfield) {
    stringstream msg;
    msg << "[ttkFiberSurface] Error2: Could not find data array '"
        << DataVcomponent << "'!" << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -2;
  }

  if(PolygonUcomponent.length()) {
    polygonUfield = polygon->GetPointData()->GetArray(PolygonUcomponent.data());
  } else {
    // default
    polygonUfield = polygon->GetPointData()->GetArray(0);
  }
  if(!polygonUfield) {
    stringstream msg;
    msg << "[ttkFiberSurface] Error3: Could not find data array '"
        << PolygonUcomponent << "'!" << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -3;
  }

  if(PolygonVcomponent.length()) {
    polygonVfield = polygon->GetPointData()->GetArray(PolygonVcomponent.data());
  } else {
    // default
    polygonVfield = polygon->GetPointData()->GetArray(0);
  }
  if(!polygonVfield) {
    stringstream msg;
    msg << "[ttkFiberSurface] Error4: Could not find data array '"
        << PolygonVcomponent << "'!" << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -4;
  }

  if(!((input->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
       //     ||(input->GetDataObjectType() == TTK_UNSTRUCTURED_GRID)
       || (input->GetDataObjectType() == VTK_IMAGE_DATA))) {
    //     ||(input->GetDataObjectType() == TTK_IMAGE_DATA))){
    stringstream msg;
    msg << "[ttkFiberSurface] Error5: Unsupported VTK data-structure ("
        << input->GetDataObjectType() << ")" << endl;
    dMsg(cerr, msg.str(), fatalMsg);
    return -5;
  }

  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

  if(!triangulation)
    return -1;
  triangulation->setWrapper(this);
  fiberSurface_.setupTriangulation(triangulation);
  fiberSurface_.setWrapper(this);
  outputVertexList_.clear();
  fiberSurface_.setGlobalVertexList(&outputVertexList_);
  fiberSurface_.setInputField(
    dataUfield->GetVoidPointer(0), dataVfield->GetVoidPointer(0));
  fiberSurface_.setPolygonEdgeNumber(polygon->GetNumberOfCells());
  threadedTriangleList_.resize(polygon->GetNumberOfCells());
  threadedVertexList_.resize(polygon->GetNumberOfCells());
  fiberSurface_.setPolygon(&inputPolygon_);

  fiberSurface_.setPointMerging(PointMerge);
  fiberSurface_.setPointMergingThreshold(PointMergeDistanceThreshold);

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
  if((!RangeOctree) || (dataUfield->GetMTime() > GetMTime())
     || (dataVfield->GetMTime() > GetMTime())) {

    {
      stringstream msg;
      msg << "[ttkFiberSurface] Resetting octree..." << endl;
      dMsg(cout, msg.str(), infoMsg);
    }

    fiberSurface_.flushOctree();
    Modified();
  }
#endif

  inputPolygon_.clear();

#if !defined(_WIN32) || defined(_WIN32) && defined(VTK_USE_64BIT_IDS)
  const long long int *cellArray = polygon->GetCells()->GetPointer();
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

  for(SimplexId i = 0; i < (SimplexId)threadedTriangleList_.size(); i++) {
    threadedTriangleList_[i].clear();
    fiberSurface_.setTriangleList(i, &(threadedTriangleList_[i]));
    threadedVertexList_[i].clear();
    fiberSurface_.setVertexList(i, &(threadedVertexList_[i]));
  }

#ifdef TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE
  switch(vtkTemplate2PackMacro(
    dataUfield->GetDataType(), dataVfield->GetDataType())) {
    ttkTemplate2Macro({
      if(RangeOctree)
        fiberSurface_.buildOctree<VTK_T1 TTK_COMMA VTK_T2>();
      fiberSurface_.computeSurface<VTK_T1 TTK_COMMA VTK_T2>();
    });
  }
#else
  switch(vtkTemplate2PackMacro(
    dataUfield->GetDataType(), dataVfield->GetDataType())) {
    ttkTemplate2Macro(fiberSurface_.computeSurface<VTK_T1 TTK_COMMA VTK_T2>());
  }
#endif

  // prepare the VTK output
  // NOTE: right now, there is a copy of the output data. this is no good.
  // to fix.

  SimplexId triangleNumber = 0;

  for(SimplexId i = 0; i < (SimplexId)threadedTriangleList_.size(); i++) {
    triangleNumber += threadedTriangleList_[i].size();
  }

  vtkSmartPointer<vtkPoints> outputVertexList
    = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkDoubleArray> outputU
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> outputV
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> outputParameterization
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkCellArray> outputTriangleList
    = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<ttkSimplexIdTypeArray> outputEdgeIds
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vtkSmartPointer<ttkSimplexIdTypeArray> outputTetIds
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vtkSmartPointer<ttkSimplexIdTypeArray> outputCaseIds
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();

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
  for(SimplexId i = 0; i < (SimplexId)outputVertexList_.size(); i++) {
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

  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
  idList->SetNumberOfIds(3);

  triangleNumber = 0;
  for(SimplexId i = 0; i < (SimplexId)threadedTriangleList_.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)threadedTriangleList_[i].size(); j++) {
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

  {
    stringstream msg;
    msg << "[ttkFiberSurface] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
