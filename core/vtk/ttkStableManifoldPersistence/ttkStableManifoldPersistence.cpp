#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <ttkMacros.h>
#include <ttkStableManifoldPersistence.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkStableManifoldPersistence);

ttkStableManifoldPersistence::ttkStableManifoldPersistence() {

  this->setDebugMsgPrefix("StableManifoldPersistence");
  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(1);
}

int ttkStableManifoldPersistence::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }else if(port == 2){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkStableManifoldPersistence::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkStableManifoldPersistence::AttachPersistence(vtkDataSet *output) const {

  ttk::Timer t;

  printMsg("Attaching persistence...", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  auto ascendingManifoldArray
    = output->GetPointData()->GetArray(ttk::MorseSmaleAscendingName);
  auto descendingManifoldArray
    = output->GetPointData()->GetArray(ttk::MorseSmaleDescendingName);

  auto sourceArray
    = output->GetCellData()->GetArray(ttk::MorseSmaleSourceIdName);
  auto destinationArray
    = output->GetCellData()->GetArray(ttk::MorseSmaleDestinationIdName);

  if((!ascendingManifoldArray) && (!descendingManifoldArray)
     && (!sourceArray)) {
    printErr("The input #0 is not a valid stable manifold.");
    return -3;
  }

  bool isSegmentation = false;

  if((ascendingManifoldArray) || (descendingManifoldArray))
    isSegmentation = true;

  vtkSmartPointer<vtkDoubleArray> persistenceArray
    = vtkSmartPointer<vtkDoubleArray>::New();
  persistenceArray->SetName(ttk::PersistenceName);

  if(!isSegmentation) {

    int cellNumber = output->GetNumberOfCells();

    persistenceArray->SetNumberOfTuples(cellNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < cellNumber; i++) {
      double cellId = -1;
      if((destinationArray) && (IsUnstable)) {
        destinationArray->GetTuple(i, &cellId);
      } else {
        sourceArray->GetTuple(i, &cellId);
      }
      double persistence = simplex2persistence[(int)cellId];
      persistenceArray->SetTuple(i, &persistence);
    }

    output->GetCellData()->AddArray(persistenceArray);
  } else {
    int vertexNumber = output->GetNumberOfPoints();

    persistenceArray->SetNumberOfTuples(vertexNumber);

    // TODO: parallel
    for(int i = 0; i < vertexNumber; i++) {
      int cellId = -1;
      double extremumId = -1;
      if(IsUnstable) {
        descendingManifoldArray->GetTuple(i, &extremumId);
        cellId = min2simplex[(int)extremumId];
      } else {
        ascendingManifoldArray->GetTuple(i, &extremumId);
        cellId = max2simplex[(int)extremumId];
      }
      double persistence = simplex2persistence[(int)cellId];
      persistenceArray->SetTuple(i, &persistence);
    }

    output->GetPointData()->AddArray(persistenceArray);
  }

  // TODO:
  // size=number of vertices for d-dimensional things
  //   persistenceArray->SetNumberOfTuples(stableManifold->GetNumberOfCells());

  //   vtkDataArray *ascendingArray =
  //     stableManifold->GetPointData();

  printMsg("Persistence attached!", 1, t.getElapsedTime(), threadNumber_);

  return 0;
}

int ttkStableManifoldPersistence::BuildSimplex2PersistenceMap(
  vtkDataSet *stableManifold,
  vtkPolyData *criticalPoints,
  vtkUnstructuredGrid *persistenceDiagram) {

  ttk::Timer t;

  printMsg("Building critical persistence map...", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  vtkDataArray *criticalPointVertexIdArray
    = criticalPoints->GetPointData()->GetArray(ttk::VertexScalarFieldName);
  vtkDataArray *criticalCellIdArray
    = criticalPoints->GetPointData()->GetArray(ttk::MorseSmaleCellIdName);
  vtkDataArray *criticalCellDimensionArray
    = criticalPoints->GetPointData()->GetArray(
      ttk::MorseSmaleCellDimensionName);
  if((!criticalPointVertexIdArray) || (!criticalCellIdArray)
     || (!criticalCellDimensionArray)) {
    printErr("The input #1 is not a valid set of critical points.");
    return -1;
  }

  vtkDataArray *vertexIdArray
    = persistenceDiagram->GetPointData()->GetArray(ttk::VertexScalarFieldName);
  vtkDataArray *persistenceArray
    = persistenceDiagram->GetCellData()->GetArray(ttk::PersistenceName);
  vtkDataArray *persistencePairIdArray
    = persistenceDiagram->GetCellData()->GetArray(
      ttk::PersistencePairIdentifierName);

  if((!vertexIdArray) || (!persistenceArray) || (!persistencePairIdArray)) {
    printErr("The input #2 is not a valid persistence diagram.");
    return -2;
  }

  int maximumVertexId = vertexIdArray->GetMaxNorm();
  std::vector<double> vertex2persistence(maximumVertexId + 1, -1);

  int cellNumber = persistenceDiagram->GetNumberOfCells();

  // NOTE: multi-saddle prevent a parallel loop here.
  for(int i = 0; i < cellNumber; i++) {

    double pairId = -1;
    persistencePairIdArray->GetTuple(i, &pairId);

    if(pairId != -1) {
      // not the diagonal

      vtkNew<vtkGenericCell> c;
      persistenceDiagram->GetCell(i, c);
      int pointId0 = c->GetPointId(0);
      int pointId1 = c->GetPointId(1);

      double persistence = -1;
      persistenceArray->GetTuple(i, &persistence);

      double vertexId0 = -1, vertexId1 = -1;
      vertexIdArray->GetTuple(pointId0, &vertexId0);
      vertexIdArray->GetTuple(pointId1, &vertexId1);

      if(vertex2persistence[(int)vertexId0] < persistence)
        vertex2persistence[(int)vertexId0] = persistence;
      if(vertex2persistence[(int)vertexId1] < persistence)
        vertex2persistence[(int)vertexId1] = persistence;
    }
  }

  int maximumSimplexId = criticalCellIdArray->GetMaxNorm();
  simplex2persistence.resize(maximumSimplexId + 1, -1);

  int criticalPointNumber = criticalCellIdArray->GetNumberOfTuples();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < criticalPointNumber; i++) {
    double cellId = -1;
    criticalCellIdArray->GetTuple(i, &cellId);

    double vertexId = -1;
    criticalPointVertexIdArray->GetTuple(i, &vertexId);

    simplex2persistence[(int)cellId] = vertex2persistence[(int)vertexId];
  }

  // taking care of the maximum to simplex map
  // get dimensionality from the stable manifold (the following maps are only
  // useful if the module is called on the segmentation output of the
  // morse-smale complex; then, we'll get the dimension right).
  if(stableManifold->GetNumberOfCells()) {
    int cellType = stableManifold->GetCellType(0);
    int dimension = -1;
    if((cellType == VTK_TETRA) || (cellType == VTK_VOXEL))
      dimension = 3;
    else if((cellType == VTK_TRIANGLE) || (cellType == VTK_PIXEL))
      dimension = 2;
    else
      dimension = 1;

    int minimumNumber = 0, maximumNumber = 0;
    min2simplex.resize(maximumSimplexId + 1, -1);
    max2simplex.resize(maximumSimplexId + 1, -1);
    for(int i = 0; i < criticalPointNumber; i++) {
      double cellId = -1;
      double criticalSimplexDimension = -1;
      criticalCellDimensionArray->GetTuple(i, &criticalSimplexDimension);
      criticalCellIdArray->GetTuple(i, &cellId);
      if(criticalSimplexDimension == 0) {
        // minimum
        min2simplex[(int)cellId] = minimumNumber;
        minimumNumber++;
      } else if(criticalSimplexDimension == dimension) {
        // maximum
        max2simplex[(int)cellId] = maximumNumber;
        maximumNumber++;
        // TODO: debug here
      }
    }
  }

  printMsg(
    "Critical persistence map built!", 1, t.getElapsedTime(), threadNumber_);

  return 0;
}

int ttkStableManifoldPersistence::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  ttk::Timer t;

  const auto stableManifold = vtkDataSet::GetData(inputVector[0]);
  const auto criticalPoints = vtkPolyData::GetData(inputVector[1]);
  const auto persistenceDiagram = vtkUnstructuredGrid::GetData(inputVector[2]);

  int ret = this->BuildSimplex2PersistenceMap(
    stableManifold, criticalPoints, persistenceDiagram);

  if(ret)
    return ret;

  auto output = vtkDataSet::GetData(outputVector);
  output->ShallowCopy(stableManifold);

  ret = this->AttachPersistence(output);

  if(ret)
    return ret;

  // TODO: parallelism
  // TODO: improve above messages (with 100% 0 to 100)
  // TODO: documentation

  printMsg("Stable manifold total time", 1, t.getElapsedTime(), threadNumber_);

  return 1;
}
