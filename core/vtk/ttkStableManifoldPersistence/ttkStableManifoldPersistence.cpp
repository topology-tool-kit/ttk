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

int ttkStableManifoldPersistence::AttachPersistence(
  const std::vector<double> &simplex2persistence, vtkDataSet *output) const {

  ttk::Timer t;

  printMsg("Attaching persistence...", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  auto ascendingManifold
    = output->GetPointData()->GetArray(ttk::MorseSmaleAscendingName);
  auto descendingManifold
    = output->GetPointData()->GetArray(ttk::MorseSmaleDescendingName);

  auto sourceArray
    = output->GetCellData()->GetArray(ttk::MorseSmaleSourceIdName);
  auto destinationArray
    = output->GetCellData()->GetArray(ttk::MorseSmaleDestinationIdName);

  if((!ascendingManifold) && (!descendingManifold) && (!sourceArray)) {
    printErr("The input #0 is not a valid stable manifold.");
    return -3;
  }

  bool isSegmentation = false;

  if((ascendingManifold) || (descendingManifold))
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
  vtkPolyData *criticalPoints,
  vtkUnstructuredGrid *persistenceDiagram,
  std::vector<double> &simplex2persistence) const {

  ttk::Timer t;

  printMsg("Building critical persistence map...", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  vtkDataArray *criticalPointVertexIdArray
    = criticalPoints->GetPointData()->GetArray(ttk::VertexScalarFieldName);
  vtkDataArray *criticalCellIdArray
    = criticalPoints->GetPointData()->GetArray(ttk::MorseSmaleCellIdName);
  if((!criticalPointVertexIdArray) || (!criticalCellIdArray)) {
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

  std::vector<double> simplex2persistence;
  int ret = this->BuildSimplex2PersistenceMap(
    criticalPoints, persistenceDiagram, simplex2persistence);

  if(ret)
    return ret;

  auto output = vtkDataSet::GetData(outputVector);
  output->ShallowCopy(stableManifold);

  ret = this->AttachPersistence(simplex2persistence, output);

  if(ret)
    return ret;

  // TODO: parallelism
  // TODO: improve above messages (with 100% 0 to 100)
  // TODO: documentation

  printMsg("Stable manifold total time", 1, t.getElapsedTime(), threadNumber_);

  return 1;
}
