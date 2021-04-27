#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
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
  } else if(port == 2) {
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

  vtkSmartPointer<vtkIntArray> pairTypeArray
    = vtkSmartPointer<vtkIntArray>::New();
  pairTypeArray->SetName(ttk::PersistencePairTypeName);

  if(!isSegmentation) {

    int cellNumber = output->GetNumberOfCells();

    persistenceArray->SetNumberOfTuples(cellNumber);
    pairTypeArray->SetNumberOfTuples(cellNumber);

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
      double persistence = simplex2persistence_[(int)cellId];
      double pairType = simplex2pairType_[(int)cellId];
      persistenceArray->SetTuple(i, &persistence);
      pairTypeArray->SetTuple(i, &pairType);
    }

    output->GetCellData()->AddArray(persistenceArray);
    output->GetCellData()->AddArray(pairTypeArray);
  } else {

    if((IsUnstable && descendingManifoldArray == nullptr)
       || ascendingManifoldArray == nullptr) {
      this->printErr("Missing array");
      return -4;
    }

    int vertexNumber = output->GetNumberOfPoints();

    persistenceArray->SetNumberOfTuples(vertexNumber);
    pairTypeArray->SetNumberOfTuples(vertexNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < vertexNumber; i++) {
      int cellId = -1;
      double extremumId = -1;
      if(IsUnstable) {
        descendingManifoldArray->GetTuple(i, &extremumId);
        cellId = min2simplex_[(int)extremumId];
      } else {
        ascendingManifoldArray->GetTuple(i, &extremumId);
        cellId = max2simplex_[(int)extremumId];
      }
      double persistence = simplex2persistence_[cellId];
      double pairType = simplex2pairType_[cellId];
      persistenceArray->SetTuple(i, &persistence);
      pairTypeArray->SetTuple(i, &pairType);
    }

    output->GetPointData()->AddArray(persistenceArray);
    output->GetPointData()->AddArray(pairTypeArray);
  }

  printMsg("Persistence attached!", 1, t.getElapsedTime(), threadNumber_);

  return 0;
}

int ttkStableManifoldPersistence::BuildSimplex2PersistenceMap(
  vtkDataSet *stableManifold,
  vtkPolyData *criticalPoints,
  vtkUnstructuredGrid *persistenceDiagram) {

  ttk::Timer t;

  printMsg("Building critical persistence map.", 0, 0, threadNumber_,
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
  vtkDataArray *persistencePairTypeArray
    = persistenceDiagram->GetCellData()->GetArray(ttk::PersistencePairTypeName);

  if((!vertexIdArray) || (!persistenceArray) || (!persistencePairIdArray)
     || (!persistencePairTypeArray)) {
    printErr("The input #2 is not a valid persistence diagram.");
    return -2;
  }

  int maximumVertexId = vertexIdArray->GetMaxNorm();
  std::vector<double> vertex2persistence(maximumVertexId + 1, -1);
  std::vector<int> vertex2pairType(maximumVertexId + 1, -1);

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

      double pairType = -1;
      persistencePairTypeArray->GetTuple(i, &pairType);

      double vertexId0 = -1, vertexId1 = -1;
      vertexIdArray->GetTuple(pointId0, &vertexId0);
      vertexIdArray->GetTuple(pointId1, &vertexId1);

      if(vertex2persistence[(int)vertexId0] < persistence) {
        vertex2persistence[(int)vertexId0] = persistence;
        vertex2pairType[(int)vertexId0] = pairType;
      }
      if(vertex2persistence[(int)vertexId1] < persistence) {
        vertex2persistence[(int)vertexId1] = persistence;
        vertex2pairType[(int)vertexId1] = pairType;
      }
    }
  }

  int maximumSimplexId = criticalCellIdArray->GetMaxNorm();
  simplex2persistence_.resize(maximumSimplexId + 1, -1);
  simplex2pairType_.resize(maximumSimplexId + 1, -1);

  int criticalPointNumber = criticalCellIdArray->GetNumberOfTuples();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < criticalPointNumber; i++) {
    double cellId = -1;
    criticalCellIdArray->GetTuple(i, &cellId);

    double vertexId = -1;
    criticalPointVertexIdArray->GetTuple(i, &vertexId);

    simplex2persistence_[(int)cellId] = vertex2persistence[(int)vertexId];
    simplex2pairType_[(int)cellId] = vertex2pairType[(int)vertexId];
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

    min2simplex_.clear();
    max2simplex_.clear();
    for(int i = 0; i < criticalPointNumber; i++) {
      double cellId = -1;
      double criticalSimplexDimension = -1;
      criticalCellDimensionArray->GetTuple(i, &criticalSimplexDimension);
      criticalCellIdArray->GetTuple(i, &cellId);
      if(criticalSimplexDimension == 0) {
        // minimum
        min2simplex_.push_back(cellId);
      } else if(criticalSimplexDimension == dimension) {
        // maximum
        max2simplex_.push_back(cellId);
      }
    }
  }

  printMsg(
    "Critical persistence map built!", 1, t.getElapsedTime(), threadNumber_);

  return 0;
}

int ttkStableManifoldPersistence::RequestData(
  vtkInformation *ttkNotUsed(request),
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

  // TODO: make an absolute version
  //  - color all separatrices on the path to the paired critical point

  printMsg("Stable manifold total time", 1, t.getElapsedTime(), threadNumber_);

  return 1;
}
