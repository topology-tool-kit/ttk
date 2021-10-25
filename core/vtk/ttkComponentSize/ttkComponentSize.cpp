#include <ttkComponentSize.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkConnectivityFilter.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkComponentSize);

ttkComponentSize::ttkComponentSize() {
  this->setDebugMsgPrefix("ComponentSize");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  vtkWarningMacro("`TTK ComponentSize' is now deprecated. Please use "
                  "`Connectivity' instead.");
}

ttkComponentSize::~ttkComponentSize() {
}

int ttkComponentSize::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkComponentSize::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkComponentSize::RequestData(vtkInformation *ttkNotUsed(request),
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  ttk::Timer t;
  size_t threadNumber = this->getThreadNumber();

  this->printMsg(
    "Computing connected components", 0, 0, ttk::debug::LineMode::REPLACE);

  auto connectivityFilter = vtkSmartPointer<vtkConnectivityFilter>::New();
  connectivityFilter->SetInputData(vtkDataSet::GetData(inputVector[0]));
  connectivityFilter->SetExtractionModeToAllRegions();
  connectivityFilter->ColorRegionsOn();
  connectivityFilter->Update();

  size_t nRegions = connectivityFilter->GetNumberOfExtractedRegions();
  if(nRegions < 1) {
    this->printErr("Unable to compute connected components.");
    return 0;
  }

  this->printMsg(
    "Computing connected components (" + std::to_string(nRegions) + ")", 1,
    t.getElapsedTime());

  t.reStart();
  this->printMsg("Computing component sizes", 0, 0, threadNumber,
                 ttk::debug::LineMode::REPLACE);

  auto output = vtkDataSet::GetData(outputVector);
  output->ShallowCopy(connectivityFilter->GetOutput());

  size_t nVertices = output->GetNumberOfPoints();
  size_t nCells = output->GetNumberOfCells();

  auto vertexIds = (vtkIdType *)ttkUtils::GetVoidPointer(
    output->GetPointData()->GetArray("RegionId"));
  auto cellIds = (vtkIdType *)ttkUtils::GetVoidPointer(
    output->GetCellData()->GetArray("RegionId"));

  if(!vertexIds || !cellIds) {
    this->printErr("Unable to retrieve vertex and cell Identifiers.");
    return 0;
  }

  this->printMsg("Computing component sizes", 0.1, t.getElapsedTime(),
                 threadNumber, ttk::debug::LineMode::REPLACE);

  // count vertices per region
  std::vector<double> regionIdToVertexCountMap(nRegions, 0);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif
  for(size_t i = 0; i < nVertices; i++) {
    auto regionId = vertexIds[i];

#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
    regionIdToVertexCountMap[regionId]++;
  }

  // count cells per region
  std::vector<double> regionIdToCellCountMap(nRegions, 0);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif
  for(size_t i = 0; i < nCells; i++) {
    auto regionId = cellIds[i];

#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
    regionIdToCellCountMap[regionId]++;
  }

  // generate vertex number fields
  {
    auto vertexNumbers = vtkSmartPointer<vtkDoubleArray>::New();
    vertexNumbers->SetNumberOfComponents(1);
    vertexNumbers->SetNumberOfTuples(nVertices);
    vertexNumbers->SetName("VertexNumber");
    auto vertexNumbersData = (double *)ttkUtils::GetVoidPointer(vertexNumbers);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif
    for(size_t i = 0; i < nVertices; i++) {
      vertexNumbersData[i] = regionIdToVertexCountMap[vertexIds[i]];
    }

    output->GetPointData()->AddArray(vertexNumbers);
  }

  // generate cell number fields
  {
    auto cellNumbers = vtkSmartPointer<vtkDoubleArray>::New();
    cellNumbers->SetNumberOfComponents(1);
    cellNumbers->SetNumberOfTuples(nCells);
    cellNumbers->SetName("CellNumber");
    auto cellNumbersData = (double *)ttkUtils::GetVoidPointer(cellNumbers);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif
    for(size_t i = 0; i < nCells; i++) {
      cellNumbersData[i] = regionIdToCellCountMap[cellIds[i]];
    }

    output->GetCellData()->AddArray(cellNumbers);
  }

  this->printMsg(
    "Computing component sizes", 1, t.getElapsedTime(), threadNumber);

  return 1;
}
