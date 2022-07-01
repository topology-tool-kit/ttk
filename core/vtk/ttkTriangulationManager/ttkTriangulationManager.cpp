#include <ttkTriangulationManager.h>

#include <vtkDataObject.h> // For port information
#include <vtkObjectFactory.h> // for new macro

#include <vtkCellData.h>
#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>

#include <CompactTriangulationPreconditioning.h>
#include <Triangulation.h>
#include <ttkUtils.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkTriangulationManager);

ttkTriangulationManager::ttkTriangulationManager() {
  this->setDebugMsgPrefix("TriangulationManager");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  // ensure that modifying the selection re-triggers the filter
  // (c.f. vtkPassSelectedArrays.cxx)
  this->ArraySelection = vtkSmartPointer<vtkDataArraySelection>::New();
  this->ArraySelection->AddObserver(
    vtkCommand::ModifiedEvent, this, &ttkTriangulationManager::Modified);
}

int ttkTriangulationManager::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTriangulationManager::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

void ttkTriangulationManager::processImplicit(
  ttk::Triangulation &triangulation) const {

  const bool prevPeriodic = triangulation.hasPeriodicBoundaries();
  if(prevPeriodic != this->Periodicity) {
    triangulation.setPeriodicBoundaryConditions(Periodicity);
    printMsg("Switching regular grid periodicity from "
             + (prevPeriodic ? std::string("ON") : std::string("OFF")) + " to "
             + (Periodicity ? std::string("ON") : std::string("OFF")));
  }

  const bool prevPreconditions = triangulation.hasImplicitPreconditions();
  if((this->PreconditioningStrategy == STRATEGY::NO_PRECONDITIONS
      && !prevPreconditions)
     || (this->PreconditioningStrategy == STRATEGY::WITH_PRECONDITIONS
         && prevPreconditions)) {
    return;
  }

  triangulation.setImplicitPreconditions(this->PreconditioningStrategy);
  const auto newPreconditions = triangulation.hasImplicitPreconditions();
  if(prevPreconditions != newPreconditions) {
    printMsg("Switching regular grid preconditions from "
             + (prevPreconditions ? std::string("ON") : std::string("OFF"))
             + " to "
             + (newPreconditions ? std::string("ON") : std::string("OFF")));
  }
}

int ttkTriangulationManager::processExplicit(
  vtkUnstructuredGrid *const output,
  vtkPointSet *const input,
  ttk::Triangulation &triangulation) const {

  if(output == nullptr || input == nullptr) {
    this->printErr("Empty data-sets");
    return 0;
  }
  ttk::Timer tm{};

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Compact explicit triangulation...");
  int status = 0; // this integer checks if the base code returns an error
  ttk::CompactTriangulationPreconditioning worker{};
  ttkTemplateMacro(
    triangulation.getType(),
    (status = worker.execute(
       static_cast<TTK_TT *>(triangulation.getData()), this->Threshold)));

  // On error cancel filter execution
  if(status != 1)
    return 0;

  // Get input data array selection
  std::vector<vtkDataArray *> pointDataArrays{};
  std::vector<vtkDataArray *> cellDataArrays{};

  vtkPointData *inPointData = input->GetPointData();
  for(int i = 0; i < inPointData->GetNumberOfArrays(); i++) {
    vtkDataArray *curArray = inPointData->GetArray(i);
    if(curArray != nullptr && curArray->GetName() != nullptr
       && ArraySelection->ArrayIsEnabled(curArray->GetName())) {
      pointDataArrays.emplace_back(curArray);
    }
  }

  vtkCellData *inCellData = input->GetCellData();
  for(int i = 0; i < inCellData->GetNumberOfArrays(); i++) {
    vtkDataArray *curArray = inCellData->GetArray(i);
    if(curArray != nullptr && curArray->GetName() != nullptr
       && ArraySelection->ArrayIsEnabled(curArray->GetName())) {
      cellDataArrays.emplace_back(curArray);
    }
  }

  output->Initialize();

  const auto &vertices{worker.getVertices()};
  const auto &nodes{worker.getNodes()};
  const auto &cells{worker.getCells()};
  std::vector<ttk::SimplexId> vertexMap(vertices.size());

  vtkNew<vtkPoints> points{};
  vtkNew<vtkIntArray> indices{};

  indices->SetNumberOfComponents(1);
  indices->SetName(ttk::compactTriangulationIndex);

  // insert vertices in the output mesh
  for(size_t i = 0; i < vertices.size(); i++) {
    float x, y, z;
    triangulation.getVertexPoint(vertices[i], x, y, z);
    points->InsertNextPoint(x, y, z);
    vertexMap[vertices[i]] = i;
    indices->InsertNextTuple1(nodes[i]);
  }
  output->SetPoints(points);

  output->GetPointData()->AddArray(indices);

  // insert cells in the output mesh
  output->Allocate(cells.size());
  const size_t dimension = triangulation.getCellVertexNumber(0);

  for(unsigned int i = 0; i < cells.size(); i++) {
    std::array<vtkIdType, 4> cell{};
    for(size_t j = 0; j < dimension; j++) {
      ttk::SimplexId vertexId;
      triangulation.getCellVertex(cells[i], j, vertexId);
      cell[j] = vertexMap[vertexId];
    }
    std::sort(cell.begin(), cell.begin() + dimension);
    if(dimension == 2) {
      output->InsertNextCell(VTK_LINE, 2, cell.data());
    } else if(dimension == 3) {
      output->InsertNextCell(VTK_TRIANGLE, 3, cell.data());
    } else if(dimension == 4) {
      output->InsertNextCell(VTK_TETRA, 4, cell.data());
    } else {
      this->printErr("Should not get here!");
    }
  }

  // Modify the selected point data arrays with new indices
  for(vtkDataArray *scalarArray : pointDataArrays) {
    vtkSmartPointer<vtkDataArray> updatedField(scalarArray->NewInstance());
    updatedField->SetName(scalarArray->GetName());
    for(size_t j = 0; j < vertices.size(); j++) {
      updatedField->InsertTuple(j, scalarArray->GetTuple(vertices[j]));
    }

    output->GetPointData()->AddArray(updatedField);
  }

  // Modify the selected cell data arrays with new indices
  for(vtkDataArray *scalarArray : cellDataArrays) {
    vtkSmartPointer<vtkDataArray> updatedField(scalarArray->NewInstance());
    updatedField->SetName(scalarArray->GetName());
    for(size_t j = 0; j < cells.size(); j++) {
      updatedField->InsertTuple(j, scalarArray->GetTuple(cells[j]));
    }

    output->GetCellData()->AddArray(updatedField);
  }

  this->printMsg("Done!", 1.0, tm.getElapsedTime(), 1);

  // return success
  return 1;
}

int ttkTriangulationManager::RequestData(vtkInformation *ttkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  auto *inputDataSet = vtkDataSet::GetData(inputVector[0]);
  auto *triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);

  if(triangulation == nullptr) {
    this->printErr("Triangulation is NULL");
    return 0;
  }

  if(inputDataSet->IsA("vtkImageData")) {
    this->processImplicit(*triangulation);
    auto *outputDataSet = vtkDataSet::GetData(outputVector, 0);
    outputDataSet->ShallowCopy(inputDataSet);
    return 1;
  } else if(inputDataSet->IsA("vtkUnstructuredGrid")
            || inputDataSet->IsA("vtkPolyData")) {
    return this->processExplicit(vtkUnstructuredGrid::GetData(outputVector, 0),
                                 vtkPointSet::SafeDownCast(inputDataSet),
                                 *triangulation);
  }

  return 1;
}
