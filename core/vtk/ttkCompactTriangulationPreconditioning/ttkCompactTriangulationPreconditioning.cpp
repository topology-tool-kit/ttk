#include <ttkCompactTriangulationPreconditioning.h>

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <vtkCellData.h>
#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

using namespace std;
using namespace ttk;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkCompactTriangulationPreconditioning);

ttkCompactTriangulationPreconditioning::
  ttkCompactTriangulationPreconditioning() {
  this->Threshold = 1000;
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  // ensure that modifying the selection re-triggers the filter
  // (c.f. vtkPassSelectedArrays.cxx)
  this->ArraySelection = vtkSmartPointer<vtkDataArraySelection>::New();
  this->ArraySelection->AddObserver(
    vtkCommand::ModifiedEvent, this,
    &ttkCompactTriangulationPreconditioning::Modified);
}

ttkCompactTriangulationPreconditioning::
  ~ttkCompactTriangulationPreconditioning() {
}

int ttkCompactTriangulationPreconditioning::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkCompactTriangulationPreconditioning::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkCompactTriangulationPreconditioning::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  // clear state
  this->clear();

  // Get input object from input vector
  // Note: has to be a vtkPointSet as required by FillInputPortInformation
  vtkPointSet *inputPointSet = vtkPointSet::GetData(inputVector[0]);
  if(!inputPointSet)
    return 0;

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");

  Triangulation *triangulation = ttkAlgorithm::GetTriangulation(inputPointSet);
  if(!triangulation)
    return 0;

  // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
  // this->preconditionTriangulation(triangulation); // implemented in base
  // class

  // Templatize over the different input array data types and call the base code
  int status = 0; // this integer checks if the base code returns an error
  ttkTemplateMacro(
    triangulation->getType(),
    (status = this->execute((TTK_TT *)triangulation->getData(), Threshold)));

  // On error cancel filter execution
  if(status != 1)
    return 0;

  // Get input data array selection
  vector<vtkDataArray *> pointDataArrays{};
  vector<vtkDataArray *> cellDataArrays{};

  vtkPointData *inPointData = inputPointSet->GetPointData();
  for(int i = 0; i < inPointData->GetNumberOfArrays(); i++) {
    vtkDataArray *curArray = inPointData->GetArray(i);
    if(curArray != nullptr && curArray->GetName() != nullptr
       && ArraySelection->ArrayIsEnabled(curArray->GetName())) {
      pointDataArrays.emplace_back(curArray);
    }
  }

  vtkCellData *inCellData = inputPointSet->GetCellData();
  for(int i = 0; i < inCellData->GetNumberOfArrays(); i++) {
    vtkDataArray *curArray = inCellData->GetArray(i);
    if(curArray != nullptr && curArray->GetName() != nullptr
       && ArraySelection->ArrayIsEnabled(curArray->GetName())) {
      cellDataArrays.emplace_back(curArray);
    }
  }

  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkPointSet *outputDataSet = vtkPointSet::GetData(outputVector, 0);
  outputDataSet->Initialize();

  vector<SimplexId> vertexMap(this->vertices.size());
  vtkUnstructuredGrid *outputMesh
    = vtkUnstructuredGrid::SafeDownCast(outputDataSet);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkIntArray> indices = vtkSmartPointer<vtkIntArray>::New();

  indices->SetNumberOfComponents(1);
  indices->SetName(compactTriangulationIndex);

  // insert vertices in the output mesh
  for(size_t i = 0; i < this->vertices.size(); i++) {
    float x, y, z;
    triangulation->getVertexPoint(this->vertices.at(i), x, y, z);
    points->InsertNextPoint(x, y, z);
    vertexMap[this->vertices.at(i)] = i;
    indices->InsertNextTuple1(this->nodes.at(i));
  }
  outputMesh->SetPoints(points);

  // vtkPointData *pointData = outputMesh->GetPointData();
  outputDataSet->GetPointData()->AddArray(indices);

  // insert cells in the output mesh
  outputMesh->Allocate(this->cells.size());
  int dimension = triangulation->getCellVertexNumber(0);

  for(unsigned int i = 0; i < this->cells.size(); i++) {
    std::array<vtkIdType, 4> cell{};
    for(int j = 0; j < dimension; j++) {
      SimplexId vertexId;
      triangulation->getCellVertex(this->cells.at(i), j, vertexId);
      cell[j] = vertexMap[vertexId];
    }
    std::sort(cell.begin(), cell.begin() + dimension);
    if(dimension == 2) {
      outputMesh->InsertNextCell(VTK_LINE, 2, cell.data());
    } else if(dimension == 3) {
      outputMesh->InsertNextCell(VTK_TRIANGLE, 3, cell.data());
    } else if(dimension == 4) {
      outputMesh->InsertNextCell(VTK_TETRA, 4, cell.data());
    } else {
      this->printErr("Should not get here!");
    }
  }

  // Modify the selected point data arrays with new indices
  for(vtkDataArray *scalarArray : pointDataArrays) {
    vtkSmartPointer<vtkDataArray> updatedField(scalarArray->NewInstance());
    updatedField->SetName(scalarArray->GetName());
    for(size_t j = 0; j < this->vertices.size(); j++) {
      updatedField->InsertTuple(j, scalarArray->GetTuple(this->vertices.at(j)));
    }

    outputDataSet->GetPointData()->AddArray(updatedField);
  }

  // Modify the selected cell data arrays with new indices
  for(vtkDataArray *scalarArray : cellDataArrays) {
    vtkSmartPointer<vtkDataArray> updatedField(scalarArray->NewInstance());
    updatedField->SetName(scalarArray->GetName());
    for(size_t j = 0; j < this->cells.size(); j++) {
      updatedField->InsertTuple(j, scalarArray->GetTuple(this->cells.at(j)));
    }

    outputDataSet->GetCellData()->AddArray(updatedField);
  }

  // return success
  return 1;
}
