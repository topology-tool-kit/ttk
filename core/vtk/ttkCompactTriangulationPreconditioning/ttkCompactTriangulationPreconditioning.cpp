#include <ttkCompactTriangulationPreconditioning.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkCompactTriangulationPreconditioning);

ttkCompactTriangulationPreconditioning::
  ttkCompactTriangulationPreconditioning() {
  this->Threshold = 1000;
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
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

  // Get input object from input vector
  // Note: has to be a vtkPointSet as required by FillInputPortInformation
  vtkPointSet *inputDataSet = vtkPointSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");

  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
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

  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkPointSet *outputDataSet = vtkPointSet::GetData(outputVector, 0);
  vector<SimplexId> vertexMap(this->vertices.size());
  vtkUnstructuredGrid *outputMesh
    = vtkUnstructuredGrid::SafeDownCast(outputDataSet);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkIntArray> indices = vtkSmartPointer<vtkIntArray>::New();

  indices->SetNumberOfComponents(1);
  indices->SetName(compactTriangulationIndex);

  // insert the vertices in the output mesh
  for(size_t i = 0; i < this->vertices.size(); i++) {
    float x, y, z;
    triangulation->getVertexPoint(this->vertices.at(i), x, y, z);
    points->InsertNextPoint(x, y, z);
    vertexMap[this->vertices.at(i)] = i;
    indices->InsertNextTuple1(this->nodes.at(i));
  }
  outputMesh->SetPoints(points);

  vtkPointData *pointData = outputMesh->GetPointData();
  pointData->AddArray(indices);

  std::vector<std::string>::iterator iter = scalarFields.begin();
  while(iter != scalarFields.end()) {
    vtkDataArray *inputArray
      = inputDataSet->GetPointData()->GetArray(iter->data());
    if(inputArray == nullptr) {
      iter++;
      continue;
    }

    iter = scalarFields.erase(iter);
    vtkDataArray *newField = nullptr;

    // // To make sure that the selected array can be processed by this filter,
    // // one should also check that the array association and format is
    // correct. if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    //   this->printErr("Input array needs to be a point data array.");
    //   return 0;
    // }

    switch(inputArray->GetDataType()) {

      case VTK_CHAR:
        newField = vtkCharArray::New();
        break;

      case VTK_DOUBLE:
        newField = vtkDoubleArray::New();
        break;

      case VTK_FLOAT:
        newField = vtkFloatArray::New();
        break;

      case VTK_INT:
        newField = vtkIntArray::New();
        break;

      case VTK_ID_TYPE:
        newField = vtkIdTypeArray::New();
        break;
    }

    newField->DeepCopy(inputArray);
    for(size_t j = 0; j < this->vertices.size(); j++) {
      newField->SetTuple(j, inputArray->GetTuple(this->vertices.at(j)));
    }

    pointData->AddArray(newField);
  }

  // insert the cells in the output mesh
  outputMesh->Allocate(this->cells.size());
  int dimension = triangulation->getCellVertexNumber(0);

  for(unsigned int i = 0; i < this->cells.size(); i++) {
    vtkIdType cell[dimension];
    for(int j = 0; j < dimension; j++) {
      SimplexId vertexId;
      triangulation->getCellVertex(this->cells.at(i), j, vertexId);
      cell[j] = vertexMap[vertexId];
    }
    sort(cell, cell + dimension);
    if(dimension == 2) {
      outputMesh->InsertNextCell(VTK_LINE, 2, cell);
    } else if(dimension == 3) {
      outputMesh->InsertNextCell(VTK_TRIANGLE, 3, cell);
    } else if(dimension == 4) {
      outputMesh->InsertNextCell(VTK_TETRA, 4, cell);
    } else {
      this->printErr("Should not get here!");
    }
  }

  // the cell is the same order as the input mesh
  vtkCellData *cellData = outputMesh->GetCellData();
  iter = scalarFields.begin();
  while(iter != scalarFields.end()) {
    vtkDataArray *inputArray
      = inputDataSet->GetCellData()->GetArray(iter->data());
    if(inputArray == nullptr) {
      iter++;
      continue;
    }

    iter = scalarFields.erase(iter);
    cellData->AddArray(inputArray);
  }

  if(!scalarFields.empty()) {
    this->printErr("Some scalar fields are not processed!");
    return 0;
  }

  // scalarFields.clear();

  // return success
  return 1;
}
