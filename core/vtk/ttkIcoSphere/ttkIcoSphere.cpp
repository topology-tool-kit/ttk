#include <ttkIcoSphere.h>
#include <ttkUtils.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkCellArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkIcoSphere);

ttkIcoSphere::ttkIcoSphere() {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}
ttkIcoSphere::~ttkIcoSphere() {
}

int ttkIcoSphere::FillInputPortInformation(int port, vtkInformation *info) {
  return 0;
}

int ttkIcoSphere::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  else
    return 0;
  return 1;
}

int ttkIcoSphere::RequestData(vtkInformation *request,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {
  // get parameter
  int nSpheres = this->GetNumberOfIcoSpheres();
  int nSubdivisions = this->GetNumberOfSubdivisions();
  float radius = this->GetRadius();
  float *centers = this->Centers ? this->Centers : this->Center;

  // prepare the output buffers
  size_t nVertices = 0;
  size_t nTriangles = 0;
  this->computeNumberOfVerticesAndTriangles(
    nSubdivisions, nVertices, nTriangles);

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(nSpheres * nVertices);

  auto cells = vtkSmartPointer<vtkCellArray>::New();

  // execute base code
  if(!this->computeIcoSpheres<vtkIdType>(
       nSpheres, nSubdivisions, radius, centers,
       (float *)points->GetVoidPointer(0),
       cells->WritePointer(nSpheres * nTriangles, nSpheres * nTriangles * 4)))
    // (vtkIdType*)ttkUtils::WritePointer(cells->GetData(),nSpheres *
    // nTriangles, nSpheres * nTriangles * 4)))
    return 0;

  // finalize output
  {
    auto output = vtkUnstructuredGrid::GetData(outputVector);
    output->SetPoints(points);
    output->SetCells(VTK_TRIANGLE, cells);
  }

  return 1;
}
