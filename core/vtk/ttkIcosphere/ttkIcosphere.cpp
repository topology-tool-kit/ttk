#include <ttkIcosphere.h>

#include <ttkUtils.h>

#include <vtkInformation.h>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkIcosphere);

ttkIcosphere::ttkIcosphere() {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}
ttkIcosphere::~ttkIcosphere() {
}

int ttkIcosphere::FillInputPortInformation(int port, vtkInformation *info) {
  return 0;
}

int ttkIcosphere::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  else
    return 0;
  return 1;
}

int ttkIcosphere::RequestData(vtkInformation *request,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {
  // get parameter
  int nSpheres = this->Centers ? this->Centers->GetNumberOfPoints() : 1;
  int nSubdivisions = this->GetNumberOfSubdivisions();
  double radius = this->GetRadius();

  bool useDoublePrecision
    = this->Centers ? this->Centers->GetDataType() == VTK_DOUBLE : false;

  // prepare the output buffers
  size_t nVertices = 0;
  size_t nTriangles = 0;
  this->computeNumberOfVerticesAndTriangles(
    nVertices, nTriangles, nSubdivisions);

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataType(useDoublePrecision ? VTK_DOUBLE : VTK_FLOAT);
  points->SetNumberOfPoints(nSpheres * nVertices);

  auto normals = vtkSmartPointer<vtkFloatArray>::New();
  if(this->ComputeNormals) {
    normals->SetName("Normals");
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(nVertices * nSpheres);
  }

  auto cells = vtkSmartPointer<vtkCellArray>::New();

  // execute base code
  // TODO: Improve here
#ifdef TTK_CELL_ARRAY_NEW
  std::vector<vtkIdType> cellArray;
  cellArray.resize(nSpheres * nTriangles * 4);
  if(useDoublePrecision) {
    if(!this->computeIcospheres<double, vtkIdType>(
         (double *)ttkUtils::GetVoidPointer(points), cellArray.data(),

         nSpheres, nSubdivisions, radius,
         this->Centers ? (double *)ttkUtils::GetVoidPointer(this->Centers)
                       : this->Center,
         this->ComputeNormals ? (float *)ttkUtils::GetVoidPointer(normals)
                              : nullptr)) {
      return 0;
    }
  } else {
    float centerFloat[3]{
      (float)this->Center[0], (float)this->Center[1], (float)this->Center[2]};
    if(!this->computeIcospheres<float, vtkIdType>(
         (float *)ttkUtils::GetVoidPointer(points), cellArray.data(),

         nSpheres, nSubdivisions, radius,
         this->Centers ? (float *)ttkUtils::GetVoidPointer(this->Centers)
                       : centerFloat,
         this->ComputeNormals ? (float *)ttkUtils::GetVoidPointer(normals)
                              : nullptr))
      return 0;
  }
  ttkUtils::FillCellArrayFromSingle(
    cellArray.data(), nSpheres * nTriangles, cells);
#else
  if(useDoublePrecision) {
    if(!this->computeIcospheres<double, vtkIdType>(
         (double *)ttkUtils::GetVoidPointer(points),
         cells->WritePointer(nSpheres * nTriangles, nSpheres * nTriangles * 4),

         nSpheres, nSubdivisions, radius,
         this->Centers ? (double *)ttkUtils::GetVoidPointer(this->Centers)
                       : this->Center,
         this->ComputeNormals ? (float *)ttkUtils::GetVoidPointer(normals)
                              : nullptr))
      return 0;
  } else {
    float centerFloat[3]{
      (float)this->Center[0], (float)this->Center[1], (float)this->Center[2]};
    if(!this->computeIcospheres<float, vtkIdType>(
         (float *)ttkUtils::GetVoidPointer(points),
         cells->WritePointer(nSpheres * nTriangles, nSpheres * nTriangles * 4),

         nSpheres, nSubdivisions, radius,
         this->Centers ? (float *)ttkUtils::GetVoidPointer(this->Centers)
                       : centerFloat,
         this->ComputeNormals ? (float *)ttkUtils::GetVoidPointer(normals)
                              : nullptr))
      return 0;
  }
#endif

  // module easily finalize output
  {
    auto output = vtkUnstructuredGrid::GetData(outputVector);
    output->SetPoints(points);

    output->SetCells(VTK_TRIANGLE, cells);
    if(this->ComputeNormals)
      output->GetPointData()->SetNormals(normals);
  }

  return 1;
}
