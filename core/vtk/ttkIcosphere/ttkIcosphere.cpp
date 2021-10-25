#include <ttkIcosphere.h>

#include <ttkUtils.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

vtkStandardNewMacro(ttkIcosphere);

ttkIcosphere::ttkIcosphere() {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}
ttkIcosphere::~ttkIcosphere() {
}

int ttkIcosphere::FillInputPortInformation(int ttkNotUsed(port),
                                           vtkInformation *ttkNotUsed(info)) {
  return 0;
}

int ttkIcosphere::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  else
    return 0;
  return 1;
}

int ttkIcosphere::RequestData(vtkInformation *ttkNotUsed(request),
                              vtkInformationVector **ttkNotUsed(inputVector),
                              vtkInformationVector *outputVector) {
  // get parameter
  size_t nSpheres = this->Centers ? this->Centers->GetNumberOfTuples() : 1;
  size_t nSubdivisions = this->GetNumberOfSubdivisions();
  double radius = this->GetRadius();

  bool useDoublePrecision
    = this->Centers ? this->Centers->GetDataType() == VTK_DOUBLE : false;

  // prepare the output buffers
  size_t nVertices = 0;
  size_t nTriangles = 0;
  this->computeNumberOfVerticesAndTriangles(
    nVertices, nTriangles, nSubdivisions);

  const size_t nTotalVertices = nSpheres * nVertices;
  const size_t nTotalTriangles = nSpheres * nTriangles;

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataType(useDoublePrecision ? VTK_DOUBLE : VTK_FLOAT);
  points->SetNumberOfPoints(nTotalVertices);

  vtkSmartPointer<vtkDataArray> normals;
  if(this->ComputeNormals) {
    if(useDoublePrecision)
      normals = vtkSmartPointer<vtkDoubleArray>::New();
    else
      normals = vtkSmartPointer<vtkFloatArray>::New();

    normals->SetName("Normals");
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(nTotalVertices);
  }

  auto offsets = vtkSmartPointer<vtkIdTypeArray>::New();
  offsets->SetNumberOfTuples(nTotalTriangles + 1);
  auto offsetsData
    = static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(offsets));
  for(size_t i = 0; i <= nTotalTriangles; i++)
    offsetsData[i] = i * 3;

  auto connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
  connectivity->SetNumberOfTuples(nTotalTriangles * 3);

  int status = 0;
  if(useDoublePrecision) {
    typedef double DT;
    status = this->computeIcospheres<DT, vtkIdType>(
      ttkUtils::GetPointer<DT>(points->GetData()),
      ttkUtils::GetPointer<vtkIdType>(connectivity),

      nSpheres, nSubdivisions, radius,
      this->Centers ? ttkUtils::GetPointer<DT>(this->Centers) : this->Center,
      this->ComputeNormals ? ttkUtils::GetPointer<DT>(normals) : nullptr);
  } else {
    typedef float DT;
    DT centerFloat[3]{
      (DT)this->Center[0], (DT)this->Center[1], (DT)this->Center[2]};
    status = this->computeIcospheres<DT, vtkIdType>(
      ttkUtils::GetPointer<DT>(points->GetData()),
      ttkUtils::GetPointer<vtkIdType>(connectivity),

      nSpheres, nSubdivisions, radius,
      this->Centers ? ttkUtils::GetPointer<DT>(this->Centers) : centerFloat,
      this->ComputeNormals ? ttkUtils::GetPointer<DT>(normals) : nullptr);
  }

  if(!status)
    return 0;

  // finalize output
  {
    auto output = vtkPolyData::GetData(outputVector);
    output->SetPoints(points);

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetData(offsets, connectivity);
    output->SetPolys(cells);

    if(this->ComputeNormals)
      output->GetPointData()->SetNormals(normals);
  }

  return 1;
}
