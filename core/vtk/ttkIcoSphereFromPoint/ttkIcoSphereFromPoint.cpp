#include <ttkIcoSphereFromPoint.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkIcoSphereFromPoint);

ttkIcoSphereFromPoint::ttkIcoSphereFromPoint() : ttkIcoSphere() {
  this->setDebugMsgPrefix("IcoSphereFromPoint");
  this->SetNumberOfInputPorts(1);
}
ttkIcoSphereFromPoint::~ttkIcoSphereFromPoint() {
}

int ttkIcoSphereFromPoint::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  else
    return 0;
  return 1;
}

template <typename VTK_TT>
int copyArrayData(vtkDataArray *oldArray,
                  vtkDataArray *newArray,
                  const size_t &nSpheres,
                  const size_t &nVerticesPerSphere,
                  const size_t &nComponents) {
  auto oldData = (VTK_TT *)oldArray->GetVoidPointer(0);
  auto newData = (VTK_TT *)newArray->GetVoidPointer(0);

  for(size_t i = 0; i < nSpheres; i++) {
    size_t sphereIndex = i * nVerticesPerSphere * nComponents;
    for(size_t j = 0; j < nComponents; j++) {
      const auto &value = oldData[i * nComponents + j];
      size_t newIndex = sphereIndex + j;
      for(size_t k = 0; k < nVerticesPerSphere; k++) {
        newData[newIndex] = value;
        newIndex += nComponents;
      }
    }
  }
  return 1;
}

int ttkIcoSphereFromPoint::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {
  auto input = vtkPointSet::GetData(inputVector[0], 0);
  size_t nPoints = input->GetNumberOfPoints();
  this->SetNumberOfIcoSpheres(nPoints);
  this->SetCenters((float *)(input->GetPoints()->GetVoidPointer(0)));

  // compute spheres
  int status
    = this->ttkIcoSphere::RequestData(request, inputVector, outputVector);
  if(!status)
    return 0;

  size_t nVertices = 0;
  size_t nTriangles = 0;
  this->computeNumberOfVerticesAndTriangles(
    this->GetNumberOfSubdivisions(), nVertices, nTriangles);

  auto output = vtkUnstructuredGrid::GetData(outputVector);
  auto outputPD = output->GetPointData();

  if(this->CopyPointData) {
    // copy point data
    auto inputPD = input->GetPointData();
    for(size_t i = 0, n = inputPD->GetNumberOfArrays(); i < n; i++) {
      auto oldArray = vtkDataArray::SafeDownCast(inputPD->GetAbstractArray(i));
      if(!oldArray)
        continue;
      auto newArray
        = vtkSmartPointer<vtkDataArray>::Take(oldArray->NewInstance());
      size_t nComponents = oldArray->GetNumberOfComponents();
      newArray->SetName(oldArray->GetName());
      newArray->SetNumberOfComponents(nComponents);
      newArray->SetNumberOfTuples(nVertices * nPoints);

      switch(newArray->GetDataType()) {
        vtkTemplateMacro((copyArrayData<VTK_TT>(
          oldArray, newArray, nPoints, nVertices, nComponents)));
      }
      outputPD->AddArray(newArray);
    }
  }

  return 1;
}