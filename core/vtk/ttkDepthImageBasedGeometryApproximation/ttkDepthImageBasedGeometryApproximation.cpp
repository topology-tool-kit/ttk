#include <ttkDepthImageBasedGeometryApproximation.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

#include <vtkInformation.h>
// #include <vtkInformationVector.h>

vtkStandardNewMacro(ttkDepthImageBasedGeometryApproximation);

ttkDepthImageBasedGeometryApproximation::
  ttkDepthImageBasedGeometryApproximation() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}
ttkDepthImageBasedGeometryApproximation::
  ~ttkDepthImageBasedGeometryApproximation() {
}

int ttkDepthImageBasedGeometryApproximation::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
  else
    return 0;
  return 1;
}

int ttkDepthImageBasedGeometryApproximation::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;
  return 1;
}

int ttkDepthImageBasedGeometryApproximation::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  ttk::Timer globalTimer;

  // Prepare input and output
  auto inputMBD = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto outputMBD = vtkMultiBlockDataSet::GetData(outputVector);

  // Process each depth image individually
  size_t nBlocks = inputMBD->GetNumberOfBlocks();
  for(size_t i = 0; i < nBlocks; i++) {
    // Get vtkImageData
    auto inputImage = vtkImageData::SafeDownCast(inputMBD->GetBlock(i));

    // Get depth array
    auto depthArray = this->GetInputArrayToProcess(0, inputImage);

    // Get camera paramters from field data
    auto inputFD = inputImage->GetFieldData();

    auto camHeight
      = vtkDoubleArray::SafeDownCast(inputFD->GetAbstractArray("CamHeight"));
    auto camPosition
      = vtkDoubleArray::SafeDownCast(inputFD->GetAbstractArray("CamPosition"));
    auto camDirection
      = vtkDoubleArray::SafeDownCast(inputFD->GetAbstractArray("CamDirection"));
    auto camUp
      = vtkDoubleArray::SafeDownCast(inputFD->GetAbstractArray("CamUp"));
    auto camNearFar
      = vtkDoubleArray::SafeDownCast(inputFD->GetAbstractArray("CamNearFar"));
    auto resolution
      = vtkDoubleArray::SafeDownCast(inputFD->GetAbstractArray("Resolution"));

    // Check if all parameters are present
    if(depthArray == nullptr || camHeight == nullptr || camPosition == nullptr
       || camDirection == nullptr || camUp == nullptr || camNearFar == nullptr
       || resolution == nullptr) {
      this->printErr("Input depth image does not have all of the required "
                     "field data arrays (see Cinema Spec D - Data Product "
                     "Specification)");
      return 0;
    }

    size_t resX = resolution->GetValue(0);
    size_t resY = resolution->GetValue(1);

    // Initialize output point buffer that holds one point for each input depth
    // image pixel
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(resX * resY);

    // Prepare output cell buffer that holds two triangles for every quad
    // consisting of 4 vertices
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    size_t nTriangles = 2 * (resX - 1) * (resY - 1);
    // Vectors to hold raw approximated geometry

    auto triangleDistortions = vtkSmartPointer<vtkDoubleArray>::New();
    triangleDistortions->SetName("TriangleDistortion");
    triangleDistortions->SetNumberOfComponents(1);
    triangleDistortions->SetNumberOfTuples(nTriangles);

    // retrieve the numerous input array of the method.
    // TODO: use setter getter instead, 10 parameters is too much.
    // Input
    double *camPositionPtr
      = static_cast<double *>(ttkUtils::GetVoidPointer(camPosition));
    double *camDirectionPtr
      = static_cast<double *>(ttkUtils::GetVoidPointer(camDirection));
    double *camUpPtr = static_cast<double *>(ttkUtils::GetVoidPointer(camUp));
    double *camNearFarPtr
      = static_cast<double *>(ttkUtils::GetVoidPointer(camNearFar));
    double *camHeightPtr
      = static_cast<double *>(ttkUtils::GetVoidPointer(camHeight));
    double *resolutionPtr
      = static_cast<double *>(ttkUtils::GetVoidPointer(resolution));
    // Output
    float *pointsPtr // may be double, error then
      = static_cast<float *>(ttkUtils::GetVoidPointer(points));
    double *triangleDistortionsPtr
      = static_cast<double *>(ttkUtils::GetVoidPointer(triangleDistortions));
    // Cells output
    std::vector<vtkIdType> cellsCo;
    cellsCo.resize(nTriangles * 3);
    std::vector<vtkIdType> cellsOff;
    cellsOff.resize(nTriangles + 1);

    // Approximate geometry
    switch(depthArray->GetDataType()) {
      vtkTemplateMacro({
        VTK_TT *depthArrayPtr
          = static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(depthArray));
        this->execute(
          // Input
          depthArrayPtr, camPositionPtr, camDirectionPtr, camUpPtr,
          camNearFarPtr, camHeightPtr, resolutionPtr,

          // Output
          pointsPtr, cellsCo.data(), cellsOff.data(), triangleDistortionsPtr);
      });
    }

    ttkUtils::FillCellArrayFromDual(
      cellsCo.data(), cellsOff.data(), nTriangles, cells);

    // Represent approximated geometry via VTK
    auto mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    {
      mesh->SetPoints(points);
      mesh->SetCells(VTK_TRIANGLE, cells);
      mesh->GetCellData()->AddArray(triangleDistortions);
    }

    // copy image point data to point data of approximation
    {
      auto meshPD = mesh->GetPointData();
      auto inputPD = inputImage->GetPointData();
      for(int p = 0; p < inputPD->GetNumberOfArrays(); p++)
        meshPD->AddArray(inputPD->GetAbstractArray(p));

      outputMBD->SetBlock(i, mesh);
    }
  }

  // Print status
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete (#images: " + std::to_string(nBlocks) + ")", 1,
                 globalTimer.getElapsedTime(), this->threadNumber_);
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
