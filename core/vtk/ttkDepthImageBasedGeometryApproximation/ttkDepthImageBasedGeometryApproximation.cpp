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
  vtkInformation *ttkNotUsed(request),
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

    if(inputImage == nullptr) {
      continue;
    }

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

    const size_t resX = resolution->GetValue(0);
    const size_t resY = resolution->GetValue(1);

    // Initialize output point buffer that holds one point for each input depth
    // image pixel
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(resX * resY);

    // Prepare output cell buffer that holds two triangles for every quad
    // consisting of 4 vertices
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    const size_t nTriangles = 2 * (resX - 1) * (resY - 1);

    auto triangleDistortions = vtkSmartPointer<vtkDoubleArray>::New();
    triangleDistortions->SetName("TriangleDistortion");
    triangleDistortions->SetNumberOfComponents(1);
    triangleDistortions->SetNumberOfTuples(nTriangles);

    auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
    connectivityArray->SetNumberOfTuples(3 * nTriangles);

    auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
    offsetArray->SetNumberOfTuples(nTriangles + 1);

    // Approximate geometry
    int status = 0;
    switch(depthArray->GetDataType()) {
      vtkTemplateMacro(
        (status = this->execute<VTK_TT, int>(
           // Output
           static_cast<float *>(ttkUtils::GetVoidPointer(points)),
           static_cast<double *>(ttkUtils::GetVoidPointer(triangleDistortions)),
           static_cast<int *>(ttkUtils::GetVoidPointer(connectivityArray)),
           static_cast<int *>(ttkUtils::GetVoidPointer(offsetArray)),

           // Input
           static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(depthArray)),
           static_cast<double *>(ttkUtils::GetVoidPointer(camPosition)),
           static_cast<double *>(ttkUtils::GetVoidPointer(camDirection)),
           static_cast<double *>(ttkUtils::GetVoidPointer(camUp)),
           static_cast<double *>(ttkUtils::GetVoidPointer(camNearFar)),
           static_cast<double *>(ttkUtils::GetVoidPointer(camHeight)),
           static_cast<double *>(ttkUtils::GetVoidPointer(resolution)))));
    }
    if(!status)
      return 0;

    // Represent approximated geometry via VTK
    auto mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    {
      auto cellArray = vtkSmartPointer<vtkCellArray>::New();
      cellArray->SetData(offsetArray, connectivityArray);

      mesh->SetPoints(points);
      mesh->SetCells(VTK_TRIANGLE, cellArray);
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
