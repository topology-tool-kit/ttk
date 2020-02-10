#include <ttkDepthImageBasedGeometryApproximation.h>
#include <ttkUtils.h>

#include <vtkObjectFactory.h> // for new macro

#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

#include <vtkInformationVector.h>

#include <vtkSetGet.h> // vtkTemplateMacro

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

    // Approximate geometry
    switch(depthArray->GetDataType()) {
      vtkTemplateMacro({
        this->execute(
          // Input
          (VTK_TT *)depthArray->GetVoidPointer(0),
          (double *)camPosition->GetVoidPointer(0),
          (double *)camDirection->GetVoidPointer(0),
          (double *)camUp->GetVoidPointer(0),
          (double *)camNearFar->GetVoidPointer(0),
          (double *)camHeight->GetVoidPointer(0),
          (double *)resolution->GetVoidPointer(0),

          // Output
          (float *)points->GetVoidPointer(0),
          cells->WritePointer(nTriangles, nTriangles * 4),
          (double *)triangleDistortions->GetVoidPointer(0));
      });
    }

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
      for(size_t p = 0; p < inputPD->GetNumberOfArrays(); p++)
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
