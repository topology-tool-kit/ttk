#include <ttkCinemaImagingNative.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkUnsignedIntArray.h>

#include <ttkCinemaImaging.h>
#include <ttkUtils.h>

#include <BoundingVolumeHierarchy.h>

ttk::ttkCinemaImagingNative::ttkCinemaImagingNative() {
  this->setDebugMsgPrefix("CinemaImaging(Native)");
};

ttk::ttkCinemaImagingNative::~ttkCinemaImagingNative() {
}

int ttk::ttkCinemaImagingNative::RenderVTKObject(
  vtkMultiBlockDataSet *outputImages,

  vtkPointSet *inputObject,
  vtkPointSet *inputGrid) const {
  int status = 0;

  auto inputObjectCells = ttkCinemaImaging::GetCells(inputObject);

  // ---------------------------------------------------------------------------
  // Prepare Field Data for Depth Values
  // ---------------------------------------------------------------------------

  // iterate over sampling locations
  this->printMsg(ttk::debug::Separator::L2);
  float *samplingPositions
    = static_cast<float *>(ttkUtils::GetVoidPointer(inputGrid->GetPoints()));
  int nSamplingPositions = inputGrid->GetNumberOfPoints();
  auto camParameters = inputGrid->GetPointData();
  auto camUp = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamUp")));
  auto camNearFar = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamNearFar")));
  auto camDir = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamDirection")));
  auto camHeight = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamHeight")));
  auto camAngle = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamAngle")));
  auto resolution = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("Resolution")));
  auto projectionMode = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("ProjectionMode")));

  auto inputObjectConnectivityList = static_cast<vtkIdType *>(
    ttkUtils::GetVoidPointer(inputObjectCells->GetConnectivityArray()));

  ttk::Timer test;
  BoundingVolumeHierarchy<vtkIdType> bvh(
    static_cast<float *>(ttkUtils::GetVoidPointer(inputObject->GetPoints())),
    inputObjectConnectivityList, inputObjectCells->GetNumberOfCells());

  this->printMsg("BVH", 1, test.getElapsedTime(), 1);

  for(int i = 0; i < nSamplingPositions; i++) {

    double camPos[3]{samplingPositions[i * 3], samplingPositions[i * 3 + 1],
                     samplingPositions[i * 3 + 2]};

    // Initialize Output
    auto outputImage = vtkSmartPointer<vtkImageData>::New();
    outputImage->SetDimensions(resolution[0], resolution[1], 1);
    outputImage->SetSpacing(1, 1, 1);
    outputImage->SetOrigin(0, 0, 0);
    outputImage->AllocateScalars(VTK_FLOAT, 1);

    size_t nPixels = resolution[i * 2] * resolution[i * 2 + 1];
    auto outputImagePD = outputImage->GetPointData();

    auto depthBuffer = outputImagePD->GetArray(0);
    depthBuffer->SetName("Depth");

    auto primitiveIdArray = vtkSmartPointer<vtkUnsignedIntArray>::New();
    primitiveIdArray->SetName("PrimitiveId");
    primitiveIdArray->SetNumberOfComponents(1);
    primitiveIdArray->SetNumberOfTuples(nPixels);
    auto primitiveIdArrayData
      = static_cast<unsigned int *>(ttkUtils::GetVoidPointer(primitiveIdArray));
    outputImagePD->AddArray(primitiveIdArray);

    auto barycentricCoordinates = vtkSmartPointer<vtkFloatArray>::New();
    barycentricCoordinates->SetName("BarycentricCoordinates");
    barycentricCoordinates->SetNumberOfComponents(2);
    barycentricCoordinates->SetNumberOfTuples(nPixels);
    auto barycentricCoordinatesData
      = static_cast<float *>(ttkUtils::GetVoidPointer(barycentricCoordinates));
    outputImagePD->AddArray(barycentricCoordinates);

    // Render Object
    status = this->renderImage(
      static_cast<float *>(ttkUtils::GetVoidPointer(depthBuffer)),
      static_cast<unsigned int *>(ttkUtils::GetVoidPointer(primitiveIdArray)),
      static_cast<float *>(ttkUtils::GetVoidPointer(barycentricCoordinates)),
      inputObject->GetNumberOfPoints(),
      static_cast<float *>(ttkUtils::GetVoidPointer(inputObject->GetPoints())),
      inputObjectCells->GetNumberOfCells(), inputObjectConnectivityList, bvh,
      &resolution[i * 2], camPos, &camDir[i * 3], &camUp[i * 3], camHeight[i],
      projectionMode[i] == 0, camAngle[i]);
    if(!status)
      return 0;

    ttkCinemaImaging::Normalize(depthBuffer, &camNearFar[i * 2]);

    status = ttkCinemaImaging::MapPointAndCellData(
      outputImage,

      inputObject, this, primitiveIdArrayData, barycentricCoordinatesData,
      inputObjectConnectivityList);
    if(!status)
      return 0;

    ttkCinemaImaging::AddAllFieldDataArrays(inputGrid, outputImage, i);

    outputImages->SetBlock(i, outputImage);
  }
  this->printMsg(ttk::debug::Separator::L2);

  return status;
};
