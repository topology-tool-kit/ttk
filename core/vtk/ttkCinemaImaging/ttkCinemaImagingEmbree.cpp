#include <ttkCinemaImagingEmbree.h>

#include <ttkCinemaImaging.h>
#include <ttkUtils.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>

#include <vtkFloatArray.h>
#include <vtkUnsignedIntArray.h>

ttk::ttkCinemaImagingEmbree::ttkCinemaImagingEmbree() {
  this->setDebugMsgPrefix("CinemaImaging(Embree)");
};
ttk::ttkCinemaImagingEmbree::~ttkCinemaImagingEmbree() {
}

int ttk::ttkCinemaImagingEmbree::RenderVTKObject(
  vtkMultiBlockDataSet *outputImages,

  vtkPointSet *inputObject,
  vtkPointSet *inputGrid) const {
  int status = 0;

#if TTK_ENABLE_EMBREE

  auto inputObjectCells = ttkCinemaImaging::GetCells(inputObject);

  RTCDevice device;
  status = this->initializeDevice(device);
  if(!status)
    return 0;

  auto inputObjectConnectivityList = static_cast<vtkIdType *>(
    ttkUtils::GetVoidPointer(inputObjectCells->GetConnectivityArray()));

  RTCScene scene;
  status = this->initializeScene<vtkIdType>(
    scene,

    device, inputObject->GetNumberOfPoints(),
    static_cast<float *>(ttkUtils::GetVoidPointer(inputObject->GetPoints())),
    inputObjectCells->GetNumberOfCells(), inputObjectConnectivityList);
  if(!status)
    return 0;

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
  auto camDir = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamDirection")));
  auto camHeight = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamHeight")));
  auto camNearFar = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamNearFar")));
  auto camAngle = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamAngle")));
  auto resolution = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("Resolution")));
  auto projectionMode = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("ProjectionMode")));

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
      primitiveIdArrayData, barycentricCoordinatesData,

      scene, &resolution[i * 2], camPos, &camDir[i * 3], &camUp[i * 3],
      projectionMode[i] == 0 ? camHeight[i] : camAngle[i],
      projectionMode[i] == 0);
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

  this->deallocateScene(device, scene);

#else
  TTK_FORCE_USE(outputImages);
  TTK_FORCE_USE(inputObject);
  TTK_FORCE_USE(inputGrid);

  this->printErr("TTK was build without EMBREE backend.");
#endif // TTK_ENABLE_EMBREE

  return status;
};
