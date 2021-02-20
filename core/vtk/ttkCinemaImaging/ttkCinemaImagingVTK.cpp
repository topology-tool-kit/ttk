#include <ttkCinemaImagingVTK.h>

#include <ttkCinemaImaging.h>

#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCameraPass.h>
#include <vtkOpenGLRenderer.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderPassCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSequencePass.h>
#include <vtkValuePass.h>
#include <vtkWindowToImageFilter.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkSignedCharArray.h>

#include <ttkUtils.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

ttk::ttkCinemaImagingVTK::ttkCinemaImagingVTK() {
  this->setDebugMsgPrefix("CinemaImaging(VTK)");
};
ttk::ttkCinemaImagingVTK::~ttkCinemaImagingVTK() {
}

int ttk::ttkCinemaImagingVTK::setupRenderer(vtkRenderer *renderer,
                                            vtkPointSet *object,
                                            vtkCamera *camera) const {
  auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputDataObject(object);

  auto actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  renderer->SetBackground(0.0, 0.0, 0.0);
  renderer->GradientBackgroundOff();
  renderer->AddActor(actor);
  renderer->SetActiveCamera(camera);

  return 1;
}

int ttk::ttkCinemaImagingVTK::setupWindow(vtkRenderWindow *window,
                                          vtkRenderer *renderer,
                                          const double resolution[2]) const {
  window->SetSize(resolution[0], resolution[1]);
  window->SetMultiSamples(0); // Disable AA
  window->OffScreenRenderingOn();
  window->AddRenderer(renderer);

  return 1;
};

int ttk::ttkCinemaImagingVTK::addValuePass(
  vtkPointSet *object,
  int fieldType,
  vtkRenderPassCollection *valuePassCollection,
  std::vector<std::string> &valuePassNames) const {
  auto fd = fieldType == 0 ? static_cast<vtkFieldData *>(object->GetPointData())
                           : static_cast<vtkFieldData *>(object->GetCellData());

  for(size_t i = 0, j = fd->GetNumberOfArrays(); i < j; i++) {
    auto array = fd->GetArray(i);
    if(!array)
      continue;

    std::string name(array->GetName());

    double minmax[2];
    array->GetRange(minmax);

    size_t nComponents = array->GetNumberOfComponents();
    for(size_t c = 0; c < nComponents; c++) {
      auto valuePass = vtkSmartPointer<vtkValuePass>::New();
      valuePass->SetInputArrayToProcess(fieldType == 0
                                          ? VTK_SCALAR_MODE_USE_POINT_FIELD_DATA
                                          : VTK_SCALAR_MODE_USE_CELL_FIELD_DATA,
                                        name.data());
      valuePass->SetInputComponentToProcess(c);

      valuePassCollection->AddItem(valuePass);
      valuePassNames.push_back(
        nComponents == 1 ? name : name + "_" + std::to_string(c));
    }
  }

  return 1;
};

int ttk::ttkCinemaImagingVTK::RenderVTKObject(
  vtkMultiBlockDataSet *outputImages,

  vtkPointSet *inputObject,
  vtkPointSet *inputGrid) const {

  auto inputObjectAsPD = vtkSmartPointer<vtkPolyData>::New();
  if(inputObject->IsA("vtkPolyData")) {
    inputObjectAsPD->ShallowCopy(inputObject);
  } else {
    auto surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputDataObject(inputObject);
    surfaceFilter->Update();
    inputObjectAsPD->ShallowCopy(surfaceFilter->GetOutput());
  }

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
  auto camFocalPoint = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamFocalPoint")));
  auto camHeight = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamHeight")));
  auto camAngle = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamAngle")));
  auto camNearFar = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("CamNearFar")));
  auto resolution = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("Resolution")));
  auto projectionMode = static_cast<double *>(
    ttkUtils::GetVoidPointer(camParameters->GetArray("ProjectionMode")));

  // Camera
  auto camera = vtkSmartPointer<vtkCamera>::New();

  // Depth Pass Elements
  auto rendererDepth = vtkSmartPointer<vtkRenderer>::New();
  this->setupRenderer(rendererDepth, inputObjectAsPD, camera);

  auto windowDepth = vtkSmartPointer<vtkRenderWindow>::New();
  this->setupWindow(windowDepth, rendererDepth, resolution);

  auto windowDepthToImageFilter
    = vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowDepthToImageFilter->SetInput(windowDepth);
  windowDepthToImageFilter->SetInputBufferTypeToZBuffer();

  // Value passes Elements
  size_t nValuePasses = 0;

  auto rendererScalars = vtkSmartPointer<vtkRenderer>::New();
  this->setupRenderer(rendererScalars, inputObjectAsPD, camera);

  auto windowScalars = vtkSmartPointer<vtkRenderWindow>::New();
  this->setupWindow(windowScalars, rendererScalars, resolution);

  auto valuePassCollection = vtkSmartPointer<vtkRenderPassCollection>::New();
  std::vector<std::string> valuePassNames;
  size_t firstValuePassIndex = 0;

  auto preventVTKBug = [](vtkPointSet *object) {
    auto pd = object->GetPointData();
    auto cd = object->GetCellData();

    if(pd->GetNumberOfArrays() < 1 && cd->GetNumberOfArrays() > 0) {
      size_t nP = object->GetNumberOfPoints();

      auto fakeArray = vtkSmartPointer<vtkSignedCharArray>::New();
      fakeArray->SetName("Fake");
      fakeArray->SetNumberOfComponents(1);
      fakeArray->SetNumberOfTuples(nP);
      auto fakeArrayData = (signed char *)fakeArray->GetVoidPointer(0);
      for(size_t i = 0; i < nP; i++)
        fakeArrayData[i] = 0;
      pd->AddArray(fakeArray);
      return 1;
    }
    return 0;
  };

  if(preventVTKBug(inputObjectAsPD)) {
    firstValuePassIndex = 1;
  };

  this->addValuePass(inputObjectAsPD, 0, valuePassCollection, valuePassNames);
  this->addValuePass(inputObjectAsPD, 1, valuePassCollection, valuePassNames);
  nValuePasses = valuePassNames.size();

  auto sequence = vtkSmartPointer<vtkSequencePass>::New();
  sequence->SetPasses(valuePassCollection);

  auto cameraPass = vtkSmartPointer<vtkCameraPass>::New();
  cameraPass->SetDelegatePass(sequence);

  auto glRenderer = vtkOpenGLRenderer::SafeDownCast(rendererScalars);
  glRenderer->SetPass(cameraPass);

  // First pass to setup everything
  windowScalars->Render();

  for(int i = 0; i < nSamplingPositions; i++) {
    ttk::Timer timer;
    int resX = resolution[i * 2];
    int resY = resolution[i * 2 + 1];

    this->printMsg("Rendering Image ("
                     + std::string(projectionMode[i] ? "P" : "O") + "|"
                     + std::to_string(resX) + "x" + std::to_string(resY) + ")",
                   0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

    double camPos[3]{samplingPositions[i * 3], samplingPositions[i * 3 + 1],
                     samplingPositions[i * 3 + 2]};

    if(projectionMode[i] == 0) {
      camera->SetParallelProjection(true);
      camera->SetParallelScale(
        camHeight[i]
        * 0.5); // *0.5 to convert CamHeight to weird VTK convention
    } else {
      camera->SetParallelProjection(false);
      camera->SetViewAngle(camAngle[i] * resolution[i * 2 + 1]
                           / resolution[i * 2]);
    }

    camera->SetPosition(camPos);
    camera->SetViewUp(&camUp[i * 3]);
    camera->SetFocalPoint(&camFocalPoint[i * 3]);
    camera->SetClippingRange(&camNearFar[i * 2]);

    // Initialize Output
    auto outputImage = vtkSmartPointer<vtkImageData>::New();
    outputImage->SetDimensions(resolution[i * 2], resolution[i * 2 + 1], 1);
    outputImage->SetSpacing(1, 1, 1);
    outputImage->SetOrigin(0, 0, 0);
    outputImage->AllocateScalars(VTK_FLOAT, 1);

    auto outputImagePD = outputImage->GetPointData();

    // Initialize as depth image
    windowDepthToImageFilter->Modified();
    windowDepthToImageFilter->Update();
    outputImage->DeepCopy(windowDepthToImageFilter->GetOutput());
    outputImagePD->GetAbstractArray(0)->SetName("Depth");

    ttkCinemaImaging::AddAllFieldDataArrays(inputGrid, outputImage, i);

    // Render Scalar Images
    if(nValuePasses > firstValuePassIndex) {
      windowScalars->Render();

      for(size_t j = firstValuePassIndex; j < nValuePasses; j++) {
        auto valuePass
          = vtkValuePass::SafeDownCast(valuePassCollection->GetItemAsObject(j));
        auto newValueArray = vtkSmartPointer<vtkFloatArray>::New();
        newValueArray->DeepCopy(
          valuePass->GetFloatImageDataArray(rendererScalars));
        newValueArray->SetName(valuePassNames[j].data());
        outputImagePD->AddArray(newValueArray);
      }
    }

    // Add to output collection
    outputImages->SetBlock(i, outputImage);

    this->printMsg("Rendering Image ("
                     + std::string(projectionMode[i] ? "P" : "O") + "|"
                     + std::to_string(resX) + "x" + std::to_string(resY) + ")",
                   1, timer.getElapsedTime(), this->threadNumber_);
  }
  this->printMsg(ttk::debug::Separator::L2);

  return 1;
};
