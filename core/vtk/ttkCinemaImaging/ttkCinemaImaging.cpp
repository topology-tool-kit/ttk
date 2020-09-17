#include <ttkCinemaImaging.h>
#include <ttkUtils.h>

#include <vtkInformation.h>
#include <vtkVersion.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSignedCharArray.h>
#include <vtkSmartPointer.h>

// Render Dependencies
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCompositeDataGeometryFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkWindowToImageFilter.h>

// Value Pass Dependencies
#if VTK_MAJOR_VERSION >= 7
#include <vtkCameraPass.h>
#include <vtkOpenGLRenderer.h>
#include <vtkRenderPassCollection.h>
#include <vtkSequencePass.h>
#include <vtkValuePass.h>
#endif

vtkStandardNewMacro(ttkCinemaImaging);

ttkCinemaImaging::ttkCinemaImaging() {
  this->setDebugMsgPrefix("CinemaImaging");
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaImaging::~ttkCinemaImaging() {
}

int ttkCinemaImaging::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  else if(port == 1)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  else
    return 0;
  return 1;
}

int ttkCinemaImaging::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  else
    return 0;
  return 1;
}

void setupRenderer(vtkRenderer *renderer,
                   vtkDataObject *object,
                   vtkCamera *camera) {
  auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputDataObject(object);

  auto actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  renderer->SetBackground(0.0, 0.0, 0.0);
  renderer->GradientBackgroundOff();
  renderer->AddActor(actor);
  renderer->SetActiveCamera(camera);
}

void setupWindow(vtkRenderWindow *window,
                 vtkRenderer *renderer,
                 double resolution[2]) {
  window->SetSize(resolution[0], resolution[1]);
  window->SetMultiSamples(0); // Disable AA
  window->OffScreenRenderingOn();
  window->AddRenderer(renderer);
};

void addValuePass(vtkDataSet *object,
                  int fieldType,
                  vtkRenderPassCollection *valuePassCollection,
                  std::vector<std::string> &valuePassNames) {
  vtkFieldData *fd = fieldType == 0 ? (vtkFieldData *)object->GetPointData()
                                    : (vtkFieldData *)object->GetCellData();

  for(size_t i = 0, j = fd->GetNumberOfArrays(); i < j; i++) {
    auto field = vtkDataArray::SafeDownCast(fd->GetAbstractArray(i));
    if(!field)
      continue;

    std::string name(field->GetName());

    size_t nComponents = field->GetNumberOfComponents();
    for(size_t c = 0; c < nComponents; c++) {
      auto valuePass = vtkSmartPointer<vtkValuePass>::New();
#if VTK_MAJOR_VERSION < 8 || (VTK_MAJOR_VERSION == 8 && VTK_MINOR_VERSION < 90)
      valuePass->SetRenderingMode(vtkValuePass::FLOATING_POINT);
#endif
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
};

void addGobalArray(std::vector<vtkSmartPointer<vtkAbstractArray>> &globalArrays,
                   std::string name,
                   size_t nValues,
                   double *values) {
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetName(name.data());
  array->SetNumberOfComponents(nValues);
  array->SetNumberOfTuples(1);
  for(size_t i = 0; i < nValues; i++)
    array->SetValue(i, values[i]);
  globalArrays.push_back(array);
};

int ttkCinemaImaging::RequestData(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  ttk::Timer timer;
  double t0 = 0;

  // ---------------------------------------------------------------------------
  // Get Input / Output
  // ---------------------------------------------------------------------------
  auto inputObject = vtkDataObject::GetData(inputVector[0]);
  auto inputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  inputAsMB->SetBlock(0, inputObject);

  auto toPolyFilter = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();
  toPolyFilter->SetInputData(inputAsMB);
  toPolyFilter->Update();
  auto inputAsPD = vtkPolyData::SafeDownCast(toPolyFilter->GetOutput());

  auto inputGrid = vtkPointSet::GetData(inputVector[1]);
  auto outputImages = vtkMultiBlockDataSet::GetData(outputVector);

  // ---------------------------------------------------------------------------
  // Parameters
  // ---------------------------------------------------------------------------
  double resolution[2];
  resolution[0] = this->Resolution[0];
  resolution[1] = this->Resolution[1];
  double camFocus[3];
  std::copy(this->CamFocus, this->CamFocus + 3, camFocus);
  double camNearFar[2];
  std::copy(this->CamNearFar, this->CamNearFar + 3, camNearFar);
  double camHeight = this->GetCamHeight();
  double camDir[3];
  double camUp[3] = {0, 0, 1};

  {
    if(this->GetCamFocusAuto() || this->GetCamNearFarAuto()
       || this->GetCamHeightAuto()) {
      double dxO, dyO, dzO, dO; // object
      double dxG, dyG, dzG, dG; // object

      auto getBoxDiameter = [](const double *bounds, double &dx, double &dy,
                               double &dz, double &d) {
        dx = bounds[1] - bounds[0];
        dy = bounds[3] - bounds[2];
        dz = bounds[5] - bounds[4];
        d = std::sqrt(dx * dx + dy * dy + dz * dz);
        return 1;
      };

      double gridBounds[6];
      inputGrid->GetBounds(gridBounds);
      getBoxDiameter(gridBounds, dxG, dyG, dzG, dG);

      double objectBounds[6];
      inputAsMB->GetBounds(objectBounds);
      getBoxDiameter(objectBounds, dxO, dyO, dzO, dO);

      if(this->GetCamFocusAuto()) {
        camFocus[0] = objectBounds[0] + dxO * 0.5;
        camFocus[1] = objectBounds[2] + dyO * 0.5;
        camFocus[2] = objectBounds[4] + dzO * 0.5;
      }

      if(this->GetCamNearFarAuto()) {
        double gDiameter = std::max(dxG, std::max(dyG, dzG));
        camNearFar[0] = gDiameter / 2 - dO / 2;
        camNearFar[1] = camNearFar[0] + dO;
      }

      if(this->GetCamHeightAuto()) {
        camHeight = dO;
      }
    }

    // Ensure that camNear is not 0
    camNearFar[0] = std::max(camNearFar[0], 0.01);

    // print parameters
    {
      std::string projectionMode
        = this->GetCamProjectionMode() ? "Perspective" : "Orthographic";
      std::string camTypeBasedParameterName
        = this->GetCamProjectionMode() ? "CamAngle" : "CamHeight";
      std::string camTypeBasedParameterValue
        = this->GetCamProjectionMode() ? std::to_string(this->GetCamAngle())
                                       : std::to_string(camHeight);
      this->printMsg({
        {"Resolution", std::to_string((int)resolution[0]) + " x "
                         + std::to_string((int)resolution[1])},
        {"Projection", projectionMode},
        {"CamNearFar",
         std::to_string(camNearFar[0]) + " / " + std::to_string(camNearFar[1])},
        {"CamFocus", "(" + std::to_string(camFocus[0]) + ", "
                       + std::to_string(camFocus[1]) + ", "
                       + std::to_string(camFocus[2]) + ")"},
      });
      this->printMsg(ttk::debug::Separator::L1);
    }
  }

  // ---------------------------------------------------------------------------
  // Camera
  // ---------------------------------------------------------------------------
  auto camera = vtkSmartPointer<vtkCamera>::New();
  camera->SetClippingRange(camNearFar);
  camera->SetFocalPoint(camFocus);
  if(this->GetCamProjectionMode() == 0) {
    camera->SetParallelProjection(true);
    camera->SetParallelScale(
      camHeight * 0.5); // *0.5 to convert CamHeight to weird VTK convention
  } else {
    camera->SetParallelProjection(false);
    camera->SetViewAngle(this->CamAngle);
  }

  // ---------------------------------------------------------------------------
  // Initialize Depth Renderer and Components
  // ---------------------------------------------------------------------------
  t0 = timer.getElapsedTime();
  this->printMsg(
    "Initializing Rendering Pipeline", 0, ttk::debug::LineMode::REPLACE);

  // Depth Pass Elements
  auto rendererDepth = vtkSmartPointer<vtkRenderer>::New();
  setupRenderer(rendererDepth, inputAsPD, camera);
  auto windowDepth = vtkSmartPointer<vtkRenderWindow>::New();
  setupWindow(windowDepth, rendererDepth, resolution);
  auto windowDepthToImageFilter
    = vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowDepthToImageFilter->SetInput(windowDepth);
  windowDepthToImageFilter->SetInputBufferTypeToZBuffer();

  // Value passes Elements
  size_t nValuePasses = 0;

#if VTK_MAJOR_VERSION >= 7
  auto rendererScalars = vtkSmartPointer<vtkRenderer>::New();
  setupRenderer(rendererScalars, inputAsPD, camera);
  auto windowScalars = vtkSmartPointer<vtkRenderWindow>::New();
  setupWindow(windowScalars, rendererScalars, resolution);

  auto valuePassCollection = vtkSmartPointer<vtkRenderPassCollection>::New();
  std::vector<std::string> valuePassNames;
  size_t firstValuePassIndex = 0;
  if(this->GetGenerateScalarImages()) {
    auto preventVTKBug = [](vtkDataSet *object) {
      auto pd = object->GetPointData();

      if(pd->GetNumberOfArrays() < 1) {
        size_t nP = object->GetNumberOfPoints();

        auto fakeArray = vtkSmartPointer<vtkSignedCharArray>::New();
        fakeArray->SetName("Fake");
        fakeArray->SetNumberOfComponents(1);
        fakeArray->SetNumberOfTuples(nP);
        auto fakeArrayData = (signed char *)ttkUtils::GetVoidPointer(fakeArray);
        for(size_t i = 0; i < nP; i++)
          fakeArrayData[i] = 0;
        pd->AddArray(fakeArray);
        return 1;
      }
      return 0;
    };

    if(preventVTKBug(inputAsPD)) {
      firstValuePassIndex = 1;
    };

    addValuePass(inputAsPD, 0, valuePassCollection, valuePassNames);
    addValuePass(inputAsPD, 1, valuePassCollection, valuePassNames);
    nValuePasses = valuePassNames.size();

    auto sequence = vtkSmartPointer<vtkSequencePass>::New();
    sequence->SetPasses(valuePassCollection);

    auto cameraPass = vtkSmartPointer<vtkCameraPass>::New();
    cameraPass->SetDelegatePass(sequence);

    auto glRenderer = vtkOpenGLRenderer::SafeDownCast(rendererScalars);
    glRenderer->SetPass(cameraPass);

    // First pass to setup everything
    windowScalars->Render();
  }
#else
  if(this->GetGenerateScalarImages()) {
    this->printErr("Rendering scalar images requires VTK  7.0 or higher");
    return 0;
  }
#endif

  this->printMsg(
    "Initializing Rendering Pipeline", 1, timer.getElapsedTime() - t0);

  // ---------------------------------------------------------------------------
  // Prepare Field Data for Depth Values
  // ---------------------------------------------------------------------------
  std::vector<vtkSmartPointer<vtkAbstractArray>> globalArrays;
  {
    addGobalArray(globalArrays, "CamHeight", 1, &camHeight);
    addGobalArray(globalArrays, "CamNearFar", 2, camNearFar);
    addGobalArray(globalArrays, "Resolution", 2, resolution);

    auto inputObjectFD = inputObject->GetFieldData();
    for(size_t i = 0, j = inputObjectFD->GetNumberOfArrays(); i < j; i++) {
      auto array = inputObjectFD->GetAbstractArray(i);
      auto copy = vtkSmartPointer<vtkAbstractArray>::Take(array->NewInstance());
      copy->DeepCopy(array);
      globalArrays.push_back(array);
    }
  }
  size_t nGlobalArrays = globalArrays.size();

  // ---------------------------------------------------------------------------
  // Iterate over Locations
  // ---------------------------------------------------------------------------
  auto inputGridPD = inputGrid->GetPointData();
  size_t nInputGridPD = inputGridPD->GetNumberOfArrays();

  double *camFocusData = nullptr;
  double *camDirData = nullptr;
  double *camUpData = nullptr;

  {
    auto checkGridArray
      = [](vtkFieldData *fd, std::string name, std::string &gridFieldNames,
           double *&data, ttkCinemaImaging *caller) {
          if(fd->HasArray(name.data())) {
            auto field
              = vtkDoubleArray::SafeDownCast(fd->GetAbstractArray(name.data()));
            if(!field) {
              caller->printErr(
                "Sampling grid field '" + name
                + "' is not a vtkDoubleArray with three components.");
              return 0;
            }
            data = (double *)ttkUtils::GetVoidPointer(field);
            gridFieldNames += " '" + name + "'";
          }
          return 1;
        };

    std::string gridFieldNames = "";
    if(!checkGridArray(
         inputGridPD, "CamFocus", gridFieldNames, camFocusData, this))
      return 0;
    if(!checkGridArray(
         inputGridPD, "CamDirection", gridFieldNames, camDirData, this))
      return 0;
    if(!checkGridArray(inputGridPD, "CamUp", gridFieldNames, camUpData, this))
      return 0;

    if(gridFieldNames.compare("") != 0)
      this->printMsg("Sampling grid has field(s): " + gridFieldNames);
  }

  // -------------------------------------------------------------------------
  // Render Images for all Camera Locations
  // -------------------------------------------------------------------------
  {
    size_t n = inputGrid->GetNumberOfPoints();
    t0 = timer.getElapsedTime();
    this->printMsg("Rendering " + std::to_string(n) + " images with "
                     + std::to_string(nValuePasses + 1) + " fields",
                   0, ttk::debug::LineMode::REPLACE);

    double camPosition[3] = {0, 0, 0};
    auto readCameraData = [](double target[3], double *src, int index) {
      target[0] = src[index];
      target[1] = src[index + 1];
      target[2] = src[index + 2];
    };

    auto addCamFieldData
      = [](vtkFieldData *fd, std::string name, double *data) {
          auto array = vtkSmartPointer<vtkDoubleArray>::New();
          array->SetName(name.data());
          array->SetNumberOfComponents(3);
          array->SetNumberOfTuples(1);
          array->SetValue(0, data[0]);
          array->SetValue(1, data[1]);
          array->SetValue(2, data[2]);
          fd->AddArray(array);
        };

    for(size_t i = 0; i < n; i++) {
      // Set Camera Position
      inputGrid->GetPoint(i, camPosition);

      // Cam Up Fix
      if(camPosition[0] == 0 && camPosition[2] == 0) {
        camPosition[0] = 0.00000000001;
        camPosition[2] = 0.00000000001;
      }
      camera->SetPosition(camPosition);

      if(camUpData != nullptr) {
        readCameraData(camUp, camUpData, i * 3);
        vtkMath::Normalize(camUp);
      }
      camera->SetViewUp(camUp);

      if(camFocusData != nullptr) {
        readCameraData(camFocus, camFocusData, i * 3);
        camera->SetFocalPoint(camFocus);
      }
      if(camDirData != nullptr) {
        readCameraData(camDir, camDirData, i * 3);
        camFocus[0] = camPosition[0] + camDir[0];
        camFocus[1] = camPosition[1] + camDir[1];
        camFocus[2] = camPosition[2] + camDir[2];
        camera->SetFocalPoint(camFocus);
      } else {
        camDir[0] = camFocus[0] - camPosition[0];
        camDir[1] = camFocus[1] - camPosition[1];
        camDir[2] = camFocus[2] - camPosition[2];
        vtkMath::Normalize(camDir);
      }

      // Initialize Output Image
      auto outputImage = vtkSmartPointer<vtkImageData>::New();
      auto outputImagePD = outputImage->GetPointData();

      // Initialize as depth image
      {
        windowDepthToImageFilter->Modified();
        windowDepthToImageFilter->Update();
        outputImage->DeepCopy(windowDepthToImageFilter->GetOutput());
        outputImagePD->GetAbstractArray(0)->SetName("Depth");
      }

// Render Scalar Images
#if VTK_MAJOR_VERSION >= 7
      if(nValuePasses > firstValuePassIndex) {
        windowScalars->Render();

        for(size_t j = firstValuePassIndex; j < nValuePasses; j++) {
          auto valuePass = vtkValuePass::SafeDownCast(
            valuePassCollection->GetItemAsObject(j));
          auto newValueArray = vtkSmartPointer<vtkFloatArray>::New();
          newValueArray->DeepCopy(
            valuePass->GetFloatImageDataArray(rendererScalars));
          newValueArray->SetName(valuePassNames[j].data());
          outputImagePD->AddArray(newValueArray);
        }
      }
#endif

      // Add Field Data
      {
        auto outputImageFD = outputImage->GetFieldData();

        // Global Arrays
        for(size_t j = 0; j < nGlobalArrays; j++) {
          outputImageFD->AddArray(globalArrays[j]);
        }

        // Specific Arrays
        addCamFieldData(outputImageFD, "CamPosition", camPosition);
        addCamFieldData(outputImageFD, "CamDirection", camDir);
        addCamFieldData(outputImageFD, "CamUp", camera->GetViewUp());

        for(size_t j = 0; j < nInputGridPD; j++) {
          auto array = inputGridPD->GetAbstractArray(j);
          auto newArray
            = vtkSmartPointer<vtkAbstractArray>::Take(array->NewInstance());
          std::string name(array->GetName());
          if(outputImageFD->HasArray(name.data()))
            name += "FromGrid";
          newArray->SetName(name.data());
          newArray->SetNumberOfComponents(array->GetNumberOfComponents());
          newArray->SetNumberOfTuples(1);
          newArray->SetTuple(0, i, array);

          outputImageFD->AddArray(newArray);
        }
      }

      // Add Image to MultiBlock
      outputImages->SetBlock(i, outputImage);
    }

    this->printMsg("Rendering " + std::to_string(n) + " images with "
                     + std::to_string(nValuePasses + 1) + " fields",
                   1, timer.getElapsedTime() - t0);
  }

  // print stats
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg(
    "Complete (#images: "
      + std::to_string(inputGrid->GetNumberOfPoints() * (nValuePasses + 1))
      + ")",
    1, timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
