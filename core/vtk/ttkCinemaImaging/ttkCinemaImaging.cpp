#include <ttkCinemaImaging.h>

#include <vtkVersion.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>

// Render Dependencies
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCompositeDataGeometryFilter.h>
#include <vtkCompositePolyDataMapper2.h>
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

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaImaging)

  int ttkCinemaImaging::RequestData(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  // Print Status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl;
    msg << "[ttkCinemaImaging] RequestData" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  Memory mem;
  Timer t;
  double t0 = 0;

  // Get Input / Output
  vtkInformation *inputObjectInfo = inputVector[0]->GetInformationObject(0);
  auto inputObject = inputObjectInfo->Get(vtkDataObject::DATA_OBJECT());

  vtkInformation *inGridInfo = inputVector[1]->GetInformationObject(0);
  auto inputGrid
    = vtkPointSet::SafeDownCast(inGridInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto outputImages = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // -------------------------------------------------------------------------
  // Initialize Shared Render Objects
  // -------------------------------------------------------------------------

  // Insert InputDataObject into MultiBlockDataSet
  auto inputMultiBlock = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  inputMultiBlock->SetBlock(0, inputObject);

  // Convert MultiBlock to PolyData
  auto toPoly = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();
  toPoly->SetInputData(inputMultiBlock);
  toPoly->Update();
  auto poly = toPoly->GetOutput();

  // Camera
  auto camera = vtkSmartPointer<vtkCamera>::New();
  camera->SetParallelProjection(true);
  camera->SetClippingRange(this->CamNearFar);
  camera->SetFocalPoint(this->CamFocus);
  camera->SetParallelScale(
    this->CamHeight * 0.5); // *0.5 to convert CamHeight to weird VTK convention

  // -------------------------------------------------------------------------
  // Initialize Depth Renderer and Components
  // -------------------------------------------------------------------------

  // Mapper
  auto mapper0 = vtkSmartPointer<vtkCompositePolyDataMapper2>::New();
  mapper0->SetInputConnection(toPoly->GetOutputPort());

  // Actor
  auto actor0 = vtkSmartPointer<vtkActor>::New();
  actor0->SetMapper(mapper0);

  // Renderer
  auto renderer0 = vtkSmartPointer<vtkRenderer>::New();
  renderer0->SetBackground(0, 0, 0);
  renderer0->AddActor(actor0);
  renderer0->SetActiveCamera(camera);

  // Window
  auto renderWindow0 = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow0->SetSize(this->Resolution);
  renderWindow0->SetMultiSamples(0); // Disable AA
  renderWindow0->OffScreenRenderingOn();
  renderWindow0->AddRenderer(renderer0);

// -------------------------------------------------------------------------
// Initialize Value Renderer and Components
// -------------------------------------------------------------------------
#if VTK_MAJOR_VERSION >= 7

  // Mapper
  auto mapper1 = vtkSmartPointer<vtkCompositePolyDataMapper2>::New();
  mapper1->SetInputConnection(toPoly->GetOutputPort());

  // Actor
  auto actor1 = vtkSmartPointer<vtkActor>::New();
  actor1->SetMapper(mapper1);

  // Renderer
  auto renderer1 = vtkSmartPointer<vtkRenderer>::New();
  renderer1->SetBackground(0, 0, 0);
  renderer1->AddActor(actor1);
  renderer1->SetActiveCamera(camera);

  // Window
  auto renderWindow1 = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow1->SetSize(this->Resolution);
  renderWindow1->SetMultiSamples(0); // Disable AA
  renderWindow1->OffScreenRenderingOn();
  renderWindow1->AddRenderer(renderer1);

  // Value Passes
  vector<pair<vtkValuePass *, string>> valuePassList;
  {
    auto valuePassCollection = vtkSmartPointer<vtkRenderPassCollection>::New();

    // Lambda function that generates vtkValuePasses for Point or Cell Data
    auto addValuePasses = [&](vtkFieldData *data,
                              int pointDataFlag // 0: Point Data, 1: Cell Data
                          ) {
      double minmax[2];
      size_t n = data->GetNumberOfArrays();

      for(size_t i = 0; i < n; i++) {
        auto values = data->GetArray(i);
        values->GetRange(minmax);

        size_t m = values->GetNumberOfComponents();
        for(size_t j = 0; j < m; j++) {
          auto valuePass = vtkSmartPointer<vtkValuePass>::New();
          valuePass->SetInputArrayToProcess(
            pointDataFlag == 0 ? VTK_SCALAR_MODE_USE_POINT_FIELD_DATA
                               : VTK_SCALAR_MODE_USE_CELL_FIELD_DATA,
            i);
          valuePass->SetRenderingMode(2);
          valuePass->SetInputComponentToProcess(j);
          valuePass->SetScalarRange(minmax[0], minmax[1]);

          valuePassCollection->AddItem(valuePass);
          valuePassList.push_back(make_pair(
            valuePass, m < 2 ? values->GetName()
                             : string(values->GetName()) + "_" + to_string(j)));
        }
      }
    };

    // Add Point Data Passes
    addValuePasses(poly->GetPointData(), 0);

    // Add Cell Data Passes
    addValuePasses(poly->GetCellData(), 1);

    // Build Render Sequence
    {
      auto sequence = vtkSmartPointer<vtkSequencePass>::New();
      sequence->SetPasses(valuePassCollection);

      auto cameraPass = vtkSmartPointer<vtkCameraPass>::New();
      cameraPass->SetDelegatePass(sequence);

      auto glRenderer = vtkOpenGLRenderer::SafeDownCast(renderer1);
      glRenderer->SetPass(cameraPass);
    }
  }
  bool renderValuePasses = valuePassList.size() > 0;
  // Render for the first time to initialize everything
  if(renderValuePasses)
    renderWindow1->Render();

#else
  {
    stringstream msg;
    msg << "[ttkCinemaImaging] ERROR: VTK version too old." << endl;
    msg << "[ttkCinemaImaging]        Value images requires VTK 7.0 or higher"
        << endl;
    dMsg(cout, msg.str(), fatalMsg);
  }
#endif

  // Render for the first time to initialize everything
  renderWindow0->Render();

  // Print Status
  {
    stringstream msg;
    t0 = t.getElapsedTime();
    msg << "[ttkCinemaImaging] VTK Rendering Pipeline initialized in " << t0
        << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  // -------------------------------------------------------------------------
  // Render Images for all Camera Locations
  // -------------------------------------------------------------------------

  // Create vtkWindowToImageFilter
  auto windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowToImageFilter->SetInput(renderWindow0);
  windowToImageFilter
    ->SetInputBufferTypeToZBuffer(); // Set output to depth buffer

  // Prepare Field Data for Depth Values
  auto ch = vtkSmartPointer<vtkDoubleArray>::New();
  ch->SetName("CamHeight");
  ch->SetNumberOfValues(1);
  ch->SetValue(0, this->CamHeight);

  auto cnf = vtkSmartPointer<vtkDoubleArray>::New();
  cnf->SetName("CamNearFar");
  cnf->SetNumberOfValues(2);
  cnf->SetValue(0, this->CamNearFar[0]);
  cnf->SetValue(1, this->CamNearFar[1]);

  auto cr = vtkSmartPointer<vtkDoubleArray>::New();
  cr->SetName("CamRes");
  cr->SetNumberOfValues(2);
  cr->SetValue(0, this->Resolution[0]);
  cr->SetValue(1, this->Resolution[1]);

  // Iterate over Locations
  double camPosition[3] = {0, 0, 0};
  size_t n = inputGrid->GetNumberOfPoints();
  auto inputGridPointData = inputGrid->GetPointData();
  size_t nInputGridPointData = inputGridPointData->GetNumberOfArrays();

  for(size_t i = 0; i < n; i++) {
    // Set Camera Position
    // TODO: In the future it shoud be possible to override focus, res,...
    inputGrid->GetPoint(i, camPosition);

    // Cam Up Fix
    if(camPosition[0] == 0 && camPosition[2] == 0) {
      camPosition[0] = 0.00000000001;
      camPosition[2] = 0.00000000001;
    }

    camera->SetPosition(camPosition);

    // Render Depth Image
    renderWindow0->Render();

    // Initialize Output Image
    auto outputImage = vtkSmartPointer<vtkImageData>::New();
    {
      windowToImageFilter->Modified();
      windowToImageFilter->Update();
      outputImage->DeepCopy(windowToImageFilter->GetOutput());
      auto depthValues
        = outputImage->GetPointData()->GetAbstractArray("ImageScalars");
      depthValues->SetName("Depth");
    }

    // Add Field Data
    {
      auto outputImageFD = outputImage->GetFieldData();

      // Camera Parameters
      outputImageFD->AddArray(ch);
      outputImageFD->AddArray(cnf);
      outputImageFD->AddArray(cr);

      // Position
      auto cp = vtkSmartPointer<vtkDoubleArray>::New();
      cp->SetName("CamPosition");
      cp->SetNumberOfValues(3);
      cp->SetValue(0, camPosition[0]);
      cp->SetValue(1, camPosition[1]);
      cp->SetValue(2, camPosition[2]);
      outputImageFD->AddArray(cp);

      // Dir
      auto cd = vtkSmartPointer<vtkDoubleArray>::New();
      cd->SetName("CamDirection");
      cd->SetNumberOfValues(3);
      double tempCD[3] = {this->CamFocus[0] - camPosition[0],
                          this->CamFocus[1] - camPosition[1],
                          this->CamFocus[2] - camPosition[2]};
      vtkMath::Normalize(tempCD);
      cd->SetValue(0, tempCD[0]);
      cd->SetValue(1, tempCD[1]);
      cd->SetValue(2, tempCD[2]);
      outputImageFD->AddArray(cd);

      // Up
      auto upd = camera->GetViewUp();
      auto cu = vtkSmartPointer<vtkDoubleArray>::New();
      cu->SetName("CamUp");
      cu->SetNumberOfValues(3);
      cu->SetValue(0, upd[0]);
      cu->SetValue(1, upd[1]);
      cu->SetValue(2, upd[2]);
      outputImageFD->AddArray(cu);

      for(size_t j = 0; j < nInputGridPointData; j++) {
        auto array = inputGridPointData->GetAbstractArray(j);
        auto newArray
          = vtkSmartPointer<vtkAbstractArray>::Take(array->NewInstance());
        newArray->SetName(array->GetName());
        newArray->SetNumberOfTuples(1);
        newArray->SetNumberOfComponents(array->GetNumberOfComponents());
        array->GetTuples(i, i, newArray);

        outputImageFD->AddArray(newArray);
      }
    }

// Add Point Data
#if VTK_MAJOR_VERSION >= 7
    if(renderValuePasses) {
      // Render Value Passes
      renderWindow1->Render();

      auto outputImagePD = outputImage->GetPointData();
      for(auto &passData : valuePassList) {
        auto data = vtkSmartPointer<vtkFloatArray>::New();
        data->DeepCopy(passData.first->GetFloatImageDataArray(renderer1));
        data->SetName(passData.second.data());
        outputImagePD->AddArray(data);
      }
    }
#endif

    // Add Image to MultiBlock
    outputImages->SetBlock(i, outputImage);

    this->updateProgress(((float)i) / ((float)(n - 1)));
  }

  // Output Performance
  {
    stringstream msg;
    msg << "[ttkCinemaImaging] "
           "-------------------------------------------------------------"
        << endl;
    msg << "[ttkCinemaImaging] " << n << " Images rendered" << endl;
    msg << "[ttkCinemaImaging]   time: " << (t.getElapsedTime() - t0) << " s"
        << endl;
    msg << "[ttkCinemaImaging] memory: " << mem.getElapsedUsage() << " MB"
        << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}
