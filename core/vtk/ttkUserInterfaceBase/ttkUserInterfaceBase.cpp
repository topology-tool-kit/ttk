// local includes
#include <ttkUserInterfaceBase.h>

#include <vtkPointData.h>
#include <vtkTexture.h>

#ifndef TTK_INSTALL_ASSETS_DIR
#define TTK_INSTALL_ASSETS_DIR "."
#endif // TTK_INSTALL_ASSETS_DIR

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCustomInteractor);

void ttkCustomInteractor::OnKeyPress() {

  vtkRenderWindowInteractor *interactor = this->Interactor;
  string key = interactor->GetKeySym();

  if(interactor->GetControlKey()) {
    // control key pressed
    if(key == "e") {
      // export
      userInterface_->exportScene("output.wrl");
    }
    if(key == "s") {
      // save in vtk file format
      userInterface_->save();
    }
  } else {

    if(key == "Return") {
      // do something if the user hits Enter
      userInterface_->execute();
    } else if(key == "Escape") {
      exit(0);
    }

    if((key == "0") || (key == "1") || (key == "2") || (key == "3")
       || (key == "4") || (key == "5") || (key == "6") || (key == "7")
       || (key == "8") || (key == "9")) {

      int outputId = -1;
      stringstream stream;

      stream << key;
      stream >> outputId;

      userInterface_->switchOutput(outputId);
    }

    if(key == "t") {
      userInterface_->switchTransparency();
    }
  }

  if(userInterface_->getKeyHandler()) {
    userInterface_->getKeyHandler()->OnKeyPress(interactor, key);
  }

  userInterface_->refresh();
}

ttkUserInterfaceBase::ttkUserInterfaceBase() {

  keyHandler_ = NULL;
  vtkWrapper_ = NULL;
  isUp_ = false;
  repeat_ = false;
  transparency_ = false;
  fullscreen_ = false;

  customInteractor_ = vtkSmartPointer<ttkCustomInteractor>::New();
  pngReader_ = vtkSmartPointer<vtkPNGReader>::New();
  renderer_ = vtkSmartPointer<vtkRenderer>::New();
  renderWindow_ = vtkSmartPointer<vtkRenderWindow>::New();
  interactor_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  texture_ = vtkSmartPointer<vtkTexture>::New();

  pngReader_->SetFileName(
    TTK_INSTALL_ASSETS_DIR
    "/textures/png/scalarFieldTexturePaleInterleavedRules.png");
  pngReader_->Update();
  hasTexture_ = !((pngReader_->GetOutput()->GetNumberOfPoints() == 1)
                  && (pngReader_->GetOutput()->GetNumberOfCells() == 1));
}

ttkUserInterfaceBase::~ttkUserInterfaceBase() {
}

int ttkUserInterfaceBase::exportScene(const string &fileName) const {

  vtkVRMLExporter *exporter = vtkVRMLExporter::New();

  exporter->SetInput(renderWindow_);
  exporter->SetFileName(fileName.data());
  exporter->Write();

  exporter->Delete();

  return 0;
}

int ttkUserInterfaceBase::init(int &argc, char **argv) {

  parser_.setOption("R", &repeat_, "Repeat the program when hitting `Return'");
  parser_.setOption(
    "fullscreen", &fullscreen_, "Maximize the window at launch");

  return ProgramBase::init(argc, argv);
}

int ttkUserInterfaceBase::refresh() {

  // collect the output and update the rendering
  int outputPortNumber = vtkWrapper_->GetNumberOfOutputPorts();

  if((int)visibleOutputs_.size() != outputPortNumber) {
    visibleOutputs_.resize(outputPortNumber, true);

    for(int i = 0; i < (int)hiddenOutputs_.size(); i++) {
      if((hiddenOutputs_[i] >= 0)
         && (hiddenOutputs_[i] < (int)visibleOutputs_.size())) {
        visibleOutputs_[hiddenOutputs_[i]] = false;
      }
    }
  }

  if((int)surfaces_.size() != outputPortNumber) {
    surfaces_.resize(outputPortNumber, NULL);
    mainActors_.resize(outputPortNumber);
    boundaryFilters_.resize(outputPortNumber);
    boundaryMappers_.resize(outputPortNumber);
    textureMapFromFields_.resize(outputPortNumber);

    for(int i = 0; i < outputPortNumber; i++) {
      mainActors_[i] = vtkSmartPointer<vtkActor>::New();
      boundaryFilters_[i] = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
      boundaryMappers_[i] = vtkSmartPointer<vtkPolyDataMapper>::New();
      textureMapFromFields_[i] = vtkSmartPointer<ttkTextureMapFromField>::New();
    }
  }

  for(int i = 0; i < outputPortNumber; i++) {

    if(visibleOutputs_[i]) {

      if(hasTexture_) {
        vtkWrapper_->GetOutput(i)->GetPointData()->SetActiveScalars(NULL);
      }

      if((repeat_) && (i < vtkWrapper_->GetNumberOfInputPorts())) {
        inputs_[i]->DeepCopy(vtkWrapper_->GetOutput(i));
        vtkWrapper_->SetInputData(i, inputs_[i]);
      }

      // extract the boundary surface of the tet-mesh
      boundaryFilters_[i]->SetInputData(vtkWrapper_->GetOutput(i));
      boundaryFilters_[i]->Update();
    }
  }

  // update the scalar field texture
  updateScalarFieldTexture();

  for(int i = 0; i < (int)mainActors_.size(); i++) {
    if(visibleOutputs_[i]) {
      mainActors_[i]->SetMapper(boundaryMappers_[i]);
    } else {
      mainActors_[i]->SetMapper(NULL);
    }
    if(transparency_) {
      mainActors_[i]->GetProperty()->SetOpacity(0.3);
    } else {
      mainActors_[i]->GetProperty()->SetOpacity(1);
    }
    renderer_->AddActor(mainActors_[i]);
  }

  if(isUp_)
    renderWindow_->Render();

  return 0;
}

int ttkUserInterfaceBase::run() {

  execute();

  {
    stringstream msg;
    msg << "[UserInterace] Initializing user interface..." << endl;
    printMsg(msg.str());
  }

  renderWindow_->AddRenderer(renderer_);
  if(fullscreen_) {
    renderWindow_->SetFullScreen(fullscreen_);
  } else {
    renderWindow_->SetSize(1920, 1080);
  }
  renderWindow_->SetWindowName("TTK - The Topology ToolKit");

  interactor_->SetRenderWindow(renderWindow_);
  interactor_->SetInteractorStyle(customInteractor_);
  customInteractor_->SetCurrentRenderer(renderer_);
  customInteractor_->setUserInterface(this);

  refresh();

  renderer_->SetBackground(44 / 255.0, 44 / 255.0, 44 / 255.0);

  interactor_->Initialize();

  {
    stringstream msg;
    msg << "[ttkUserInterfaceBase] Running user interface!" << endl;
    printMsg(msg.str());
  }

  isUp_ = true;

  interactor_->Start();

  return 0;
}

int ttkUserInterfaceBase::switchOutput(const int &outputId) {

  if(!vtkWrapper_)
    return -1;

  if((outputId < 0) || (outputId >= vtkWrapper_->GetNumberOfOutputPorts()))
    return -2;

  stringstream msg;
  msg << "[ttkUserInterfaceBase] Turning output #" << outputId << " ";

  if(visibleOutputs_[outputId]) {
    msg << "off";
  } else {
    msg << "on";
  }
  msg << endl;
  printMsg(msg.str());

  visibleOutputs_[outputId] = !visibleOutputs_[outputId];

  return 0;
}

int ttkUserInterfaceBase::switchTransparency() {

  stringstream msg;
  msg << "[ttkUserInterfaceBase] Switching transparency ";

  if(transparency_) {
    msg << "off";
  } else {
    msg << "on";
  }
  msg << endl;
  printMsg(msg.str());

  transparency_ = !transparency_;

  return 0;
}

int ttkUserInterfaceBase::updateScalarFieldTexture() {

  for(int i = 0; i < (int)boundaryFilters_.size(); i++) {

    // use a texture for the color-value visualization
    if((boundaryFilters_[i]->GetOutput()->GetPointData())
       && (boundaryFilters_[i]->GetOutput()->GetPointData()->GetArray(0))) {

      textureMapFromFields_[i]->SetInputDataObject(
        0, boundaryFilters_[i]->GetOutput());
      textureMapFromFields_[i]->Update();
      surfaces_[i] = vtkPolyData::SafeDownCast(
        textureMapFromFields_[i]->GetOutputDataObject(0));

      texture_->SetInputConnection(pngReader_->GetOutputPort());

      if(hasTexture_) {
        mainActors_[i]->SetTexture(texture_);
      }
    } else {
      surfaces_[i] = boundaryFilters_[i]->GetOutput();
    }

    boundaryMappers_[i]->SetInputData(surfaces_[i]);
  }

  return 0;
}
