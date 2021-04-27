#include <ttkCinemaDarkroomShader.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>

#include <vtkCamera.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkOpenGLRenderer.h>
#include <vtkPlaneSource.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>

#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>

#include <vtkOpenGLTexture.h>
#include <vtkProperty.h>
#include <vtkShaderProperty.h>
#include <vtkTextureObject.h>

#include <vtkFramebufferPass.h>
#include <vtkRenderStepsPass.h>

#include <boost/algorithm/string.hpp>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkCinemaDarkroomShader);

ttkCinemaDarkroomShader::ttkCinemaDarkroomShader() {
  this->CreateFullScreenQuad();
  this->CreateRenderer();

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaDarkroomShader::~ttkCinemaDarkroomShader() {
}

int ttkCinemaDarkroomShader::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomShader::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

std::string ttkCinemaDarkroomShader::GetVertexShaderCode() {
  return std::string(R"(

//VTK::System::Dec  // always start with this line in your VS

attribute vec4 vertexMC;
varying vec4 vPos;

void main () {
    vPos = vertexMC/2. + vec4(0.5,0.5,0.5,0);
    gl_Position = vertexMC;
}

  )");
}

std::string ttkCinemaDarkroomShader::GetFragmentShaderCode() {
  return std::string(R"(

//VTK::System::Dec // always start with these lines in your FS
//VTK::Output::Dec // always start with these lines in your FS

varying vec4 vPos;

void main(void) {
    gl_FragData[0] = vec4(1,0,0,1);
}

  )");
}

std::string
  ttkCinemaDarkroomShader::PerformReplacements(const std::string &source) {
  std::string result = source;

  for(const auto &it : this->Replacements)
    boost::replace_all(result, it.first, it.second.toString());

  return result;
}

int ttkCinemaDarkroomShader::AddReplacement(const std::string &name,
                                            const std::vector<double> &values,
                                            const bool &isInt) {
  auto it = this->Replacements.find(name);
  if(it == this->Replacements.end()) {
    this->Replacements.emplace(std::piecewise_construct,
                               std::forward_as_tuple(name),
                               std::forward_as_tuple(values, isInt));
  } else {
    it->second.values = values;
  }

  return 1;
}

int ttkCinemaDarkroomShader::CreateFullScreenQuad() {
  auto ps = vtkSmartPointer<vtkPlaneSource>::New();
  ps->SetOrigin(-1, -1, 0);
  ps->SetPoint1(1, -1, 0);
  ps->SetPoint2(-1, 1, 0);
  ps->Update();

  this->FullScreenQuad = vtkSmartPointer<vtkPolyData>::New();
  this->FullScreenQuad->ShallowCopy(ps->GetOutput());
  this->FullScreenQuadActor = vtkSmartPointer<vtkActor>::New();

  auto mapper = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
  mapper->SetInputData(this->FullScreenQuad);
  this->FullScreenQuadActor->SetMapper(mapper);

  return 1;
}

int ttkCinemaDarkroomShader::CreateRenderer() {
  // Renderer
  this->Renderer = vtkSmartPointer<vtkOpenGLRenderer>::New();
  this->Renderer->AddActor(this->FullScreenQuadActor);
  this->Renderer->SetBackground(0, 0, 0);

  // Camera
  auto camera = vtkSmartPointer<vtkCamera>::New();
  camera->SetParallelProjection(true);
  camera->SetClippingRange(0, 2);
  camera->SetPosition(0, 0, 1);
  camera->SetFocalPoint(0, 0, 0);
  camera->SetParallelScale(
    1); // Will be ignored because quad positions are fixed
  this->Renderer->SetActiveCamera(camera);

  this->RenderWindow = vtkSmartPointer<vtkRenderWindow>::New();

  return 1;
}

int ttkCinemaDarkroomShader::InitRenderer(vtkImageData *outputImage) {
  int dim[3];
  outputImage->GetDimensions(dim);

  this->AddReplacement("cResolution", {(double)dim[0], (double)dim[1]});

  // Window
  int *size = this->RenderWindow->GetSize();
  if(size[0] != dim[0] || size[1] != dim[1]) {
    ttk::Timer timer;
    this->printMsg("Initializing Renderer (" + std::to_string(dim[0]) + "x"
                     + std::to_string(dim[1]) + ")",
                   0, 0, 1, ttk::debug::LineMode::REPLACE);

    this->RenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    this->RenderWindow->AddRenderer(this->Renderer);
    this->RenderWindow->SetMultiSamples(0); // Disable AA
    this->RenderWindow->OffScreenRenderingOn();

    auto windowAsOGL = vtkOpenGLRenderWindow::SafeDownCast(this->RenderWindow);
    if(windowAsOGL == nullptr) {
      return 0;
    }
    windowAsOGL->SetSize(dim[0], dim[1]);
    windowAsOGL->Initialize();

    this->printMsg("Initializing Renderer (" + std::to_string(dim[0]) + "x"
                     + std::to_string(dim[1]) + ")",
                   1, timer.getElapsedTime(), 1);
  }

  return 1;
}

int ttkCinemaDarkroomShader::AddTexture(vtkImageData *image,
                                        int arrayIdx,
                                        int textureIdx) {
  int dim[3];
  image->GetDimensions(dim);

  auto inputArray = this->GetInputArrayToProcess(arrayIdx, image);
  if(!inputArray || this->GetInputArrayAssociation(arrayIdx, image) != 0) {
    this->printErr("Unable to retrieve input point data array "
                   + std::to_string(arrayIdx) + ".");
    return 0;
  }

  std::string textureName = "tex" + std::to_string(textureIdx);

  auto properties = this->FullScreenQuadActor->GetProperty();
  // if texture already exists remove it
  if(properties->GetTexture(textureName.data()))
    properties->RemoveTexture(textureName.data());

  // create texture
  auto texture = vtkSmartPointer<vtkOpenGLTexture>::New();

  auto textureObj = vtkSmartPointer<vtkTextureObject>::New();
  textureObj->SetContext(
    vtkOpenGLRenderWindow::SafeDownCast(this->RenderWindow));
  textureObj->SetWrapT(vtkTextureObject::ClampToEdge);
  textureObj->SetWrapS(vtkTextureObject::ClampToEdge);
  textureObj->Create2DFromRaw(
    dim[0], dim[1], inputArray->GetNumberOfComponents(),
    inputArray->GetDataType(), ttkUtils::GetVoidPointer(inputArray));
  texture->SetTextureObject(textureObj);
  texture->InterpolateOn();

  // add texture to properties
  properties->SetTexture(textureName.data(), texture);

  return 1;
}

int ttkCinemaDarkroomShader::Render(vtkImageData *image,
                                    const std::string &name) {
  ttk::Timer timer;
  int dim[3];
  image->GetDimensions(dim);
  this->printMsg(
    "Rendering (" + std::to_string(dim[0]) + "x" + std::to_string(dim[1]) + ")",
    0, 0, 1, ttk::debug::LineMode::REPLACE);

  this->FullScreenQuadActor->GetShaderProperty()->SetVertexShaderCode(
    this->PerformReplacements(this->GetVertexShaderCode()).data());
  this->FullScreenQuadActor->GetShaderProperty()->SetFragmentShaderCode(
    this->PerformReplacements(this->GetFragmentShaderCode()).data());

  auto buffer = vtkSmartPointer<vtkUnsignedCharArray>::New();
  buffer->SetName(name.data());
  buffer->SetNumberOfComponents(4);
  buffer->SetNumberOfTuples(dim[0] * dim[1]);

  this->RenderWindow->Render();
  this->RenderWindow->GetRGBACharPixelData(
    0, 0, dim[0] - 1, dim[1] - 1, 1, buffer);

  image->GetPointData()->AddArray(buffer);
  image->GetPointData()->SetActiveScalars(buffer->GetName());

  this->printMsg(
    "Rendering (" + std::to_string(dim[0]) + "x" + std::to_string(dim[1]) + ")",
    1, timer.getElapsedTime(), 1);

  return 1;
}
