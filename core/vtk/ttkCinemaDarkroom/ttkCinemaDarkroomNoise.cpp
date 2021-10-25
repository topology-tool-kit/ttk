#include <ttkCinemaDarkroomNoise.h>

#include <ttkUtils.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkCinemaDarkroomNoise);

ttkCinemaDarkroomNoise::ttkCinemaDarkroomNoise() : ttkCinemaDarkroomShader() {
  this->setDebugMsgPrefix("CinemaDarkroomNoise");
}
ttkCinemaDarkroomNoise::~ttkCinemaDarkroomNoise() {
}

int ttkCinemaDarkroomNoise::RequestData(vtkInformation *ttkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  auto inputImage = vtkImageData::GetData(inputVector[0]);
  auto outputImage = vtkImageData::GetData(outputVector);
  outputImage->ShallowCopy(inputImage);

  size_t nPoints = outputImage->GetNumberOfPoints();

  auto noise = vtkSmartPointer<vtkFloatArray>::New();
  noise->SetName("Noise");
  noise->SetNumberOfTuples(nPoints);
  auto noiseData = static_cast<float *>(ttkUtils::GetVoidPointer(noise));

  int dim[3];
  outputImage->GetDimensions(dim);

  ttk::Timer timer;
  const std::string msg = "Computing Noise (" + std::to_string(dim[0]) + "x"
                          + std::to_string(dim[1]) + "x"
                          + std::to_string(dim[2]) + ")";

  this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE);

  for(size_t i = 0; i < nPoints; i++)
    noiseData[i] = ((float)std::rand()) / ((float)RAND_MAX);

  this->printMsg(msg, 0, 0, 1);

  outputImage->GetPointData()->AddArray(noise);

  return 1;
}