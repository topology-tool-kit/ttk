#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>

#include <array>

#include <ttkTextureMapFromField.h>

vtkStandardNewMacro(ttkTextureMapFromField);

ttkTextureMapFromField::ttkTextureMapFromField() {
  this->setDebugMsgPrefix("TextureMapFromField");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkTextureMapFromField::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
};
int ttkTextureMapFromField::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
};

int ttkTextureMapFromField::RequestData(vtkInformation *ttkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

  ttk::Timer t;

  output->ShallowCopy(input);

  const auto inputScalarFieldU = this->GetInputArrayToProcess(0, input);
  const auto inputScalarFieldV = this->GetInputArrayToProcess(1, input);

  if(inputScalarFieldU == nullptr || inputScalarFieldV == nullptr) {
    return -1;
  }

  vtkNew<vtkFloatArray> textureCoordinates{};
  textureCoordinates->SetNumberOfComponents(2);
  textureCoordinates->SetName("UV coordinates from field");

  if(textureCoordinates->GetNumberOfTuples() != output->GetNumberOfPoints()) {
    textureCoordinates->SetNumberOfTuples(output->GetNumberOfPoints());
  }

  double uRange[2], vRange[2];
  inputScalarFieldU->GetRange(uRange);
  inputScalarFieldV->GetRange(vRange);

  std::vector<std::array<double, 2>> coordinates(threadNumber_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < output->GetNumberOfPoints(); i++) {

    ttk::ThreadId threadId = 0;

#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    coordinates[threadId][0] = coordinates[threadId][1] = 0;

    if(!OnlyVComponent) {
      inputScalarFieldU->GetTuple(i, &(coordinates[threadId][0]));
      if(!RepeatUTexture) {
        coordinates[threadId][0]
          = (coordinates[threadId][0] - uRange[0]) / (uRange[1] - uRange[0]);
      }
    }

    if(!OnlyUComponent) {
      inputScalarFieldV->GetTuple(i, &(coordinates[threadId][1]));
      if(!RepeatVTexture) {
        coordinates[threadId][1]
          = (coordinates[threadId][1] - vRange[0]) / (vRange[1] - vRange[0]);
      }
    }

    textureCoordinates->SetTuple(i, coordinates[threadId].data());
  }

  output->GetPointData()->SetTCoords(textureCoordinates);

  this->printMsg(
    "Computed texture map", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 1;
}
