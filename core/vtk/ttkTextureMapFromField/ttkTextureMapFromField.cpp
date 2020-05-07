#include <ttkTextureMapFromField.h>

#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <array>

vtkStandardNewMacro(ttkTextureMapFromField);

ttkTextureMapFromField::ttkTextureMapFromField() {
  this->setDebugMsgPrefix("TextureMapFromField");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkTextureMapFromField::~ttkTextureMapFromField() {
}

int ttkTextureMapFromField::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTextureMapFromField::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkTextureMapFromField::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  ttk::Timer t;

  this->printMsg("Computing Texture Map", 0, 0, ttk::debug::LineMode::REPLACE);

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

  output->ShallowCopy(input);

  auto inputScalarFieldU = this->GetInputArrayToProcess(0, inputVector);
  if(!inputScalarFieldU)
    return -1;
  auto inputScalarFieldV = this->GetInputArrayToProcess(1, inputVector);
  if(!inputScalarFieldV)
    return -2;

  auto textureCoordinates = vtkSmartPointer<vtkFloatArray>::New();
  textureCoordinates->SetNumberOfComponents(2);
  textureCoordinates->SetName("UV coordinates from field");
  textureCoordinates->SetNumberOfTuples(output->GetNumberOfPoints());

  double uRange[2], vRange[2];
  inputScalarFieldU->GetRange(uRange);
  inputScalarFieldV->GetRange(vRange);

  std::vector<std::array<double, 2>> coordinates(threadNumber_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(ttk::SimplexId i = 0; i < output->GetNumberOfPoints(); i++) {

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

  this->printMsg("Computing Texture Map", 1, t.getElapsedTime());

  return 1;
}