#include <ttkProjectionFromField.h>

#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkProjectionFromField);

ttkProjectionFromField::ttkProjectionFromField() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  this->setDebugMsgPrefix("ProjectionFromField");
}

int ttkProjectionFromField::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkProjectionFromField::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkProjectionFromField::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  ttk::Timer t;

  vtkPointSet *input = vtkPointSet::GetData(inputVector[0]);
  vtkPointSet *output = vtkPointSet::GetData(outputVector, 0);

  output->ShallowCopy(input);

  vtkDataArray *inputScalarFieldU = nullptr;
  vtkDataArray *inputScalarFieldV = nullptr;
  vtkDataArray *textureCoordinates = nullptr;

  if(UseTextureCoordinates) {
    textureCoordinates = input->GetPointData()->GetTCoords();

    if(!textureCoordinates)
      return -3;

    printMsg("Starting computation with texture coordinates...");
  } else {

    inputScalarFieldU = this->GetInputArrayToProcess(0, inputVector);
    inputScalarFieldV = this->GetInputArrayToProcess(1, inputVector);

    if(!inputScalarFieldU)
      return -1;

    if(!inputScalarFieldV)
      return -2;

    printMsg("Starting computation...");
    printMsg({{"  U-component", inputScalarFieldU->GetName()},
              {"  V-component", inputScalarFieldV->GetName()}});
  }

  vtkNew<vtkPoints> pointSet{};

  if(pointSet->GetNumberOfPoints() != input->GetNumberOfPoints()) {
    pointSet->SetNumberOfPoints(input->GetNumberOfPoints());
  }

  std::vector<std::vector<double>> points(threadNumber_);

  for(int i = 0; i < threadNumber_; i++) {
    points[i].resize(3);
    points[i][2] = 0;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < input->GetNumberOfPoints(); i++) {

#ifdef TTK_ENABLE_OPENMP
    const int threadId = omp_get_thread_num();
#else
    const int threadId = 0;
#endif

    if(UseTextureCoordinates) {
      textureCoordinates->GetTuple(i, points[threadId].data());
    } else {
      points[threadId][0] = inputScalarFieldU->GetComponent(i, 0);
      points[threadId][1] = inputScalarFieldV->GetComponent(i, 0);
    }

    pointSet->SetPoint(
      i, points[threadId][0], points[threadId][1], points[threadId][2]);
  }

  output->SetPoints(pointSet);

  printMsg(std::to_string(input->GetNumberOfPoints()) + " points projected", 1,
           t.getElapsedTime(), threadNumber_);

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
