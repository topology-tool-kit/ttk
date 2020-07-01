#include <ttkProjectionFromField.h>

#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkProjectionFromField);

ttkProjectionFromField::ttkProjectionFromField() {
  // init
  pointSet_ = vtkSmartPointer<vtkPoints>::New();

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  setDebugMsgPrefix("ProjectionFromField");
}

ttkProjectionFromField::~ttkProjectionFromField() {
}

int ttkProjectionFromField::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  else
    return 0;

  return 1;
}

int ttkProjectionFromField::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;

  return 1;
}

int ttkProjectionFromField::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  Timer t;

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

  if(pointSet_->GetNumberOfPoints() != input->GetNumberOfPoints()) {
    pointSet_->SetNumberOfPoints(input->GetNumberOfPoints());
  }

  vector<vector<double>> points(threadNumber_);
  for(ThreadId i = 0; i < threadNumber_; i++) {
    points[i].resize(3);
    points[i][2] = 0;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < input->GetNumberOfPoints(); i++) {

    ThreadId threadId = 0;

#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    if(UseTextureCoordinates) {
      textureCoordinates->GetTuple(i, points[threadId].data());
    } else {
      points[threadId][0] = inputScalarFieldU->GetComponent(i, 0);
      points[threadId][1] = inputScalarFieldV->GetComponent(i, 0);
    }

    pointSet_->SetPoint(
      i, points[threadId][0], points[threadId][1], points[threadId][2]);
  }

  output->SetPoints(pointSet_);

  printMsg(std::to_string(input->GetNumberOfPoints()) + " points projected", 1,
           t.getElapsedTime(), threadNumber_);

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
