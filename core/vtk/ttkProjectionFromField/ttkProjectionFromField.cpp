#include <ttkProjectionFromField.h>

#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPoints.h>

vtkStandardNewMacro(ttkProjectionFromField);

ttkProjectionFromField::ttkProjectionFromField() {
  this->setDebugMsgPrefix("ProjectionFromField");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkProjectionFromField::~ttkProjectionFromField() {
}

int ttkProjectionFromField::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkProjectionFromField::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkProjectionFromField::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {
  ttk::Timer t;

  this->printMsg(
    "Projecting Points",
    0, 0,
    this->threadNumber_, ttk::debug::LineMode::REPLACE
  );

  auto input = vtkPointSet::GetData( inputVector[0] );
  auto output = vtkPointSet::GetData( outputVector );
  output->ShallowCopy(input);

  vtkDataArray *inputScalarFieldU = NULL;
  vtkDataArray *inputScalarFieldV = NULL;
  vtkDataArray *textureCoordinates = NULL;

  if(UseTextureCoordinates) {
    textureCoordinates = input->GetPointData()->GetTCoords();
    if(!textureCoordinates){
        this->printErr("Unable to retrieve texture coordinates.");
      return 0;
    }
  } else {
    inputScalarFieldU = this->GetInputArrayToProcess(0, inputVector);
    if(!inputScalarFieldU){
        this->printErr("Unable to retrieve input array 0 (U).");
      return 0;
    }

    inputScalarFieldV = this->GetInputArrayToProcess(1, inputVector);
    if(!inputScalarFieldV){
        this->printErr("Unable to retrieve input array 1 (V).");
      return 0;
    }
  }

  auto points_ = vtkSmartPointer<vtkPoints>::New();
  if(points_->GetNumberOfPoints() != input->GetNumberOfPoints()) {
    points_->SetNumberOfPoints(input->GetNumberOfPoints());
  }

  std::vector<std::vector<double>> points(threadNumber_);
  for(ttk::ThreadId i = 0; i < threadNumber_; i++) {
    points[i].resize(3);
    points[i][2] = 0;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(ttk::SimplexId i = 0; i < input->GetNumberOfPoints(); i++) {

    ttk::ThreadId threadId = 0;

#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    if(UseTextureCoordinates) {
      textureCoordinates->GetTuple(i, points[threadId].data());
    } else {
      points[threadId][0] = inputScalarFieldU->GetComponent(i, 0);
      points[threadId][1] = inputScalarFieldV->GetComponent(i, 0);
    }

    points_->SetPoint(
      i, points[threadId][0], points[threadId][1], points[threadId][2]);
  }

  output->SetPoints(points_);

  this->printMsg(
    "Projecting Points",
    1, t.getElapsedTime(), this->threadNumber_
  );

  return 1;
}
