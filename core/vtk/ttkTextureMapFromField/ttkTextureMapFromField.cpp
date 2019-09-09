#include <array>
#include <ttkTextureMapFromField.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkTextureMapFromField)

  ttkTextureMapFromField::ttkTextureMapFromField() {

  // init
  OnlyUComponent = true;
  OnlyVComponent = false;

  RepeatUTexture = RepeatVTexture = false;
  textureCoordinates_ = NULL;

  UseAllCores = true;
}

ttkTextureMapFromField::~ttkTextureMapFromField() {

  if(textureCoordinates_)
    textureCoordinates_->Delete();
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkTextureMapFromField::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkTextureMapFromField::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkTextureMapFromField] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkTextureMapFromField::doIt(vtkDataSet *input, vtkDataSet *output) {

  Timer t;

  output->ShallowCopy(input);

  vtkDataArray *inputScalarFieldU = NULL;
  vtkDataArray *inputScalarFieldV = NULL;

  if(UComponent.length()) {
    inputScalarFieldU = input->GetPointData()->GetArray(UComponent.data());
  } else {
    inputScalarFieldU = input->GetPointData()->GetArray(0);
  }

  if(!inputScalarFieldU)
    return -1;

  if(VComponent.length()) {
    inputScalarFieldV = input->GetPointData()->GetArray(VComponent.data());
  } else {
    inputScalarFieldV = input->GetPointData()->GetArray(0);
  }

  if(!inputScalarFieldV)
    return -2;

  if(!textureCoordinates_) {
    textureCoordinates_ = vtkFloatArray::New();
    textureCoordinates_->SetNumberOfComponents(2);
    textureCoordinates_->SetName("UV coordinates from field");
  }

  if(textureCoordinates_->GetNumberOfTuples() != output->GetNumberOfPoints()) {
    textureCoordinates_->SetNumberOfTuples(output->GetNumberOfPoints());
  }

  double uRange[2], vRange[2];
  inputScalarFieldU->GetRange(uRange);
  inputScalarFieldV->GetRange(vRange);

  std::vector<std::array<double, 2>> coordinates(threadNumber_);

  SimplexId count = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < output->GetNumberOfPoints(); i++) {

    ThreadId threadId = 0;

#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    if(!needsToAbort()) {

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

      textureCoordinates_->SetTuple(i, coordinates[threadId].data());

      if(debugLevel_ > 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
        {
          if(!(count % (output->GetNumberOfPoints() / 10))) {
            updateProgress((count + 1.0) / output->GetNumberOfPoints());
          }

          count++;
        }
      }
    }
  }

  output->GetPointData()->SetTCoords(textureCoordinates_);

  {
    stringstream msg;
    msg << "[ttkTextureMapFromField] Texture map computed in "
        << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkTextureMapFromField::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  Memory m;

  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkTextureMapFromField] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
