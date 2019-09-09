#include <ttkProjectionFromField.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkProjectionFromField)

  ttkProjectionFromField::ttkProjectionFromField() {
  UseAllCores = true;

  // init
  pointSet_ = vtkSmartPointer<vtkPoints>::New();
}

ttkProjectionFromField::~ttkProjectionFromField() {
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkProjectionFromField::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkProjectionFromField::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkProjectionFromField] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkProjectionFromField::doIt(vtkPointSet *input, vtkPointSet *output) {

  Timer t;

  output->ShallowCopy(input);

  vtkDataArray *inputScalarFieldU = NULL;
  vtkDataArray *inputScalarFieldV = NULL;
  vtkDataArray *textureCoordinates = NULL;

  if(UseTextureCoordinates) {
    textureCoordinates = input->GetPointData()->GetTCoords();

    if(!textureCoordinates)
      return -1;
  } else {

    if(UComponent.length()) {
      inputScalarFieldU = input->GetPointData()->GetArray(UComponent.data());
    } else {
      inputScalarFieldU = input->GetPointData()->GetArray(0);
    }

    if(!inputScalarFieldU)
      return -2;

    if(VComponent.length()) {
      inputScalarFieldV = input->GetPointData()->GetArray(VComponent.data());
    } else {
      inputScalarFieldV = input->GetPointData()->GetArray(0);
    }

    if(!inputScalarFieldV)
      return -3;
  }

  if(pointSet_->GetNumberOfPoints() != input->GetNumberOfPoints()) {
    pointSet_->SetNumberOfPoints(input->GetNumberOfPoints());
  }

  vector<vector<double>> points(threadNumber_);
  for(ThreadId i = 0; i < threadNumber_; i++) {
    points[i].resize(3);
    points[i][2] = 0;
  }

  SimplexId count = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < input->GetNumberOfPoints(); i++) {

    ThreadId threadId = 0;

#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    if(!needsToAbort()) {

      if(UseTextureCoordinates) {
        textureCoordinates->GetTuple(i, points[threadId].data());
      } else {
        points[threadId][0] = inputScalarFieldU->GetComponent(i, 0);
        points[threadId][1] = inputScalarFieldV->GetComponent(i, 0);
      }

      pointSet_->SetPoint(
        i, points[threadId][0], points[threadId][1], points[threadId][2]);

      if(debugLevel_ > Debug::advancedInfoMsg) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
        {
          if(!(count % (input->GetNumberOfPoints() / 10))) {
            updateProgress((count + 1.0) / input->GetNumberOfPoints());
          }

          count++;
        }
      }
    }
  }

  output->SetPoints(pointSet_);

  {
    stringstream msg;
    msg << "[ttkProjectionFromField] Data-set projected in "
        << t.getElapsedTime() << " s. (" << input->GetNumberOfPoints()
        << " points)." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkProjectionFromField::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  Memory m;

  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkPointSet *input = vtkPointSet::GetData(inputVector[0]);
  vtkPointSet *output = vtkPointSet::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkProjectionFromField] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
