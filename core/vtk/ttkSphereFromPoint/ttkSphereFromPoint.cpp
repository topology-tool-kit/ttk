#include <ttkSphereFromPoint.h>

#include <vtkAppendPolyData.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkSphereFromPoint);

ttkSphereFromPoint::ttkSphereFromPoint() {

  masterAppender_ = NULL;

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  setDebugMsgPrefix("SphereFromPoint");

  vtkWarningMacro("`TTK SphereFromPoint' is now deprecated. Please use "
                  "`TTK IcospheresFromPoint' instead.");
}

ttkSphereFromPoint::~ttkSphereFromPoint() {

  if(masterAppender_)
    masterAppender_->Delete();

  for(SimplexId i = 0; i < (SimplexId)appenderList_.size(); i++) {
    appenderList_[i]->Delete();
  }

  for(SimplexId i = 0; i < (SimplexId)sphereList_.size(); i++) {
    sphereList_[i]->Delete();
  }

  for(SimplexId i = 0; i < (SimplexId)dataArrayList_.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)dataArrayList_[i].size(); j++) {
      dataArrayList_[i][j]->Delete();
    }
  }
}

int ttkSphereFromPoint::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkSphereFromPoint::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkSphereFromPoint::RequestData(vtkInformation *ttkNotUsed(request),
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {

  Timer t;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  if(masterAppender_) {
    masterAppender_->Delete();
    masterAppender_ = NULL;
  }

  for(SimplexId i = 0; i < (SimplexId)appenderList_.size(); i++) {
    appenderList_[i]->Delete();
  }

  masterAppender_ = vtkAppendPolyData::New();
  appenderList_.resize(threadNumber_);
  for(SimplexId i = 0; i < (SimplexId)appenderList_.size(); i++) {
    appenderList_[i] = vtkAppendPolyData::New();
  }

  if((SimplexId)sphereList_.size() > input->GetNumberOfPoints()) {
    for(SimplexId i = input->GetNumberOfPoints();
        i < (SimplexId)sphereList_.size(); i++) {
      sphereList_[i]->Delete();
    }

    sphereList_.resize(input->GetNumberOfPoints());
  } else if((SimplexId)sphereList_.size() < input->GetNumberOfPoints()) {
    SimplexId oldSize = sphereList_.size();

    sphereList_.resize(input->GetNumberOfPoints());

    for(SimplexId i = oldSize; i < (SimplexId)sphereList_.size(); i++) {
      sphereList_[i] = vtkSphereSource::New();
    }
  }

  if(dataArrayList_.size()) {
    for(SimplexId i = 0; i < (SimplexId)dataArrayList_.size(); i++) {
      for(SimplexId j = 0; j < (SimplexId)dataArrayList_[i].size(); j++) {
        dataArrayList_[i][j]->Delete();
      }
      dataArrayList_[i].clear();
    }
  }
  dataArrayList_.resize(threadNumber_);

  vector<vector<double>> p(threadNumber_);

  for(SimplexId i = 0; i < (SimplexId)p.size(); i++) {
    p[i].resize(3);
  }

  if(!input->GetNumberOfPoints())
    return -1;
  // make this function thread-safe
  input->GetPoint(0, p[0].data());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < input->GetNumberOfPoints(); i++) {

    ThreadId threadId = 0;

#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    input->GetPoint(i, p[threadId].data());

    sphereList_[i]->SetCenter(p[threadId][0], p[threadId][1], p[threadId][2]);

    sphereList_[i]->SetRadius(Radius);

    sphereList_[i]->SetThetaResolution(ThetaResolution);
    sphereList_[i]->SetStartTheta(StartTheta);
    sphereList_[i]->SetEndTheta(EndTheta);

    sphereList_[i]->SetPhiResolution(PhiResolution);
    sphereList_[i]->SetStartPhi(StartPhi);
    sphereList_[i]->SetEndPhi(EndPhi);

    sphereList_[i]->Update();

    vtkPolyData *sphereSurface = sphereList_[i]->GetOutput();

    // copy the data values
    for(SimplexId j = 0; j < input->GetPointData()->GetNumberOfArrays(); j++) {

      vtkDataArray *array = input->GetPointData()->GetArray(j);

      if(array != nullptr && array->GetNumberOfComponents() == 1) {

        double value = 0;
        array->GetTuple(i, &value);

        vtkDataArray *dataArray = array->NewInstance();
        if(dataArray != nullptr) {
          dataArray->SetName(array->GetName());
          dataArray->SetNumberOfTuples(sphereSurface->GetNumberOfPoints());
          for(SimplexId k = 0; k < sphereSurface->GetNumberOfPoints(); k++) {
            dataArray->SetTuple(k, &value);
          }
        }
        sphereSurface->GetPointData()->AddArray(dataArray);
        dataArrayList_[threadId].emplace_back(dataArray);
      } else {
        printMsg(
          "Unsupported number of components :(", debug::Priority::DETAIL);
      }
    }

    appenderList_[threadId]->AddInputConnection(
      sphereList_[i]->GetOutputPort());
  }

  for(SimplexId i = 0; i < (SimplexId)appenderList_.size(); i++) {
    if(appenderList_[i]->GetInput())
      masterAppender_->AddInputConnection(appenderList_[i]->GetOutputPort());
  }
  masterAppender_->Update();

  output->ShallowCopy(masterAppender_->GetOutput());

  printMsg(std::to_string(input->GetNumberOfPoints()) + " spheres generated", 1,
           t.getElapsedTime(), threadNumber_);

  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
