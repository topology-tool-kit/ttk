#include <ttkCinemaDarkroomCompositing.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkCinemaDarkroomCompositing);

ttkCinemaDarkroomCompositing::ttkCinemaDarkroomCompositing() {
  this->setDebugMsgPrefix("CinemaDarkroomCompositing");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaDarkroomCompositing::~ttkCinemaDarkroomCompositing() {
}

int ttkCinemaDarkroomCompositing::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

template <typename DT>
int computeMask(unsigned char *mask,
                const DT *depth0,
                const DT *depth1,
                const size_t nPixels) {
  for(size_t i = 0; i < nPixels; i++) {
    mask[i] = depth0[i] > depth1[i] ? 1 : 0;
  }
  return 1;
};

template <typename DT>
int compositeArray(DT *array0,
                   const unsigned char *mask,
                   const DT *array1,
                   const size_t nPixels,
                   const size_t nComponents) {
  for(size_t i = 0; i < nPixels; i++) {
    if(mask[i] == 1) {
      for(size_t j = i * nComponents, k = j + nComponents; j < k; j++)
        array0[j] = array1[j];
    }
  }
  return 1;
};

int ttkCinemaDarkroomCompositing::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  auto output = vtkImageData::GetData(outputVector);
  size_t nInputs = inputVector[0]->GetNumberOfInformationObjects();

  ttk::Timer timer;
  this->printMsg("Compositing " + std::to_string(nInputs) + " Images", 0, 0,
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  for(size_t i = 0; i < nInputs; i++) {
    auto input = vtkImageData::GetData(inputVector[0], i);
    if(i == 0) {
      output->DeepCopy(input);
      continue;
    }

    // get depth arrays
    auto depth0Array = this->GetInputArrayToProcess(0, output);
    auto depth1Array = this->GetInputArrayToProcess(0, input);
    size_t nPixels = depth0Array->GetNumberOfTuples();

    // compute mask
    auto mask = vtkSmartPointer<vtkUnsignedCharArray>::New();
    mask->SetNumberOfTuples(nPixels);
    switch(depth0Array->GetDataType()) {
      vtkTemplateMacro(computeMask<VTK_TT>(
        static_cast<unsigned char *>(ttkUtils::GetVoidPointer(mask)),
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(depth0Array)),
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(depth1Array)), nPixels));
    }

    auto outputPD = output->GetPointData();
    auto inputPD = input->GetPointData();

    for(int a = 0; a < outputPD->GetNumberOfArrays(); a++) {
      auto outputArray = outputPD->GetArray(a);
      if(!outputArray)
        continue;
      auto inputArray = inputPD->GetArray(outputArray->GetName());
      if(!inputArray)
        continue;

      switch(outputArray->GetDataType()) {
        vtkTemplateMacro(compositeArray<VTK_TT>(
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputArray)),
          static_cast<unsigned char *>(ttkUtils::GetVoidPointer(mask)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputArray)), nPixels,
          outputArray->GetNumberOfComponents()));
      }
    }

    // perform compositing
    this->printMsg("Compositing " + std::to_string(nInputs) + " Images",
                   (float)(i + 1) / (float)nInputs, timer.getElapsedTime(),
                   this->threadNumber_, ttk::debug::LineMode::REPLACE);
  }

  this->printMsg("Compositing " + std::to_string(nInputs) + " Images", 1,
                 timer.getElapsedTime(), this->threadNumber_);

  return 1;
}
