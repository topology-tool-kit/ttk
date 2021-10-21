#include <ttkCinemaDarkroomCompositing.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
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
    info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
    info->Append(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
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
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  const size_t nInputs = inputVector[0]->GetNumberOfInformationObjects();

  auto inputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  auto outputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  if(nInputs < 1) {
    this->printWrn("Empty Input");
    return 1;
  }

  auto firstInput = vtkDataObject::GetData(inputVector[0], 0);

  if(firstInput->IsA("vtkImageData")) {
    // if input are repeated vtkImageData objects
    for(size_t i = 0; i < nInputs; i++) {
      auto collection = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      auto image = vtkImageData::GetData(inputVector[0], i);
      if(!image) {
        this->printErr("Input is not a list of vtkImageData objects.");
        return 0;
      }
      collection->SetBlock(0, image);
      inputAsMB->SetBlock(i, collection);
    }
  } else if(nInputs == 1 && firstInput->IsA("vtkMultiBlockDataSet")) {
    // if input is a single list of vtkImageData objects
    auto temp = vtkMultiBlockDataSet::SafeDownCast(firstInput);
    auto nBlocks = temp->GetNumberOfBlocks();

    if(nBlocks > 0) {
      if(vtkMultiBlockDataSet::SafeDownCast(temp->GetBlock(0))) {
        inputAsMB->ShallowCopy(firstInput);
      } else {
        // check if all blocks are images
        for(size_t i = 0; i < nBlocks; i++) {
          auto collection = vtkSmartPointer<vtkMultiBlockDataSet>::New();
          auto image = vtkImageData::SafeDownCast(temp->GetBlock(i));
          if(!image) {
            this->printErr(
              "Input is not a single list of vtkImageData objects.");
            return 0;
          }
          collection->SetBlock(0, image);
          inputAsMB->SetBlock(i, collection);
        }
      }
    }
  } else if(firstInput->IsA("vtkMultiBlockDataSet")) {
    // inputs are multiple lists of vtkImageData objects
    for(size_t i = 0; i < nInputs; i++) {
      auto list = vtkMultiBlockDataSet::GetData(inputVector[0], i);
      if(!list) {
        this->printErr(
          "Inputs are not multiple lists of vtkImageData objects.");
        return 0;
      }
      inputAsMB->SetBlock(i, list);
    }
  } else {
    this->printErr("Unsupported input data structure.");
    return 0;
  }

  const size_t nImagesPerCompositingStep = inputAsMB->GetNumberOfBlocks();
  const size_t nCompositingSteps
    = vtkMultiBlockDataSet::SafeDownCast(inputAsMB->GetBlock(0))
        ->GetNumberOfBlocks();

  ttk::Timer timer;
  this->printMsg("Compositing " + std::to_string(nCompositingSteps) + "x"
                   + std::to_string(nImagesPerCompositingStep) + " Images",
                 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

  auto getImage = [](vtkMultiBlockDataSet *root, int i, int j) {
    return (vtkImageData *)((vtkMultiBlockDataSet *)root->GetBlock(i))
      ->GetBlock(j);
  };

  for(size_t i = 0; i < nCompositingSteps; i++) {

    auto firstImage = getImage(inputAsMB, 0, i);

    auto output = vtkSmartPointer<vtkImageData>::New();
    output->DeepCopy(firstImage);

    for(size_t j = 1; j < nImagesPerCompositingStep; j++) {
      auto image = getImage(inputAsMB, j, i);
      if(!image) {
        this->printErr("Unsupported input data structure.");
        return 0;
      }

      // get depth arrays
      auto depth0Array = this->GetInputArrayToProcess(0, output);
      auto depth1Array = this->GetInputArrayToProcess(0, image);
      if(!depth0Array || !depth1Array) {
        this->printErr("Unable to retrieve depth arrays.");
        return 0;
      }

      const size_t nPixels = depth0Array->GetNumberOfTuples();

      // compute mask
      auto mask = vtkSmartPointer<vtkUnsignedCharArray>::New();
      mask->SetNumberOfTuples(nPixels);
      switch(depth0Array->GetDataType()) {
        vtkTemplateMacro(computeMask<VTK_TT>(
          static_cast<unsigned char *>(ttkUtils::GetVoidPointer(mask)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(depth0Array)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(depth1Array)),
          nPixels));
      }

      auto outputPD = output->GetPointData();
      auto inputPD = image->GetPointData();

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
            static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputArray)),
            nPixels, outputArray->GetNumberOfComponents()));
        }
      }
    }

    outputAsMB->SetBlock(i, output);

    // perform compositing
    this->printMsg("Compositing " + std::to_string(nCompositingSteps) + "x"
                     + std::to_string(nImagesPerCompositingStep) + " Images",
                   (float)(i + 1) / (float)nInputs, timer.getElapsedTime(),
                   this->threadNumber_, ttk::debug::LineMode::REPLACE);
  }

  this->printMsg("Compositing " + std::to_string(nCompositingSteps) + "x"
                   + std::to_string(nImagesPerCompositingStep) + " Images",
                 1, timer.getElapsedTime(), this->threadNumber_);

  auto output = vtkDataObject::GetData(outputVector);
  if(firstInput->IsA("vtkImageData")) {
    output->ShallowCopy(outputAsMB->GetBlock(0));
  } else {
    output->ShallowCopy(outputAsMB);
  }

  return 1;
}
