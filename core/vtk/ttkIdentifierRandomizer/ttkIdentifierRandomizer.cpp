#include <Shuffle.h>
#include <ttkIdentifierRandomizer.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <map>
#include <numeric>
#include <random>

vtkStandardNewMacro(ttkIdentifierRandomizer);

ttkIdentifierRandomizer::ttkIdentifierRandomizer() {

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  this->setDebugMsgPrefix("IdentifierRandomizer");
}

int ttkIdentifierRandomizer::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
    return 1;
  }
  return 0;
}

int ttkIdentifierRandomizer::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <typename T>
int shuffleScalarFieldValues(const T *const inputField,
                             T *const outputField,
                             const int nValues,
                             const int seed,
                             const bool compactRange,
                             const int nThreads = 1) {

  // copy input field into vector
  std::vector<T> inputValues(inputField, inputField + nValues);

  // reduce the copy
  TTK_PSORT(nThreads, inputValues.begin(), inputValues.end());
  const auto last = std::unique(inputValues.begin(), inputValues.end());
  inputValues.erase(last, inputValues.end());

  // copy the range of values
  std::vector<T> shuffledValues(inputValues.size());
  if(compactRange) {
    std::iota(shuffledValues.begin(), shuffledValues.end(), T{});
  } else {
    std::copy(inputValues.begin(), inputValues.end(), shuffledValues.begin());
  }

  // shuffle them using the seed
  std::mt19937 random_engine{};
  random_engine.seed(seed);
  // use the Fisher-Yates algorithm instead of std::shuffle, whose
  // results are platform-dependent
  ttk::shuffle(shuffledValues, random_engine);

  // link original value to shuffled value correspondence
  std::map<T, T> originalToShuffledValues{};
  for(size_t i = 0; i < inputValues.size(); ++i) {
    originalToShuffledValues[inputValues[i]] = shuffledValues[i];
  }

// write shuffled values inside the output scalar field
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < nValues; ++i) {
    outputField[i] = originalToShuffledValues[inputField[i]];
  }

  TTK_FORCE_USE(nThreads);
  return 1;
}

template <typename T>
int ttkIdentifierRandomizer::shuffleScalarFieldValuesMultiBlock(
  vtkMultiBlockDataSet *input,
  vtkMultiBlockDataSet *output,
  const int nThreads) {

  int n_blocks = input->GetNumberOfBlocks();
  std::vector<T> inputValues;
  for(int i = 0; i < n_blocks; i++) {
    vtkDataSet *block = vtkDataSet::SafeDownCast(input->GetBlock(i));
    vtkDataArray *inputScalarField = this->GetInputArrayToProcess(0, block);
    int nValues = inputScalarField->GetNumberOfTuples();
    const T *const inputScalarFieldPtr
      = static_cast<T *>(ttkUtils::GetVoidPointer(inputScalarField));
    inputValues.insert(
      inputValues.end(), inputScalarFieldPtr, inputScalarFieldPtr + nValues);
  }

  TTK_PSORT(nThreads, inputValues.begin(), inputValues.end());
  const auto last = std::unique(inputValues.begin(), inputValues.end());
  inputValues.erase(last, inputValues.end());
  std::vector<T> shuffledValues(inputValues.size());
  if(CompactRange) {
    std::iota(shuffledValues.begin(), shuffledValues.end(), T{});
  } else {
    std::copy(inputValues.begin(), inputValues.end(), shuffledValues.begin());
  }
  // shuffle them using the seed
  std::mt19937 random_engine{};
  random_engine.seed(RandomSeed);
  // use the Fisher-Yates algorithm instead of std::shuffle, whose
  // results are platform-dependent
  ttk::shuffle(shuffledValues, random_engine);

  // link original value to shuffled value correspondence
  std::map<T, T> originalToShuffledValues{};
  for(size_t i = 0; i < inputValues.size(); ++i) {
    originalToShuffledValues[inputValues[i]] = shuffledValues[i];
  }

// write shuffled values inside the output scalar field
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < n_blocks; ++i) {
    vtkDataSet *block = vtkDataSet::SafeDownCast(input->GetBlock(i));
    vtkDataArray *inputScalarField = this->GetInputArrayToProcess(0, block);
    int nValues = inputScalarField->GetNumberOfTuples();

    vtkSmartPointer<vtkDataArray> const outputArray
      = vtkSmartPointer<vtkDataArray>::Take(inputScalarField->NewInstance());
    outputArray->SetName(inputScalarField->GetName());
    outputArray->SetNumberOfComponents(1);
    outputArray->SetNumberOfTuples(inputScalarField->GetNumberOfTuples());
    T *const inputArrayPtr
      = static_cast<T *>(ttkUtils::GetVoidPointer(inputScalarField));
    T *const outputArrayPtr
      = static_cast<T *>(ttkUtils::GetVoidPointer(outputArray));
    for(int j = 0; j < nValues; ++j) {
      T newValue = originalToShuffledValues[inputArrayPtr[j]];
      outputArrayPtr[j] = newValue;
    }

    vtkDataSet *outputBlock = vtkDataSet::SafeDownCast(output->GetBlock(i));
    if(block->GetPointData()->GetArray(inputScalarField->GetName())) {
      outputBlock->GetPointData()->AddArray(outputArray);
    } else {
      outputBlock->GetCellData()->AddArray(outputArray);
    }
  }
  TTK_FORCE_USE(nThreads);
  return 1;
}

int ttkIdentifierRandomizer::RequestData(vtkInformation *ttkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  ttk::Timer t;

  vtkDataSet *input
    = vtkDataSet::SafeDownCast(vtkDataSet::GetData(inputVector[0]));
  if(input) {
    vtkDataSet *output = vtkDataSet::GetData(outputVector);
    output->ShallowCopy(input);
    // use a pointer-base copy for the input data -- to adapt if your wrapper
    // does not produce an output of the type of the input.

    // in the following, the target scalar field of the input is replaced in the
    // variable 'output' with the result of the computation.
    // if your wrapper produces an output of the same type of the input, you
    // should proceed in the same way.
    vtkDataArray *inputScalarField
      = this->GetInputArrayToProcess(0, inputVector);

    if(!inputScalarField) {
      printErr("Could not retrieve mandatory input array :(");
      return 0;
    }

    bool isPointData = false;
    if(input->GetPointData()->GetArray(inputScalarField->GetName())
       == inputScalarField) {
      isPointData = true;
    }

    this->printMsg("Shuffling " + std::string{isPointData ? "vertex" : "cell"}
                   + " field `" + std::string{inputScalarField->GetName()}
                   + "'...");

    // allocate the memory for the output scalar field
    vtkSmartPointer<vtkDataArray> const outputArray
      = vtkSmartPointer<vtkDataArray>::Take(inputScalarField->NewInstance());
    outputArray->SetName(inputScalarField->GetName());
    outputArray->SetNumberOfComponents(1);
    outputArray->SetNumberOfTuples(inputScalarField->GetNumberOfTuples());

    switch(outputArray->GetDataType()) {
      vtkTemplateMacro(shuffleScalarFieldValues(
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalarField)),
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputArray)),
        outputArray->GetNumberOfTuples(), this->RandomSeed, this->CompactRange,
        this->threadNumber_));
    }

    if(isPointData)
      output->GetPointData()->AddArray(outputArray);
    else
      output->GetCellData()->AddArray(outputArray);

    printMsg("Processed " + std::to_string(outputArray->GetNumberOfTuples())
               + (isPointData ? " vertices." : " cells."),
             1, t.getElapsedTime(), 1);

    printMsg(ttk::debug::Separator::L1);

  } else {
    vtkMultiBlockDataSet *input_mb = vtkMultiBlockDataSet::SafeDownCast(
      vtkMultiBlockDataSet::GetData(inputVector[0]));
    vtkMultiBlockDataSet *output_mb
      = vtkMultiBlockDataSet::GetData(outputVector);
    output_mb->ShallowCopy(input_mb);
    if(!input_mb) {
      printMsg("Invalid input.");
    }
    int n_blocks = input_mb->GetNumberOfBlocks();
    int currentType = 0;
    // check if the multiblock input and the input scalar field are valid and
    // retrieve the data type
    for(int i = 0; i < n_blocks; i++) {
      vtkDataSet *block = vtkDataSet::SafeDownCast(input_mb->GetBlock(i));
      if(!block)
        printMsg("Block " + std::to_string(i) + " invalid.");
      vtkDataArray *inputScalarField = this->GetInputArrayToProcess(0, block);

      if(!inputScalarField) {
        printWrn(
          "Block " + std::to_string(i)
          + " does not have the required input scalar field as data array.");
        continue;
      }

      if(i == 0) {
        currentType = inputScalarField->GetDataType();
        continue;
      } else {
        if(currentType != inputScalarField->GetDataType()) {
          printErr("All block's input scalar field must have the same type.");
          return 0;
        }
      }
    }
    switch(currentType) {
      vtkTemplateMacro(shuffleScalarFieldValuesMultiBlock<VTK_TT>(
        input_mb, output_mb, this->threadNumber_));
    }
  }

  return 1;
}
