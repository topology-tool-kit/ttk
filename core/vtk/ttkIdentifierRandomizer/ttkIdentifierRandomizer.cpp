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
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
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

  // link original value to shuffled value correspondance
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

int ttkIdentifierRandomizer::RequestData(vtkInformation *ttkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  ttk::Timer t;

  bool isPointData = false;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  // in the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  vtkDataArray *inputScalarField = this->GetInputArrayToProcess(0, inputVector);

  if(!inputScalarField) {
    printErr("Could not retrieve mandatory input array :(");
    return 0;
  }

  if(input->GetPointData()->GetArray(inputScalarField->GetName())
     == inputScalarField) {
    isPointData = true;
  }

  this->printMsg("Shuffling " + std::string{isPointData ? "vertex" : "cell"}
                 + " field `" + std::string{inputScalarField->GetName()}
                 + "'...");

  // allocate the memory for the output scalar field
  vtkSmartPointer<vtkDataArray> outputArray
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

  return 1;
}
