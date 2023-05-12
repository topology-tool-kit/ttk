#include <ttkMacros.h>
#include <ttkPathCompression.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>
#include <vtkUnsignedCharArray.h>

vtkStandardNewMacro(ttkPathCompression);

ttkPathCompression::ttkPathCompression() {
  this->setDebugMsgPrefix("PathCompression");
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkPathCompression::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPathCompression::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <typename triangulationType>
int ttkPathCompression::dispatch(const SimplexId *const inputOrderArray,
                                 const triangulationType &triangulation) {

  const int ret = this->execute<triangulationType>(
    segmentations_, inputOrderArray, triangulation);

  if(ret != 0)
    return !this->printErr("PathCompression.execute() error");

  return ret;
}

int ttkPathCompression::RequestData(vtkInformation *ttkNotUsed(request),
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputMorseComplexes = vtkDataSet::GetData(outputVector, 0);

  if(!input)
    return !this->printErr("Input pointer is NULL.");

  if(input->GetNumberOfPoints() == 0)
    return !this->printErr("Input has no point.");

  if(!outputMorseComplexes)
    return !this->printErr("Output pointers are NULL.");

  const auto triangulation = ttkAlgorithm::GetTriangulation(input);

  if(triangulation == nullptr)
    return !this->printErr("Triangulation is null");

  this->preconditionTriangulation(triangulation);

  const auto inputScalars = this->GetInputArrayToProcess(0, inputVector);

  if(inputScalars == nullptr)
    return !this->printErr("No input scalars");

  auto inputOrderArray = ttkAlgorithm::GetOrderArray(
    input, 0, 1, this->ForceInputOffsetScalarField);

  if(inputOrderArray == nullptr)
    return !this->printErr("No order array");

  if(inputOrderArray->GetDataType() != VTK_INT
     && inputOrderArray->GetDataType() != VTK_ID_TYPE)
    return !this->printErr("input offset field type not supported.");

  this->printMsg("Launching computation on field `"
                 + std::string(inputScalars->GetName()) + "'...");

  const SimplexId numberOfVertices = triangulation->getNumberOfVertices();

  if(!numberOfVertices)
    return !this->printErr("Input has no vertices.");

  vtkNew<ttkSimplexIdTypeArray> ascendingSegmentation{};
  vtkNew<ttkSimplexIdTypeArray> descendingSegmentation{};
  vtkNew<ttkSimplexIdTypeArray> morseSmaleSegmentation{};

  if(!ascendingSegmentation || !descendingSegmentation
     || !morseSmaleSegmentation)
    return !this->printErr("Segmentation vtkDataArray allocation problem.");

  ascendingSegmentation->SetNumberOfComponents(1);
  ascendingSegmentation->SetNumberOfTuples(numberOfVertices);
  ascendingSegmentation->SetName("AscendingSegmentation");

  descendingSegmentation->SetNumberOfComponents(1);
  descendingSegmentation->SetNumberOfTuples(numberOfVertices);
  descendingSegmentation->SetName("DescendingSegmentation");

  morseSmaleSegmentation->SetNumberOfComponents(1);
  morseSmaleSegmentation->SetNumberOfTuples(numberOfVertices);
  morseSmaleSegmentation->SetName("MorseSmaleSegmentation");

  this->segmentations_
    = {ttkUtils::GetPointer<SimplexId>(ascendingSegmentation),
       ttkUtils::GetPointer<SimplexId>(descendingSegmentation),
       ttkUtils::GetPointer<SimplexId>(morseSmaleSegmentation)};

  int ret{};

  ttkTemplateMacro(
    triangulation->getType(),
    (ret = dispatch<TTK_TT>(ttkUtils::GetPointer<SimplexId>(inputOrderArray),
                            *static_cast<TTK_TT *>(triangulation->getData()))));

  if(ret != 0)
    return -1;

  outputMorseComplexes->ShallowCopy(input);

  if(ComputeAscendingSegmentation || ComputeDescendingSegmentation
     || ComputeMSSegmentationHash) {
    vtkPointData *pointData = outputMorseComplexes->GetPointData();

    if(!pointData)
      return !this->printErr("outputMorseComplexes has no point data.");

    if(ComputeDescendingSegmentation || ComputeMSSegmentationHash)
      pointData->AddArray(descendingSegmentation);
    if(ComputeAscendingSegmentation || ComputeMSSegmentationHash)
      pointData->AddArray(ascendingSegmentation);
    if(ComputeMSSegmentationHash)
      pointData->AddArray(morseSmaleSegmentation);
  }

  return 1;
}
