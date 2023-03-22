#include <ttkClusteringMetrics.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkVariantArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkClusteringMetrics);

ttkClusteringMetrics::ttkClusteringMetrics() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int ttkClusteringMetrics::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkClusteringMetrics::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkClusteringMetrics::RequestData(vtkInformation *ttkNotUsed(request),
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  // Get input object from input vector
  vtkTable *input1 = vtkTable::GetData(inputVector[0]);
  vtkTable *input2 = vtkTable::GetData(inputVector[1]);
  vtkTable *output = vtkTable::GetData(outputVector);
  if(!input1 || !input2 || !output)
    return 0;
  // We do not use the GetInputArrayToProcess but a version of it allowing us to
  // retrieve an array of int instead of an array of double.
  vtkDataArray *inputClustering1 = this->GetInputArrayToProcess(0, input1);
  vtkDataArray *inputClustering2 = this->GetInputArrayToProcess(1, input2);

  if(!inputClustering1) {
    this->printErr("Unable to retrieve input array for first clustering.");
    return 1;
  }
  if(!inputClustering2) {
    this->printErr("Unable to retrieve input array for second clustering.");
    return 1;
  }

  vtkIntArray *intArray1 = vtkIntArray::SafeDownCast(inputClustering1);
  vtkIntArray *intArray2 = vtkIntArray::SafeDownCast(inputClustering2);

  // To avoid making a copy of the data.
  const int *values1 = ttkUtils::GetPointer<int>(inputClustering1);
  const int *values2 = ttkUtils::GetPointer<int>(inputClustering2);

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");
  this->printMsg("  First clustering column: "
                 + std::string(intArray1->GetName()));
  this->printMsg("  Second clustering column: "
                 + std::string(intArray2->GetName()));

  size_t nbVal1 = inputClustering1->GetNumberOfTuples();
  size_t nbVal2 = inputClustering2->GetNumberOfTuples();

  if(nbVal1 != nbVal2) {
    this->printMsg("Error : the two clusterings must have the same size\n");
    return 0;
  }
  size_t nbVal = nbVal1;

  double nmiValue = 0, ariValue = 0;
  this->execute(values1, values2, nbVal, nmiValue, ariValue);
  vtkNew<vtkDoubleArray> nmiValArray{}, ariValArray{};
  output->SetNumberOfRows(1);

  nmiValArray->SetName("NMIValue");
  nmiValArray->SetNumberOfTuples(1);
  nmiValArray->SetTuple1(0, nmiValue);
  output->AddColumn(nmiValArray);

  ariValArray->SetName("ARIValue");
  ariValArray->SetNumberOfTuples(1);
  ariValArray->SetTuple1(0, ariValue);
  output->AddColumn(ariValArray);

  // return success
  return 1;
}
