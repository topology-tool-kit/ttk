#include <ttkMorphologicalOperators.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkMorphologicalOperators);

ttkMorphologicalOperators::ttkMorphologicalOperators() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkMorphologicalOperators::~ttkMorphologicalOperators() {
}

int ttkMorphologicalOperators::FillInputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkMorphologicalOperators::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkMorphologicalOperators::RequestData(vtkInformation *ttkNotUsed(request),
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector) {
  // get input and output
  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);
  if(!input || !output)
    return 0;

  // copy input to output
  output->ShallowCopy(input);

  // get input label array
  auto inputLabels = this->GetInputArrayToProcess(0, inputVector);
  if(!inputLabels)
    return 0;
  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input labels needs to be a point data array.");
    return 0;
  }
  if(inputLabels->GetNumberOfComponents() != 1) {
    this->printErr("Input labels needs to be a scalar array.");
    return 0;
  }

  // -------------------------------------------------------------------------
  // Replace Variables in PivotLabel (e.g. {time[2]})
  // -------------------------------------------------------------------------
  double pivotLabel = 0;
  {
    std::string temp, errorMsg;
    if(!ttkUtils::replaceVariables(
         this->GetPivotLabel(), input->GetFieldData(), temp, errorMsg)) {
      this->printErr(errorMsg);
      return 0;
    }

    std::vector<double> values;
    ttkUtils::stringListToDoubleVector(temp, values);

    if(values.size() < 1) {
      this->printErr("Unable to parse pivot label as double.");
      return 0;
    }
    pivotLabel = values[0];
  }

  // create output labels
  auto outputLabels
    = vtkSmartPointer<vtkDataArray>::Take(inputLabels->NewInstance());
  outputLabels->DeepCopy(inputLabels);
  output->GetPointData()->AddArray(outputLabels);

  // get triangulation
  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  // precondition triangulation
  this->preconditionTriangulation(triangulation);

  int status = 0;

  ttkVtkTemplateMacro(
    inputLabels->GetDataType(), triangulation->getType(),
    (status = this->execute<VTK_TT, TTK_TT>(
       // Output
       static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputLabels)),

       // Input
       this->Mode, this->Iterations, this->Grayscale,
       static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputLabels)), pivotLabel,
       static_cast<TTK_TT *>(triangulation->getData()))));

  if(!status)
    return 0;

  return 1;
}
