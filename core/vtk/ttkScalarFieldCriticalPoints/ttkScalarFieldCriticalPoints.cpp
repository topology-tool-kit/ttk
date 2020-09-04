#include <ttkScalarFieldCriticalPoints.h>

#include <vtkInformation.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkSignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkScalarFieldCriticalPoints);

ttkScalarFieldCriticalPoints::ttkScalarFieldCriticalPoints() {

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkScalarFieldCriticalPoints::~ttkScalarFieldCriticalPoints() {
}

int ttkScalarFieldCriticalPoints::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  else
    return 0;

  return 1;
}

int ttkScalarFieldCriticalPoints::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  else
    return 0;

  return 1;
}

int ttkScalarFieldCriticalPoints::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector, 0);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  if(VertexBoundary)
    triangulation->preconditionBoundaryVertices();

  // in the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.

  vtkDataArray *inputScalarField = this->GetInputArrayToProcess(0, inputVector);
  if(!inputScalarField)
    return 0;

  vtkDataArray *offsetField
    = this->GetOrderArray(input, 0, 1, ForceInputOffsetScalarField);

  // setting up the base layer
  this->preconditionTriangulation(triangulation);
  this->setOutput(&criticalPoints_);

  printMsg("Starting computation...");
  printMsg({{"  Scalar Array", inputScalarField->GetName()},
            {"  Offset Array", offsetField ? offsetField->GetName() : "None"}});

  int status = 0;
  if(offsetField->GetDataType() == VTK_INT) {
    ttkTemplateMacro(
      triangulation->getType(),
      (status = this->execute(
         static_cast<int *>(ttkUtils::GetVoidPointer(offsetField)),
         (TTK_TT *)triangulation->getData())));
  } else if(offsetField->GetDataType() == VTK_ID_TYPE) {
    ttkTemplateMacro(
      triangulation->getType(),
      (status = this->execute(
         static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(offsetField)),
         (TTK_TT *)triangulation->getData())));
  } else {
    this->printErr("Wrong offset field type");
    return 0;
  }
  if(status < 0)
    return 0;

  // allocate the output
  vtkNew<vtkSignedCharArray> vertexTypes{};
  vertexTypes->SetNumberOfComponents(1);
  vertexTypes->SetNumberOfTuples(criticalPoints_.size());
  vertexTypes->SetName("CriticalType");

  vtkNew<vtkPoints> pointSet{};
  pointSet->SetNumberOfPoints(criticalPoints_.size());
  double p[3];
  for(size_t i = 0; i < criticalPoints_.size(); i++) {
    input->GetPoint(criticalPoints_[i].first, p);
    pointSet->SetPoint(i, p);
    vertexTypes->SetTuple1(i, (float)criticalPoints_[i].second);
  }
  output->SetPoints(pointSet);
  output->GetPointData()->AddArray(vertexTypes);

  if(VertexBoundary) {
    vtkNew<vtkSignedCharArray> vertexBoundary{};
    vertexBoundary->SetNumberOfComponents(1);
    vertexBoundary->SetNumberOfTuples(criticalPoints_.size());
    vertexBoundary->SetName("IsOnBoundary");

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t i = 0; i < criticalPoints_.size(); i++) {
      vertexBoundary->SetTuple1(
        i, (signed char)triangulation->isVertexOnBoundary(
             criticalPoints_[i].first));
    }

    output->GetPointData()->AddArray(vertexBoundary);
  } else {
    output->GetPointData()->RemoveArray("IsOnBoundary");
  }

  if(VertexIds) {
    vtkNew<ttkSimplexIdTypeArray> vertexIds{};
    vertexIds->SetNumberOfComponents(1);
    vertexIds->SetNumberOfTuples(criticalPoints_.size());
    vertexIds->SetName(ttk::VertexScalarFieldName);

    for(size_t i = 0; i < criticalPoints_.size(); i++) {
      vertexIds->SetTuple1(i, criticalPoints_[i].first);
    }

    output->GetPointData()->AddArray(vertexIds);
  } else {
    output->GetPointData()->RemoveArray(ttk::VertexScalarFieldName);
  }

  if(VertexScalars) {
    for(SimplexId i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {

      vtkDataArray *scalarField = input->GetPointData()->GetArray(i);
      vtkSmartPointer<vtkDataArray> scalarArray{scalarField->NewInstance()};

      scalarArray->SetNumberOfComponents(scalarField->GetNumberOfComponents());
      scalarArray->SetNumberOfTuples(criticalPoints_.size());
      scalarArray->SetName(scalarField->GetName());
      std::vector<double> value(scalarField->GetNumberOfComponents());
      for(size_t j = 0; j < criticalPoints_.size(); j++) {
        scalarField->GetTuple(criticalPoints_[j].first, value.data());
        scalarArray->SetTuple(j, value.data());
      }
      output->GetPointData()->AddArray(scalarArray);
    }
  } else {
    for(SimplexId i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {
      output->GetPointData()->RemoveArray(
        input->GetPointData()->GetArray(i)->GetName());
    }
  }

  return 1;
}
