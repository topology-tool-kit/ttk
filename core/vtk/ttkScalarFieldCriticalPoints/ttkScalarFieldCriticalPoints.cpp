#include <ttkScalarFieldCriticalPoints.h>

#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
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

int ttkScalarFieldCriticalPoints::RequestData(vtkInformation 
*request,vtkInformationVector **inputVector,vtkInformationVector *outputVector) 
{
  
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *output = 
    vtkUnstructuredGrid::GetData(outputVector, 0);
  
  Timer t;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    printErr("Input pointer is null :(");
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    printErr("Input has no points :(");
    return -1;
  }
#endif

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    printErr("Input triangulation is nullptr :(");
    return -1;
  }
#endif

  if(VertexBoundary)
    triangulation->preconditionBoundaryVertices();

  // in the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  

  vtkDataArray *inputScalarField = 
    this->GetInputArrayToProcess(0, inputVector);
  if(!inputScalarField) return 0;

  vtkDataArray *offsetField = nullptr;

  if(this->GetInputArrayInformation(1))
    offsetField = this->GetInputArrayToProcess(1, inputVector);

  if((!offsetField) || (!ForceInputOffsetScalarField)) {
    offsetField = input->GetPointData()->GetArray(ttk::OffsetScalarFieldName);
  }

  sosOffsets_.resize(inputScalarField->GetNumberOfTuples());
  for(SimplexId i = 0; i < inputScalarField->GetNumberOfTuples(); i++) {
    SimplexId offset = i;
    if(offsetField){
      offset = offsetField->GetTuple1(i);
    }
    sosOffsets_[i] = offset;
  }
    
  // setting up the base layer
  this->setupTriangulation(triangulation);
  this->setSosOffsets(&sosOffsets_);
  this->setOutput(&criticalPoints_);
    
  printMsg("Starting computation on array `" 
    + string(inputScalarField->GetName()) + "'...");
  if(offsetField) {
    printMsg("  Using offset array `" + string(offsetField->GetName())
             + "'...");
  }

  int status = 0;
  ttkVtkTemplateMacro(
    triangulation->getType(),
    inputScalarField->GetDataType(),
    (status = this->execute<VTK_TT, TTK_TT>(
      (TTK_TT *) triangulation->getData(),
      (VTK_TT *) ttkUtils::GetVoidPointer(inputScalarField)
    ))
  );
  if(status < 0)
    return 0;

  // allocate the output
  vtkSmartPointer<vtkCharArray> vertexTypes
    = vtkSmartPointer<vtkCharArray>::New();

  vertexTypes->SetNumberOfComponents(1);
  vertexTypes->SetNumberOfTuples(criticalPoints_.size());
  vertexTypes->SetName("CriticalType");

  vtkSmartPointer<vtkPoints> pointSet = vtkSmartPointer<vtkPoints>::New();
  pointSet->SetNumberOfPoints(criticalPoints_.size());
  double p[3];
  for(SimplexId i = 0; i < (SimplexId)criticalPoints_.size(); i++) {
    input->GetPoint(criticalPoints_[i].first, p);
    pointSet->SetPoint(i, p);
    vertexTypes->SetTuple1(i, (float)criticalPoints_[i].second);
  }
  output->SetPoints(pointSet);
  output->GetPointData()->AddArray(vertexTypes);

  if(VertexBoundary) {
    vtkSmartPointer<vtkCharArray> vertexBoundary
      = vtkSmartPointer<vtkCharArray>::New();
    vertexBoundary->SetNumberOfComponents(1);
    vertexBoundary->SetNumberOfTuples(criticalPoints_.size());
    vertexBoundary->SetName("IsOnBoundary");

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)criticalPoints_.size(); i++) {
      vertexBoundary->SetTuple1(
        i, (char)triangulation->isVertexOnBoundary(criticalPoints_[i].first));
    }

    output->GetPointData()->AddArray(vertexBoundary);
  } else {
    output->GetPointData()->RemoveArray("IsOnBoundary");
  }

  if(VertexIds) {
    vtkSmartPointer<ttkSimplexIdTypeArray> vertexIds
      = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
    vertexIds->SetNumberOfComponents(1);
    vertexIds->SetNumberOfTuples(criticalPoints_.size());
    vertexIds->SetName(ttk::VertexScalarFieldName);

    for(SimplexId i = 0; i < (SimplexId)criticalPoints_.size(); i++) {
      vertexIds->SetTuple1(i, criticalPoints_[i].first);
    }

    output->GetPointData()->AddArray(vertexIds);
  } else {
    output->GetPointData()->RemoveArray(ttk::VertexScalarFieldName);
  }

  if(VertexScalars) {
    for(SimplexId i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {

      vtkDataArray *scalarField = input->GetPointData()->GetArray(i);
      vtkSmartPointer<vtkDataArray> scalarArray;

      auto copyToScalarArray = [&]() {
        scalarArray->SetNumberOfComponents(
          scalarField->GetNumberOfComponents());
        scalarArray->SetNumberOfTuples(criticalPoints_.size());
        scalarArray->SetName(scalarField->GetName());
        std::vector<double> value(scalarField->GetNumberOfComponents());
        for(SimplexId j = 0; j < (SimplexId)criticalPoints_.size(); j++) {
          scalarField->GetTuple(criticalPoints_[j].first, value.data());
          scalarArray->SetTuple(j, value.data());
        }
        output->GetPointData()->AddArray(scalarArray);
      };

      switch(scalarField->GetDataType()) {
        case VTK_CHAR:
          scalarArray = vtkSmartPointer<vtkCharArray>::New();
          copyToScalarArray();
          break;
        case VTK_DOUBLE:
          scalarArray = vtkSmartPointer<vtkDoubleArray>::New();
          copyToScalarArray();
          break;
        case VTK_FLOAT:
          scalarArray = vtkSmartPointer<vtkFloatArray>::New();
          copyToScalarArray();
          break;
        case VTK_INT:
          scalarArray = vtkSmartPointer<vtkIntArray>::New();
          copyToScalarArray();
          break;
        case VTK_ID_TYPE:
          scalarArray = vtkSmartPointer<vtkIdTypeArray>::New();
          copyToScalarArray();
          break;
        default: {
          stringstream msg;
          msg << "[ttkScalarFieldCriticalPoints] Scalar attachment: "
              << "unsupported data type :(" << endl;
          dMsg(cerr, msg.str(), detailedInfoMsg);
        } break;
      }
    }
  } else {
    for(SimplexId i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {
      output->GetPointData()->RemoveArray(
        input->GetPointData()->GetArray(i)->GetName());
    }
  }

  return 1;
}
