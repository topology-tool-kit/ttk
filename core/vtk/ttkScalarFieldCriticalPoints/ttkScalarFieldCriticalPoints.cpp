#include <ttkScalarFieldCriticalPoints.h>

#include <vtkInformation.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkScalarFieldCriticalPoints);

ttkScalarFieldCriticalPoints::ttkScalarFieldCriticalPoints() {

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkScalarFieldCriticalPoints::~ttkScalarFieldCriticalPoints() = default;

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
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  else
    return 0;

  return 1;
}

int ttkScalarFieldCriticalPoints::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector, 0);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);

  int const keepGoing = checkEmptyMPIInput<Triangulation>(triangulation);
  if(keepGoing < 2) {
    return keepGoing;
  }

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

#ifdef TTK_ENABLE_MPI_TIME
  ttk::Timer t_mpi;
  ttk::startMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
#endif

  printMsg("Starting computation...");
  printMsg({{"  Scalar Array", inputScalarField->GetName()},
            {"  Offset Array", offsetField ? offsetField->GetName() : "None"}});

  int status = 0;
  ttkTemplateMacro(
    triangulation->getType(),
    (status = this->execute(
       static_cast<SimplexId *>(ttkUtils::GetVoidPointer(offsetField)),
       (TTK_TT *)triangulation->getData())));

  if(status < 0)
    return 0;

#ifdef TTK_ENABLE_MPI_TIME
  double elapsedTime = ttk::endMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
  if(ttk::MPIrank_ == 0) {
    printMsg("Computation performed using " + std::to_string(ttk::MPIsize_)
             + " MPI processes lasted :" + std::to_string(elapsedTime));
  }
#endif

  // allocate the output
  vtkNew<vtkSignedCharArray> vertexTypes{};
  vertexTypes->SetNumberOfComponents(1);
  vertexTypes->SetNumberOfTuples(criticalPoints_.size());
  vertexTypes->SetName("CriticalType");

  vtkNew<vtkPoints> pointSet{};
  pointSet->SetNumberOfPoints(criticalPoints_.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < criticalPoints_.size(); i++) {
    std::array<double, 3> p{};
    input->GetPoint(criticalPoints_[i].first, p.data());
    pointSet->SetPoint(i, p.data());
    vertexTypes->SetTuple1(i, (float)criticalPoints_[i].second);
  }

  ttkUtils::CellVertexFromPoints(output, pointSet);
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
#ifdef TTK_ENABLE_MPI
      if(hasInitializedMPI()) {
        vertexIds->SetTuple1(
          i, triangulation->getVertexGlobalId(criticalPoints_[i].first));
      } else
#endif // TTK_ENABLE_MPI
      {
        vertexIds->SetTuple1(i, criticalPoints_[i].first);
      }
    }

    output->GetPointData()->AddArray(vertexIds);
  } else {
    output->GetPointData()->RemoveArray(ttk::VertexScalarFieldName);
  }

  if(VertexScalars) {
    for(SimplexId i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {

      vtkDataArray *scalarField = input->GetPointData()->GetArray(i);
      vtkSmartPointer<vtkDataArray> const scalarArray{
        scalarField->NewInstance()};

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
