#include <ttkIdentifiers.h>

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIdentifiers);

ttkIdentifiers::ttkIdentifiers() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkIdentifiers::~ttkIdentifiers() = default;

int ttkIdentifiers::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkIdentifiers::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkIdentifiers::RequestData(vtkInformation *ttkNotUsed(request),
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  Timer t;

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    this->printErr("Triangulation is NULL");
    return 0;
  }
  output->ShallowCopy(input);

#ifndef TTK_ENABLE_MPI
  // In case TTK is not compiled with MPI, the user can still create Global Ids
  int vertexNumber = input->GetNumberOfPoints();
  int cellNumber = input->GetNumberOfCells();

  vtkSmartPointer<vtkIdTypeArray> vtkVertexIdentifiers
    = vtkSmartPointer<vtkIdTypeArray>::New();
  vtkSmartPointer<vtkIdTypeArray> vtkCellIdentifiers
    = vtkSmartPointer<vtkIdTypeArray>::New();
  vtkVertexIdentifiers->SetName("GlobalPointIds");
  vtkVertexIdentifiers->SetNumberOfComponents(1);
  vtkVertexIdentifiers->SetNumberOfTuples(input->GetNumberOfPoints());

  vtkCellIdentifiers->SetName("GlobalCellIds");
  vtkCellIdentifiers->SetNumberOfComponents(1);
  vtkCellIdentifiers->SetNumberOfTuples(input->GetNumberOfCells());

#pragma omp parallel for
  for(ttk::SimplexId i = 0; i < vertexNumber; i++) {
    vtkVertexIdentifiers->SetTuple1(i, i);
  }

#pragma omp parallel for
  for(ttk::SimplexId i = 0; i < cellNumber; i++) {
    vtkCellIdentifiers->SetTuple1(i, i);
  }

  output->GetPointData()->AddArray(vtkVertexIdentifiers);
  output->GetCellData()->AddArray(vtkCellIdentifiers);
#endif

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
