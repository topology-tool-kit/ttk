#include <ttkIdentifiers.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIdentifiers);

ttkIdentifiers::ttkIdentifiers() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  setDebugMsgPrefix("Identifiers");

  vtkWarningMacro("`TTK Identifiers' is now deprecated. Please use "
                  "`Generate Ids' instead.");
}

ttkIdentifiers::~ttkIdentifiers() {
}

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

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifiers
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vtkSmartPointer<ttkSimplexIdTypeArray> cellIdentifiers
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();

  vertexIdentifiers->SetName(VertexFieldName.data());
  vertexIdentifiers->SetNumberOfComponents(1);
  vertexIdentifiers->SetNumberOfTuples(input->GetNumberOfPoints());

  cellIdentifiers->SetName(CellFieldName.data());
  cellIdentifiers->SetNumberOfComponents(1);
  cellIdentifiers->SetNumberOfTuples(input->GetNumberOfCells());

  SimplexId vertexNumber = input->GetNumberOfPoints();
  SimplexId cellNumber = input->GetNumberOfCells();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < vertexNumber; i++) {
    // avoid any processing if the abort signal is sent
    vertexIdentifiers->SetTuple1(i, i);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < cellNumber; i++) {
    // avoid any processing if the abort signal is sent
    cellIdentifiers->SetTuple1(i, i);
  }

  output->GetPointData()->AddArray(vertexIdentifiers);
  output->GetCellData()->AddArray(cellIdentifiers);

  printMsg("Processed " + std::to_string(vertexNumber) + " vertices and "
             + std::to_string(cellNumber) + " cells",
           1, t.getElapsedTime(), threadNumber_);

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
