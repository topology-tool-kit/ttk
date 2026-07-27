#include <ttkIdentifiers.h>

#include <Triangulation.h>
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

  Timer const t;

  int const keepGoing = checkEmptyMPIInput<vtkDataSet>(input);
  if(keepGoing < 2) {
    return keepGoing;
  }

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    this->printErr("Triangulation is NULL");
    return 0;
  }

  // The following is reserved to vtkImageData. For UnstructuredGrid and
  // vtkPolyData, the global identifiers are always stored as vtkDataArrays
  // automatically during preconditioning

  ttk::SimplexId const numberOfVertices = triangulation->getNumberOfVertices();
  vtkNew<vtkIdTypeArray> globalPointIds;
  globalPointIds->SetNumberOfTuples(numberOfVertices);
  globalPointIds->SetNumberOfComponents(1);
  globalPointIds->SetName("GlobalPointIds");
  for(int i = 0; i < numberOfVertices; i++) {
#ifdef TTK_ENABLE_MPI
    if(input->GetDataObjectType() == VTK_IMAGE_DATA) {
      globalPointIds->SetTuple1(i, triangulation->getVertexGlobalId(i));
    }
#else
    globalPointIds->SetTuple1(i, i);
#endif
  }

#ifdef TTK_ENABLE_MPI
  if(input->GetDataObjectType() == VTK_IMAGE_DATA) {
#endif
    input->GetPointData()->AddArray(globalPointIds);
#ifdef TTK_ENABLE_MPI
  }
#endif

  output->ShallowCopy(input);

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
