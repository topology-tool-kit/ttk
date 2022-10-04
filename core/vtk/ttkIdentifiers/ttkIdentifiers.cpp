#include <ttkIdentifiers.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <map>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIdentifiers);

ttkIdentifiers::ttkIdentifiers() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  setDebugMsgPrefix("Identifiers");

  //   vtkWarningMacro("`TTK Identifiers' is now deprecated. Please use "
  //                   "`Generate Ids' instead.");
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

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;
  this->preconditionTriangulation(triangulation);
  // vtkImageData* data = vtkImageData::SafeDownCast(input);
  // data->ComputeBounds();
  // double* bounds = data->GetBounds();
  //   cout << "bounds: " << bounds[0] << "," << bounds[1] << "  " << bounds[2]
  //        << "," << bounds[3] << " " << bounds[4] << "," <<
  //        bounds[5] << endl;
  // double* point = data->GetPoint(0);
  // printMsg("Point: "+std::to_string(point[0])+" "+std::to_string(point[1])+"
  // "+std::to_string(point[2])); data->SetExtent(wholeExtent); data->Update();

  Timer t;

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  std::vector<ttk::SimplexId> vertexIdentifiers(input->GetNumberOfPoints(), -1);
  this->setVertexIdentifiers(&vertexIdentifiers);

  std::vector<ttk::SimplexId> cellIdentifiers(input->GetNumberOfCells(), -1);
  this->setCellIdentifiers(&cellIdentifiers);

  this->setVertexNumber(input->GetNumberOfPoints());
  this->setCellNumber(input->GetNumberOfCells());

  this->setVertRankArray(ttkUtils::GetPointer<ttk::SimplexId>(
    input->GetPointData()->GetArray("RankArray")));
  this->setCellRankArray(ttkUtils::GetPointer<ttk::SimplexId>(
    input->GetCellData()->GetArray("RankArray")));

  this->setVertGhost(ttkUtils::GetPointer<unsigned char>(
    input->GetPointData()->GetArray("vtkGhostType")));
  this->setCellGhost(ttkUtils::GetPointer<unsigned char>(
    input->GetCellData()->GetArray("vtkGhostType")));

  double *boundingBox = input->GetBounds();
  initializeNeighbors(boundingBox);

  this->initializeMPITypes();

  vtkIdList *points = vtkIdList::New();
  input->GetCellPoints(0, points);
  this->setNbPoints(points->GetNumberOfIds());
  this->setBounds(boundingBox);

  ttkTemplateMacro(triangulation->getType(),
                   this->execute((TTK_TT *)triangulation->getData()));

  // #ifdef TTK_ENABLE_OPENMP
  // #pragma omp parallel for num_threads(threadNumber_)
  // #endif
  //   for(SimplexId i = 0; i < vertexNumber; i++) {
  //     // avoid any processing if the abort signal is sent
  //     vertexIdentifiers->SetTuple1(i, i);
  //   }

  // #ifdef TTK_ENABLE_OPENMP
  // #pragma omp parallel for num_threads(threadNumber_)
  // #endif
  //   for(SimplexId i = 0; i < cellNumber; i++) {
  //     // avoid any processing if the abort signal is sent
  //     cellIdentifiers->SetTuple1(i, i);
  //   }
  vtkSmartPointer<ttkSimplexIdTypeArray> vtkVertexIdentifiers
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vtkSmartPointer<ttkSimplexIdTypeArray> vtkCellIdentifiers
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vtkVertexIdentifiers->SetName(VertexFieldName.data());
  vtkVertexIdentifiers->SetNumberOfComponents(1);
  vtkVertexIdentifiers->SetNumberOfTuples(input->GetNumberOfPoints());

  vtkCellIdentifiers->SetName(CellFieldName.data());
  vtkCellIdentifiers->SetNumberOfComponents(1);
  vtkCellIdentifiers->SetNumberOfTuples(input->GetNumberOfCells());

  for(SimplexId i = 0; i < vertexNumber_; i++) {
    vtkVertexIdentifiers->SetTuple1(i, vertexIdentifiers_->at(i));
  }

  for(SimplexId i = 0; i < cellNumber_; i++) {
    vtkCellIdentifiers->SetTuple1(i, cellIdentifiers_->at(i));
  }

  output->GetPointData()->AddArray(vtkVertexIdentifiers);
  output->GetCellData()->AddArray(vtkCellIdentifiers);

  printMsg("Processed " + std::to_string(vertexNumber_) + " vertices and "
             + std::to_string(cellNumber_) + " cells",
           1, t.getElapsedTime(), threadNumber_);

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
