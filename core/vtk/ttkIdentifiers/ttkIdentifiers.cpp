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
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

vtkCellArray *GetCells(vtkDataSet *dataSet) {
  switch(dataSet->GetDataObjectType()) {
    case VTK_UNSTRUCTURED_GRID: {
      auto dataSetAsUG = static_cast<vtkUnstructuredGrid *>(dataSet);
      return dataSetAsUG->GetCells();
    }
    case VTK_POLY_DATA: {
      auto dataSetAsPD = static_cast<vtkPolyData *>(dataSet);
      return dataSetAsPD->GetNumberOfPolys() > 0   ? dataSetAsPD->GetPolys()
             : dataSetAsPD->GetNumberOfLines() > 0 ? dataSetAsPD->GetLines()
                                                   : dataSetAsPD->GetVerts();
    }
  }
  return nullptr;
}

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

  std::vector<ttk::SimplexId> vertexIdentifiers(input->GetNumberOfPoints(), -1);
  this->setVertexIdentifiers(&vertexIdentifiers);

  std::vector<ttk::SimplexId> cellIdentifiers(input->GetNumberOfCells(), -1);
  this->setCellIdentifiers(&cellIdentifiers);

  this->setVertexNumber(input->GetNumberOfPoints());
  this->setCellNumber(input->GetNumberOfCells());

  Timer t;

#ifdef TTK_ENABLE_MPI
  double *boundingBox = input->GetBounds();
  this->setBounds(boundingBox);

  switch(input->GetDataObjectType()) {
    case VTK_UNSTRUCTURED_GRID:
    case VTK_POLY_DATA: {

      this->setVertRankArray(ttkUtils::GetPointer<ttk::SimplexId>(
        input->GetPointData()->GetArray("RankArray")));
      this->setCellRankArray(ttkUtils::GetPointer<ttk::SimplexId>(
        input->GetCellData()->GetArray("RankArray")));

      this->setVertGhost(ttkUtils::GetPointer<unsigned char>(
        input->GetPointData()->GetArray("vtkGhostType")));
      this->setCellGhost(ttkUtils::GetPointer<unsigned char>(
        input->GetCellData()->GetArray("vtkGhostType")));
      vtkPointSet *pointSet = vtkPointSet::GetData(inputVector[0]);
      setPointSet(
        static_cast<float *>(ttkUtils::GetVoidPointer(pointSet->GetPoints())));
      vtkCellArray *cells = GetCells(pointSet);

      if(!cells->IsStorage64Bit()) {
        if(cells->CanConvertTo64BitStorage()) {
          this->printWrn("Converting the cell array to 64-bit storage");
          bool success = cells->ConvertTo64BitStorage();
          if(!success) {
            this->printErr(
              "Error converting the provided cell array to 64-bit storage");
            return {};
          }
        } else {
          this->printErr(
            "Cannot convert the provided cell array to 64-bit storage");
          return {};
        }
      }

      setConnectivity(static_cast<ttk::LongSimplexId *>(
        ttkUtils::GetVoidPointer(cells->GetConnectivityArray())));

      std::vector<std::vector<ttk::SimplexId>> pointsToCells(vertexNumber_);
      vtkIdList *cellList = vtkIdList::New();
      for(ttk::SimplexId i = 0; i < vertexNumber_; i++) {
        input->GetPointCells(i, cellList);
        for(int j = 0; j < cellList->GetNumberOfIds(); j++) {
          if(cellGhost_[cellList->GetId(j)] == 0) {
            pointsToCells[i].push_back(cellList->GetId(j));
          }
        }
      }

      this->setPointsToCells(pointsToCells);

      initializeNeighbors(boundingBox);

      this->initializeMPITypes();

      vtkIdList *pointCell = vtkIdList::New();
      input->GetCellPoints(0, pointCell);
      this->setNbPoints(pointCell->GetNumberOfIds());
      setDomainDimension(nbPoints_ - 1);
      this->executePolyData();
      break;
    }
    case VTK_IMAGE_DATA: {
      vtkImageData *data = vtkImageData::SafeDownCast(input);
      setDims(data->GetDimensions());
      setSpacing(data->GetSpacing());

      this->executeImageData();

      break;
    }
    default: {
      this->printErr("Unable to triangulate `"
                     + std::string(input->GetClassName()) + "`");
    }
  }

#else
  this->execute();
#endif
  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

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
