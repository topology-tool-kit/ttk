#include <ttkGhostCellPreconditioning.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkGhostCellPreconditioning);

ttkGhostCellPreconditioning::ttkGhostCellPreconditioning() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->setDebugMsgPrefix("GhostCellPreconditioning");

  // ensure that modifying the selection re-triggers the filter
  // (c.f. vtkPassSelectedArrays.cxx)
  this->ArraySelection->AddObserver(
    vtkCommand::ModifiedEvent, this, &ttkGhostCellPreconditioning::Modified);
}

int ttkGhostCellPreconditioning::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkGhostCellPreconditioning::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkGhostCellPreconditioning::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);
  ttk::Timer tm{};

  if(input == nullptr || output == nullptr) {
    return 0;
  }
  output->ShallowCopy(input);

  auto pointData = input->GetPointData();
  auto cellData = input->GetCellData();
  ttk::SimplexId nVertices = input->GetNumberOfPoints();
  ttk::SimplexId nCells = input->GetNumberOfCells();

  this->printMsg("#Points: " + std::to_string(nVertices));
  this->printMsg("#Cells: " + std::to_string(nCells));

  auto *verticesGlobalIds
    = ttkUtils::GetPointer<ttk::LongSimplexId>(pointData->GetGlobalIds());
  auto *verticesGhostCells
    = ttkUtils::GetPointer<unsigned char>(pointData->GetArray("vtkGhostType"));
  auto *cellsGlobalIds
    = ttkUtils::GetPointer<ttk::LongSimplexId>(cellData->GetGlobalIds());
  auto *cellsGhostCells
    = ttkUtils::GetPointer<unsigned char>(cellData->GetArray("vtkGhostType"));

  if(verticesGlobalIds != nullptr && verticesGhostCells != nullptr
     && cellsGlobalIds != nullptr && cellsGhostCells != nullptr) {
#ifdef TTK_ENABLE_MPI
    if(ttk::hasInitializedMPI()) {
      if(ttk::MPIrank_ == 0)
        this->printMsg(
          "Global Point Ids and Ghost Cells exist, therefore we can continue!");
      this->printMsg("#Ranks " + std::to_string(ttk::MPIsize_)
                     + ", this is rank " + std::to_string(ttk::MPIrank_));
      std::vector<int> verticesRankArray(nVertices, 0);
      std::vector<int> cellsRankArray(nCells, 0);
      double *boundingBox = input->GetBounds();
      std::vector<int> neighbors{};
      ttk::preconditionNeighborsUsingBoundingBox(boundingBox, neighbors);
      ttk::produceRankArray(verticesRankArray, verticesGlobalIds,
                            verticesGhostCells, nVertices, boundingBox,
                            neighbors);
      ttk::produceRankArray(cellsRankArray, cellsGlobalIds, cellsGhostCells,
                            nCells, boundingBox, neighbors);

      vtkNew<vtkIntArray> vtkVerticesRankArray{};
      vtkVerticesRankArray->SetName("RankArray");
      vtkVerticesRankArray->SetNumberOfComponents(1);
      vtkVerticesRankArray->SetNumberOfTuples(nVertices);

      vtkNew<vtkIntArray> vtkCellsRankArray{};
      vtkCellsRankArray->SetName("RankArray");
      vtkCellsRankArray->SetNumberOfComponents(1);
      vtkCellsRankArray->SetNumberOfTuples(nCells);

      for(int i = 0; i < nVertices; i++) {
        vtkVerticesRankArray->SetComponent(i, 0, verticesRankArray[i]);
      }
      for(int i = 0; i < nCells; i++) {
        vtkCellsRankArray->SetComponent(i, 0, cellsRankArray[i]);
      }

      output->GetPointData()->AddArray(vtkVerticesRankArray);
      output->GetCellData()->AddArray(vtkCellsRankArray);

      this->printMsg("Preprocessed RankArray", 1.0, tm.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    } else {
      this->printMsg("Necessary arrays are present,  TTK is built with MPI "
                     "support, but not run with mpirun. Running sequentially.");
      return 0;
    }
#else
    this->printMsg(
      "Necessary arrays are present, but TTK is not built with MPI support");
    return 0;

#endif
  } else {
    this->printMsg("Necessary arrays are not present.");
    return 0;
  }
}
