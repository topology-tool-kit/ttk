#include <ttkMacros.h>
#include <ttkPersistenceCurve.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkUtils.h>

#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkPersistenceCurve);

ttkPersistenceCurve::ttkPersistenceCurve() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(4);
}

vtkTable *ttkPersistenceCurve::GetOutput() {
  return this->GetOutput(0);
}

vtkTable *ttkPersistenceCurve::GetOutput(int port) {
  return vtkTable::SafeDownCast(this->GetOutputDataObject(port));
}

int ttkPersistenceCurve::FillInputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkPersistenceCurve::FillOutputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2 || port == 3) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

static void getPersistenceCurve(vtkTable *outputCurve,
                                const ttk::PersistenceCurve::PlotType &plot,
                                const std::string &pairsType) {

  if(plot.empty()) {
    return;
  }

  vtkNew<vtkDoubleArray> persistenceScalars{};
  const auto persQualifier = "Persistence (" + pairsType + ")";
  persistenceScalars->SetName(persQualifier.c_str());
  persistenceScalars->SetNumberOfTuples(plot.size());

  vtkNew<ttkSimplexIdTypeArray> numberOfPairsScalars{};
  const auto nPairsQualifier = "Number Of Pairs (" + pairsType + ")";
  numberOfPairsScalars->SetName(nPairsQualifier.c_str());
  numberOfPairsScalars->SetNumberOfTuples(plot.size());

  for(size_t i = 0; i < plot.size(); ++i) {
    persistenceScalars->SetTuple1(i, plot[i].first);
    numberOfPairsScalars->SetTuple1(i, plot[i].second);
  }

  vtkNew<vtkTable> persistenceCurve{};
  persistenceCurve->AddColumn(persistenceScalars);
  persistenceCurve->AddColumn(numberOfPairsScalars);

  outputCurve->ShallowCopy(persistenceCurve);
}

int ttkPersistenceCurve::RequestData(vtkInformation *ttkNotUsed(request),
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {

  ttk::Timer timer;

  auto *input = vtkUnstructuredGrid::GetData(inputVector[0]);

  auto *minSadTable = vtkTable::GetData(outputVector, 0);
  auto *sadSadTable = vtkTable::GetData(outputVector, 1);
  auto *sadMaxTable = vtkTable::GetData(outputVector, 2);
  auto *allPairsTable = vtkTable::GetData(outputVector, 3);

  ttk::DiagramType diagram{};
  int ret = VTUToDiagram(diagram, input, *this);
  if(ret != 0) {
    this->printErr("Could not read input Persistence Diagram");
    return 0;
  }

  std::array<PlotType, 4> plots{};
  ret = this->execute(plots, diagram);
  if(ret != 0) {
    this->printErr("Could not execute base code");
    return -1;
  }

  getPersistenceCurve(minSadTable, plots[0], "minimum-saddle pairs");
  getPersistenceCurve(sadSadTable, plots[1], "saddle-saddle pairs");
  getPersistenceCurve(sadMaxTable, plots[2], "maximum-saddle pairs");
  getPersistenceCurve(allPairsTable, plots[3], "all pairs");

  this->printMsg("Completed", 1, timer.getElapsedTime(), threadNumber_);
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
