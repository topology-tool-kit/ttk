#include <VectorSimplification.h>
#include <ttkMacros.h>
#include <ttkUtils.h>
#include <ttkVectorWeightCurve.h>

#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkVectorWeightCurve);

ttkVectorWeightCurve::ttkVectorWeightCurve() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

vtkTable *ttkVectorWeightCurve::GetOutput() {
  return this->GetOutput(0);
}

vtkTable *ttkVectorWeightCurve::GetOutput(int port) {
  return vtkTable::SafeDownCast(this->GetOutputDataObject(port));
}

int ttkVectorWeightCurve::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkVectorWeightCurve::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

static void getVectorWeightCurve(
  vtkTable *outputCurve,
  std::vector<ttk::VectorSimplification::PlotPoint> &plotPoints,
  bool displayCycles,
  bool displayCancelTypes) {

  if(plotPoints.empty()) {
    return;
  }

  vtkNew<vtkDoubleArray> persistenceScalars{};
  persistenceScalars->SetName("Weight Value");
  persistenceScalars->SetNumberOfTuples(plotPoints.size());

  vtkNew<ttkSimplexIdTypeArray> numberOfCPScalars{};
  numberOfCPScalars->SetName("Number of Critical Points");
  numberOfCPScalars->SetNumberOfTuples(plotPoints.size());

  vtkNew<ttkSimplexIdTypeArray> numberOfSinksScalars{};
  vtkNew<ttkSimplexIdTypeArray> numberOfSourcesScalars{};
  if(displayCancelTypes) {
    numberOfSinksScalars->SetName("Number of Sinks");
    numberOfSinksScalars->SetNumberOfTuples(plotPoints.size());

    numberOfSourcesScalars->SetName("Number of Sources");
    numberOfSourcesScalars->SetNumberOfTuples(plotPoints.size());
  }

  vtkNew<ttkSimplexIdTypeArray> numberOfOrbitsAdded{};
  vtkNew<ttkSimplexIdTypeArray> numberOfOrbitsRemoved{};
  if(displayCycles) {
    numberOfOrbitsAdded->SetName("Number of Orbits Generated");
    numberOfOrbitsAdded->SetNumberOfTuples(plotPoints.size());

    numberOfOrbitsRemoved->SetName("Number of Orbits Cancelled");
    numberOfOrbitsRemoved->SetNumberOfTuples(plotPoints.size());
  }

  for(size_t i = 0; i < plotPoints.size(); ++i) {
    persistenceScalars->SetTuple1(i, plotPoints[i].weight);
    numberOfCPScalars->SetTuple1(i, plotPoints[i].numCP);
    if(displayCancelTypes) {
      numberOfSinksScalars->SetTuple1(i, plotPoints[i].numSinks);
      numberOfSourcesScalars->SetTuple1(i, plotPoints[i].numSources);
    }
    if(displayCycles) {
      numberOfOrbitsAdded->SetTuple1(i, plotPoints[i].orbitsAdded);
      numberOfOrbitsRemoved->SetTuple1(i, plotPoints[i].orbitsRemoved);
    }
  }

  vtkNew<vtkTable> vectorWeightCurve{};
  vectorWeightCurve->AddColumn(persistenceScalars);
  vectorWeightCurve->AddColumn(numberOfCPScalars);
  if(displayCancelTypes) {
    vectorWeightCurve->AddColumn(numberOfSinksScalars);
    vectorWeightCurve->AddColumn(numberOfSourcesScalars);
  }
  if(displayCycles) {
    vectorWeightCurve->AddColumn(numberOfOrbitsAdded);
    vectorWeightCurve->AddColumn(numberOfOrbitsRemoved);
  }

  outputCurve->ShallowCopy(vectorWeightCurve);
}

int ttkVectorWeightCurve::RequestData(vtkInformation *ttkNotUsed(request),
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  ttk::Timer timer;

  auto input = vtkDataSet::GetData(inputVector[0]);

  auto triangulation = ttkAlgorithm::GetTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }
  if(!input->GetNumberOfPoints()) {
    this->printErr("Input has no point.");
    return -1;
  }
#endif

  this->vs_.preconditionTriangulation(triangulation);

  const auto inputVectors = this->GetInputArrayToProcess(0, input);

  if(inputVectors == nullptr) {
    this->printErr("Input vector array is NULL");
    return 0;
  }
  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  if(inputVectors->GetNumberOfComponents() != 3) {
    this->printErr("Input array needs to be a vector array(3D).");
    return 0;
  }

  auto *allPairsTable = vtkTable::GetData(outputVector, 0);

  int ret{};

  // Generate and simplify the discrete field
  ttkVtkTemplateMacro(
    inputVectors->GetDataType(), triangulation->getType(),
    (ret = this->vs_.buildField<VTK_TT, TTK_TT>(
       ttkUtils::GetVoidPointer(inputVectors), inputVectors->GetMTime(),
       *static_cast<TTK_TT *>(triangulation->getData()))));
  if(ret != 0) {
    this->printErr("Could not generate field");
    return -1;
  }
  this->printMsg("Generated the discrete field");
  std::vector<ttk::VectorSimplification::PlotPoint> vs_pairs;
  this->UpdateTraceFullOrbits();

  ttkVtkTemplateMacro(
    inputVectors->GetDataType(), triangulation->getType(),
    (ret = this->vs_.performSimplification<VTK_TT, TTK_TT>(
       0, true, vs_pairs, *static_cast<TTK_TT *>(triangulation->getData()))));
  if(ret != 0) {
    this->printErr("Could not simplify field");
    return -1;
  }
  this->printMsg("Simplified the discete field");

  getVectorWeightCurve(
    allPairsTable, vs_pairs, this->DisplayOrbits, this->DisplayExtrema);

  this->printMsg("Completed", 1, timer.getElapsedTime(), threadNumber_);
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
