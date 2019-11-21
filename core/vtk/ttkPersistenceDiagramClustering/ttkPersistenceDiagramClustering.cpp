#include <ttkPersistenceDiagramClustering.h>

#ifndef macroDiagramTuple
#define macroDiagramTuple                                                     \
  std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,               \
             ttk::CriticalType, VTK_TT, ttk::SimplexId, VTK_TT, float, float, \
             float, VTK_TT, float, float, float>
#endif

#ifndef macroMatchingTuple
#define macroMatchingTuple std::tuple<ttk::SimplexId, ttk::SimplexId, VTK_TT>
#endif

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceDiagramClustering)

  ttkPersistenceDiagramClustering::ttkPersistenceDiagramClustering() {
  UseAllCores = false;
  WassersteinMetric = "2";
  UseOutputMatching = true;
  TimeLimit = 9999999;
  NumberOfClusters = 1;
  Deterministic = 1;
  ThreadNumber = 1;
  PairTypeClustering = -1;
  numberOfInputsFromCommandLine = 1;
  UseProgressive = 1;
  UseAccelerated = 0;
  UseKmeansppInit = 0;
  Alpha = 1;
  DeltaLim = 0.01;
  Lambda = 1;
  Spacing = 1;
  oldSpacing = 1;
  Method = 0;
  needUpdate_ = true;
  UseInterruptible = true;
  UseAdditionalPrecision = false;
  ForceUseOfAlgorithm = false;
  DistanceWritingOptions = 0;
  DisplayMethod = 0;

  final_centroids_ = NULL;
  intermediateDiagrams_ = NULL;
  all_matchings_ = NULL;

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
}

ttkPersistenceDiagramClustering::~ttkPersistenceDiagramClustering() {
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkPersistenceDiagramClustering::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkPersistenceDiagramClustering::updateProgress(const float &progress) {
  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramClustering] " << progress * 100
        << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

template <typename VTK_TT>
int ttkPersistenceDiagramClustering::dispatch(
  int numInputs,
  std::vector<vtkUnstructuredGrid *> &inputDiagram,
  vtkUnstructuredGrid *outputClusters,
  vtkUnstructuredGrid *outputCentroids,
  vtkUnstructuredGrid *outputMatchings) {

  int ret{};
  vector<vector<macroDiagramTuple>> *intermediateDiagrams;
  vector<vector<vector<macroMatchingTuple>>> *all_matchings;
  vector<vector<macroDiagramTuple>> *final_centroids;
  if(needUpdate_) {
    if(intermediateDiagrams_) {
      vector<vector<macroDiagramTuple>> *tmpPTR
        = (vector<vector<macroDiagramTuple>> *)intermediateDiagrams_;
      delete tmpPTR;
    }
    intermediateDiagrams_ = new vector<vector<macroDiagramTuple>>(numInputs);

    if(final_centroids_) {
      vector<vector<macroDiagramTuple>> *tmpPTR
        = (vector<vector<macroDiagramTuple>> *)final_centroids_;
      delete tmpPTR;
    }
    final_centroids_ = new vector<vector<macroDiagramTuple>>;

    if(all_matchings_) {
      vector<vector<vector<macroMatchingTuple>>> *tmpPTR
        = (vector<vector<vector<macroMatchingTuple>>> *)all_matchings_;
      delete tmpPTR;
    }
    all_matchings_ = new vector<vector<vector<macroMatchingTuple>>>(3);
  }

  final_centroids = (vector<vector<macroDiagramTuple>> *)final_centroids_;
  intermediateDiagrams
    = (vector<vector<macroDiagramTuple>> *)intermediateDiagrams_;
  all_matchings = (vector<vector<vector<macroMatchingTuple>>> *)all_matchings_;

  if(needUpdate_) {

    max_dimension_total_ = 0;
    for(int i = 0; i < numInputs; i++) {
      double max_dimension = getPersistenceDiagram<VTK_TT>(
        &(intermediateDiagrams->at(i)), inputDiagram[i], Spacing, 0);
      if(max_dimension_total_ < max_dimension) {
        max_dimension_total_ = max_dimension;
      }
    }

    if(Method == 0) {
      // Progressive approach
      PersistenceDiagramClustering<VTK_TT> persistenceDiagramsClustering;
      persistenceDiagramsClustering.setWrapper(this);

      string wassersteinMetric = WassersteinMetric;

      if(!UseInterruptible) {
        TimeLimit = 999999999;
      }
      persistenceDiagramsClustering.setWasserstein(wassersteinMetric);
      persistenceDiagramsClustering.setDeterministic(Deterministic);
      persistenceDiagramsClustering.setForceUseOfAlgorithm(ForceUseOfAlgorithm);
      persistenceDiagramsClustering.setPairTypeClustering(PairTypeClustering);
      persistenceDiagramsClustering.setNumberOfInputs(numInputs);
      persistenceDiagramsClustering.setDebugLevel(debugLevel_);
      persistenceDiagramsClustering.setTimeLimit(TimeLimit);
      persistenceDiagramsClustering.setUseProgressive(UseProgressive);
      persistenceDiagramsClustering.setThreadNumber(threadNumber_);
      persistenceDiagramsClustering.setAlpha(Alpha);
      persistenceDiagramsClustering.setDeltaLim(DeltaLim);
      persistenceDiagramsClustering.setUseDeltaLim(UseAdditionalPrecision);
      persistenceDiagramsClustering.setLambda(Lambda);
      persistenceDiagramsClustering.setNumberOfClusters(NumberOfClusters);
      persistenceDiagramsClustering.setUseAccelerated(UseAccelerated);
      persistenceDiagramsClustering.setUseKmeansppInit(UseKmeansppInit);
      persistenceDiagramsClustering.setDistanceWritingOptions(
        DistanceWritingOptions);

      persistenceDiagramsClustering.setDiagrams((void *)intermediateDiagrams);
      inv_clustering_
        = persistenceDiagramsClustering.execute(final_centroids, all_matchings);

      needUpdate_ = false;
    }

    else {
      // AUCTION APPROACH
      final_centroids->resize(1);
      inv_clustering_.resize(numInputs);
      for(int i_input = 0; i_input < numInputs; i_input++) {
        inv_clustering_[i_input] = 0;
      }
      PersistenceDiagramBarycenter<VTK_TT> persistenceDiagramsBarycenter;
      persistenceDiagramsBarycenter.setWrapper(this);

      string wassersteinMetric = WassersteinMetric;
      persistenceDiagramsBarycenter.setWasserstein(wassersteinMetric);
      persistenceDiagramsBarycenter.setMethod(2);
      persistenceDiagramsBarycenter.setNumberOfInputs(numInputs);
      persistenceDiagramsBarycenter.setTimeLimit(TimeLimit);
      persistenceDiagramsBarycenter.setDeterministic(Deterministic);
      persistenceDiagramsBarycenter.setUseProgressive(UseProgressive);
      persistenceDiagramsBarycenter.setDebugLevel(debugLevel_);
      persistenceDiagramsBarycenter.setThreadNumber(threadNumber_);
      persistenceDiagramsBarycenter.setAlpha(Alpha);
      persistenceDiagramsBarycenter.setLambda(Lambda);
      // persistenceDiagramsBarycenter.setReinitPrices(ReinitPrices);
      // persistenceDiagramsBarycenter.setEpsilonDecreases(EpsilonDecreases);
      // persistenceDiagramsBarycenter.setEarlyStoppage(EarlyStoppage);

      persistenceDiagramsBarycenter.setDiagrams((void *)intermediateDiagrams);

      persistenceDiagramsBarycenter.execute(
        &(final_centroids->at(0)), all_matchings);

      needUpdate_ = false;
    }
  }

  outputMatchings->ShallowCopy(
    createMatchings(final_centroids, inv_clustering_, *intermediateDiagrams,
                    all_matchings, max_dimension_total_, Spacing));
  outputClusters->ShallowCopy(createOutputClusteredDiagrams(
    *intermediateDiagrams, inv_clustering_, max_dimension_total_, Spacing));
  outputCentroids->ShallowCopy(createOutputCentroids<VTK_TT>(
    final_centroids, inv_clustering_, max_dimension_total_, Spacing));

  return ret;
}

int ttkPersistenceDiagramClustering::doIt(vtkDataSet **input,
                                          vtkUnstructuredGrid *outputClusters,
                                          vtkUnstructuredGrid *outputCentroids,
                                          vtkUnstructuredGrid *outputMatchings,
                                          int numInputs) {
  // Get arrays from input datas
  // vtkDataArray* inputDiagram[numInputs] = { NULL };
  //
  //
  vector<vtkUnstructuredGrid *> inputDiagram(numInputs);
  for(int i = 0; i < numInputs; ++i) {
    inputDiagram[i] = vtkUnstructuredGrid::SafeDownCast(input[i]);
  }
  // Calling the executing package

  int dataType
    = inputDiagram[0]->GetCellData()->GetArray("Persistence")->GetDataType();

  switch(dataType) {
    vtkTemplateMacro(dispatch<VTK_TT>(numInputs, inputDiagram, outputClusters,
                                      outputCentroids, outputMatchings));
  }

  return 0;
}

int ttkPersistenceDiagramClustering::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(!this->Superclass::FillInputPortInformation(port, info)) {
    return 0;
  }
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}

int ttkPersistenceDiagramClustering::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(!this->Superclass::FillOutputPortInformation(port, info)) {
    return 0;
  }
  if(port == 0 || port == 1 || port == 3)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramClustering::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  Memory m;

  // Number of input files
  int numInputs = numberOfInputsFromCommandLine;
  if(numInputs == 1) {
    numInputs = inputVector[0]->GetNumberOfInformationObjects();
  }
  {
    stringstream msg;
    dMsg(cout, msg.str(), infoMsg);
  }
  // Get input datas
  vtkDataSet **input = new vtkDataSet *[numInputs];
  for(int i = 0; i < numInputs; ++i) {
    if(numberOfInputsFromCommandLine > 1) {
      input[i] = vtkDataSet::GetData(inputVector[i], 0);
    } else {
      input[i] = vtkDataSet::GetData(inputVector[0], i);
    }
    if(!input[i]) {
      std::cout << "No data in input[" << i << "]" << std::endl;
    } else if(this->GetMTime() < input[i]->GetMTime()) {
      needUpdate_ = true;
    }
  }
  // TODO Set output
  vtkInformation *outInfo1;
  outInfo1 = outputVector->GetInformationObject(0);
  vtkDataSet *output1
    = vtkDataSet::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *output_clusters
    = vtkUnstructuredGrid::SafeDownCast(output1);

  vtkInformation *outInfo2;
  outInfo2 = outputVector->GetInformationObject(1);
  vtkDataSet *output2
    = vtkDataSet::SafeDownCast(outInfo2->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *output_centroids
    = vtkUnstructuredGrid::SafeDownCast(output2);

  vtkInformation *outInfo3;
  outInfo3 = outputVector->GetInformationObject(2);
  vtkDataSet *output3
    = vtkDataSet::SafeDownCast(outInfo3->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *output_matchings
    = vtkUnstructuredGrid::SafeDownCast(output3);

  doIt(input, output_clusters, output_centroids, output_matchings, numInputs);
  delete[] input;

  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramClustering] Memory usage: "
        << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
