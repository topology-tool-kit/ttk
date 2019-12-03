#include <ttkPersistenceDiagramClustering.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceDiagramClustering)

ttkPersistenceDiagramClustering::ttkPersistenceDiagramClustering() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(5);
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
  const std::vector<vtkUnstructuredGrid *> &inputDiagram,
  vtkUnstructuredGrid *outputClusters,
  vtkUnstructuredGrid *outputCentroids,
  vtkUnstructuredGrid *outputMatchings,
  vtkTable *diagramsDistTable,
  vtkTable *centroidsDistTable) {

  using macroDiagramTuple
    = std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,
                 ttk::CriticalType, VTK_TT, ttk::SimplexId, VTK_TT, float,
                 float, float, VTK_TT, float, float, float>;
  using macroMatchingTuple = std::tuple<ttk::SimplexId, ttk::SimplexId, VTK_TT>;

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
        intermediateDiagrams->at(i), inputDiagram[i], Spacing, 0);
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
      persistenceDiagramsClustering.setOutputDistanceMatrix(
        OutputDistanceMatrix);
      persistenceDiagramsClustering.setPerClusterDistanceMatrix(
        PerClusterDistanceMatrix);

      inv_clustering_ = persistenceDiagramsClustering.execute(
        *intermediateDiagrams, *final_centroids, *all_matchings);

      diagramsDistMat = persistenceDiagramsClustering.getDiagramsDistMat();
      distanceToCentroid
        = persistenceDiagramsClustering.getDistanceToCentroid();
      centroidsDistMat = persistenceDiagramsClustering.getCentroidsDistMat();

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

      persistenceDiagramsBarycenter.execute(
        *intermediateDiagrams, final_centroids->at(0), *all_matchings);

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

  if(!OutputDistanceMatrix) {
    // early return
    return ret;
  }

  // zero-padd column name to keep Row Data columns ordered
  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  // copy diagrams distance matrix to output
  for(size_t i = 0; i < diagramsDistMat.size(); ++i) {
    std::string name{"Diagram"};
    zeroPad(name, diagramsDistMat.size(), i);

    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(numInputs);
    col->SetName(name.c_str());
    for(size_t j = 0; j < diagramsDistMat[i].size(); ++j) {
      col->SetTuple1(j, diagramsDistMat[i][j]);
    }
    diagramsDistTable->AddColumn(col);
  }

  // add clusters id to diagrams distance matrix output
  vtkNew<vtkIntArray> clusterId{};
  clusterId->SetName("ClusterId");
  clusterId->SetArray(inv_clustering_.data(), inv_clustering_.size(), 1);
  diagramsDistTable->AddColumn(clusterId);

  // add distance to cluster centroid to diagrams distance matrix output
  vtkNew<vtkDoubleArray> dCentroid{};
  dCentroid->SetName("DistanceToCentroid");
  dCentroid->SetArray(distanceToCentroid.data(), distanceToCentroid.size(), 1);
  diagramsDistTable->AddColumn(dCentroid);

  // aggregate input field data into diagrams distance matrix output field data
  auto fd = diagramsDistTable->GetFieldData();
  fd->CopyStructure(inputDiagram[0]->GetFieldData());
  fd->SetNumberOfTuples(inputDiagram.size());
  for(size_t i = 0; i < inputDiagram.size(); ++i) {
    fd->SetTuple(i, 0, inputDiagram[i]->GetFieldData());
  }

  // also copy field data arrays to row data
  for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
    diagramsDistTable->AddColumn(fd->GetAbstractArray(i));
  }

  // copy centroids distance matrix to output
  for(size_t i = 0; i < centroidsDistMat.size(); ++i) {
    std::string name{"Cluster"};
    zeroPad(name, centroidsDistMat.size(), i);

    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(centroidsDistMat.size());
    col->SetName(name.c_str());
    for(size_t j = 0; j < centroidsDistMat[i].size(); ++j) {
      col->SetTuple1(j, centroidsDistMat[i][j]);
    }
    centroidsDistTable->AddColumn(col);
  }

  return ret;
}

int ttkPersistenceDiagramClustering::doIt(
  const std::vector<vtkUnstructuredGrid *> &input,
  vtkUnstructuredGrid *outputClusters,
  vtkUnstructuredGrid *outputCentroids,
  vtkUnstructuredGrid *outputMatchings,
  vtkTable *diagramsDistTable,
  vtkTable *centroidsDistTable,
  int numInputs) {

  // Calling the executing package
  int dataType
    = input[0]->GetCellData()->GetArray("Persistence")->GetDataType();

  switch(dataType) {
    vtkTemplateMacro(dispatch<VTK_TT>(numInputs, input, outputClusters,
                                      outputCentroids, outputMatchings,
                                      diagramsDistTable, centroidsDistTable));
  }

  return 0;
}

int ttkPersistenceDiagramClustering::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(!this->Superclass::FillInputPortInformation(port, info)) {
    return 0;
  }
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
}

int ttkPersistenceDiagramClustering::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(!this->Superclass::FillOutputPortInformation(port, info)) {
    return 0;
  }
  if(port == 0 || port == 1 || port == 2)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  if(port == 3 || port == 4) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
  }
  return 1;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramClustering::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  Memory m;

  // Number of input files
  int numInputs = numberOfInputsFromCommandLine;

  // Get input data
  std::vector<vtkUnstructuredGrid *> input;

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);

  if(blocks != nullptr) {
    numInputs = blocks->GetNumberOfBlocks();
    input.resize(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      input[i] = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
      if(this->GetMTime() < input[i]->GetMTime()) {
        needUpdate_ = true;
      }
    }
  }

  // Set outputs
  auto output_clusters = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  auto output_centroids = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));
  auto output_matchings = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(2)->Get(vtkDataObject::DATA_OBJECT()));
  auto output_diagrams_matrix = vtkTable::SafeDownCast(
    outputVector->GetInformationObject(3)->Get(vtkDataObject::DATA_OBJECT()));
  auto output_centroids_matrix = vtkTable::SafeDownCast(
    outputVector->GetInformationObject(4)->Get(vtkDataObject::DATA_OBJECT()));

  doIt(input, output_clusters, output_centroids, output_matchings,
       output_diagrams_matrix, output_centroids_matrix, numInputs);

  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramClustering] Memory usage: "
        << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
