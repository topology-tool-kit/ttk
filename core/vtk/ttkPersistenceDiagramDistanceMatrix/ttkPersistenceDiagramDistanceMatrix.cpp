#include <ttkPersistenceDiagramDistanceMatrix.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceDiagramDistanceMatrix);

ttkPersistenceDiagramDistanceMatrix::ttkPersistenceDiagramDistanceMatrix() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkPersistenceDiagramDistanceMatrix::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkPersistenceDiagramDistanceMatrix::updateProgress(const float &progress) {
  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramDistanceMatrix] " << progress * 100
        << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkPersistenceDiagramDistanceMatrix::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagramDistanceMatrix::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramDistanceMatrix::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  Memory m;

  // hard-coded number of clusters
  this->NumberOfClusters = 4;

  // Number of input files
  int numInputs = numberOfInputsFromCommandLine;

  // Get input data
  std::vector<vtkUnstructuredGrid *> inputDiagrams;

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);

  if(blocks != nullptr) {
    numInputs = blocks->GetNumberOfBlocks();
    inputDiagrams.resize(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      inputDiagrams[i] = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
      if(this->GetMTime() < inputDiagrams[i]->GetMTime()) {
        needUpdate_ = true;
      }
    }
  }

  // Set output
  auto diagramsDistTable = vtkTable::GetData(outputVector);

  if(needUpdate_) {
    intermediateDiagrams_.resize(numInputs);
    all_matchings_.resize(3);

    max_dimension_total_ = 0;
    for(int i = 0; i < numInputs; i++) {
      double max_dimension
        = getPersistenceDiagram(intermediateDiagrams_[i], inputDiagrams[i]);
      if(max_dimension_total_ < max_dimension) {
        max_dimension_total_ = max_dimension;
      }
    }

    if(Method == 0) {
      // Progressive approach
      PersistenceDiagramDistanceMatrix<double> persistenceDiagramsClustering;
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
      persistenceDiagramsClustering.setUseFullDiagrams(UseFullDiagrams);
      persistenceDiagramsClustering.setPerClusterDistanceMatrix(
        PerClusterDistanceMatrix);

      inv_clustering_ = persistenceDiagramsClustering.execute(
        intermediateDiagrams_, final_centroids_, all_matchings_);

      diagramsDistMat = persistenceDiagramsClustering.getDiagramsDistMat();
      distanceToCentroid
        = persistenceDiagramsClustering.getDistanceToCentroid();

      needUpdate_ = false;
    }

    else {
      // AUCTION APPROACH
      final_centroids_.resize(1);
      inv_clustering_.resize(numInputs);
      for(int i_input = 0; i_input < numInputs; i_input++) {
        inv_clustering_[i_input] = 0;
      }
      PersistenceDiagramBarycenter<double> persistenceDiagramsBarycenter;
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
        intermediateDiagrams_, final_centroids_[0], all_matchings_);

      needUpdate_ = false;
    }
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

  // aggregate input field data
  vtkNew<vtkFieldData> fd{};
  fd->CopyStructure(inputDiagrams[0]->GetFieldData());
  fd->SetNumberOfTuples(inputDiagrams.size());
  for(size_t i = 0; i < inputDiagrams.size(); ++i) {
    fd->SetTuple(i, 0, inputDiagrams[i]->GetFieldData());
  }

  // copy input field data to output row data
  for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
    diagramsDistTable->AddColumn(fd->GetAbstractArray(i));
  }

  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramDistanceMatrix] Memory usage: "
        << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}

double ttkPersistenceDiagramDistanceMatrix::getPersistenceDiagram(
  std::vector<diagramType> &diagram,
  vtkUnstructuredGrid *CTPersistenceDiagram_) {
  vtkIntArray *vertexIdentifierScalars
    = vtkIntArray::SafeDownCast(CTPersistenceDiagram_->GetPointData()->GetArray(
      ttk::VertexScalarFieldName));

  vtkIntArray *nodeTypeScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("CriticalType"));

  vtkIntArray *pairIdentifierScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("PairIdentifier"));

  vtkIntArray *extremumIndexScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("PairType"));

  const auto persistenceScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("Persistence"));
  const auto birthScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Birth"));
  const auto deathScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Death"));

  const auto critCoordinates = vtkFloatArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Coordinates"));
  const auto points = CTPersistenceDiagram_->GetPoints();

  const bool embed = birthScalars != nullptr && deathScalars != nullptr;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!embed && critCoordinates == nullptr) {
    // missing data
    return -2;
  }
#endif // TTK_ENABLE_KAMIKAZE

  int pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();
  // FIX : no more missed pairs
  for(int pair_index = 0; pair_index < pairingsSize; pair_index++) {
    const float index_of_pair = pair_index;
    if(*pairIdentifierScalars->GetTuple(pair_index) != -1)
      pairIdentifierScalars->SetTuple(pair_index, &index_of_pair);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(pairingsSize < 1 || !vertexIdentifierScalars || !pairIdentifierScalars
     || !nodeTypeScalars || !persistenceScalars || !extremumIndexScalars
     || !points) {
    return -2;
  }
#endif // TTK_ENABLE_KAMIKAZE

  diagram.resize(pairingsSize + 1);
  int nbNonCompact = 0;
  double max_dimension = 0;

  // skip diagonal cell (corresponding points already dealt with)
  for(int i = 0; i < pairingsSize - 1; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
    int nodeType1 = nodeTypeScalars->GetValue(2 * i);
    int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    std::array<double, 3> coordsBirth{}, coordsDeath{};

    const auto i0 = 2 * i;
    const auto i1 = 2 * i + 1;

    double birth, death;

    if(embed) {
      points->GetPoint(i0, coordsBirth.data());
      points->GetPoint(i1, coordsDeath.data());
      birth = birthScalars->GetValue(i0);
      death = deathScalars->GetValue(i1);
    } else {
      critCoordinates->GetTuple(i0, coordsBirth.data());
      critCoordinates->GetTuple(i1, coordsDeath.data());
      birth = points->GetPoint(i0)[0];
      death = points->GetPoint(i1)[1];
    }

    if(pairIdentifier != -1 && pairIdentifier < pairingsSize) {
      if(pairIdentifier == 0) {
        max_dimension = persistence;

        diagram[0] = std::make_tuple(
          vertexId1, CriticalType::Local_minimum, vertexId2,
          CriticalType::Saddle1, persistence, pairType, birth, coordsBirth[0],
          coordsBirth[1], coordsBirth[2], death, coordsDeath[0], coordsDeath[1],
          coordsDeath[2]);
        diagram[pairingsSize] = std::make_tuple(
          vertexId1, CriticalType::Saddle1, vertexId2,
          CriticalType::Local_maximum, persistence, pairType, birth,
          coordsBirth[0], coordsBirth[1], coordsBirth[2], death, coordsDeath[0],
          coordsDeath[1], coordsDeath[2]);

      } else {
        diagram[pairIdentifier] = std::make_tuple(
          vertexId1, (BNodeType)nodeType1, vertexId2, (BNodeType)nodeType2,
          persistence, pairType, birth, coordsBirth[0], coordsBirth[1],
          coordsBirth[2], death, coordsDeath[0], coordsDeath[1],
          coordsDeath[2]);
      }
    }
    if(pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if(nbNonCompact == 0) {
        std::stringstream msg;
        msg << "[TTKPersistenceDiagramDistanceMatrix] Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). " << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }
    }
  }

  if(nbNonCompact > 0) {
    {
      std::stringstream msg;
      msg << "[TTKPersistenceDiagramDistanceMatrix] Missed " << nbNonCompact
          << " pairs due to non-compactness." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  return max_dimension;
}
