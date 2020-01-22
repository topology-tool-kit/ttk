#include <ttkPersistenceDiagramClustering.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceDiagramClustering)

  ttkPersistenceDiagramClustering::ttkPersistenceDiagramClustering() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(4);
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

int ttkPersistenceDiagramClustering::doIt(
  const std::vector<vtkUnstructuredGrid *> &inputDiagram,
  vtkUnstructuredGrid *outputClusters,
  vtkUnstructuredGrid *outputCentroids,
  vtkUnstructuredGrid *outputMatchings,
  vtkTable *diagramsDistTable,
  int numInputs) {

  int ret{};
  if(needUpdate_) {
    intermediateDiagrams_.resize(numInputs);
    all_matchings_.resize(3);

    max_dimension_total_ = 0;
    for(int i = 0; i < numInputs; i++) {
      double max_dimension
        = getPersistenceDiagram(intermediateDiagrams_[i], inputDiagram[i]);
      if(max_dimension_total_ < max_dimension) {
        max_dimension_total_ = max_dimension;
      }
    }

    if(Method == 0) {
      // Progressive approach
      PersistenceDiagramClustering<double> persistenceDiagramsClustering;
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

  outputMatchings->ShallowCopy(createMatchings());
  outputClusters->ShallowCopy(createOutputClusteredDiagrams());
  outputCentroids->ShallowCopy(createOutputCentroids());

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

  return ret;
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
  if(port == 3) {
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

  doIt(input, output_clusters, output_centroids, output_matchings,
       output_diagrams_matrix, numInputs);

  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramClustering] Memory usage: "
        << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}

double ttkPersistenceDiagramClustering::getPersistenceDiagram(
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

  if(NumberOfClusters == 1) {
    diagram.resize(pairingsSize);
  } else {
    diagram.resize(pairingsSize + 1);
  }
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

        if(NumberOfClusters == 1) {
          diagram[0] = std::make_tuple(
            vertexId1, CriticalType::Local_minimum, vertexId2,
            CriticalType::Local_maximum, persistence, pairType, birth,
            coordsBirth[0], coordsBirth[1], coordsBirth[2], death,
            coordsDeath[0], coordsDeath[1], coordsDeath[2]);
        } else {
          diagram[0] = std::make_tuple(
            vertexId1, CriticalType::Local_minimum, vertexId2,
            CriticalType::Saddle1, persistence, pairType, birth, coordsBirth[0],
            coordsBirth[1], coordsBirth[2], death, coordsDeath[0],
            coordsDeath[1], coordsDeath[2]);
          diagram[pairingsSize] = std::make_tuple(
            vertexId1, CriticalType::Saddle1, vertexId2,
            CriticalType::Local_maximum, persistence, pairType, birth,
            coordsBirth[0], coordsBirth[1], coordsBirth[2], death,
            coordsDeath[0], coordsDeath[1], coordsDeath[2]);
        }

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
        msg << "[TTKPersistenceDiagramClustering] Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). " << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }
    }
  }

  if(nbNonCompact > 0) {
    {
      std::stringstream msg;
      msg << "[TTKPersistenceDiagramClustering] Missed " << nbNonCompact
          << " pairs due to non-compactness." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  return max_dimension;
}

vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceDiagramClustering::createOutputCentroids() {
  if(debugLevel_ > 5) {
    std::cout << "[ttkPersistenceDiagramClustering] Creating vtk diagrams"
              << std::endl;
  }
  vtkNew<vtkPoints> points{};

  vtkNew<vtkUnstructuredGrid> persistenceDiagram{};

  vtkNew<vtkIntArray> nodeType{};
  nodeType->SetName("CriticalType");

  vtkNew<vtkDoubleArray> persistenceScalars{};
  persistenceScalars->SetName("Persistence");

  vtkNew<vtkIntArray> idOfPair{};
  idOfPair->SetName("PairID");

  vtkNew<vtkDoubleArray> persistenceScalarsPoint{};
  persistenceScalarsPoint->SetName("Persistence");

  vtkNew<vtkIntArray> idOfDiagramPoint{};
  idOfDiagramPoint->SetName("ClusterID");

  vtkNew<vtkIntArray> pairType{};
  pairType->SetName("PairType");

  vtkNew<vtkFloatArray> coordsScalars{};
  coordsScalars->SetNumberOfComponents(3);
  coordsScalars->SetName("Coordinates");

  int count = 0;
  for(unsigned int j = 0; j < final_centroids_.size(); ++j) {
    const std::vector<diagramType> &diagram = final_centroids_[j];

    // First, add diagram points to the global input diagram
    for(unsigned int i = 0; i < diagram.size(); ++i) {
      vtkIdType ids[2];
      const diagramType &t = diagram[i];
      double x1 = std::get<6>(t);
      double y1 = x1;
      if(DisplayMethod == 1 && Spacing != 0) {
        x1 += 3 * (abs(Spacing) + 0.2) * max_dimension_total_ * j;
      }
      double z1 = 0; // Change 1 to j if you want to isolate the diagrams

      float coords1[3];
      coords1[0] = std::get<7>(t);
      coords1[1] = std::get<8>(t);
      coords1[2] = std::get<9>(t);

      double x2 = std::get<6>(t);
      double y2 = std::get<10>(t);
      double z2 = 0; // Change 1 to j if you want to isolate the

      float coords2[3];
      coords2[0] = std::get<11>(t);
      coords2[1] = std::get<12>(t);
      coords2[2] = std::get<13>(t);

      idOfPair->InsertTuple1(count, i);

      points->InsertNextPoint(x1, y1, z1);
      coordsScalars->InsertTuple3(
        2 * count, coords1[0], coords1[1], coords1[2]);
      idOfDiagramPoint->InsertTuple1(2 * count, j);
      const ttk::CriticalType n1Type = std::get<1>(t);
      switch(n1Type) {
        case BLocalMin:
          nodeType->InsertTuple1(2 * count, 0);
          break;

        case BSaddle1:
          nodeType->InsertTuple1(2 * count, 1);
          break;

        case BSaddle2:
          nodeType->InsertTuple1(2 * count, 2);
          break;

        case BLocalMax:
          nodeType->InsertTuple1(2 * count, 3);
          break;
        default:
          nodeType->InsertTuple1(2 * count, 0);
      }
      if(DisplayMethod == 1 && Spacing != 0) {
        points->InsertNextPoint(
          x2 + 3 * (abs(Spacing) + 0.2) * max_dimension_total_ * j, y2, z2);
      } else {
        points->InsertNextPoint(x2, y2, z2);
      }
      coordsScalars->InsertTuple3(
        2 * count + 1, coords2[0], coords2[1], coords2[2]);
      idOfDiagramPoint->InsertTuple1(2 * count + 1, j);
      const ttk::CriticalType n2Type = std::get<3>(t);
      switch(n2Type) {
        case BLocalMin:
          nodeType->InsertTuple1(2 * count + 1, 0);
          break;

        case BSaddle1:
          nodeType->InsertTuple1(2 * count + 1, 1);
          break;

        case BSaddle2:
          nodeType->InsertTuple1(2 * count + 1, 2);
          break;

        case BLocalMax:
          nodeType->InsertTuple1(2 * count + 1, 3);
          break;
        default:
          nodeType->InsertTuple1(2 * count + 1, 0);
      }

      ids[0] = 2 * count;
      ids[1] = 2 * count + 1;

      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      persistenceScalars->InsertTuple1(count, y2 - x2);
      persistenceScalarsPoint->InsertTuple1(2 * count, y2 - x2);
      persistenceScalarsPoint->InsertTuple1(2 * count + 1, y2 - x2);
      const ttk::SimplexId type = std::get<5>(t);
      switch(type) {
        case 0:
          pairType->InsertTuple1(count, 0);
          break;

        case 1:
          pairType->InsertTuple1(count, 1);
          break;

        case 2:
          pairType->InsertTuple1(count, 2);
          break;
        default:
          pairType->InsertTuple1(count, 0);
      }
      count++;
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(pairType);
  persistenceDiagram->GetCellData()->AddArray(idOfPair);
  persistenceDiagram->GetPointData()->AddArray(nodeType);
  persistenceDiagram->GetPointData()->AddArray(coordsScalars);
  persistenceDiagram->GetPointData()->AddArray(idOfDiagramPoint);
  persistenceDiagram->GetPointData()->AddArray(persistenceScalarsPoint);

  return persistenceDiagram;
}

vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceDiagramClustering::createOutputClusteredDiagrams() {
  if(debugLevel_ > 5) {
    std::cout << "[ttkPersistenceDiagramClustering] Creating vtk Outputs"
              << std::endl;
  }
  vtkNew<vtkPoints> points{};

  vtkNew<vtkUnstructuredGrid> persistenceDiagram{};

  vtkNew<vtkIntArray> nodeType{};
  nodeType->SetName("CriticalType");

  vtkNew<vtkDoubleArray> persistenceScalars{};
  persistenceScalars->SetName("Persistence");

  vtkNew<vtkIntArray> idOfPair{};
  idOfPair->SetName("PairID");

  vtkNew<vtkDoubleArray> persistenceScalarsPoint{};
  persistenceScalarsPoint->SetName("Persistence");

  vtkNew<vtkIntArray> idOfDiagramPoint{};
  idOfDiagramPoint->SetName("DiagramID");

  vtkNew<vtkIntArray> idOfCluster{};
  idOfCluster->SetName("ClusterID");

  vtkNew<vtkIntArray> pairType{};
  pairType->SetName("PairType");

  vtkNew<vtkFloatArray> coordsScalars{};
  coordsScalars->SetNumberOfComponents(3);
  coordsScalars->SetName("Coordinates");

  vtkNew<ttkSimplexIdTypeArray> vertexSField{};
  vertexSField->SetName(ttk::VertexScalarFieldName);
  vertexSField->SetNumberOfComponents(1);

  std::vector<int> cluster_size;
  std::vector<int> idxInCluster(intermediateDiagrams_.size());
  for(unsigned int j = 0; j < intermediateDiagrams_.size(); ++j) {
    idxInCluster[j] = 0;
  }
  // RE-Invert clusters
  if(Spacing > 0) {
    for(unsigned int j = 0; j < intermediateDiagrams_.size(); ++j) {
      unsigned int c = inv_clustering_[j];
      if(c + 1 > cluster_size.size()) {
        cluster_size.resize(c + 1);
        cluster_size[c] = 1;
        idxInCluster[j] = 0;
      } else {
        cluster_size[c]++;
        idxInCluster[j] = cluster_size[c] - 1;
      }
    }
  }
  int count = 0;
  for(unsigned int j = 0; j < intermediateDiagrams_.size(); ++j) {
    const std::vector<diagramType> &diagram = intermediateDiagrams_[j];

    unsigned int c = inv_clustering_[j];
    // First, add diagram points to the global input diagram
    for(unsigned int i = 0; i < diagram.size(); ++i) {
      vtkIdType ids[2];
      const diagramType t = diagram[i];
      double x1 = std::get<6>(t);
      double y1 = x1;
      double z1 = 0;

      float coords1[3];
      coords1[0] = std::get<7>(t);
      coords1[1] = std::get<8>(t);
      coords1[2] = std::get<9>(t);
      double x2 = std::get<6>(t);
      double y2 = std::get<10>(t);
      double z2 = 0;
      if(DisplayMethod == 1 && Spacing > 0) {
        // cout<<"j "<<j<<" size "<<cluster_size[inv_clustering[j]]<<endl;
        // cout<<"count "<<count_diagram<<endl;
        double angle = 2 * 3.1415926 * (double)(idxInCluster[j])
                       / cluster_size[inv_clustering_[j]];
        x1 += (abs(Spacing) + .2) * 3 * max_dimension_total_ * c
              + Spacing * max_dimension_total_ * cos(angle);
        x2 += (abs(Spacing) + .2) * 3 * max_dimension_total_ * c
              + Spacing * max_dimension_total_ * cos(angle);
        y1 += Spacing * max_dimension_total_ * sin(angle);
        y2 += Spacing * max_dimension_total_ * sin(angle);
      } else if(DisplayMethod == 2) {
        z2 = Spacing;
        z1 = Spacing;
        if(j == 0) {
          z2 = -Spacing;
          z1 = -Spacing;
        }
      }

      float coords2[3];
      coords2[0] = std::get<11>(t);
      coords2[1] = std::get<12>(t);
      coords2[2] = std::get<13>(t);

      idOfPair->InsertTuple1(count, i);

      points->InsertNextPoint(x1, y1, z1);
      coordsScalars->InsertTuple3(
        2 * count, coords1[0], coords1[1], coords1[2]);
      idOfDiagramPoint->InsertTuple1(2 * count, j);
      // std::cout<<"\nMAX DIM \n"<<max_dimension<<std::endl;
      idOfCluster->InsertTuple1(2 * count, c);
      const ttk::CriticalType n1Type = std::get<1>(t);
      switch(n1Type) {
        case BLocalMin:
          nodeType->InsertTuple1(2 * count, 0);
          break;

        case BSaddle1:
          nodeType->InsertTuple1(2 * count, 1);
          break;

        case BSaddle2:
          nodeType->InsertTuple1(2 * count, 2);
          break;

        case BLocalMax:
          nodeType->InsertTuple1(2 * count, 3);
          break;
        default:
          nodeType->InsertTuple1(2 * count, 0);
      }

      points->InsertNextPoint(x2, y2, z2);
      coordsScalars->InsertTuple3(
        2 * count + 1, coords2[0], coords2[1], coords2[2]);
      idOfDiagramPoint->InsertTuple1(2 * count + 1, j);
      idOfCluster->InsertTuple1(2 * count + 1, c);
      const ttk::CriticalType n2Type = std::get<3>(t);
      switch(n2Type) {
        case BLocalMin:
          nodeType->InsertTuple1(2 * count + 1, 0);
          break;

        case BSaddle1:
          nodeType->InsertTuple1(2 * count + 1, 1);
          break;

        case BSaddle2:
          nodeType->InsertTuple1(2 * count + 1, 2);
          break;

        case BLocalMax:
          nodeType->InsertTuple1(2 * count + 1, 3);
          break;
        default:
          nodeType->InsertTuple1(2 * count + 1, 0);
      }

      ids[0] = 2 * count;
      ids[1] = 2 * count + 1;

      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      persistenceScalars->InsertTuple1(count, y2 - x2);
      persistenceScalarsPoint->InsertTuple1(2 * count, y2 - x2);
      persistenceScalarsPoint->InsertTuple1(2 * count + 1, y2 - x2);
      const ttk::SimplexId type = std::get<5>(t);
      switch(type) {
        case 0:
          pairType->InsertTuple1(count, 0);
          break;

        case 1:
          pairType->InsertTuple1(count, 1);
          break;

        case 2:
          pairType->InsertTuple1(count, 2);
          break;
        default:
          pairType->InsertTuple1(count, 0);
      }
      vertexSField->InsertTuple1(2 * count, std::get<0>(t));
      vertexSField->InsertTuple1(2 * count + 1, std::get<2>(t));

      count++;
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(pairType);
  persistenceDiagram->GetCellData()->AddArray(idOfPair);
  persistenceDiagram->GetPointData()->AddArray(nodeType);
  persistenceDiagram->GetPointData()->AddArray(coordsScalars);
  persistenceDiagram->GetPointData()->AddArray(idOfDiagramPoint);
  persistenceDiagram->GetPointData()->AddArray(idOfCluster);
  persistenceDiagram->GetPointData()->AddArray(persistenceScalarsPoint);
  persistenceDiagram->GetPointData()->AddArray(vertexSField);

  return persistenceDiagram;
}

vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceDiagramClustering::createMatchings() {
  if(debugLevel_ > 5) {
    std::cout << "[ttkPersistenceDiagramClustering] Creating vtk Matchings"
              << std::endl;
  }
  vtkNew<vtkPoints> matchingPoints{};

  vtkNew<vtkUnstructuredGrid> matchingMesh{};

  vtkNew<vtkIntArray> idOfDiagramMatchingPoint{};
  idOfDiagramMatchingPoint->SetName("DiagramID");

  vtkNew<vtkIntArray> idOfPoint{};
  idOfPoint->SetName("PointID");

  vtkNew<vtkIntArray> idOfDiagramMatching{};
  idOfDiagramMatching->SetName("DiagramID");

  vtkNew<vtkIntArray> idOfCluster{};
  idOfCluster->SetName("ClusterID");

  vtkNew<vtkDoubleArray> cost{};
  cost->SetName("Cost");

  vtkNew<vtkIntArray> pairType{};
  pairType->SetName("PairType");

  vtkNew<vtkIntArray> matchingCount{};
  matchingCount->SetName("MatchNumber");

  std::vector<int> cluster_size;
  std::vector<int> idxInCluster(intermediateDiagrams_.size());

  std::vector<int> matchings_count(final_centroids_[0].size(), 0);
  std::vector<int> count_to_good;

  for(unsigned int j = 0; j < intermediateDiagrams_.size(); ++j) {
    idxInCluster[j] = 0;
  }
  // RE-Invert clusters
  if(DisplayMethod == 1 && Spacing > 0) {
    for(unsigned int j = 0; j < intermediateDiagrams_.size(); ++j) {
      unsigned int c = inv_clustering_[j];
      if(c + 1 > cluster_size.size()) {
        cluster_size.resize(c + 1);
        cluster_size[c] = 1;
        idxInCluster[j] = 0;
      } else {
        cluster_size[c]++;
        idxInCluster[j] = cluster_size[c] - 1;
      }
    }
  }
  int count = 0;
  for(unsigned int j = 0; j < intermediateDiagrams_.size(); ++j) {
    int c = inv_clustering_[j];
    const auto &diagram = intermediateDiagrams_[j];
    std::vector<matchingType> matchings_j
      = all_matchings_[inv_clustering_[j]][j];
    for(unsigned int i = 0; i < matchings_j.size(); ++i) {

      vtkIdType ids[2];
      ids[0] = 2 * count;
      ids[1] = 2 * count + 1;
      matchingType m = matchings_j[i];
      int bidder_id = std::get<0>(m);
      int good_id = std::get<1>(m);
      if(NumberOfClusters == 1 && good_id > -1) {
        matchings_count[good_id] += 1;
        count_to_good.push_back(good_id);
      }

      // avoid a weird out-of-bound access
      good_id
        = std::min(good_id, static_cast<int>(final_centroids_[c].size() - 1));

      diagramType t1 = final_centroids_[c][good_id];
      double x1 = std::get<6>(t1);
      double y1 = std::get<10>(t1);
      double z1 = 0;

      // endl;
      if(bidder_id < (int)diagram.size()) {
        diagramType t2 = diagram[bidder_id];
        double x2 = std::get<6>(t2);
        double y2 = std::get<10>(t2);
        double z2 = 0; // Change 1 to j if you want to isolate the diagrams

        if(DisplayMethod == 1 && Spacing > 0) {
          double angle
            = 2 * 3.1415926 * (double)(idxInCluster[j]) / cluster_size[c];
          x1 += (abs(Spacing) + .2) * 3 * max_dimension_total_ * c;
          x2 += (abs(Spacing) + .2) * 3 * max_dimension_total_ * c
                + Spacing * max_dimension_total_ * cos(angle);
          y2 += Spacing * max_dimension_total_ * sin(angle);
        } else if(DisplayMethod == 2) {
          z2 = Spacing;
          if(intermediateDiagrams_.size() == 2 and j == 0) {
            z2 = -Spacing;
          }
        }
        if(good_id > -1) {
          matchingPoints->InsertNextPoint(x1, y1, z1);
          matchingPoints->InsertNextPoint(x2, y2, z2);
          matchingMesh->InsertNextCell(VTK_LINE, 2, ids);
          idOfDiagramMatching->InsertTuple1(count, j);
          idOfCluster->InsertTuple1(count, inv_clustering_[j]);
          cost->InsertTuple1(count, std::get<2>(m));
          idOfDiagramMatchingPoint->InsertTuple1(2 * count, j);
          idOfDiagramMatchingPoint->InsertTuple1(2 * count + 1, j);
          idOfPoint->InsertTuple1(2 * count, good_id);
          idOfPoint->InsertTuple1(2 * count + 1, bidder_id);

          const ttk::SimplexId type = std::get<5>(t2);
          switch(type) {
            case 0:
              pairType->InsertTuple1(count, 0);
              break;

            case 1:
              pairType->InsertTuple1(count, 1);
              break;

            case 2:
              pairType->InsertTuple1(count, 2);
              break;
            default:
              pairType->InsertTuple1(count, 0);
          }
          count++;
        }
      }
    }
  }
  if(NumberOfClusters == 1 and intermediateDiagrams_.size() == 2) {
    for(int i = 0; i < count; i++) {
      matchingCount->InsertTuple1(i, matchings_count[count_to_good[i]]);
    }
  }

  matchingMesh->SetPoints(matchingPoints);
  matchingMesh->GetPointData()->AddArray(idOfDiagramMatchingPoint);
  matchingMesh->GetPointData()->AddArray(idOfPoint);
  matchingMesh->GetCellData()->AddArray(idOfDiagramMatching);
  matchingMesh->GetCellData()->AddArray(idOfCluster);
  matchingMesh->GetCellData()->AddArray(pairType);
  matchingMesh->GetCellData()->AddArray(cost);
  if(NumberOfClusters == 1 and intermediateDiagrams_.size() == 2) {
    matchingMesh->GetCellData()->AddArray(matchingCount);
  }

  return matchingMesh;
}
