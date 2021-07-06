#include <ttkMacros.h>
#include <ttkPersistenceDiagramClustering.h>
#include <ttkUtils.h>
#include <vtkFieldData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

using namespace ttk;

vtkStandardNewMacro(ttkPersistenceDiagramClustering);

ttkPersistenceDiagramClustering::ttkPersistenceDiagramClustering() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
}

int ttkPersistenceDiagramClustering::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  } else {
    return 0;
  }
  return 1;
}

int ttkPersistenceDiagramClustering::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  } else if(port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  } else {
    return 0;
  }
  return 1;
}

void ttkPersistenceDiagramClustering::Modified() {
  needUpdate_ = true;
  ttkAlgorithm::Modified();
}

void outputClusteredDiagrams(
  vtkMultiBlockDataSet *output,
  const std::vector<vtkUnstructuredGrid *> &inputDiags,
  const std::vector<int> &inv_clustering,
  const enum DisplayMethodType dm,
  const double spacing,
  const double max_persistence) {

  // index of diagram in its cluster
  std::vector<int> diagIdInClust{};
  // number of diagrams per cluster
  std::vector<int> clustSize{};

  // prep work for displaying diagrams as cluster stars
  if(dm == DisplayMethodType::STARS) {
    // total number of clusters
    const auto nClusters
      = 1 + *std::max_element(inv_clustering.begin(), inv_clustering.end());
    clustSize.resize(nClusters, 0);
    diagIdInClust.resize(inputDiags.size());
    for(size_t i = 0; i < inv_clustering.size(); ++i) {
      auto &diagsInClust = clustSize[inv_clustering[i]];
      diagIdInClust[i] = diagsInClust;
      diagsInClust++;
    }
  }

  output->SetNumberOfBlocks(inputDiags.size());
  for(size_t i = 0; i < inputDiags.size(); ++i) {
    const auto diag = inputDiags[i];
    // add clustering info to diagram
    vtkNew<vtkIntArray> clusterId{};
    clusterId->SetName("ClusterID");
    clusterId->SetNumberOfComponents(1);
    clusterId->SetNumberOfTuples(diag->GetNumberOfPoints());
    clusterId->Fill(inv_clustering[i]);
    diag->GetPointData()->AddArray(clusterId);

    if(dm == DisplayMethodType::MATCHINGS && spacing > 0) {
      // translate diagrams along the Z axis
      vtkNew<vtkTransform> tr{};
      tr->Translate(0, 0, i == 0 ? -spacing : spacing);
      vtkNew<vtkTransformFilter> trf{};
      trf->SetTransform(tr);
      trf->SetInputData(diag);
      trf->Update();
      output->SetBlock(i, trf->GetOutputDataObject(0));
    } else if(dm == DisplayMethodType::STARS && spacing > 0) {
      const auto c = inv_clustering[i];
      const auto angle = 2.0 * M_PI * static_cast<double>(diagIdInClust[i])
                         / static_cast<double>(clustSize[c]);
      // translate diagrams in the XY plane
      vtkNew<vtkTransform> tr{};
      tr->Translate(3.0 * (spacing + 0.2) * max_persistence * c
                      + spacing * max_persistence * std::cos(angle) + 0.2,
                    spacing * max_persistence * std::sin(angle), 0);
      vtkNew<vtkTransformFilter> trf{};
      trf->SetTransform(tr);
      trf->SetInputData(diag);
      trf->Update();
      output->SetBlock(i, trf->GetOutputDataObject(0));

    } else {
      // add diagram to output multi-block dataset
      output->SetBlock(i, diag);
    }
  }
}

void outputCentroids(
  vtkMultiBlockDataSet *output,
  std::vector<ttkPersistenceDiagramClustering::diagramType> &final_centroids,
  const DisplayMethodType dm,
  const double spacing,
  const double max_persistence) {

  for(size_t i = 0; i < final_centroids.size(); ++i) {
    const auto &bary{final_centroids[i]};

    const auto nPoints = 2 * bary.size();
    vtkNew<vtkPoints> points{};
    points->SetNumberOfPoints(nPoints);

    vtkNew<vtkUnstructuredGrid> diagram{};
    diagram->SetPoints(points);

    // point data
    vtkNew<vtkIntArray> critType{};
    critType->SetName("CriticalType");
    critType->SetNumberOfTuples(nPoints);
    diagram->GetPointData()->AddArray(critType);

    vtkNew<vtkIntArray> clusterId{};
    clusterId->SetName("ClusterID");
    clusterId->SetNumberOfTuples(nPoints);
    clusterId->Fill(i);
    diagram->GetPointData()->AddArray(clusterId);

    vtkNew<vtkFloatArray> coords{};
    coords->SetNumberOfComponents(3);
    coords->SetName("Coordinates");
    coords->SetNumberOfTuples(nPoints);
    diagram->GetPointData()->AddArray(coords);

    vtkNew<vtkDoubleArray> pointPers{};
    pointPers->SetName("Persistence");
    pointPers->SetNumberOfTuples(nPoints);
    diagram->GetPointData()->AddArray(pointPers);

    // cell data
    vtkNew<vtkIntArray> pairId{};
    pairId->SetName("PairID");
    pairId->SetNumberOfTuples(bary.size());
    diagram->GetCellData()->AddArray(pairId);

    vtkNew<vtkIntArray> pairType{};
    pairType->SetName("PairType");
    pairType->SetNumberOfTuples(bary.size());
    diagram->GetCellData()->AddArray(pairType);

    vtkNew<vtkDoubleArray> pairPers{};
    pairPers->SetName("Persistence");
    pairPers->SetNumberOfTuples(bary.size());
    diagram->GetCellData()->AddArray(pairPers);

    for(size_t j = 0; j < bary.size(); ++j) {
      const auto &pair{bary[j]};
      const auto birth{std::get<6>(pair)};
      const auto death{std::get<10>(pair)};
      const auto pType{std::get<5>(pair)};
      const auto birthType{std::get<1>(pair)};
      const auto deathType{std::get<3>(pair)};
      std::array<float, 3> coordsBirth{
        std::get<7>(pair), std::get<8>(pair), std::get<9>(pair)};
      std::array<float, 3> coordsDeath{
        std::get<11>(pair), std::get<12>(pair), std::get<13>(pair)};

      // cell data
      pairId->SetTuple1(j, j);
      pairType->SetTuple1(j, pType);
      pairPers->SetTuple1(j, death - birth);

      // point data
      coords->SetTuple(2 * j + 0, coordsBirth.data());
      coords->SetTuple(2 * j + 1, coordsDeath.data());
      pointPers->SetTuple1(2 * j + 0, death - birth);
      pointPers->SetTuple1(2 * j + 1, death - birth);
      critType->SetTuple1(2 * j + 0, static_cast<int>(birthType));
      critType->SetTuple1(2 * j + 1, static_cast<int>(deathType));

      auto birthShift{birth};
      if(dm == DisplayMethodType::STARS && spacing > 0) {
        // shift diagram along the X axis
        birthShift += 3.0 * (spacing + 0.2) * max_persistence * i;
      }
      points->SetPoint(2 * j + 0, birthShift, birth, 0);
      points->SetPoint(2 * j + 1, birthShift, death, 0);

      const std::array<vtkIdType, 2> ids{
        2 * static_cast<vtkIdType>(j) + 0,
        2 * static_cast<vtkIdType>(j) + 1,
      };
      diagram->InsertNextCell(VTK_LINE, 2, ids.data());
    }

    // attach barycenter diagram to multi-block output
    output->SetBlock(i, diagram);
  }
}

int ttkPersistenceDiagramClustering::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  Memory m;

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);

  // Flat storage for diagrams extracted from blocks
  std::vector<vtkUnstructuredGrid *> input;

  // Number of input diagrams
  int numInputs = 0;

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

  // Get output pointers
  auto output_clusters = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_centroids = vtkMultiBlockDataSet::GetData(outputVector, 1);
  auto output_matchings = vtkUnstructuredGrid::GetData(outputVector, 2);

  if(needUpdate_) {
    // clear data before computation
    intermediateDiagrams_ = {};
    all_matchings_ = {};
    final_centroids_ = {};

    intermediateDiagrams_.resize(numInputs);
    all_matchings_.resize(3);

    max_dimension_total_ = 0;
    for(int i = 0; i < numInputs; i++) {
      double max_dimension
        = getPersistenceDiagram(intermediateDiagrams_[i], input[i]);
      if(max_dimension_total_ < max_dimension) {
        max_dimension_total_ = max_dimension;
      }
    }

    if(Method == 0) {

      // Progressive approach
      if(!UseInterruptible) {
        TimeLimit = 999999999;
      }

      inv_clustering_ = this->execute<double>(
        intermediateDiagrams_, final_centroids_, all_matchings_);
      needUpdate_ = false;

    } else {

      // AUCTION APPROACH
      final_centroids_.resize(1);
      inv_clustering_.resize(numInputs);
      for(int i_input = 0; i_input < numInputs; i_input++) {
        inv_clustering_[i_input] = 0;
      }
      PersistenceDiagramBarycenter<double> pdBarycenter;

      const auto wassersteinMetric = std::to_string(WassersteinMetric);
      pdBarycenter.setWasserstein(wassersteinMetric);
      pdBarycenter.setMethod(2);
      pdBarycenter.setNumberOfInputs(numInputs);
      pdBarycenter.setTimeLimit(TimeLimit);
      pdBarycenter.setDeterministic(Deterministic);
      pdBarycenter.setUseProgressive(UseProgressive);
      pdBarycenter.setDebugLevel(debugLevel_);
      pdBarycenter.setThreadNumber(threadNumber_);
      pdBarycenter.setAlpha(Alpha);
      pdBarycenter.setLambda(Lambda);
      pdBarycenter.execute(
        intermediateDiagrams_, final_centroids_[0], all_matchings_);

      needUpdate_ = false;
    }
  }

  output_matchings->ShallowCopy(createMatchings());

  outputClusteredDiagrams(output_clusters, input, this->inv_clustering_,
                          static_cast<enum DisplayMethodType>(this->DisplayMethod),
                          this->Spacing, this->max_dimension_total_);
  outputCentroids(output_centroids, this->final_centroids_,
                  static_cast<enum DisplayMethodType>(this->DisplayMethod),
                  this->Spacing, this->max_dimension_total_);

  // collect input FieldData to annotate outputClusters
  auto fd = output_clusters->GetFieldData();
  bool hasStructure{false};
  for(const auto diag : input) {
    if(diag->GetFieldData() != nullptr) {
      // ensure fd has the same structure as the first non-null input
      // FieldData
      if(!hasStructure) {
        fd->CopyStructure(diag->GetFieldData());
        hasStructure = true;
      }
      // copy data
      fd->InsertNextTuple(0, diag->GetFieldData());
    }
  }

  // add clusterId to outputClusters FieldData
  vtkNew<vtkIntArray> cid{};
  cid->SetName("ClusterId");
  for(const auto c : this->inv_clustering_) {
    cid->InsertNextValue(c);
  }
  fd->AddArray(cid);
  return 1;
}

double ttkPersistenceDiagramClustering::getPersistenceDiagram(
  diagramType &diagram, vtkUnstructuredGrid *CTPersistenceDiagram_) {
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

  // If diagram has the diagonal (we assume it is last)
  if(*pairIdentifierScalars->GetTuple(pairingsSize - 1) == -1)
    pairingsSize -= 1;

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
  for(int i = 0; i < pairingsSize; ++i) {

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
        msg << "Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). " << std::endl;
        this->printWrn(msg.str());
      }
    }
  }

  if(nbNonCompact > 0) {
    {
      std::stringstream msg;
      msg << "Missed " << nbNonCompact << " pairs due to non-compactness."
          << std::endl;
      this->printWrn(msg.str());
    }
  }

  return max_dimension;
}

vtkNew<vtkUnstructuredGrid> ttkPersistenceDiagramClustering::createMatchings() {
  this->printMsg("Creating vtk matchings", debug::Priority::VERBOSE);
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
      const size_t bidder_id = std::get<0>(m);
      const size_t good_id = std::get<1>(m);

      // avoid out-of-bound accesses
      if(good_id >= matchings_count.size() || bidder_id >= diagram.size()) {
        continue;
      }

      if(NumberOfClusters == 1) {
        matchings_count[good_id] += 1;
        count_to_good.push_back(good_id);
      }

      pairTuple t1 = final_centroids_[c][good_id];
      double x1 = std::get<6>(t1);
      double y1 = std::get<10>(t1);
      double z1 = 0;

      pairTuple t2 = diagram[bidder_id];
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
