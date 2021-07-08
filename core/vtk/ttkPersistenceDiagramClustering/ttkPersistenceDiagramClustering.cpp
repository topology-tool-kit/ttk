#include <ttkMacros.h>
#include <ttkPersistenceDiagramClustering.h>
#include <ttkUtils.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>

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

void ttkPersistenceDiagramClustering::diagramToVTU(
  vtkUnstructuredGrid *output,
  const diagramType &diagram,
  const int cid,
  const double max_persistence) const {

  const auto nPoints = 2 * diagram.size();
  if(nPoints == 0) {
    this->printWrn("Diagram with no points");
    return;
  }

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(nPoints);
  output->SetPoints(points);

  // point data
  vtkNew<vtkIntArray> critType{};
  critType->SetName("CriticalType");
  critType->SetNumberOfTuples(nPoints);
  output->GetPointData()->AddArray(critType);

  vtkNew<vtkIntArray> clusterId{};
  clusterId->SetName("ClusterID");
  clusterId->SetNumberOfComponents(1);
  clusterId->SetNumberOfTuples(nPoints);
  clusterId->Fill(cid);
  output->GetPointData()->AddArray(clusterId);

  vtkNew<vtkFloatArray> coords{};
  coords->SetNumberOfComponents(3);
  coords->SetName("Coordinates");
  coords->SetNumberOfTuples(nPoints);
  output->GetPointData()->AddArray(coords);

  vtkNew<vtkDoubleArray> pointPers{};
  pointPers->SetName("Persistence");
  pointPers->SetNumberOfTuples(nPoints);
  output->GetPointData()->AddArray(pointPers);

  // cell data
  vtkNew<vtkIntArray> pairId{};
  pairId->SetName("PairID");
  pairId->SetNumberOfTuples(diagram.size() + 1);
  output->GetCellData()->AddArray(pairId);

  vtkNew<vtkIntArray> pairType{};
  pairType->SetName("PairType");
  pairType->SetNumberOfTuples(diagram.size() + 1);
  output->GetCellData()->AddArray(pairType);

  vtkNew<vtkDoubleArray> pairPers{};
  pairPers->SetName("Persistence");
  pairPers->SetNumberOfTuples(diagram.size() + 1);
  output->GetCellData()->AddArray(pairPers);

  for(size_t j = 0; j < diagram.size(); ++j) {
    const auto &pair{diagram[j]};
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

    points->SetPoint(2 * j + 0, birth, birth, 0);
    points->SetPoint(2 * j + 1, birth, death, 0);

    const std::array<vtkIdType, 2> ids{
      2 * static_cast<vtkIdType>(j) + 0,
      2 * static_cast<vtkIdType>(j) + 1,
    };
    output->InsertNextCell(VTK_LINE, 2, ids.data());
  }

  // add diagonal
  const auto minmax_birth = std::minmax_element(
    diagram.begin(), diagram.end(), [](const pairTuple &a, const pairTuple &b) {
      return std::get<6>(a) < std::get<6>(b);
    });
  const std::array<vtkIdType, 2> ids{
    2 * (minmax_birth.first - diagram.begin()),
    2 * (minmax_birth.second - diagram.begin()),
  };
  output->InsertNextCell(VTK_LINE, 2, ids.data());
  pairId->SetTuple1(diagram.size(), diagram.size());
  pairType->SetTuple1(diagram.size(), -1);
  // use the max persistence of all input diagrams...
  pairPers->SetTuple1(diagram.size(), max_persistence);
}

void ttkPersistenceDiagramClustering::outputClusteredDiagrams(
  vtkMultiBlockDataSet *output,
  const std::vector<diagramType> &diags,
  const std::vector<int> &inv_clustering,
  const enum DISPLAY dm,
  const double spacing,
  const double max_persistence) const {

  // index of diagram in its cluster
  std::vector<int> diagIdInClust{};
  // number of diagrams per cluster
  std::vector<int> clustSize{};

  // prep work for displaying diagrams as cluster stars
  if(dm == DISPLAY::STARS) {
    // total number of clusters
    const auto nClusters
      = 1 + *std::max_element(inv_clustering.begin(), inv_clustering.end());
    clustSize.resize(nClusters, 0);
    diagIdInClust.resize(diags.size());
    for(size_t i = 0; i < inv_clustering.size(); ++i) {
      auto &diagsInClust = clustSize[inv_clustering[i]];
      diagIdInClust[i] = diagsInClust;
      diagsInClust++;
    }
  }

  output->SetNumberOfBlocks(diags.size());

  for(size_t i = 0; i < diags.size(); ++i) {
    vtkNew<vtkUnstructuredGrid> vtu{};
    this->diagramToVTU(
      vtu, diags[i], this->inv_clustering_[i], this->max_dimension_total_);

    if(dm == DISPLAY::MATCHINGS && spacing > 0) {
      // translate diagrams along the Z axis
      vtkNew<vtkTransform> tr{};
      tr->Translate(0, 0, i == 0 ? -spacing : spacing);
      vtkNew<vtkTransformFilter> trf{};
      trf->SetTransform(tr);
      trf->SetInputData(vtu);
      trf->Update();
      output->SetBlock(i, trf->GetOutputDataObject(0));
    } else if(dm == DISPLAY::STARS && spacing > 0) {
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
      trf->SetInputData(vtu);
      trf->Update();
      output->SetBlock(i, trf->GetOutputDataObject(0));

    } else {
      // add diagram to output multi-block dataset
      output->SetBlock(i, vtu);
    }
  }
}

void ttkPersistenceDiagramClustering::outputCentroids(
  vtkMultiBlockDataSet *output,
  std::vector<diagramType> &final_centroids,
  const DISPLAY dm,
  const double spacing,
  const double max_persistence) const {

  for(size_t i = 0; i < final_centroids.size(); ++i) {
    vtkNew<vtkUnstructuredGrid> vtu{};
    this->diagramToVTU(vtu, final_centroids[i], i, this->max_dimension_total_);

    if(dm == DISPLAY::STARS && spacing > 0) {
      // shift centroid along the X axis
      vtkNew<vtkTransform> tr{};
      tr->Translate(3.0 * (spacing + 0.2) * max_persistence * i, 0, 0);
      vtkNew<vtkTransformFilter> trf{};
      trf->SetTransform(tr);
      trf->SetInputData(vtu);
      trf->Update();
      output->SetBlock(i, trf->GetOutputDataObject(0));

    } else {
      // add centroid to output multi-block dataset
      output->SetBlock(i, vtu);
    }
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

    // store the persistence of every min-max global pair
    std::vector<double> max_persistences(numInputs);

    for(int i = 0; i < numInputs; i++) {
      this->VTUToDiagram(this->intermediateDiagrams_[i], input[i]);
      max_persistences[i] = std::get<4>(intermediateDiagrams_[i][0]);
    }

    this->max_dimension_total_
      = *std::max_element(max_persistences.begin(), max_persistences.end());

    if(this->Method == METHOD::PROGRESSIVE) {

      if(!UseInterruptible) {
        TimeLimit = 999999999;
      }

      inv_clustering_ = this->execute<double>(
        intermediateDiagrams_, final_centroids_, all_matchings_);
      needUpdate_ = false;

    } else if(this->Method == METHOD::AUCTION) {

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

  outputClusteredDiagrams(output_clusters, this->intermediateDiagrams_,
                          this->inv_clustering_, this->DisplayMethod,
                          this->Spacing, this->max_dimension_total_);
  outputCentroids(output_centroids, this->final_centroids_, this->DisplayMethod,
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

void ttkPersistenceDiagramClustering::VTUToDiagram(
  diagramType &diagram, vtkUnstructuredGrid *vtu) const {

  const auto pd = vtu->GetPointData();
  const auto cd = vtu->GetCellData();

  if(pd == nullptr) {
    this->printErr("VTU diagram with NULL Point Data");
    return;
  }
  if(cd == nullptr) {
    this->printErr("VTU diagram with NULL Cell Data");
    return;
  }

  // cell data
  const auto pairId = vtkIntArray::SafeDownCast(cd->GetArray("PairIdentifier"));
  const auto pairType = vtkIntArray::SafeDownCast(cd->GetArray("PairType"));
  const auto pairPers
    = vtkDoubleArray::SafeDownCast(cd->GetArray("Persistence"));

  // point data
  const auto vertexId
    = vtkIntArray::SafeDownCast(pd->GetArray(ttk::VertexScalarFieldName));
  const auto critType = vtkIntArray::SafeDownCast(pd->GetArray("CriticalType"));
  const auto birthScalars = vtkDoubleArray::SafeDownCast(pd->GetArray("Birth"));
  const auto deathScalars = vtkDoubleArray::SafeDownCast(pd->GetArray("Death"));
  const auto coords = vtkFloatArray::SafeDownCast(pd->GetArray("Coordinates"));

  const auto points = vtu->GetPoints();

  const bool embed = birthScalars != nullptr && deathScalars != nullptr;

  if(!embed && coords == nullptr) {
    this->printErr("Missing coordinates array on non-embedded diagram");
    return;
  }

  int nPairs = pairId->GetNumberOfTuples();

  // FIX : no more missed pairs
  for(int pair_index = 0; pair_index < nPairs; pair_index++) {
    const float index_of_pair = pair_index;
    if(pairId->GetTuple1(pair_index) != -1)
      pairId->SetTuple(pair_index, &index_of_pair);
  }

  // skip diagram diagonal if present (assuming it's the last pair in the
  // diagram)
  if(pairId->GetTuple1(nPairs - 1) == -1)
    nPairs -= 1;

  if(nPairs < 1 || vertexId == nullptr || pairId == nullptr
     || critType == nullptr || pairPers == nullptr || pairType == nullptr
     || points == nullptr) {
    this->printErr("Either no pairs in diagram or some array is NULL");
    return;
  }

  if(NumberOfClusters == 1) {
    diagram.resize(nPairs);
  } else {
    diagram.resize(nPairs + 1);
  }

  // count the number of pairs whose index is >= nPairs
  int nbNonCompact = 0;

  // skip diagonal cell (corresponding points already dealt with)
  for(int i = 0; i < nPairs; ++i) {

    const int v0 = vertexId->GetValue(2 * i);
    const int v1 = vertexId->GetValue(2 * i + 1);
    const int ct0 = critType->GetValue(2 * i);
    const int ct1 = critType->GetValue(2 * i + 1);

    const int pId = pairId->GetValue(i);
    const int pType = pairType->GetValue(i);
    const double pers = pairPers->GetValue(i);

    std::array<double, 3> coordsBirth{}, coordsDeath{};
    double birth, death;

    if(embed) {
      points->GetPoint(2 * i + 0, coordsBirth.data());
      points->GetPoint(2 * i + 1, coordsDeath.data());
      birth = birthScalars->GetValue(2 * i + 0);
      death = deathScalars->GetValue(2 * i + 1);
    } else {
      coords->GetTuple(2 * i + 0, coordsBirth.data());
      coords->GetTuple(2 * i + 1, coordsDeath.data());
      birth = points->GetPoint(2 * i + 0)[0];
      death = points->GetPoint(2 * i + 1)[1];
    }

    if(pId != -1 && pId < nPairs) {

      if(pId == 0) {
        // deal with the global min-max pair separately
        // (what do we do with other infinite pairs?)

        if(NumberOfClusters == 1) {
          diagram[0] = std::make_tuple(
            v0, CriticalType::Local_minimum, v1, CriticalType::Local_maximum,
            pers, pType, birth, coordsBirth[0], coordsBirth[1], coordsBirth[2],
            death, coordsDeath[0], coordsDeath[1], coordsDeath[2]);
        } else {
          // duplicate the global min-max pair into two: one min-saddle pair and
          // one saddle-max pair
          diagram[0] = std::make_tuple(
            v0, CriticalType::Local_minimum, v1, CriticalType::Saddle1, pers,
            pType, birth, coordsBirth[0], coordsBirth[1], coordsBirth[2], death,
            coordsDeath[0], coordsDeath[1], coordsDeath[2]);
          // store the saddle max pair at the vector end
          diagram[nPairs] = std::make_tuple(
            v0, CriticalType::Saddle1, v1, CriticalType::Local_maximum, pers,
            pType, birth, coordsBirth[0], coordsBirth[1], coordsBirth[2], death,
            coordsDeath[0], coordsDeath[1], coordsDeath[2]);
        }

      } else {
        // all other pairs
        diagram[pId] = std::make_tuple(
          v0, static_cast<CriticalType>(ct0), v1,
          static_cast<CriticalType>(ct1), pers, pType, birth, coordsBirth[0],
          coordsBirth[1], coordsBirth[2], death, coordsDeath[0], coordsDeath[1],
          coordsDeath[2]);
      }
    }

    if(pId >= nPairs) {
      nbNonCompact++;
    }
  }

  if(nbNonCompact > 0) {
    this->printWrn("Missed " + std::to_string(nbNonCompact)
                   + " pairs due to non-compactness.");
  }
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
  if(DisplayMethod == DISPLAY::STARS && Spacing > 0) {
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

      if(DisplayMethod == DISPLAY::STARS && Spacing > 0) {
        double angle
          = 2 * 3.1415926 * (double)(idxInCluster[j]) / cluster_size[c];
        x1 += (abs(Spacing) + .2) * 3 * max_dimension_total_ * c;
        x2 += (abs(Spacing) + .2) * 3 * max_dimension_total_ * c
              + Spacing * max_dimension_total_ * cos(angle);
        y2 += Spacing * max_dimension_total_ * sin(angle);
      } else if(DisplayMethod == DISPLAY::MATCHINGS) {
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
