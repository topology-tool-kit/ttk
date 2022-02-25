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
  if(port == 0 || port == 1 || port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  } else {
    return 0;
  }
  return 1;
}

void ttkPersistenceDiagramClustering::Modified() {
  needUpdate_ = true;
  ttkAlgorithm::Modified();
}

int ttkPersistenceDiagramClustering::RequestData(
  vtkInformation *ttkNotUsed(request),
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
  auto output_matchings = vtkMultiBlockDataSet::GetData(outputVector, 2);

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

  outputClusteredDiagrams(output_clusters, input, this->inv_clustering_,
                          this->DisplayMethod, this->Spacing,
                          this->max_dimension_total_);
  outputCentroids(output_centroids, this->final_centroids_, this->DisplayMethod,
                  this->Spacing, this->max_dimension_total_);
  outputMatchings(
    output_matchings, this->NumberOfClusters, this->intermediateDiagrams_,
    this->all_matchings_, this->final_centroids_, this->inv_clustering_,
    this->DisplayMethod, this->Spacing, this->max_dimension_total_);

  // forward input diagrams FieldData to output_clusters blocks
  for(size_t i = 0; i < input.size(); ++i) {
    const auto diag{input[i]};
    const auto block{
      vtkUnstructuredGrid::SafeDownCast(output_clusters->GetBlock(i))};
    if(block != nullptr && block->GetFieldData() != nullptr
       && diag->GetFieldData() != nullptr) {
      block->GetFieldData()->ShallowCopy(diag->GetFieldData());
      // add clusterId to FieldData
      vtkNew<vtkIntArray> cid{};
      cid->SetName("ClusterId");
      cid->SetNumberOfTuples(1);
      cid->SetTuple1(0, this->inv_clustering_[i]);
      block->GetFieldData()->AddArray(cid);
    }
  }

  // add distance results to output_matchings FieldData
  vtkNew<vtkDoubleArray> minSad{};
  minSad->SetName("MinSaddleCost");
  minSad->SetNumberOfTuples(1);
  minSad->SetTuple1(0, this->distances[0]);

  vtkNew<vtkDoubleArray> sadSad{};
  sadSad->SetName("SaddleSaddleCost");
  sadSad->SetNumberOfTuples(1);
  sadSad->SetTuple1(0, this->distances[1]);

  vtkNew<vtkDoubleArray> sadMax{};
  sadMax->SetName("SaddleMaxCost");
  sadMax->SetNumberOfTuples(1);
  sadMax->SetTuple1(0, this->distances[2]);

  vtkNew<vtkDoubleArray> wass{};
  wass->SetName("WassersteinDistance");
  wass->SetNumberOfTuples(1);
  wass->SetTuple1(
    0, std::accumulate(this->distances.begin(), this->distances.end(), 0.0));

  for(size_t i = 0; i < output_matchings->GetNumberOfBlocks(); ++i) {
    const auto block{
      vtkUnstructuredGrid::SafeDownCast(output_matchings->GetBlock(i))};
    if(block != nullptr && block->GetFieldData() != nullptr) {
      block->GetFieldData()->AddArray(minSad);
      block->GetFieldData()->AddArray(sadSad);
      block->GetFieldData()->AddArray(sadMax);
      block->GetFieldData()->AddArray(wass);
    }
  }

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

  // compact pairIds in [0, nPairs - 1] (diagonal excepted)
  for(int i = 0; i < nPairs; i++) {
    if(pairId->GetTuple1(i) != -1) {
      pairId->SetTuple1(i, i);
    }
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

  vtkNew<ttkSimplexIdTypeArray> vsf{};
  vsf->SetName(ttk::VertexScalarFieldName);
  vsf->SetNumberOfTuples(nPoints);
  output->GetPointData()->AddArray(vsf);

  // cell data
  vtkNew<vtkIntArray> pairId{};
  pairId->SetName("PairIdentifier");
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
    const auto birtVertId{std::get<0>(pair)};
    const auto deathVertId{std::get<2>(pair)};
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
    vsf->SetTuple1(2 * j + 0, birtVertId);
    vsf->SetTuple1(2 * j + 1, deathVertId);

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
  pairId->SetTuple1(diagram.size(), -1);
  pairType->SetTuple1(diagram.size(), -1);
  // use twice the max persistence of all input diagrams...
  pairPers->SetTuple1(diagram.size(), 2.0 * max_persistence);
}

void ttkPersistenceDiagramClustering::outputClusteredDiagrams(
  vtkMultiBlockDataSet *output,
  const std::vector<vtkUnstructuredGrid *> &diags,
  const std::vector<int> &inv_clustering,
  const DISPLAY dm,
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
    vtu->ShallowCopy(diags[i]);

    vtkNew<vtkIntArray> clusterId{};
    clusterId->SetName("ClusterID");
    clusterId->SetNumberOfComponents(1);
    clusterId->SetNumberOfTuples(vtu->GetNumberOfPoints());
    clusterId->Fill(inv_clustering[i]);
    vtu->GetPointData()->AddArray(clusterId);

    // add Persistence data array on vertices
    vtkNew<vtkDoubleArray> pointPers{};
    pointPers->SetName("Persistence");
    pointPers->SetNumberOfTuples(vtu->GetNumberOfPoints());
    vtu->GetPointData()->AddArray(pointPers);

    // diagonal uses two existing points
    for(int j = 0; j < vtu->GetNumberOfCells() - 1; ++j) {
      const auto persArray = vtu->GetCellData()->GetArray("Persistence");
      const auto pers = persArray->GetTuple1(j);
      pointPers->SetTuple1(2 * j + 0, pers);
      pointPers->SetTuple1(2 * j + 1, pers);
    }

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
  const std::vector<diagramType> &final_centroids,
  const DISPLAY dm,
  const double spacing,
  const double max_persistence) const {

  for(size_t i = 0; i < final_centroids.size(); ++i) {
    vtkNew<vtkUnstructuredGrid> vtu{};
    this->diagramToVTU(vtu, final_centroids[i], i, this->max_dimension_total_);

    vtkNew<vtkIntArray> cid{};
    cid->SetName("ClusterId");
    cid->SetNumberOfTuples(1);
    cid->SetTuple1(0, i);
    vtu->GetFieldData()->AddArray(cid);

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

void ttkPersistenceDiagramClustering::outputMatchings(
  vtkMultiBlockDataSet *output,
  const size_t nClusters,
  const std::vector<diagramType> &diags,
  const std::vector<std::vector<std::vector<matchingType>>>
    &matchingsPerCluster,
  const std::vector<diagramType> &centroids,
  const std::vector<int> &inv_clustering,
  const DISPLAY dm,
  const double spacing,
  const double max_persistence) const {

  // index of diagram in its cluster
  std::vector<int> diagIdInClust{};
  // number of diagrams per cluster
  std::vector<int> clustSize{};

  // prep work for displaying diagrams as cluster stars
  if(dm == DISPLAY::STARS) {
    clustSize.resize(nClusters, 0);
    diagIdInClust.resize(diags.size());
    for(size_t i = 0; i < inv_clustering.size(); ++i) {
      auto &diagsInClust = clustSize[inv_clustering[i]];
      diagIdInClust[i] = diagsInClust;
      diagsInClust++;
    }
  }

  // count the number of bidders per centroid pair
  // (when with only 1 cluster and 2 diagrams)
  std::vector<int> matchings_count(centroids[0].size());
  std::vector<int> count_to_good{};

  for(size_t i = 0; i < diags.size(); ++i) {
    const auto cid = inv_clustering[i];
    const auto &diag{diags[i]};
    const auto &matchings{matchingsPerCluster[cid][i]};

    vtkNew<vtkUnstructuredGrid> matchingsGrid{};

    const auto nCells{matchings.size()};
    const auto nPoints{2 * matchings.size()};

    vtkNew<vtkPoints> points{};
    points->SetNumberOfPoints(nPoints);
    matchingsGrid->SetPoints(points);

    // point data
    vtkNew<vtkIntArray> diagIdVerts{};
    diagIdVerts->SetName("DiagramID");
    diagIdVerts->SetNumberOfTuples(nPoints);
    matchingsGrid->GetPointData()->AddArray(diagIdVerts);

    vtkNew<vtkIntArray> pointId{};
    pointId->SetName("PointID");
    pointId->SetNumberOfTuples(nPoints);
    matchingsGrid->GetPointData()->AddArray(pointId);

    // cell data
    vtkNew<vtkIntArray> diagIdCells{};
    diagIdCells->SetName("DiagramID");
    diagIdCells->SetNumberOfTuples(nCells);
    matchingsGrid->GetCellData()->AddArray(diagIdCells);

    vtkNew<vtkIntArray> clusterId{};
    clusterId->SetName("ClusterID");
    clusterId->SetNumberOfTuples(nCells);
    clusterId->Fill(cid);
    matchingsGrid->GetCellData()->AddArray(clusterId);

    vtkNew<vtkDoubleArray> matchCost{};
    matchCost->SetName("Cost");
    matchCost->SetNumberOfTuples(nCells);
    matchingsGrid->GetCellData()->AddArray(matchCost);

    vtkNew<vtkIntArray> pairType{};
    pairType->SetName("PairType");
    pairType->SetNumberOfTuples(nCells);
    matchingsGrid->GetCellData()->AddArray(pairType);

    vtkNew<vtkIntArray> isDiagonal{};
    isDiagonal->SetName("IsDiagonal");
    isDiagonal->SetNumberOfTuples(nCells);
    matchingsGrid->GetCellData()->AddArray(isDiagonal);

    for(size_t j = 0; j < matchings.size(); ++j) {
      const auto &m{matchings[j]};
      const auto bidderId{std::get<0>(m)};
      const auto goodId{std::get<1>(m)};

      // avoid out-of-bound accesses
      if(goodId >= static_cast<ttk::SimplexId>(centroids[cid].size())
         || bidderId >= static_cast<ttk::SimplexId>(diag.size())) {
        this->printWrn("Out-of-bounds access averted");
        continue;
      }

      if(nClusters == 1) {
        matchings_count[goodId] += 1;
        count_to_good.push_back(goodId);
      }

      const auto &p0{centroids[cid][goodId]};
      std::array<double, 3> coords0{std::get<6>(p0), std::get<10>(p0), 0};

      std::array<double, 3> coords1{};

      if(bidderId >= 0) {
        const auto &p1{diag[bidderId]};
        coords1[0] = std::get<6>(p1);
        coords1[1] = std::get<10>(p1);
        coords1[2] = 0;
        isDiagonal->SetTuple1(j, 0);
      } else {
        double diagonal_projection = (std::get<6>(p0) + std::get<10>(p0)) / 2;
        coords1[0] = diagonal_projection;
        coords1[1] = diagonal_projection;
        coords1[2] = 0;
        isDiagonal->SetTuple1(j, 1);
      }

      if(dm == DISPLAY::STARS && spacing > 0) {
        const auto angle = 2.0 * M_PI * static_cast<double>(diagIdInClust[i])
                           / static_cast<double>(clustSize[cid]);
        const auto shift
          = 3.0 * (std::abs(spacing) + 0.2) * max_persistence * cid;
        coords0[0] += shift;
        coords1[0] += shift + spacing * max_persistence * std::cos(angle);
        coords1[1] += spacing * max_persistence * std::sin(angle);

      } else if(dm == DISPLAY::MATCHINGS) {
        coords1[2] = (diags.size() == 2 && i == 0) ? -spacing : spacing;
      }

      points->SetPoint(2 * j + 0, coords0.data());
      points->SetPoint(2 * j + 1, coords1.data());
      std::array<vtkIdType, 2> ids{
        2 * static_cast<vtkIdType>(j) + 0,
        2 * static_cast<vtkIdType>(j) + 1,
      };
      matchingsGrid->InsertNextCell(VTK_LINE, 2, ids.data());

      diagIdCells->SetTuple1(j, i);
      matchCost->SetTuple1(j, std::get<2>(m));
      diagIdVerts->SetTuple1(2 * j + 0, i);
      diagIdVerts->SetTuple1(2 * j + 1, i);
      pointId->SetTuple1(2 * j + 0, goodId);
      pointId->SetTuple1(2 * j + 1, bidderId);
      pairType->SetTuple1(j, std::get<5>(p0));
    }

    output->SetBlock(i, matchingsGrid);
  }

  // add matchings number
  if(nClusters == 1 && diags.size() == 2) {
    size_t nPrevCells{};
    for(size_t i = 0; i < diags.size(); ++i) {
      vtkNew<vtkIntArray> matchNumber{};
      matchNumber->SetName("MatchNumber");
      const auto matchings
        = vtkUnstructuredGrid::SafeDownCast(output->GetBlock(i));
      const auto nCells = matchings->GetNumberOfCells();
      matchNumber->SetNumberOfTuples(nCells);
      for(int j = 0; j < nCells; j++) {
        const auto goodId = count_to_good[j + (i == 0 ? 0 : nPrevCells)];
        matchNumber->SetTuple1(j, matchings_count[goodId]);
      }
      matchings->GetCellData()->AddArray(matchNumber);
      nPrevCells = nCells;
    }
  }
}
