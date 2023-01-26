#include "BaseClass.h"
#include "DataTypes.h"
#include <string>
#include <ttkIntegralLines.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <ArrayLinkedList.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkUnstructuredGrid.h>

#include <array>
#include <numeric>

vtkStandardNewMacro(ttkIntegralLines);

ttkIntegralLines::ttkIntegralLines() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkIntegralLines::~ttkIntegralLines() = default;

int ttkIntegralLines::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");

  return 1;
}

int ttkIntegralLines::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");

  return 1;
}

template <typename triangulationType>
int ttkIntegralLines::getTrajectories(
  vtkDataSet *input,
  triangulationType *triangulation,
  std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
    &trajectories,
  std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
    &forkIdentifiers,
  std::vector<ttk::ArrayLinkedList<std::vector<double>, TABULAR_SIZE>>
    &distancesFromSeed,
  std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
    &seedIdentifiers,
#ifdef TTK_ENABLE_MPI
  std::vector<ttk::SimplexId> &globalVertexId,
  std::vector<ttk::SimplexId> &globalCellId,
#endif
  vtkUnstructuredGrid *output) {
  if(input == nullptr || output == nullptr
     || input->GetPointData() == nullptr) {
    this->printErr("Null pointers in getTrajectories parameters");
    return 0;
  }
  ttk::SimplexId threadNumber = trajectories.size();
  vtkNew<vtkUnstructuredGrid> ug{};
  vtkNew<vtkPoints> pts{};
  vtkNew<vtkDoubleArray> dist{};
  vtkNew<vtkIdTypeArray> identifier{};
#ifdef TTK_ENABLE_MPI
  vtkNew<vtkIdTypeArray> vtkEdgeIdentifiers{};
  vtkNew<vtkIdTypeArray> vtkVertexGlobalIdArray{};
  vtkNew<vtkIntArray> vtkVertexRankArray{};
  vtkNew<vtkIntArray> vtkEdgeRankArray{};
#endif
  vtkNew<vtkIdTypeArray> vtkForkIdentifiers{};
  vtkNew<vtkUnsignedCharArray> outputMaskField{};

  outputMaskField->SetNumberOfComponents(1);
  outputMaskField->SetName(ttk::MaskScalarFieldName);

  dist->SetNumberOfComponents(1);
  dist->SetName("DistanceFromSeed");
  identifier->SetNumberOfComponents(1);
  identifier->SetName("SeedIdentifier");
  vtkForkIdentifiers->SetNumberOfComponents(1);
  vtkForkIdentifiers->SetName("ForkIdentifiers");

#ifdef TTK_ENABLE_MPI
  vtkVertexGlobalIdArray->SetNumberOfComponents(1);
  vtkVertexGlobalIdArray->SetName("GlobalPointIds");
  vtkEdgeIdentifiers->SetNumberOfComponents(1);
  vtkEdgeIdentifiers->SetName("GlobalCellIds");
  vtkEdgeRankArray->SetNumberOfComponents(1);
  vtkEdgeRankArray->SetName("RankArray");
  vtkVertexRankArray->SetNumberOfComponents(1);
  vtkVertexRankArray->SetName("RankArray");
#endif
  const auto numberOfArrays = input->GetPointData()->GetNumberOfArrays();

  std::vector<vtkDataArray *> scalarArrays{};
  scalarArrays.reserve(numberOfArrays);
  for(int k = 0; k < numberOfArrays; ++k) {
    const auto a = input->GetPointData()->GetArray(k);
    if(a->GetNumberOfComponents() == 1) {
      // only keep scalar arrays
      scalarArrays.push_back(a);
    }
  }

  std::vector<vtkSmartPointer<vtkDataArray>> inputScalars(scalarArrays.size());
  for(size_t k = 0; k < scalarArrays.size(); ++k) {
    inputScalars[k]
      = vtkSmartPointer<vtkDataArray>::Take(scalarArrays[k]->NewInstance());
    inputScalars[k]->SetNumberOfComponents(1);
    inputScalars[k]->SetName(scalarArrays[k]->GetName());
  }
  std::array<float, 3> p;
  std::array<vtkIdType, 2> ids;
#ifdef TTK_ENABLE_MPI
  ttk::SimplexId vertexCounter = 0;
  ttk::SimplexId edgeCounter = 0;
#endif
  for(int thread = 0; thread < threadNumber; thread++) {
    std::list<std::array<std::vector<ttk::SimplexId>, TABULAR_SIZE>>::iterator
      trajectory
      = trajectories[thread].list.begin();
    std::list<std::array<ttk::SimplexId, TABULAR_SIZE>>::iterator forkIdentifier
      = forkIdentifiers[thread].list.begin();
    std::list<std::array<std::vector<double>, TABULAR_SIZE>>::iterator
      distanceFromSeed
      = distancesFromSeed[thread].list.begin();
    std::list<std::array<ttk::SimplexId, TABULAR_SIZE>>::iterator seedIdentifier
      = seedIdentifiers[thread].list.begin();
    while(trajectory != trajectories[thread].list.end()) {
      for(int i = 0; i < TABULAR_SIZE; i++) {
        if((*trajectory)[i].size() > 0) {
          ttk::SimplexId vertex = (*trajectory)[i].at(0);
          triangulation->getVertexPoint(vertex, p[0], p[1], p[2]);
          ids[0] = pts->InsertNextPoint(p.data());
          // distanceScalars
          dist->InsertNextTuple1((*distanceFromSeed)[i].at(0));
#ifdef TTK_ENABLE_MPI
          if(triangulation->getVertexRank((*trajectory)[i].back())
             == ttk::MPIrank_) {
            outputMaskField->InsertNextTuple1(0);
          } else {
            outputMaskField->InsertNextTuple1(1);
          }
          vtkVertexGlobalIdArray->InsertNextTuple1(
            globalVertexId.at(vertexCounter));
          vertexCounter++;
          vtkVertexRankArray->InsertNextTuple1(
            triangulation->getVertexRank(vertex));
#else
          outputMaskField->InsertNextTuple1(0);
#endif
          identifier->InsertNextTuple1((*seedIdentifier)[i]);
          vtkForkIdentifiers->InsertNextTuple1((*forkIdentifier)[i]);
          // inputScalars
          for(size_t k = 0; k < scalarArrays.size(); ++k) {
            inputScalars[k]->InsertNextTuple1(
              scalarArrays[k]->GetTuple1(vertex));
          }
          for(size_t j = 1; j < (*trajectory)[i].size(); ++j) {
            vertex = (*trajectory)[i].at(j);
#ifdef TTK_ENABLE_MPI
            vtkVertexGlobalIdArray->InsertNextTuple1(
              globalVertexId.at(vertexCounter));
            vertexCounter++;
            vtkEdgeIdentifiers->InsertNextTuple1(globalCellId.at(edgeCounter));
            edgeCounter++;
            vtkEdgeRankArray->InsertNextTuple1(
              triangulation->getVertexRank((*trajectory)[i].at(j - 1)));
            vtkVertexRankArray->InsertNextTuple1(
              triangulation->getVertexRank(vertex));
#endif
            outputMaskField->InsertNextTuple1(1);
            vtkForkIdentifiers->InsertNextTuple1((*forkIdentifier)[i]);
            triangulation->getVertexPoint(vertex, p[0], p[1], p[2]);
            ids[1] = pts->InsertNextPoint(p.data());
            // distanceScalars
            dist->InsertNextTuple1((*distanceFromSeed)[i].at(j));
            identifier->InsertNextTuple1((*seedIdentifier)[i]);
            // inputScalars
            for(unsigned int k = 0; k < scalarArrays.size(); ++k)
              inputScalars[k]->InsertNextTuple1(
                scalarArrays[k]->GetTuple1(vertex));
            ug->InsertNextCell(VTK_LINE, 2, ids.data());
            // iteration
            ids[0] = ids[1];
          }
#ifdef TTK_ENABLE_MPI
          if(triangulation->getVertexRank((*trajectory)[i].back())
             == ttk::MPIrank_) {
            outputMaskField->SetTuple1(
              outputMaskField->GetNumberOfComponents() - 1, 0);
          }
#endif
        } else {
          break;
        }
      }
      trajectory++;
      distanceFromSeed++;
      seedIdentifier++;
      forkIdentifier++;
    }
  }

  ug->SetPoints(pts);
  ug->GetPointData()->AddArray(dist);
  ug->GetPointData()->AddArray(identifier);
  ug->GetPointData()->AddArray(vtkForkIdentifiers);
  ug->GetPointData()->AddArray(outputMaskField);

  for(unsigned int k = 0; k < scalarArrays.size(); ++k) {
    ug->GetPointData()->AddArray(inputScalars[k]);
  }
#ifdef TTK_ENABLE_MPI
  ug->GetPointData()->AddArray(vtkVertexRankArray);
  ug->GetCellData()->AddArray(vtkEdgeRankArray);
  ug->GetCellData()->SetGlobalIds(vtkEdgeIdentifiers);
  ug->GetPointData()->SetGlobalIds(vtkVertexGlobalIdArray);
#endif
  output->ShallowCopy(ug);

  return 1;
}

#ifdef TTK_ENABLE_MPI
template <typename triangulationType>
int ttkIntegralLines::getGlobalIdentifiers(
  std::vector<ttk::SimplexId> &globalVertexId,
  std::vector<ttk::SimplexId> &globalCellId,
  std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
    &seedIdentifiers,
  std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
    &forkIdentifiers,
  std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
    &trajectories,
  std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
    &localVertexIdentifiers,
  triangulationType *triangulation) {
  ttk::SimplexId outputVertexNumber = 0;
  ttk::SimplexId outputCellNumber = 0;
  ttk::SimplexId realVertexNumber = 0;
  ttk::SimplexId realCellNumber = 0;
  ttk::SimplexId intervalSize;
  // Counts vertices and edges number (with and without ghosts)
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for reduction(+:outputVertexNumber,outputCellNumber,realCellNumber,realVertexNumber) schedule(static,1) private(intervalSize)
#endif
  for(int thread = 0; thread < threadNumber_; thread++) {
    std::list<std::array<std::vector<ttk::SimplexId>, TABULAR_SIZE>>::iterator
      trajectory
      = trajectories[thread].list.begin();
    while(trajectory != trajectories[thread].list.end()) {
      for(int i = 0; i < TABULAR_SIZE; i++) {
        if((*trajectory)[i].size() > 0) {
          intervalSize = static_cast<ttk::SimplexId>((*trajectory)[i].size());
          outputVertexNumber += intervalSize;
          outputCellNumber += intervalSize - 1;
          if(triangulation->getVertexRank((*trajectory)[i].at(0))
             != ttk::MPIrank_) {
            intervalSize--;
          }
          if((*trajectory)[i].size() > 1) {
            realCellNumber += intervalSize - 1;
            if(triangulation->getVertexRank((*trajectory)[i].back())
               != ttk::MPIrank_) {
              intervalSize--;
            }
          }
          realVertexNumber += intervalSize;
        } else {
          break;
        }
      }
      trajectory++;
    }
  }

  ttk::SimplexId vertIndex;
  ttk::SimplexId cellIndex;

  // Perform exclusive prefix sum to find local offset for vertices and
  // cells
  MPI_Exscan(&realVertexNumber, &vertIndex, 1, ttk::getMPIType(vertIndex),
             MPI_SUM, ttk::MPIcomm_);
  MPI_Exscan(&realCellNumber, &cellIndex, 1, ttk::getMPIType(cellIndex),
             MPI_SUM, ttk::MPIcomm_);

  // Rank 0 received garbage values, it is replaced by the correct offset
  // (always 0)
  if(ttk::MPIrank_ == 0) {
    vertIndex = 0;
    cellIndex = 0;
  }

  // Generate Global identifers and package ghost data for the process that owns
  // the ghost data
  ttk::SimplexId startCellId = 0;
  ttk::SimplexId startVertexId = 0;
  globalVertexId.resize(outputVertexNumber);
  globalCellId.resize(outputCellNumber);
  std::vector<std::vector<ttk::intgl::GhostElementsToSort>> unmatchedGhosts(
    neighborNumber_);
  std::vector<std::vector<ttk::SimplexId>> unmatchedGhostsVertexLocalId(
    neighborNumber_);
  std::vector<std::vector<ttk::SimplexId>> unmatchedGhostsEdgeLocalId(
    neighborNumber_);
  for(int thread = 0; thread < threadNumber_; thread++) {
    std::list<std::array<ttk::SimplexId, TABULAR_SIZE>>::iterator forkIdentifier
      = forkIdentifiers[thread].list.begin();
    std::list<std::array<std::vector<ttk::SimplexId>, TABULAR_SIZE>>::iterator
      localVertexIdentifier
      = localVertexIdentifiers[thread].list.begin();
    std::list<std::array<std::vector<ttk::SimplexId>, TABULAR_SIZE>>::iterator
      trajectory
      = trajectories[thread].list.begin();
    std::list<std::array<ttk::SimplexId, TABULAR_SIZE>>::iterator seedIdentifier
      = seedIdentifiers[thread].list.begin();
    while(trajectory != trajectories[thread].list.end()) {
      for(int i = 0; i < TABULAR_SIZE; i++) {
        if((*trajectory)[i].size() > 0) {
          if(triangulation->getVertexRank((*trajectory)[i].at(0))
             != ttk::MPIrank_) {
            globalVertexId.at(startVertexId) = -1;
          } else {
            globalVertexId.at(startVertexId) = vertIndex;
            vertIndex++;
          }
          startVertexId++;
          if((*trajectory)[i].size() > 1) {
            if(triangulation->getVertexRank((*trajectory)[i].at(0))
               != ttk::MPIrank_) {
              globalCellId.at(startCellId) = -1;
              unmatchedGhosts
                .at(neighborsToId_[triangulation->getVertexRank(
                  (*trajectory)[i].at(0))])
                .push_back(ttk::intgl::GhostElementsToSort{
                  vertIndex, (*localVertexIdentifier)[i].at(0),
                  (*seedIdentifier)[i], (*forkIdentifier)[i], -1,
                  startVertexId - 1, startCellId});
            } else {
              globalCellId.at(startCellId) = cellIndex;
              cellIndex++;
            }
            startCellId++;
            if((*trajectory)[i].size() > 2) {
              std::iota(globalCellId.begin() + startCellId,
                        globalCellId.begin() + startCellId
                          + (*trajectory)[i].size() - 2,
                        cellIndex);
              std::iota(globalVertexId.begin() + startVertexId,
                        globalVertexId.begin() + startVertexId
                          + (*trajectory)[i].size() - 2,
                        vertIndex);
              startCellId += (*trajectory)[i].size() - 2;
              startVertexId += (*trajectory)[i].size() - 2;
              vertIndex += (*trajectory)[i].size() - 2;
              cellIndex += (*trajectory)[i].size() - 2;
            }
            if(triangulation->getVertexRank((*trajectory)[i].back())
               != ttk::MPIrank_) {
              globalVertexId.at(startVertexId) = -1;
              unmatchedGhosts
                .at(neighborsToId_[triangulation->getVertexRank(
                  (*trajectory)[i].back())])
                .push_back(ttk::intgl::GhostElementsToSort{
                  globalVertexId.at(startVertexId - 1),
                  (*localVertexIdentifier)[i].at(
                    (*localVertexIdentifier)[i].size() - 2),
                  (*seedIdentifier)[i], (*forkIdentifier)[i],
                  globalCellId.at(startCellId - 1), startVertexId,
                  startCellId - 1});
            } else {
              globalVertexId.at(startVertexId) = vertIndex;
              vertIndex++;
            }
            startVertexId++;
          }
        } else {
          break;
        }
      }
      seedIdentifier++;
      forkIdentifier++;
      localVertexIdentifier++;
      trajectory++;
    }
  }
  // unmatchedGhosts contains the ghost data for each process.
  // For two processes i and j that are neighbors, their respective
  // unmatchedGhosts[k_i] and unmatchedGhosts[k_j] have the same
  // number of elements and each element represents a segment of
  // integral line. This means that if unmatchedGhosts[k_i] and
  // unmatchedGhosts[k_j] are sorted using the same comparator,
  // element l of unmatchedGhosts[k_i] and element l of unmatchedGhosts[k_j]
  // represent the same segment of integral line. Therefore, once
  // unmatchedGhosts vectors are sorted for each process, all that
  // is left to do is exchange the global identifiers of vertices
  // and edges, as the receiving process already knows to which local
  // vertices in its domain the ghosts corresponds to.
  if(ttk::MPIsize_ > 1) {
    for(int i = 0; i < neighborNumber_; i++) {
      TTK_PSORT(threadNumber_, unmatchedGhosts.at(i).begin(),
                unmatchedGhosts.at(i).end());
    }
    this->exchangeGhosts(unmatchedGhosts, globalVertexId, globalCellId);
  }

  return 0;
}

int ttkIntegralLines::exchangeGhosts(
  std::vector<std::vector<ttk::intgl::GhostElementsToSort>> &unmatchedGhosts,
  std::vector<ttk::SimplexId> &globalVertexId,
  std::vector<ttk::SimplexId> &globalCellId) {
  ttk::SimplexId id{0};
  std::vector<std::vector<ttk::SimplexId>> globalIdsToSend(neighborNumber_);
  std::vector<std::vector<ttk::SimplexId>> globalIdsToReceive(neighborNumber_);
  for(int i = 0; i < neighborNumber_; i++) {
    globalIdsToSend[i].resize(2 * unmatchedGhosts[i].size());
    globalIdsToReceive[i].resize(2 * unmatchedGhosts[i].size());
  }
  // The sending buffer is prepared
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < neighborNumber_; i++) {
    for(size_t j = 0; j < unmatchedGhosts[i].size(); j++) {
      globalIdsToSend[i][2 * j] = unmatchedGhosts[i][j].ownedGlobalId;
      globalIdsToSend[i][2 * j + 1] = unmatchedGhosts[i][j].globalEdgeId;
    }
  }
  // Each process sends and receives global identifiers of ghosts
  // to its neighbors
  for(int neighbor : neighbors_) {
    MPI_Sendrecv(globalIdsToSend[neighborsToId_[neighbor]].data(),
                 globalIdsToSend[neighborsToId_[neighbor]].size(),
                 ttk::getMPIType(id), neighbor, ttk::MPIrank_,
                 globalIdsToReceive[neighborsToId_[neighbor]].data(),
                 globalIdsToSend[neighborsToId_[neighbor]].size(),
                 ttk::getMPIType(id), neighbor, neighbor, ttk::MPIcomm_,
                 MPI_STATUS_IGNORE);
  }
  // The global identifiers of ghosts are inserted in the globalVertexId
  // and globalCellId vectors
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < neighborNumber_; i++) {
    for(size_t j = 0; j < unmatchedGhosts[i].size(); j++) {
      globalVertexId.at(unmatchedGhosts[i][j].ghostVertexLocalId)
        = globalIdsToReceive[i][2 * j];
      if(globalCellId.at(unmatchedGhosts[i][j].ghostEdgeLocalId) == -1) {
        globalCellId.at(unmatchedGhosts[i][j].ghostEdgeLocalId)
          = globalIdsToReceive[i][j * 2 + 1];
      }
    }
  }
  return 0;
}

#endif // TTK_ENABLE_MPI

int ttkIntegralLines::RequestData(vtkInformation *ttkNotUsed(request),
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {

  vtkDataSet *domain = vtkDataSet::GetData(inputVector[0], 0);
  vtkPointSet *seeds = vtkPointSet::GetData(inputVector[1], 0);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector, 0);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(domain);
  vtkDataArray *inputScalars = this->GetInputArrayToProcess(0, domain);

  vtkDataArray *inputOffsets
    = this->GetOrderArray(domain, 0, 1, ForceInputOffsetScalarField);

  const ttk::SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();
  this->setVertexNumber(numberOfPointsInDomain);
  int numberOfPointsInSeeds = seeds->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  int totalSeeds;
#else
#ifdef TTK_ENABLE_MPI
  int totalSeeds;
#endif
#endif

#ifdef TTK_ENABLE_MPI_TIME
  ttk::Timer t_mpi;
  ttk::startMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
#endif
#ifdef TTK_ENABLE_MPI
  // Necessary when using MPI
  std::vector<ttk::SimplexId> inputIdentifiers{};
  if(ttk::MPIsize_ > 1) {
    MPI_Reduce(&numberOfPointsInSeeds, &totalSeeds, 1, MPI_INTEGER, MPI_SUM, 0,
               ttk::MPIcomm_);
    int isDistributed;

    if(ttk::MPIrank_ == 0) {
      isDistributed = numberOfPointsInSeeds != totalSeeds;
    }
    MPI_Bcast(&isDistributed, 1, MPI_INTEGER, 0, ttk::MPIcomm_);
    MPI_Bcast(&totalSeeds, 1, MPI_INTEGER, 0, ttk::MPIcomm_);

    if(!isDistributed) {
      this->setGlobalElementCounter(totalSeeds);
      vtkDataArray *globalSeedsId;
      if(ttk::MPIrank_ == 0) {
        globalSeedsId = seeds->GetPointData()->GetArray("GlobalPointIds");
      } else {
        globalSeedsId = vtkDataArray::CreateDataArray(VTK_ID_TYPE);
      }

      if(ttk::MPIrank_ != 0) {
        globalSeedsId->SetNumberOfComponents(1);
        globalSeedsId->SetNumberOfTuples(totalSeeds);
      }
      ttk::LongSimplexId id = 0;
      MPI_Bcast(ttkUtils::GetPointer<ttk::LongSimplexId>(globalSeedsId),
                totalSeeds, ttk::getMPIType(id), 0, ttk::MPIcomm_);
      ttk::SimplexId localId = -1;
      for(int i = 0; i < totalSeeds; i++) {
        localId = triangulation->getVertexLocalIdIfExists(
          globalSeedsId->GetTuple1(i));
        if(localId != -1) {
          inputIdentifiers.push_back(localId);
        }
      }
      numberOfPointsInSeeds = inputIdentifiers.size();
    } else {
      std::vector<ttk::SimplexId> idSpareStorage{};
      ttk::SimplexId *inputIdentifierGlobalId;
      inputIdentifierGlobalId = this->GetIdentifierArrayPtr(
        ForceInputVertexScalarField, 2, ttk::VertexScalarFieldName, seeds,
        idSpareStorage);
      ttk::SimplexId localId = 0;
      for(int i = 0; i < numberOfPointsInSeeds; i++) {
        localId = triangulation->getVertexLocalId(inputIdentifierGlobalId[i]);
        if(localId != -1) {
          inputIdentifiers.push_back(localId);
        }
      }
      numberOfPointsInSeeds = inputIdentifiers.size();
      MPI_Allreduce(&numberOfPointsInSeeds, &totalSeeds, 1, MPI_INTEGER,
                    MPI_SUM, ttk::MPIcomm_);
      this->setGlobalElementCounter(totalSeeds);
    }
  } else {
    this->setGlobalElementCounter(numberOfPointsInSeeds);
    inputIdentifiers.resize(numberOfPointsInSeeds);
    totalSeeds = numberOfPointsInSeeds;
    std::vector<ttk::SimplexId> idSpareStorage{};
    ttk::SimplexId *inputIdentifierGlobalId;
    inputIdentifierGlobalId = this->GetIdentifierArrayPtr(
      ForceInputVertexScalarField, 2, ttk::VertexScalarFieldName, seeds,
      idSpareStorage);
    for(int i = 0; i < numberOfPointsInSeeds; i++) {
      inputIdentifiers.at(i)
        = triangulation->getVertexLocalId(inputIdentifierGlobalId[i]);
    }
  }
#else
  std::vector<ttk::SimplexId> idSpareStorage{};
  ttk::SimplexId *identifiers = this->GetIdentifierArrayPtr(
    ForceInputVertexScalarField, 2, ttk::VertexScalarFieldName, seeds,
    idSpareStorage);
  std::unordered_set<ttk::SimplexId> isSeed;
  for(ttk::SimplexId k = 0; k < numberOfPointsInSeeds; ++k) {
    isSeed.insert(identifiers[k]);
  }
  std::vector<ttk::SimplexId> inputIdentifiers(isSeed.begin(), isSeed.end());
#ifndef TTK_ENABLE_KAMIKAZE
  totalSeeds = inputIdentifiers.size();
#endif
  isSeed.clear();
#endif

  std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
    trajectories(
      threadNumber_,
      ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>());
  std::vector<ttk::ArrayLinkedList<std::vector<double>, TABULAR_SIZE>>
    distancesFromSeed(
      threadNumber_, ttk::ArrayLinkedList<std::vector<double>, TABULAR_SIZE>());
  std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
    seedIdentifiers(
      threadNumber_, ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>());
  std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
    forkIdentifiers(
      threadNumber_, ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>());
  std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
    localVertexIdentifiers(
      threadNumber_,
      ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>());

  this->setVertexNumber(numberOfPointsInDomain);
  this->setSeedNumber(numberOfPointsInSeeds);
  this->setDirection(Direction);
  this->setInputScalarField(inputScalars->GetVoidPointer(0));
  this->setInputOffsets(ttkUtils::GetPointer<ttk::SimplexId>(inputOffsets));
  this->setVertexIdentifierScalarField(&inputIdentifiers);
  this->setOutputTrajectories(&trajectories);
  this->setOutputDistancesFromSeed(&distancesFromSeed);
  this->setOutputSeedIdentifiers(&seedIdentifiers);
  this->setOutputForkIdentifiers(&forkIdentifiers);
  this->setOutputLocalVertexIdentifiers(&localVertexIdentifiers);
  this->preconditionTriangulation(triangulation);
  this->setChunkSize(
    std::max(std::max(std::min(1000, (int)numberOfPointsInSeeds),
                      (int)numberOfPointsInSeeds / (threadNumber_ * 100)),
             1));
#ifdef TTK_ENABLE_MPI
  std::vector<std::vector<std::vector<ttk::ElementToBeSent>>> toSend(
    ttk::MPIsize_);
  this->setNeighbors(triangulation->getNeighborRanks());
  if(ttk::MPIsize_ > 1) {
    toSend.resize(this->neighborNumber_);
    for(int i = 0; i < this->neighborNumber_; i++) {
      toSend[i].resize(this->threadNumber_);
      for(int j = 0; j < this->threadNumber_; j++) {
        toSend[i][j].reserve((int)numberOfPointsInSeeds * 0.005
                             / this->threadNumber_);
      }
    }
  }
  this->setToSend(&toSend);
  this->createMessageType();
#endif
#ifdef TTK_ENABLE_MPI_TIME
  double elapsedTime = ttk::endMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
  if(ttk::MPIrank_ == 0) {
    printMsg("Preparation performed using " + std::to_string(ttk::MPIsize_)
             + " MPI processes lasted :" + std::to_string(elapsedTime));
  }
#endif
#ifndef TTK_ENABLE_KAMIKAZE
  // triangulation problem
  if(!triangulation) {
    this->printErr("wrong triangulation.");
    return -1;
  }
  // field problem
  if(!inputScalars) {
    this->printErr("wrong scalar field.");
    return -1;
  }
  // field problem
  if(inputOffsets->GetDataType() != VTK_INT
     and inputOffsets->GetDataType() != VTK_ID_TYPE) {
    this->printErr("input offset field type not supported.");
    return -1;
  }

  // no points.
  if(numberOfPointsInDomain <= 0) {
    this->printErr("domain has no points.");
    return -1;
  }
  // no points.
  if(totalSeeds <= 0) {
    this->printErr("seeds have no points.");
    return -1;
  }
#endif
#ifdef TTK_ENABLE_MPI_TIME
  ttk::startMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
#endif
  int status = 0;
  ttkTemplateMacro(triangulation->getType(),
                   (status = this->execute<TTK_TT>(
                      static_cast<TTK_TT *>(triangulation->getData()))));

#ifdef TTK_ENABLE_MPI_TIME
  elapsedTime = ttk::endMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
  if(ttk::MPIrank_ == 0) {
    printMsg("Computation performed using " + std::to_string(ttk::MPIsize_)
             + " MPI processes lasted :" + std::to_string(elapsedTime));
  }
#endif
#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(status != 0) {
    this->printErr("IntegralLines.execute() error code : "
                   + std::to_string(status));
    return 0;
  }
#endif
#ifdef TTK_ENABLE_MPI
  std::vector<ttk::SimplexId> globalVertexId;
  std::vector<ttk::SimplexId> globalCellId;
  ttkTemplateMacro(triangulation->getType(),
                   (getGlobalIdentifiers<TTK_TT>(
                     globalVertexId, globalCellId, seedIdentifiers,
                     forkIdentifiers, trajectories, localVertexIdentifiers,
                     static_cast<TTK_TT *>(triangulation->getData()))));
#endif
  // make the vtk trajectories
  ttkTemplateMacro(
    triangulation->getType(),
    (getTrajectories<TTK_TT>(
      domain, static_cast<TTK_TT *>(triangulation->getData()), trajectories,
      forkIdentifiers, distancesFromSeed, seedIdentifiers,
#ifdef TTK_ENABLE_MPI
      globalVertexId, globalCellId,
#endif
      output)));

  return (int)(status == 0);
}
