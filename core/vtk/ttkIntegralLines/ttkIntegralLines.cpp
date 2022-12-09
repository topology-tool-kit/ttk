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
  std::vector<ttk::ArrayLinkedList<std::vector<double>, TABULAR_SIZE>>
    &distancesFromSeed,
  std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
    &seedIdentifiers,
  std::vector<ttk::SimplexId> &offsets,
  std::vector<ttk::SimplexId> &edgeOffsets,
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
  vtkNew<vtkIdTypeArray> vtkEdgeIdentifiers{};
  vtkNew<vtkIntArray> vtkRankArray{};
  vtkNew<vtkIdTypeArray> vtkVertexGlobalIdArray{};

  dist->SetNumberOfComponents(1);
  dist->SetName("DistanceFromSeed");
  identifier->SetNumberOfComponents(1);
  identifier->SetName("SeedIdentifier");
  vtkVertexGlobalIdArray->SetNumberOfComponents(1);
  vtkVertexGlobalIdArray->SetName("GlobalPointIds");

#ifdef TTK_ENABLE_MPI
  vtkEdgeIdentifiers->SetNumberOfComponents(1);
  vtkEdgeIdentifiers->SetName("GlobalCellIds");
  vtkRankArray->SetNumberOfComponents(1);
  vtkRankArray->SetName("RankArray");
#else
  vtkEdgeIdentifiers->SetNumberOfComponents(1);
  vtkEdgeIdentifiers->SetName("EdgeIdentifiers");
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
#ifdef TTK_ENABLE_MPI
  const int *vertexRankArray = triangulation->getVertexRankArray();
#endif
  std::array<float, 3> p;
  std::array<vtkIdType, 2> ids;
  ttk::SimplexId globalVertexId;
  ttk::SimplexId counter = 0;
  ttk::SimplexId edgeOffset;
  for(int thread = 0; thread < threadNumber; thread++) {
    std::list<std::array<std::vector<ttk::SimplexId>, TABULAR_SIZE>>::iterator
      trajectory
      = trajectories[thread].list.begin();
    std::list<std::array<std::vector<double>, TABULAR_SIZE>>::iterator
      distanceFromSeed
      = distancesFromSeed[thread].list.begin();
    std::list<std::array<ttk::SimplexId, TABULAR_SIZE>>::iterator seedIdentifier
      = seedIdentifiers[thread].list.begin();
    while(trajectory != trajectories[thread].list.end()) {
      for(int i = 0; i < TABULAR_SIZE; i++) {
        if((*trajectory)[i].size() > 0) {
          ttk::SimplexId vertex = (*trajectory)[i].at(0);
          // init
          globalVertexId = offsets.at(counter);
          edgeOffset = edgeOffsets.at(counter);
          counter++;
          triangulation->getVertexPoint(vertex, p[0], p[1], p[2]);
          ids[0] = pts->InsertNextPoint(p.data());
          // distanceScalars
          dist->InsertNextTuple1((*distanceFromSeed)[i].at(0));
          identifier->InsertNextTuple1((*seedIdentifier)[i]);
          if(vertexRankArray[vertex] != ttk::MPIrank_) {
            globalVertexId--;
          }
          vtkVertexGlobalIdArray->InsertNextTuple1(globalVertexId);
          // inputScalars
          for(size_t k = 0; k < scalarArrays.size(); ++k) {
            inputScalars[k]->InsertNextTuple1(
              scalarArrays[k]->GetTuple1(vertex));
          }
          for(size_t j = 1; j < (*trajectory)[i].size(); ++j) {
            vtkEdgeIdentifiers->InsertNextTuple1(globalVertexId - edgeOffset);
            globalVertexId++;
            vertex = (*trajectory)[i].at(j);
            triangulation->getVertexPoint(vertex, p[0], p[1], p[2]);
            ids[1] = pts->InsertNextPoint(p.data());
            vtkVertexGlobalIdArray->InsertNextTuple1(globalVertexId);
            // distanceScalars
            dist->InsertNextTuple1((*distanceFromSeed)[i].at(j));
            identifier->InsertNextTuple1((*seedIdentifier)[i]);
            // inputScalars
            for(unsigned int k = 0; k < scalarArrays.size(); ++k)
              inputScalars[k]->InsertNextTuple1(
                scalarArrays[k]->GetTuple1(vertex));
            ug->InsertNextCell(VTK_LINE, 2, ids.data());
            // Computation of RankArray of the edge
            // WARNING: this computation is different from the
            // rest of TTK
            // If both vertices of the edge belong to the same process
            // the edge belong to that process
            // If vertices belong to different processes, the edge belongs to
            // the same process as its vertex of smallest global id.
            vtkRankArray->InsertNextTuple1(
              vertexRankArray[(*trajectory)[i].at(j - 1)]);
            // iteration
            ids[0] = ids[1];
          }
        } else {
          break;
        }
      }
      trajectory++;
      distanceFromSeed++;
      seedIdentifier++;
    }
  }

  ug->SetPoints(pts);
  ug->GetPointData()->AddArray(dist);
  ug->GetPointData()->AddArray(identifier);
#ifdef TTK_ENABLE_MPI
  ug->GetCellData()->AddArray(vtkRankArray);
#endif
  ug->GetCellData()->SetGlobalIds(vtkEdgeIdentifiers);
  for(unsigned int k = 0; k < scalarArrays.size(); ++k) {
    ug->GetPointData()->AddArray(inputScalars[k]);
  }
  ug->GetPointData()->SetGlobalIds(vtkVertexGlobalIdArray);
  output->ShallowCopy(ug);

  return 1;
}

int ttkIntegralLines::getGlobalIdentifiers(
  std::vector<ttk::SimplexId> &offsets,
  std::vector<ttk::SimplexId> &edgeOffsets,
  std::vector<ttk::ArrayLinkedList<ttk::SimplexId, TABULAR_SIZE>>
    &seedIdentifiers,
  std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
    &trajectories,
  std::vector<ttk::ArrayLinkedList<std::vector<ttk::SimplexId>, TABULAR_SIZE>>
    &localVertexIdentifiers,
  const int *vertexRankArray) {
  std::vector<IntervalToSend> intervals;
  ttk::SimplexId threadNumber = localVertexIdentifiers.size();
  ttk::SimplexId nbTrajectories = 0;
  ttk::SimplexId intervalSize = 0;
  unsigned char isFirstGhost;
  for(int thread = 0; thread < threadNumber; thread++) {
    std::list<std::array<ttk::SimplexId, TABULAR_SIZE>>::iterator seedIdentifier
      = seedIdentifiers[thread].list.begin();
    std::list<std::array<std::vector<ttk::SimplexId>, TABULAR_SIZE>>::iterator
      localVertexIdentifier
      = localVertexIdentifiers[thread].list.begin();
    std::list<std::array<std::vector<ttk::SimplexId>, TABULAR_SIZE>>::iterator
      trajectory
      = trajectories[thread].list.begin();
    while(localVertexIdentifier != localVertexIdentifiers[thread].list.end()) {
      for(int i = 0; i < TABULAR_SIZE; i++) {
        if((*localVertexIdentifier)[i].size() > 0) {
          intervalSize
            = static_cast<ttk::SimplexId>((*localVertexIdentifier)[i].size());
          if(vertexRankArray[(*trajectory)[i].at(0)] != ttk::MPIrank_) {
            intervalSize--;
            isFirstGhost = 1;
          } else {
            isFirstGhost = 0;
          }
          if((*trajectory)[i].size() > 1
             && vertexRankArray[(*trajectory)[i].back()] != ttk::MPIrank_) {
            intervalSize--;
          }
          intervals.push_back(IntervalToSend{(*localVertexIdentifier)[i].at(0),
                                             (*seedIdentifier)[i], intervalSize,
                                             isFirstGhost, ttk::MPIrank_});
          nbTrajectories++;
        } else {
          break;
        }
      }
      seedIdentifier++;
      localVertexIdentifier++;
      trajectory++;
    }
  }

  ttk::SimplexId id = 0;
  MPI_Datatype MPI_INTERVAL;
  MPI_Datatype types[]
    = {ttk::getMPIType(id), ttk::getMPIType(id), ttk::getMPIType(id),
       ttk::getMPIType(static_cast<unsigned char>(0)), MPI_INTEGER};
  int lengths[] = {1, 1, 1, 1, 1};
  const long int mpi_offsets[] = {offsetof(IntervalToSend, localVertexId),
                                  offsetof(IntervalToSend, seedIdentifier),
                                  offsetof(IntervalToSend, intervalSize),
                                  offsetof(IntervalToSend, isFirstGhost),
                                  offsetof(IntervalToSend, rankArray)};
  MPI_Type_create_struct(5, lengths, mpi_offsets, types, &MPI_INTERVAL);
  MPI_Type_commit(&MPI_INTERVAL);

  MPI_Datatype MPI_INTERVAL_OFFSET;
  MPI_Datatype typesWithOffset[]
    = {ttk::getMPIType(id), ttk::getMPIType(id), ttk::getMPIType(id),
       ttk::getMPIType(id), ttk::getMPIType(id)};
  int lengthsWithOffset[] = {1, 1, 1, 1, 1};
  const long int mpi_offsetsWithOffset[]
    = {offsetof(IntervalToSendWithOffset, localVertexId),
       offsetof(IntervalToSendWithOffset, seedIdentifier),
       offsetof(IntervalToSendWithOffset, intervalSize),
       offsetof(IntervalToSendWithOffset, offset),
       offsetof(IntervalToSendWithOffset, intervalId)};
  MPI_Type_create_struct(5, lengthsWithOffset, mpi_offsetsWithOffset,
                         typesWithOffset, &MPI_INTERVAL_OFFSET);
  MPI_Type_commit(&MPI_INTERVAL_OFFSET);

  int intervalNumber = intervals.size();
  std::vector<IntervalToSendWithOffset> intervalsWithOffset(intervalNumber);
  std::vector<int> recvcounts(ttk::MPIsize_);
  std::vector<int> displs(ttk::MPIsize_, 0);

  MPI_Gather(&intervalNumber, 1, MPI_INTEGER, recvcounts.data(), 1, MPI_INTEGER,
             0, ttk::MPIcomm_);
  for(int i = 1; i < ttk::MPIsize_; i++) {
    displs[i] = displs[i - 1] + recvcounts[i - 1];
  }
  std::vector<IntervalToSend> allIntervals(displs.back() + recvcounts.back());
  MPI_Gatherv(intervals.data(), intervalNumber, MPI_INTERVAL,
              allIntervals.data(), recvcounts.data(), displs.data(),
              MPI_INTERVAL, 0, ttk::MPIcomm_);
  ttk::SimplexId allIntervalsNumber = allIntervals.size();
  std::vector<IntervalToSendWithOffset> allIntervalsWithOffset(
    displs.back() + recvcounts.back());
  std::vector<ttk::SimplexId> currentCount(ttk::MPIsize_, 0);
  ttk::SimplexId offset = 0;
  ttk::SimplexId intervalId = 0;
  if(ttk::MPIrank_ == 0) {
    // Sort
    const auto comp = [this](const IntervalToSend i1, const IntervalToSend i2) {
      if(i1.seedIdentifier != i2.seedIdentifier) {
        return i1.seedIdentifier < i2.seedIdentifier;
      }
      if(i1.localVertexId == i2.localVertexId
         && i1.isFirstGhost != i2.isFirstGhost) {
        if((i1.isFirstGhost == 0 && i1.intervalSize == 1)
           || (i2.isFirstGhost == 0 && i2.intervalSize == 1)) {
          return (i1.isFirstGhost == 0 && i1.intervalSize == 1);
        } else {
          return i1.isFirstGhost > i2.isFirstGhost;
        }
      }
      return i1.localVertexId <= i2.localVertexId;
    };

    std::sort(allIntervals.begin(), allIntervals.end(), comp);

    allIntervalsWithOffset.at(displs[allIntervals.at(0).rankArray]
                              + currentCount[allIntervals.at(0).rankArray])
      = IntervalToSendWithOffset{allIntervals.at(0).localVertexId,
                                 allIntervals.at(0).seedIdentifier,
                                 allIntervals.at(0).intervalSize, offset, 0};
    offset += allIntervals.at(0).intervalSize;
    currentCount[allIntervals.at(0).rankArray]++;
    // Compute offset
    for(int i = 1; i < allIntervalsNumber; i++) {
      if(!(allIntervals.at(i - 1).seedIdentifier
             == allIntervals.at(i).seedIdentifier
           && allIntervals.at(i - 1).rankArray
                != allIntervals.at(i).rankArray)) {
        intervalId++;
      }
      allIntervalsWithOffset.at(displs[allIntervals.at(i).rankArray]
                                + currentCount[allIntervals.at(i).rankArray])
        = IntervalToSendWithOffset{
          allIntervals.at(i).localVertexId, allIntervals.at(i).seedIdentifier,
          allIntervals.at(i).intervalSize, offset, intervalId};
      offset += allIntervals.at(i).intervalSize;
      currentCount[allIntervals.at(i).rankArray]++;
    }
  }

  MPI_Scatterv(allIntervalsWithOffset.data(), recvcounts.data(), displs.data(),
               MPI_INTERVAL_OFFSET, intervalsWithOffset.data(), intervalNumber,
               MPI_INTERVAL_OFFSET, 0, ttk::MPIcomm_);

  // Reception of data
  offsets.reserve(nbTrajectories);
  edgeOffsets.reserve(nbTrajectories);
  std::vector<IntervalToSendWithOffset>::iterator it;
  for(int thread = 0; thread < threadNumber; thread++) {
    std::list<std::array<ttk::SimplexId, TABULAR_SIZE>>::iterator seedIdentifier
      = seedIdentifiers[thread].list.begin();
    std::list<std::array<std::vector<ttk::SimplexId>, TABULAR_SIZE>>::iterator
      localVertexIdentifier
      = localVertexIdentifiers[thread].list.begin();
    std::list<std::array<std::vector<ttk::SimplexId>, TABULAR_SIZE>>::iterator
      trajectory
      = trajectories[thread].list.begin();
    while(localVertexIdentifier != localVertexIdentifiers[thread].list.end()) {
      for(int i = 0; i < TABULAR_SIZE; i++) {
        if((*localVertexIdentifier)[i].size() > 0) {
          intervalSize
            = static_cast<ttk::SimplexId>((*localVertexIdentifier)[i].size());
          if(vertexRankArray[(*trajectory)[i].at(0)] != ttk::MPIrank_) {
            intervalSize--;
          }
          if((*trajectory)[i].size() > 1
             && vertexRankArray[(*trajectory)[i].back()] != ttk::MPIrank_) {
            intervalSize--;
          }
          it = std::find(
            intervalsWithOffset.begin(), intervalsWithOffset.end(),
            IntervalToSendWithOffset{(*localVertexIdentifier)[i].at(0),
                                     (*seedIdentifier)[i], intervalSize, 0, 0});
          if(it != intervalsWithOffset.end()) {
            offsets.push_back(
              intervalsWithOffset.at(it - intervalsWithOffset.begin()).offset);
            edgeOffsets.push_back(
              intervalsWithOffset.at(it - intervalsWithOffset.begin())
                .intervalId);
            intervalsWithOffset.erase(it);
          }
        } else {
          break;
        }
      }
      seedIdentifier++;
      localVertexIdentifier++;
      trajectory++;
    }
  }
  return 1;
}

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
#if TTK_ENABLE_MPI
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
  vertexRankArray_ = triangulation->getVertexRankArray();
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
        if(localId != -1 && vertexRankArray_[localId] == ttk::MPIrank_) {
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
        localId
          = triangulation->getVertexLocalIdIfExists(inputIdentifierGlobalId[i]);
        if(vertexRankArray_[localId] == ttk::MPIrank_) {
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
  this->setOutputLocalVertexIdentifiers(&localVertexIdentifiers);
  this->preconditionTriangulation(triangulation);
  this->setChunkSize(
    std::max(std::min(1000, (int)numberOfPointsInSeeds),
             (int)numberOfPointsInSeeds / (threadNumber_ * 100)));
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
  // field problem
  if(!inputIdentifiers.size()) {
    this->printErr("wrong identifiers.");
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
  std::vector<ttk::SimplexId> offsets;
  std::vector<ttk::SimplexId> edgeOffsets;
  this->getGlobalIdentifiers(offsets, edgeOffsets, seedIdentifiers,
                             trajectories, localVertexIdentifiers,
                             triangulation->getVertexRankArray());

  // make the vtk trajectories
  ttkTemplateMacro(
    triangulation->getType(),
    (getTrajectories<TTK_TT>(
      domain, static_cast<TTK_TT *>(triangulation->getData()), trajectories,
      distancesFromSeed, seedIdentifiers, offsets, edgeOffsets, output)));

  // // Write data to csv
  // std::ofstream myfile;
  // myfile.open("/home/eveleguillou/experiment/IntegralLines/Benchmark/"
  //             + std::to_string(ttk::MPIsize_) + "_proc_integraLines_"
  //             + std::to_string(ttk::MPIrank_) + ".csv");
  // myfile << "DistanceFromSeed,SeedIdentifier,GlobalPointIds,vtkGhostType\n";
  // vtkDataArray *ghostArray =
  // output->GetPointData()->GetArray("vtkGhostType"); vtkDataArray
  // *seedIdentifier
  //   = output->GetPointData()->GetArray("SeedIdentifier");
  // vtkDataArray *globalIdsForCsv
  //   = output->GetPointData()->GetArray("GlobalPointIds");
  // vtkDataArray *distance =
  // output->GetPointData()->GetArray("DistanceFromSeed"); for(int i = 0; i <
  // ghostArray->GetNumberOfTuples(); i++) {
  //   myfile << std::to_string(distance->GetTuple1(i)) + ","
  //               + std::to_string(seedIdentifier->GetTuple1(i)) + ","
  //               + std::to_string(globalIdsForCsv->GetTuple1(i)) + ","
  //               + std::to_string(ghostArray->GetTuple1(i)) + "\n";
  // }
  // myfile.close();

  // myfile.open("/home/eveleguillou/experiment/IntegralLines/Benchmark/"
  //             + std::to_string(ttk::MPIsize_) + "_proc_integraLinesCellData_"
  //             + std::to_string(ttk::MPIrank_) + ".csv");
  // myfile << "GlobalCellIds,RankArray\n";
  // vtkDataArray *edgeId = output->GetCellData()->GetArray("GlobalCellIds");
  // vtkDataArray *rankArray = output->GetCellData()->GetArray("RankArray");
  // for(int i = 0; i < edgeId->GetNumberOfTuples(); i++) {
  //   myfile << std::to_string(edgeId->GetTuple1(i)) + ","
  //               + std::to_string(rankArray->GetTuple1(i)) + "\n";
  // }
  // myfile.close();

  return (int)(status == 0);
}
