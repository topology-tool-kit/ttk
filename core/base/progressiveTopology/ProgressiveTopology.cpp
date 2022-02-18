#include <ProgressiveTopology.h>

int ttk::ProgressiveTopology::computeProgressivePD(
  std::vector<PersistencePair> &CTDiagram, const SimplexId *offsets) {
  int ret = -1;
  printMsg("Progressive Persistence Diagram computation");
  ret = executeCPProgressive(1, offsets);
  CTDiagram = std::move(CTDiagram_);
  CTDiagram_.clear();

  return ret;
}

int ttk::ProgressiveTopology::executeCPProgressive(
  int computePersistenceDiagram, const SimplexId *offsets) {

  printMsg(ttk::debug::Separator::L1);

  if(resumeProgressive_) {
    resumeProgressive(computePersistenceDiagram, offsets);
    return 0;
  }

  Timer timer;

  decimationLevel_ = startingDecimationLevel_;
  multiresTriangulation_.setDecimationLevel(0);
  const SimplexId vertexNumber = multiresTriangulation_.getVertexNumber();

#ifdef TTK_ENABLE_KAMIKAZE
  if(vertexNumber == 0) {
    this->printErr("No points in triangulation");
    return 1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  double tm_allocation = timer.getElapsedTime();

  // clean state (in case previous operation was a resume)
  clearResumableState();

  const auto dim = multiresTriangulation_.getDimensionality();
  const size_t maxNeigh = dim == 3 ? 14 : (dim == 2 ? 6 : 0);

  std::vector<std::vector<SimplexId>> saddleCCMin{}, saddleCCMax{};
  std::vector<std::vector<SimplexId>> vertexRepresentativesMin{},
    vertexRepresentativesMax{};
  std::vector<polarity> toPropageMin{}, toPropageMax{};
  std::vector<polarity> isUpToDateMin{}, isUpToDateMax{};
  std::vector<char> vertexTypes{};

  if(computePersistenceDiagram) {
    saddleCCMin.resize(vertexNumber);
    saddleCCMax.resize(vertexNumber);
    vertexRepresentativesMin.resize(vertexNumber);
    vertexRepresentativesMax.resize(vertexNumber);
    toPropageMin.resize(vertexNumber, 0);
    toPropageMax.resize(vertexNumber, 0);
    isUpToDateMin.resize(vertexNumber, 0);
    isUpToDateMax.resize(vertexNumber, 0);
  } else {
    vertexTypes.resize(vertexNumber, static_cast<char>(CriticalType::Regular));
  }

  // index in vertexLinkByBoundaryType
  std::vector<uint8_t> vertexLink(vertexNumber);
  VLBoundaryType vertexLinkByBoundaryType{};
  std::vector<DynamicTree> link(vertexNumber);
  std::vector<polarity> isNew(vertexNumber, 255);
  std::vector<std::vector<std::pair<polarity, polarity>>> vertexLinkPolarity(
    vertexNumber);
  std::vector<polarity> toProcess(vertexNumber, 0), toReprocess{};

  if(this->startingDecimationLevel_ > this->stoppingDecimationLevel_
     || this->isResumable_) {
    // only needed for progressive computation
    toReprocess.resize(vertexNumber, 0);
  }

  // lock vertex thread access for firstPropage
  std::vector<Lock> vertLockMin(vertexNumber), vertLockMax(vertexNumber);

  // pre-allocate memory
  if(preallocateMemory_) {
    double tm_prealloc = timer.getElapsedTime();
    printMsg("Pre-allocating data structures", 0, 0, threadNumber_,
             ttk::debug::LineMode::REPLACE);
    for(SimplexId i = 0; i < vertexNumber; ++i) {
      vertexLinkPolarity[i].reserve(maxNeigh);
      link[i].alloc(maxNeigh);
    }
    printMsg("Pre-allocating data structures", 1,
             timer.getElapsedTime() - tm_prealloc, threadNumber_);
  }

  tm_allocation = timer.getElapsedTime() - tm_allocation;
  printMsg("Total memory allocation", 1, tm_allocation, threadNumber_);

  // computation of implicit link
  std::vector<SimplexId> boundReps{};
  multiresTriangulation_.findBoundaryRepresentatives(boundReps);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < boundReps.size(); i++) {
    if(boundReps[i] != -1) {
      buildVertexLinkByBoundary(boundReps[i], vertexLinkByBoundaryType);
    }
  }

  printMsg(this->resolutionInfoString(), 0,
           timer.getElapsedTime() - tm_allocation, this->threadNumber_,
           ttk::debug::LineMode::REPLACE);
  multiresTriangulation_.setDecimationLevel(decimationLevel_);

  if(computePersistenceDiagram) {
    initSaddleSeeds(isNew, vertexLinkPolarity, toPropageMin, toPropageMax,
                    toProcess, link, vertexLink, vertexLinkByBoundaryType,
                    saddleCCMin, saddleCCMax, offsets);
    initPropagation(toPropageMin, toPropageMax, vertexRepresentativesMin,
                    vertexRepresentativesMax, saddleCCMin, saddleCCMax,
                    vertLockMin, vertLockMax, isUpToDateMin, isUpToDateMax,
                    offsets);

    // compute pairs in non-progressive mode
    computePersistencePairsFromSaddles(
      CTDiagram_, offsets, vertexRepresentativesMin, vertexRepresentativesMax,
      toPropageMin, toPropageMax);
  } else {
    initCriticalPoints(isNew, vertexLinkPolarity, toProcess, link, vertexLink,
                       vertexLinkByBoundaryType, vertexTypes, offsets);
  }

  printMsg(this->resolutionInfoString(), 1,
           timer.getElapsedTime() - tm_allocation, this->threadNumber_);

  // skip subsequent propagations if time limit is exceeded
  stopComputationIf(timer.getElapsedTime() - tm_allocation
                    > 0.9 * this->timeLimit_);

  while(decimationLevel_ > stoppingDecimationLevel_) {
    Timer tmIter{};
    decimationLevel_--;
    multiresTriangulation_.setDecimationLevel(decimationLevel_);

    printMsg(this->resolutionInfoString(), 0,
             timer.getElapsedTime() - tm_allocation, this->threadNumber_,
             ttk::debug::LineMode::REPLACE);

    if(computePersistenceDiagram) {
      updateSaddleSeeds(isNew, vertexLinkPolarity, toPropageMin, toPropageMax,
                        toProcess, toReprocess, link, vertexLink,
                        vertexLinkByBoundaryType, saddleCCMin, saddleCCMax,
                        isUpToDateMin, isUpToDateMax, offsets);
      updatePropagation(toPropageMin, toPropageMax, vertexRepresentativesMin,
                        vertexRepresentativesMax, saddleCCMin, saddleCCMax,
                        vertLockMin, vertLockMax, isUpToDateMin, isUpToDateMax,
                        offsets);
      computePersistencePairsFromSaddles(
        CTDiagram_, offsets, vertexRepresentativesMin, vertexRepresentativesMax,
        toPropageMin, toPropageMax);
    } else {
      updateCriticalPoints(isNew, vertexLinkPolarity, toProcess, toReprocess,
                           link, vertexLink, vertexLinkByBoundaryType,
                           vertexTypes, offsets);
    }

    const auto itDuration = tmIter.getElapsedTime();
    const auto nextItDuration
      = predictNextIterationDuration(itDuration, CTDiagram_.size() + 1);

    printMsg(this->resolutionInfoString(), 1,
             timer.getElapsedTime() - tm_allocation, this->threadNumber_);

    // skip subsequent propagations if time limit is exceeded
    stopComputationIf(timer.getElapsedTime() + nextItDuration - tm_allocation
                      > this->timeLimit_);

    this->printMsg("current iteration", 1.0, itDuration, 1,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
    this->printMsg("next iteration duration prediction", 1.0, nextItDuration, 1,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  // ADD GLOBAL MIN-MAX PAIR
  if(computePersistenceDiagram) {
    CTDiagram_.emplace_back(this->globalMin_, this->globalMax_, -1);
  }

  // store state for resuming computation
  if(this->isResumable_ and decimationLevel_ > 0) {
    if(computePersistenceDiagram) {
      this->vertexRepresentativesMax_ = std::move(vertexRepresentativesMax);
      this->vertexRepresentativesMin_ = std::move(vertexRepresentativesMin);
      this->saddleCCMin_ = std::move(saddleCCMin);
      this->saddleCCMax_ = std::move(saddleCCMax);
    }
    this->isNew_ = std::move(isNew);
    this->toProcess_ = std::move(toProcess);
    this->toReprocess_ = std::move(toReprocess);
    this->vertexLinkPolarity_ = std::move(vertexLinkPolarity);
    this->link_ = std::move(link);
    this->vertexLink_ = std::move(vertexLink);
    this->vertexLinkByBoundaryType_ = std::move(vertexLinkByBoundaryType);
    this->resumeProgressive_ = true;
  }

  // prepare outputs
  vertexTypes_ = std::move(vertexTypes);

  // finally sort the diagram
  sortPersistenceDiagram2(CTDiagram_, offsets);
  // sortPersistenceDiagram(CTDiagram,  offsets);
  this->printMsg(
    "Total", 1.0, timer.getElapsedTime() - tm_allocation, this->threadNumber_);
  return 0;
}

int ttk::ProgressiveTopology::resumeProgressive(int computePersistenceDiagram,
                                                const SimplexId *offsets) {

  const auto vertexNumber = multiresTriangulation_.getVertexNumber();

  this->printMsg(
    "Resuming computation from resolution level "
    + std::to_string(multiresTriangulation_.DL_to_RL(decimationLevel_))
    + " to level "
    + std::to_string(
      multiresTriangulation_.DL_to_RL(stoppingDecimationLevel_)));

  // lock vertex thread access for firstPropage
  std::vector<Lock> vertLockMin(vertexNumber), vertLockMax(vertexNumber);
  // propagation markers
  std::vector<polarity> toPropageMin{}, toPropageMax{};
  std::vector<polarity> isUpToDateMin{}, isUpToDateMax{};

  if(computePersistenceDiagram) {
    toPropageMin.resize(vertexNumber, 0);
    toPropageMax.resize(vertexNumber, 0);
    isUpToDateMin.resize(vertexNumber, 0);
    isUpToDateMax.resize(vertexNumber, 0);
  }

  Timer timer;

  while(this->decimationLevel_ > this->stoppingDecimationLevel_) {
    Timer tmIter{};
    this->decimationLevel_--;
    multiresTriangulation_.setDecimationLevel(this->decimationLevel_);
    printMsg(this->resolutionInfoString(), 0, timer.getElapsedTime(),
             this->threadNumber_, ttk::debug::LineMode::REPLACE);
    if(computePersistenceDiagram) {
      updateSaddleSeeds(isNew_, vertexLinkPolarity_, toPropageMin, toPropageMax,
                        toProcess_, toReprocess_, link_, vertexLink_,
                        vertexLinkByBoundaryType_, saddleCCMin_, saddleCCMax_,
                        isUpToDateMin, isUpToDateMax, offsets);
      updatePropagation(toPropageMin, toPropageMax, vertexRepresentativesMin_,
                        vertexRepresentativesMax_, saddleCCMin_, saddleCCMax_,
                        vertLockMin, vertLockMax, isUpToDateMin, isUpToDateMax,
                        offsets);
      computePersistencePairsFromSaddles(
        CTDiagram_, offsets, vertexRepresentativesMin_,
        vertexRepresentativesMax_, toPropageMin, toPropageMax);
    } else {
      updateCriticalPoints(isNew_, vertexLinkPolarity_, toProcess_,
                           toReprocess_, link_, vertexLink_,
                           vertexLinkByBoundaryType_, vertexTypes_, offsets);
    }

    const auto itDuration = tmIter.getElapsedTime();
    const auto nextItDuration
      = predictNextIterationDuration(itDuration, CTDiagram_.size() + 1);

    printMsg(this->resolutionInfoString(), 1, timer.getElapsedTime(),
             this->threadNumber_);

    // skip subsequent propagations if time limit is exceeded
    stopComputationIf(timer.getElapsedTime() + nextItDuration
                      > this->timeLimit_);
  }

  // ADD GLOBAL MIN-MAX PAIR
  if(computePersistenceDiagram) {
    CTDiagram_.emplace_back(this->globalMin_, this->globalMax_, -1);
  }
  // finally sort the diagram
  sortPersistenceDiagram2(CTDiagram_, offsets);
  // sortPersistenceDiagram(CTDiagram,  offsets);
  this->printMsg("Total", 1.0, timer.getElapsedTime(), this->threadNumber_);

  // clean state (we don't need it anymore)
  if(this->decimationLevel_ == 0) {
    clearResumableState();
  }

  return 0;
}

void ttk::ProgressiveTopology::computePersistencePairsFromSaddles(
  std::vector<PersistencePair> &CTDiagram,
  const SimplexId *const offsets,
  std::vector<std::vector<SimplexId>> &vertexRepresentativesMin,
  std::vector<std::vector<SimplexId>> &vertexRepresentativesMax,
  const std::vector<polarity> &toPropageMin,
  const std::vector<polarity> &toPropageMax) const {

  Timer timer{};
  std::vector<triplet> tripletsMax{}, tripletsMin{};
  const SimplexId nbDecVert = multiresTriangulation_.getDecimatedVertexNumber();

  for(SimplexId localId = 0; localId < nbDecVert; localId++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(localId);
    if(toPropageMin[globalId]) {
      getTripletsFromSaddles(globalId, tripletsMin, vertexRepresentativesMin);
    }
    if(toPropageMax[globalId]) {
      getTripletsFromSaddles(globalId, tripletsMax, vertexRepresentativesMax);
    }
  }

  this->printMsg("TRIPLETS", 1.0, timer.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  double tm_pairs = timer.getElapsedTime();

  sortTriplets(tripletsMax, offsets, true);
  sortTriplets(tripletsMin, offsets, false);

  const auto tm_sort = timer.getElapsedTime();
  this->printMsg("TRIPLETS SORT", 1.0, tm_sort - tm_pairs, this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  typename std::remove_reference<decltype(CTDiagram)>::type CTDiagramMin{},
    CTDiagramMax{};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif // TTK_ENABLE_OPENMP
    tripletsToPersistencePairs(
      CTDiagramMin, vertexRepresentativesMax, tripletsMax, offsets, true);
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif // TTK_ENABLE_OPENMP
    tripletsToPersistencePairs(
      CTDiagramMax, vertexRepresentativesMin, tripletsMin, offsets, false);
  }
  CTDiagram = std::move(CTDiagramMin);
  CTDiagram.insert(CTDiagram.end(), CTDiagramMax.begin(), CTDiagramMax.end());

  this->printMsg("PAIRS", 1.0, timer.getElapsedTime() - tm_sort,
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
}

void ttk::ProgressiveTopology::sortTriplets(std::vector<triplet> &triplets,
                                            const SimplexId *const offsets,
                                            const bool splitTree) const {
  if(triplets.empty())
    return;

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return offsets[a] < offsets[b];
  };

  // Sorting step
  const auto cmp = [=](const triplet &t1, const triplet &t2) {
    const SimplexId s1 = std::get<0>(t1);
    const SimplexId s2 = std::get<0>(t2);
    const SimplexId m1 = std::get<2>(t1);
    const SimplexId m2 = std::get<2>(t2);
    if(s1 != s2)
      return lt(s1, s2) != splitTree;
    else // s1 == s2
      return lt(m1, m2) == splitTree;
  };

  TTK_PSORT(this->threadNumber_, triplets.begin(), triplets.end(), cmp);
}

void ttk::ProgressiveTopology::tripletsToPersistencePairs(
  std::vector<PersistencePair> &pairs,
  std::vector<std::vector<SimplexId>> &vertexRepresentatives,
  std::vector<triplet> &triplets,
  const SimplexId *const offsets,
  const bool splitTree) const {

  Timer tm;
  if(triplets.empty())
    return;
  size_t numberOfPairs = 0;

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return offsets[a] < offsets[b];
  };

  const auto getRep = [&](SimplexId v) -> SimplexId {
    auto r = vertexRepresentatives[v][0];
    while(r != v) {
      v = r;
      r = vertexRepresentatives[v][0];
    }
    return r;
  };

  for(const auto &t : triplets) {
    SimplexId r1 = getRep(std::get<1>(t));
    SimplexId r2 = getRep(std::get<2>(t));
    if(r1 != r2) {
      SimplexId s = std::get<0>(t);
      numberOfPairs++;

      // Add pair
      if(splitTree) {
        // r1 = min(r1, r2), r2 = max(r1, r2)
        if(lt(r2, r1)) {
          std::swap(r1, r2);
        }
        // pair saddle-max: s -> min(r1, r2);
        pairs.emplace_back(s, r1, 2);

      } else {
        // r1 = max(r1, r2), r2 = min(r1, r2)
        if(lt(r1, r2)) {
          std::swap(r1, r2);
        }
        // pair min-saddle: max(r1, r2) -> s;
        pairs.emplace_back(r1, s, 0);
      }

      vertexRepresentatives[std::get<1>(t)][0] = r2;
      vertexRepresentatives[r1][0] = r2;
    }
  }

  const std::string ptype = splitTree ? "sad-max" : "min-sad";
  this->printMsg("found all " + ptype + " pairs", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
}

void ttk::ProgressiveTopology::buildVertexLinkPolarity(
  const SimplexId vertexId,
  std::vector<std::pair<polarity, polarity>> &vlp,
  const SimplexId *const offsets) const {

  const SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  vlp.resize(neighborNumber);
  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId0 = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId0);
    const bool lower0 = offsets[neighborId0] < offsets[vertexId];
    const polarity isUpper0 = static_cast<polarity>(!lower0) * 255;
    vlp[i] = std::make_pair(isUpper0, 0);
  }
}

void ttk::ProgressiveTopology::initDynamicLink(
  const SimplexId &vertexId,
  std::vector<std::pair<polarity, polarity>> &vlp,
  uint8_t &vertexLink,
  DynamicTree &link,
  VLBoundaryType &vlbt,
  const SimplexId *const offsets) const {

  if(vlp.empty()) {
    buildVertexLinkPolarity(vertexId, vlp, offsets);
  }

  const SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  link.alloc(neighborNumber);

  // associate vertex link boundary
  vertexLink = multiresTriangulation_.getVertexBoundaryIndex(vertexId);

  // update the link polarity for old points that are processed for
  // the first time
  const auto &vl = vlbt[vertexLink];
  for(size_t edgeId = 0; edgeId < vl.size(); edgeId++) {
    const SimplexId n0 = vl[edgeId].first;
    const SimplexId n1 = vl[edgeId].second;
    if(vlp[n0].first == vlp[n1].first) {
      // the smallest id (n0) becomes the parent of n1
      link.insertEdge(n1, n0);
    }
  }
}

void ttk::ProgressiveTopology::updateCriticalPoints(
  std::vector<polarity> &isNew,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  std::vector<polarity> &toProcess,
  std::vector<polarity> &toReprocess,
  std::vector<DynamicTree> &link,
  std::vector<uint8_t> &vertexLink,
  VLBoundaryType &vertexLinkByBoundaryType,
  std::vector<char> &vertexTypes,
  const SimplexId *const offsets) const {

  Timer tm;
  const auto nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  // find breaking edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId localId = 0; localId < nDecVerts; localId++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(localId);
    if(isNew[globalId]) {
      if(decimationLevel_ > stoppingDecimationLevel_ || isResumable_) {
        buildVertexLinkPolarity(
          globalId, vertexLinkPolarity[globalId], offsets);
      }
    } else {
      getMonotonyChangeByOldPointCP(globalId, isNew, toProcess, toReprocess,
                                    vertexLinkPolarity[globalId], offsets);
    }
  }

  this->printMsg("MONOTONY", 1.0, tm.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  double t_critical = tm.getElapsedTime();
  // second Loop  process or reprocess
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < nDecVerts; i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);

    if(isNew[globalId]) { // new point
      if(toProcess[globalId]) {
        initDynamicLink(globalId, vertexLinkPolarity[globalId],
                        vertexLink[globalId], link[globalId],
                        vertexLinkByBoundaryType, offsets);
        vertexTypes[globalId] = getCriticalTypeFromLink(
          vertexLinkPolarity[globalId], link[globalId]);
      }
      isNew[globalId] = false;

    } else { // old point
      if(toReprocess[globalId]) {
        if(toProcess[globalId]) { // was already processed : need to reprocess
          updateDynamicLink(link[globalId], vertexLinkPolarity[globalId],
                            vertexLinkByBoundaryType[vertexLink[globalId]]);
        } else { // first processing
          updateLinkPolarity(globalId, vertexLinkPolarity[globalId], offsets);
          initDynamicLink(globalId, vertexLinkPolarity[globalId],
                          vertexLink[globalId], link[globalId],
                          vertexLinkByBoundaryType, offsets);
          toProcess[globalId] = 255; // mark as processed
        }
        vertexTypes[globalId] = getCriticalTypeFromLink(
          vertexLinkPolarity[globalId], link[globalId]);
        toReprocess[globalId] = 0;
      }
    }

  } // end for openmp

  this->printMsg("CRITICAL POINTS UPDATE", 1.0,
                 tm.getElapsedTime() - t_critical, this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);
}

void ttk::ProgressiveTopology::updateSaddleSeeds(
  std::vector<polarity> &isNew,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  std::vector<polarity> &toPropageMin,
  std::vector<polarity> &toPropageMax,
  std::vector<polarity> &toProcess,
  std::vector<polarity> &toReprocess,
  std::vector<DynamicTree> &link,
  std::vector<uint8_t> &vertexLink,
  VLBoundaryType &vertexLinkByBoundaryType,
  std::vector<std::vector<SimplexId>> &saddleCCMin,
  std::vector<std::vector<SimplexId>> &saddleCCMax,
  std::vector<polarity> &isUpdatedMin,
  std::vector<polarity> &isUpdatedMax,
  const SimplexId *const offsets) const {

  Timer tm;
  const auto nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  // find breaking edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId localId = 0; localId < nDecVerts; localId++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(localId);
    if(isNew[globalId]) {
      if(decimationLevel_ > stoppingDecimationLevel_ || isResumable_) {
        buildVertexLinkPolarity(
          globalId, vertexLinkPolarity[globalId], offsets);
      }
    } else {
      getMonotonyChangeByOldPointCP(globalId, isNew, toProcess, toReprocess,
                                    vertexLinkPolarity[globalId], offsets);
    }
  }

  this->printMsg("MONOTONY", 1.0, tm.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  double t_critical = tm.getElapsedTime();
  // second Loop  process or reprocess
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < nDecVerts; i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);

    if(isNew[globalId]) { // new point
      if(toProcess[globalId]) {
        initDynamicLink(globalId, vertexLinkPolarity[globalId],
                        vertexLink[globalId], link[globalId],
                        vertexLinkByBoundaryType, offsets);
        getValencesFromLink(globalId, vertexLinkPolarity[globalId],
                            link[globalId], toPropageMin, toPropageMax,
                            saddleCCMin, saddleCCMax);
      }
      isNew[globalId] = false;

    } else { // old point
      if(toReprocess[globalId]) {
        if(toProcess[globalId]) { // was already processed : need to reprocess
          updateDynamicLink(link[globalId], vertexLinkPolarity[globalId],
                            vertexLinkByBoundaryType[vertexLink[globalId]]);
        } else { // first processing
          updateLinkPolarity(globalId, vertexLinkPolarity[globalId], offsets);
          initDynamicLink(globalId, vertexLinkPolarity[globalId],
                          vertexLink[globalId], link[globalId],
                          vertexLinkByBoundaryType, offsets);
          toProcess[globalId] = 255; // mark as processed
        }
        getValencesFromLink(globalId, vertexLinkPolarity[globalId],
                            link[globalId], toPropageMin, toPropageMax,
                            saddleCCMin, saddleCCMax);
        toReprocess[globalId] = 0;
      }
    }

    // reset updated flag
    isUpdatedMin[globalId] = 0;
    isUpdatedMax[globalId] = 0;
  } // end for openmp

  this->printMsg("CRITICAL POINTS UPDATE", 1.0,
                 tm.getElapsedTime() - t_critical, this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);
}

bool ttk::ProgressiveTopology::getMonotonyChangeByOldPointCP(
  const SimplexId vertexId,
  const std::vector<polarity> &isNew,
  std::vector<polarity> &toProcess,
  std::vector<polarity> &toReprocess,
  std::vector<std::pair<polarity, polarity>> &vlp,
  const SimplexId *const offsets) const {

  bool hasMonotonyChanged = false;
  const SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);

    // check for monotony changes
    const bool lower = offsets[neighborId] < offsets[vertexId];
    const polarity isUpper = lower ? 0 : 255;
    const polarity isUpperOld = vlp[i].first;

    if(isUpper != isUpperOld) { // change of monotony
      hasMonotonyChanged = true;

      toReprocess[vertexId] = 255;
      toProcess[neighborId] = 255;
      const SimplexId neighborNumberNew
        = multiresTriangulation_.getVertexNeighborNumber(neighborId);
      for(SimplexId j = 0; j < neighborNumberNew; j++) {
        SimplexId neighborIdNew = -1;
        multiresTriangulation_.getVertexNeighbor(neighborId, j, neighborIdNew);
        if(isNew[neighborIdNew])
          toProcess[neighborIdNew] = 255;
      }

      vlp[i].second = 255;
    }
  }
  return hasMonotonyChanged;
}

ttk::SimplexId ttk::ProgressiveTopology::propageFromSaddles(
  const SimplexId vertexId,
  std::vector<Lock> &vertLock,
  std::vector<polarity> &toPropage,
  std::vector<std::vector<SimplexId>> &vertexRepresentatives,
  std::vector<std::vector<SimplexId>> &saddleCC,
  std::vector<polarity> &isUpdated,
  std::vector<SimplexId> &globalExtremum,
  const SimplexId *const offsets,
  const bool splitTree) const {

  auto &toProp = toPropage[vertexId];
  auto &reps = vertexRepresentatives[vertexId];
  auto &updated = isUpdated[vertexId];

  const auto gt = [=](const SimplexId a, const SimplexId b) {
    return (offsets[a] > offsets[b]) == splitTree;
  };

  if(updated) {
    return reps[0];
  }

  if(this->threadNumber_ > 1) {
    vertLock[vertexId].lock();
  }

  if(toProp) { // SADDLE POINT
    const auto &CC = saddleCC[vertexId];
    reps.clear();
    reps.reserve(CC.size());
    for(size_t r = 0; r < CC.size(); r++) {
      SimplexId neighborId = -1;
      SimplexId localId = CC[r];
      multiresTriangulation_.getVertexNeighbor(vertexId, localId, neighborId);
      SimplexId ret = propageFromSaddles(
        neighborId, vertLock, toPropage, vertexRepresentatives, saddleCC,
        isUpdated, globalExtremum, offsets, splitTree);
      reps.emplace_back(ret);
    }

    if(reps.size() > 1) {
      // sort & remove duplicate elements
      std::sort(reps.begin(), reps.end(), gt);
      const auto last = std::unique(reps.begin(), reps.end());
      reps.erase(last, reps.end());
    }

    updated = 255;
    if(this->threadNumber_ > 1) {
      vertLock[vertexId].unlock();
    }

    return reps[0];

  } else {

    SimplexId ret = vertexId;
    SimplexId neighborNumber
      = multiresTriangulation_.getVertexNeighborNumber(vertexId);
    SimplexId maxNeighbor = vertexId;
    for(SimplexId i = 0; i < neighborNumber; i++) {
      SimplexId neighborId = -1;
      multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);
      if(gt(neighborId, maxNeighbor)) {
        maxNeighbor = neighborId;
      }
    }
    if(maxNeighbor != vertexId) { // not an extremum
      ret = propageFromSaddles(maxNeighbor, vertLock, toPropage,
                               vertexRepresentatives, saddleCC, isUpdated,
                               globalExtremum, offsets, splitTree);
    } else {
#ifdef TTK_ENABLE_OPENMP
      const auto tid = omp_get_thread_num();
#else
      const auto tid = 0;
#endif // TTK_ENABLE_OPENMP
      if(gt(vertexId, globalExtremum[tid])) {
        globalExtremum[tid] = vertexId;
      }
    }
    reps.resize(1);
    reps[0] = ret;
    updated = 255;
    if(this->threadNumber_ > 1) {
      vertLock[vertexId].unlock();
    }
    return ret;
  }
}

void ttk::ProgressiveTopology::initCriticalPoints(
  std::vector<polarity> &isNew,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  std::vector<polarity> &toProcess,
  std::vector<DynamicTree> &link,
  std::vector<uint8_t> &vertexLink,
  VLBoundaryType &vertexLinkByBoundaryType,
  std::vector<char> &vertexTypes,
  const SimplexId *const offsets) const {

  Timer timer{};
  const size_t nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  // computes the critical types of all points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
    buildVertexLinkPolarity(globalId, vertexLinkPolarity[globalId], offsets);
    initDynamicLink(globalId, vertexLinkPolarity[globalId],
                    vertexLink[globalId], link[globalId],
                    vertexLinkByBoundaryType, offsets);
    vertexTypes[globalId]
      = getCriticalTypeFromLink(vertexLinkPolarity[globalId], link[globalId]);
    toProcess[globalId] = 255;
    isNew[globalId] = 0;
  }

  this->printMsg("initial critical types", 1.0, timer.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
}

void ttk::ProgressiveTopology::initSaddleSeeds(
  std::vector<polarity> &isNew,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  std::vector<polarity> &toPropageMin,
  std::vector<polarity> &toPropageMax,
  std::vector<polarity> &toProcess,
  std::vector<DynamicTree> &link,
  std::vector<uint8_t> &vertexLink,
  VLBoundaryType &vertexLinkByBoundaryType,
  std::vector<std::vector<SimplexId>> &saddleCCMin,
  std::vector<std::vector<SimplexId>> &saddleCCMax,
  const SimplexId *const offsets) const {

  Timer timer{};
  const size_t nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  // computes the critical types of all points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
    buildVertexLinkPolarity(globalId, vertexLinkPolarity[globalId], offsets);
    initDynamicLink(globalId, vertexLinkPolarity[globalId],
                    vertexLink[globalId], link[globalId],
                    vertexLinkByBoundaryType, offsets);
    getValencesFromLink(globalId, vertexLinkPolarity[globalId], link[globalId],
                        toPropageMin, toPropageMax, saddleCCMin, saddleCCMax);
    toProcess[globalId] = 255;
    isNew[globalId] = 0;
  }

  this->printMsg("initial critical types", 1.0, timer.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
}

void ttk::ProgressiveTopology::initPropagation(
  std::vector<polarity> &toPropageMin,
  std::vector<polarity> &toPropageMax,
  std::vector<std::vector<SimplexId>> &vertexRepresentativesMin,
  std::vector<std::vector<SimplexId>> &vertexRepresentativesMax,
  std::vector<std::vector<SimplexId>> &saddleCCMin,
  std::vector<std::vector<SimplexId>> &saddleCCMax,
  std::vector<Lock> &vertLockMin,
  std::vector<Lock> &vertLockMax,
  std::vector<polarity> &isUpdatedMin,
  std::vector<polarity> &isUpdatedMax,
  const SimplexId *const offsets) const {

  Timer timer{};
  const size_t nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  std::vector<SimplexId> globalMaxThr(threadNumber_, 0);
  std::vector<SimplexId> globalMinThr(threadNumber_, 0);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId v = multiresTriangulation_.localToGlobalVertexId(i);
    if(toPropageMin[v]) {
      propageFromSaddles(v, vertLockMin, toPropageMin, vertexRepresentativesMin,
                         saddleCCMin, isUpdatedMin, globalMinThr, offsets,
                         false);
    }
    if(toPropageMax[v]) {
      propageFromSaddles(v, vertLockMax, toPropageMax, vertexRepresentativesMax,
                         saddleCCMax, isUpdatedMax, globalMaxThr, offsets,
                         true);
    }
  }

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return offsets[a] < offsets[b];
  };

  globalMin_ = *std::min_element(globalMinThr.begin(), globalMinThr.end(), lt);
  globalMax_ = *std::max_element(globalMaxThr.begin(), globalMaxThr.end(), lt);

  this->printMsg("FIRSTPROPAGATION", 1.0, timer.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
}

void ttk::ProgressiveTopology::updatePropagation(
  std::vector<polarity> &toPropageMin,
  std::vector<polarity> &toPropageMax,
  std::vector<std::vector<SimplexId>> &vertexRepresentativesMin,
  std::vector<std::vector<SimplexId>> &vertexRepresentativesMax,
  std::vector<std::vector<SimplexId>> &saddleCCMin,
  std::vector<std::vector<SimplexId>> &saddleCCMax,
  std::vector<Lock> &vertLockMin,
  std::vector<Lock> &vertLockMax,
  std::vector<polarity> &isUpdatedMin,
  std::vector<polarity> &isUpdatedMax,
  const SimplexId *const offsets) const {

  Timer tm{};
  const size_t nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  std::vector<SimplexId> globalMaxThr(threadNumber_, 0);
  std::vector<SimplexId> globalMinThr(threadNumber_, 0);

  // propage along split tree
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId v = multiresTriangulation_.localToGlobalVertexId(i);
    if(toPropageMin[v]) {
      propageFromSaddles(v, vertLockMin, toPropageMin, vertexRepresentativesMin,
                         saddleCCMin, isUpdatedMin, globalMinThr, offsets,
                         false);
    }
    if(toPropageMax[v]) {
      propageFromSaddles(v, vertLockMax, toPropageMax, vertexRepresentativesMax,
                         saddleCCMax, isUpdatedMax, globalMaxThr, offsets,
                         true);
    }
  }

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return offsets[a] < offsets[b];
  };

  globalMin_ = *std::min_element(globalMinThr.begin(), globalMinThr.end(), lt);
  globalMax_ = *std::max_element(globalMaxThr.begin(), globalMaxThr.end(), lt);

  this->printMsg("PROPAGATION UPDATE", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
}

void ttk::ProgressiveTopology::updateLinkPolarity(
  const SimplexId vertexId,
  std::vector<std::pair<polarity, polarity>> &vlp,
  const SimplexId *const offsets) const {

  for(size_t i = 0; i < vlp.size(); i++) {
    SimplexId neighborId = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);
    const bool lower = offsets[neighborId] < offsets[vertexId];
    const polarity isUpper = lower ? 0 : 255;
    vlp[i] = std::make_pair(isUpper, 0);
  }
}

void ttk::ProgressiveTopology::updateDynamicLink(
  DynamicTree &link,
  std::vector<std::pair<polarity, polarity>> &vlp,
  std::vector<std::pair<SimplexId, SimplexId>> &vl) const {

  std::vector<SimplexId> monotony_changes_list{};

  for(size_t neighborId = 0; neighborId < vlp.size(); neighborId++) {
    const polarity isBroken = vlp[neighborId].second;
    if(isBroken != 0) {
      monotony_changes_list.emplace_back(neighborId);
    }
  }

  // loop on the link
  //   for each edge that shares n0
  //      if only one break and different polarity : remove
  //      else if only one break and same polarity : insert
  //      else : do nothing
  std::vector<std::vector<std::pair<SimplexId, SimplexId>>> edgesToInsertLater(
    vlp.size());
  std::vector<std::vector<std::pair<SimplexId, SimplexId>>> edgesToRemoveLater(
    vlp.size());

  for(size_t e = 0; e < vl.size(); e++) {
    const SimplexId n0 = vl[e].first;
    const SimplexId n1 = vl[e].second;
    const polarity isBroken0 = vlp[n0].second;
    const polarity isBroken1 = vlp[n1].second;

    if(isBroken0 != 0 and isBroken1 == 0) {
      if(vlp[n0].first != vlp[n1].first) {
        edgesToInsertLater[n0].emplace_back(n0, n1);
      } else {
        edgesToRemoveLater[n0].emplace_back(n0, n1);
      }
    } else if(isBroken0 == 0 and isBroken1 != 0) {
      if(vlp[n0].first != vlp[n1].first) {
        edgesToInsertLater[n1].emplace_back(n1, n0);
      } else {
        edgesToRemoveLater[n1].emplace_back(n1, n0);
      }
    }
  }

  // REMOVE EDGES:
  for(const auto brokenNode : monotony_changes_list) {
    vlp[brokenNode].first = ~vlp[brokenNode].first;
    vlp[brokenNode].second = 0;

    link.getNode(brokenNode)->evert();
    for(const auto &edge : edgesToRemoveLater[brokenNode]) {
      link.removeEdge(edge.first, edge.second);
    }
  }
  if(!edgesToRemoveLater.empty()) {
    // reconnect link
    for(const auto &edge : vl) {
      if(vlp[edge.first].first == vlp[edge.second].first) {
        link.insertEdge(edge.first, edge.second);
      }
    }
  }

  // INSERT EDGES
  for(const auto brokenNode : monotony_changes_list) {
    for(const auto &edge : edgesToInsertLater[brokenNode]) {
      link.insertEdge(edge.first, edge.second);
    }
  }
}

double ttk::ProgressiveTopology::predictNextIterationDuration(
  const double currItDuration, const size_t nCurrPairs) const {

  // number of vertices at current iteration
  const double nCurrVerts = multiresTriangulation_.getDecimatedVertexNumber();
  // prediction of duration at iteration n + 1 from iteration n
  // (linear regression, R^2 = 0.994)
  return -0.21 + 0.77 / (decimationLevel_ + 1) - 4.0 * nCurrPairs / nCurrVerts
         + currItDuration
             * (3.3 - 2.32 / nCurrPairs + 32.3 * nCurrPairs / nCurrVerts);
}

void ttk::ProgressiveTopology::stopComputationIf(const bool b) {
  if(b) {
    if(this->decimationLevel_ > this->stoppingDecimationLevel_) {
      this->printMsg("Computation stopped at resolution level "
                     + std::to_string(multiresTriangulation_.DL_to_RL(
                       this->decimationLevel_)));
    }
    this->stoppingDecimationLevel_ = this->decimationLevel_;
  }
}

void ttk::ProgressiveTopology::clearResumableState() {
  // force de-allocation
  vertexRepresentativesMin_ = {};
  vertexRepresentativesMax_ = {};
  vertexLinkPolarity_ = {};
  isNew_ = {};
  vertexLink_ = {};
  link_ = {};
  toProcess_ = {};
  toReprocess_ = {};
  saddleCCMin_ = {};
  saddleCCMax_ = {};
  vertexTypes_ = {};
}

char ttk::ProgressiveTopology::getCriticalTypeFromLink(
  const std::vector<std::pair<polarity, polarity>> &vlp,
  DynamicTree &link) const {

  const auto nbCC = link.getNbCC();

  int dimensionality = multiresTriangulation_.getDimensionality();
  SimplexId downValence = 0, upValence = 0;

  std::vector<size_t> CCIds;
  CCIds.reserve(nbCC);
  link.retrieveNbCC(CCIds);
  for(size_t i = 0; i < CCIds.size(); i++) {
    const SimplexId neighbor = CCIds[i];
    const polarity isUpper = vlp[neighbor].first;
    if(isUpper) {
      upValence++;
    } else {
      downValence++;
    }
  }

  if(downValence == -1 && upValence == -1) {
    return -1;
  } else if(downValence == 0 && upValence == 1) {
    return static_cast<char>(CriticalType::Local_minimum);
  } else if(downValence == 1 && upValence == 0) {
    return static_cast<char>(CriticalType::Local_maximum);
  } else if(downValence == 1 && upValence == 1) {
    // regular point
    return static_cast<char>(CriticalType::Regular);
  } else {
    // saddles
    if(dimensionality == 2) {
      if((downValence == 2 && upValence == 1)
         || (downValence == 1 && upValence == 2)
         || (downValence == 2 && upValence == 2)) {
        // regular saddle
        return static_cast<char>(CriticalType::Saddle1);
      } else {
        // monkey saddle, saddle + extremum
        return static_cast<char>(CriticalType::Degenerate);
        // NOTE: you may have multi-saddles on the boundary in that
        // configuration
        // to make this computation 100% correct, one would need to
        // disambiguate boundary from interior vertices
      }
    } else if(dimensionality == 3) {
      if(downValence == 2 && upValence == 1) {
        return static_cast<char>(CriticalType::Saddle1);
      } else if(downValence == 1 && upValence == 2) {
        return static_cast<char>(CriticalType::Saddle2);
      } else {
        // monkey saddle, saddle + extremum
        return static_cast<char>(CriticalType::Degenerate);
        // NOTE: we may have a similar effect in 3D (TODO)
      }
    }
  }

  // -2: regular points
  return static_cast<char>(CriticalType::Regular);
}

void ttk::ProgressiveTopology::sortPersistenceDiagram2(
  std::vector<PersistencePair> &diagram, const SimplexId *const offsets) const {

  auto cmp = [offsets](const PersistencePair &a, const PersistencePair &b) {
    return offsets[a.birth] < offsets[b.birth];
  };

  std::sort(diagram.begin(), diagram.end(), cmp);
}

int ttk::ProgressiveTopology::computeProgressiveCP(
  std::vector<std::pair<SimplexId, char>> *criticalPoints,
  const SimplexId *offsets) {

  printMsg("Progressive Critical Points computation");
  int ret = -1;
  ret = executeCPProgressive(0, offsets);

  SimplexId vertexNumber = multiresTriangulation_.getVertexNumber();
  criticalPoints->clear();
  criticalPoints->reserve(vertexNumber);
  for(SimplexId i = 0; i < vertexNumber; i++) {
    if(vertexTypes_[i] != (char)(CriticalType::Regular)) {
      criticalPoints->emplace_back(i, vertexTypes_[i]);
    }
  }
  return ret;
}
