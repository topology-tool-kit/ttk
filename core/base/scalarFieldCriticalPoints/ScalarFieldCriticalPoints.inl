#ifndef SCALARFIELDCRITICALPOINTS_INL
#define SCALARFIELDCRITICALPOINTS_INL

#include <ScalarFieldCriticalPoints.h>

template <class dataType>
ttk::ScalarFieldCriticalPoints<dataType>::ScalarFieldCriticalPoints() {

  dimension_ = 0;
  vertexNumber_ = 0;
  scalarValues_ = NULL;
  vertexLinkEdgeLists_ = NULL;
  criticalPoints_ = NULL;
  sosOffsets_ = NULL;
  triangulation_ = NULL;

  forceNonManifoldCheck = false;

  //   threadNumber_ = 1;
}

template <class dataType>
ttk::ScalarFieldCriticalPoints<dataType>::~ScalarFieldCriticalPoints() {
}

template <class dataType>
int ttk::ScalarFieldCriticalPoints<dataType>::execute() {

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if((!dimension_) && ((!triangulation_) || (triangulation_->isEmpty())))
    return -1;
  if((!vertexNumber_) && ((!triangulation_) || (triangulation_->isEmpty())))
    return -2;
  if(!scalarValues_)
    return -3;
  if((!vertexLinkEdgeLists_)
     && ((!triangulation_) || (triangulation_->isEmpty())))
    return -4;
  if(!criticalPoints_)
    return -5;
#endif

  if(triangulation_) {
    vertexNumber_ = triangulation_->getNumberOfVertices();
    dimension_ = triangulation_->getCellVertexNumber(0) - 1;
  }

  if(!sosOffsets_) {
    // let's use our own local copy
    sosOffsets_ = &localSosOffSets_;
  }
  if((SimplexId)sosOffsets_->size() != vertexNumber_) {
    Timer preProcess;
    sosOffsets_->resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; i++)
      (*sosOffsets_)[i] = i;

    {
      std::stringstream msg;
      msg << "[ScalarFieldCriticalPoints] Offset pre-processing done in "
          << preProcess.getElapsedTime() << " s. Go!" << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  Timer t;

  std::vector<char> vertexTypes(vertexNumber_);

  if(triangulation_) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)vertexNumber_; i++) {

      vertexTypes[i] = getCriticalType(i, triangulation_);
    }
  } else if(vertexLinkEdgeLists_) {
    // legacy implementation
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)vertexNumber_; i++) {

      vertexTypes[i] = getCriticalType(i, (*vertexLinkEdgeLists_)[i]);
    }
  }

  SimplexId minimumNumber = 0, maximumNumber = 0, saddleNumber = 0,
            oneSaddleNumber = 0, twoSaddleNumber = 0, monkeySaddleNumber = 0;

  // debug msg
  if(debugLevel_ >= Debug::infoMsg) {
    if(dimension_ == 3) {
      for(SimplexId i = 0; i < vertexNumber_; i++) {
        switch(vertexTypes[i]) {

          case static_cast<char>(CriticalType::Local_minimum):
            minimumNumber++;
            break;

          case static_cast<char>(CriticalType::Saddle1):
            oneSaddleNumber++;
            break;

          case static_cast<char>(CriticalType::Saddle2):
            twoSaddleNumber++;
            break;

          case static_cast<char>(CriticalType::Local_maximum):
            maximumNumber++;
            break;

          case static_cast<char>(CriticalType::Degenerate):
            monkeySaddleNumber++;
            break;
        }
      }
    } else if(dimension_ == 2) {
      for(SimplexId i = 0; i < vertexNumber_; i++) {
        switch(vertexTypes[i]) {

          case static_cast<char>(CriticalType::Local_minimum):
            minimumNumber++;
            break;

          case static_cast<char>(CriticalType::Saddle1):
            saddleNumber++;
            break;

          case static_cast<char>(CriticalType::Local_maximum):
            maximumNumber++;
            break;

          case static_cast<char>(CriticalType::Degenerate):
            monkeySaddleNumber++;
            break;
        }
      }
    }

    {
      std::stringstream msg;
      msg << "[ScalarFieldCriticalPoints] " << minimumNumber << " minima."
          << std::endl;
      if(dimension_ == 3) {
        msg << "[ScalarFieldCriticalPoints] " << oneSaddleNumber
            << " 1-saddle(s)." << std::endl;
        msg << "[ScalarFieldCriticalPoints] " << twoSaddleNumber
            << " 2-saddle(s)." << std::endl;
      }
      if(dimension_ == 2) {
        msg << "[ScalarFieldCriticalPoints] " << saddleNumber << " saddle(s)."
            << std::endl;
      }
      msg << "[ScalarFieldCriticalPoints] " << monkeySaddleNumber
          << " multi-saddle(s)." << std::endl;
      msg << "[ScalarFieldCriticalPoints] " << maximumNumber << " maxima."
          << std::endl;

      //       msg << "[ScalarFieldCriticalPoints] Euler characteristic
      // approximation:";
      //       if(monkeySaddleNumber){
      //         msg << " approximation";
      //       }
      //       msg << ": ";
      //       if(dimension_ == 3){
      //         msg
      //           << minimumNumber - oneSaddleNumber + twoSaddleNumber -
      // maximumNumber;
      //       }
      //       if(dimension_ == 2){
      //         msg << minimumNumber - saddleNumber + maximumNumber;
      //       }
      //       msg << std::endl;
      dMsg(std::cout, msg.str(), Debug::infoMsg);
    }
  }

  // prepare the output
  criticalPoints_->clear();
  criticalPoints_->reserve(vertexNumber_);
  for(SimplexId i = 0; i < vertexNumber_; i++) {
    if(vertexTypes[i] != static_cast<char>(CriticalType::Regular)) {
      criticalPoints_->emplace_back(i, vertexTypes[i]);
    }
  }

  {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Data-set (" << vertexNumber_
        << " vertices) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), 2);
  }

  return 0;
}

template <class dataType>
std::pair<ttk::SimplexId, ttk::SimplexId>
  ttk::ScalarFieldCriticalPoints<dataType>::getNumberOfLowerUpperComponents(
    const SimplexId vertexId, Triangulation *triangulation) const {

  SimplexId neighborNumber = triangulation->getVertexNeighborNumber(vertexId);
  std::vector<SimplexId> lowerNeighbors, upperNeighbors;

  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId = 0;
    triangulation->getVertexNeighbor(vertexId, i, neighborId);

    if(isSosLowerThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                      (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      lowerNeighbors.push_back(neighborId);
    }

    // upper link
    if(isSosHigherThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      upperNeighbors.push_back(neighborId);
    }
  }

  // shortcut, if min or max do not construct the complete star
  if(!forceNonManifoldCheck && lowerNeighbors.empty()) {
    // minimum
    return std::make_pair(0, 1);
  }

  if(!forceNonManifoldCheck && upperNeighbors.empty()) {
    // maximum
    return std::make_pair(1, 0);
  }

  // now do the actual work
  std::vector<UnionFind> lowerSeeds(lowerNeighbors.size());
  std::vector<UnionFind *> lowerList(lowerNeighbors.size());
  std::vector<UnionFind> upperSeeds(upperNeighbors.size());
  std::vector<UnionFind *> upperList(upperNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)lowerSeeds.size(); i++) {
    lowerList[i] = &(lowerSeeds[i]);
  }
  for(SimplexId i = 0; i < (SimplexId)upperSeeds.size(); i++) {
    upperList[i] = &(upperSeeds[i]);
  }

  SimplexId vertexStarSize = triangulation->getVertexStarNumber(vertexId);

  for(SimplexId i = 0; i < vertexStarSize; i++) {
    SimplexId cellId = 0;
    triangulation->getVertexStar(vertexId, i, cellId);

    SimplexId cellSize = triangulation->getCellVertexNumber(cellId);
    for(SimplexId j = 0; j < cellSize; j++) {
      SimplexId neighborId0 = -1;
      triangulation->getCellVertex(cellId, j, neighborId0);

      if(neighborId0 != vertexId) {
        // we are on the link

        bool lower0 = isSosLowerThan(
          (*sosOffsets_)[neighborId0], scalarValues_[neighborId0],
          (*sosOffsets_)[vertexId], scalarValues_[vertexId]);

        // connect it to everybody except himself and vertexId
        for(SimplexId k = j + 1; k < cellSize; k++) {

          SimplexId neighborId1 = -1;
          triangulation->getCellVertex(cellId, k, neighborId1);

          if((neighborId1 != neighborId0) && (neighborId1 != vertexId)) {

            bool lower1 = isSosLowerThan(
              (*sosOffsets_)[neighborId1], scalarValues_[neighborId1],
              (*sosOffsets_)[vertexId], scalarValues_[vertexId]);

            std::vector<SimplexId> *neighbors = &lowerNeighbors;
            std::vector<UnionFind *> *seeds = &lowerList;

            if(!lower0) {
              neighbors = &upperNeighbors;
              seeds = &upperList;
            }

            if(lower0 == lower1) {
              // connect their union-find sets!
              SimplexId lowerId0 = -1, lowerId1 = -1;
              for(SimplexId l = 0; l < (SimplexId)neighbors->size(); l++) {
                if((*neighbors)[l] == neighborId0) {
                  lowerId0 = l;
                }
                if((*neighbors)[l] == neighborId1) {
                  lowerId1 = l;
                }
              }
              if((lowerId0 != -1) && (lowerId1 != -1)) {
                (*seeds)[lowerId0]
                  = makeUnion((*seeds)[lowerId0], (*seeds)[lowerId1]);
                (*seeds)[lowerId1] = (*seeds)[lowerId0];
              }
            }
          }
        }
      }
    }
  }

  // let's remove duplicates now

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = upperList[i]->find();

  std::vector<UnionFind *>::iterator it;
  std::sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));

  std::sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId
        << ": lowerLink-#CC=" << lowerList.size()
        << " upperLink-#CC=" << upperList.size() << std::endl;

    dMsg(std::cout, msg.str(), Debug::advancedInfoMsg);
  }

  return std::make_pair(lowerList.size(), upperList.size());
}

template <class dataType>
char ttk::ScalarFieldCriticalPoints<dataType>::getCriticalType(
  const SimplexId &vertexId, Triangulation *triangulation) const {

  SimplexId downValence, upValence;
  std::tie(downValence, upValence)
    = getNumberOfLowerUpperComponents(vertexId, triangulation);

  if(downValence == 0 && upValence == 1) {
    return static_cast<char>(CriticalType::Local_minimum);
  } else if(downValence == 1 && upValence == 0) {
    return static_cast<char>(CriticalType::Local_maximum);
  } else if(downValence == 1 && upValence == 1) {
    // regular point
    return static_cast<char>(CriticalType::Regular);
  } else {
    // saddles
    if(dimension_ == 2) {
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
    } else if(dimension_ == 3) {
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

template <class dataType>
char ttk::ScalarFieldCriticalPoints<dataType>::getCriticalType(
  const SimplexId &vertexId,
  const std::vector<std::pair<SimplexId, SimplexId>> &vertexLink) const {

  std::map<SimplexId, SimplexId> global2LowerLink, global2UpperLink;
  std::map<SimplexId, SimplexId>::iterator neighborIt;

  SimplexId lowerCount = 0, upperCount = 0;

  for(SimplexId i = 0; i < (SimplexId)vertexLink.size(); i++) {

    SimplexId neighborId = vertexLink[i].first;

    // first vertex
    // lower link search
    if(isSosLowerThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                      (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      neighborIt = global2LowerLink.find(neighborId);
      if(neighborIt == global2LowerLink.end()) {
        // not in there, add it
        global2LowerLink[neighborId] = lowerCount;
        lowerCount++;
      }
    }

    // upper link
    if(isSosHigherThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      neighborIt = global2UpperLink.find(neighborId);
      if(neighborIt == global2UpperLink.end()) {
        // not in there, add it
        global2UpperLink[neighborId] = upperCount;
        upperCount++;
      }
    }

    // second vertex
    neighborId = vertexLink[i].second;

    // lower link search
    if(isSosLowerThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                      (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      neighborIt = global2LowerLink.find(neighborId);
      if(neighborIt == global2LowerLink.end()) {
        // not in there, add it
        global2LowerLink[neighborId] = lowerCount;
        lowerCount++;
      }
    }

    // upper link
    if(isSosHigherThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      neighborIt = global2UpperLink.find(neighborId);
      if(neighborIt == global2UpperLink.end()) {
        // not in there, add it
        global2UpperLink[neighborId] = upperCount;
        upperCount++;
      }
    }
  }

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId << " lower link ("
        << lowerCount << " vertices)" << std::endl;

    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId << " upper link ("
        << upperCount << " vertices)" << std::endl;

    dMsg(std::cout, msg.str(), Debug::advancedInfoMsg);
  }

  if(!lowerCount) {
    // minimum
    return static_cast<char>(CriticalType::Local_minimum);
  }
  if(!upperCount) {
    // maximum
    return static_cast<char>(CriticalType::Local_maximum);
  }

  // so far 40% of the computation, that's ok.

  // now enumerate the connected components of the lower and upper links
  // NOTE: a breadth first search might be faster than a UF
  // if so, one would need the one-skeleton data structure, not the edge list
  std::vector<UnionFind> lowerSeeds(lowerCount);
  std::vector<UnionFind> upperSeeds(upperCount);
  std::vector<UnionFind *> lowerList(lowerCount);
  std::vector<UnionFind *> upperList(upperCount);
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = &(lowerSeeds[i]);
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = &(upperSeeds[i]);

  for(SimplexId i = 0; i < (SimplexId)vertexLink.size(); i++) {

    SimplexId neighborId0 = vertexLink[i].first;
    SimplexId neighborId1 = vertexLink[i].second;

    // process the lower link
    if((isSosLowerThan((*sosOffsets_)[neighborId0], scalarValues_[neighborId0],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId]))
       && (isSosLowerThan((*sosOffsets_)[neighborId1],
                          scalarValues_[neighborId1], (*sosOffsets_)[vertexId],
                          scalarValues_[vertexId]))) {

      // both vertices are lower, let's add that edge and update the UF
      std::map<SimplexId, SimplexId>::iterator n0It
        = global2LowerLink.find(neighborId0);
      std::map<SimplexId, SimplexId>::iterator n1It
        = global2LowerLink.find(neighborId1);

      lowerList[n0It->second]
        = makeUnion(lowerList[n0It->second], lowerList[n1It->second]);
      lowerList[n1It->second] = lowerList[n0It->second];
    }

    // process the upper link
    if((isSosHigherThan((*sosOffsets_)[neighborId0], scalarValues_[neighborId0],
                        (*sosOffsets_)[vertexId], scalarValues_[vertexId]))
       && (isSosHigherThan((*sosOffsets_)[neighborId1],
                           scalarValues_[neighborId1], (*sosOffsets_)[vertexId],
                           scalarValues_[vertexId]))) {

      // both vertices are lower, let's add that edge and update the UF
      std::map<SimplexId, SimplexId>::iterator n0It
        = global2UpperLink.find(neighborId0);
      std::map<SimplexId, SimplexId>::iterator n1It
        = global2UpperLink.find(neighborId1);

      upperList[n0It->second]
        = makeUnion(upperList[n0It->second], upperList[n1It->second]);
      upperList[n1It->second] = upperList[n0It->second];
    }
  }

  // let's remove duplicates
  std::vector<UnionFind *>::iterator it;
  // update the UFs if necessary
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = upperList[i]->find();

  sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));

  sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId
        << ": lowerLink-#CC=" << lowerList.size()
        << " upperLink-#CC=" << upperList.size() << std::endl;

    dMsg(std::cout, msg.str(), Debug::advancedInfoMsg);
  }

  if((lowerList.size() == 1) && (upperList.size() == 1))
    // regular point
    return static_cast<char>(CriticalType::Regular);
  else {
    // saddles
    if(dimension_ == 2) {
      if((lowerList.size() > 2) || (upperList.size() > 2)) {
        // monkey saddle
        return static_cast<char>(CriticalType::Degenerate);
      } else {
        // regular saddle
        return static_cast<char>(CriticalType::Saddle1);
        // NOTE: you may have multi-saddles on the boundary in that
        // configuration
        // to make this computation 100% correct, one would need to disambiguate
        // boundary from interior vertices
      }
    } else if(dimension_ == 3) {
      if((lowerList.size() == 2) && (upperList.size() == 1)) {
        return static_cast<char>(CriticalType::Saddle1);
      } else if((lowerList.size() == 1) && (upperList.size() == 2)) {
        return static_cast<char>(CriticalType::Saddle2);
      } else {
        // monkey saddle
        return static_cast<char>(CriticalType::Degenerate);
        // NOTE: we may have a similar effect in 3D (TODO)
      }
    }
  }

  // -2: regular points
  return static_cast<char>(CriticalType::Regular);
}

#endif /* end of include guard: SCALARFIELDCRITICALPOINTS_INL */
