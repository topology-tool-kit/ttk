#include <ScalarFieldCriticalPoints.h>

ttk::ScalarFieldCriticalPoints::ScalarFieldCriticalPoints() {
  this->setDebugMsgPrefix("ScalarFieldCriticalPoints");
}

char ttk::ScalarFieldCriticalPoints::getCriticalType(
  const SimplexId &vertexId,
  const SimplexId *const offsets,
  const std::vector<std::pair<SimplexId, SimplexId>> &vertexLink) const {

  std::map<SimplexId, SimplexId> global2LowerLink, global2UpperLink;
  std::map<SimplexId, SimplexId>::iterator neighborIt;

  SimplexId lowerCount = 0, upperCount = 0;

  for(SimplexId i = 0; i < (SimplexId)vertexLink.size(); i++) {

    SimplexId neighborId = vertexLink[i].first;

    // first vertex
    // lower link search
    if(offsets[neighborId] < offsets[vertexId]) {

      neighborIt = global2LowerLink.find(neighborId);
      if(neighborIt == global2LowerLink.end()) {
        // not in there, add it
        global2LowerLink[neighborId] = lowerCount;
        lowerCount++;
      }
    }

    // upper link
    if(offsets[neighborId] > offsets[vertexId]) {

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
    if(offsets[neighborId] < offsets[vertexId]) {

      neighborIt = global2LowerLink.find(neighborId);
      if(neighborIt == global2LowerLink.end()) {
        // not in there, add it
        global2LowerLink[neighborId] = lowerCount;
        lowerCount++;
      }
    }

    // upper link
    if(offsets[neighborId] > offsets[vertexId]) {

      neighborIt = global2UpperLink.find(neighborId);
      if(neighborIt == global2UpperLink.end()) {
        // not in there, add it
        global2UpperLink[neighborId] = upperCount;
        upperCount++;
      }
    }
  }

  if(debugLevel_ >= (int)(debug::Priority::VERBOSE)) {
    printMsg("Vertex #" + std::to_string(vertexId) + " lower link ("
               + std::to_string(lowerCount) + " vertices)",
             debug::Priority::VERBOSE);
    printMsg("Vertex #" + std::to_string(vertexId) + " upper link ("
               + std::to_string(upperCount) + " vertices)",
             debug::Priority::VERBOSE);
  }

  if(!lowerCount) {
    // minimum
    return (char)(CriticalType::Local_minimum);
  }
  if(!upperCount) {
    // maximum
    return (char)(CriticalType::Local_maximum);
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
    if(offsets[neighborId0] < offsets[vertexId]
       && offsets[neighborId1] < offsets[vertexId]) {

      // both vertices are lower, let's add that edge and update the UF
      std::map<SimplexId, SimplexId>::iterator n0It
        = global2LowerLink.find(neighborId0);
      std::map<SimplexId, SimplexId>::iterator n1It
        = global2LowerLink.find(neighborId1);

      lowerList[n0It->second] = UnionFind::makeUnion(
        lowerList[n0It->second], lowerList[n1It->second]);
      lowerList[n1It->second] = lowerList[n0It->second];
    }

    // process the upper link
    if(offsets[neighborId0] > offsets[vertexId]
       && offsets[neighborId1] > offsets[vertexId]) {

      // both vertices are lower, let's add that edge and update the UF
      std::map<SimplexId, SimplexId>::iterator n0It
        = global2UpperLink.find(neighborId0);
      std::map<SimplexId, SimplexId>::iterator n1It
        = global2UpperLink.find(neighborId1);

      upperList[n0It->second] = UnionFind::makeUnion(
        upperList[n0It->second], upperList[n1It->second]);
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

  if(debugLevel_ >= (int)(debug::Priority::VERBOSE)) {
    printMsg("Vertex #" + std::to_string(vertexId)
               + ": lowerLink-#CC=" + std::to_string(lowerList.size())
               + " upperLink-#CC=" + std::to_string(upperList.size()),
             debug::Priority::VERBOSE);
  }

  if((lowerList.size() == 1) && (upperList.size() == 1))
    // regular point
    return (char)(CriticalType::Regular);
  else {
    // saddles
    if(dimension_ == 2) {
      if((lowerList.size() > 2) || (upperList.size() > 2)) {
        // monkey saddle
        return (char)(CriticalType::Degenerate);
      } else {
        // regular saddle
        return (char)(CriticalType::Saddle1);
        // NOTE: you may have multi-saddles on the boundary in that
        // configuration
        // to make this computation 100% correct, one would need to disambiguate
        // boundary from interior vertices
      }
    } else if(dimension_ == 3) {
      if((lowerList.size() == 2) && (upperList.size() == 1)) {
        return (char)(CriticalType::Saddle1);
      } else if((lowerList.size() == 1) && (upperList.size() == 2)) {
        return (char)(CriticalType::Saddle2);
      } else {
        // monkey saddle
        return (char)(CriticalType::Degenerate);
        // NOTE: we may have a similar effect in 3D (TODO)
      }
    }
  }

  // -2: regular points
  return (char)(CriticalType::Regular);
}
