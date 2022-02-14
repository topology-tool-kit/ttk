#include <MultiresTopology.h>

void ttk::MultiresTopology::getValencesFromLink(
  const SimplexId vertexId,
  const std::vector<std::pair<polarity, polarity>> &vlp,
  DynamicTree &link,
  std::vector<polarity> &toPropageMin,
  std::vector<polarity> &toPropageMax,
  std::vector<std::vector<SimplexId>> &saddleCCMin,
  std::vector<std::vector<SimplexId>> &saddleCCMax) const {

  const auto nbCC = link.getNbCC();

  SimplexId downValence = 0, upValence = 0;
  saddleCCMin[vertexId].clear();
  saddleCCMax[vertexId].clear();

  if(nbCC > 2) {
    std::vector<size_t> CCIds;
    CCIds.reserve(nbCC);
    link.retrieveNbCC(CCIds);
    for(size_t i = 0; i < CCIds.size(); i++) {
      const SimplexId neighbor = CCIds[i];
      const polarity isUpper = vlp[neighbor].first;
      if(isUpper) {
        saddleCCMax[vertexId].emplace_back(neighbor);
        upValence++;
      } else {
        saddleCCMin[vertexId].emplace_back(neighbor);
        downValence++;
      }
    }

    if(downValence > 1) {
      toPropageMin[vertexId] = 255;
    } else {
      saddleCCMin[vertexId].clear();
      toPropageMin[vertexId] = 0;
    }
    if(upValence > 1) {
      toPropageMax[vertexId] = 255;
    } else {
      saddleCCMax[vertexId].clear();
      toPropageMax[vertexId] = 0;
    }

  } else { // not a saddle
    toPropageMax[vertexId] = 0;
    toPropageMin[vertexId] = 0;
  }
}

void ttk::MultiresTopology::buildVertexLinkByBoundary(
  const SimplexId vertexId, VLBoundaryType &vlbt) const {

  const auto bid = multiresTriangulation_.getVertexBoundaryIndex(vertexId);
  const auto nneigh = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  vlbt[bid].reserve(nneigh);

  for(SimplexId i = 0; i < nneigh; i++) {
    SimplexId n0 = 0;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, n0);
    for(SimplexId j = i + 1; j < nneigh; j++) {
      SimplexId n1 = 0;
      multiresTriangulation_.getVertexNeighbor(vertexId, j, n1);
      if(multiresTriangulation_.areVerticesNeighbors(n0, n1)) {
        vlbt[bid].emplace_back(i, j);
      }
    }
  }
}

void ttk::MultiresTopology::getTripletsFromSaddles(
  const SimplexId vertexId,
  std::vector<triplet> &triplets,
  const std::vector<std::vector<SimplexId>> &vertexReps) const {

  const auto &reps = vertexReps[vertexId];
  const SimplexId m = reps[0];
#ifndef TTK_ENABLE_KAMIKAZE
  const auto &repsm = vertexReps[m];
  if(m == -1 || repsm.empty() || repsm[0] != m) {
    this->printErr("HERE PROBLEM");
  }
#endif // TTK_ENABLE_KAMIKAZE
  for(size_t i = 1; i < reps.size(); i++) {
    const SimplexId n = reps[i];
#ifndef TTK_ENABLE_KAMIKAZE
    const auto &repsn = vertexReps[n];
    if(n == -1 || repsn.empty() || repsn[0] != n) {
      this->printErr("HERE2 PROBLEM");
    }
#endif // TTK_ENABLE_KAMIKAZE
    triplets.emplace_back(vertexId, m, n);
  }
}

char ttk::MultiresTopology::getCriticalTypeFromLink(
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

std::string ttk::MultiresTopology::resolutionInfoString() {
  std::stringstream res;
  res << "Resolution level "
      << multiresTriangulation_.DL_to_RL(decimationLevel_);
  if(decimationLevel_ == 0) {
    res << " (final)";
  }
  return res.str();
}

void ttk::MultiresTopology::updateLinkPolarityPonctual(
  std::vector<std::pair<polarity, polarity>> &vlp) const {

  for(size_t i = 0; i < vlp.size(); i++) {
    if(vlp[i].second) {
      vlp[i].first = ~vlp[i].first;
      vlp[i].second = 0;
    }
  }
}
