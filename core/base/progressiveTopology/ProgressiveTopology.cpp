#include <ProgressiveTopology.h>

void ttk::ProgressiveTopology::getValencesFromLink(
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

void ttk::ProgressiveTopology::buildVertexLinkByBoundary(
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

void ttk::ProgressiveTopology::getTripletsFromSaddles(
  const SimplexId vertexId,
  std::vector<triplet> &triplets,
  const std::vector<std::vector<SimplexId>> &vertexReps) const {

  const auto &reps = vertexReps[vertexId];
  const SimplexId m = reps[0];
#ifndef TTK_ENABLE_KAMIKAZE
  const auto &repsm = vertexReps[m];
  if(m == -1 || repsm.empty() || repsm[0] != m) {
    std::cout << "HERE PROBLEM" << std::endl;
  }
#endif // TTK_ENABLE_KAMIKAZE
  for(size_t i = 1; i < reps.size(); i++) {
    const SimplexId n = reps[i];
#ifndef TTK_ENABLE_KAMIKAZE
    const auto &repsn = vertexReps[n];
    if(n == -1 || repsn.empty() || repsn[0] != n) {
      std::cout << "HERE2 PROBLEM" << std::endl;
    }
#endif // TTK_ENABLE_KAMIKAZE
    triplets.emplace_back(vertexId, m, n);
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
      this->printMsg("Computation stopped at decimation level "
                     + std::to_string(this->decimationLevel_));
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
  const SimplexId vertexId,
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
  const double *dummyScalars{};
  ret = executeCPProgressive<double>(0, dummyScalars, offsets);

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

std::string ttk::ProgressiveTopology::resolutionInfoString() {
  std::stringstream res;
  res << "Resolution level "
      << multiresTriangulation_.DL_to_RL(decimationLevel_);
  if(decimationLevel_ == 0) {
    res << " (final)";
  }
  return res.str();
}
