#include <SurfaceQuadrangulation.h>

ttk::SurfaceQuadrangulation::SurfaceQuadrangulation()
  : criticalPointsNumber_{}, criticalPointsCellIds_{}, separatriceNumber_{},
    sepSourceId_{}, sepDestId_{}, outputCells_{} {
}

bool ttk::SurfaceQuadrangulation::hasCommonManifold(
  const std::vector<size_t> &verts) const {

  // get the manifold of every neighbors of our input vector points
  std::vector<std::set<SimplexId>> common_manifolds(verts.size());

  // TTK identifiers of input vertices
  std::vector<SimplexId> vertsId(verts.size());

  std::transform(verts.begin(), verts.end(), vertsId.begin(),
                 [&](size_t a) { return criticalPointsIdentifier_[a]; });

  // iterate around a neighborhood around the critical points to store
  // their manifold index
  for(size_t i = 0; i < verts.size(); ++i) {

    std::set<SimplexId> processedNeighbors;
    std::queue<SimplexId> neighborsToProcess;

    neighborsToProcess.push(vertsId[i]);

    while(!neighborsToProcess.empty()) {
      auto curr = neighborsToProcess.front();
      neighborsToProcess.pop();
      common_manifolds[i].insert(segmentation_[curr]);
      processedNeighbors.insert(curr);
      auto nneigh = triangulation_->getVertexNeighborNumber(vertsId[i]);
      for(SimplexId j = 0; j < nneigh; ++j) {
        SimplexId next;
        triangulation_->getVertexNeighbor(curr, j, next);
        if(processedNeighbors.find(next) == processedNeighbors.end()) {
          neighborsToProcess.push(next);
        }
      }
      if(processedNeighbors.size() > 20) {
        break;
      }
    }
  }

  // intersect every set to get a unique common manifold
  std::set<SimplexId> curr;
  auto last = common_manifolds[0];

  for(size_t i = 1; i < common_manifolds.size(); ++i) {
    std::set_intersection(last.begin(), last.end(), common_manifolds[i].begin(),
                          common_manifolds[i].end(),
                          std::inserter(curr, curr.begin()));
    std::swap(last, curr);
    curr.clear();
  }

  return !last.empty();
}

// main routine
int ttk::SurfaceQuadrangulation::execute() const {

  using std::cout;
  using std::endl;
  using std::vector;

  Timer t;

  // clear output
  outputCells_->clear();

  // unique separatrices source -> destination mapping
  std::set<std::pair<SimplexId, SimplexId>> sepMappingSet;

  for(SimplexId i = 0; i < separatriceNumber_; i++) {
    sepMappingSet.insert(std::make_pair(sepSourceId_[i], sepDestId_[i]));
  }

  // separatrices: destinations (extrema) -> sources (saddle points)
  vector<std::set<SimplexId>> sepMappingDests(sepMappingSet.size());

  for(auto &p : sepMappingSet) {
    for(SimplexId i = 0; i < criticalPointsNumber_; i++) {
      if(p.first == criticalPointsCellIds_[i]) {
        for(SimplexId j = 0; j < criticalPointsNumber_; j++) {
          if(p.second == criticalPointsCellIds_[j]) {
            sepMappingDests[j].insert(i);
          }
        }
      }
    }
  }

  // quadrangle vertices: dest, source, dest, source
  size_t i, j, k, l;

  // iterate twice over dests
  for(i = 0; i < sepMappingDests.size(); i++) {
    // skip if no sources
    if(sepMappingDests[i].empty()) {
      continue;
    }
    // begin second loop at i+1 to avoid duplicates and improve
    // performance
    for(k = i + 1; k < sepMappingDests.size(); k++) {
      // skip same dest or if no sources
      if(k == i || sepMappingDests[k].empty()) {
        continue;
      }

      // skip if same type (we are looking for quadrangles containing
      // one minimum, one maximum and two saddle points)
      if(criticalPointsType_[k] == criticalPointsType_[i]) {
        continue;
      }

      // list of common sources to i and k
      vector<SimplexId> common_dests;
      std::set_intersection(
        sepMappingDests[i].begin(), sepMappingDests[i].end(),
        sepMappingDests[k].begin(), sepMappingDests[k].end(),
        std::back_inserter(common_dests));

      // find at least two common sources: j and l
      if(common_dests.size() >= 2) {

        // iterate over all possible common sources
        for(size_t m = 0; m < common_dests.size(); m++) {
          // avoid duplicates by beginning at m+1
          for(size_t n = m + 1; n < common_dests.size(); n++) {

            // gotcha!
            j = common_dests[m];
            l = common_dests[n];

            // check for a common shared manifold (looking around
            // saddle points only seems to be sufficient)
            std::vector<size_t> verts{j, l};

            bool validQuad = hasCommonManifold(verts);

            // fill output vector
            if(validQuad) {
              outputCells_->emplace_back(4);
              outputCells_->emplace_back(i);
              outputCells_->emplace_back(j);
              outputCells_->emplace_back(k);
              outputCells_->emplace_back(l);
            }
          }
        }
      } else if(common_dests.size() == 1
                && (sepMappingDests[i].size() == 1
                    || sepMappingDests[k].size() == 1)) {
        // we have degenerate quadrangles: i, j, k, j
        j = common_dests[0];

        // fill output vector
        outputCells_->emplace_back(4);
        outputCells_->emplace_back(i);
        outputCells_->emplace_back(j);
        outputCells_->emplace_back(k);
        outputCells_->emplace_back(j);
      }
    }
  }

  // minimun manifold id
  int minSegId = 0;
  // maximum manifold id
  int maxSegId = 0;
  // compute max + 1 of segmentation indices
  for(size_t i = 0; i < segmentationNumber_; ++i) {
    if(segmentation_[i] > maxSegId) {
      maxSegId = segmentation_[i];
    }
    if(segmentation_[i] < minSegId) {
      minSegId = segmentation_[i];
    }
  }
  // total number of manifolds
  int nseg = maxSegId - minSegId + 1;

  // number of produced quads
  size_t quadNumber = outputCells_->size() / 5;

  // check that we got the right number of quadrangles
  assert(quadNumber == static_cast<size_t>(nseg));

  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] Produced " << quadNumber << " after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
