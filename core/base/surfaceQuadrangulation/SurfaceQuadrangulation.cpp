#include <SurfaceQuadrangulation.h>

ttk::SurfaceQuadrangulation::SurfaceQuadrangulation()
  : criticalPointsNumber_{}, criticalPointsCellIds_{}, separatriceNumber_{},
    sepSourceId_{}, sepDestId_{}, outputCells_{} {
}

// main routine
int ttk::SurfaceQuadrangulation::execute() const {

  using std::cout;
  using std::endl;
  using std::vector;

  Timer t;

  // unique separatrices source -> destination mapping
  std::set<std::pair<SimplexId, SimplexId>> sepMappingSet;

  for(SimplexId i = 0; i < separatriceNumber_; i++) {
    sepMappingSet.insert(std::make_pair(sepSourceId_[i], sepDestId_[i]));
  }

  // separatrices sources -> list of destinations mapping
  vector<std::set<SimplexId>> sepMappingSources(sepMappingSet.size());

  for(auto &p : sepMappingSet) {
    for(SimplexId i = 0; i < criticalPointsNumber_; i++) {
      if(p.first == criticalPointsCellIds_[i]) {
        for(SimplexId j = 0; j < criticalPointsNumber_; j++) {
          if(p.second == criticalPointsCellIds_[j]) {
            sepMappingSources[i].insert(j);
          }
        }
      }
    }
  }

  // quadrangle vertices: source, dest, source, dest
  size_t i, j, k, l;

  // iterate twice over sources
  for(i = 0; i < sepMappingSources.size(); i++) {
    // skip if no dests
    if(sepMappingSources[i].empty()) {
      continue;
    }
    for(k = 0; k < sepMappingSources.size(); k++) {
      // skip same source or if no dests
      if(k == i || sepMappingSources[k].empty()) {
        continue;
      }

      // list of common dests to i and k
      vector<SimplexId> common_dests;
      std::set_intersection(
        sepMappingSources[i].begin(), sepMappingSources[i].end(),
        sepMappingSources[k].begin(), sepMappingSources[k].end(),
        std::back_inserter(common_dests));

      // find at least two common dests: j and l
      if(common_dests.size() >= 2) {
        // gotcha!
        j = common_dests[0];
        l = common_dests[1];

        // fill output vector
        outputCells_->emplace_back(4);
        outputCells_->emplace_back(i);
        outputCells_->emplace_back(j);
        outputCells_->emplace_back(k);
        outputCells_->emplace_back(l);
      }
    }
  }

  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
