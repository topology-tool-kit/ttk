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

            // fill output vector
            outputCells_->emplace_back(4);
            outputCells_->emplace_back(i);
            outputCells_->emplace_back(j);
            outputCells_->emplace_back(k);
            outputCells_->emplace_back(l);
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

  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
