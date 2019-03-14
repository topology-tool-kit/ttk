#include <SurfaceQuadrangulation.h>

ttk::SurfaceQuadrangulation::SurfaceQuadrangulation()
  : vertexNumber_{}, subdivisionLevel_{5}, relaxationIterations_{100},
    criticalPointsNumber_{}, criticalPoints_{}, criticalPointsIdentifiers_{},
    separatriceNumber_{}, sepId_{}, sepSourceId_{}, sepDestId_{} {
}

// main routine
int ttk::SurfaceQuadrangulation::execute() const {

  using std::cout;
  using std::endl;
  using std::map;
  using std::pair;
  using std::set;
  using std::stringstream;
  using std::vector;

  Timer t;

  // unique separatrices source -> destination mapping
  set<pair<SimplexId, SimplexId>> sepMappingSet;

  for(SimplexId i = 0; i < separatriceNumber_; i++) {
    sepMappingSet.insert(std::make_pair(sepSourceId_[i], sepDestId_[i]));
  }

  // separatrices sources -> list of destinations mapping
  vector<vector<SimplexId>> sepMappingSources(sepMappingSet.size());

  for(auto &p : sepMappingSet) {
    for(SimplexId i = 0; i < criticalPointsNumber_; i++) {
      if(p.first == criticalPointsCellIds_[i]) {
        for(SimplexId j = 0; j < criticalPointsNumber_; j++) {
          if(p.second == criticalPointsCellIds_[j]) {
            sepMappingSources[i].emplace_back(j);
          }
        }
      }
    }
  }

  // for each source, iterate over each dest
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sepMappingSources.size(); i++) {
    for(auto &j : sepMappingSources[i]) {
      // find in sources the ones that have j as dest
      size_t k;
      for(k = 0; k < sepMappingSources.size(); k++) {
        if(k == i) {
          continue;
        }
        auto beg = sepMappingSources[k].begin();
        auto end = sepMappingSources[k].end();
        if(std::find(beg, end, j) != end) {
          break;
        }
      }
      // skip this iteration if no k found
      if(k == sepMappingSources.size()) {
        continue;
      }
      // list of common dests to i and k
      vector<SimplexId> common_dests;
      std::set_intersection(
        sepMappingSources[i].begin(), sepMappingSources[i].end(),
        sepMappingSources[k].begin(), sepMappingSources[k].end(),
        std::back_inserter(common_dests));

      // find a fourth vertex that is a common dest to i and k but not j
      if(common_dests.size() > 1) {
        // gotcha!
        SimplexId l;
        for(auto &other_dest : common_dests) {
          if(other_dest == j) {
            // skip second vertex
            continue;
          }
          l = other_dest;
          // we don't need more than four vertices
          break;
        }

        // fill output vector one thread at a time
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif // TTK_ENABLE_OPENMP
        {
          outputCells_->emplace_back(4);
          outputCells_->emplace_back(i);
          outputCells_->emplace_back(j);
          outputCells_->emplace_back(k);
          outputCells_->emplace_back(l);
        }
      }
    }
  }

#if 0
  // print connexions between critical points
  for(SimplexId i = 0; i < criticalPointsNumber_; i++) {
    cout << criticalPointsCellIds_[i] << " ";
    for(auto &v : sepMappingSources[i]) {
      cout << v << " ";
    }
    cout << endl;
  }
  cout << endl;

  // output only segments linking two critical points
  for(auto &p : sepMappingSet) {
    outputCells_->emplace_back(2);
    for(SimplexId i = 0; i < criticalPointsNumber_; i++) {
      if(p.first == criticalPointsCellIds_[i]) {
        outputCells_->emplace_back(i);
        break;
      }
    }
    for(SimplexId i = 0; i < criticalPointsNumber_; i++) {
      if(p.second == criticalPointsCellIds_[i]) {
        outputCells_->emplace_back(i);
        break;
      }
    }
  }
#endif

  {
    stringstream msg;
    msg << "[SurfaceQuadrangulation] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
