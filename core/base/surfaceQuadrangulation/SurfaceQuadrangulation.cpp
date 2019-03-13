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
  // separatrices destinations -> list of sources mapping
  vector<vector<SimplexId>> sepMappingDests(sepMappingSet.size());

  for(auto &p : sepMappingSet) {
    for(SimplexId i = 0; i < criticalPointsNumber_; i++) {
      if(p.first == criticalPointsCellIds_[i]) {
        sepMappingSources[i].emplace_back(p.second);
      }
      if(p.second == criticalPointsCellIds_[i]) {
        sepMappingDests[i].emplace_back(p.first);
      }
    }
  }

  // for each source, iterate over each dest
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sepMappingSources.size(); i++) {
    for(auto &dest : sepMappingSources[i]) {
      // index of dest in the critical points array
      SimplexId j;
      for(j = 0; j < criticalPointsNumber_; j++) {
        if(dest == criticalPointsCellIds_[j]) {
          break;
        }
      }
      // all the sources that lead to dest
      auto dest_sources = sepMappingDests[j];
      for(auto &src : dest_sources) {
        // filter out outer loop source
        if(src == criticalPointsCellIds_[i]) {
          continue;
        }
        // index of the current src
        SimplexId k;
        for(k = 0; k < criticalPointsNumber_; k++) {
          if(src == criticalPointsCellIds_[k]) {
            break;
          }
        }
        // find a second common dest of i and k (fourth vertex)
        vector<SimplexId> common_dests;
        std::set_intersection(
          sepMappingSources[i].begin(), sepMappingSources[i].end(),
          sepMappingSources[k].begin(), sepMappingSources[k].end(),
          std::back_inserter(common_dests));
        if(common_dests.size() > 1) {
          // gotcha!
          SimplexId l;
          for(auto &other_dest : common_dests) {
            if(other_dest == dest) {
              // skip second vertex
              continue;
            }
            // find index of fourth vertex
            for(l = 0; l < criticalPointsNumber_; l++) {
              if(other_dest == criticalPointsCellIds_[l]) {
                break;
              }
            }
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

  for(SimplexId i = 0; i < criticalPointsNumber_; i++) {
    cout << criticalPointsCellIds_[i] << " ";
    for(auto &v : sepMappingDests[i]) {
      cout << v << " ";
    }
    cout << endl;
  }

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
