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

#if 0
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
#endif

  // TODO find quadrangles with separatrices sources <-> dest mappings

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

  {
    stringstream msg;
    msg << "[SurfaceQuadrangulation] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
