#include <Geometry.h>
#include <SurfaceQuadrangulation.h>
#include <cmath>
#include <queue>
#include <set>

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

  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] Common manifolds between vertices";
    for(auto &id : vertsId) {
      msg << " " << id;
    }
    msg << ":";
    for(auto &elem : last) {
      msg << " " << elem;
    }
    msg << std::endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);
  }

  return !last.empty();
}

int ttk::SurfaceQuadrangulation::dualQuadrangulate(
  const std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &sepEdges)
  const {

  // quadrangles vertices are only extrema

  // maps sources (saddle points) to vector of their destinations (extrema)
  std::map<SimplexId, std::vector<SimplexId>> sourceDests{};

  for(auto &p : sepEdges) {
    SimplexId i;
    for(i = 0; i < criticalPointsNumber_; i++) {
      if(p.first == criticalPointsCellIds_[i]) {
        break;
      }
    }
    SimplexId j;
    for(j = 0; j < criticalPointsNumber_; j++) {
      if(p.second == criticalPointsCellIds_[j]) {
        break;
      }
    }

    auto &v = sourceDests[i];
    v.emplace_back(j);
  }

  for(auto &elt : sourceDests) {
    auto extrema = elt.second;
    if(extrema.size() == 4) {
      auto i = extrema[0];
      auto j = i;
      auto k = i;
      auto l = i;
      // filter extrema by nature (minimum: 0 or maximum: 2)
      if(criticalPointsType_[extrema[1]] == criticalPointsType_[i]) {
        j = extrema[2];
        k = extrema[1];
        l = extrema[3];
      } else if(criticalPointsType_[extrema[2]] == criticalPointsType_[i]) {
        j = extrema[1];
        k = extrema[2];
        l = extrema[3];
      } else if(criticalPointsType_[extrema[3]] == criticalPointsType_[i]) {
        j = extrema[2];
        k = extrema[3];
        l = extrema[1];
      }
      outputCells_->emplace_back(4);
      outputCells_->emplace_back(i);
      outputCells_->emplace_back(j);
      outputCells_->emplace_back(k);
      outputCells_->emplace_back(l);
    }
  }

  return 0;
}

int ttk::SurfaceQuadrangulation::quadrangulate(
  const std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &sepEdges,
  size_t &ndegen) const {
  // quadrangle vertices are either extrema or saddle points

  // separatrices: destinations (extrema) -> sources (saddle points)
  std::vector<std::set<SimplexId>> sepMappingDests(sepEdges.size());

  // number of separatrices coming out of this edge
  std::vector<size_t> pointSepNumber(criticalPointsNumber_);

  for(auto &p : sepEdges) {
    for(SimplexId i = 0; i < criticalPointsNumber_; i++) {
      if(p.first == criticalPointsCellIds_[i]) {
        for(SimplexId j = 0; j < criticalPointsNumber_; j++) {
          if(p.second == criticalPointsCellIds_[j]) {
            sepMappingDests[j].insert(i);
            pointSepNumber[j]++;
            // should we also use the saddle points valence?
            pointSepNumber[i]++;
          }
        }
      }
    }
  }

  // iterate twice over dests
  for(size_t i = 0; i < sepMappingDests.size(); i++) {
    // skip if no sources
    if(sepMappingDests[i].empty()) {
      continue;
    }
    // begin second loop at i+1 to avoid duplicates and improve
    // performance
    for(size_t k = i + 1; k < sepMappingDests.size(); k++) {
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
      std::vector<SimplexId> common_dests;
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
            size_t j = common_dests[m];
            size_t l = common_dests[n];

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
        size_t j = common_dests[0];
        ndegen++;

        // fill output vector
        outputCells_->emplace_back(4);
        outputCells_->emplace_back(i);
        outputCells_->emplace_back(j);
        outputCells_->emplace_back(k);
        outputCells_->emplace_back(j);
      }
    }
  }

  return 0;
}

int ttk::SurfaceQuadrangulation::postProcess() const {
  // post-processing: try to detect missing or extra quadrangles by
  // comparing separatrices number coming out of extrema

  // TODO remove extra quadrangles?

  // separatrices bounds indices and cell ids
  std::vector<std::pair<size_t, SimplexId>> sepFlatEdgesPos{};

  for(SimplexId i = 0; i < separatriceNumber_; ++i) {
    if(sepMask_[i] == 1) {
      continue;
    }
    sepFlatEdgesPos.emplace_back(std::make_pair(i, sepCellIds_[i]));
  }

  // duplicate separatrices bounds
  std::vector<
    std::pair<std::pair<size_t, SimplexId>, std::pair<size_t, SimplexId>>>
    dupSep{};

  for(size_t i = 0; i < sepFlatEdgesPos.size() / 2; ++i) {
    for(size_t j = i + 1; j < sepFlatEdgesPos.size() / 2; ++j) {
      if(sepFlatEdgesPos[2 * i].second == sepFlatEdgesPos[2 * j].second
         && sepFlatEdgesPos[2 * i + 1].second
              == sepFlatEdgesPos[2 * j + 1].second) {
        dupSep.emplace_back(
          std::make_pair(sepFlatEdgesPos[2 * i], sepFlatEdgesPos[2 * i + 1]));
        dupSep.emplace_back(
          std::make_pair(sepFlatEdgesPos[2 * j], sepFlatEdgesPos[2 * j + 1]));
      }
    }
  }

  // store sep bounds -> middle id in output points
  std::map<std::pair<size_t, size_t>, size_t> sepMiddles{};

  // generate new points at the middle of duplicate separatrices
  for(auto &p : dupSep) {
    auto a = p.first.first;
    auto b = p.second.first;
    findSeparatrixMiddle(a, b);
    // store separatrix bounds and middle id
    sepMiddles.insert(
      std::make_pair(std::make_pair(a, b), outputPoints_->size() / 3 - 1));
  }

  // for every pair of critical points belonging to a duplicate edge,
  // the indices of the corresponding separatrices in the separatrices array
  std::map<std::pair<SimplexId, SimplexId>,
           std::vector<std::pair<size_t, size_t>>>
    points2Seps{};

  for(const auto &p : sepMiddles) {
    points2Seps[std::make_pair(
                  sepCellIds_[p.first.first], sepCellIds_[p.first.second])]
      .emplace_back(p.first);
  }

  // ad-hoc quad data structure
  struct Quad {
    long long i;
    long long j;
    long long k;
    long long l;
  };

  auto generateNewQuads
    = [&](Quad *const q, const SimplexId v0, const SimplexId v1,
          const SimplexId v2, const SimplexId v3,
          const std::pair<SimplexId, SimplexId> &edge0,
          const std::pair<SimplexId, SimplexId> &edge1) {
        std::vector<SimplexId> mids0{}, mids1{};

        for(const auto &p : points2Seps[edge0]) {
          mids0.emplace_back(sepMiddles[p]);
        }
        for(const auto &p : points2Seps[edge1]) {
          mids1.emplace_back(sepMiddles[p]);
        }

        float *pt0 = &outputPoints_->at(3 * mids0[0]);
        float *pt10 = &outputPoints_->at(3 * mids1[0]);
        float *pt11 = &outputPoints_->at(3 * mids1[1]);

        auto dist00 = Geometry::distance(pt0, pt10);
        auto dist01 = Geometry::distance(pt0, pt11);

        if(dist01 < dist00) {
          std::swap(mids1[0], mids1[1]);
        }

        // first new quad replaces the old one
        q->i = v0;
        q->j = mids0[0];
        q->k = mids1[0];
        q->l = v3;

        // other new quads are placed at the end of the output vector
        outputCells_->emplace_back(4);
        outputCells_->emplace_back(v1);
        outputCells_->emplace_back(mids0[0]);
        outputCells_->emplace_back(mids1[0]);
        outputCells_->emplace_back(v2);
        outputCells_->emplace_back(4);
        outputCells_->emplace_back(v0);
        outputCells_->emplace_back(mids0[1]);
        outputCells_->emplace_back(mids1[1]);
        outputCells_->emplace_back(v3);
        outputCells_->emplace_back(4);
        outputCells_->emplace_back(v1);
        outputCells_->emplace_back(mids0[1]);
        outputCells_->emplace_back(mids1[1]);
        outputCells_->emplace_back(v2);
      };

  auto fixCurrentQuad
    = [&](Quad *const q, const SimplexId v0, const SimplexId v1,
          const SimplexId v2, const SimplexId v3,
          const std::pair<SimplexId, SimplexId> &edge) {
        std::vector<SimplexId> mids{};

        for(const auto &p : points2Seps[edge]) {
          mids.emplace_back(sepMiddles[p]);
        }

        if(outputPoints_ == nullptr) {
          // avoid potential null pointer dereference
          return;
        }

        float *pt0 = &outputPoints_->at(3 * mids[0]);
        float *pt1 = &outputPoints_->at(3 * mids[1]);

        // idea: sort [v0, v1, mid[0] and mid[1]] according to
        // distance to v2 and v3 to choose the best quadrangle with
        // the shortest distance between its vertices

        // distances from v2 to v0, v1, mid[0], mid[1]
        std::vector<float> dists0(mids.size() + 2);
        // distances from v3 to v0, v1, mid[0], mid[1]
        std::vector<float> dists1(mids.size() + 2);

        float *pv0 = &criticalPoints_[3 * v0];
        float *pv1 = &criticalPoints_[3 * v1];
        float *pv2 = &criticalPoints_[3 * v2];
        float *pv3 = &criticalPoints_[3 * v3];

        dists0[0] = Geometry::distance(pv2, pv0);
        dists0[1] = Geometry::distance(pv2, pv1);
        dists0[2] = Geometry::distance(pv2, pt0);
        dists0[3] = Geometry::distance(pv2, pt1);

        dists1[0] = Geometry::distance(pv3, pv0);
        dists1[1] = Geometry::distance(pv3, pv1);
        dists1[2] = Geometry::distance(pv3, pt0);
        dists1[3] = Geometry::distance(pv3, pt1);

        std::vector<float> distssum(dists0.size());
        std::transform(dists0.begin(), dists0.end(), dists1.begin(),
                       distssum.begin(), std::plus<float>());

        // pick one critical point and one sep middle nearest to v2 and v3
        if(distssum[0] < distssum[1]) {
          q->i = v0;
          if(distssum[2] < distssum[3]) {
            q->j = mids[0];
          } else {
            q->j = mids[1];
          }
        } else {
          q->j = v1;
          if(distssum[2] < distssum[3]) {
            q->i = mids[0];
          } else {
            q->i = mids[1];
          }
        }
        q->k = v2;
        q->l = v3;
      };

  // store current number of quads
  auto nquads = outputCells_->size();

  // subdivise quadrangles
  for(size_t i = 0; i < nquads / 5; ++i) {
    auto q = reinterpret_cast<Quad *>(&outputCells_->at(5 * i + 1));
    // bounds indices in separatrices array
    auto qi = criticalPointsCellIds_[q->i];
    auto qj = criticalPointsCellIds_[q->j];
    auto qk = criticalPointsCellIds_[q->k];
    auto ql = criticalPointsCellIds_[q->l];

    // edges (sep sources before dests)
    auto ij = std::make_pair(qj, qi);
    auto jk = std::make_pair(qj, qk);
    auto kl = std::make_pair(ql, qk);
    auto li = std::make_pair(ql, qi);

    auto pse = points2Seps.end();
    if(points2Seps.find(ij) != pse && points2Seps.find(kl) != pse) {
      // from i to l
      generateNewQuads(q, q->i, q->j, q->k, q->l, ij, kl);
    } else if(points2Seps.find(ij) != pse) {
      fixCurrentQuad(q, q->i, q->j, q->k, q->l, ij);
    } else if(points2Seps.find(kl) != pse) {
      fixCurrentQuad(q, q->k, q->l, q->i, q->j, kl);
    }
    if(points2Seps.find(jk) != pse && points2Seps.find(li) != pse) {
      // rotated: from j to i
      generateNewQuads(q, q->j, q->k, q->l, q->i, jk, li);
    } else if(points2Seps.find(jk) != pse) {
      fixCurrentQuad(q, q->j, q->k, q->l, q->i, jk);
    } else if(points2Seps.find(li) != pse) {
      fixCurrentQuad(q, q->l, q->i, q->j, q->k, li);
    }
  }

  return 0;
}

size_t ttk::SurfaceQuadrangulation::findSeparatrixMiddle(const size_t a,
                                                         const size_t b) const {

  const int dim = 3;

  std::vector<float> distFromA(b - a + 1);
  std::array<float, dim> prev{}, curr{};

  curr[0] = sepPoints_[dim * a];
  curr[1] = sepPoints_[dim * a + 1];
  curr[2] = sepPoints_[dim * a + 2];

  // integrate distances at every point of this separatrix
  for(size_t i = 1; i < b - a + 1; ++i) {
    prev = std::move(curr);
    curr[0] = sepPoints_[dim * (a + i)];
    curr[1] = sepPoints_[dim * (a + i) + 1];
    curr[2] = sepPoints_[dim * (a + i) + 2];
    distFromA[i]
      = distFromA[i - 1] + ttk::Geometry::distance(&curr[0], &prev[0]);
  }

  auto distAB = distFromA.back();
  for(auto &el : distFromA) {
    el = std::abs(el - distAB / 2.0);
  }

  // index in separatrices point data array of separatrix middle
  auto id = a + std::min_element(distFromA.begin(), distFromA.end())
            - distFromA.begin();

  // new point!
  outputPoints_->emplace_back(sepPoints_[dim * id]);
  outputPoints_->emplace_back(sepPoints_[dim * id + 1]);
  outputPoints_->emplace_back(sepPoints_[dim * id + 2]);

  // new point identifier (on the triangular mesh)
  switch(sepCellDims_[id]) {
    case 0:
      outputPointsIds_->emplace_back(sepCellIds_[id]);
      break;
    case 1: {
      // take the first vertex of the edge
      SimplexId pos;
      triangulation_->getEdgeVertex(sepCellIds_[id], 0, pos);
      outputPointsIds_->emplace_back(pos);
      break;
    }
    case 2: {
      // take the first vertex of the triangle
      SimplexId pos;
      triangulation_->getTriangleVertex(sepCellIds_[id], 0, pos);
      outputPointsIds_->emplace_back(pos);
      break;
    }
    default:
      break;
  }

  return id;
}

// main routine
int ttk::SurfaceQuadrangulation::execute() const {

  using std::cout;
  using std::endl;
  using std::vector;

  Timer t;

  // clear output
  outputCells_->clear();
  outputPoints_->clear();
  outputPointsIds_->clear();
  outputPoints_->resize(3 * criticalPointsNumber_);
  outputPointsIds_->resize(criticalPointsNumber_);

  // fill in critical points 3d coordinates and identifiers
  for(SimplexId i = 0; i < criticalPointsNumber_; ++i) {
    outputPoints_->at(3 * i) = criticalPoints_[3 * i];
    outputPoints_->at(3 * i + 1) = criticalPoints_[3 * i + 1];
    outputPoints_->at(3 * i + 2) = criticalPoints_[3 * i + 2];
    outputPointsIds_->at(i) = criticalPointsIdentifier_[i];
  }

  // filter sepCellIds_ array according to sepMask_
  std::vector<SimplexId> sepFlatEdges{};

  for(SimplexId i = 0; i < separatriceNumber_; ++i) {
    if(sepMask_[i] == 1) {
      continue;
    }
    sepFlatEdges.emplace_back(sepCellIds_[i]);
  }

  if(sepFlatEdges.size() % 2 != 0) {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] Error: odd number of separatrices edges"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
    return -1;
  }

  // holds separatrices edges for every separatrix
  std::vector<std::pair<SimplexId, SimplexId>> sepEdges{};

  for(size_t i = 0; i < sepFlatEdges.size() / 2; ++i) {
    sepEdges.emplace_back(
      std::make_pair(sepFlatEdges[2 * i], sepFlatEdges[2 * i + 1]));
  }

  // number of degenerate quadrangles
  size_t ndegen = 0;

  if(dualQuadrangulation_) {
    dualQuadrangulate(sepEdges);
  } else {
    // direct quadrangulation with saddle points
    quadrangulate(sepEdges, ndegen);
    postProcess();
  }

  // maximum manifold id
  int maxSegId = 0;
  // compute max + 1 of segmentation indices
  for(size_t a = 0; a < segmentationNumber_; ++a) {
    if(segmentation_[a] > maxSegId) {
      maxSegId = segmentation_[a];
    }
  }
  // total number of manifolds
  int nseg = maxSegId + 1;

  // number of produced quads
  size_t quadNumber = outputCells_->size() / 5;

  // print number of quadrangles wrt number of MSC segmentation
  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] " << quadNumber << " quads (" << ndegen
        << " degenerated, " << nseg << " manifolds)" << endl;
    dMsg(cout, msg.str(), detailedInfoMsg);
  }

  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] Produced " << quadNumber
        << " quadrangles after " << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
