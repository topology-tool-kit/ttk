#include <algorithm>

#include <PersistenceDiagramDistanceMatrix.h>

using namespace ttk;

std::vector<std::vector<double>> PersistenceDiagramDistanceMatrix::execute(
  std::vector<std::vector<DiagramTuple>> &intermediateDiagrams) const {

  Timer tm{};

  const auto nInputs = intermediateDiagrams.size();

  std::vector<std::vector<DiagramTuple>> inputDiagramsMin(nInputs);
  std::vector<std::vector<DiagramTuple>> inputDiagramsSad(nInputs);
  std::vector<std::vector<DiagramTuple>> inputDiagramsMax(nInputs);

  std::vector<BidderDiagram<double>> bidder_diagrams_min{};
  std::vector<BidderDiagram<double>> bidder_diagrams_sad{};
  std::vector<BidderDiagram<double>> bidder_diagrams_max{};
  std::vector<BidderDiagram<double>> current_bidder_diagrams_min{};
  std::vector<BidderDiagram<double>> current_bidder_diagrams_sad{};
  std::vector<BidderDiagram<double>> current_bidder_diagrams_max{};

  // Create diagrams for min, saddle and max persistence pairs
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; i++) {
    std::vector<DiagramTuple> &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      const DiagramTuple t = CTDiagram[j];
      const ttk::CriticalType nt1 = std::get<1>(t);
      const ttk::CriticalType nt2 = std::get<3>(t);
      const double dt = std::get<4>(t);

      if(dt > 0) {
        if(nt1 == CriticalType::Local_minimum
           && nt2 == CriticalType::Local_maximum) {
          inputDiagramsMax[i].emplace_back(t);
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax[i].emplace_back(t);
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            inputDiagramsMin[i].emplace_back(t);
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            inputDiagramsSad[i].emplace_back(t);
          }
        }
      }
    }
  }

  if(do_min_) {
    setBidderDiagrams(nInputs, inputDiagramsMin, bidder_diagrams_min,
                      current_bidder_diagrams_min);
  }
  if(do_sad_) {
    setBidderDiagrams(nInputs, inputDiagramsSad, bidder_diagrams_sad,
                      current_bidder_diagrams_sad);
  }
  if(do_max_) {
    setBidderDiagrams(nInputs, inputDiagramsMax, bidder_diagrams_max,
                      current_bidder_diagrams_max);
  }

  // Getting current diagrams (with only at most min_points_to_add points)
  std::array<double, 3> max_persistence{
    2 * getMostPersistent(bidder_diagrams_min),
    2 * getMostPersistent(bidder_diagrams_sad),
    2 * getMostPersistent(bidder_diagrams_max),
  };

  if(do_min_) {
    enrichCurrentBidderDiagrams(max_persistence[0], min_points_to_add_,
                                bidder_diagrams_min,
                                current_bidder_diagrams_min);
  }
  if(do_sad_) {
    enrichCurrentBidderDiagrams(max_persistence[1], min_points_to_add_,
                                bidder_diagrams_sad,
                                current_bidder_diagrams_sad);
  }
  if(do_max_) {
    enrichCurrentBidderDiagrams(max_persistence[2], min_points_to_add_,
                                bidder_diagrams_max,
                                current_bidder_diagrams_max);
  }

  std::vector<std::vector<double>> distMat(nInputs);
  if(useFullDiagrams_) {
    getDiagramsDistMat(nInputs, distMat, bidder_diagrams_min,
                       bidder_diagrams_sad, bidder_diagrams_max);
  } else {
    getDiagramsDistMat(nInputs, distMat, current_bidder_diagrams_min,
                       current_bidder_diagrams_sad,
                       current_bidder_diagrams_max);
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);

  return distMat;
}

double PersistenceDiagramDistanceMatrix::getMostPersistent(
  const std::vector<BidderDiagram<double>> &bidder_diags) const {

  double max_persistence = 0;

  for(unsigned int i = 0; i < bidder_diags.size(); ++i) {
    for(int j = 0; j < bidder_diags[i].size(); ++j) {
      const double persistence = bidder_diags[i].get(j).getPersistence();
      if(persistence > max_persistence) {
        max_persistence = persistence;
      }
    }
  }

  return max_persistence;
}

double PersistenceDiagramDistanceMatrix::computeDistance(
  const BidderDiagram<double> &D1,
  const BidderDiagram<double> &D2,
  const double delta_lim) const {

  GoodDiagram<double> D2_bis{};
  for(int i = 0; i < D2.size(); i++) {
    const Bidder<double> &b = D2.get(i);
    Good<double> g(b.x_, b.y_, b.isDiagonal(), D2_bis.size());
    g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
    g.setPrice(0);
    D2_bis.addGood(g);
  }

  Auction<double> auction(
    wasserstein_, geometrical_factor_, lambda_, delta_lim, use_kdtree_);
  auction.BuildAuctionDiagrams(&D1, &D2_bis);

  std::vector<MatchingTuple> matchings;
  double cost = auction.run(&matchings);

  return cost;
}

void PersistenceDiagramDistanceMatrix::getDiagramsDistMat(
  const size_t nInputs,
  std::vector<std::vector<double>> &distanceMatrix,
  const std::vector<BidderDiagram<double>> &diags_min,
  const std::vector<BidderDiagram<double>> &diags_sad,
  const std::vector<BidderDiagram<double>> &diags_max) const {

  distanceMatrix.resize(nInputs);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; ++i) {
    distanceMatrix[i].resize(nInputs, std::numeric_limits<float>::max());

    // matrix diagonal
    distanceMatrix[i][i] = 0.0;

    for(size_t j = i + 1; j < nInputs; ++j) {
      double distance{};

      if(do_min_) {
        auto &dimin = diags_min[i];
        auto &djmin = diags_min[j];
        distance += computeDistance(dimin, djmin, deltaLim_);
      }
      if(do_sad_) {
        auto &disad = diags_sad[i];
        auto &djsad = diags_sad[j];
        distance += computeDistance(disad, djsad, deltaLim_);
      }
      if(do_max_) {
        auto &dimax = diags_max[i];
        auto &djmax = diags_max[j];
        distance += computeDistance(dimax, djmax, deltaLim_);
      }

      distanceMatrix[i][j] = distance;
    }
  }

  // distance matrix is symmetric
  for(size_t i = 0; i < nInputs; ++i) {
    for(size_t j = i + 1; j < nInputs; ++j) {
      distanceMatrix[j][i] = distanceMatrix[i][j];
    }
  }
}

void PersistenceDiagramDistanceMatrix::setBidderDiagrams(
  const size_t nInputs,
  std::vector<std::vector<DiagramTuple>> &inputDiagrams,
  std::vector<BidderDiagram<double>> &bidder_diags,
  std::vector<BidderDiagram<double>> &current_bidder_diags) const {

  for(size_t i = 0; i < nInputs; i++) {
    bidder_diags.emplace_back();
    current_bidder_diags.emplace_back();

    auto &diag = inputDiagrams[i];
    auto &bidders = bidder_diags.back();

    for(size_t j = 0; j < diag.size(); j++) {
      // Add bidder to bidders
      Bidder<double> b(diag[j], j, lambda_);
      b.setPositionInAuction(bidders.size());
      bidders.addBidder(b);
      if(b.isDiagonal() || b.x_ == b.y_) {
        this->printWrn("Diagonal point in diagram!");
      }
    }
  }
}

double PersistenceDiagramDistanceMatrix::enrichCurrentBidderDiagrams(
  const double prev_min_persistence,
  const size_t min_points_to_add,
  const std::vector<BidderDiagram<double>> &bidder_diags,
  std::vector<BidderDiagram<double>> &current_bidder_diags) const {

  double new_min_persistence = 0.0;

  // 1. Get size of the largest current diagram, deduce the maximal number
  // of points to append
  const auto nInputs = current_bidder_diags.size();
  size_t max_diagram_size = 0;
  for(const auto &diag : current_bidder_diags) {
    max_diagram_size
      = std::max(static_cast<size_t>(diag.size()), max_diagram_size);
  }
  size_t max_points_to_add
    = std::max(min_points_to_add, min_points_to_add + max_diagram_size / 10);
  // 2. Get which points can be added, deduce the new minimal persistence
  std::vector<std::vector<int>> candidates_to_be_added(nInputs);
  std::vector<std::vector<size_t>> idx(nInputs);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; i++) {
    double local_min_persistence = std::numeric_limits<double>::min();
    std::vector<double> persistences;
    for(int j = 0; j < bidder_diags[i].size(); j++) {
      Bidder<double> b = bidder_diags[i].get(j);
      double persistence = b.getPersistence();
      if(persistence >= 0.0 && persistence <= prev_min_persistence) {
        candidates_to_be_added[i].emplace_back(j);
        idx[i].emplace_back(idx[i].size());
        persistences.emplace_back(persistence);
      }
    }
    const auto cmp = [&persistences](const size_t a, const size_t b) {
      return ((persistences[a] > persistences[b])
              || ((persistences[a] == persistences[b]) && (a > b)));
    };
    std::sort(idx[i].begin(), idx[i].end(), cmp);
    const auto size = candidates_to_be_added[i].size();
    if(size >= max_points_to_add) {
      double last_persistence_added
        = persistences[idx[i][max_points_to_add - 1]];
      if(last_persistence_added > local_min_persistence) {
        local_min_persistence = last_persistence_added;
      }
    }
    if(i == 0) {
      new_min_persistence = local_min_persistence;
    } else {
      if(local_min_persistence < new_min_persistence) {
        new_min_persistence = local_min_persistence;
      }
    }
    // 3. Add the points to the current diagrams
    const auto s = candidates_to_be_added[i].size();
    for(size_t j = 0; j < std::min(max_points_to_add, s); j++) {
      Bidder<double> b
        = bidder_diags[i].get(candidates_to_be_added[i][idx[i][j]]);
      const double persistence = b.getPersistence();
      if(persistence >= new_min_persistence) {
        b.id_ = current_bidder_diags[i].size();
        b.setPositionInAuction(current_bidder_diags[i].size());
        current_bidder_diags[i].addBidder(b);
      }
    }
  }
  return new_min_persistence;
}
