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
  for(size_t i = 0; i < nInputs; i++) {
    std::vector<DiagramTuple> &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      DiagramTuple t = CTDiagram[j];

      ttk::CriticalType nt1 = std::get<1>(t);
      ttk::CriticalType nt2 = std::get<3>(t);

      double dt = std::get<4>(t);
      if(dt > 0) {
        if(nt1 == CriticalType::Local_minimum
           && nt2 == CriticalType::Local_maximum) {
          inputDiagramsMax[i].push_back(t);
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax[i].push_back(t);
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            inputDiagramsMin[i].push_back(t);
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            inputDiagramsSad[i].push_back(t);
          }
        }
      }
    }
  }

  setBidderDiagrams(nInputs, inputDiagramsMin, inputDiagramsSad,
                    inputDiagramsMax, bidder_diagrams_min, bidder_diagrams_sad,
                    bidder_diagrams_max, current_bidder_diagrams_min,
                    current_bidder_diagrams_sad, current_bidder_diagrams_max);

  // Getting current diagrams (with only at most min_points_to_add points)
  std::array<double, 3> max_persistence{};
  std::array<double, 3> lowest_persistence{};
  std::array<double, 3> min_persistence{};

  for(int i_crit = 0; i_crit < 3; i_crit++) {
    max_persistence[i_crit]
      = 2
        * getMostPersistent(i_crit, bidder_diagrams_min, bidder_diagrams_sad,
                            bidder_diagrams_max);
    lowest_persistence[i_crit] = getLessPersistent(
      i_crit, bidder_diagrams_min, bidder_diagrams_sad, bidder_diagrams_max);
    min_persistence[i_crit] = 0;
  }

  std::array<int, 3> min_points_to_add{10, 10, 10};
  std::array<std::vector<double>, 3> min_diag_price{};
  for(auto &arr : min_diag_price) {
    arr.resize(nInputs, 0);
  }

  min_persistence = enrichCurrentBidderDiagrams(
    max_persistence, min_persistence, min_diag_price, min_points_to_add,
    bidder_diagrams_min, bidder_diagrams_sad, bidder_diagrams_max,
    current_bidder_diagrams_min, current_bidder_diagrams_sad,
    current_bidder_diagrams_max);

  const auto distMat = getDiagramsDistMat(
    nInputs, useFullDiagrams_, bidder_diagrams_min, bidder_diagrams_sad,
    bidder_diagrams_max, current_bidder_diagrams_min,
    current_bidder_diagrams_sad, current_bidder_diagrams_max);

  {
    std::stringstream msg;
    msg << "[PersistenceDiagramDistanceMatrix] Processed in "
        << tm.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  return distMat;
}

double PersistenceDiagramDistanceMatrix::getMostPersistent(
  const int type,
  const std::vector<BidderDiagram<double>> &bidder_diags_min,
  const std::vector<BidderDiagram<double>> &bidder_diags_sad,
  const std::vector<BidderDiagram<double>> &bidder_diags_max) const {
  double max_persistence = 0;

  const auto maxPersistence
    = [&](const std::vector<BidderDiagram<double>> bidder_diags) {
        for(unsigned int i = 0; i < bidder_diags.size(); ++i) {
          for(int j = 0; j < bidder_diags[i].size(); ++j) {
            const double persistence = bidder_diags[i].get(j).getPersistence();
            if(persistence > max_persistence) {
              max_persistence = persistence;
            }
          }
        }
      };

  if(do_min_ && (type == -1 || type == 0)) {
    maxPersistence(bidder_diags_min);
  }
  if(do_sad_ && (type == -1 || type == 1)) {
    maxPersistence(bidder_diags_sad);
  }
  if(do_max_ && (type == -1 || type == 2)) {
    maxPersistence(bidder_diags_max);
  }
  return max_persistence;
}

double PersistenceDiagramDistanceMatrix::getLessPersistent(
  const int type,
  const std::vector<BidderDiagram<double>> &bidder_diags_min,
  const std::vector<BidderDiagram<double>> &bidder_diags_sad,
  const std::vector<BidderDiagram<double>> &bidder_diags_max) const {
  // type == -1 : query the min of all the types of diagrams.
  // type = 0 : min,  1 : sad,   2 : max
  // std::cout << "type = " << type << std::endl;

  double min_persistence = std::numeric_limits<double>::max();

  const auto minPersistence
    = [&](const std::vector<BidderDiagram<double>> bidder_diags) {
        for(unsigned int i = 0; i < bidder_diags.size(); ++i) {
          for(int j = 0; j < bidder_diags[i].size(); ++j) {
            const double persistence = bidder_diags[i].get(j).getPersistence();
            if(persistence < min_persistence) {
              min_persistence = persistence;
            }
          }
        }
      };

  if(do_min_ && (type == -1 || type == 0)) {
    minPersistence(bidder_diags_min);
  }
  if(do_sad_ && (type == -1 || type == 1)) {
    minPersistence(bidder_diags_sad);
  }
  if(do_max_ && (type == -1 || type == 2)) {
    minPersistence(bidder_diags_max);
  }
  return min_persistence;
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

std::vector<std::vector<double>>
  PersistenceDiagramDistanceMatrix::getDiagramsDistMat(
    const size_t nInputs,
    const bool useFullDiagrams,
    const std::vector<BidderDiagram<double>> &bidder_diags_min,
    const std::vector<BidderDiagram<double>> &bidder_diags_sad,
    const std::vector<BidderDiagram<double>> &bidder_diags_max,
    const std::vector<BidderDiagram<double>> &current_bidder_diags_min,
    const std::vector<BidderDiagram<double>> &current_bidder_diags_sad,
    const std::vector<BidderDiagram<double>> &current_bidder_diags_max) const {

  std::vector<std::vector<double>> distanceMatrix(nInputs);

  const auto &diags_min
    = useFullDiagrams ? bidder_diags_min : current_bidder_diags_min;
  const auto &diags_saddle
    = useFullDiagrams ? bidder_diags_sad : current_bidder_diags_sad;
  const auto &diags_max
    = useFullDiagrams ? bidder_diags_max : current_bidder_diags_max;

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
        auto &disad = diags_saddle[i];
        auto &djsad = diags_saddle[j];
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

  return distanceMatrix;
}

void PersistenceDiagramDistanceMatrix::setBidderDiagrams(
  const size_t nInputs,
  std::vector<std::vector<DiagramTuple>> &inputDiagramsMin,
  std::vector<std::vector<DiagramTuple>> &inputDiagramsSad,
  std::vector<std::vector<DiagramTuple>> &inputDiagramsMax,
  std::vector<BidderDiagram<double>> &bidder_diags_min,
  std::vector<BidderDiagram<double>> &bidder_diags_sad,
  std::vector<BidderDiagram<double>> &bidder_diags_max,
  std::vector<BidderDiagram<double>> &current_bidder_diags_min,
  std::vector<BidderDiagram<double>> &current_bidder_diags_sad,
  std::vector<BidderDiagram<double>> &current_bidder_diags_max) const {

  const auto setBidderDiag
    = [&](const size_t i, std::vector<std::vector<DiagramTuple>> &inputDiagrams,
          std::vector<BidderDiagram<double>> &bidder_diags,
          std::vector<BidderDiagram<double>> &current_bidder_diags) {
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
            std::cout << "Diagonal point in diagram !!!" << std::endl;
          }
        }
      };

  for(size_t i = 0; i < nInputs; i++) {
    if(do_min_) {
      setBidderDiag(
        i, inputDiagramsMin, bidder_diags_min, current_bidder_diags_min);
    }
    if(do_sad_) {
      setBidderDiag(
        i, inputDiagramsSad, bidder_diags_sad, current_bidder_diags_sad);
    }
    if(do_max_) {
      setBidderDiag(
        i, inputDiagramsMax, bidder_diags_max, current_bidder_diags_max);
    }
  }
}

std::array<double, 3>
  PersistenceDiagramDistanceMatrix::enrichCurrentBidderDiagrams(
    const std::array<double, 3> &previous_min_persistence,
    const std::array<double, 3> &min_persistence,
    const std::array<std::vector<double>, 3> initial_diagonal_prices,
    const std::array<int, 3> min_points_to_add,
    const std::vector<BidderDiagram<double>> &bidder_diags_min,
    const std::vector<BidderDiagram<double>> &bidder_diags_sad,
    const std::vector<BidderDiagram<double>> &bidder_diags_max,
    std::vector<BidderDiagram<double>> &current_bidder_diags_min,
    std::vector<BidderDiagram<double>> &current_bidder_diags_sad,
    std::vector<BidderDiagram<double>> &current_bidder_diags_max) const {

  const auto enrich
    = [&](const double curr_min_persistence, const double prev_min_persistence,
          const std::vector<double> &initial_diag_prices,
          const size_t min_pts_to_add,
          const std::vector<BidderDiagram<double>> &bidder_diags,
          std::vector<BidderDiagram<double>> &current_bidder_diags) {
        auto new_min_persistence = curr_min_persistence;

        // 1. Get size of the largest current diagram, deduce the maximal number
        // of points to append
        const auto nInputs = current_bidder_diags.size();
        size_t max_diagram_size = 0;
        for(const auto &diag : current_bidder_diags) {
          max_diagram_size
            = std::max(static_cast<size_t>(diag.size()), max_diagram_size);
        }
        size_t max_points_to_add
          = std::max(min_pts_to_add, min_pts_to_add + max_diagram_size / 10);
        // 2. Get which points can be added, deduce the new minimal persistence
        std::vector<std::vector<int>> candidates_to_be_added(nInputs);
        std::vector<std::vector<size_t>> idx(nInputs);

        for(size_t i = 0; i < nInputs; i++) {
          double local_min_persistence = std::numeric_limits<double>::min();
          std::vector<double> persistences;
          for(int j = 0; j < bidder_diags[i].size(); j++) {
            Bidder<double> b = bidder_diags[i].get(j);
            double persistence = b.getPersistence();
            if(persistence >= curr_min_persistence
               && persistence <= prev_min_persistence) {
              candidates_to_be_added[i].push_back(j);
              idx[i].push_back(idx[i].size());
              persistences.push_back(persistence);
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
        }
        // 3. Add the points to the current diagrams
        for(size_t i = 0; i < nInputs; i++) {
          const auto s = candidates_to_be_added[i].size();
          for(size_t j = 0; j < std::min(max_points_to_add, s); j++) {
            Bidder<double> b
              = bidder_diags[i].get(candidates_to_be_added[i][idx[i][j]]);
            const double persistence = b.getPersistence();
            if(persistence >= new_min_persistence) {
              b.id_ = current_bidder_diags[i].size();
              b.setPositionInAuction(current_bidder_diags[i].size());
              b.setDiagonalPrice(initial_diag_prices[i]);
              current_bidder_diags[i].addBidder(b);
            }
          }
        }
        return new_min_persistence;
      };

  const std::array<double, 3> new_min_persistence
    = {do_min_ ? enrich(min_persistence[0], previous_min_persistence[0],
                        initial_diagonal_prices[0], min_points_to_add[0],
                        bidder_diags_min, current_bidder_diags_min)
               : previous_min_persistence[0],
       do_sad_ ? enrich(min_persistence[1], previous_min_persistence[1],
                        initial_diagonal_prices[1], min_points_to_add[1],
                        bidder_diags_sad, current_bidder_diags_sad)
               : previous_min_persistence[1],
       do_max_ ? enrich(min_persistence[2], previous_min_persistence[2],
                        initial_diagonal_prices[2], min_points_to_add[2],
                        bidder_diags_max, current_bidder_diags_max)
               : previous_min_persistence[2]};

  return new_min_persistence;
}
