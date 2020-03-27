#include <algorithm>

#include <PersistenceDiagramDistanceMatrix.h>

using namespace ttk;

std::vector<std::vector<double>> PersistenceDiagramDistanceMatrix::execute(
  std::vector<std::vector<DiagramTuple>> &intermediateDiagrams) {

  Timer tm{};

  inputDiagramsMin_.resize(numberOfInputs_);
  inputDiagramsMax_.resize(numberOfInputs_);
  inputDiagramsSaddle_.resize(numberOfInputs_);

  // Create diagrams for min, saddle and max persistence pairs
  for(int i = 0; i < numberOfInputs_; i++) {
    std::vector<DiagramTuple> &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      DiagramTuple t = CTDiagram[j];

      ttk::CriticalType nt1 = std::get<1>(t);
      ttk::CriticalType nt2 = std::get<3>(t);

      double dt = std::get<4>(t);
      // if (abs<double>(dt) < zeroThresh) continue;
      if(dt > 0) {
        if(nt1 == CriticalType::Local_minimum
           && nt2 == CriticalType::Local_maximum) {
          inputDiagramsMax_[i].push_back(t);
          do_max_ = true;
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax_[i].push_back(t);
            do_max_ = true;
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            inputDiagramsMin_[i].push_back(t);
            do_min_ = true;
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            inputDiagramsSaddle_[i].push_back(t);
            do_sad_ = true;
          }
        }
      }
    }
  }

  {
    std::stringstream msg;
    switch(pairTypeClustering_) {
      case(0):
        msg << "[PersistenceDiagramDistanceMatrix] Only MIN-SAD Pairs";
        do_max_ = false;
        do_sad_ = false;
        break;
      case(1):
        msg << "[PersistenceDiagramDistanceMatrix] Only SAD-SAD Pairs";
        do_max_ = false;
        do_min_ = false;
        break;
      case(2):
        msg << "[PersistenceDiagramDistanceMatrix] Only SAD-MAX Pairs";
        do_min_ = false;
        do_sad_ = false;
        break;
      default:
        msg << "[PersistenceDiagramDistanceMatrix] All critical pairs: "
               "global clustering";
        break;
    }
    msg << std::endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);
  }

  original_dos[0] = do_min_;
  original_dos[1] = do_sad_;
  original_dos[2] = do_max_;

  std::array<bool *, 3> current_dos{&do_min_, &do_sad_, &do_max_};

  bool converged = false;
  std::array<bool, 3> diagrams_complete{};
  for(int c = 0; c < 3; c++) {
    diagrams_complete[c] = !original_dos[c];
  }
  bool all_diagrams_complete
    = diagrams_complete[0] && diagrams_complete[1] && diagrams_complete[2];
  n_iterations_ = 0;
  double total_time = 0;

  // double cost = std::numeric_limits<double>::max();
  setBidderDiagrams(inputDiagramsMin_, inputDiagramsSaddle_, inputDiagramsMax_,
                    bidder_diagrams_min_, bidder_diagrams_saddle_,
                    bidder_diagrams_max_, current_bidder_diagrams_min_,
                    current_bidder_diagrams_saddle_,
                    current_bidder_diagrams_max_);

  cost_ = std::numeric_limits<double>::max();
  double min_cost_min = std::numeric_limits<double>::max();
  double min_cost_max = std::numeric_limits<double>::max();
  double min_cost_sad = std::numeric_limits<double>::max();
  // double last_min_cost_obtained = -1;
  double last_min_cost_obtained_min = -1;
  double last_min_cost_obtained_sad = -1;
  double last_min_cost_obtained_max = -1;
  std::array<double, 3> epsilon0{};
  std::array<double, 3> epsilon_candidate{};
  std::array<double, 3> rho{};

  // std::cout<<"checkpoint"<<std::endl;
  // Getting current diagrams (with only at most min_points_to_add points)
  std::array<double, 3> max_persistence{};
  std::array<double, 3> lowest_persistence{};
  std::array<double, 3> min_persistence{};

  // std::cout<<"checkpoint"<<std::endl;
  for(int i_crit = 0; i_crit < 3; i_crit++) {
    // std::cout<<"checkpoint"<<i_crit<<std::endl;
    max_persistence[i_crit]
      = 2
        * getMostPersistent(i_crit, bidder_diagrams_min_,
                            bidder_diagrams_saddle_, bidder_diagrams_max_);
    lowest_persistence[i_crit]
      = getLessPersistent(i_crit, bidder_diagrams_min_, bidder_diagrams_saddle_,
                          bidder_diagrams_max_);
    min_persistence[i_crit] = 0;
    // std::cout<<"size eps "<<epsilon_.size()<<std::endl;
    epsilon_[i_crit] = pow(0.5 * max_persistence[i_crit], 2)
                       / 8.; // max_persistence actually holds 2 times the
                             // highest persistence
    epsilon0[i_crit] = epsilon_[i_crit];
  }
  // std::cout<<"checkpoint"<<std::endl;
  std::array<int, 3> min_points_to_add{};
  min_points_to_add[0] = 10;
  min_points_to_add[1] = 10;
  min_points_to_add[2] = 10;

  std::array<std::vector<double>, 3> min_diag_price{};
  for(auto &arr : min_diag_price) {
    arr.resize(numberOfInputs_, 0);
  }
  // std::cout << "firstinrich" << std::endl;
  if(debugLevel_ > 5) {
    cout << "enrich with rho : " << min_persistence[2]
         << " and initial epsilon " << epsilon_[2]
         << endl; //"  and barycenter size : "<<centroids_max_[0].size()<<endl;
  }
  min_persistence = enrichCurrentBidderDiagrams(
    max_persistence, min_persistence, min_diag_price, min_points_to_add,
    bidder_diagrams_min_, bidder_diagrams_saddle_, bidder_diagrams_max_,
    current_bidder_diagrams_min_, current_bidder_diagrams_saddle_,
    current_bidder_diagrams_max_);

  if(debugLevel_ > 5) {
    std::cout << "first enrich done" << std::endl;
  }

  for(int c = 0; c < 3; c++) {
    if(min_persistence[c] <= lowest_persistence[c]) {
      diagrams_complete[c] = true;
    } // max_persistence actually holds 2 times the highest persistence
  }
  all_diagrams_complete
    = diagrams_complete[0] && diagrams_complete[1] && diagrams_complete[2];

  while(!converged || !all_diagrams_complete) {
    Timer t_inside;

    n_iterations_++;

    for(int i_crit = 0; i_crit < 3; i_crit++) {
      if(*(current_dos[i_crit])) {
        rho[i_crit] = min_persistence[i_crit] > 0
                        ? std::sqrt(8.0 * epsilon_[i_crit])
                        : -1;
      }
    }

    if(n_iterations_ > 1) {
      do_min_ = do_min_ && (min_persistence[0] > rho[0]);
      do_sad_ = do_sad_ && (min_persistence[1] > rho[1]);
      do_max_ = do_max_ && (min_persistence[2] > rho[2]);

      for(int i_crit = 0; i_crit < 3; i_crit++) {
        if(*(current_dos[i_crit])) {
          epsilon_candidate[i_crit] = pow(min_persistence[i_crit], 2) / 8.;
          if(epsilon_candidate[i_crit] > epsilon_[i_crit]) {
            // Should always be the case except if min_persistence is
            // equal to zero
            epsilon_[i_crit] = epsilon_candidate[i_crit];
          }
        }
      }

      if(epsilon_[0] < 5e-5) {
        // Add all remaining points for final convergence.
        min_persistence[0] = 0;
        min_points_to_add[0] = std::numeric_limits<int>::max();
      }
      if(epsilon_[1] < 5e-5) {
        // Add all remaining points for final convergence.
        min_persistence[1] = 0;
        min_points_to_add[1] = std::numeric_limits<int>::max();
      }
      if(epsilon_[2] < 5e-5) {
        // Add all remaining points for final convergence.
        min_persistence[2] = 0;
        min_points_to_add[2] = std::numeric_limits<int>::max();
      }

      if(do_min_ || do_sad_ || do_max_) {
        min_persistence = enrichCurrentBidderDiagrams(
          min_persistence, rho, min_diag_price, min_points_to_add,
          bidder_diagrams_min_, bidder_diagrams_saddle_, bidder_diagrams_max_,
          current_bidder_diagrams_min_, current_bidder_diagrams_saddle_,
          current_bidder_diagrams_max_);
      }

      for(int i_crit = 0; i_crit < 3; i_crit++) {
        if(*(current_dos[i_crit])) {
          if(min_persistence[i_crit] <= lowest_persistence[i_crit]) {
            diagrams_complete[i_crit] = true;
          }
        }
      }

      if(diagrams_complete[0] && diagrams_complete[1] && diagrams_complete[2]) {
        all_diagrams_complete = true;
      }

      resetDosToOriginalValues();
    }

    if(do_min_ && !UseDeltaLim_) {
      precision_min_ = (epsilon_[0] < epsilon0[0] / 500.);
    }
    if(do_sad_ && !UseDeltaLim_) {
      precision_sad_ = (epsilon_[1] < epsilon0[1] / 500.);
    }
    if(do_max_ && !UseDeltaLim_) {
      precision_max_ = (epsilon_[2] < epsilon0[2] / 500.);
    }

    if(epsilon_[0] < epsilon_min_ /*&& diagrams_complete[0]*/) {
      if(debugLevel_ > 4) {
        cout << "[min barycenter] epsilon under minimal value " << endl;
      }
      do_min_ = false;
      epsilon_[0] = epsilon_min_;
      diagrams_complete[0] = true;
    }
    if(epsilon_[1] < epsilon_min_ /*&& diagrams_complete[1]*/) {
      if(debugLevel_ > 4) {
        cout << "[sad barycenter] epsilon under minimal value " << endl;
      }
      do_sad_ = false;
      epsilon_[1] = epsilon_min_;
      diagrams_complete[1] = true;
    }
    if(epsilon_[2] < epsilon_min_ /*&& diagrams_complete[2]*/) {
      if(debugLevel_ > 4) {
        cout << "[max barycenter] epsilon under minimal value " << endl;
      }
      do_max_ = false;
      epsilon_[2] = epsilon_min_;
      diagrams_complete[2] = true;
    }

    if(diagrams_complete[0] && diagrams_complete[1] && diagrams_complete[2]) {
      all_diagrams_complete = true;
    }

    precision_criterion_ = precision_min_ && precision_sad_ && precision_max_;
    bool precision_criterion_reached = precision_criterion_;

    if(cost_min_ < min_cost_min && n_iterations_ > 2
       && diagrams_complete[0] /*&& precision_min_*/) {
      min_cost_min = cost_min_;
      last_min_cost_obtained_min = 0;
    } else if(n_iterations_ > 2 && precision_min_ && diagrams_complete[0]) {
      last_min_cost_obtained_min += 1;
      if(last_min_cost_obtained_min > 1) {
        do_min_ = false;
      }
    }

    if(cost_sad_ < min_cost_sad && n_iterations_ > 2
       && diagrams_complete[1] /*&& precision_sad_*/) {
      min_cost_sad = cost_sad_;
      last_min_cost_obtained_sad = 0;
    } else if(n_iterations_ > 2 && precision_sad_ && diagrams_complete[1]) {
      last_min_cost_obtained_sad += 1;
      if(last_min_cost_obtained_sad > 1 && diagrams_complete[1]) {
        do_sad_ = false;
      }
    }

    if(cost_max_ < min_cost_max && n_iterations_ > 2
       && diagrams_complete[2] /*&& precision_max_*/) {
      min_cost_max = cost_max_;
      last_min_cost_obtained_max = 0;
    } else if(n_iterations_ > 2 && precision_max_ && diagrams_complete[2]) {
      last_min_cost_obtained_max += 1;
      if(last_min_cost_obtained_max > 1 && diagrams_complete[2]) {
        do_max_ = false;
      }
    }

    converged = converged
                || (all_diagrams_complete && !do_min_ && !do_sad_ && !do_max_
                    && (precision_criterion_reached));

    total_time += t_inside.getElapsedTime(); // - t_real_cost.getElapsedTime();

    if(total_time + t_inside.getElapsedTime() > 0.9 * time_limit_) {
      min_cost_min = cost_min_;
      min_cost_sad = cost_sad_;
      min_cost_max = cost_max_;
      converged = true;
    }
    if(total_time > 0.1 * time_limit_) {
      all_diagrams_complete = true;
      diagrams_complete[0] = true;
      diagrams_complete[1] = true;
      diagrams_complete[2] = true;
    }
  }

  resetDosToOriginalValues();

  const auto distMat = getDiagramsDistMat(
    numberOfInputs_, useFullDiagrams_, bidder_diagrams_min_,
    bidder_diagrams_saddle_, bidder_diagrams_max_, current_bidder_diagrams_min_,
    current_bidder_diagrams_saddle_, current_bidder_diagrams_max_);

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
  double delta_lim{0.01};

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

      if(original_dos[0]) {
        auto &dimin = diags_min[i];
        auto &djmin = diags_min[j];
        distance += computeDistance(dimin, djmin, delta_lim);
      }
      if(original_dos[1]) {
        auto &disad = diags_saddle[i];
        auto &djsad = diags_saddle[j];
        distance += computeDistance(disad, djsad, delta_lim);
      }
      if(original_dos[2]) {
        auto &dimax = diags_max[i];
        auto &djmax = diags_max[j];
        distance += computeDistance(dimax, djmax, delta_lim);
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

  for(int i = 0; i < numberOfInputs_; i++) {
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
