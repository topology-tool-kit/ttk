#include <algorithm>
#include <cmath>
#include <cstdlib> /* srand, rand */
#include <iostream>
#include <iterator>
#include <random>

#include <KDTree.h>
#include <PersistenceDiagramDistanceMatrix.h>

using namespace ttk;

void PersistenceDiagramDistanceMatrix::execute(
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

  std::vector<bool *> current_prec{
    &precision_min_, &precision_sad_, &precision_max_};
  std::vector<bool *> current_dos{&do_min_, &do_sad_, &do_max_};

  bool converged = false;
  std::vector<bool> diagrams_complete(3);
  for(int c = 0; c < 3; c++) {
    diagrams_complete[c] = !original_dos[c];
  }
  bool all_diagrams_complete
    = diagrams_complete[0] && diagrams_complete[1] && diagrams_complete[2];
  n_iterations_ = 0;
  double total_time = 0;

  // double cost = std::numeric_limits<double>::max();
  setBidderDiagrams();
  cost_ = std::numeric_limits<double>::max();
  double min_cost_min = std::numeric_limits<double>::max();
  double min_cost_max = std::numeric_limits<double>::max();
  double min_cost_sad = std::numeric_limits<double>::max();
  // double last_min_cost_obtained = -1;
  double last_min_cost_obtained_min = -1;
  double last_min_cost_obtained_sad = -1;
  double last_min_cost_obtained_max = -1;
  std::vector<double> epsilon0(3);
  std::vector<double> epsilon_candidate(3);
  std::vector<double> rho(3);

  // std::cout<<"checkpoint"<<std::endl;
  // Getting current diagrams (with only at most min_points_to_add points)
  std::vector<double> max_persistence(3);
  std::vector<double> lowest_persistence(3);
  std::vector<double> min_persistence(3);

  // std::cout<<"checkpoint"<<std::endl;
  for(int i_crit = 0; i_crit < 3; i_crit++) {
    // std::cout<<"checkpoint"<<i_crit<<std::endl;
    max_persistence[i_crit] = 2 * getMostPersistent(i_crit);
    lowest_persistence[i_crit] = getLessPersistent(i_crit);
    min_persistence[i_crit] = 0;
    // std::cout<<"size eps "<<epsilon_.size()<<std::endl;
    epsilon_[i_crit] = pow(0.5 * max_persistence[i_crit], 2)
                       / 8.; // max_persistence actually holds 2 times the
                             // highest persistence
    epsilon0[i_crit] = epsilon_[i_crit];
  }
  // std::cout<<"checkpoint"<<std::endl;
  std::vector<int> min_points_to_add(3);
  min_points_to_add[0] = 10;
  min_points_to_add[1] = 10;
  min_points_to_add[2] = 10;

  std::vector<std::vector<double>> min_diag_price(3);
  for(int c = 0; c < 3; ++c) {
    for(int i = 0; i < numberOfInputs_; i++) {
      min_diag_price[c].push_back(0);
    }
  }
  // std::cout << "firstinrich" << std::endl;
  if(debugLevel_ > 5) {
    cout << "enrich with rho : " << min_persistence[2]
         << " and initial epsilon " << epsilon_[2]
         << endl; //"  and barycenter size : "<<centroids_max_[0].size()<<endl;
  }
  min_persistence = enrichCurrentBidderDiagrams(
    max_persistence, min_persistence, min_diag_price, min_points_to_add);
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
          min_persistence, rho, min_diag_price, min_points_to_add);
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

  {
    std::stringstream msg;
    msg << "[PersistenceDiagramDistanceMatrix] Processed in "
        << tm.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }
}

double PersistenceDiagramDistanceMatrix::getMostPersistent(int type) {
  double max_persistence = 0;
  // std::cout << "type = " << type << std::endl;
  if(do_min_ && (type == -1 || type == 0)) {
    for(unsigned int i = 0; i < bidder_diagrams_min_.size(); ++i) {
      for(int j = 0; j < bidder_diagrams_min_[i].size(); ++j) {
        Bidder<double> b = bidder_diagrams_min_[i].get(j);
        double persistence = b.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }
    }
  }

  if(do_sad_ && (type == -1 || type == 1)) {
    for(unsigned int i = 0; i < bidder_diagrams_saddle_.size(); ++i) {
      for(int j = 0; j < bidder_diagrams_saddle_[i].size(); ++j) {
        Bidder<double> b = bidder_diagrams_saddle_[i].get(j);
        double persistence = b.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }
    }
  }

  if(do_max_ && (type == -1 || type == 2)) {
    for(unsigned int i = 0; i < bidder_diagrams_max_.size(); ++i) {
      for(int j = 0; j < bidder_diagrams_max_[i].size(); ++j) {
        Bidder<double> b = bidder_diagrams_max_[i].get(j);
        double persistence = b.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }
    }
  }
  return max_persistence;
}

double PersistenceDiagramDistanceMatrix::getLessPersistent(int type) {
  // type == -1 : query the min of all the types of diagrams.
  // type = 0 : min,  1 : sad,   2 : max
  // std::cout << "type = " << type << std::endl;
  double min_persistence = std::numeric_limits<double>::max();
  if(do_min_ && (type == -1 || type == 0)) {
    for(unsigned int i = 0; i < bidder_diagrams_min_.size(); ++i) {
      for(int j = 0; j < bidder_diagrams_min_[i].size(); ++j) {
        Bidder<double> b = bidder_diagrams_min_[i].get(j);
        double persistence = b.getPersistence();
        if(persistence < min_persistence) {
          min_persistence = persistence;
        }
      }
    }
  }

  if(do_sad_ && (type == -1 || type == 1)) {
    for(unsigned int i = 0; i < bidder_diagrams_saddle_.size(); ++i) {
      for(int j = 0; j < bidder_diagrams_saddle_[i].size(); ++j) {
        Bidder<double> b = bidder_diagrams_saddle_[i].get(j);
        double persistence = b.getPersistence();
        if(persistence < min_persistence) {
          min_persistence = persistence;
        }
      }
    }
  }

  if(do_max_ && (type == -1 || type == 2)) {
    for(unsigned int i = 0; i < bidder_diagrams_max_.size(); ++i) {
      for(int j = 0; j < bidder_diagrams_max_[i].size(); ++j) {
        Bidder<double> b = bidder_diagrams_max_[i].get(j);
        double persistence = b.getPersistence();
        if(persistence < min_persistence) {
          min_persistence = persistence;
        }
      }
    }
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
  PersistenceDiagramDistanceMatrix::getDiagramsDistMat() {

  std::vector<std::vector<double>> diagramsDistanceMatrix(numberOfInputs_);
  double delta_lim{0.01};

  const auto &diags_min
    = useFullDiagrams_ ? bidder_diagrams_min_ : current_bidder_diagrams_min_;
  const auto &diags_saddle = useFullDiagrams_ ? bidder_diagrams_saddle_
                                              : current_bidder_diagrams_saddle_;
  const auto &diags_max
    = useFullDiagrams_ ? bidder_diagrams_max_ : current_bidder_diagrams_max_;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < numberOfInputs_; ++i) {
    diagramsDistanceMatrix[i].resize(
      numberOfInputs_, std::numeric_limits<float>::max());

    // matrix diagonal
    diagramsDistanceMatrix[i][i] = 0.0;

    for(int j = i + 1; j < numberOfInputs_; ++j) {
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

      diagramsDistanceMatrix[i][j] = distance;
    }
  }

  // distance matrix is symmetric
  for(int i = 0; i < numberOfInputs_; ++i) {
    for(int j = i + 1; j < numberOfInputs_; ++j) {
      diagramsDistanceMatrix[j][i] = diagramsDistanceMatrix[i][j];
    }
  }

  return diagramsDistanceMatrix;
}

void PersistenceDiagramDistanceMatrix::setBidderDiagrams() {
  for(int i = 0; i < numberOfInputs_; i++) {
    if(do_min_) {
      std::vector<DiagramTuple> *CTDiagram = &(inputDiagramsMin_[i]);
      BidderDiagram<double> bidders;
      for(unsigned int j = 0; j < CTDiagram->size(); j++) {
        // Add bidder to bidders
        Bidder<double> b((*CTDiagram)[j], j, lambda_);

        b.setPositionInAuction(bidders.size());
        bidders.addBidder(b);
        if(b.isDiagonal() || b.x_ == b.y_) {
          std::cout << "Diagonal point in diagram !!!" << std::endl;
        }
      }
      bidder_diagrams_min_.push_back(bidders);
      current_bidder_diagrams_min_.push_back(BidderDiagram<double>());
      std::vector<int> ids(bidders.size());
      for(unsigned int j = 0; j < ids.size(); j++) {
        ids[j] = -1;
      }
    }

    if(do_sad_) {
      std::vector<DiagramTuple> *CTDiagram = &(inputDiagramsSaddle_[i]);

      BidderDiagram<double> bidders;
      for(unsigned int j = 0; j < CTDiagram->size(); j++) {
        // Add bidder to bidders
        Bidder<double> b((*CTDiagram)[j], j, lambda_);

        b.setPositionInAuction(bidders.size());
        bidders.addBidder(b);
        if(b.isDiagonal() || b.x_ == b.y_) {
          std::cout << "Diagonal point in diagram !!!" << std::endl;
        }
      }
      bidder_diagrams_saddle_.push_back(bidders);
      current_bidder_diagrams_saddle_.push_back(BidderDiagram<double>());
      std::vector<int> ids(bidders.size());
      for(unsigned int j = 0; j < ids.size(); j++) {
        ids[j] = -1;
      }
    }

    if(do_max_) {
      std::vector<DiagramTuple> *CTDiagram = &(inputDiagramsMax_[i]);

      BidderDiagram<double> bidders;
      for(unsigned int j = 0; j < CTDiagram->size(); j++) {
        // Add bidder to bidders
        Bidder<double> b((*CTDiagram)[j], j, lambda_);

        b.setPositionInAuction(bidders.size());
        bidders.addBidder(b);
        if(b.isDiagonal() || b.x_ == b.y_) {
          std::cout << "Diagonal point in diagram !!!" << std::endl;
        }
      }
      bidder_diagrams_max_.push_back(bidders);
      current_bidder_diagrams_max_.push_back(BidderDiagram<double>());
      std::vector<int> ids(bidders.size());
      for(unsigned int j = 0; j < ids.size(); j++) {
        ids[j] = -1;
      }
    }
  }
  return;
}

std::vector<double>
  PersistenceDiagramDistanceMatrix::enrichCurrentBidderDiagrams(
    std::vector<double> previous_min_persistence,
    std::vector<double> min_persistence,
    std::vector<std::vector<double>> initial_diagonal_prices,
    std::vector<int> min_points_to_add) {

  std::vector<double> new_min_persistence = min_persistence;

  if(!do_min_) {
    new_min_persistence[0] = previous_min_persistence[0];
  }
  if(!do_sad_) {
    new_min_persistence[1] = previous_min_persistence[1];
  }
  if(!do_max_) {
    new_min_persistence[2] = previous_min_persistence[2];
  }

  // 1. Get size of the largest current diagram, deduce the maximal number of
  // points to append
  int max_diagram_size_min = 0;
  int max_diagram_size_sad = 0;
  int max_diagram_size_max = 0;
  if(do_min_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      if(current_bidder_diagrams_min_[i].size() > max_diagram_size_min) {
        max_diagram_size_min = current_bidder_diagrams_min_[i].size();
      }
    }
  }
  if(do_sad_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      if(current_bidder_diagrams_saddle_[i].size() > max_diagram_size_sad) {
        max_diagram_size_sad = current_bidder_diagrams_saddle_[i].size();
      }
    }
  }
  if(do_max_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      if(current_bidder_diagrams_max_[i].size() > max_diagram_size_max) {
        max_diagram_size_max = current_bidder_diagrams_max_[i].size();
      }
    }
  }
  int max_points_to_add_min = std::max(
    min_points_to_add[0], min_points_to_add[0] + max_diagram_size_min / 10);
  int max_points_to_add_sad = std::max(
    min_points_to_add[1], min_points_to_add[1] + max_diagram_size_sad / 10);
  int max_points_to_add_max = std::max(
    min_points_to_add[2], min_points_to_add[2] + max_diagram_size_max / 10);
  // 2. Get which points can be added, deduce the new minimal persistence
  std::vector<std::vector<int>> candidates_to_be_added_min(numberOfInputs_);
  std::vector<std::vector<int>> candidates_to_be_added_sad(numberOfInputs_);
  std::vector<std::vector<int>> candidates_to_be_added_max(numberOfInputs_);
  std::vector<std::vector<int>> idx_min(numberOfInputs_);
  std::vector<std::vector<int>> idx_sad(numberOfInputs_);
  std::vector<std::vector<int>> idx_max(numberOfInputs_);

  if(do_min_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      double local_min_persistence = std::numeric_limits<double>::min();
      std::vector<double> persistences;
      for(int j = 0; j < bidder_diagrams_min_[i].size(); j++) {
        Bidder<double> b = bidder_diagrams_min_[i].get(j);
        double persistence = b.getPersistence();
        if(persistence >= min_persistence[0]
           && persistence <= previous_min_persistence[0]) {
          candidates_to_be_added_min[i].push_back(j);
          idx_min[i].push_back(idx_min[i].size());
          persistences.push_back(persistence);
        }
      }
      sort(
        idx_min[i].begin(), idx_min[i].end(), [&persistences](int &a, int &b) {
          return ((persistences[a] > persistences[b])
                  || ((persistences[a] == persistences[b]) && (a > b)));
        });
      int size = candidates_to_be_added_min[i].size();
      if(size >= max_points_to_add_min) {
        double last_persistence_added_min
          = persistences[idx_min[i][max_points_to_add_min - 1]];
        if(last_persistence_added_min > local_min_persistence) {
          local_min_persistence = last_persistence_added_min;
        }
      }
      if(i == 0) {
        new_min_persistence[0] = local_min_persistence;
      } else {
        if(local_min_persistence < new_min_persistence[0]) {
          new_min_persistence[0] = local_min_persistence;
        }
      }
    }
  }
  if(do_sad_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      double local_min_persistence = std::numeric_limits<double>::min();
      std::vector<double> persistences;
      for(int j = 0; j < bidder_diagrams_saddle_[i].size(); j++) {
        Bidder<double> b = bidder_diagrams_saddle_[i].get(j);
        double persistence = b.getPersistence();
        if(persistence >= min_persistence[1]
           && persistence <= previous_min_persistence[1]) {
          candidates_to_be_added_sad[i].push_back(j);
          idx_sad[i].push_back(idx_sad[i].size());
          persistences.push_back(persistence);
        }
      }
      sort(
        idx_sad[i].begin(), idx_sad[i].end(), [&persistences](int &a, int &b) {
          return ((persistences[a] > persistences[b])
                  || ((persistences[a] == persistences[b]) && (a > b)));
        });
      int size = candidates_to_be_added_sad[i].size();
      if(size >= max_points_to_add_sad) {
        double last_persistence_added_sad
          = persistences[idx_sad[i][max_points_to_add_sad - 1]];
        if(last_persistence_added_sad > local_min_persistence) {
          local_min_persistence = last_persistence_added_sad;
        }
      }
      if(i == 0) {
        new_min_persistence[1] = local_min_persistence;
      } else {
        if(local_min_persistence < new_min_persistence[1]) {
          new_min_persistence[1] = local_min_persistence;
        }
      }
    }
  }
  if(do_max_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      double local_min_persistence = std::numeric_limits<double>::min();
      std::vector<double> persistences;
      for(int j = 0; j < bidder_diagrams_max_[i].size(); j++) {
        Bidder<double> b = bidder_diagrams_max_[i].get(j);
        double persistence = b.getPersistence();
        if(persistence >= min_persistence[2]
           && persistence <= previous_min_persistence[2]) {
          candidates_to_be_added_max[i].push_back(j);
          idx_max[i].push_back(idx_max[i].size());
          persistences.push_back(persistence);
        }
      }
      sort(
        idx_max[i].begin(), idx_max[i].end(), [&persistences](int &a, int &b) {
          return ((persistences[a] > persistences[b])
                  || ((persistences[a] == persistences[b]) && (a > b)));
        });
      int size = candidates_to_be_added_max[i].size();
      if(size >= max_points_to_add_max) {
        double last_persistence_added_max
          = persistences[idx_max[i][max_points_to_add_max - 1]];
        if(last_persistence_added_max > local_min_persistence) {
          local_min_persistence = last_persistence_added_max;
        }
      }
      if(i == 0) {
        new_min_persistence[2] = local_min_persistence;
      } else {
        if(local_min_persistence < new_min_persistence[2]) {
          new_min_persistence[2] = local_min_persistence;
        }
      }
    }
  }

  // 3. Add the points to the current diagrams
  if(do_min_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      int size = candidates_to_be_added_min[i].size();
      for(int j = 0; j < std::min(max_points_to_add_min, size); j++) {
        Bidder<double> b = bidder_diagrams_min_[i].get(
          candidates_to_be_added_min[i][idx_min[i][j]]);
        double persistence = b.getPersistence();
        if(persistence >= new_min_persistence[0]) {
          b.id_ = current_bidder_diagrams_min_[i].size();
          b.setPositionInAuction(current_bidder_diagrams_min_[i].size());
          b.setDiagonalPrice(initial_diagonal_prices[0][i]);
          current_bidder_diagrams_min_[i].addBidder(b);
        }
      }
    }
  }
  if(do_sad_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      int size = candidates_to_be_added_sad[i].size();
      for(int j = 0; j < std::min(max_points_to_add_sad, size); j++) {
        Bidder<double> b = bidder_diagrams_saddle_[i].get(
          candidates_to_be_added_sad[i][idx_sad[i][j]]);
        double persistence = b.getPersistence();
        if(persistence >= new_min_persistence[1]) {
          b.id_ = current_bidder_diagrams_saddle_[i].size();
          b.setPositionInAuction(current_bidder_diagrams_saddle_[i].size());
          b.setDiagonalPrice(initial_diagonal_prices[1][i]);
          current_bidder_diagrams_saddle_[i].addBidder(b);
        }
      }
    }
  }
  if(do_max_) {
    for(int i = 0; i < numberOfInputs_; i++) {
      int size = candidates_to_be_added_max[i].size();
      for(int j = 0; j < std::min(max_points_to_add_max, size); j++) {
        Bidder<double> b = bidder_diagrams_max_[i].get(
          candidates_to_be_added_max[i][idx_max[i][j]]);
        double persistence = b.getPersistence();
        if(persistence >= new_min_persistence[2]) {
          b.id_ = current_bidder_diagrams_max_[i].size();
          b.setPositionInAuction(current_bidder_diagrams_max_[i].size());
          b.setDiagonalPrice(initial_diagonal_prices[2][i]);
          current_bidder_diagrams_max_[i].addBidder(b);
        }
      }
    }
  }

  return new_min_persistence;
}
