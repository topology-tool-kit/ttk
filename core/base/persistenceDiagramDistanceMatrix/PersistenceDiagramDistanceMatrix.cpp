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

  {
    std::stringstream msg;
    msg << "[PersistenceDiagramDistanceMatrix] Clustering " << numberOfInputs_
        << " diagrams in " << k_ << " cluster(s)." << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

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
    diagrams_complete[c] = (!use_progressive_) || (!original_dos[c]);
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

  if(use_progressive_) {
    // min_persistence = max_persistence/2.;
    // min_persistence = 0;
  } else {
    min_points_to_add[0] = std::numeric_limits<int>::max();
    min_points_to_add[1] = std::numeric_limits<int>::max();
    min_points_to_add[2] = std::numeric_limits<int>::max();
  }
  std::vector<std::vector<double>> min_diag_price(3);
  std::vector<std::vector<double>> min_off_diag_price(3);
  for(int c = 0; c < 3; ++c) {
    for(int i = 0; i < numberOfInputs_; i++) {
      min_diag_price[c].push_back(0);
      min_off_diag_price[c].push_back(0);
    }
  }
  // std::cout << "firstinrich" << std::endl;
  if(debugLevel_ > 5) {
    cout << "enrich with rho : " << min_persistence[2]
         << " and initial epsilon " << epsilon_[2]
         << endl; //"  and barycenter size : "<<centroids_max_[0].size()<<endl;
  }
  min_persistence = enrichCurrentBidderDiagrams(
    max_persistence, min_persistence, min_diag_price, min_off_diag_price,
    min_points_to_add, false);
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
  if(all_diagrams_complete) {
    use_progressive_ = false;
  }

  // Initializing centroids and clusters
  if(use_kmeanspp_) {
    initializeCentroidsKMeanspp();
  } else {
    initializeCentroids();
  }

  initializeEmptyClusters();

  if(use_accelerated_) {
    initializeAcceleratedKMeans();
    getCentroidDistanceMatrix();
    acceleratedUpdateClusters();
  } else {
    updateClusters();
    old_clustering_ = clustering_;
  }

  while(!converged || (!all_diagrams_complete && use_progressive_)) {
    Timer t_inside;

    n_iterations_++;

    for(int i_crit = 0; i_crit < 3; i_crit++) {
      if(*(current_dos[i_crit])) {
        rho[i_crit] = min_persistence[i_crit] > 0
                        ? std::sqrt(8.0 * epsilon_[i_crit])
                        : -1;
      }
    }

    if(use_progressive_ && n_iterations_ > 1) {
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
          min_persistence, rho, min_diag_price, min_off_diag_price,
          min_points_to_add, true);
      }
      barycenter_inputs_reset_flag = true;

      for(int i_crit = 0; i_crit < 3; i_crit++) {
        if(*(current_dos[i_crit])) {
          if(min_persistence[i_crit] <= lowest_persistence[i_crit]) {
            diagrams_complete[i_crit] = true;
          }
        }
      }

      if(diagrams_complete[0] && diagrams_complete[1] && diagrams_complete[2]) {
        use_progressive_ = false;
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
      use_progressive_ = false;
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
      use_progressive_ = false;
    }
  }

  resetDosToOriginalValues();

  computeDiagramsDistanceMatrix();

  {
    std::stringstream msg;
    msg << "[PersistenceDiagramDistanceMatrix] Processed in "
        << tm.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }
}

vector<vector<int>> PersistenceDiagramDistanceMatrix::get_centroids_sizes() {
  return centroids_sizes_;
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

std::vector<std::vector<double>>
  PersistenceDiagramDistanceMatrix::getMinPrices() {
  std::vector<std::vector<double>> min_prices(3);
  // cout<<"dos : "<<original_dos[0]<<original_dos[1]<<original_dos[2]<<endl;
  if(original_dos[0]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[0].push_back(std::numeric_limits<double>::max());
      for(int j = 0; j < centroids_with_price_min_[i].size(); ++j) {
        Good<double> g = centroids_with_price_min_[i].get(j);
        double price = g.getPrice();
        if(price < min_prices[0][i]) {
          min_prices[0][i] = price;
        }
      }
    }
  }

  if(original_dos[1]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[1].push_back(std::numeric_limits<double>::max());
      for(int j = 0; j < centroids_with_price_saddle_[i].size(); ++j) {
        Good<double> g = centroids_with_price_saddle_[i].get(j);
        double price = g.getPrice();
        if(price < min_prices[1][i]) {
          min_prices[1][i] = price;
        }
      }
    }
  }

  if(original_dos[2]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[2].push_back(std::numeric_limits<double>::max());
      for(int j = 0; j < centroids_with_price_max_[i].size(); ++j) {
        Good<double> g = centroids_with_price_max_[i].get(j);
        double price = g.getPrice();
        if(price < min_prices[2][i]) {
          min_prices[2][i] = price;
        }
      }
    }
  }

  return min_prices;
}

std::vector<std::vector<double>>
  PersistenceDiagramDistanceMatrix::getMinDiagonalPrices() {
  std::vector<std::vector<double>> min_prices(3);
  if(original_dos[0]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[0].push_back(std::numeric_limits<double>::max());
      for(int j = 0; j < current_bidder_diagrams_min_[i].size(); ++j) {
        Bidder<double> b = current_bidder_diagrams_min_[i].get(j);
        double price = b.diagonal_price_;
        if(price < min_prices[0][i]) {
          min_prices[0][i] = price;
        }
      }
      if(min_prices[0][i] >= std::numeric_limits<double>::max() / 2.) {
        min_prices[0][i] = 0;
      }
    }
  }

  if(original_dos[1]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[1].push_back(std::numeric_limits<double>::max());
      for(int j = 0; j < current_bidder_diagrams_saddle_[i].size(); ++j) {
        Bidder<double> b = current_bidder_diagrams_saddle_[i].get(j);
        double price = b.diagonal_price_;
        if(price < min_prices[1][i]) {
          min_prices[1][i] = price;
        }
      }
      if(min_prices[1][i] >= std::numeric_limits<double>::max() / 2.) {
        min_prices[1][i] = 0;
      }
    }
  }

  if(original_dos[2]) {
    for(int i = 0; i < numberOfInputs_; ++i) {
      min_prices[2].push_back(std::numeric_limits<double>::max());
      for(int j = 0; j < current_bidder_diagrams_max_[i].size(); ++j) {
        Bidder<double> b = current_bidder_diagrams_max_[i].get(j);
        double price = b.diagonal_price_;
        if(price < min_prices[2][i]) {
          min_prices[2][i] = price;
        }
      }
      if(min_prices[2][i] >= std::numeric_limits<double>::max() / 2.) {
        min_prices[2][i] = 0;
      }
    }
  }
  return min_prices;
}

double PersistenceDiagramDistanceMatrix::computeDistance(
  const BidderDiagram<double> &D1,
  const BidderDiagram<double> &D2,
  const double delta_lim) {
  GoodDiagram<double> D2_bis = diagramToCentroid(D2);
  return computeDistance(D1, D2_bis, delta_lim);
}

double PersistenceDiagramDistanceMatrix::computeDistance(
  const BidderDiagram<double> D1,
  const GoodDiagram<double> D2,
  const double delta_lim) {
  std::vector<MatchingTuple> matchings;
  const auto D2_bis = centroidWithZeroPrices(D2);
  Auction<double> auction(
    wasserstein_, geometrical_factor_, lambda_, delta_lim, use_kdtree_);
  auction.BuildAuctionDiagrams(&D1, &D2_bis);
  double cost = auction.run(&matchings);
  return cost;
}

double PersistenceDiagramDistanceMatrix::computeDistance(
  BidderDiagram<double> *const D1,
  const GoodDiagram<double> *const D2,
  const double delta_lim) {
  std::vector<MatchingTuple> matchings;
  Auction<double> auction(
    wasserstein_, geometrical_factor_, lambda_, delta_lim, use_kdtree_);
  int size1 = D1->size();
  auction.BuildAuctionDiagrams(D1, D2);
  double cost = auction.run(&matchings);
  // Diagonal Points were added in the original diagram. The following line
  // removes them.
  D1->bidders_.resize(size1);
  return cost;
}

double PersistenceDiagramDistanceMatrix::computeDistance(
  const GoodDiagram<double> &D1,
  const GoodDiagram<double> &D2,
  const double delta_lim) {
  BidderDiagram<double> D1_bis = centroidToDiagram(D1);
  return computeDistance(D1_bis, D2, delta_lim);
}

GoodDiagram<double> PersistenceDiagramDistanceMatrix::centroidWithZeroPrices(
  const GoodDiagram<double> centroid) {
  GoodDiagram<double> GD = GoodDiagram<double>();
  for(int i = 0; i < centroid.size(); i++) {
    Good<double> g = centroid.get(i);
    g.setPrice(0);
    GD.addGood(g);
  }
  return GD;
}

BidderDiagram<double> PersistenceDiagramDistanceMatrix::diagramWithZeroPrices(
  const BidderDiagram<double> diagram) {
  BidderDiagram<double> BD = BidderDiagram<double>();
  for(int i = 0; i < diagram.size(); i++) {
    Bidder<double> b = diagram.get(i);
    b.setDiagonalPrice(0);
    BD.addBidder(b);
  }
  return BD;
}

BidderDiagram<double> PersistenceDiagramDistanceMatrix::centroidToDiagram(
  const GoodDiagram<double> centroid) {
  BidderDiagram<double> BD = BidderDiagram<double>();
  for(int i = 0; i < centroid.size(); i++) {
    Good<double> g = centroid.get(i);

    Bidder<double> b = Bidder<double>(g.x_, g.y_, g.isDiagonal(), BD.size());
    b.SetCriticalCoordinates(g.coords_x_, g.coords_y_, g.coords_z_);
    b.setPositionInAuction(BD.size());
    BD.addBidder(b);
  }
  return BD;
}

GoodDiagram<double> PersistenceDiagramDistanceMatrix::diagramToCentroid(
  const BidderDiagram<double> diagram) {
  GoodDiagram<double> GD = GoodDiagram<double>();
  for(int i = 0; i < diagram.size(); i++) {
    Bidder<double> b = diagram.get(i);

    Good<double> g = Good<double>(b.x_, b.y_, b.isDiagonal(), GD.size());
    g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
    GD.addGood(g);
  }
  return GD;
}

void PersistenceDiagramDistanceMatrix::initializeEmptyClusters() {
  clustering_ = std::vector<std::vector<int>>(k_);
}

void PersistenceDiagramDistanceMatrix::initializeCentroids() {
  std::vector<int> idx(numberOfInputs_);
  // To perform a random draw with replacement, the vector {1, 2, ...,
  // numberOfInputs_} is shuffled, and we consider its k_ first elements to be
  // the initial centroids.
  for(int i = 0; i < numberOfInputs_; i++) {
    idx[i] = i;
  }
  if(!deterministic_)
    std::random_shuffle(idx.begin(), idx.end());

  for(int c = 0; c < k_; c++) {
    if(do_min_) {
      GoodDiagram<double> centroid_min
        = diagramToCentroid(current_bidder_diagrams_min_[idx[c]]);
      centroids_min_.push_back(centroid_min);
    }
    if(do_sad_) {
      GoodDiagram<double> centroid_sad
        = diagramToCentroid(current_bidder_diagrams_saddle_[idx[c]]);
      centroids_saddle_.push_back(centroid_sad);
    }
    if(do_max_) {
      GoodDiagram<double> centroid_max
        = diagramToCentroid(current_bidder_diagrams_max_[idx[c]]);
      centroids_max_.push_back(centroid_max);
    }
  }
}

void PersistenceDiagramDistanceMatrix::initializeCentroidsKMeanspp() {
  std::vector<int> indexes_clusters;
  int random_idx = deterministic_ ? 0 : rand() % numberOfInputs_;
  indexes_clusters.push_back(random_idx);

  if(do_min_) {
    GoodDiagram<double> centroid_min
      = diagramToCentroid(current_bidder_diagrams_min_[random_idx]);
    centroids_min_.push_back(centroid_min);
  }
  if(do_sad_) {
    GoodDiagram<double> centroid_sad
      = diagramToCentroid(current_bidder_diagrams_saddle_[random_idx]);
    centroids_saddle_.push_back(centroid_sad);
  }
  if(do_max_) {
    GoodDiagram<double> centroid_max
      = diagramToCentroid(current_bidder_diagrams_max_[random_idx]);
    centroids_max_.push_back(centroid_max);
  }
  // cout<<"CP 0. sizes of bidders : "<<current_bidder_diagrams_min_.size()<<"
  // "<<current_bidder_diagrams_max_.size()<<endl;
  while((int)indexes_clusters.size() < k_) {
    std::vector<double> min_distance_to_centroid(numberOfInputs_);
    std::vector<double> probabilities(numberOfInputs_);

    // Uncomment for a deterministic algorithm
    double maximal_distance = 0;
    int candidate_centroid = 0;

    for(int i = 0; i < numberOfInputs_; i++) {
      // cout<<"test1"<<i<<endl;
      min_distance_to_centroid[i] = std::numeric_limits<double>::max();
      if(std::find(indexes_clusters.begin(), indexes_clusters.end(), i)
         != indexes_clusters.end()) {
        // cout<<"go 0"<<endl;
        min_distance_to_centroid[i] = 0;
      } else {
        // cout<<"go 1"<<endl;
        for(unsigned int j = 0; j < indexes_clusters.size(); ++j) {
          // cout<<"test "<<j<<" sizes :
          // "<<current_bidder_diagrams_min_.size()<<"
          // "<<centroids_min_.size()<<endl;
          double distance = 0;
          if(do_min_) {
            // cout<<"1"<<endl;
            GoodDiagram<double> centroid_min
              = centroidWithZeroPrices(centroids_min_[j]);
            // cout<<"2"<<endl;
            distance += computeDistance(
              current_bidder_diagrams_min_[i], centroid_min, 0.01);
            // cout<<"3"<<endl;
          }
          // cout<<"test "<<j<<endl;
          if(do_sad_) {
            GoodDiagram<double> centroid_saddle
              = centroidWithZeroPrices(centroids_saddle_[j]);
            distance += computeDistance(
              current_bidder_diagrams_saddle_[i], centroid_saddle, 0.01);
          }
          // cout<<"test "<<j<<endl;
          if(do_max_) {
            GoodDiagram<double> centroid_max
              = centroidWithZeroPrices(centroids_max_[j]);
            distance += computeDistance(
              current_bidder_diagrams_max_[i], centroid_max, 0.01);
          }
          // cout<<"test "<<j<<endl;
          if(distance < min_distance_to_centroid[i]) {
            min_distance_to_centroid[i] = distance;
          }
          // cout<<"test "<<j<<endl;
        }
      }
      probabilities[i] = pow(min_distance_to_centroid[i], 2);

      // The following block is useful in case of need for a deterministic
      // algoritm
      if(deterministic_ && min_distance_to_centroid[i] > maximal_distance) {
        maximal_distance = min_distance_to_centroid[i];
        candidate_centroid = i;
      }
    }
    // cout<<"CP 1 "<<candidate_centroid<<endl;
    // Comment the following four lines to make it deterministic
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<int> distribution(
      probabilities.begin(), probabilities.end());

    if(!deterministic_) {
      candidate_centroid = distribution(gen);
    }

    indexes_clusters.push_back(candidate_centroid);
    if(do_min_) {
      GoodDiagram<double> centroid_min
        = diagramToCentroid(current_bidder_diagrams_min_[candidate_centroid]);
      centroids_min_.push_back(centroid_min);
    }
    if(do_sad_) {
      GoodDiagram<double> centroid_sad = diagramToCentroid(
        current_bidder_diagrams_saddle_[candidate_centroid]);
      centroids_saddle_.push_back(centroid_sad);
    }
    if(do_max_) {
      GoodDiagram<double> centroid_max
        = diagramToCentroid(current_bidder_diagrams_max_[candidate_centroid]);
      centroids_max_.push_back(centroid_max);
    }
  }
}

void PersistenceDiagramDistanceMatrix::initializeAcceleratedKMeans() {
  // r_ is a vector stating for each diagram if its distance to its centroid is
  // up to date (false) or needs to be recomputed (true)
  r_ = std::vector<bool>(numberOfInputs_);
  // u_ is a vector of upper bounds of the distance of each diagram to its
  // closest centroid
  u_ = std::vector<double>(numberOfInputs_);
  inv_clustering_ = std::vector<int>(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; i++) {
    r_[i] = true;
    u_[i] = std::numeric_limits<double>::max();
    inv_clustering_[i] = -1;
  }
  // l_ is the matrix of lower bounds for the distance from each diagram
  // to each centroid
  l_ = std::vector<std::vector<double>>(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; ++i) {
    l_[i] = std::vector<double>(k_);
    for(int c = 0; c < k_; ++c) {
      l_[i][c] = 0;
    }
  }

  // And d_ is a (K x K) matrix storing the distances between each pair of
  // centroids
  centroidsDistanceMatrix_.resize(k_);
  for(int i = 0; i < k_; ++i) {
    centroidsDistanceMatrix_[i].resize(k_, 0.0);
  }
  return;
}

std::vector<std::vector<double>>
  PersistenceDiagramDistanceMatrix::getDistanceMatrix() {
  std::vector<std::vector<double>> D(numberOfInputs_);

  for(int i = 0; i < numberOfInputs_; ++i) {
    BidderDiagram<double> D1_min, D1_sad, D1_max;
    if(do_min_) {
      D1_min = diagramWithZeroPrices(current_bidder_diagrams_min_[i]);
    }
    if(do_sad_) {
      D1_sad = diagramWithZeroPrices(current_bidder_diagrams_saddle_[i]);
    }
    if(do_max_) {
      D1_max = diagramWithZeroPrices(current_bidder_diagrams_max_[i]);
    }
    for(int c = 0; c < k_; ++c) {
      GoodDiagram<double> D2_min, D2_sad, D2_max;
      double distance = 0;
      if(do_min_) {
        D2_min = centroids_min_[c];
        distance += computeDistance(D1_min, D2_min, 0.01);
      }
      if(do_sad_) {
        D2_sad = centroids_saddle_[c];
        distance += computeDistance(D1_sad, D2_sad, 0.01);
      }
      if(do_max_) {
        D2_max = centroids_max_[c];
        distance += computeDistance(D1_max, D2_max, 0.01);
      }
      D[i].push_back(distance);
    }
  }
  return D;
}

void PersistenceDiagramDistanceMatrix::getCentroidDistanceMatrix() {
  for(int i = 0; i < k_; ++i) {
    GoodDiagram<double> D1_min, D1_sad, D1_max;
    if(do_min_) {
      D1_min = centroidWithZeroPrices(centroids_min_[i]);
    }
    if(do_sad_) {
      D1_sad = centroidWithZeroPrices(centroids_saddle_[i]);
    }
    if(do_max_) {
      D1_max = centroidWithZeroPrices(centroids_max_[i]);
    }
    for(int j = i + 1; j < k_; ++j) {
      double distance{};
      GoodDiagram<double> D2_min, D2_sad, D2_max;
      if(do_min_) {
        D2_min = centroidWithZeroPrices(centroids_min_[j]);
        distance += computeDistance(D1_min, D2_min, 0.01);
      }
      if(do_sad_) {
        D2_sad = centroidWithZeroPrices(centroids_saddle_[j]);
        distance += computeDistance(D1_sad, D2_sad, 0.01);
      }
      if(do_max_) {
        D2_max = centroidWithZeroPrices(centroids_max_[j]);
        distance += computeDistance(D1_max, D2_max, 0.01);
      }

      centroidsDistanceMatrix_[i][j] = distance;
      centroidsDistanceMatrix_[j][i] = distance;
    }
  }
  return;
}

void PersistenceDiagramDistanceMatrix::computeDiagramsDistanceMatrix() {

  diagramsDistanceMatrix_.resize(numberOfInputs_);
  double delta_lim{0.01};

  const auto &diags_min
    = useFullDiagrams_ ? bidder_diagrams_min_ : current_bidder_diagrams_min_;
  const auto &diags_saddle = useFullDiagrams_ ? bidder_diagrams_saddle_
                                              : current_bidder_diagrams_saddle_;
  const auto &diags_max
    = useFullDiagrams_ ? bidder_diagrams_max_ : current_bidder_diagrams_max_;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < numberOfInputs_; ++i) {
    diagramsDistanceMatrix_[i].resize(
      numberOfInputs_, std::numeric_limits<float>::max());

    // matrix diagonal
    diagramsDistanceMatrix_[i][i] = 0.0;

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

      diagramsDistanceMatrix_[i][j] = distance;
    }
  }

  // distance matrix is symmetric
  for(int i = 0; i < numberOfInputs_; ++i) {
    for(int j = i + 1; j < numberOfInputs_; ++j) {
      diagramsDistanceMatrix_[j][i] = diagramsDistanceMatrix_[i][j];
    }
  }
}

void PersistenceDiagramDistanceMatrix::updateClusters() {
  if(k_ > 1) {
    std::vector<std::vector<double>> distance_matrix = getDistanceMatrix();
    old_clustering_ = clustering_;
    invertClusters();
    initializeEmptyClusters();

    for(int i = 0; i < numberOfInputs_; ++i) {
      double min_distance_to_centroid = std::numeric_limits<double>::max();
      int cluster = -1;
      for(int c = 0; c < k_; ++c) {
        if(distance_matrix[i][c] < min_distance_to_centroid) {
          min_distance_to_centroid = distance_matrix[i][c];
          cluster = c;
        }
      }

      clustering_[cluster].push_back(i);
      if(cluster != inv_clustering_[i]) {
        // New centroid attributed to this diagram
        resetDosToOriginalValues();
        barycenter_inputs_reset_flag = true;
        if(do_min_) {
          centroids_with_price_min_[i]
            = centroidWithZeroPrices(centroids_min_[cluster]);
        }
        if(do_sad_) {
          centroids_with_price_saddle_[i]
            = centroidWithZeroPrices(centroids_saddle_[cluster]);
        }
        if(do_max_) {
          centroids_with_price_max_[i]
            = centroidWithZeroPrices(centroids_max_[cluster]);
        }
        inv_clustering_[i] = cluster;
      }
    }
  } else {
    old_clustering_ = clustering_;
    invertClusters();
    initializeEmptyClusters();

    for(int i = 0; i < numberOfInputs_; i++) {
      clustering_[0].push_back(i);
      if(n_iterations_ < 1) {
        if(do_min_) {
          centroids_with_price_min_[i]
            = centroidWithZeroPrices(centroids_min_[0]);
        }
        if(do_sad_) {
          centroids_with_price_saddle_[i]
            = centroidWithZeroPrices(centroids_saddle_[0]);
        }
        if(do_max_) {
          centroids_with_price_max_[i]
            = centroidWithZeroPrices(centroids_max_[0]);
        }
      }
      inv_clustering_[i] = 0;
    }
  }
  return;
}

void PersistenceDiagramDistanceMatrix::invertClusters() {
  /// Converts the clustering (vector of vector of diagram's id) into
  /// a vector of size numberOfInputs_ containg the cluster of each input
  /// diagram.

  // Initializes clusters with -1
  inv_clustering_ = std::vector<int>(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; ++i) {
    inv_clustering_[i] = -1;
  }

  // Fill in the clusters
  for(int c = 0; c < k_; ++c) {
    for(unsigned int j = 0; j < clustering_[c].size(); ++j) {
      int idx = clustering_[c][j];
      inv_clustering_[idx] = c;
    }
  }
}

void PersistenceDiagramDistanceMatrix::invertInverseClusters() {
  clustering_ = std::vector<std::vector<int>>(k_);
  for(int i = 0; i < numberOfInputs_; ++i) {
    clustering_[inv_clustering_[i]].push_back(i);
  }

  // Check if a cluster was left without diagram
  for(int c = 0; c < k_; ++c) {
    if(clustering_[c].size() == 0) {
      std::cout << "Problem in invertInverseClusters()... \nCluster " << c
                << " was left with no diagram attached to it... " << std::endl;
    }
  }
}

void PersistenceDiagramDistanceMatrix::acceleratedUpdateClusters() {
  // Step 1
  getCentroidDistanceMatrix();
  old_clustering_ = clustering_;
  // self.old_clusters = copy.copy(self.clusters)
  invertClusters();
  initializeEmptyClusters();
  bool do_min = original_dos[0];
  bool do_sad = original_dos[1];
  bool do_max = original_dos[2];

  for(int i = 0; i < numberOfInputs_; ++i) {
    // Step 3 find potential changes of clusters
    BidderDiagram<double> D1_min, D1_sad, D1_max;
    if(do_min) {
      D1_min = diagramWithZeroPrices(current_bidder_diagrams_min_[i]);
    }
    if(do_sad) {
      D1_sad = diagramWithZeroPrices(current_bidder_diagrams_saddle_[i]);
    }
    if(do_max) {
      D1_max = diagramWithZeroPrices(current_bidder_diagrams_max_[i]);
    }

    for(int c = 0; c < k_; ++c) {
      if(inv_clustering_[i] == -1) {
        // If not yet assigned, assign it first to a random cluster

        if(deterministic_) {
          inv_clustering_[i] = i % k_;
        } else {
          std::cout << " - ASSIGNED TO A RANDOM CLUSTER " << '\n';
          inv_clustering_[i] = rand() % (k_);
        }

        r_[i] = true;
        if(do_min) {
          centroids_with_price_min_[i]
            = centroidWithZeroPrices(centroids_min_[inv_clustering_[i]]);
        }
        if(do_sad) {
          centroids_with_price_saddle_[i]
            = centroidWithZeroPrices(centroids_saddle_[inv_clustering_[i]]);
        }
        if(do_max) {
          centroids_with_price_max_[i]
            = centroidWithZeroPrices(centroids_max_[inv_clustering_[i]]);
        }
      }

      if(c != inv_clustering_[i] && u_[i] > l_[i][c]
         && u_[i] > 0.5 * centroidsDistanceMatrix_[inv_clustering_[i]][c]) {
        // Step 3a, If necessary, recompute the distance to centroid
        if(r_[i]) {
          double distance = 0;
          GoodDiagram<double> centroid_min, centroid_sad, centroid_max;
          if(do_min) {
            centroid_min
              = centroidWithZeroPrices(centroids_min_[inv_clustering_[i]]);
            distance += computeDistance(D1_min, centroid_min, 0.01);
          }
          if(do_sad) {
            centroid_sad
              = centroidWithZeroPrices(centroids_saddle_[inv_clustering_[i]]);
            distance += computeDistance(D1_sad, centroid_sad, 0.01);
          }
          if(do_max) {
            centroid_max
              = centroidWithZeroPrices(centroids_max_[inv_clustering_[i]]);
            distance += computeDistance(D1_max, centroid_max, 0.01);
          }
          r_[i] = false;
          u_[i] = distance;
          l_[i][inv_clustering_[i]] = distance;
        }
        // Step 3b, check if still potential change of clusters
        if((n_iterations_ > 2 || n_iterations_ < 1)
           && (u_[i] > l_[i][c]
               || u_[i]
                    > 0.5 * centroidsDistanceMatrix_[inv_clustering_[i]][c])) {
          BidderDiagram<double> diagram_min, diagram_sad, diagram_max;
          GoodDiagram<double> centroid_min, centroid_sad, centroid_max;
          double distance = 0;

          if(do_min) {
            centroid_min = centroidWithZeroPrices(centroids_min_[c]);
            diagram_min
              = diagramWithZeroPrices(current_bidder_diagrams_min_[i]);
            distance += computeDistance(diagram_min, centroid_min, 0.01);
          }
          if(do_sad) {
            centroid_sad = centroidWithZeroPrices(centroids_saddle_[c]);
            diagram_sad
              = diagramWithZeroPrices(current_bidder_diagrams_saddle_[i]);
            distance += computeDistance(diagram_sad, centroid_sad, 0.01);
          }
          if(do_max) {
            centroid_max = centroidWithZeroPrices(centroids_max_[c]);
            diagram_max
              = diagramWithZeroPrices(current_bidder_diagrams_max_[i]);
            distance += computeDistance(diagram_max, centroid_max, 0.01);
          }
          l_[i][c] = distance;
          // TODO Prices are lost here... If distance<self.u[i], we should keep
          // the prices
          if(distance < u_[i]) {
            // Changing cluster
            resetDosToOriginalValues();
            barycenter_inputs_reset_flag = true;
            u_[i] = distance;
            inv_clustering_[i] = c;

            if(do_min) {
              centroids_with_price_min_[i]
                = centroidWithZeroPrices(centroids_min_[c]);
            }
            if(do_sad) {
              centroids_with_price_saddle_[i]
                = centroidWithZeroPrices(centroids_saddle_[c]);
            }
            if(do_max) {
              centroids_with_price_max_[i]
                = centroidWithZeroPrices(centroids_max_[c]);
            }
          }
        }
      }
    }
  }
  invertInverseClusters();
  for(int c = 0; c < k_; ++c) {
    if(clustering_[c].size() == 0) {
      std::cout << "Adding artificial centroid because a cluster was empty"
                << std::endl;
      bool idx_acceptable = false;
      int idx = 0;
      int increment = 0;

      // std::cout<< " u_ : [ ";
      // for(int i=0; i<u_.size(); i++){
      //         std::cout<<" "<<u_[i];
      //         }
      //         std::cout<<" ] "<<std::endl;
      std::vector<double> copy_of_u(u_.size());
      copy_of_u = u_;
      while(!idx_acceptable) {
        auto argMax = std::max_element(copy_of_u.begin(), copy_of_u.end());
        idx = std::distance(copy_of_u.begin(), argMax);
        if(inv_clustering_[idx] < k_ && inv_clustering_[idx] >= 0
           && clustering_[inv_clustering_[idx]].size() > 1) {
          idx_acceptable = true;
          int cluster_removal = inv_clustering_[idx];
          clustering_[cluster_removal].erase(
            std::remove(clustering_[cluster_removal].begin(),
                        clustering_[cluster_removal].end(), idx),
            clustering_[cluster_removal].end());
        } else {
          if(copy_of_u.size() > (size_t)idx) {
            copy_of_u.erase(argMax);
          } else {
            idx_acceptable = true;
            int cluster_max = 0;
            if(clustering_[cluster_max].size() > 0) {
              idx = clustering_[cluster_max][0];
            }
            for(int i_test = 1; i_test < k_; i_test++) {
              if(clustering_[i_test].size() > clustering_[cluster_max].size()) {
                cluster_max = i_test;
                idx = clustering_[cluster_max][0];
              }
            }
            int cluster_removal = inv_clustering_[idx];
            clustering_[cluster_removal].erase(
              std::remove(clustering_[cluster_removal].begin(),
                          clustering_[cluster_removal].end(), idx),
              clustering_[cluster_removal].end());
          }
        }
        increment += 1;
      }

      clustering_[c].push_back(idx);
      inv_clustering_[idx] = c;

      if(do_min) {
        centroids_min_[c]
          = diagramToCentroid(current_bidder_diagrams_min_[idx]);
        centroids_with_price_min_[idx]
          = centroidWithZeroPrices(centroids_min_[c]);
      }
      if(do_sad) {
        centroids_saddle_[c]
          = diagramToCentroid(current_bidder_diagrams_saddle_[idx]);
        centroids_with_price_saddle_[idx]
          = centroidWithZeroPrices(centroids_saddle_[c]);
      }
      if(do_max) {
        centroids_max_[c]
          = diagramToCentroid(current_bidder_diagrams_max_[idx]);
        centroids_with_price_max_[idx]
          = centroidWithZeroPrices(centroids_max_[c]);
      }
      resetDosToOriginalValues();
      barycenter_inputs_reset_flag = true;
    }
  }
  return;
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
      centroids_with_price_min_.push_back(GoodDiagram<double>());
      std::vector<int> ids(bidders.size());
      for(unsigned int j = 0; j < ids.size(); j++) {
        ids[j] = -1;
      }
      current_bidder_ids_min_.push_back(ids);
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
      centroids_with_price_saddle_.push_back(GoodDiagram<double>());
      std::vector<int> ids(bidders.size());
      for(unsigned int j = 0; j < ids.size(); j++) {
        ids[j] = -1;
      }
      current_bidder_ids_sad_.push_back(ids);
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
      centroids_with_price_max_.push_back(GoodDiagram<double>());
      std::vector<int> ids(bidders.size());
      for(unsigned int j = 0; j < ids.size(); j++) {
        ids[j] = -1;
      }
      current_bidder_ids_max_.push_back(ids);
    }
  }
  return;
}

std::vector<double>
  PersistenceDiagramDistanceMatrix::enrichCurrentBidderDiagrams(
    std::vector<double> previous_min_persistence,
    std::vector<double> min_persistence,
    std::vector<std::vector<double>> initial_diagonal_prices,
    std::vector<std::vector<double>> initial_off_diagonal_prices,
    std::vector<int> min_points_to_add,
    bool add_points_to_barycenter) {

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
  // cout<<"\n max points to add for first min bidder
  // "<<max_points_to_add_max<<endl;
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
  // cout<<"candidates to be added to first diag :
  // "<<candidates_to_be_added_min[0].size()<<endl;

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
    int compteur_for_adding_points = 0;
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
          // cout<<"\n   size before adding to "<<i<<"th bidder:
          // "<<current_bidder_diagrams_min_[i].size()<<endl;
          current_bidder_diagrams_min_[i].addBidder(b);
          current_bidder_ids_min_[i]
                                 [candidates_to_be_added_min[i][idx_min[i][j]]]
            = current_bidder_diagrams_min_[i].size() - 1;
          // cout<<"   size after adding
          // "<<current_bidder_diagrams_min_[i].size()<<endl;

          if(use_accelerated_ && n_iterations_ > 0) {
            for(int c = 0; c < k_; ++c) {
              // Step 5 of Accelerated KMeans: Update the lower bound on
              // distance thanks to the triangular inequality
              l_[i][c] = pow(
                pow(l_[i][c], 1. / wasserstein_) - persistence / pow(2, 0.5),
                wasserstein_);
              if(l_[i][c] < 0) {
                l_[i][c] = 0;
              }
            }
            // Step 6, update the upper bound on the distance to the centroid
            // thanks to the triangle inequality
            u_[i]
              = pow(pow(u_[i], 1. / wasserstein_) + persistence / pow(2, 0.5),
                    wasserstein_);
            r_[i] = true;
          }
          int to_be_added_to_barycenter
            = deterministic_ ? compteur_for_adding_points % numberOfInputs_
                             : rand() % numberOfInputs_;
          if(to_be_added_to_barycenter == 0 && add_points_to_barycenter) {
            // std::cout << "here we are adding points to the centroid_min" <<
            // std::endl;
            for(int k = 0; k < numberOfInputs_; k++) {
              if(inv_clustering_[i] == inv_clustering_[k]) {
                // std::cout<< "index
                // "<<centroids_with_price_min_[k].size()<<std::endl;
                Good<double> g = Good<double>(
                  b.x_, b.y_, false, centroids_with_price_min_[k].size());
                g.setPrice(initial_off_diagonal_prices[0][k]);
                g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
                centroids_with_price_min_[k].addGood(g);
                // std::cout<<"added for "<<k<<std::endl;
              }
            }
            // std::cout<<"size of centroid
            // "<<centroids_min_[inv_clustering_[i]].size()<<std::endl;
            Good<double> g = Good<double>(
              b.x_, b.y_, false, centroids_min_[inv_clustering_[i]].size());
            g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
            centroids_min_[inv_clustering_[i]].addGood(g);
            // std::cout << "all added" << std::endl;
          }
        }
        compteur_for_adding_points++;
      }
      if(debugLevel_ > 5)
        std::cout << " Diagram " << i
                  << " size : " << current_bidder_diagrams_min_[i].size()
                  << std::endl;
    }
    if(debugLevel_ > 3) {
      // cout<<" Added "<<compteur_for_adding_points<<" in min-sad
      // diagram"<<endl;
    }
  }
  if(do_sad_) {
    int compteur_for_adding_points = 0;
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
          current_bidder_ids_sad_[i]
                                 [candidates_to_be_added_sad[i][idx_sad[i][j]]]
            = current_bidder_diagrams_saddle_[i].size() - 1;

          if(use_accelerated_ && n_iterations_ > 0) {
            for(int c = 0; c < k_; ++c) {
              // Step 5 of Accelerated KMeans: Update the lower bound on
              // distance thanks to the triangular inequality
              l_[i][c] = pow(
                pow(l_[i][c], 1. / wasserstein_) - persistence / pow(2, 0.5),
                wasserstein_);
              if(l_[i][c] < 0) {
                l_[i][c] = 0;
              }
            }
            // Step 6, update the upper bound on the distance to the centroid
            // thanks to the triangle inequality
            u_[i]
              = pow(pow(u_[i], 1. / wasserstein_) + persistence / pow(2, 0.5),
                    wasserstein_);
            r_[i] = true;
          }
          int to_be_added_to_barycenter
            = deterministic_ ? compteur_for_adding_points % numberOfInputs_
                             : rand() % numberOfInputs_;
          if(to_be_added_to_barycenter == 0 && add_points_to_barycenter) {
            // std::cout << "here we are adding points to the centroid_min" <<
            // std::endl;
            for(int k = 0; k < numberOfInputs_; k++) {
              if(inv_clustering_[i] == inv_clustering_[k]) {
                Good<double> g = Good<double>(
                  b.x_, b.y_, false, centroids_with_price_saddle_[k].size());
                g.setPrice(initial_off_diagonal_prices[1][k]);
                g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
                centroids_with_price_saddle_[k].addGood(g);
                // std::cout<<"added for "<<k<<std::endl;
              }
            }
            // std::cout << "all added" << std::endl;
          }
        }
        compteur_for_adding_points++;
      }
      if(debugLevel_ > 5)
        std::cout << " Diagram " << i
                  << " size : " << current_bidder_diagrams_saddle_[i].size()
                  << std::endl;
    }
    if(debugLevel_ > 3) {
      cout << " Added " << compteur_for_adding_points << " in sad-sad diagram"
           << endl;
    }
  }
  if(do_max_) {
    int compteur_for_adding_points = 0;
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
          current_bidder_ids_max_[i]
                                 [candidates_to_be_added_max[i][idx_max[i][j]]]
            = current_bidder_diagrams_max_[i].size() - 1;

          if(use_accelerated_ && n_iterations_ > 0) {
            for(int c = 0; c < k_; ++c) {
              // Step 5 of Accelerated KMeans: Update the lower bound on
              // distance thanks to the triangular inequality
              l_[i][c] = pow(
                pow(l_[i][c], 1. / wasserstein_) - persistence / pow(2, 0.5),
                wasserstein_);
              if(l_[i][c] < 0) {
                l_[i][c] = 0;
              }
            }
            // Step 6, update the upper bound on the distance to the centroid
            // thanks to the triangle inequality
            u_[i]
              = pow(pow(u_[i], 1. / wasserstein_) + persistence / pow(2, 0.5),
                    wasserstein_);
            r_[i] = true;
          }
          int to_be_added_to_barycenter
            = deterministic_ ? compteur_for_adding_points % numberOfInputs_
                             : rand() % numberOfInputs_;
          if(to_be_added_to_barycenter == 0 && add_points_to_barycenter) {
            // std::cout << "here we are adding points to the centroid_max" <<
            // std::endl;
            for(int k = 0; k < numberOfInputs_; k++) {
              if(inv_clustering_[i] == inv_clustering_[k]) {
                Good<double> g = Good<double>(
                  b.x_, b.y_, false, centroids_with_price_max_[k].size());
                g.setPrice(initial_off_diagonal_prices[2][k]);
                g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
                centroids_with_price_max_[k].addGood(g);
                // std::cout<<"added for "<<k<<std::endl;
              }
            }
            Good<double> g = Good<double>(
              b.x_, b.y_, false, centroids_max_[inv_clustering_[i]].size());
            g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
            centroids_max_[inv_clustering_[i]].addGood(g);
            // std::cout << "all added" << std::endl;
          }
        }
        compteur_for_adding_points++;
      }
      if(debugLevel_ > 5)
        std::cout << " Diagram " << i
                  << " size : " << current_bidder_diagrams_max_[i].size()
                  << std::endl;
    }
    if(debugLevel_ > 3) {
      cout << " Added " << compteur_for_adding_points << " in sad-max diagram"
           << endl;
    }
  }

  return new_min_persistence;
}
