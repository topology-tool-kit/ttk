#include "PersistenceDiagramAuction.h"

void ttk::PersistenceDiagramAuction::runAuctionRound(int &n_biddings,
                                                     const int kdt_index) {
  double const max_price = getMaximalPrice();
  double epsilon = epsilon_;
  if(epsilon_ < 1e-6 * max_price) {
    // Risks of floating point limits reached...
    epsilon = 1e-6 * max_price;
  }
  while(unassignedBidders_.size() > 0) {
    n_biddings++;
    int const pos = unassignedBidders_.front();
    Bidder &b = this->bidders_[pos];
    unassignedBidders_.pop();

    GoodDiagram &all_goods = b.isDiagonal() ? diagonal_goods_ : goods_;
    Good &twin_good = b.id_ >= 0 ? diagonal_goods_[b.id_] : goods_[-b.id_ - 1];
    // double eps = epsilon_*(1+0.05*n_biddings/bidders_.size());
    int idx_reassigned;
    if(b.isDiagonal()) {
      if(use_kdt_) {
        idx_reassigned = b.runDiagonalKDTBidding(
          &all_goods, twin_good, wasserstein_, epsilon, geometricalFactor_,
          correspondence_kdt_map_, diagonal_queue_, kdt_index);
      } else {
        idx_reassigned
          = b.runDiagonalBidding(&all_goods, twin_good, wasserstein_, epsilon,
                                 geometricalFactor_, diagonal_queue_);
      }
    } else {
      if(use_kdt_) {
        // We can use the kd-tree to speed up the search
        idx_reassigned
          = b.runKDTBidding(&all_goods, twin_good, wasserstein_, epsilon,
                            geometricalFactor_, &kdt_, kdt_index);
      } else {
        idx_reassigned = b.runBidding(
          &all_goods, twin_good, wasserstein_, epsilon, geometricalFactor_);
      }
    }
    if(idx_reassigned >= 0) {
      Bidder &reassigned = bidders_[idx_reassigned];
      reassigned.resetProperty();
      unassignedBidders_.push(idx_reassigned);
    }
  }
}

double ttk::PersistenceDiagramAuction::getMaximalPrice() {
  double max_price = 0;
  for(size_t i = 0; i < goods_.size(); ++i) {
    Good const &g = goods_[i];
    double const price = g.getPrice();
    if(price > max_price) {
      max_price = price;
    }
  }

  for(size_t i = 0; i < diagonal_goods_.size(); ++i) {
    Good const &g = diagonal_goods_[i];
    double const price = g.getPrice();
    if(price > max_price) {
      max_price = price;
    }
  }

  return max_price;
}

double ttk::PersistenceDiagramAuction::getMatchingsAndDistance(
  std::vector<MatchingType> &matchings, bool get_diagonal_matches) {
  double wassersteinDistance = 0;
  for(size_t i = 0; i < bidders_.size(); i++) {
    Bidder const &b = bidders_[i];
    if(!b.isDiagonal()) {
      int const good_id = b.getProperty().id_;
      double cost;

      if(good_id > -1) {
        // good is not diagonal
        cost = b.cost(b.getProperty(), wasserstein_, geometricalFactor_);
        MatchingType const t = std::make_tuple(i, good_id, cost);
        matchings.emplace_back(t);
      } else {
        if(!b.getProperty().isDiagonal()) {
          this->printWrn("both the bidder and the good are diagonal points");
        }
        cost = 2.0 * Geometry::pow(std::abs((b.y_ - b.x_) / 2.0), wasserstein_);
        if(get_diagonal_matches) {
          MatchingType const t = std::make_tuple(i, good_id, cost);
          matchings.emplace_back(t);
        }
      }
      wassersteinDistance += cost;
    } else {
      // b is diagonal
      const Good &g = b.getProperty();
      double const cost
        = 2.0 * Geometry::pow(std::abs((g.y_ - g.x_) / 2.0), wasserstein_);
      if(get_diagonal_matches) {
        MatchingType const t = std::make_tuple(b.id_, g.id_, cost);
        matchings.emplace_back(t);
      }
      wassersteinDistance += cost;
    }
  }
  return wassersteinDistance;
}

double ttk::PersistenceDiagramAuction::initLowerBoundCost(const int kdt_index) {
  lowerBoundCost_ = 0;
  for(unsigned int i = 0; i < bidders_.size(); ++i) {
    if(bidders_[i].isDiagonal())
      continue;

    // Get closest good
    double bestCost = std::numeric_limits<double>::max();
    std::vector<KDT *> neighbours;
    std::vector<double> costs;
    if(use_kdt_) {
      std::array<double, 5> coordinates;
      bidders_[i].GetKDTCoordinates(geometricalFactor_, coordinates);
      kdt_.getKClosest(1, coordinates, neighbours, costs, kdt_index);
      int const bestIndex = neighbours[0]->id_;
      bestCost
        = bidders_[i].cost(goods_[bestIndex], wasserstein_, geometricalFactor_);
    } else {
      for(unsigned int j = 0; j < goods_.size(); ++j) {
        double const cost
          = bidders_[i].cost(goods_[j], wasserstein_, geometricalFactor_);
        if(cost < bestCost)
          bestCost = cost;
      }
    }

    // Compare with diagonal good
    Good g{bidders_[i].x_, bidders_[i].y_, true, -bidders_[i].id_ - 1};
    g.projectOnDiagonal();
    double const cost = bidders_[i].cost(g, wasserstein_, geometricalFactor_);
    if(cost < bestCost)
      bestCost = cost;

    // Update lower bound
    lowerBoundCost_ += bestCost;
  }
  return lowerBoundCost_;
}

double ttk::PersistenceDiagramAuction::run(std::vector<MatchingType> &matchings,
                                           const int kdt_index) {
  initLowerBoundCostWeight(delta_lim_);
  initLowerBoundCost(kdt_index);
  initializeEpsilon();
  int n_biddings = 0;
  double delta = 5;
  while(delta > delta_lim_) {
    epsilon_ /= 5;
    this->buildUnassignedBidders();
    this->reinitializeGoods();
    this->runAuctionRound(n_biddings, kdt_index);
    delta = this->getRelativePrecision();
  }
  double const wassersteinDistance
    = this->getMatchingsAndDistance(matchings, true);
  return wassersteinDistance;
}

double ttk::PersistenceDiagramAuctionActor::cost(
  const PersistenceDiagramAuctionActor &g,
  const int wasserstein,
  const double geometricalFactor) const {

  if(is_diagonal_ && g.isDiagonal()) {
    return 0;
  } else if(is_diagonal_) {
    return geometricalFactor
             * (2 * Geometry::pow(std::abs(g.y_ / 2 - g.x_ / 2), wasserstein))
           + (1 - geometricalFactor) * getPairGeometricalLength(wasserstein);
  } else if(g.isDiagonal()) {
    return geometricalFactor
             * (2 * Geometry::pow(std::abs(y_ / 2 - x_ / 2), wasserstein))
           + (1 - geometricalFactor) * g.getPairGeometricalLength(wasserstein);
  } else {
    return geometricalFactor
             * (Geometry::pow(std::abs(x_ - g.x_), wasserstein)
                + Geometry::pow(std::abs(y_ - g.y_), wasserstein))
           + (1 - geometricalFactor)
               * (Geometry::pow(
                    std::abs(coords_[0] - g.coords_[0]), wasserstein)
                  + Geometry::pow(
                    std::abs(coords_[1] - g.coords_[1]), wasserstein)
                  + Geometry::pow(
                    std::abs(coords_[2] - g.coords_[2]), wasserstein));
  }
}

int ttk::Bidder::runBidding(GoodDiagram *goods,
                            Good &twinGood,
                            int wasserstein,
                            double epsilon,
                            double geometricalFactor) {
  double best_val = std::numeric_limits<double>::lowest();
  double second_val = std::numeric_limits<double>::lowest();
  Good *best_good{};
  for(size_t i = 0; i < goods->size(); i++) {
    Good &g = (*goods)[i];
    double val = -this->cost(g, wasserstein, geometricalFactor);
    val -= g.getPrice();
    if(val > best_val) {
      second_val = best_val;
      best_val = val;
      best_good = &g;
    } else if(val > second_val) {
      second_val = val;
    }
  }
  // And now check for the corresponding twin bidder
  Good &g = twinGood;
  double val = -this->cost(g, wasserstein, geometricalFactor);
  val -= g.getPrice();
  if(val > best_val) {
    second_val = best_val;
    best_val = val;
    best_good = &g;
  } else if(val > second_val) {
    second_val = val;
  }

  if(second_val == std::numeric_limits<double>::lowest()) {
    // There is only one acceptable good for the bidder
    second_val = best_val;
  }
  if(best_good == nullptr) {
    return -1;
  }
  double const old_price = best_good->getPrice();
  double new_price = old_price + best_val - second_val + epsilon;
  if(new_price > std::numeric_limits<double>::max() / 2) {
    new_price = old_price + epsilon;
  }
  // Assign bidder to best_good
  this->setProperty(*best_good);
  this->setPricePaid(new_price);

  // Assign best_good to bidder and unassign the previous owner of best_good
  // if need be
  int const idx_reassigned = best_good->getOwner();
  best_good->assign(this->position_in_auction_, new_price);
  return idx_reassigned;
}

int ttk::Bidder::runDiagonalBidding(
  GoodDiagram *goods,
  Good &twinGood,
  int wasserstein,
  double epsilon,
  double geometricalFactor,
  std::priority_queue<std::pair<int, double>,
                      std::vector<std::pair<int, double>>,
                      Compare> &diagonal_queue) {

  // First, find the lowest and second lowest weights for diagonal goods
  // Take this time to update the weights in the priority queue of the goods
  // tested. It is not necessary to update weights for all diagonal bidders in
  // the pririty queue since weights can only increase and we are interested
  // only in the lowest ones
  bool updated_top_pair
    = false; // Boolean which equals true iff the top pair in the priority
             // queue is given the good price
  bool const non_empty_goods = goods->size() > 0;
  std::pair<int, double> best_pair;
  if(non_empty_goods) {
    while(!updated_top_pair) {
      std::pair<int, double> top_pair = diagonal_queue.top();
      diagonal_queue.pop();

      double const queue_weight = top_pair.second;
      const auto &good = (*goods)[top_pair.first];
      if(good.getPrice() > queue_weight) {
        // If the weight in the priority queue is not the good one, update
        std::get<1>(top_pair) = good.getPrice();
        diagonal_queue.push(top_pair);
      } else {
        updated_top_pair = true;
        best_pair = top_pair;
      }
    }
  }

  std::pair<int, double> second_pair;
  if(non_empty_goods && diagonal_queue.size() > 0) {
    bool updated_second_pair = false;
    while(!updated_second_pair) {
      second_pair = diagonal_queue.top();
      double const queue_weight = second_pair.second;
      const auto &good = (*goods)[second_pair.first];
      if(good.getPrice() != queue_weight) {
        // If the weight in the priority queue is not the good one, update it
        diagonal_queue.pop();
        std::get<1>(second_pair) = good.getPrice();
        diagonal_queue.push(second_pair);
      } else {
        updated_second_pair = true;
      }
    }
  }

  Good *best_good{};
  double best_val = 0;
  double second_val = 0;
  if(non_empty_goods) {
    best_val = -best_pair.second;
    second_val = diagonal_queue.size() > 0 ? -second_pair.second : best_val;
    best_good = &(*goods)[best_pair.first];
  }

  // And now check for the corresponding twin bidder
  bool is_twin = false;
  Good &g = twinGood;
  double val = -this->cost(g, wasserstein, geometricalFactor);
  val -= g.getPrice();
  if(non_empty_goods) {
    if(val > best_val) {
      second_val = best_val;
      best_val = val;
      best_good = &g;
      is_twin = true;
    } else if(val > second_val || diagonal_queue.empty()) {
      second_val = val;
    }
  } else {
    best_val = val;
    best_good = &g;
    is_twin = true;
    second_val = best_val;
  }

  if(second_val == std::numeric_limits<double>::lowest()) {
    // There is only one acceptable good for the bidder
    second_val = best_val;
  }
  double const old_price = best_good->getPrice();
  double new_price = old_price + best_val - second_val + epsilon;
  if(new_price > std::numeric_limits<double>::max() / 2) {
    new_price = old_price + epsilon;
  }

  // Assign bidder to best_good
  this->setProperty(*best_good);
  this->setPricePaid(new_price);

  // Assign best_good to bidder and unassign the previous owner of best_good
  // if need be
  int const idx_reassigned = best_good->getOwner();
  best_good->assign(this->position_in_auction_, new_price);
  if(!is_twin) {
    std::get<1>(best_pair) = new_price;
  }
  if(non_empty_goods) {
    diagonal_queue.push(best_pair);
  }
  return idx_reassigned;
}

int ttk::Bidder::runDiagonalKDTBidding(
  GoodDiagram *goods,
  Good &twinGood,
  int wasserstein,
  double epsilon,
  double geometricalFactor,
  std::vector<KDT *> &correspondence_kdt_map,
  std::priority_queue<std::pair<int, double>,
                      std::vector<std::pair<int, double>>,
                      Compare> &diagonal_queue,
  const int kdt_index) {
  /// Runs bidding of a diagonal bidder

  // First, find the lowest and second lowest weights for diagonal goods
  // Take this time to update the weights in the priority queue of the goods
  // tested. It is not necessary to update weights for all diagonal bidders in
  // the pririty queue since weights can only increase and we are interested
  // only in the lowest ones
  bool updated_top_pair
    = false; // Boolean which equals true iff the top pair in the priority
             // queue is given the good price
  bool const non_empty_goods = goods->size() > 0;
  std::pair<int, double> best_pair;
  if(non_empty_goods) {
    while(!updated_top_pair) {
      std::pair<int, double> top_pair = diagonal_queue.top();
      double const queue_weight = top_pair.second;
      const auto &good = (*goods)[top_pair.first];

      diagonal_queue.pop();
      if(good.getPrice() > queue_weight) {
        // If the weight in the priority queue is not the good one, update it
        std::get<1>(top_pair) = good.getPrice();
        diagonal_queue.push(top_pair);
      } else {
        updated_top_pair = true;
        best_pair = top_pair;
      }
    }
  }
  std::pair<int, double> second_pair;
  if(non_empty_goods && diagonal_queue.size() > 0) {
    bool updated_second_pair = false;
    while(!updated_second_pair) {
      second_pair = diagonal_queue.top();
      double const queue_weight = second_pair.second;
      const auto &good = (*goods)[second_pair.first];
      if(good.getPrice() > queue_weight) {
        // If the weight in the priority queue is not the good one, update it
        diagonal_queue.pop();
        std::get<1>(second_pair) = good.getPrice();
        diagonal_queue.push(second_pair);
      } else {
        updated_second_pair = true;
      }
    }
  }

  Good *best_good{};
  double best_val = 0;
  double second_val = 0;
  if(non_empty_goods) {
    best_val = -best_pair.second;
    second_val = diagonal_queue.size() > 0 ? -second_pair.second : best_val;
    best_good = &(*goods)[best_pair.first];
  }

  // And now check for the corresponding twin bidder
  bool is_twin = false;
  Good &g = twinGood;
  double val = -this->cost(g, wasserstein, geometricalFactor);
  val -= g.getPrice();
  if(non_empty_goods) {
    if(val > best_val) {
      second_val = best_val;
      best_val = val;
      best_good = &g;
      is_twin = true;
    } else if(val > second_val || diagonal_queue.empty()) {
      second_val = val;
    }
  } else {
    best_val = val;
    best_good = &g;
    is_twin = true;
    second_val = best_val;
  }

  if(second_val == std::numeric_limits<double>::lowest()) {
    // There is only one acceptable good for the bidder
    second_val = best_val;
  }
  double const old_price = best_good->getPrice();
  double new_price = old_price + best_val - second_val + epsilon;
  if(new_price > std::numeric_limits<double>::max() / 2) {
    new_price = old_price + epsilon;
  }

  // Assign bidder to best_good
  this->setProperty(*best_good);
  this->setPricePaid(new_price);

  // Assign best_good to bidder and unassign the previous owner of best_good
  // if need be
  int const idx_reassigned = best_good->getOwner();
  best_good->assign(this->position_in_auction_, new_price);
  if(is_twin) {
    // Update weight in KDTree if the closest good is in it
    correspondence_kdt_map[best_good->id_]->updateWeight(new_price, kdt_index);
    if(non_empty_goods) {
      diagonal_queue.push(best_pair);
    }
  } else {
    if(non_empty_goods) {
      std::get<1>(best_pair) = new_price;
      diagonal_queue.push(best_pair);
    }
  }
  return idx_reassigned;
}

int ttk::Bidder::runKDTBidding(GoodDiagram *goods,
                               Good &twinGood,
                               int wasserstein,
                               double epsilon,
                               double geometricalFactor,
                               KDT *kdt,
                               const int kdt_index) {

  /// Runs bidding of a non-diagonal bidder
  std::vector<KDT *> neighbours;
  std::vector<double> costs;

  std::array<double, 5> coordinates;
  GetKDTCoordinates(geometricalFactor, coordinates);

  kdt->getKClosest(2, coordinates, neighbours, costs, kdt_index);
  double best_val, second_val;
  KDT *closest_kdt;
  Good *best_good{};
  if(costs.size() == 2) {
    std::array<int, 2> idx{0, 1};
    std::sort(idx.begin(), idx.end(),
              [&costs](int &a, int &b) { return costs[a] < costs[b]; });

    closest_kdt = neighbours[idx[0]];
    best_good = &(*goods)[closest_kdt->id_];
    // Value is defined as the opposite of cost (each bidder aims at
    // maximizing it)
    best_val = -costs[idx[0]];
    second_val = -costs[idx[1]];
  } else {
    // If the kdtree contains only one point
    closest_kdt = neighbours[0];
    best_good = &(*goods)[closest_kdt->id_];
    best_val = -costs[0];
    second_val = best_val;
  }
  // And now check for the corresponding twin bidder
  bool twin_chosen = false;
  Good &g = twinGood;
  double val = -this->cost(g, wasserstein, geometricalFactor);
  val -= g.getPrice();
  if(val > best_val) {
    second_val = best_val;
    best_val = val;
    best_good = &g;
    twin_chosen = true;
  } else if(val > second_val) {
    second_val = val;
  }

  if(second_val == std::numeric_limits<double>::lowest()) {
    // There is only one acceptable good for the bidder
    second_val = best_val;
  }
  double const old_price = best_good->getPrice();
  double new_price = old_price + best_val - second_val + epsilon;
  if(new_price > std::numeric_limits<double>::max() / 2) {
    new_price = old_price + epsilon;
  }
  // Assign bidder to best_good
  this->setProperty(*best_good);
  this->setPricePaid(new_price);
  // Assign best_good to bidder and unassign the previous owner of best_good
  // if need be
  int const idx_reassigned = best_good->getOwner();

  best_good->assign(this->position_in_auction_, new_price);
  // Update the price in the KDTree
  if(!twin_chosen) {
    closest_kdt->updateWeight(new_price, kdt_index);
  }
  return idx_reassigned;
}
