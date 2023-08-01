/// \ingroup base
/// \class ttk::PDBarycenter
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramClustering

#include <PersistenceDiagramBarycenter.h>

#include <cmath>
#include <cstdlib> /* srand, rand */
#include <numeric>

std::vector<std::vector<ttk::MatchingType>>
  ttk::PDBarycenter::execute(DiagramType &barycenter) {
  return executeAuctionBarycenter(barycenter);
}

void ttk::PDBarycenter::runMatching(
  double *total_cost,
  double epsilon,
  std::vector<int> &sizes,
  KDT &kdt,
  std::vector<KDT *> &correspondence_kdt_map,
  std::vector<double> *min_diag_price,
  std::vector<double> *min_price,
  std::vector<std::vector<MatchingType>> *all_matchings,
  bool use_kdt,
  bool actual_distance) {
  Timer const time_matchings;

  double local_cost = *total_cost;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic, 1) reduction(+:local_cost)
#endif
  for(int i = 0; i < numberOfInputs_; i++) {
    double const delta_lim = 0.01;
    PersistenceDiagramAuction auction(
      current_bidder_diagrams_[i], barycenter_goods_[i], wasserstein_,
      geometrical_factor_, lambda_, delta_lim, kdt, correspondence_kdt_map,
      epsilon, min_diag_price->at(i), use_kdt);
    int n_biddings = 0;
    auction.initLowerBoundCostWeight(delta_lim);
    auction.initLowerBoundCost(i);
    auction.buildUnassignedBidders();
    auction.reinitializeGoods();
    auction.runAuctionRound(n_biddings, i);
    auction.updateDiagonalPrices();
    min_diag_price->at(i) = auction.getMinimalDiagonalPrice();
    min_price->at(i) = getMinimalPrice(i);
    std::vector<MatchingType> matchings;
    double const cost = auction.getMatchingsAndDistance(matchings, true);
    all_matchings->at(i) = matchings;
    if(actual_distance) {
      local_cost += sqrt(cost);
    } else {
      local_cost += cost;
    }

    double const quotient
      = epsilon * auction.getAugmentedNumberOfBidders() / cost;
    precision_[i] = quotient < 1 ? 1. / sqrt(1 - quotient) - 1 : 10;
    if(auction.getRelativePrecision() == 0)
      precision_[i] = 0;
    // Resizes the diagram which was enrich with diagonal bidders
    // during the auction
    // TODO do this inside the auction !
    current_bidder_diagrams_[i].resize(sizes[i]);
  }
  *total_cost = local_cost;
}

void ttk::PDBarycenter::runMatchingAuction(
  double *total_cost,
  std::vector<int> &sizes,
  KDT &kdt,
  std::vector<KDT *> &correspondence_kdt_map,
  std::vector<double> *min_diag_price,
  std::vector<std::vector<MatchingType>> *all_matchings,
  bool use_kdt,
  bool actual_distance) {
  double local_cost = *total_cost;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic, 1) reduction(+:local_cost)
#endif
  for(int i = 0; i < numberOfInputs_; i++) {
    PersistenceDiagramAuction auction(
      current_bidder_diagrams_[i], barycenter_goods_[i], wasserstein_,
      geometrical_factor_, lambda_, 0.01, kdt, correspondence_kdt_map, 0,
      (*min_diag_price)[i], use_kdt);
    std::vector<MatchingType> matchings;
    double const cost = auction.run(matchings, i);
    all_matchings->at(i) = matchings;
    if(actual_distance) {
      local_cost += sqrt(cost);
    } else {
      local_cost += cost;
    }

    // Resizes the diagram which was enrich with diagonal bidders
    // during the auction
    // TODO do this inside the auction !
    current_bidder_diagrams_[i].resize(sizes[i]);
  }
  *total_cost = local_cost;
}

bool ttk::PDBarycenter::hasBarycenterConverged(
  std::vector<std::vector<MatchingType>> &matchings,
  std::vector<std::vector<MatchingType>> &previous_matchings) {
  if(points_added_ > 0 || points_deleted_ > 0 || previous_matchings.empty()) {
    return false;
  }

  for(size_t j = 0; j < matchings.size(); j++) {
    for(size_t i = 0; i < matchings[j].size(); i++) {
      MatchingType t = matchings[j][i];
      MatchingType previous_t = previous_matchings[j][i];

      if(std::get<1>(t) != std::get<1>(previous_t)
         && (std::get<0>(t) >= 0 && std::get<0>(previous_t) >= 0)) {
        return false;
      }
    }
  }
  return true;
}

std::vector<std::vector<ttk::MatchingType>> ttk::PDBarycenter::correctMatchings(
  std::vector<std::vector<MatchingType>> &previous_matchings) {

  std::vector<std::vector<MatchingType>> corrected_matchings(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; i++) {
    // 1. Invert the current_bidder_ids_ vector
    std::vector<int> new_to_old_id(current_bidder_diagrams_[i].size());
    for(size_t j = 0; j < current_bidder_ids_[i].size(); j++) {
      int const new_id = current_bidder_ids_[i][j];
      if(new_id >= 0) {
        new_to_old_id[new_id] = j;
      }
    }
    // 2. Reconstruct the matchings
    std::vector<MatchingType> matchings_diagram_i;
    for(size_t j = 0; j < previous_matchings[i].size(); j++) {
      MatchingType m = previous_matchings[i][j];
      int const new_id = std::get<0>(m);
      if(new_id >= 0 && std::get<1>(m) >= 0) {
        std::get<0>(m) = new_to_old_id[new_id];
        matchings_diagram_i.push_back(m);
      }
    }
    corrected_matchings[i] = matchings_diagram_i;
  }
  return corrected_matchings;
}

double ttk::PDBarycenter::updateBarycenter(
  std::vector<std::vector<MatchingType>> &matchings) {
  // 1. Initialize variables used in the sequel
  Timer const t_update;
  size_t const n_goods = barycenter_goods_[0].size();

  size_t const n_diagrams = current_bidder_diagrams_.size();
  points_added_ = 0;
  points_deleted_ = 0;
  double max_shift = 0;

  std::vector<size_t> count_diag_matchings(
    n_goods); // Number of diagonal matchings for each point of the barycenter
  std::vector<double> x(n_goods);
  std::vector<double> y(n_goods);
  std::vector<double> crit_coords_x(n_goods);
  std::vector<double> crit_coords_y(n_goods);
  std::vector<double> crit_coords_z(n_goods);
  for(size_t i = 0; i < n_goods; i++) {
    count_diag_matchings[i] = 0;
    x[i] = 0;
    y[i] = 0;
    crit_coords_x[i] = 0;
    crit_coords_y[i] = 0;
    crit_coords_z[i] = 0;
  }
  std::vector<double> min_prices(n_diagrams);
  for(size_t j = 0; j < n_diagrams; j++) {
    min_prices[j] = std::numeric_limits<double>::max();
  }

  std::vector<Bidder *>
    points_to_append; // Will collect bidders linked to diagonal
  // 2. Preprocess the matchings
  for(size_t j = 0; j < matchings.size(); j++) {
    for(size_t i = 0; i < matchings[j].size(); i++) {
      int const bidder_id = std::get<0>(matchings[j][i]);
      int const good_id = std::get<1>(matchings[j][i]);
      if(good_id < 0 && bidder_id >= 0) {
        // Future new barycenter point
        points_to_append.push_back(&current_bidder_diagrams_[j].at(bidder_id));
      }

      else if(good_id >= 0 && bidder_id >= 0) {
        // Update coordinates (to be divided by the number of diagrams later on)
        x[good_id] += current_bidder_diagrams_[j].at(bidder_id).x_;
        y[good_id] += current_bidder_diagrams_[j].at(bidder_id).y_;
        if(geometrical_factor_ < 1) {
          const auto &critical_coordinates = current_bidder_diagrams_[j]
                                               .at(bidder_id)
                                               .GetCriticalCoordinates();
          crit_coords_x[good_id] += critical_coordinates[0];
          crit_coords_y[good_id] += critical_coordinates[1];
          crit_coords_z[good_id] += critical_coordinates[2];
        }
      } else if(good_id >= 0 && bidder_id < 0) {
        // Counting the number of times this barycenter point is linked to the
        // diagonal
        count_diag_matchings[good_id] = count_diag_matchings[good_id] + 1;
      }
    }
  }

  // 3. Update the previous points of the barycenter
  for(size_t i = 0; i < n_goods; i++) {
    if(count_diag_matchings[i] < n_diagrams) {
      // Barycenter point i is matched at least to one off-diagonal bidder
      // 3.1 Compute the arithmetic mean of off-diagonal bidders linked to it
      double const x_bar
        = x[i] / (double)(n_diagrams - count_diag_matchings[i]);
      double const y_bar
        = y[i] / (double)(n_diagrams - count_diag_matchings[i]);
      // 3.2 Compute the new coordinates of the point (the more linked to the
      // diagonal it was, the closer to the diagonal it'll be)
      double const new_x
        = ((double)(n_diagrams - count_diag_matchings[i]) * x_bar
           + (double)count_diag_matchings[i] * (x_bar + y_bar) / 2.)
          / (double)n_diagrams;
      double const new_y
        = ((double)(n_diagrams - count_diag_matchings[i]) * y_bar
           + (double)count_diag_matchings[i] * (x_bar + y_bar) / 2.)
          / (double)n_diagrams;
      // TODO Weight by persistence
      double const new_crit_coord_x
        = crit_coords_x[i] / (double)(n_diagrams - count_diag_matchings[i]);
      double const new_crit_coord_y
        = crit_coords_y[i] / (double)(n_diagrams - count_diag_matchings[i]);
      double const new_crit_coord_z
        = crit_coords_z[i] / (double)(n_diagrams - count_diag_matchings[i]);

      // 3.3 Compute and store how much the point has shifted
      // TODO adjust shift with geometrical_factor_
      double const dx = barycenter_goods_[0].at(i).x_ - new_x;
      double const dy = barycenter_goods_[0].at(i).y_ - new_y;
      double const shift = Geometry::pow(std::abs(dx), wasserstein_)
                           + Geometry::pow(std::abs(dy), wasserstein_);
      if(shift > max_shift) {
        max_shift = shift;
      }
      // 3.4 Update the position of the point
      for(size_t j = 0; j < n_diagrams; j++) {
        barycenter_goods_[j].at(i).SetCoordinates(new_x, new_y);
        if(geometrical_factor_ < 1) {
          barycenter_goods_[j].at(i).SetCriticalCoordinates(
            new_crit_coord_x, new_crit_coord_y, new_crit_coord_z);
        }
        if(barycenter_goods_[j].at(i).getPrice() < min_prices[j]) {
          min_prices[j] = barycenter_goods_[j].at(i).getPrice();
        }
      }
      // TODO Reinitialize/play with prices here if you wish
    }
  }
  for(size_t j = 0; j < n_diagrams; j++) {
    if(min_prices[j] >= std::numeric_limits<double>::max() / 2.) {
      min_prices[j] = 0;
    }
  }

  // 4. Delete off-diagonal barycenter points not linked to any
  // off-diagonal bidder
  for(size_t i = 0; i < n_goods; i++) {
    if(count_diag_matchings[i] == n_diagrams) {
      points_deleted_ += 1;
      double const shift
        = 2
          * Geometry::pow(
            barycenter_goods_[0].at(i).getPersistence() / 2., wasserstein_);
      if(shift > max_shift) {
        max_shift = shift;
      }
      for(size_t j = 0; j < n_diagrams; j++) {
        barycenter_goods_[j].at(i).id_ = -1;
      }
    }
  }

  // 5. Append the new points to the barycenter
  for(size_t k = 0; k < points_to_append.size(); k++) {
    points_added_ += 1;
    Bidder *b = points_to_append[k];
    double const gx
      = (b->x_ + (n_diagrams - 1) * (b->x_ + b->y_) / 2.) / (n_diagrams);
    double const gy
      = (b->y_ + (n_diagrams - 1) * (b->x_ + b->y_) / 2.) / (n_diagrams);
    const auto &critical_coordinates = b->GetCriticalCoordinates();
    for(size_t j = 0; j < n_diagrams; j++) {
      Good g = Good(gx, gy, false, barycenter_goods_[j].size());
      g.setPrice(min_prices[j]);
      if(geometrical_factor_ < 1) {
        g.SetCriticalCoordinates(std::get<0>(critical_coordinates),
                                 std::get<1>(critical_coordinates),
                                 std::get<2>(critical_coordinates));
      }
      barycenter_goods_[j].emplace_back(g);
      double const shift
        = 2
          * Geometry::pow(
            barycenter_goods_[j].at(g.id_).getPersistence() / 2., wasserstein_);
      if(shift > max_shift) {
        max_shift = shift;
      }
    }
  }

  // 6. Finally, recreate barycenter_goods
  for(size_t j = 0; j < n_diagrams; j++) {
    int count = 0;
    GoodDiagram new_barycenter;
    for(size_t i = 0; i < barycenter_goods_[j].size(); i++) {
      Good g = barycenter_goods_[j].at(i);
      if(g.id_ != -1) {
        g.id_ = count;
        new_barycenter.emplace_back(g);
        count++;
      }
    }
    barycenter_goods_[j] = new_barycenter;
  }

  return max_shift;
}

double ttk::PDBarycenter::getEpsilon(double rho) {
  return rho * rho / 8.0;
}

double ttk::PDBarycenter::getRho(double epsilon) {
  return std::sqrt(8.0 * epsilon);
}

void ttk::PDBarycenter::setBidderDiagrams() {

  for(int i = 0; i < numberOfInputs_; i++) {
    DiagramType *CTDiagram = &((*inputDiagrams_)[i]);

    BidderDiagram bidders;
    for(size_t j = 0; j < CTDiagram->size(); j++) {
      // Add bidder to bidders
      Bidder b((*CTDiagram)[j], j, lambda_);

      b.setPositionInAuction(bidders.size());
      bidders.emplace_back(b);
      if(b.isDiagonal() || b.x_ == b.y_) {
        this->printWrn("Diagonal point in diagram !!!");
      }
    }
    bidder_diagrams_.push_back(bidders);
    current_bidder_diagrams_.emplace_back();
    std::vector<int> ids(bidders.size());
    for(size_t j = 0; j < ids.size(); j++) {
      ids[j] = -1;
    }
    current_bidder_ids_.push_back(ids);
  }
}

double ttk::PDBarycenter::enrichCurrentBidderDiagrams(
  double previous_min_persistence,
  double min_persistence,
  std::vector<double> &initial_diagonal_prices,
  std::vector<double> &initial_off_diagonal_prices,
  int min_points_to_add,
  bool add_points_to_barycenter) {

  double new_min_persistence = min_persistence;

  // 1. Get size of the largest current diagram, deduce the maximal number of
  // points to append
  size_t max_diagram_size{};
  for(int i = 0; i < numberOfInputs_; i++) {
    if(current_bidder_diagrams_[i].size() > max_diagram_size) {
      max_diagram_size = current_bidder_diagrams_[i].size();
    }
  }
  int const max_points_to_add = std::max(
    min_points_to_add, min_points_to_add + (int)(max_diagram_size / 10));

  // 2. Get which points can be added, deduce the new minimal persistence
  std::vector<std::vector<int>> candidates_to_be_added(numberOfInputs_);
  std::vector<std::vector<int>> idx(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; i++) {

    std::vector<double> persistences;
    for(size_t j = 0; j < bidder_diagrams_[i].size(); j++) {
      Bidder const b = bidder_diagrams_[i].at(j);
      double const persistence = b.getPersistence();
      if(persistence >= min_persistence
         && persistence < previous_min_persistence) {
        candidates_to_be_added[i].push_back(j);
        idx[i].push_back(idx[i].size());
        persistences.push_back(persistence);
      }
    }
    sort(idx[i].begin(), idx[i].end(), [&persistences](int &a, int &b) {
      return ((persistences[a] > persistences[b])
              || ((persistences[a] == persistences[b]) && (a > b)));
    });
    int const size = candidates_to_be_added[i].size();
    if(size >= max_points_to_add) {
      double const last_persistence_added
        = persistences[idx[i][max_points_to_add - 1]];
      if(last_persistence_added > new_min_persistence) {
        new_min_persistence = last_persistence_added;
      }
    }
  }

  // 3. Add the points to the current diagrams

  // only to give determinism
  int compteur_for_adding_points = 0;

  for(int i = 0; i < numberOfInputs_; i++) {
    int const size = candidates_to_be_added[i].size();
    for(int j = 0; j < std::min(max_points_to_add, size); j++) {
      Bidder b = bidder_diagrams_[i].at(candidates_to_be_added[i][idx[i][j]]);
      if(b.getPersistence() >= new_min_persistence) {
        b.id_ = current_bidder_diagrams_[i].size();
        b.setPositionInAuction(current_bidder_diagrams_[i].size());
        b.setDiagonalPrice(initial_diagonal_prices[i]);
        current_bidder_diagrams_[i].emplace_back(b);
        // b.id_ --> position of b in current_bidder_diagrams_[i]
        current_bidder_ids_[i][candidates_to_be_added[i][idx[i][j]]]
          = current_bidder_diagrams_[i].size() - 1;

        int const to_be_added_to_barycenter
          = deterministic_ ? compteur_for_adding_points % numberOfInputs_
                           : rand() % numberOfInputs_;
        // We add the bidder as a good with probability 1/n_diagrams
        if(to_be_added_to_barycenter == 0 && add_points_to_barycenter) {
          for(int k = 0; k < numberOfInputs_; k++) {
            Good g = Good(b.x_, b.y_, false, barycenter_goods_[k].size());
            g.setPrice(initial_off_diagonal_prices[k]);
            g.SetCriticalCoordinates(b.coords_[0], b.coords_[1], b.coords_[2]);
            barycenter_goods_[k].emplace_back(g);
          }
        }
      }
      compteur_for_adding_points++;
    }
  }
  return new_min_persistence;
}

double ttk::PDBarycenter::getMaxPersistence() {
  double max_persistence = 0;
  for(int i = 0; i < numberOfInputs_; i++) {
    BidderDiagram &D = bidder_diagrams_[i];
    for(size_t j = 0; j < D.size(); j++) {
      // Add bidder to bidders
      Bidder const &b = D.at(j);
      double const persistence = b.getPersistence();
      if(persistence > max_persistence) {
        max_persistence = persistence;
      }
    }
  }
  return max_persistence;
}

double ttk::PDBarycenter::getMinimalPrice(int i) {
  double min_price = std::numeric_limits<double>::max();

  GoodDiagram &D = barycenter_goods_[i];
  if(D.empty()) {
    return 0;
  }
  for(size_t j = 0; j < D.size(); j++) {
    Good const &b = D.at(j);
    double const price = b.getPrice();
    if(price < min_price) {
      min_price = price;
    }
  }
  if(min_price >= std::numeric_limits<double>::max() / 2.) {
    return 0;
  }
  return min_price;
}

double ttk::PDBarycenter::getLowestPersistence() {
  double lowest_persistence = std::numeric_limits<double>::max();
  for(int i = 0; i < numberOfInputs_; i++) {
    BidderDiagram &D = bidder_diagrams_[i];
    for(size_t j = 0; j < D.size(); j++) {
      // Add bidder to bidders
      Bidder const &b = D.at(j);
      double const persistence = b.getPersistence();
      if(persistence < lowest_persistence && persistence > 0) {
        lowest_persistence = persistence;
      }
    }
  }
  if(lowest_persistence >= std::numeric_limits<double>::max() / 2.) {
    return 0;
  }
  return lowest_persistence;
}

void ttk::PDBarycenter::setInitialBarycenter(double min_persistence) {
  int size = 0;
  int random_idx;
  DiagramType *CTDiagram;
  int iter = 0;
  while(size == 0) {
    random_idx
      = deterministic_ ? iter % numberOfInputs_ : rand() % numberOfInputs_;
    CTDiagram = &((*inputDiagrams_)[random_idx]);
    for(int i = 0; i < numberOfInputs_; i++) {
      GoodDiagram goods;
      int count = 0;
      for(size_t j = 0; j < CTDiagram->size(); j++) {
        // Add bidder to bidders
        Good const g = Good((*CTDiagram)[j], count, lambda_);
        if(g.getPersistence() >= min_persistence) {
          goods.emplace_back(g);
          count++;
        }
      }
      if(static_cast<int>(barycenter_goods_.size()) < (i + 1)) {
        barycenter_goods_.push_back(goods);
      } else {
        barycenter_goods_[i] = goods;
      }
    }
    size = barycenter_goods_[0].size();
    iter++;
  }
}

typename ttk::PDBarycenter::KDTreePair ttk::PDBarycenter::getKDTree() const {
  Timer tm;
  auto kdt = std::make_unique<KDT>(true, wasserstein_);

  const int dimension = geometrical_factor_ >= 1 ? 2 : 5;

  std::vector<double> coordinates;
  std::vector<std::vector<double>> weights;

  for(size_t i = 0; i < barycenter_goods_[0].size(); i++) {
    const Good &g = barycenter_goods_[0].at(i);
    coordinates.push_back(geometrical_factor_ * g.x_);
    coordinates.push_back(geometrical_factor_ * g.y_);
    if(geometrical_factor_ < 1) {
      coordinates.push_back((1 - geometrical_factor_) * g.coords_[0]);
      coordinates.push_back((1 - geometrical_factor_) * g.coords_[1]);
      coordinates.push_back((1 - geometrical_factor_) * g.coords_[2]);
    }
  }

  for(size_t idx = 0; idx < barycenter_goods_.size(); idx++) {
    std::vector<double> const empty_weights;
    weights.push_back(empty_weights);
    for(size_t i = 0; i < barycenter_goods_[idx].size(); i++) {
      const Good &g = barycenter_goods_[idx].at(i);
      weights[idx].push_back(g.getPrice());
    }
  }
  // Correspondence map : position in barycenter_goods_ --> KDT node

  auto correspondence_kdt_map
    = kdt->build(coordinates.data(), barycenter_goods_[0].size(), dimension,
                 weights, barycenter_goods_.size());
  this->printMsg(" Building KDTree", 1, tm.getElapsedTime(),
                 debug::LineMode::NEW, debug::Priority::VERBOSE);
  return std::make_pair(std::move(kdt), correspondence_kdt_map);
}

std::vector<std::vector<ttk::MatchingType>>
  ttk::PDBarycenter::executeAuctionBarycenter(DiagramType &barycenter) {

  std::vector<std::vector<MatchingType>> previous_matchings;
  double const min_persistence = 0;
  double min_cost = std::numeric_limits<double>::max();
  int last_min_cost_obtained = 0;

  this->setBidderDiagrams();
  this->setInitialBarycenter(
    min_persistence); // false for a determinist initialization

  double const max_persistence = getMaxPersistence();

  std::vector<double> min_diag_price(numberOfInputs_);
  std::vector<double> min_price(numberOfInputs_);
  for(int i = 0; i < numberOfInputs_; i++) {
    min_diag_price[i] = 0;
    min_price[i] = 0;
  }

  int const min_points_to_add = std::numeric_limits<int>::max();
  this->enrichCurrentBidderDiagrams(2 * max_persistence, min_persistence,
                                    min_diag_price, min_price,
                                    min_points_to_add, false);

  bool converged = false;
  bool finished = false;
  double total_cost;

  while(!finished) {
    Timer const tm;

    std::pair<std::unique_ptr<KDT>, std::vector<KDT *>> pair;
    bool use_kdt = false;
    // If the barycenter is empty, do not compute the kdt (or it will crash :/)
    // TODO Fix KDTree to handle empty inputs...
    if(!barycenter_goods_[0].empty()) {
      pair = this->getKDTree();
      use_kdt = true;
    }

    std::vector<std::vector<MatchingType>> all_matchings(numberOfInputs_);
    std::vector<int> sizes(numberOfInputs_);
    for(int i = 0; i < numberOfInputs_; i++) {
      sizes[i] = current_bidder_diagrams_[i].size();
    }

    total_cost = 0;

    barycenter.clear();
    for(size_t j = 0; j < barycenter_goods_[0].size(); j++) {
      Good const &g = barycenter_goods_[0].at(j);
      barycenter.emplace_back(PersistencePair{CriticalVertex{0, nt1_, g.x_, {}},
                                              CriticalVertex{0, nt2_, g.y_, {}},
                                              diagramType_, true});
    }

    bool const actual_distance = (numberOfInputs_ == 2);
    runMatchingAuction(&total_cost, sizes, *pair.first, pair.second,
                       &min_diag_price, &all_matchings, use_kdt,
                       actual_distance);

    this->printMsg("Barycenter cost : " + std::to_string(total_cost),
                   debug::Priority::DETAIL);

    if(converged) {
      finished = true;
    }

    if(!finished) {
      updateBarycenter(all_matchings);

      if(min_cost > total_cost) {
        min_cost = total_cost;
        last_min_cost_obtained = 0;
      } else {
        last_min_cost_obtained += 1;
      }

      converged = converged || last_min_cost_obtained > 1;
    }

    previous_matchings = std::move(all_matchings);
    // END OF TIMER

    for(size_t i = 0; i < barycenter_goods_.size(); ++i) {
      for(size_t j = 0; j < barycenter_goods_[i].size(); ++j) {
        barycenter_goods_[i].at(j).setPrice(0);
      }
    }
    for(size_t i = 0; i < current_bidder_diagrams_.size(); ++i) {
      for(size_t j = 0; j < current_bidder_diagrams_[i].size(); ++j) {
        current_bidder_diagrams_[i].at(j).setDiagonalPrice(0);
      }
    }
    for(int i = 0; i < numberOfInputs_; i++) {
      min_diag_price[i] = 0;
      min_price[i] = 0;
    }
  }
  barycenter.resize(0);
  for(size_t j = 0; j < barycenter_goods_[0].size(); j++) {
    Good const &g = barycenter_goods_[0].at(j);
    barycenter.emplace_back(PersistencePair{CriticalVertex{0, nt1_, g.x_, {}},
                                            CriticalVertex{0, nt2_, g.y_, {}},
                                            diagramType_, true});
  }

  cost_ = total_cost;
  std::vector<std::vector<MatchingType>> corrected_matchings
    = correctMatchings(previous_matchings);
  return corrected_matchings;
}

double ttk::PDBarycenter::computeRealCost() {
  double total_real_cost = 0;
  std::vector<MatchingType> fake_matchings;
  for(int i = 0; i < numberOfInputs_; i++) {
    PersistenceDiagramAuction auction(
      wasserstein_, geometrical_factor_, lambda_, 0.01, true);
    GoodDiagram const current_barycenter = barycenter_goods_[0];
    BidderDiagram const current_bidder_diagram = bidder_diagrams_[i];
    auction.BuildAuctionDiagrams(current_bidder_diagram, current_barycenter);
    double const cost = auction.run(fake_matchings);
    total_real_cost += cost * cost;
  }
  return sqrt(total_real_cost);
}

bool ttk::PDBarycenter::isPrecisionObjectiveMet(double precision_objective,
                                                int mode) {
  if(mode == 0) { // ABSOLUTE PRECISION
    for(int i_input = 0; i_input < numberOfInputs_; i_input++) {
      if(precision_[i_input] > precision_objective) {
        return false;
      }
    }
  } else if(mode == 1) { // AVERAGE PRECISION
    double const average_precision
      = std::accumulate(precision_.begin(), precision_.end(), 0.0)
        / numberOfInputs_;
    if(average_precision > precision_objective) {
      return false;
    }
  }
  return true;
}
