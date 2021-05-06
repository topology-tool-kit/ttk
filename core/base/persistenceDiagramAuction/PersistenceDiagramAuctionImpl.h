#ifndef _PERSISTENCEDIAGRAMAUCTIONIMPL_H
#define _PERSISTENCEDIAGRAMAUCTIONIMPL_H

#ifndef matchingTuple
#define matchingTuple std::tuple<SimplexId, SimplexId, dataType>
#endif

template <typename dataType>
void ttk::PersistenceDiagramAuction<dataType>::runAuctionRound(
  int &n_biddings, const int kdt_index) {
  dataType max_price = getMaximalPrice();
  dataType epsilon = epsilon_;
  if(epsilon_ < 1e-6 * max_price) {
    // Risks of floating point limits reached...
    epsilon = 1e-6 * max_price;
  }
  while(unassignedBidders_.size() > 0) {
    n_biddings++;
    int pos = unassignedBidders_.front();
    Bidder<dataType> &b = bidders_.get(pos);
    unassignedBidders_.pop();

    GoodDiagram<dataType> &all_goods
      = b.isDiagonal() ? diagonal_goods_ : goods_;
    Good<dataType> &twin_good
      = b.id_ >= 0 ? diagonal_goods_.get(b.id_) : goods_.get(-b.id_ - 1);
    // dataType eps = epsilon_*(1+0.05*n_biddings/bidders_.size());
    int idx_reassigned;
    if(b.isDiagonal()) {
      if(use_kdt_) {
        idx_reassigned = b.runDiagonalKDTBidding(
          &all_goods, twin_good, wasserstein_, epsilon, geometricalFactor_,
          correspondance_kdt_map_, diagonal_queue_, kdt_index);
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
    /*if(n_biddings>-1){
      std::cout.precision(10);
      std::cout<< n_biddings << ", out : " <<idx_reassigned<< ", " << b.id_ <<
    "->" << b.getProperty()->id_ << ", distance : "<< b.cost(b.getProperty(),
    wasserstein_, geometricalFactor_) << ", price : "
    <<b.getProperty()->getPrice() << ", diagonal ? " << b.isDiagonal() << ", "<<
    b.getProperty()->isDiagonal() << std::endl;
    }*/
    if(idx_reassigned >= 0) {
      Bidder<dataType> &reassigned = bidders_.get(idx_reassigned);
      reassigned.resetProperty();
      unassignedBidders_.push(idx_reassigned);
    }
  }
}

template <typename dataType>
dataType ttk::PersistenceDiagramAuction<dataType>::getMaximalPrice() {
  dataType max_price = 0;
  for(int i = 0; i < goods_.size(); ++i) {
    Good<dataType> &g = goods_.get(i);
    dataType price = g.getPrice();
    if(price > max_price) {
      max_price = price;
    }
  }

  for(int i = 0; i < diagonal_goods_.size(); ++i) {
    Good<dataType> &g = diagonal_goods_.get(i);
    dataType price = g.getPrice();
    if(price > max_price) {
      max_price = price;
    }
  }

  return max_price;
}

template <typename dataType>
dataType ttk::PersistenceDiagramAuction<dataType>::getMatchingsAndDistance(
  std::vector<matchingTuple> *matchings, bool get_diagonal_matches) {
  dataType wassersteinDistance = 0;
  for(int i = 0; i < bidders_.size(); i++) {
    Bidder<dataType> &b = bidders_.get(i);
    if(!b.isDiagonal()) {
      int good_id = b.getProperty().id_;
      dataType cost;

      if(good_id > -1) {
        // good is not diagonal
        cost = b.cost(b.getProperty(), wasserstein_, geometricalFactor_);
        matchingTuple t = std::make_tuple(i, good_id, cost);
        matchings->push_back(t);
      } else {
        if(!b.getProperty().isDiagonal()) {
          std::cout << "[PersistenceDiagramAuctionImpl] Huho : both the bidder "
                       "and the good are "
                       "diagonal points"
                    << std::endl;
        }
        cost
          = 2 * Geometry::pow(abs<dataType>((b.y_ - b.x_) / 2), wasserstein_);
        if(get_diagonal_matches) {
          matchingTuple t = std::make_tuple(i, good_id, cost);
          matchings->push_back(t);
        }
      }
      wassersteinDistance += cost;
    } else {
      // b is diagonal
      const Good<dataType> &g = b.getProperty();
      dataType cost
        = 2 * Geometry::pow(abs<dataType>((g.y_ - g.x_) / 2), wasserstein_);
      if(get_diagonal_matches) {
        matchingTuple t = std::make_tuple(b.id_, g.id_, cost);
        matchings->push_back(t);
      }
      wassersteinDistance += cost;
    }
  }
  return wassersteinDistance;
}

template <typename dataType>
dataType ttk::PersistenceDiagramAuction<dataType>::run(
  std::vector<matchingTuple> *matchings) {
  initializeEpsilon();
  int n_biddings = 0;
  dataType delta = 5;
  while(delta > delta_lim_) {
    epsilon_ /= 5;
    this->buildUnassignedBidders();
    this->reinitializeGoods();
    this->runAuctionRound(n_biddings);
    delta = this->getRelativePrecision();
  }
  dataType wassersteinDistance = this->getMatchingsAndDistance(matchings, true);
  return wassersteinDistance;
}
#endif
