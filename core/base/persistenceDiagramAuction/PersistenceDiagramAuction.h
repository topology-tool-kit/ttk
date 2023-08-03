#pragma once

#include <PersistenceDiagramAuctionActor.h>

#include <limits>

namespace ttk {
  class PersistenceDiagramAuction : public Debug {

  public:
    inline int getAugmentedNumberOfBidders() {
      return bidders_.size();
    }

    using KDT = Bidder::KDT;

    KDT default_kdt_{};
    KDT &kdt_{default_kdt_};
    std::vector<KDT *> default_correspondence_kdt_map_{};
    std::vector<KDT *> &correspondence_kdt_map_{
      default_correspondence_kdt_map_};

    inline void initLowerBoundCostWeight(double delta_lim) {
      lowerBoundCostWeight_ = 1 + delta_lim;
    }

    double initLowerBoundCost(const int kdt_index = 0);

    PersistenceDiagramAuction(int wasserstein,
                              double geometricalFactor,
                              double lambda,
                              double delta_lim,
                              bool use_kdTree)
      : wasserstein_{wasserstein}, geometricalFactor_{geometricalFactor},
        lambda_{lambda}, delta_lim_{delta_lim}, use_kdt_{use_kdTree} {
    }

    PersistenceDiagramAuction(BidderDiagram &bidders,
                              GoodDiagram &goods,
                              int wasserstein,
                              double geometricalFactor,
                              double lambda,
                              double delta_lim,
                              KDT &kdt,
                              std::vector<KDT *> &correspondence_kdt_map,
                              double epsilon = {},
                              double initial_diag_price = {},
                              bool use_kdTree = true)
      : kdt_{kdt}, correspondence_kdt_map_{correspondence_kdt_map},
        bidders_{bidders}, goods_{goods} {

      n_bidders_ = bidders.size();
      n_goods_ = goods.size();

      for(int i = 0; i < n_bidders_; i++) {
        // Add diagonal goods
        Bidder const &b = bidders_[i];
        Good g{b.x_, b.y_, true, -b.id_ - 1};
        g.projectOnDiagonal();
        if(b.diagonal_price_ > 0) {
          g.setPrice(b.diagonal_price_);
        } else {
          g.setPrice(initial_diag_price);
        }
        diagonal_goods_.emplace_back(g);
        std::pair<int, double> const pair = std::make_pair(i, g.getPrice());
        diagonal_queue_.push(pair);
      }
      for(int i = 0; i < n_goods_; i++) {
        // Add diagonal bidders
        Good const &g = goods_[i];
        Bidder b{g.x_, g.y_, true, -g.id_ - 1};
        b.projectOnDiagonal();
        b.setPositionInAuction(bidders_.size());
        bidders_.emplace_back(b);
      }

      epsilon_ = epsilon;
      wasserstein_ = wasserstein;
      delta_lim_ = delta_lim;
      geometricalFactor_ = geometricalFactor;
      lambda_ = lambda;

      use_kdt_ = (use_kdTree && !goods_.empty());
    }

    void runAuctionRound(int &n_biddings, const int kdt_index = 0);
    double getMatchingsAndDistance(std::vector<MatchingType> &matchings,
                                   bool get_diagonal_matches = false);
    double run(std::vector<MatchingType> &matchings, const int kdt_index = 0);
    double run() {
      std::vector<MatchingType> matchings{};
      return this->run(matchings);
    }
    double getMaximalPrice();

    void BuildAuctionDiagrams(const BidderDiagram &BD, const GoodDiagram &GD) {
      this->n_bidders_ = BD.size();
      this->n_goods_ = GD.size();
      // delete_kdTree_ = false;
      this->bidders_ = BD;
      this->goods_ = GD;

      for(int i = 0; i < n_bidders_; i++) {
        // Add diagonal goods
        Bidder const &b = bidders_[i];
        Good g{b.x_, b.y_, true, -b.id_ - 1};
        g.projectOnDiagonal();
        diagonal_goods_.emplace_back(g);
        std::pair<int, double> const pair = std::make_pair(i, g.getPrice());
        diagonal_queue_.push(pair);
      }
      for(int i = 0; i < n_goods_; i++) {
        // Add diagonal bidders
        Good const &g = goods_[i];
        Bidder b{g.x_, g.y_, true, -g.id_ - 1};
        b.projectOnDiagonal();
        b.setPositionInAuction(bidders_.size());
        bidders_.emplace_back(b);
      }
      if(!goods_.empty()) {
        // use_kdt_ = use_kdt_;
        this->buildKDTree();
      } else {
        use_kdt_ = false;
      }
    }

    void BuildAuctionDiagrams(const DiagramType &diagram1,
                              const DiagramType &diagram2) {
      n_bidders_ = diagram1.size();
      n_goods_ = diagram2.size();
      this->setBidders(diagram1);
      this->setGoods(diagram2);
      for(int i = 0; i < n_bidders_; i++) {
        // Add diagonal goods
        Bidder const &b = bidders_[i];
        Good g{b.x_, b.y_, true, -b.id_ - 1};
        g.projectOnDiagonal();
        diagonal_goods_.emplace_back(g);
        std::pair<int, double> const pair = std::make_pair(i, g.getPrice());
        diagonal_queue_.push(pair);
      }
      for(int i = 0; i < n_goods_; i++) {
        // Add diagonal bidders
        Good const &g = goods_[i];
        Bidder b{g.x_, g.y_, true, -g.id_ - 1};
        b.projectOnDiagonal();
        b.setPositionInAuction(bidders_.size());
        bidders_.emplace_back(b);
      }
      if(!bidders_.empty()) {
        use_kdt_ = true;
        this->buildKDTree();
      } else {
        use_kdt_ = false;
      }
    }

    void setBidders(const DiagramType &diagram1) {
      for(size_t i = 0; i < diagram1.size(); i++) {
        // Add bidder to bidders
        Bidder b{diagram1[i], static_cast<int>(i), lambda_};
        b.setPositionInAuction(bidders_.size());
        bidders_.emplace_back(b);
      }
      n_bidders_ = bidders_.size();
    }

    void setGoods(const DiagramType &diagram2) {
      for(size_t i = 0; i < diagram2.size(); i++) {
        // Add bidder to bidders
        Good const g{diagram2[i], static_cast<int>(i), lambda_};
        goods_.emplace_back(g);
      }
      n_goods_ = goods_.size();
    }

    void buildKDTree() {
      Timer const t;
      default_kdt_ = KDT{true, wasserstein_};
      const int dimension
        = geometricalFactor_ >= 1 ? (geometricalFactor_ <= 0 ? 3 : 2) : 5;
      std::vector<double> coordinates;
      for(size_t i = 0; i < goods_.size(); i++) {
        Good &g = goods_[i];
        if(geometricalFactor_ > 0) {
          coordinates.push_back(geometricalFactor_ * g.x_);
          coordinates.push_back(geometricalFactor_ * g.y_);
        }
        if(geometricalFactor_ < 1) {
          coordinates.push_back((1 - geometricalFactor_) * g.coords_[0]);
          coordinates.push_back((1 - geometricalFactor_) * g.coords_[1]);
          coordinates.push_back((1 - geometricalFactor_) * g.coords_[2]);
        }
      }
      correspondence_kdt_map_
        = kdt_.build(coordinates.data(), goods_.size(), dimension);
    }

    void setEpsilon(const double epsilon) {
      epsilon_ = epsilon;
    }

    void initializeEpsilon() {
      double max_persistence = 0;
      for(const auto &b : this->bidders_) {
        double const persistence = b.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }
      for(const auto &g : this->goods_) {
        double const persistence = g.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }
      this->epsilon_ = 5.0 / 4 * Geometry::pow(max_persistence, wasserstein_);
    }

    void buildUnassignedBidders() {
      for(size_t i = 0; i < bidders_.size(); i++) {
        Bidder &b = bidders_[i];
        b.resetProperty();
        unassignedBidders_.push(i);
      }
    }

    void reinitializeGoods() {
      for(auto &g : this->goods_) {
        g.setOwner(-1);
      }
      for(auto &g : this->diagonal_goods_) {
        g.setOwner(-1);
      }
    }

    double getMatchingDistance() {
      double d = 0;
      for(size_t i = 0; i < bidders_.size(); i++) {
        Bidder const &b = bidders_[i];
        d += b.cost(b.getProperty(), wasserstein_, geometricalFactor_);
      }
      return d;
    }

    double getRelativePrecision() {
      double const d = this->getMatchingDistance();
      if(d < 1e-6 or d <= (lowerBoundCost_ * lowerBoundCostWeight_)) {
        return 0;
      }
      double const denominator = d - bidders_.size() * epsilon_;
      if(denominator <= 0) {
        return 1;
      } else {
        return Geometry::pow(d / denominator, 1 / ((float)wasserstein_)) - 1;
      }
    }

    void updateDiagonalPrices() {
      for(size_t i = 0; i < diagonal_goods_.size(); i++) {
        Bidder &b = bidders_[i];
        b.setDiagonalPrice(diagonal_goods_[i].getPrice());
      }
    }

    double getMinimalDiagonalPrice() {
      if(diagonal_goods_.empty()) {
        return 0;
      }
      double min_price = std::numeric_limits<double>::max();
      for(size_t i = 0; i < diagonal_goods_.size(); i++) {
        double const price = diagonal_goods_[i].getPrice();
        if(price < min_price) {
          min_price = price;
        }
      }
      return min_price;
    }

  protected:
    int wasserstein_{2}; // Power in Wassertsein distance (by default set to 2)
    BidderDiagram default_bidders_{};
    BidderDiagram &bidders_{default_bidders_};
    GoodDiagram default_goods_{};
    GoodDiagram &goods_{default_goods_};
    GoodDiagram diagonal_goods_;
    std::priority_queue<std::pair<int, double>,
                        std::vector<std::pair<int, double>>,
                        Compare>
      diagonal_queue_{};
    std::queue<int> unassignedBidders_{};

    int n_bidders_{0};
    int n_goods_{0};

    double epsilon_{1};
    double geometricalFactor_{};
    double lambda_{};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double delta_lim_{};
    double lowerBoundCost_, lowerBoundCostWeight_;
    bool use_kdt_{true};

  }; // namespace ttk
} // namespace ttk
