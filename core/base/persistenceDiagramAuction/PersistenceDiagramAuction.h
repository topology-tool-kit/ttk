#pragma once

#include <Debug.h>
#include <KDTree.h>
#include <PersistenceDiagramAuctionActor.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>

namespace ttk {
  class PersistenceDiagramAuction : public Debug {

  public:
    inline int getAugmentedNumberOfBidders() {
      return bidders_.size();
    }

    KDTree<double> default_kdt_{};
    KDTree<double> &kdt_{default_kdt_};
    std::vector<KDTree<double> *> default_correspondance_kdt_map_{};
    std::vector<KDTree<double> *> &correspondance_kdt_map_{
      default_correspondance_kdt_map_};

    PersistenceDiagramAuction(int wasserstein,
                              double geometricalFactor,
                              double lambda,
                              double delta_lim,
                              bool use_kdTree)
      : wasserstein_{wasserstein}, geometricalFactor_{geometricalFactor},
        lambda_{lambda}, delta_lim_{delta_lim}, use_kdt_{use_kdTree} {
    }

    PersistenceDiagramAuction(
      BidderDiagram &bidders,
      GoodDiagram &goods,
      int wasserstein,
      double geometricalFactor,
      double lambda,
      double delta_lim,
      KDTree<double> &kdt,
      std::vector<KDTree<double> *> &correspondance_kdt_map,
      double epsilon = {},
      double initial_diag_price = {},
      bool use_kdTree = true)
      : kdt_{kdt}, correspondance_kdt_map_{correspondance_kdt_map},
        bidders_{bidders}, goods_{goods} {

      n_bidders_ = bidders.size();
      n_goods_ = goods.size();

      for(int i = 0; i < n_bidders_; i++) {
        // Add diagonal goods
        Bidder &b = bidders_[i];
        Good g{b.x_, b.y_, true, -b.id_ - 1};
        g.projectOnDiagonal();
        if(b.diagonal_price_ > 0) {
          g.setPrice(b.diagonal_price_);
        } else {
          g.setPrice(initial_diag_price);
        }
        diagonal_goods_.emplace_back(g);
        std::pair<int, double> pair = std::make_pair(i, g.getPrice());
        diagonal_queue_.push(pair);
      }
      for(int i = 0; i < n_goods_; i++) {
        // Add diagonal bidders
        Good &g = goods_[i];
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

      use_kdt_ = (use_kdTree && goods_.size() > 0);
    }

    void runAuctionRound(int &n_biddings, const int kdt_index = 0);
    double getMatchingsAndDistance(std::vector<MatchingType> *matchings,
                                   bool get_diagonal_matches = false);
    double run(std::vector<MatchingType> *matchings);
    double run() {
      std::vector<MatchingType> matchings{};
      return this->run(&matchings);
    }
    double getMaximalPrice();

    void BuildAuctionDiagrams(const BidderDiagram *BD, const GoodDiagram *GD) {
      n_bidders_ = BD->size();
      n_goods_ = GD->size();
      // delete_kdTree_ = false;
      bidders_ = *BD;
      goods_ = *GD;

      for(int i = 0; i < n_bidders_; i++) {
        // Add diagonal goods
        Bidder &b = bidders_[i];
        Good g{b.x_, b.y_, true, -b.id_ - 1};
        g.projectOnDiagonal();
        diagonal_goods_.emplace_back(g);
        std::pair<int, double> pair = std::make_pair(i, g.getPrice());
        diagonal_queue_.push(pair);
      }
      for(int i = 0; i < n_goods_; i++) {
        // Add diagonal bidders
        Good &g = goods_[i];
        Bidder b{g.x_, g.y_, true, -g.id_ - 1};
        b.projectOnDiagonal();
        b.setPositionInAuction(bidders_.size());
        bidders_.emplace_back(b);
      }
      if(goods_.size() > 0) {
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
        Bidder &b = bidders_[i];
        Good g{b.x_, b.y_, true, -b.id_ - 1};
        g.projectOnDiagonal();
        diagonal_goods_.emplace_back(g);
        std::pair<int, double> pair = std::make_pair(i, g.getPrice());
        diagonal_queue_.push(pair);
      }
      for(int i = 0; i < n_goods_; i++) {
        // Add diagonal bidders
        Good &g = goods_[i];
        Bidder b{g.x_, g.y_, true, -g.id_ - 1};
        b.projectOnDiagonal();
        b.setPositionInAuction(bidders_.size());
        bidders_.emplace_back(b);
      }
      if(bidders_.size() > 0) {
        use_kdt_ = true;
        this->buildKDTree();
      } else {
        use_kdt_ = false;
      }
    }

    void setBidders(const DiagramType &diagram1) {
      int d1Size = (int)diagram1.size();

      for(int i = 0; i < d1Size; i++) {
        // Add bidder to bidders
        Bidder b{diagram1[i], i, lambda_};
        b.setPositionInAuction(bidders_.size());
        bidders_.emplace_back(b);
      }
      n_bidders_ = bidders_.size();
    }

    void setGoods(const DiagramType &diagram2) {
      int d2Size = (int)diagram2.size();

      for(int i = 0; i < d2Size; i++) {
        // Add bidder to bidders
        Good g{diagram2[i], i, lambda_};
        goods_.emplace_back(g);
      }
      n_goods_ = goods_.size();
    }

    void buildKDTree() {
      Timer t;
      default_kdt_ = KDTree<double>(true, wasserstein_);
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
      correspondance_kdt_map_
        = kdt_.build(coordinates.data(), goods_.size(), dimension);
    }

    void setEpsilon(double epsilon) {
      epsilon_ = epsilon;
    }

    void initializeEpsilon() {
      double max_persistence = 0;
      for(size_t i = 0; i < bidders_.size(); i++) {
        Bidder &b = bidders_[i];
        double persistence = b.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }

      for(size_t i = 0; i < goods_.size(); i++) {
        Good &g = goods_[i];
        double persistence = g.getPersistence();
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
      for(size_t i = 0; i < goods_.size(); i++) {
        Good &g = goods_[i];
        g.setOwner(-1);
      }
      for(size_t i = 0; i < diagonal_goods_.size(); i++) {
        Good &g = diagonal_goods_[i];
        g.setOwner(-1);
      }
    }

    double getMatchingDistance() {
      double d = 0;
      for(size_t i = 0; i < bidders_.size(); i++) {
        Bidder &b = bidders_[i];
        d += b.cost(b.getProperty(), wasserstein_, geometricalFactor_);
      }
      return d;
    }

    double getRelativePrecision() {
      double d = this->getMatchingDistance();
      if(d < 1e-12) {
        return 0;
      }
      double denominator = d - bidders_.size() * epsilon_;
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
      if(diagonal_goods_.size() == 0) {
        return 0;
      }
      double min_price = std::numeric_limits<double>::max();
      for(size_t i = 0; i < diagonal_goods_.size(); i++) {
        double price = diagonal_goods_[i].getPrice();
        if(price < min_price) {
          min_price = price;
        }
      }
      return min_price;
    }

    void setEpsilonconst(double epsilon) {
      epsilon_ = epsilon;
    }

    template <typename type>
    static type abs(const type var) {
      return (var >= 0) ? var : -var;
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
    bool use_kdt_{true};

  }; // namespace ttk
} // namespace ttk
