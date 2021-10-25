#pragma once

#ifndef matchingTuple
#define matchingTuple std::tuple<ttk::SimplexId, ttk::SimplexId, dataType>
#endif

#ifndef diagramTuple
#define diagramTuple                                                       \
  std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,            \
             ttk::CriticalType, dataType, ttk::SimplexId, dataType, float, \
             float, float, dataType, float, float, float>
#endif

#ifndef BNodeType
#define BNodeType ttk::CriticalType
#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1 ttk::CriticalType::Saddle1
#define BSaddle2 ttk::CriticalType::Saddle2
#define BIdVertex int
#endif

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
  template <typename dataType>
  struct Compare {
    constexpr bool
      operator()(std::pair<int, dataType> const &a,
                 std::pair<int, dataType> const &b) const noexcept {
      return a.second > b.second;
    }
  };

  template <typename dataType>
  class PersistenceDiagramAuction : public Debug {

  public:
    inline int getAugmentedNumberOfBidders() {
      return bidders_.size();
    }

    KDTree<dataType> default_kdt_{};
    KDTree<dataType> &kdt_{default_kdt_};
    std::vector<KDTree<dataType> *> default_correspondance_kdt_map_{};
    std::vector<KDTree<dataType> *> &correspondance_kdt_map_{
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
      BidderDiagram<dataType> &bidders,
      GoodDiagram<dataType> &goods,
      int wasserstein,
      double geometricalFactor,
      double lambda,
      double delta_lim,
      KDTree<dataType> &kdt,
      std::vector<KDTree<dataType> *> &correspondance_kdt_map,
      dataType epsilon = {},
      dataType initial_diag_price = {},
      bool use_kdTree = true)
      : kdt_{kdt}, correspondance_kdt_map_{correspondance_kdt_map},
        bidders_{bidders}, goods_{goods} {

      n_bidders_ = bidders.size();
      n_goods_ = goods.size();

      for(int i = 0; i < n_bidders_; i++) {
        // Add diagonal goods
        Bidder<dataType> &b = bidders_.get(i);
        Good<dataType> g = Good<dataType>(b.x_, b.y_, true, -b.id_ - 1);
        g.projectOnDiagonal();
        if(b.diagonal_price_ > 0) {
          g.setPrice(b.diagonal_price_);
        } else {
          g.setPrice(initial_diag_price);
        }
        diagonal_goods_.addGood(g);
        std::pair<int, dataType> pair = std::make_pair(i, g.getPrice());
        diagonal_queue_.push(pair);
      }
      for(int i = 0; i < n_goods_; i++) {
        // Add diagonal bidders
        Good<dataType> &g = goods_.get(i);
        Bidder<dataType> b = Bidder<dataType>(g.x_, g.y_, true, -g.id_ - 1);
        b.projectOnDiagonal();
        b.setPositionInAuction(bidders_.size());
        bidders_.addBidder(b);
      }

      epsilon_ = epsilon;
      wasserstein_ = wasserstein;
      delta_lim_ = delta_lim;
      geometricalFactor_ = geometricalFactor;
      lambda_ = lambda;

      use_kdt_ = (use_kdTree && goods_.size() > 0);
    }

    void runAuctionRound(int &n_biddings, const int kdt_index = 0);
    dataType getMatchingsAndDistance(std::vector<matchingTuple> *matchings,
                                     bool get_diagonal_matches = false);
    dataType run(std::vector<matchingTuple> *matchings);
    dataType run() {
      std::vector<matchingTuple> matchings{};
      return this->run(&matchings);
    }
    dataType getMaximalPrice();

    void BuildAuctionDiagrams(const BidderDiagram<dataType> *BD,
                              const GoodDiagram<dataType> *GD) {
      n_bidders_ = BD->size();
      n_goods_ = GD->size();
      // delete_kdTree_ = false;
      bidders_ = *BD;
      goods_ = *GD;

      for(int i = 0; i < n_bidders_; i++) {
        // Add diagonal goods
        Bidder<dataType> &b = bidders_.get(i);
        Good<dataType> g = Good<dataType>(b.x_, b.y_, true, -b.id_ - 1);
        g.projectOnDiagonal();
        diagonal_goods_.addGood(g);
        std::pair<int, dataType> pair = std::make_pair(i, g.getPrice());
        diagonal_queue_.push(pair);
      }
      for(int i = 0; i < n_goods_; i++) {
        // Add diagonal bidders
        Good<dataType> &g = goods_.get(i);
        Bidder<dataType> b = Bidder<dataType>(g.x_, g.y_, true, -g.id_ - 1);
        b.projectOnDiagonal();
        b.setPositionInAuction(bidders_.size());
        bidders_.addBidder(b);
      }
      if(goods_.size() > 0) {
        // use_kdt_ = use_kdt_;
        this->buildKDTree();
      } else {
        use_kdt_ = false;
      }
    }

    void BuildAuctionDiagrams(std::vector<diagramTuple> diagram1,
                              std::vector<diagramTuple> diagram2) {
      n_bidders_ = diagram1.size();
      n_goods_ = diagram2.size();
      this->setBidders(diagram1);
      this->setGoods(diagram2);
      for(int i = 0; i < n_bidders_; i++) {
        // Add diagonal goods
        Bidder<dataType> &b = bidders_.get(i);
        Good<dataType> g = Good<dataType>(b.x_, b.y_, true, -b.id_ - 1);
        g.projectOnDiagonal();
        diagonal_goods_.addGood(g);
        std::pair<int, dataType> pair = std::make_pair(i, g.getPrice());
        diagonal_queue_.push(pair);
      }
      for(int i = 0; i < n_goods_; i++) {
        // Add diagonal bidders
        Good<dataType> &g = goods_.get(i);
        Bidder<dataType> b = Bidder<dataType>(g.x_, g.y_, true, -g.id_ - 1);
        b.projectOnDiagonal();
        b.setPositionInAuction(bidders_.size());
        bidders_.addBidder(b);
      }
      if(bidders_.size() > 0) {
        use_kdt_ = true;
        this->buildKDTree();
      } else {
        use_kdt_ = false;
      }
    }

    void setBidders(std::vector<diagramTuple> diagram1) {
      int d1Size = (int)diagram1.size();

      for(int i = 0; i < d1Size; i++) {
        // Add bidder to bidders
        Bidder<dataType> b = Bidder<dataType>(diagram1[i], i, lambda_);
        b.setPositionInAuction(bidders_.size());
        bidders_.addBidder(b);
      }
      n_bidders_ = bidders_.size();
    }

    void setGoods(std::vector<diagramTuple> diagram2) {
      int d2Size = (int)diagram2.size();

      for(int i = 0; i < d2Size; i++) {
        // Add bidder to bidders
        Good<dataType> g = Good<dataType>(diagram2[i], i, lambda_);
        goods_.addGood(g);
      }
      n_goods_ = goods_.size();
    }

    void buildKDTree() {
      Timer t;
      default_kdt_ = KDTree<dataType>(true, wasserstein_);
      const int dimension
        = geometricalFactor_ >= 1 ? (geometricalFactor_ <= 0 ? 3 : 2) : 5;
      std::vector<dataType> coordinates;
      for(int i = 0; i < goods_.size(); i++) {
        Good<dataType> &g = goods_.get(i);
        if(geometricalFactor_ > 0) {
          coordinates.push_back(geometricalFactor_ * g.x_);
          coordinates.push_back(geometricalFactor_ * g.y_);
        }
        if(geometricalFactor_ < 1) {
          coordinates.push_back((1 - geometricalFactor_) * g.coords_x_);
          coordinates.push_back((1 - geometricalFactor_) * g.coords_y_);
          coordinates.push_back((1 - geometricalFactor_) * g.coords_z_);
        }
      }
      correspondance_kdt_map_
        = kdt_.build(coordinates.data(), goods_.size(), dimension);
    }

    void setEpsilon(dataType epsilon) {
      epsilon_ = epsilon;
    }

    void initializeEpsilon() {
      dataType max_persistence = 0;
      for(int i = 0; i < bidders_.size(); i++) {
        Bidder<dataType> &b = bidders_.get(i);
        dataType persistence = b.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }

      for(int i = 0; i < goods_.size(); i++) {
        Good<dataType> &g = goods_.get(i);
        dataType persistence = g.getPersistence();
        if(persistence > max_persistence) {
          max_persistence = persistence;
        }
      }
      this->epsilon_ = 5.0 / 4 * Geometry::pow(max_persistence, wasserstein_);
    }

    void buildUnassignedBidders() {
      for(int i = 0; i < bidders_.size(); i++) {
        Bidder<dataType> &b = bidders_.get(i);
        b.resetProperty();
        unassignedBidders_.push(i);
      }
    }

    void reinitializeGoods() {
      for(int i = 0; i < goods_.size(); i++) {
        Good<dataType> &g = goods_.get(i);
        g.setOwner(-1);
      }
      for(int i = 0; i < diagonal_goods_.size(); i++) {
        Good<dataType> &g = diagonal_goods_.get(i);
        g.setOwner(-1);
      }
    }

    dataType getMatchingDistance() {
      dataType d = 0;
      for(int i = 0; i < bidders_.size(); i++) {
        Bidder<dataType> &b = bidders_.get(i);
        d += b.cost(b.getProperty(), wasserstein_, geometricalFactor_);
      }
      return d;
    }

    dataType getRelativePrecision() {
      dataType d = this->getMatchingDistance();
      if(d < 1e-12) {
        return 0;
      }
      dataType denominator = d - bidders_.size() * epsilon_;
      if(denominator <= 0) {
        return 1;
      } else {
        return Geometry::pow(d / denominator, 1 / ((float)wasserstein_)) - 1;
      }
    }

    void updateDiagonalPrices() {
      for(int i = 0; i < diagonal_goods_.size(); i++) {
        Bidder<dataType> &b = bidders_.get(i);
        b.setDiagonalPrice(diagonal_goods_.get(i).getPrice());
      }
    }

    dataType getMinimalDiagonalPrice() {
      if(diagonal_goods_.size() == 0) {
        return 0;
      }
      dataType min_price = std::numeric_limits<dataType>::max();
      for(int i = 0; i < diagonal_goods_.size(); i++) {
        dataType price = diagonal_goods_.get(i).getPrice();
        if(price < min_price) {
          min_price = price;
        }
      }
      return min_price;
    }

    void setEpsilonconst(dataType epsilon) {
      epsilon_ = epsilon;
    }

    template <typename type>
    static type abs(const type var) {
      return (var >= 0) ? var : -var;
    }

  protected:
    int wasserstein_{2}; // Power in Wassertsein distance (by default set to 2)
    BidderDiagram<dataType> default_bidders_{};
    BidderDiagram<dataType> &bidders_{default_bidders_};
    GoodDiagram<dataType> default_goods_{};
    GoodDiagram<dataType> &goods_{default_goods_};
    GoodDiagram<dataType> diagonal_goods_;
    std::priority_queue<std::pair<int, dataType>,
                        std::vector<std::pair<int, dataType>>,
                        Compare<dataType>>
      diagonal_queue_{};
    std::queue<int> unassignedBidders_{};

    int n_bidders_{0};
    int n_goods_{0};

    dataType epsilon_{1};
    double geometricalFactor_{};
    double lambda_{};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double delta_lim_{};
    bool use_kdt_{true};

    // KDTree<dataType>* kdt_;
  }; // namespace ttk
} // namespace ttk

#include <PersistenceDiagramAuctionImpl.h>
