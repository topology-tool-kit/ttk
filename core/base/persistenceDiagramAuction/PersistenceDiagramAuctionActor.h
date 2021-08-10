#pragma once

#include <Debug.h>
#include <KDTree.h>
#include <PersistenceDiagramUtils.h>

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>

namespace ttk {

  struct Compare {
    constexpr bool operator()(std::pair<int, double> const &a,
                              std::pair<int, double> const &b) const noexcept {
      return a.second > b.second;
    }
  };

  class PersistenceDiagramAuctionActor : public Debug {
  public:
    PersistenceDiagramAuctionActor() = default;

    PersistenceDiagramAuctionActor(double x, double y, bool is_diagonal, int id)
      : x_{x}, y_{y}, id_{id}, is_diagonal_{is_diagonal} {
    }

    void SetCoordinates(const double x, const double y) {
      this->x_ = x;
      this->y_ = y;
    }

    void SetCriticalCoordinates(const float coords_x,
                                const float coords_y,
                                const float coords_z) {
      this->coords_ = {coords_x, coords_y, coords_z};
    }

    void SetCriticalCoordinates(const std::array<float, 3> coords) {
      this->coords_ = coords;
    }

    std::array<float, 3> GetCriticalCoordinates() const {
      return this->coords_;
    }

    void projectOnDiagonal() {
      x_ = (x_ + y_) / 2;
      y_ = x_;
      is_diagonal_ = true;
    }

    int getId() const {
      return id_;
    };

    double getPersistence() const {
      return this->y_ - this->x_;
    }

    bool isDiagonal() const {
      return this->is_diagonal_;
    }

    double cost(const PersistenceDiagramAuctionActor &g,
                const int wasserstein,
                const double geometricalFactor) const;

    inline double cost(const PersistenceDiagramAuctionActor *g,
                       const int wasserstein,
                       const double geometricalFactor) const {
      return this->cost(*g, wasserstein, geometricalFactor);
    }

    double getPairGeometricalLength(const int wasserstein) const {
      return Geometry::pow(geom_pair_length_[0], wasserstein)
             + Geometry::pow(geom_pair_length_[1], wasserstein)
             + Geometry::pow(geom_pair_length_[2], wasserstein);
    }

    // data members
    double x_{}, y_{};
    int id_{};
    std::array<float, 3> coords_{};

  protected:
    bool is_diagonal_{false};
    std::array<double, 3> geom_pair_length_{};
  };

  class Good : public PersistenceDiagramAuctionActor {
  public:
    Good() = default;
    Good(double x, double y, bool is_diagonal, int id)
      : PersistenceDiagramAuctionActor(x, y, is_diagonal, id) {
    }

    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    Good(const PairTuple &tuple, const int id, const double lambda)
      : PersistenceDiagramAuctionActor{} {
      this->id_ = id;

      const auto type1 = std::get<1>(tuple);
      const auto type2 = std::get<3>(tuple);

      double x = std::get<6>(tuple);
      double y = std::get<10>(tuple);
      this->SetCoordinates(x, y);

      double coords_x;
      double coords_y;
      double coords_z;

      this->geom_pair_length_[0]
        = std::abs(std::get<7>(tuple) - std::get<11>(tuple));
      this->geom_pair_length_[1]
        = std::abs(std::get<8>(tuple) - std::get<12>(tuple));
      this->geom_pair_length_[2]
        = std::abs(std::get<9>(tuple) - std::get<13>(tuple));

      if(type2 == CriticalType::Local_maximum) {
        coords_x
          = (1 - lambda) * std::get<7>(tuple) + lambda * std::get<11>(tuple);
        coords_y
          = (1 - lambda) * std::get<8>(tuple) + lambda * std::get<12>(tuple);
        coords_z
          = (1 - lambda) * std::get<9>(tuple) + lambda * std::get<13>(tuple);
      } else if(type1 == CriticalType::Local_minimum) {
        coords_x
          = lambda * std::get<7>(tuple) + (1 - lambda) * std::get<11>(tuple);
        coords_y
          = lambda * std::get<8>(tuple) + (1 - lambda) * std::get<12>(tuple);
        coords_z
          = lambda * std::get<9>(tuple) + (1 - lambda) * std::get<13>(tuple);
      } else { // pair saddle-saddle
        std::cout << "\n\n\n SADDLE-SADDLE \n\n\n" << '\n';
        coords_x = (std::get<7>(tuple) + std::get<11>(tuple)) / 2;
        coords_y = (std::get<8>(tuple) + std::get<12>(tuple)) / 2;
        coords_z = (std::get<9>(tuple) + std::get<13>(tuple)) / 2;
      }

      this->SetCriticalCoordinates(coords_x, coords_y, coords_z);
      this->is_diagonal_ = x == y;
    }

    Good(const PairTuple &tuple, const int id)
      : PersistenceDiagramAuctionActor{} {
      double x = std::get<6>(tuple);
      double y = std::get<10>(tuple);
      this->SetCoordinates(x, y);
      this->id_ = id;
      this->is_diagonal_ = x == y;
    }

    inline void setPrice(const double price) {
      price_ = price;
    }

    inline double getPrice() const {
      return price_;
    }

    inline void assign(const int b, const double price) {
      price_ = price;
      owner_ = b;
    }

    inline int getOwner() const {
      return owner_;
    }

    inline void setOwner(const int idx) {
      owner_ = idx;
    }

    Good &operator=(const Good &) = default;

  protected:
    double price_{0};
    // Position in Auction.bidders_ of the owner of this good.
    // If the good is not owned, owner_ is set to -1
    int owner_{-1};
  };

  using GoodDiagram = std::vector<Good>;

  class Bidder : public PersistenceDiagramAuctionActor {
  public:
    double diagonal_price_{};

    Bidder() = default;
    Bidder(const double x, const double y, const bool is_diagonal, const int id)
      : PersistenceDiagramAuctionActor{x, y, is_diagonal, id} {
      price_paid_ = 0;
      diagonal_price_ = 0;
      this->is_diagonal_ = is_diagonal;
      this->id_ = id;
    }

    Bidder(const PairTuple &tuple, const int id, const double lambda)
      : PersistenceDiagramAuctionActor{} {

      this->id_ = id;
      const auto type1 = std::get<1>(tuple);
      const auto type2 = std::get<3>(tuple);

      double x = std::get<6>(tuple);
      double y = std::get<10>(tuple);
      this->SetCoordinates(x, y);

      double coords_x;
      double coords_y;
      double coords_z;

      this->geom_pair_length_[0]
        = std::abs(std::get<7>(tuple) - std::get<11>(tuple));
      this->geom_pair_length_[1]
        = std::abs(std::get<8>(tuple) - std::get<12>(tuple));
      this->geom_pair_length_[2]
        = std::abs(std::get<9>(tuple) - std::get<13>(tuple));

      if(type2 == CriticalType::Local_maximum) {
        coords_x
          = (1 - lambda) * std::get<7>(tuple) + lambda * std::get<11>(tuple);
        coords_y
          = (1 - lambda) * std::get<8>(tuple) + lambda * std::get<12>(tuple);
        coords_z
          = (1 - lambda) * std::get<9>(tuple) + lambda * std::get<13>(tuple);
      } else if(type1 == CriticalType::Local_minimum) {
        coords_x
          = lambda * std::get<7>(tuple) + (1 - lambda) * std::get<11>(tuple);
        coords_y
          = lambda * std::get<8>(tuple) + (1 - lambda) * std::get<12>(tuple);
        coords_z
          = lambda * std::get<9>(tuple) + (1 - lambda) * std::get<13>(tuple);
      } else { // pair saddle-saddle
        coords_x = (std::get<7>(tuple) + std::get<11>(tuple)) / 2;
        coords_y = (std::get<8>(tuple) + std::get<12>(tuple)) / 2;
        coords_z = (std::get<9>(tuple) + std::get<13>(tuple)) / 2;
      }

      this->SetCriticalCoordinates(coords_x, coords_y, coords_z);

      this->is_diagonal_ = std::abs(x - y) < Geometry::powIntTen<double>(-12);
      price_paid_ = 0;
      diagonal_price_ = 0;
    }

    // Off-diagonal Bidding (with or without the use of a KD-Tree
    int runBidding(GoodDiagram *goods,
                   Good &diagonalGood,
                   int wasserstein,
                   double epsilon,
                   double geometricalFactor);
    int runKDTBidding(GoodDiagram *goods,
                      Good &diagonalGood,
                      int wasserstein,
                      double epsilon,
                      double geometricalFactor,
                      KDTree<double> *kdt,
                      const int kdt_index = 0);

    // Diagonal Bidding (with or without the use of a KD-Tree
    int runDiagonalBidding(
      GoodDiagram *goods,
      Good &twinGood,
      int wasserstein,
      double epsilon,
      double geometricalFactor,
      std::priority_queue<std::pair<int, double>,
                          std::vector<std::pair<int, double>>,
                          Compare> &diagonal_queue);
    int runDiagonalKDTBidding(
      GoodDiagram *goods,
      Good &diagonalGood,
      int wasserstein,
      double epsilon,
      double geometricalFactor,
      std::vector<KDTree<double> *> &correspondance_kdt_map,
      std::priority_queue<std::pair<int, double>,
                          std::vector<std::pair<int, double>>,
                          Compare> &diagonal_queue,
      const int kdt_index = 0);

    // Utility wrapper functions
    const Good &getProperty() const {
      return this->property_;
    }
    void setDiagonalPrice(const double price) {
      this->diagonal_price_ = price;
    }
    void setPricePaid(const double price) {
      this->price_paid_ = price;
    }
    void setProperty(const Good &g) {
      this->property_ = g;
    }
    void resetProperty() {
      this->property_ = {};
    }
    void setPositionInAuction(const int pos) {
      this->position_in_auction_ = pos;
    }
    int getPositionInAuction() const {
      return this->position_in_auction_;
    }

  protected:
    bool is_diagonal_{};
    double price_paid_{};
    Good property_{};

  private:
    // Attribute stating at which position in Auction.bidders_ this bidder can
    // be found In a single Auction, this attribute could be deducted from id_,
    // but with the use of Barycenter Goods will be added and deleted, which
    // would add and delete diagonal bidders.
    int position_in_auction_{};
  };

  using BidderDiagram = std::vector<Bidder>;

} // namespace ttk
