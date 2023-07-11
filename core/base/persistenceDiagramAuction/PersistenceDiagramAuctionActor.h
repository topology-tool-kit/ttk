#pragma once

#include <Debug.h>
#include <KDTree.h>
#include <PersistenceDiagramUtils.h>

#include <array>
#include <cmath>
#include <queue>

namespace ttk {

  struct Compare {
    constexpr bool operator()(std::pair<int, double> const &a,
                              std::pair<int, double> const &b) const noexcept {
      return a.second > b.second;
    }
  };

  class PersistenceDiagramAuctionActor {
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

    void GetKDTCoordinates(double geometricalFactor,
                           std::array<double, 5> &coordinates) const {
      coordinates[0] = geometricalFactor * this->x_;
      coordinates[1] = geometricalFactor * this->y_;
      if(geometricalFactor < 1) {
        coordinates[2] = (1 - geometricalFactor) * this->coords_[0];
        coordinates[3] = (1 - geometricalFactor) * this->coords_[1];
        coordinates[4] = (1 - geometricalFactor) * this->coords_[2];
      } else {
        coordinates[2] = 0;
        coordinates[3] = 0;
        coordinates[4] = 0;
      }
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
    Good(const Good &) = default;
    Good(double x, double y, bool is_diagonal, int id)
      : PersistenceDiagramAuctionActor(x, y, is_diagonal, id) {
    }

    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    Good(const PersistencePair &pair, const int id, const double lambda) {

      this->id_ = id;
      this->SetCoordinates(pair.birth.sfValue, pair.death.sfValue);
      this->geom_pair_length_ = {
        std::abs(pair.birth.coords[0] - pair.death.coords[0]),
        std::abs(pair.birth.coords[1] - pair.death.coords[1]),
        std::abs(pair.birth.coords[2] - pair.death.coords[2]),
      };

      std::array<float, 3> coords{};
      const float lb = lambda;

      if(pair.death.type == CriticalType::Local_maximum) {
        coords = {
          (1.0f - lb) * pair.birth.coords[0] + lb * pair.death.coords[0],
          (1.0f - lb) * pair.birth.coords[1] + lb * pair.death.coords[1],
          (1.0f - lb) * pair.birth.coords[2] + lb * pair.death.coords[2],
        };
      } else if(pair.birth.type == CriticalType::Local_minimum) {
        coords = {
          lb * pair.birth.coords[0] + (1.0f - lb) * pair.death.coords[0],
          lb * pair.birth.coords[1] + (1.0f - lb) * pair.death.coords[1],
          lb * pair.birth.coords[2] + (1.0f - lb) * pair.death.coords[2],
        };
      } else { // pair saddle-saddle
        coords = {
          (pair.birth.coords[0] + pair.death.coords[0]) / 2.0f,
          (pair.birth.coords[1] + pair.death.coords[1]) / 2.0f,
          (pair.birth.coords[2] + pair.death.coords[2]) / 2.0f,
        };
      }

      this->SetCriticalCoordinates(coords);
      this->is_diagonal_ = (pair.birth.sfValue == pair.death.sfValue);
    }

    Good(const PersistencePair &pair, const int id) {
      this->SetCoordinates(pair.birth.sfValue, pair.death.sfValue);
      this->id_ = id;
      this->is_diagonal_ = (pair.birth.sfValue == pair.death.sfValue);
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

    Bidder(const PersistencePair &pair, const int id, const double lambda)
      : diagonal_price_{}, price_paid_{} {

      this->id_ = id;
      this->SetCoordinates(pair.birth.sfValue, pair.death.sfValue);
      this->geom_pair_length_ = {
        std::abs(pair.birth.coords[0] - pair.death.coords[0]),
        std::abs(pair.birth.coords[1] - pair.death.coords[1]),
        std::abs(pair.birth.coords[2] - pair.death.coords[2]),
      };

      std::array<float, 3> coords{};
      const float lb = lambda;

      if(pair.death.type == CriticalType::Local_maximum) {
        coords = {
          (1.0f - lb) * pair.birth.coords[0] + lb * pair.death.coords[0],
          (1.0f - lb) * pair.birth.coords[1] + lb * pair.death.coords[1],
          (1.0f - lb) * pair.birth.coords[2] + lb * pair.death.coords[2],
        };
      } else if(pair.birth.type == CriticalType::Local_minimum) {
        coords = {
          lb * pair.birth.coords[0] + (1.0f - lb) * pair.death.coords[0],
          lb * pair.birth.coords[1] + (1.0f - lb) * pair.death.coords[1],
          lb * pair.birth.coords[2] + (1.0f - lb) * pair.death.coords[2],
        };
      } else { // pair saddle-saddle
        coords = {
          (pair.birth.coords[0] + pair.death.coords[0]) / 2.0f,
          (pair.birth.coords[1] + pair.death.coords[1]) / 2.0f,
          (pair.birth.coords[2] + pair.death.coords[2]) / 2.0f,
        };
      }

      this->SetCriticalCoordinates(coords);
      this->is_diagonal_ = std::abs(pair.birth.sfValue - pair.death.sfValue)
                           < Geometry::powIntTen<double>(-12);
    }

    using KDT = KDTree<double, std::array<double, 5>>;

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
                      KDT *kdt,
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
      std::vector<KDT *> &correspondence_kdt_map,
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
