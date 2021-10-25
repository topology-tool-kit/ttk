#ifndef _PERSISTENCEDIAGRAMAUCTIONACTOR_H
#define _PERSISTENCEDIAGRAMAUCTIONACTOR_H

#include <Debug.h>
#include <KDTree.h>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>

namespace ttk {
  template <typename dataType>
  class Good;
  template <typename dataType>
  class Bidder;
  template <typename dataType>
  class GoodDiagram;
  template <typename dataType>
  class BidderDiagram;
  template <typename dataType>
  struct Compare;

  template <typename dataType>
  class PersistenceDiagramAuctionActor : public Debug {
  public:
    dataType x_;
    dataType y_;
    int id_;
    float coords_x_, coords_y_, coords_z_;

    PersistenceDiagramAuctionActor() {
      id_ = 0;
      x_ = 0;
      y_ = 0;
      coords_x_ = 0;
      coords_y_ = 0;
      coords_z_ = 0;
      is_diagonal_ = false;
      // geom_pair_length_.resize(3);
      // geom_pair_length_[0] = 0;
      // geom_pair_length_[1] = 0;
      // geom_pair_length_[2] = 0;
    }

    PersistenceDiagramAuctionActor(dataType x,
                                   dataType y,
                                   bool is_diagonal,
                                   int id) {
      id_ = id;
      x_ = x;
      y_ = y;
      is_diagonal_ = is_diagonal;

      coords_x_ = 0;
      coords_y_ = 0;
      coords_z_ = 0;

      // geom_pair_length_.resize(3);
      // geom_pair_length_[0] = 0;
      // geom_pair_length_[1] = 0;
      // geom_pair_length_[2] = 0;
    }
    ~PersistenceDiagramAuctionActor() = default;

    void SetCoordinates(dataType &x, dataType &y);
    void SetCriticalCoordinates(float coords_x, float coords_y, float coords_z);
    std::tuple<float, float, float> GetCriticalCoordinates() const;
    void projectOnDiagonal();
    int getId() const;
    dataType getPersistence() const;
    bool isDiagonal() const;

    template <typename type>
    inline static type abs(const type var) {
      return (var >= 0) ? var : -var;
    }

    dataType cost(const PersistenceDiagramAuctionActor &g,
                  int &wasserstein,
                  double &geometricalFactor);

    inline dataType cost(const PersistenceDiagramAuctionActor *g,
                         int &wasserstein,
                         double &geometricalFactor);

    double getPairGeometricalLength(const int wasserstein) const;

  protected:
    bool is_diagonal_;
    std::array<double, 3> geom_pair_length_{};
  };

  template <typename dataType>
  int PersistenceDiagramAuctionActor<dataType>::getId() const {
    return id_;
  }

  template <typename dataType>
  void PersistenceDiagramAuctionActor<dataType>::SetCoordinates(dataType &x,
                                                                dataType &y) {
    x_ = x;
    y_ = y;
  }

  template <typename dataType>
  void PersistenceDiagramAuctionActor<dataType>::SetCriticalCoordinates(
    float coords_x, float coords_y, float coords_z) {
    coords_x_ = coords_x;
    coords_y_ = coords_y;
    coords_z_ = coords_z;
  }

  template <typename dataType>
  std::tuple<float, float, float>
    PersistenceDiagramAuctionActor<dataType>::GetCriticalCoordinates() const {
    return std::make_tuple(coords_x_, coords_y_, coords_z_);
  }

  template <typename dataType>
  void PersistenceDiagramAuctionActor<dataType>::projectOnDiagonal() {
    x_ = (x_ + y_) / 2;
    y_ = x_;
    is_diagonal_ = true;
  }

  template <typename dataType>
  dataType PersistenceDiagramAuctionActor<dataType>::getPersistence() const {
    return y_ - x_;
  }

  template <typename dataType>
  bool PersistenceDiagramAuctionActor<dataType>::isDiagonal() const {
    return (is_diagonal_);
  }

  template <typename dataType>
  double PersistenceDiagramAuctionActor<dataType>::getPairGeometricalLength(
    const int wasserstein) const {
    return Geometry::pow(geom_pair_length_[0], wasserstein)
           + Geometry::pow(geom_pair_length_[1], wasserstein)
           + Geometry::pow(geom_pair_length_[2], wasserstein);
  }

  template <typename dataType>
  dataType PersistenceDiagramAuctionActor<dataType>::cost(
    const PersistenceDiagramAuctionActor &g,
    int &wasserstein,
    double &geometricalFactor) {
    if(is_diagonal_ && g.isDiagonal()) {
      return 0;
    } else if(is_diagonal_) {
      return geometricalFactor
               * (2
                  * Geometry::pow(
                    abs<dataType>(g.y_ / 2 - g.x_ / 2), wasserstein))
             + (1 - geometricalFactor) * getPairGeometricalLength(wasserstein);
    } else if(g.isDiagonal()) {
      return geometricalFactor
               * (2
                  * Geometry::pow(abs<dataType>(y_ / 2 - x_ / 2), wasserstein))
             + (1 - geometricalFactor)
                 * g.getPairGeometricalLength(wasserstein);
    } else {
      return geometricalFactor
               * (Geometry::pow(abs<dataType>(x_ - g.x_), wasserstein)
                  + Geometry::pow(abs<dataType>(y_ - g.y_), wasserstein))
             + (1 - geometricalFactor)
                 * (Geometry::pow(
                      abs<dataType>(coords_x_ - g.coords_x_), wasserstein)
                    + Geometry::pow(
                      abs<dataType>(coords_y_ - g.coords_y_), wasserstein)
                    + Geometry::pow(
                      abs<dataType>(coords_z_ - g.coords_z_), wasserstein));
    }
  }

  template <typename dataType>
  dataType PersistenceDiagramAuctionActor<dataType>::cost(
    const PersistenceDiagramAuctionActor *g,
    int &wasserstein,
    double &geometricalFactor) {
    return this->cost(*g, wasserstein, geometricalFactor);
  }

  template <typename dataType>
  class Good : public PersistenceDiagramAuctionActor<dataType> {

  public:
    Good() = default;
    Good(dataType x, dataType y, bool is_diagonal, int id)
      : PersistenceDiagramAuctionActor<dataType>(x, y, is_diagonal, id) {
      PersistenceDiagramAuctionActor<dataType>::is_diagonal_ = is_diagonal;
      PersistenceDiagramAuctionActor<dataType>::id_ = id;
    }

    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    Good(diagramTuple &tuple, int id, double lambda)
      : PersistenceDiagramAuctionActor<dataType>() {
      PersistenceDiagramAuctionActor<dataType>::id_ = id;

      BNodeType type1 = std::get<1>(tuple);
      BNodeType type2 = std::get<3>(tuple);

      dataType x = std::get<6>(tuple);
      dataType y = std::get<10>(tuple);
      this->SetCoordinates(x, y);

      dataType coords_x;
      dataType coords_y;
      dataType coords_z;

      this->geom_pair_length_[0]
        = std::abs(std::get<7>(tuple) - std::get<11>(tuple));
      this->geom_pair_length_[1]
        = std::abs(std::get<8>(tuple) - std::get<12>(tuple));
      this->geom_pair_length_[2]
        = std::abs(std::get<9>(tuple) - std::get<13>(tuple));

      if(type2 == BLocalMax) {
        coords_x
          = (1 - lambda) * std::get<7>(tuple) + lambda * std::get<11>(tuple);
        coords_y
          = (1 - lambda) * std::get<8>(tuple) + lambda * std::get<12>(tuple);
        coords_z
          = (1 - lambda) * std::get<9>(tuple) + lambda * std::get<13>(tuple);
      } else if(type1 == BLocalMin) {
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

      if(x == y) {
        PersistenceDiagramAuctionActor<dataType>::is_diagonal_ = true;
      } else {
        PersistenceDiagramAuctionActor<dataType>::is_diagonal_ = false;
      }
    }

    Good(diagramTuple &tuple, int id)
      : PersistenceDiagramAuctionActor<dataType>() {
      PersistenceDiagramAuctionActor<dataType>::id_ = id;
      dataType x = std::get<6>(tuple);
      dataType y = std::get<10>(tuple);
      this->SetCoordinates(x, y);

      if(x == y) {
        PersistenceDiagramAuctionActor<dataType>::is_diagonal_ = true;
      } else {
        PersistenceDiagramAuctionActor<dataType>::is_diagonal_ = false;
      }
    }

    inline void setPrice(dataType price) {
      price_ = price;
    }

    inline dataType getPrice() const {
      return price_;
    }

    inline void assign(int b, dataType price) {
      price_ = price;
      owner_ = b;
    }

    inline int getOwner() const {
      return owner_;
    }

    inline void setOwner(int idx) {
      owner_ = idx;
    }

    inline Good<dataType> copy() {
      Good<dataType> g = Good<dataType>();
      g.id_ = this->id_;

      g.SetCoordinates(this->x_, this->y_);
      g.SetCriticalCoordinates(
        this->coords_x_, this->coords_y_, this->coords_z_);

      g.is_diagonal_ = this->is_diagonal_;

      g.price_ = this->price_;
      g.owner_ = this->owner_;
      return g;
    }

  protected:
    dataType price_{0};
    // Position in Auction.bidders_ of the owner of this good.
    // If the good is not owned, owner_ is set to -1
    int owner_{-1};
  };

  template <typename dataType>
  class GoodDiagram {

  public:
    GoodDiagram() = default;

    ~GoodDiagram() = default;

    void addGood(Good<dataType> &g);
    Good<dataType> &get(int idx);
    const Good<dataType> &get(int idx) const;
    void set(Good<dataType> &g, int idx);

    inline int size() const {
      return goods_.size();
    }

  private:
    std::vector<Good<dataType>> goods_;
  };

  template <typename dataType>
  void GoodDiagram<dataType>::addGood(Good<dataType> &g) {
    goods_.push_back(g);
  }

  template <typename dataType>
  Good<dataType> &GoodDiagram<dataType>::get(int idx) {
    return goods_[idx];
  }
  template <typename dataType>
  const Good<dataType> &GoodDiagram<dataType>::get(int idx) const {
    return goods_[idx];
  }

  template <typename dataType>
  void GoodDiagram<dataType>::set(Good<dataType> &g, int idx) {
    goods_[idx] = g;
  }

  template <typename dataType>
  class Bidder : public PersistenceDiagramAuctionActor<dataType> {
  public:
    dataType diagonal_price_;

    Bidder() = default;
    Bidder(dataType x, dataType y, bool is_diagonal, int id)
      : PersistenceDiagramAuctionActor<dataType>(x, y, is_diagonal, id) {
      price_paid_ = 0;
      diagonal_price_ = 0;
      PersistenceDiagramAuctionActor<dataType>::is_diagonal_ = is_diagonal;
      PersistenceDiagramAuctionActor<dataType>::id_ = id;
    }

    Bidder(diagramTuple &tuple, int id, double lambda)
      : PersistenceDiagramAuctionActor<dataType>() {
      PersistenceDiagramAuctionActor<dataType>::id_ = id;

      BNodeType type1 = std::get<1>(tuple);
      BNodeType type2 = std::get<3>(tuple);

      dataType x = std::get<6>(tuple);
      dataType y = std::get<10>(tuple);
      this->SetCoordinates(x, y);

      dataType coords_x;
      dataType coords_y;
      dataType coords_z;

      this->geom_pair_length_[0]
        = std::abs(std::get<7>(tuple) - std::get<11>(tuple));
      this->geom_pair_length_[1]
        = std::abs(std::get<8>(tuple) - std::get<12>(tuple));
      this->geom_pair_length_[2]
        = std::abs(std::get<9>(tuple) - std::get<13>(tuple));

      if(type2 == BLocalMax) {
        coords_x
          = (1 - lambda) * std::get<7>(tuple) + lambda * std::get<11>(tuple);
        coords_y
          = (1 - lambda) * std::get<8>(tuple) + lambda * std::get<12>(tuple);
        coords_z
          = (1 - lambda) * std::get<9>(tuple) + lambda * std::get<13>(tuple);
      } else if(type1 == BLocalMin) {
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

      if(std::abs(static_cast<double>(x) - static_cast<double>(y))
         < Geometry::powIntTen<double>(-12)) {
        PersistenceDiagramAuctionActor<dataType>::is_diagonal_ = true;
      } else {
        PersistenceDiagramAuctionActor<dataType>::is_diagonal_ = false;
      }
      price_paid_ = 0;
      diagonal_price_ = 0;
    }

    ~Bidder() {
    }

    // Off-diagonal Bidding (with or without the use of a KD-Tree
    int runBidding(GoodDiagram<dataType> *goods,
                   Good<dataType> &diagonalGood,
                   int wasserstein,
                   dataType epsilon,
                   double geometricalFactor);
    int runKDTBidding(GoodDiagram<dataType> *goods,
                      Good<dataType> &diagonalGood,
                      int wasserstein,
                      dataType epsilon,
                      double geometricalFactor,
                      KDTree<dataType> *kdt,
                      const int kdt_index = 0);

    // Diagonal Bidding (with or without the use of a KD-Tree
    int runDiagonalBidding(
      GoodDiagram<dataType> *goods,
      Good<dataType> &twinGood,
      int wasserstein,
      dataType epsilon,
      double geometricalFactor,
      std::priority_queue<std::pair<int, dataType>,
                          std::vector<std::pair<int, dataType>>,
                          Compare<dataType>> &diagonal_queue);
    int runDiagonalKDTBidding(
      GoodDiagram<dataType> *goods,
      Good<dataType> &diagonalGood,
      int wasserstein,
      dataType epsilon,
      double geometricalFactor,
      std::vector<KDTree<dataType> *> &correspondance_kdt_map,
      std::priority_queue<std::pair<int, dataType>,
                          std::vector<std::pair<int, dataType>>,
                          Compare<dataType>> &diagonal_queue,
      const int kdt_index = 0);

    // Utility wrapper functions
    const Good<dataType> &getProperty();
    void setDiagonalPrice(dataType price);
    void setPricePaid(dataType price);
    void setProperty(const Good<dataType> &g);
    void resetProperty();
    void setTwin(Good<dataType> *g);
    void setPositionInAuction(int pos);
    int getPositionInAuction();

  protected:
    bool is_diagonal_;
    dataType price_paid_;
    Good<dataType> property_{};

  private:
    // Attribute stating at which position in Auction.bidders_ this bidder can
    // be found In a single Auction, this attribute could be deducted from id_,
    // but with the use of Barycenter Goods will be added and deleted, which
    // would add and delete diagonal bidders.
    int position_in_auction_;
  };

  template <typename dataType>
  const Good<dataType> &Bidder<dataType>::getProperty() {
    return property_;
  }

  template <typename dataType>
  int Bidder<dataType>::runBidding(GoodDiagram<dataType> *goods,
                                   Good<dataType> &twinGood,
                                   int wasserstein,
                                   dataType epsilon,
                                   double geometricalFactor) {
    dataType best_val = std::numeric_limits<dataType>::lowest();
    dataType second_val = std::numeric_limits<dataType>::lowest();
    Good<dataType> *best_good{};
    for(int i = 0; i < goods->size(); i++) {
      Good<dataType> &g = goods->get(i);
      dataType val = -this->cost(g, wasserstein, geometricalFactor);
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
    Good<dataType> &g = twinGood;
    dataType val = -this->cost(g, wasserstein, geometricalFactor);
    val -= g.getPrice();
    if(val > best_val) {
      second_val = best_val;
      best_val = val;
      best_good = &g;
    } else if(val > second_val) {
      second_val = val;
    }

    if(second_val == std::numeric_limits<dataType>::lowest()) {
      // There is only one acceptable good for the bidder
      second_val = best_val;
    }
    if(best_good == nullptr) {
      return -1;
    }
    dataType old_price = best_good->getPrice();
    dataType new_price = old_price + best_val - second_val + epsilon;
    if(new_price > std::numeric_limits<dataType>::max() / 2) {
      new_price = old_price + epsilon;
      std::cout << "Huho 376" << std::endl;
    }
    // Assign bidder to best_good
    this->setProperty(*best_good);
    this->setPricePaid(new_price);

    // Assign best_good to bidder and unassign the previous owner of best_good
    // if need be
    int idx_reassigned = best_good->getOwner();
    best_good->assign(this->position_in_auction_, new_price);
    return idx_reassigned;
  }

  template <typename dataType>
  int Bidder<dataType>::runDiagonalBidding(
    GoodDiagram<dataType> *goods,
    Good<dataType> &twinGood,
    int wasserstein,
    dataType epsilon,
    double geometricalFactor,
    std::priority_queue<std::pair<int, dataType>,
                        std::vector<std::pair<int, dataType>>,
                        Compare<dataType>> &diagonal_queue) {

    // First, find the lowest and second lowest weights for diagonal goods
    // Take this time to update the weights in the priority queue of the goods
    // tested. It is not necessary to update weights for all diagonal bidders in
    // the pririty queue since weights can only increase and we are interested
    // only in the lowest ones
    bool updated_top_pair
      = false; // Boolean which equals true iff the top pair in the priority
               // queue is given the good price
    bool non_empty_goods = goods->size() > 0;
    std::pair<int, dataType> best_pair;
    if(non_empty_goods) {
      while(!updated_top_pair) {
        std::pair<int, dataType> top_pair = diagonal_queue.top();
        diagonal_queue.pop();

        dataType queue_weight = top_pair.second;
        Good<dataType> *good = &(goods->get(top_pair.first));
        if(good->getPrice() > queue_weight) {
          // If the weight in the priority queue is not the good one, update
          std::get<1>(top_pair) = good->getPrice();
          diagonal_queue.push(top_pair);
        } else {
          updated_top_pair = true;
          best_pair = top_pair;
        }
      }
    }

    std::pair<int, dataType> second_pair;
    if(non_empty_goods && diagonal_queue.size() > 0) {
      bool updated_second_pair = false;
      while(!updated_second_pair) {
        second_pair = diagonal_queue.top();
        dataType queue_weight = second_pair.second;
        Good<dataType> *good = &(goods->get(second_pair.first));
        if(good->getPrice() != queue_weight) {
          // If the weight in the priority queue is not the good one, update it
          diagonal_queue.pop();
          std::get<1>(second_pair) = good->getPrice();
          diagonal_queue.push(second_pair);
        } else {
          updated_second_pair = true;
        }
      }
    }

    Good<dataType> *best_good{};
    dataType best_val = 0;
    dataType second_val = 0;
    if(non_empty_goods) {
      best_val = -best_pair.second;
      second_val = diagonal_queue.size() > 0 ? -second_pair.second : best_val;
      best_good = &(goods->get(best_pair.first));
    }

    // And now check for the corresponding twin bidder
    bool is_twin = false;
    Good<dataType> &g = twinGood;
    dataType val = -this->cost(g, wasserstein, geometricalFactor);
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

    if(second_val == std::numeric_limits<dataType>::lowest()) {
      // There is only one acceptable good for the bidder
      second_val = best_val;
      std::cout << "Huho 474" << std::endl;
    }
    dataType old_price = best_good->getPrice();
    dataType new_price = old_price + best_val - second_val + epsilon;
    if(new_price > std::numeric_limits<dataType>::max() / 2) {
      new_price = old_price + epsilon;
    }

    // Assign bidder to best_good
    this->setProperty(*best_good);
    this->setPricePaid(new_price);

    // Assign best_good to bidder and unassign the previous owner of best_good
    // if need be
    int idx_reassigned = best_good->getOwner();
    best_good->assign(this->position_in_auction_, new_price);
    if(!is_twin) {
      std::get<1>(best_pair) = new_price;
    }
    if(non_empty_goods) {
      diagonal_queue.push(best_pair);
    }
    return idx_reassigned;
  }

  template <typename dataType>
  int Bidder<dataType>::runDiagonalKDTBidding(
    GoodDiagram<dataType> *goods,
    Good<dataType> &twinGood,
    int wasserstein,
    dataType epsilon,
    double geometricalFactor,
    std::vector<KDTree<dataType> *> &correspondance_kdt_map,
    std::priority_queue<std::pair<int, dataType>,
                        std::vector<std::pair<int, dataType>>,
                        Compare<dataType>> &diagonal_queue,
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
    bool non_empty_goods = goods->size() > 0;
    std::pair<int, dataType> best_pair;
    if(non_empty_goods) {
      while(!updated_top_pair) {
        std::pair<int, dataType> top_pair = diagonal_queue.top();
        dataType queue_weight = top_pair.second;
        Good<dataType> *good = &(goods->get(top_pair.first));

        diagonal_queue.pop();
        if(good->getPrice() > queue_weight) {
          // If the weight in the priority queue is not the good one, update it
          std::get<1>(top_pair) = good->getPrice();
          diagonal_queue.push(top_pair);
        } else {
          updated_top_pair = true;
          best_pair = top_pair;
        }
      }
    }
    // std::cout << "got here 1" << std::endl;
    std::pair<int, dataType> second_pair;
    if(non_empty_goods && diagonal_queue.size() > 0) {
      bool updated_second_pair = false;
      while(!updated_second_pair) {
        second_pair = diagonal_queue.top();
        dataType queue_weight = second_pair.second;
        Good<dataType> *good = &(goods->get(second_pair.first));
        if(good->getPrice() > queue_weight) {
          // If the weight in the priority queue is not the good one, update it
          diagonal_queue.pop();
          std::get<1>(second_pair) = good->getPrice();
          diagonal_queue.push(second_pair);
        } else {
          updated_second_pair = true;
        }
      }
    }

    // std::cout << "got here 1/2" << std::endl;
    Good<dataType> *best_good{};
    dataType best_val = 0;
    dataType second_val = 0;
    if(non_empty_goods) {
      best_val = -best_pair.second;
      second_val = diagonal_queue.size() > 0 ? -second_pair.second : best_val;
      best_good = &(goods->get(best_pair.first));
    }
    // std::cout << "got here 2" << std::endl;

    // And now check for the corresponding twin bidder
    bool is_twin = false;
    Good<dataType> &g = twinGood;
    dataType val = -this->cost(g, wasserstein, geometricalFactor);
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

    // std::cout << "got here 3" << std::endl;
    if(second_val == std::numeric_limits<dataType>::lowest()) {
      // There is only one acceptable good for the bidder
      second_val = best_val;
    }
    dataType old_price = best_good->getPrice();
    dataType new_price = old_price + best_val - second_val + epsilon;
    if(new_price > std::numeric_limits<dataType>::max() / 2) {
      new_price = old_price + epsilon;
      std::cout << "Huho 592" << std::endl;
    }

    // Assign bidder to best_good
    this->setProperty(*best_good);
    this->setPricePaid(new_price);

    // std::cout << "got here 4" << std::endl;
    // Assign best_good to bidder and unassign the previous owner of best_good
    // if need be
    int idx_reassigned = best_good->getOwner();
    best_good->assign(this->position_in_auction_, new_price);
    if(is_twin) {
      // std::cout << "got here 5" << std::endl;
      // Update weight in KDTree if the closest good is in it
      correspondance_kdt_map[best_good->id_]->updateWeight(
        new_price, kdt_index);
      if(non_empty_goods) {
        diagonal_queue.push(best_pair);
      }
    } else {
      // std::cout << "got here 6" << std::endl;
      if(non_empty_goods) {
        std::get<1>(best_pair) = new_price;
        diagonal_queue.push(best_pair);
      }
    }
    // std::cout << "left that function" << std::endl;
    return idx_reassigned;
  }

  template <typename dataType>
  int Bidder<dataType>::runKDTBidding(GoodDiagram<dataType> *goods,
                                      Good<dataType> &twinGood,
                                      int wasserstein,
                                      dataType epsilon,
                                      double geometricalFactor,
                                      KDTree<dataType> *kdt,
                                      const int kdt_index) {
    /// Runs bidding of a non-diagonal bidder
    std::vector<KDTree<dataType> *> neighbours;
    std::vector<dataType> costs;

    std::vector<dataType> coordinates;
    coordinates.push_back(geometricalFactor * this->x_);
    coordinates.push_back(geometricalFactor * this->y_);
    if(geometricalFactor < 1) {
      coordinates.push_back((1 - geometricalFactor) * this->coords_x_);
      coordinates.push_back((1 - geometricalFactor) * this->coords_y_);
      coordinates.push_back((1 - geometricalFactor) * this->coords_z_);
    }

    // std::cout<<"got to 1"<<std::endl;
    kdt->getKClosest(2, coordinates, neighbours, costs, kdt_index);
    // std::cout<<"got to 2"<<std::endl;
    dataType best_val, second_val;
    KDTree<dataType> *closest_kdt;
    Good<dataType> *best_good{};
    if(costs.size() == 2) {
      // std::cout<<"got to 735"<<std::endl;
      std::vector<int> idx(2);
      idx[0] = 0;
      idx[1] = 1;
      sort(idx.begin(), idx.end(),
           [&costs](int &a, int &b) { return costs[a] < costs[b]; });

      closest_kdt = neighbours[idx[0]];
      best_good = &(goods->get(closest_kdt->id_));
      // Value is defined as the opposite of cost (each bidder aims at
      // maximizing it)
      best_val = -costs[idx[0]];
      second_val = -costs[idx[1]];
    } else {
      // std::cout<<"got to 748"<<std::endl;
      // If the kdtree contains only one point
      closest_kdt = neighbours[0];
      best_good = &(goods->get(closest_kdt->id_));
      best_val = -costs[0];
      second_val = best_val;
    }
    // std::cout<<"got to 755"<<std::endl;
    // And now check for the corresponding twin bidder
    bool twin_chosen = false;
    Good<dataType> &g = twinGood;
    dataType val = -this->cost(g, wasserstein, geometricalFactor);
    val -= g.getPrice();
    if(val > best_val) {
      second_val = best_val;
      best_val = val;
      best_good = &g;
      twin_chosen = true;
    } else if(val > second_val) {
      second_val = val;
    }

    // std::cout << "got here 2" << std::endl;
    if(second_val == std::numeric_limits<dataType>::lowest()) {
      // There is only one acceptable good for the bidder
      second_val = best_val;
    }
    dataType old_price = best_good->getPrice();
    dataType new_price = old_price + best_val - second_val + epsilon;
    if(new_price > std::numeric_limits<dataType>::max() / 2) {
      new_price = old_price + epsilon;
      std::cout << "Huho 681" << std::endl;
    }
    // Assign bidder to best_good
    this->setProperty(*best_good);
    this->setPricePaid(new_price);
    // Assign best_good to bidder and unassign the previous owner of best_good
    // if need be
    int idx_reassigned = best_good->getOwner();

    // std::cout << "got here 3" << std::endl;
    best_good->assign(this->position_in_auction_, new_price);
    // Update the price in the KDTree
    if(!twin_chosen) {
      closest_kdt->updateWeight(new_price, kdt_index);
    }
    return idx_reassigned;
  }

  template <typename dataType>
  void Bidder<dataType>::setDiagonalPrice(dataType price) {
    diagonal_price_ = price;
  }

  template <typename dataType>
  void Bidder<dataType>::setPricePaid(dataType price) {
    price_paid_ = price;
  }

  template <typename dataType>
  void Bidder<dataType>::setProperty(const Good<dataType> &g) {
    property_ = g;
  }

  template <typename dataType>
  void Bidder<dataType>::resetProperty() {
    property_ = {};
  }

  template <typename dataType>
  void Bidder<dataType>::setPositionInAuction(int pos) {
    position_in_auction_ = pos;
  }

  template <typename dataType>
  int Bidder<dataType>::getPositionInAuction() {
    return position_in_auction_;
  }

  template <typename dataType>
  class BidderDiagram {

  public:
    std::vector<Bidder<dataType>> bidders_;

    BidderDiagram() = default;

    ~BidderDiagram() = default;

    inline void addBidder(Bidder<dataType> &b) {
      bidders_.push_back(b);
    }

    inline int size() const {
      return bidders_.size();
    }

    inline Bidder<dataType> &get(int idx) {
      return bidders_[idx];
    }
    inline const Bidder<dataType> &get(int idx) const {
      return bidders_[idx];
    }
  };

} // namespace ttk

#endif
