/// \ingroup base
/// \class ttk::KDTree
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date June 2018
///
/// \brief TTK KD-Tree
///

#pragma once

// base code includes
#include <Debug.h>
#include <Geometry.h> // for pow

#include <algorithm>
#include <limits>
#include <memory>

namespace ttk {
  template <typename dataType, typename Container>
  class KDTree {

  protected:
    // Boolean indicating if the current node is a left node of its parent
    bool is_left_{};
    // Indicates according to which coordinate the tree splits its elements
    int coords_number_{0};
    // Power used for the computation of distances. p=2 yields euclidean
    // distance
    int p_{2};
    // Wether or not the KDTree should include weights that add up to distance
    // for the computation of nearest neighbours
    bool include_weights_{false};

  public:
    using KDTreeRoot = std::unique_ptr<KDTree>;
    using KDTreeMap = std::vector<KDTree *>;

    KDTreeRoot left_{}; // Lower half for the coordinate specified
    KDTreeRoot right_{}; // Higher half
    KDTree *parent_{};
    // ID of the object saved here. The whole object is not kept in the KDTree
    // Users should keep track of them in a table for instance
    int id_{};
    Container coordinates_{};
    Container coords_min_{};
    Container coords_max_{};
    int level_{};

    std::vector<dataType> weight_{};
    std::vector<dataType> min_subweights_{};

    KDTree() = default;
    KDTree(const bool include_weights, const int p)
      : p_{p}, include_weights_{include_weights} {
    }

    KDTree(KDTree *const father, const int coords_number, const bool is_left)
      : is_left_{is_left}, coords_number_{coords_number},
        include_weights_{father->include_weights_}, parent_{father} {
    }

    KDTreeMap build(dataType *data,
                    const int &ptNumber,
                    const int &dimension,
                    const std::vector<std::vector<dataType>> &weights = {},
                    const int weight_number = 1);

    void buildRecursive(dataType *data,
                        std::vector<int> &idx_side,
                        const int &ptNumber,
                        const int &dimension,
                        KDTree<dataType, Container> *parent,
                        KDTreeMap &correspondance_map,
                        const std::vector<std::vector<dataType>> &weights = {},
                        const int weight_number = 1);

    inline void updateWeight(const dataType new_weight,
                             const int weight_index = 0) {
      weight_[weight_index] = new_weight;
      updateMinSubweight(weight_index);
    }

    void updateMinSubweight(const int weight_index = 0);
    void getKClosest(const unsigned int k,
                     const Container &coordinates,
                     KDTreeMap &neighbours,
                     std::vector<dataType> &costs,
                     const int weight_index = 0);

    template <typename PowerFunc>
    void recursiveGetKClosest(const unsigned int k,
                              const Container &coordinates,
                              KDTreeMap &neighbours,
                              std::vector<dataType> &costs,
                              const int weight_index,
                              const PowerFunc &power);

    template <typename PowerFunc>
    inline dataType cost(const Container &coordinates,
                         const PowerFunc &power) const {
      dataType cost = 0;
      for(size_t i = 0; i < coordinates.size(); i++) {
        cost += power(std::abs(coordinates[i] - coordinates_[i]));
      }
      return cost;
    }

    template <typename PowerFunc>
    inline dataType distanceToBox(const KDTree<dataType, Container> &subtree,
                                  const Container &coordinates,
                                  const PowerFunc &power) const {
      dataType d_min = 0;
      for(size_t axis = 0; axis < coordinates.size(); axis++) {
        if(subtree.coords_min_[axis] > coordinates[axis]) {
          d_min += power(subtree.coords_min_[axis] - coordinates[axis]);
        } else if(subtree.coords_max_[axis] < coordinates[axis]) {
          d_min += power(coordinates[axis] - subtree.coords_max_[axis]);
        }
      }
      return d_min;
    }

    inline const std::vector<dataType> &getCoordinates() const {
      return coordinates_;
    }
    inline dataType getWeight(const int weight_index = 0) const {
      return weight_[weight_index];
    }
    inline dataType getMinSubWeight(const int weight_index = 0) const {
      return min_subweights_[weight_index];
    }
    inline bool isLeaf() const {
      return left_ == nullptr && right_ == nullptr;
    }
    inline bool isRoot() const {
      return parent_ == nullptr;
    }
  };
} // namespace ttk

template <typename dataType, typename Container>
typename ttk::KDTree<dataType, Container>::KDTreeMap
  ttk::KDTree<dataType, Container>::build(
    dataType *data,
    const int &ptNumber,
    const int &dimension,
    const std::vector<std::vector<dataType>> &weights,
    const int weight_number) {

  KDTreeMap correspondance_map(ptNumber);
  // First, perform a argsort on the data
  // initialize original index locations
  for(int axis = 0; axis < dimension; axis++) {
    coords_min_[axis] = std::numeric_limits<dataType>::lowest();
    coords_max_[axis] = std::numeric_limits<dataType>::max();
  }
  std::vector<int> idx(ptNumber);
  for(int i = 0; i < ptNumber; i++) {
    idx[i] = i;
  }
  // sort indexes based on comparing values in coordinates
  sort(idx.begin(), idx.end(), [&](int i1, int i2) {
    return data[dimension * i1 + coords_number_]
           < data[dimension * i2 + coords_number_];
  });
  int median_loc = (int)(ptNumber - 1) / 2;
  int median_idx = idx[median_loc];
  correspondance_map[median_idx] = this;

  for(int axis = 0; axis < dimension; axis++) {
    coordinates_[axis] = data[dimension * median_idx + axis];
  }

  id_ = median_idx;
  parent_ = nullptr;
  level_ = 0;

  this->weight_.clear();
  this->min_subweights_.clear();

  if(weights.empty()) {
    this->weight_.resize(weight_number);
    this->min_subweights_.resize(weight_number);
  } else {
    for(int i = 0; i < weight_number; i++) {
      weight_.push_back(weights[i][median_idx]);
      min_subweights_.push_back(weights[i][median_idx]);
    }
  }

  if(idx.size() > 2) {
    // Build left leaf
    std::vector<int> idx_left(median_loc);
    for(int i = 0; i < median_loc; i++) {
      idx_left[i] = idx[i];
    }

    this->left_
      = std::make_unique<KDTree>(this, (coords_number_ + 1) % dimension, true);
    this->left_->buildRecursive(data, idx_left, ptNumber, dimension, this,
                                correspondance_map, weights, weight_number);
  }

  if(idx.size() > 1) {
    // Build right leaf
    std::vector<int> idx_right(ptNumber - median_loc - 1);
    for(int i = 0; i < ptNumber - median_loc - 1; i++) {
      idx_right[i] = idx[i + median_loc + 1];
    }
    this->right_
      = std::make_unique<KDTree>(this, (coords_number_ + 1) % dimension, false);
    this->right_->buildRecursive(data, idx_right, ptNumber, dimension, this,
                                 correspondance_map, weights, weight_number);
  }

  return correspondance_map;
}

template <typename dataType, typename Container>
void ttk::KDTree<dataType, Container>::buildRecursive(
  dataType *data,
  std::vector<int> &idx_side,
  const int &ptNumber,
  const int &dimension,
  KDTree<dataType, Container> *parent,
  KDTreeMap &correspondance_map,
  const std::vector<std::vector<dataType>> &weights,
  const int weight_number) {

  // First, perform a argsort on the data
  sort(idx_side.begin(), idx_side.end(), [&](int i1, int i2) {
    return data[dimension * i1 + coords_number_]
           < data[dimension * i2 + coords_number_];
  });
  int median_loc = (int)(idx_side.size() - 1) / 2;
  int median_idx = idx_side[median_loc];
  correspondance_map[median_idx] = this;

  for(int axis = 0; axis < dimension; axis++) {
    coordinates_[axis] = data[dimension * median_idx + axis];
  }

  id_ = median_idx;
  parent_ = parent;
  level_ = parent->level_ + 1;

  this->weight_.clear();
  this->min_subweights_.clear();

  if(weights.empty()) {
    this->weight_.resize(weight_number);
    this->min_subweights_.resize(weight_number);
  } else {
    for(int i = 0; i < weight_number; i++) {
      weight_.push_back(weights[i][median_idx]);
      min_subweights_.push_back(weights[i][median_idx]);
    }

    if(idx_side.size() > 1) {
      // Once we get to a leaf, update min_subweights of the parents
      for(int w = 0; w < weight_number; w++) {
        this->updateMinSubweight(w);
      }
    }
  }

  // Create bounding box
  for(int axis = 0; axis < dimension; axis++) {
    coords_min_[axis] = parent_->coords_min_[axis];
    coords_max_[axis] = parent_->coords_max_[axis];
  }
  if(is_left_ && !this->isRoot()) {
    coords_max_[parent_->coords_number_]
      = parent_->coordinates_[parent_->coords_number_];
  } else if(!is_left_ && !this->isRoot()) {
    coords_min_[parent_->coords_number_]
      = parent_->coordinates_[parent_->coords_number_];
  }

  if(idx_side.size() > 2) {
    // Build left leaf
    std::vector<int> idx_left(median_loc);
    for(int i = 0; i < median_loc; i++) {
      idx_left[i] = idx_side[i];
    }

    this->left_
      = std::make_unique<KDTree>(this, (coords_number_ + 1) % dimension, true);
    this->left_->buildRecursive(data, idx_left, ptNumber, dimension, this,
                                correspondance_map, weights, weight_number);
  }

  if(idx_side.size() > 1) {
    // Build right leaf
    std::vector<int> idx_right(idx_side.size() - median_loc - 1);
    for(unsigned int i = 0; i < idx_side.size() - median_loc - 1; i++) {
      idx_right[i] = idx_side[i + median_loc + 1];
    }
    this->right_
      = std::make_unique<KDTree>(this, (coords_number_ + 1) % dimension, false);
    this->right_->buildRecursive(data, idx_right, ptNumber, dimension, this,
                                 correspondance_map, weights, weight_number);
  }
}

template <typename dataType, typename Container>
void ttk::KDTree<dataType, Container>::updateMinSubweight(
  const int weight_index) {
  dataType new_min_subweight;
  if(this->isLeaf()) {
    new_min_subweight = weight_[weight_index];
  } else if(!left_) {
    new_min_subweight
      = std::min(right_->min_subweights_[weight_index], weight_[weight_index]);
  } else if(!right_) {
    new_min_subweight
      = std::min(left_->min_subweights_[weight_index], weight_[weight_index]);
  } else {
    new_min_subweight
      = std::min(std::min(left_->min_subweights_[weight_index],
                          right_->min_subweights_[weight_index]),
                 weight_[weight_index]);
  }

  if(new_min_subweight != min_subweights_[weight_index]) {
    min_subweights_[weight_index] = new_min_subweight;
    if(!this->isRoot()) {
      parent_->updateMinSubweight(weight_index);
    }
  }
}

template <typename dataType, typename Container>
void ttk::KDTree<dataType, Container>::getKClosest(const unsigned int k,
                                                   const Container &coordinates,
                                                   KDTreeMap &neighbours,
                                                   std::vector<dataType> &costs,
                                                   const int weight_index) {

  const auto p{this->p_};

  /// Puts the k closest points to the given coordinates in the "neighbours"
  /// vector along with their costs in the "costs" vector The output is not
  /// sorted, if you are interested in the k nearest neighbours in the order,
  /// will need to sort them according to their cost.
  if(this->isLeaf()) {
    dataType cost{};
    TTK_POW_LAMBDA(cost = this->cost, dataType, p, coordinates);
    cost += weight_[weight_index];
    neighbours.push_back(this);
    costs.push_back(cost);
  } else {
    neighbours.reserve(k);
    costs.reserve(k);
    TTK_POW_LAMBDA(this->recursiveGetKClosest, dataType, p, k, coordinates,
                   neighbours, costs, weight_index);
  }
  // TODO sort neighbours and costs !
}

template <typename dataType, typename Container>
template <typename PowerFunc>
void ttk::KDTree<dataType, Container>::recursiveGetKClosest(
  const unsigned int k,
  const Container &coordinates,
  KDTreeMap &neighbours,
  std::vector<dataType> &costs,
  const int weight_index,
  const PowerFunc &power) {
  // 1- Look wether or not to include the current point in the nearest
  // neighbours
  dataType cost = this->cost(coordinates, power);
  cost += weight_[weight_index];

  if(costs.size() < k) {
    neighbours.push_back(this);
    costs.push_back(cost);
  } else {
    // 1.1- Find the most costly amongst neighbours
    const auto idx_max_cost = std::distance(
      costs.begin(), std::max_element(costs.begin(), costs.begin() + k));
    const dataType max_cost = costs[idx_max_cost];

    // 1.2- If the current KDTree is less costly, put it in the neighbours and
    // update costs.
    if(cost < max_cost) {
      costs[idx_max_cost] = cost;
      neighbours[idx_max_cost] = this;
    }
  }

  // 2- Recursively visit KDTrees that are worth it
  if(left_) {
    const dataType max_cost = *std::max_element(costs.begin(), costs.end());
    const dataType min_subweight = left_->min_subweights_[weight_index];
    const dataType d_min = this->distanceToBox(*left_, coordinates, power);
    if(costs.size() < k || d_min + min_subweight < max_cost) {
      // 2.2- It is possible that there exists a point in this subtree that is
      // less costly than max_cost
      left_->recursiveGetKClosest(
        k, coordinates, neighbours, costs, weight_index, power);
    }
  }

  if(right_) {
    const dataType max_cost = *std::max_element(costs.begin(), costs.end());
    const dataType min_subweight = right_->min_subweights_[weight_index];
    const dataType d_min = this->distanceToBox(*right_, coordinates, power);
    if(costs.size() < k || d_min + min_subweight < max_cost) {
      // 2.2- It is possible that there exists a point in this subtree that is
      // less costly than max_cost
      right_->recursiveGetKClosest(
        k, coordinates, neighbours, costs, weight_index, power);
    }
  }
}
