/// \ingroup base
/// \class ttk::KDTree
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date June 2018
///
/// \brief TTK KD-Tree
///

#pragma once

#ifndef _KDTREE_H
#define _KDTREE_H

// base code includes
#include <Debug.h>
#include <Geometry.h> // for pow
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace ttk {
  template <typename dataType>
  class KDTree : public Debug {

  protected:
    // Boolean indicating if the current node is a left node of its parent
    bool is_left_;
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
    std::vector<dataType> coordinates_{};
    std::vector<dataType> coords_min_{};
    std::vector<dataType> coords_max_{};
    int level_{};

    std::vector<dataType> weight_{};
    std::vector<dataType> min_subweights_{};

    KDTree() = default;
    KDTree(bool include_weights, int p)
      : p_{p}, include_weights_{include_weights} {
    }

    KDTree(KDTree *father, int coords_number, bool is_left)
      : is_left_{is_left}, coords_number_{coords_number},
        include_weights_{father->include_weights_}, parent_{father} {
    }

    KDTreeMap build(dataType *coordinates,
                    const int &ptNumber,
                    const int &dimension,
                    const int weight_number = 1);
    KDTreeMap build(dataType *coordinates,
                    const int &ptNumber,
                    const int &dimension,
                    std::vector<std::vector<dataType>> &weights,
                    const int weight_number = 1);

    void buildRecursive(dataType *coordinates,
                        std::vector<int> indexes,
                        const int &ptNumber,
                        const int &dimension,
                        KDTree<dataType> *parent,
                        KDTreeMap &correspondance_map,
                        std::vector<std::vector<dataType>> &weights,
                        const int weight_number = 1);
    void buildRecursive(dataType *coordinates,
                        std::vector<int> indexes,
                        const int &ptNumber,
                        const int &dimension,
                        KDTree<dataType> *parent,
                        KDTreeMap &correspondance_map,
                        const int weight_number = 1);

    void updateWeight(dataType new_weight, const int weight_index = 0);
    void updateMinSubweight(const int weight_index = 0);
    void getKClosest(const unsigned int k,
                     const std::vector<dataType> &coordinates,
                     KDTreeMap &neighbours,
                     std::vector<dataType> &costs,
                     const int weight_index = 0);
    void recursiveGetKClosest(const unsigned int k,
                              const std::vector<dataType> &coordinates,
                              KDTreeMap &neighbours,
                              std::vector<dataType> &costs,
                              const int weight_index = 0);

    dataType cost(const std::vector<dataType> &coordinates);
    dataType distanceToBox(const KDTree<dataType> &subtree,
                           const std::vector<dataType> &coordinates);
    std::vector<dataType> getCoordinates();
    dataType getWeight(const int weight_index = 0);
    dataType getMinSubWeight(const int weight_index = 0);

    bool isLeaf();
    bool isRoot();

    template <typename type>
    inline static type abs(const type var) {
      return (var > 0) ? var : -var;
    }
  };

  template <typename dataType>
  std::vector<dataType> KDTree<dataType>::getCoordinates() {
    return coordinates_;
  }

  template <typename dataType>
  dataType KDTree<dataType>::getWeight(const int weight_index) {
    return weight_[weight_index];
  }

  template <typename dataType>
  dataType KDTree<dataType>::getMinSubWeight(const int weight_index) {
    return min_subweights_[weight_index];
  }

  template <typename dataType>
  typename KDTree<dataType>::KDTreeMap
    KDTree<dataType>::build(dataType *data,
                            const int &ptNumber,
                            const int &dimension,
                            const int weight_number) {
    KDTreeMap correspondance_map(ptNumber);
    // First, perform a argsort on the data
    // initialize original index locations
    for(int axis = 0; axis < dimension; axis++) {
      coords_min_.push_back(std::numeric_limits<dataType>::lowest());
      coords_max_.push_back(std::numeric_limits<dataType>::max());
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
      coordinates_.push_back(data[dimension * median_idx + axis]);
    }
    for(int w = 0; w < weight_number; w++) {
      weight_.push_back(0);
      min_subweights_.push_back(0);
    }
    id_ = median_idx;
    parent_ = nullptr;
    level_ = 0;

    if(idx.size() > 2) {
      // Build left leaf
      std::vector<int> idx_left(median_loc);
      for(int i = 0; i < median_loc; i++) {
        idx_left[i] = idx[i];
      }

      this->left_ = std::unique_ptr<KDTree>(
        new KDTree(this, (coords_number_ + 1) % dimension, true));
      this->left_->buildRecursive(data, idx_left, ptNumber, dimension, this,
                                  correspondance_map, weight_number);
    }

    if(idx.size() > 1) {
      // Build right leaf
      std::vector<int> idx_right(ptNumber - median_loc - 1);
      for(int i = 0; i < ptNumber - median_loc - 1; i++) {
        idx_right[i] = idx[i + median_loc + 1];
      }
      this->right_ = std::unique_ptr<KDTree>(
        new KDTree(this, (coords_number_ + 1) % dimension, false));
      this->right_->buildRecursive(data, idx_right, ptNumber, dimension, this,
                                   correspondance_map, weight_number);
    }

    return correspondance_map;
  }

  template <typename dataType>
  void KDTree<dataType>::buildRecursive(dataType *data,
                                        std::vector<int> idx_side,
                                        const int &ptNumber,
                                        const int &dimension,
                                        KDTree<dataType> *parent,
                                        KDTreeMap &correspondance_map,
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
      coordinates_.push_back(data[dimension * median_idx + axis]);
    }
    for(int w = 0; w < weight_number; w++) {
      weight_.push_back(0);
      min_subweights_.push_back(0);
    }
    id_ = median_idx;
    parent_ = parent;
    level_ = parent->level_ + 1;

    // Create bounding box
    for(int axis = 0; axis < dimension; axis++) {
      coords_min_.push_back(parent_->coords_min_[axis]);
      coords_max_.push_back(parent_->coords_max_[axis]);
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

      this->left_ = std::unique_ptr<KDTree>(
        new KDTree(this, (coords_number_ + 1) % dimension, true));
      this->left_->buildRecursive(data, idx_left, ptNumber, dimension, this,
                                  correspondance_map, weight_number);
    }

    if(idx_side.size() > 1) {
      // Build right leaf
      std::vector<int> idx_right(idx_side.size() - median_loc - 1);
      for(unsigned int i = 0; i < idx_side.size() - median_loc - 1; i++) {
        idx_right[i] = idx_side[i + median_loc + 1];
      }
      this->right_ = std::unique_ptr<KDTree>(
        new KDTree(this, (coords_number_ + 1) % dimension, false));
      this->right_->buildRecursive(data, idx_right, ptNumber, dimension, this,
                                   correspondance_map, weight_number);
    }
    return;
  }

  template <typename dataType>
  void KDTree<dataType>::updateWeight(dataType new_weight,
                                      const int weight_index) {
    weight_[weight_index] = new_weight;
    updateMinSubweight(weight_index);
  }

  template <typename dataType>
  void KDTree<dataType>::updateMinSubweight(const int weight_index) {
    dataType new_min_subweight;
    if(this->isLeaf()) {
      new_min_subweight = weight_[weight_index];
    } else if(!left_) {
      new_min_subweight = std::min(
        right_->min_subweights_[weight_index], weight_[weight_index]);
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

  template <typename dataType>
  void KDTree<dataType>::getKClosest(const unsigned int k,
                                     const std::vector<dataType> &coordinates,
                                     KDTreeMap &neighbours,
                                     std::vector<dataType> &costs,
                                     const int weight_index) {
    /// Puts the k closest points to the given coordinates in the "neighbours"
    /// vector along with their costs in the "costs" vector The output is not
    /// sorted, if you are interested in the k nearest neighbours in the order,
    /// will need to sort them according to their cost.
    if(this->isLeaf()) {
      dataType cost = this->cost(coordinates);
      cost += weight_[weight_index];
      neighbours.push_back(this);
      costs.push_back(cost);
    } else {
      this->recursiveGetKClosest(
        k, coordinates, neighbours, costs, weight_index);
    }
    // TODO sort neighbours and costs !
    return;
  }

  template <typename dataType>
  void KDTree<dataType>::recursiveGetKClosest(
    const unsigned int k,
    const std::vector<dataType> &coordinates,
    KDTreeMap &neighbours,
    std::vector<dataType> &costs,
    const int weight_index) {
    // 1- Look wether or not to include the current point in the nearest
    // neighbours
    dataType cost = this->cost(coordinates);
    cost += weight_[weight_index];

    if(costs.size() < k) {
      neighbours.push_back(this);
      costs.push_back(cost);
    } else {
      // 1.1- Find the most costly amongst neighbours
      std::vector<int> idx(k);
      for(unsigned int i = 0; i < k; i++) {
        idx[i] = i;
      }
      int idx_max_cost = *std::max_element(
        idx.begin(), idx.end(),
        [&costs](int &a, int &b) { return costs[a] < costs[b]; });
      dataType max_cost = costs[idx_max_cost];

      // 1.2- If the current KDTree is less costly, put it in the neighbours and
      // update costs.
      if(cost < max_cost) {
        costs[idx_max_cost] = cost;
        neighbours[idx_max_cost] = this;
      }
    }

    // 2- Recursively visit KDTrees that are worth it
    if(left_) {
      dataType max_cost = *std::max_element(costs.begin(), costs.end());
      dataType &min_subweight = left_->min_subweights_[weight_index];
      dataType d_min = this->distanceToBox(*left_, coordinates);
      if(costs.size() < k || d_min + min_subweight < max_cost) {
        // 2.2- It is possible that there exists a point in this subtree that is
        // less costly than max_cost
        left_->recursiveGetKClosest(
          k, coordinates, neighbours, costs, weight_index);
      }
    }

    if(right_) {
      dataType max_cost = *std::max_element(costs.begin(), costs.end());
      dataType &min_subweight = right_->min_subweights_[weight_index];
      dataType d_min = this->distanceToBox(*right_, coordinates);
      if(costs.size() < k || d_min + min_subweight < max_cost) {
        // 2.2- It is possible that there exists a point in this subtree that is
        // less costly than max_cost
        right_->recursiveGetKClosest(
          k, coordinates, neighbours, costs, weight_index);
      }
    }
    return;
  }

  template <typename dataType>
  dataType KDTree<dataType>::cost(const std::vector<dataType> &coordinates) {
    dataType cost = 0;
    for(size_t i = 0; i < coordinates.size(); i++) {
      cost += Geometry::pow(abs(coordinates[i] - coordinates_[i]), p_);
    }
    return cost;
  }

  template <typename dataType>
  dataType
    KDTree<dataType>::distanceToBox(const KDTree<dataType> &subtree,
                                    const std::vector<dataType> &coordinates) {
    dataType d_min = 0;
    for(size_t axis = 0; axis < coordinates.size(); axis++) {
      if(subtree.coords_min_[axis] > coordinates[axis]) {
        d_min
          += Geometry::pow(subtree.coords_min_[axis] - coordinates[axis], p_);
      } else if(subtree.coords_max_[axis] < coordinates[axis]) {
        d_min
          += Geometry::pow(coordinates[axis] - subtree.coords_max_[axis], p_);
      }
    }
    return d_min;
  }

  template <typename dataType>
  bool KDTree<dataType>::isLeaf() {
    return left_ == nullptr && right_ == nullptr;
  }

  template <typename dataType>
  bool KDTree<dataType>::isRoot() {
    return parent_ == nullptr;
  }
} // namespace ttk
#include <buildWeights.h>

#endif
