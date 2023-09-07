/// \ingroup base
/// \class ttk::AssignmentAuction
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
///
/// Auction algorithm for Balanced and Unbalanced Assignment Problem
///
/// For the unbalanced problem:
///   The cost matrix in input has a size of (n + 1) x (m + 1)
///   - n is the number of jobs, m the number of workers
///   - the nth row contains the cost of not assigning workers
///   - the mth column is the same but with jobs
///   - the last cell (costMatrix[n][m]) is not used

#pragma once

#include <AssignmentSolver.h>

#include <limits>
#include <queue>

namespace ttk {

  template <class dataType>
  class AssignmentAuction : virtual public Debug,
                            public AssignmentSolver<dataType> {

  public:
    AssignmentAuction() = default;

    ~AssignmentAuction() override = default;

    int run(std::vector<MatchingType> &matchings) override;
    void runAuctionRound(std::vector<std::vector<dataType>> &cMatrix);

    void initFirstRound();
    void initBiddersAndGoods();
    void initEpsilon();
    void epsilonScaling();
    void makeBalancedMatrix(std::vector<std::vector<dataType>> &matrix);

    bool stoppingCriterion(std::vector<std::vector<dataType>> &cMatrix);
    dataType getRelativePrecision(std::vector<std::vector<dataType>> &cMatrix);
    dataType getMatchingDistance(std::vector<std::vector<dataType>> &cMatrix);

    inline void setBalanced(bool balanced) override {
      AssignmentSolver<dataType>::setBalanced(balanced);
      if(this->balancedAssignment)
        goodPrices.resize(this->colSize, 0);
      else
        goodPrices.resize((this->colSize - 1) + (this->rowSize - 1), 0);
    }

    inline void setNumberOfRounds(int noRounds) {
      numberOfRounds = noRounds;
    }

    inline int getIter() {
      return iter;
    }

    inline void setEpsilon(double eps) {
      epsilon = eps;
    }

    inline void setEpsilonDiviserMultiplier(double div) {
      epsilonDiviserMultiplier = div;
    }

    inline void setPrices(std::vector<double> prices) {
      goodPrices.resize(prices.size(), 0);
      for(unsigned int i = 0; i < prices.size(); ++i)
        goodPrices[i] = prices[i];
    }

    inline std::vector<double> &getPrices() {
      return goodPrices;
    }

  private:
    int numberOfRounds = -1;
    int iter = 0;
    double epsilon = -1;
    double epsilonDiviserMultiplier = 0;
    double delta_lim = 0.01;

    dataType lowerBoundCostWeight = 1 + delta_lim;

    // Filled
    std::vector<int> bidderAssignments{1, -1}, bestBidderAssignments;
    std::vector<int> goodAssignments{};
    std::vector<double> goodPrices{};
    dataType lowerBoundCost;
  }; // AssignmentAuction Class

  template <typename type>
  static type abs(const type var) {
    return (var >= 0) ? var : -var;
  }

  template <class dataType>
  dataType getMaxValue(std::vector<std::vector<dataType>> &matrix,
                       bool balancedAsgn) {
    unsigned int const nRows = matrix.size();
    unsigned int const nCols = matrix[0].size();
    dataType maxValue = std::numeric_limits<dataType>::lowest();
    for(unsigned int i = 0; i < nRows; ++i) {
      for(unsigned int j = 0; j < nCols; ++j) {
        if(not balancedAsgn and i == nRows - 1 and j == nCols - 1)
          continue;
        maxValue = (matrix[i][j] > maxValue) ? matrix[i][j] : maxValue;
      }
    }
    return maxValue;
  }

  template <typename dataType>
  void AssignmentAuction<dataType>::initEpsilon() {
    if(epsilon == -1.0) {
      dataType maxValue
        = getMaxValue<dataType>(this->costMatrix, this->balancedAssignment);
      epsilon = maxValue / 4.0;
      if(epsilon == 0.0)
        epsilon = 1.0;
      epsilon
        /= ((epsilonDiviserMultiplier == 0) ? 1 : epsilonDiviserMultiplier * 5);
    }
  }

  template <typename dataType>
  void AssignmentAuction<dataType>::epsilonScaling() {
    epsilon /= 5;
    /*if(epsilon < 1e-6)
      epsilon = 1e-6;*/
  }

  template <typename dataType>
  void AssignmentAuction<dataType>::initBiddersAndGoods() {
    bidderAssignments.clear();
    goodAssignments.clear();
    if(this->balancedAssignment) {
      bidderAssignments.resize(this->rowSize, -1);
      goodAssignments.resize(this->colSize, -1);
    } else {
      bidderAssignments.resize((this->rowSize - 1) + (this->colSize - 1), -1);
      goodAssignments.resize((this->colSize - 1) + (this->rowSize - 1), -1);
    }
  }

  template <typename dataType>
  void AssignmentAuction<dataType>::initFirstRound() {
    iter = 0;
    bidderAssignments[0] = -1;
  }

  template <typename dataType>
  void AssignmentAuction<dataType>::makeBalancedMatrix(
    std::vector<std::vector<dataType>> &matrix) {
    unsigned int const nRows = matrix.size();
    unsigned int const nCols = matrix[0].size();
    matrix[nRows - 1][nCols - 1] = 0;

    // Add rows
    for(unsigned int i = 0; i < nCols - 2; ++i) {
      std::vector<dataType> const newLine(matrix[nRows - 1]);
      matrix.push_back(newLine);
    }
    // Add columns
    for(unsigned int i = 0; i < (nRows - 1) + (nCols - 1); ++i) {
      for(unsigned int j = 0; j < nRows - 2; ++j) {
        matrix[i].push_back(matrix[i][nCols - 1]);
      }
    }
  }

  // ----------------------------------------
  // Main Functions
  // ----------------------------------------
  template <typename dataType>
  void AssignmentAuction<dataType>::runAuctionRound(
    std::vector<std::vector<dataType>> &cMatrix) {
    std::queue<int> unassignedBidders{};
    for(unsigned int i = 0; i < bidderAssignments.size(); ++i)
      unassignedBidders.push(i);

    while(!unassignedBidders.empty()) {
      int const bidderId = unassignedBidders.front();
      unassignedBidders.pop();

      // Get good with highest value
      dataType bestValue = std::numeric_limits<dataType>::lowest();
      dataType bestSecondValue = std::numeric_limits<dataType>::lowest();
      int bestGoodId = -1;
      for(unsigned int goodId = 0; goodId < goodPrices.size(); ++goodId) {
        if(cMatrix[bidderId][goodId] == static_cast<dataType>(-1))
          continue;
        dataType goodPrice = goodPrices[goodId];
        dataType value = -cMatrix[bidderId][goodId] - goodPrice;
        if(value > bestValue) {
          bestSecondValue = bestValue;
          bestValue = value;
          bestGoodId = goodId;
        } else if(value > bestSecondValue)
          bestSecondValue = value;
      }

      // Update assignments
      bidderAssignments[bidderId] = bestGoodId;
      if(goodAssignments[bestGoodId] != -1)
        unassignedBidders.push(goodAssignments[bestGoodId]);
      goodAssignments[bestGoodId] = bidderId;

      // If there is only one acceptable good for the bidder
      if(bestSecondValue == std::numeric_limits<dataType>::lowest())
        bestSecondValue = bestValue;

      // Update price
      double const delta = abs<dataType>(bestValue - bestSecondValue) + epsilon;
      double newPrice = goodPrices[bestGoodId] + delta;
      if(newPrice > std::numeric_limits<double>::max() / 2) {
        // Avoid price explosion
        newPrice = goodPrices[bestGoodId] + epsilon;
      }
      goodPrices[bestGoodId] = newPrice;
    }
  }

  template <typename dataType>
  dataType getLowerBoundCost(std::vector<std::vector<dataType>> &costMatrix) {
    std::vector<dataType> minCol(
      costMatrix[0].size(), std::numeric_limits<dataType>::max()),
      minRow(costMatrix.size(), std::numeric_limits<dataType>::max());
    for(unsigned int i = 0; i < costMatrix.size(); ++i) {
      for(unsigned int j = 0; j < costMatrix[i].size(); ++j) {
        if(costMatrix[i][j] < minCol[j])
          minCol[j] = costMatrix[i][j];
        if(costMatrix[i][j] < minRow[i])
          minRow[i] = costMatrix[i][j];
      }
    }
    dataType minColCost = 0, minRowCost = 0;
    for(unsigned int i = 0; i < minCol.size(); ++i)
      minColCost += minCol[i];
    for(unsigned int i = 0; i < minRow.size(); ++i)
      minRowCost += minRow[i];
    dataType lowerBoundCost = std::max(minColCost, minRowCost);

    return lowerBoundCost;
  }

  template <typename dataType>
  int AssignmentAuction<dataType>::run(std::vector<MatchingType> &matchings) {
    initEpsilon();
    dataType bestCost = std::numeric_limits<dataType>::max();

    // Try to avoid price war
    double const tempPrice
      = *std::max_element(goodPrices.begin(), goodPrices.end());
    std::vector<double> savedPrices;
    for(unsigned int i = 0; i < goodPrices.size(); ++i) {
      auto old = goodPrices[i];
      goodPrices[i]
        = goodPrices[i] * epsilon / ((tempPrice == 0) ? 1 : tempPrice);
      auto t = old - goodPrices[i];
      savedPrices.push_back(t);
    }

    // Make balanced cost matrix
    if(not this->balancedAssignment)
      this->makeBalancedMatrix(this->costMatrix);

    // Get lower bound cost
    lowerBoundCost = getLowerBoundCost(this->costMatrix);

    // Run auction
    initFirstRound();
    while(not stoppingCriterion(this->costMatrix)) {
      initBiddersAndGoods();
      runAuctionRound(this->costMatrix);

      dataType cost = getMatchingDistance(this->costMatrix);
      if(cost < bestCost) {
        bestCost = cost;
        bestBidderAssignments = bidderAssignments;
      }

      epsilonScaling();
      iter++;

      if(numberOfRounds != -1 and iter >= numberOfRounds)
        break;
    }
    bidderAssignments = bestBidderAssignments;

    // Create output matching
    for(unsigned int bidderId = 0; bidderId < bidderAssignments.size();
        ++bidderId) {
      int const i = bidderId;
      int const j = bidderAssignments[bidderId];
      if(this->balancedAssignment
         or (not this->balancedAssignment
             and not(i >= this->rowSize - 1 and j >= this->colSize - 1))) {
        matchings.push_back(std::make_tuple(i, j, this->costMatrix[i][j]));
      }
    }

    // Set prices as before
    for(unsigned int i = 0; i < goodPrices.size(); ++i)
      goodPrices[i] += savedPrices[i];

    return 0;
  }

  // ----------------------------------------
  // Stopping Criterion Functions
  // ----------------------------------------
  // Adapted from Persistence Diagrams Auction
  template <typename dataType>
  bool AssignmentAuction<dataType>::stoppingCriterion(
    std::vector<std::vector<dataType>> &cMatrix) {
    if(bidderAssignments[0] == -1) // Auction not started
      return false;
    dataType delta = 5;
    delta = getRelativePrecision(cMatrix);
    return not(delta > delta_lim);
  }

  // Adapted from Persistence Diagrams Auction
  template <typename dataType>
  dataType AssignmentAuction<dataType>::getRelativePrecision(
    std::vector<std::vector<dataType>> &cMatrix) {
    dataType d = this->getMatchingDistance(cMatrix);
    if(d < 1e-6 or d <= (lowerBoundCost * lowerBoundCostWeight)) {
      return 0;
    }
    dataType denominator = d - bidderAssignments.size() * epsilon;
    if(denominator <= 0) {
      return 1;
    } else {
      return d / denominator - 1;
    }
  }

  // Adapted from Persistence Diagrams Auction
  template <typename dataType>
  dataType AssignmentAuction<dataType>::getMatchingDistance(
    std::vector<std::vector<dataType>> &cMatrix) {
    dataType d = 0;
    for(unsigned int bidderId = 0; bidderId < bidderAssignments.size();
        ++bidderId) {
      int const i = bidderId;
      int const j = bidderAssignments[bidderId];
      d += cMatrix[i][j];
    }
    return d;
  }

} // namespace ttk
