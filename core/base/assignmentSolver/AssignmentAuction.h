/// \ingroup base
/// \class ttk::AssignmentAuction
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
///
/// Auction algorithm for Balanced and Unbalanced Assignement Problem
///
/// For the unbalanced problem:
///   The cost matrix in input has a size of (n + 1) x (m + 1)
///   - n is the number of jobs, m the number of workers
///   - the nth row contains the cost of not assigning workers
///   - the mth column is the same but with jobs
///   - the last cell (costMatrix[n][m]) is not used

#ifndef _ASSIGNMENTAUCTION_H
#define _ASSIGNMENTAUCTION_H

#include <Debug.h>

#include "AssignmentSolver.h"

namespace ttk {

  template <class dataType>
  class AssignmentAuction : virtual public Debug,
                            public AssignmentSolver<dataType> {

  public:
    AssignmentAuction(){};

    ~AssignmentAuction(){};

    int run(std::vector<asgnMatchingTuple> &matchings);
    void runAuctionRound(std::vector<std::vector<dataType>> &cMatrix);

    void initFirstRound();
    void initBiddersAndGoods();
    void initEpsilon();
    void epsilonScaling();
    void makeBalancedMatrix(std::vector<std::vector<dataType>> &matrix);

    bool stoppingCriterion(std::vector<std::vector<dataType>> &cMatrix);
    dataType getRelativePrecision(std::vector<std::vector<dataType>> &cMatrix);
    dataType getMatchingDistance(std::vector<std::vector<dataType>> &cMatrix);

    inline void setBalanced(bool balanced) {
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

    template <class vecType>
    void printVector(std::vector<vecType> &vec) {
      std::stringstream ss;
      for(auto valTemp : vec)
        ss << valTemp << " ";
      printMsg(ss.str(), debug::Priority::VERBOSE);
      printMsg(debug::Separator::L1, debug::Priority::VERBOSE);
    }

  private:
    int numberOfRounds = -1;
    int iter = 0;
    double epsilon = -1;
    double epsilonDiviserMultiplier = 0;
    dataType delta_lim = 0.01;

    std::vector<int> bidderAssignments{1, -1};
    std::vector<int> goodAssignments{};
    std::vector<double> goodPrices{};
  }; // AssignmentAuction Class

  template <typename type>
  static type abs(const type var) {
    return (var >= 0) ? var : -var;
  }

  template <class dataType>
  dataType getMaxValue(std::vector<std::vector<dataType>> &matrix,
                       bool balancedAsgn) {
    unsigned int nRows = matrix.size();
    unsigned int nCols = matrix[0].size();
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

  template <class dataType>
  dataType getSecondMinValueVector(std::vector<dataType> &vec) {
    dataType secondMin = std::numeric_limits<dataType>::max();
    dataType firstMin = std::numeric_limits<dataType>::max();
    for(auto elem : vec) {
      if(elem < firstMin) {
        secondMin = firstMin;
        firstMin = elem;
      } else if(elem < secondMin)
        secondMin = elem;
    }
    return secondMin;
  }

  template <typename dataType>
  void AssignmentAuction<dataType>::initEpsilon() {
    if(epsilon == -1) {
      dataType maxValue
        = getMaxValue<dataType>(this->costMatrix, this->balancedAssignment);
      // int tRowSize = this->balancedAssignment ? this->rowSize :
      // (this->rowSize-1)+(this->colSize-1); int tColSize =
      // this->balancedAssignment ? this->colSize :
      // (this->rowSize-1)+(this->colSize-1); epsilon = maxValue *
      // std::min(tRowSize, tColSize)/2;
      epsilon = maxValue / 4;
      // epsilon = std::pow(maxValue, 2)/4;
      // epsilon += *std::max_element(goodPrices.begin(), goodPrices.end());
      // epsilon += getSecondMinValueVector(goodPrices);
      if(epsilon == 0)
        epsilon = 1;
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
    /*goodPrices.clear();
    goodPrices.resize(this->colSize, 0);*/
  }

  template <typename dataType>
  void AssignmentAuction<dataType>::initFirstRound() {
    iter = 0;
    bidderAssignments[0] = -1;
    // epsilon /= ((epsilonDiviserMultiplier==0) ? 1 : epsilonDiviserMultiplier
    // * 5);
  }

  template <typename dataType>
  void AssignmentAuction<dataType>::makeBalancedMatrix(
    std::vector<std::vector<dataType>> &matrix) {
    unsigned int nRows = matrix.size();
    unsigned int nCols = matrix[0].size();
    matrix[nRows - 1][nCols - 1] = 0;

    // Add rows
    for(unsigned int i = 0; i < nCols - 2; ++i) {
      std::vector<dataType> newLine(matrix[nRows - 1]);
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
      int bidderId = unassignedBidders.front();
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

      // Update price
      double delta = abs<dataType>(bestValue - bestSecondValue) + epsilon;
      goodPrices[bestGoodId] = goodPrices[bestGoodId] + delta;

      // printVector(goodPrices);
    }
  }

  template <typename dataType>
  int AssignmentAuction<dataType>::run(
    std::vector<asgnMatchingTuple> &matchings) {
    initEpsilon();

    // Try to avoid price war
    double tempPrice = *std::max_element(goodPrices.begin(), goodPrices.end());
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
      makeBalancedMatrix(this->costMatrix);
    // printTableVector(this->costMatrix);

    // Run auction
    initFirstRound();
    // printVector(goodPrices);
    while(not stoppingCriterion(this->costMatrix)) {
      initBiddersAndGoods();
      runAuctionRound(this->costMatrix);
      // std::cout << epsilon << std::endl;
      // printVector(goodPrices);
      // printVector(bidderAssignments);
      epsilonScaling();
      iter++;
      if(numberOfRounds != -1 and iter >= numberOfRounds)
        break;
    }
    // printVector(goodPrices);

    // Create output matching
    for(unsigned int bidderId = 0; bidderId < bidderAssignments.size();
        ++bidderId) {
      /*int i = std::min(bidderId, this->rowSize-1);
      int j = std::min(bidderAssignments[bidderId], this->colSize-1);*/
      int i = bidderId;
      int j = bidderAssignments[bidderId];
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
    // std::cout << "delta = " << delta << std::endl;
    return not(delta > delta_lim);
  }

  // Adapted from Persistence Diagrams Auction
  template <typename dataType>
  dataType AssignmentAuction<dataType>::getRelativePrecision(
    std::vector<std::vector<dataType>> &cMatrix) {
    dataType d = this->getMatchingDistance(cMatrix);
    if(d < 1e-12) {
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
      int i = bidderId;
      int j = bidderAssignments[bidderId];
      d += cMatrix[i][j];
    }
    return d;
  }

} // namespace ttk

#endif
