/// \ingroup base
/// \class ttk::AssignmentExhaustive
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
///
/// Exhaustive Search for Unbalanced Assignment Problem
///   (TODO manage balanced problem)
///
/// The cost matrix in input has a size of (n + 1) x (m + 1)
/// - n is the number of jobs, m the number of workers
/// - the nth row contains the cost of not assigning workers
/// - the mth column is the same but with jobs
/// - the last cell (costMatrix[n][m]) is not used
///
/// An exhaustive search for the assignment problem has an exponential
/// complexity Use this algorithm only for small assignment problem

#pragma once

#include <AssignmentSolver.h>
#include <Debug.h>

#include <iterator>
#include <limits>
#include <map>
#include <queue>
#include <set>

namespace ttk {

  template <class dataType>
  class AssignmentExhaustive : virtual public Debug,
                               public AssignmentSolver<dataType> {

  public:
    AssignmentExhaustive() = default;

    ~AssignmentExhaustive() override = default;

    int run(std::vector<MatchingType> &matchings) override;

    void setSaveAsgn(bool doSave) {
      saveAsgn = doSave;
    }

    dataType tryAssignment(std::vector<int> &asgn,
                           std::vector<MatchingType> &matchings);

    std::string vectorToString(std::vector<dataType> &my_vector) {
      std::stringstream result;
      std::copy(my_vector.begin(), my_vector.end(),
                std::ostream_iterator<dataType>(result, ""));
      return result.str();
    }

    void enumerateAssignments(unsigned int min_dim,
                              unsigned int max_dim,
                              std::vector<std::vector<int>> &allAsgn) {
      std::queue<std::tuple<std::vector<int>, std::vector<int>, bool,
                            std::vector<bool>, int>>
        queue;

      std::vector<int> asgn_temp, unasgn_temp;
      std::vector<bool> const done_temp(max_dim, false);
      queue.emplace(asgn_temp, unasgn_temp, false, done_temp, -1);

      while(!queue.empty()) {
        auto queueElem = queue.front();
        queue.pop();
        std::vector<int> const asgn = std::get<0>(queueElem);
        std::vector<int> const unasgn = std::get<1>(queueElem);
        bool const toUnasgn = std::get<2>(queueElem);
        std::vector<bool> done = std::get<3>(queueElem);
        unsigned int const maxDone = std::get<4>(queueElem);

        unsigned int const ind = asgn.size();
        for(unsigned int j = 0; j < max_dim + 1; ++j) {
          std::vector<int> new_asgn(asgn);
          std::vector<int> new_unasgn(unasgn);
          std::vector<bool> new_done(done);
          int new_maxDone = maxDone;
          if(j < max_dim) {
            if(not done[j]) {
              if(ind >= min_dim and maxDone > j)
                continue;
            } else
              continue;

            if(not toUnasgn) {
              new_asgn.push_back(j);
              if(ind >= min_dim)
                new_maxDone = std::max(maxDone, j);
            } else {
              new_unasgn.push_back(j);
              new_maxDone = std::max(maxDone, j);
            }
            new_done[j] = true;

            unsigned int const new_ind = new_asgn.size();
            if(new_ind == max_dim) {
              // A new assignment is made here
              for(auto &new_unasgn_elem : new_unasgn)
                new_asgn.push_back(new_unasgn_elem);
              // std::sort(new_asgn.begin()+min_dim, new_asgn.end());
              allAsgn.push_back(new_asgn);
            } else
              queue.emplace(new_asgn, new_unasgn, false, new_done, new_maxDone);

          } else {
            if(not toUnasgn and ind < min_dim) {
              new_asgn.push_back(max_dim);
              queue.emplace(new_asgn, new_unasgn, true, new_done, new_maxDone);
            }
          }
        }
      }
    }

    void printAssignments(std::vector<std::vector<int>> &allAsgn) {
      printMsg(debug::Separator::L1, debug::Priority::VERBOSE);
      std::stringstream ss;
      ss << "{";
      for(auto &vecTemp : allAsgn) {
        ss << "{";
        for(unsigned int i = 0; i < vecTemp.size(); ++i) {
          auto valTemp = vecTemp[i];
          ss << valTemp;
          if(i != vecTemp.size() - 1)
            ss << ",";
        }
        ss << "}, ";
      }
      ss << "}";
      printMsg(ss.str(), debug::Priority::VERBOSE);
      printMsg(debug::Separator::L1, debug::Priority::VERBOSE);
    }

  private:
    // This attribute will store the computed assignments
    std::map<std::string, std::vector<std::vector<int>>> savedAsgn;
    bool saveAsgn = false;
  };

  template <typename dataType>
  dataType AssignmentExhaustive<dataType>::tryAssignment(
    std::vector<int> &asgn, std::vector<MatchingType> &matchings) {
    unsigned int const nRows = this->costMatrix.size() - 1;
    unsigned int const nCols = this->costMatrix[0].size() - 1;
    // int max_dim = std::max(nRows, nCols);
    unsigned int const min_dim = std::min(nRows, nCols);
    bool const transpose = nRows > nCols;

    dataType cost = 0;
    for(unsigned int ind = 0; ind < asgn.size(); ++ind) {
      int const indMatrix = std::min(ind, min_dim);
      int const i = (!transpose) ? indMatrix : asgn[ind];
      int const j = (!transpose) ? asgn[ind] : indMatrix;
      cost += this->costMatrix[i][j];
      matchings.push_back(std::make_tuple(i, j, this->costMatrix[i][j]));
    }
    return cost;
  }

  template <typename dataType>
  int AssignmentExhaustive<dataType>::run(
    std::vector<MatchingType> &matchings) {
    int const nRows = this->costMatrix.size() - 1;
    int const nCols = this->costMatrix[0].size() - 1;
    int const max_dim = std::max(nRows, nCols);
    int const min_dim = std::min(nRows, nCols);

    // --- Construct all possible assignments
    std::vector<std::vector<int>> allAsgn;
    // We hard write the basic assignments to avoid the call of
    // enumerateAssignments
    // These assignments are automatically generated by enumerateAssignments
    if(min_dim == 1 and max_dim == 1)
      allAsgn = {{0}, {1, 0}};
    else if(min_dim == 1 and max_dim == 2)
      allAsgn = {{0, 1}, {2, 0, 1}, {1, 0}};
    else if(min_dim == 1 and max_dim == 3)
      allAsgn = {{0, 1, 2}, {3, 0, 1, 2}, {1, 0, 2}, {2, 0, 1}};
    else if(min_dim == 1 and max_dim == 4)
      allAsgn = {{0, 1, 2, 3},
                 {4, 0, 1, 2, 3},
                 {1, 0, 2, 3},
                 {2, 0, 1, 3},
                 {3, 0, 1, 2}};
    else if(min_dim == 1 and max_dim == 5)
      allAsgn = {{0, 1, 2, 3, 4}, {5, 0, 1, 2, 3, 4}, {1, 0, 2, 3, 4},
                 {2, 0, 1, 3, 4}, {3, 0, 1, 2, 4},    {4, 0, 1, 2, 3}};
    else if(min_dim == 1 and max_dim == 6)
      allAsgn = {{0, 1, 2, 3, 4, 5}, {6, 0, 1, 2, 3, 4, 5}, {1, 0, 2, 3, 4, 5},
                 {2, 0, 1, 3, 4, 5}, {3, 0, 1, 2, 4, 5},    {4, 0, 1, 2, 3, 5},
                 {5, 0, 1, 2, 3, 4}};
    else if(min_dim == 2 and max_dim == 2)
      allAsgn = {{0, 1}, {0, 2, 1}, {2, 1, 0}, {2, 2, 0, 1},
                 {1, 0}, {1, 2, 0}, {2, 0, 1}};
    else {
      std::stringstream ss;
      ss << min_dim << "_" << max_dim;
      std::string const asgnName = ss.str();
      // not found
      if(not saveAsgn or savedAsgn.find(asgnName) == savedAsgn.end()) {
        if(saveAsgn)
          printMsg(asgnName, debug::Priority::VERBOSE);
        enumerateAssignments(min_dim, max_dim, allAsgn);
        if(saveAsgn) {
          savedAsgn[asgnName] = allAsgn;
          std::stringstream ss2;
          ss2 << allAsgn.size() << " done";
          printMsg(ss2.str(), debug::Priority::VERBOSE);
        }
      } else
        allAsgn = savedAsgn[asgnName];
    }

    // --- Try these assignments and get the better
    dataType bestCost = std::numeric_limits<dataType>::max();
    std::vector<MatchingType> bestMatching;
    for(std::vector<int> &asgn : allAsgn) {
      std::vector<MatchingType> tempMatching;
      dataType cost = tryAssignment(asgn, tempMatching);
      if(bestCost > cost) {
        bestCost = cost;
        bestMatching = tempMatching;
      }
    }
    matchings = bestMatching;

    return 0;
  }

} // namespace ttk
