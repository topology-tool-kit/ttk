/// \ingroup base
/// \class ttk::AssignmentExhaustive
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
///
/// Exhaustive Search for Unbalanced Assignement Problem
///   (TODO manage balanced problem)
///
/// The cost matrix in input has a size of (n + 1) x (m + 1)
/// - n is the number of jobs, m the number of workers
/// - the nth row contains the cost of not assigning workers
/// - the mth column is the same but with jobs
/// - the last cell (costMatrix[n][m]) is not used
///
/// An exhaustive search for the assignment problem has an exponential
/// complexity Use this algorithm only for small assigment problem

#ifndef _ASSIGNMENTEXHAUSTIVE_H
#define _ASSIGNMENTEXHAUSTIVE_H

#include <Debug.h>

#include "AssignmentSolver.h"

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

    ~AssignmentExhaustive() = default;

    int run(std::vector<asgnMatchingTuple> &matchings);

    dataType tryAssignment(std::vector<int> &asgn,
                           std::vector<asgnMatchingTuple> &matchings);

    std::string vectorToString(std::vector<dataType> &my_vector) {
      std::stringstream result;
      std::copy(my_vector.begin(), my_vector.end(),
                std::ostream_iterator<dataType>(result, ""));
      return result.str();
    }

    // std::vector<std::vector<int>> constructAssignments(int min_dim, int
    // max_dim); std::vector<std::vector<int>> constructAssignments2(int
    // min_dim, int max_dim); std::vector<std::vector<int>>
    // constructAssignments3(int min_dim, int max_dim);

    // This version creates a lot of duplicates (worst version)
    // std::vector<std::vector<int>>
    // AssignmentExhaustive::constructAssignments(int min_dim, int max_dim){
    std::vector<std::vector<int>> constructAssignments(unsigned int min_dim,
                                                       unsigned int max_dim) {
      std::vector<std::vector<int>> allAsgn;
      std::queue<std::vector<int>> queue;
      std::set<std::string> unasgnAdded;

      std::vector<int> asgn_temp;
      queue.push(asgn_temp);

      while(!queue.empty()) {
        std::vector<int> asgn = queue.front();
        queue.pop();
        unsigned int ind = asgn.size();
        for(unsigned int j = 0; j < max_dim; ++j) {
          std::vector<int> new_asgn(asgn);
          // TODO avoid high complexity of find in vector
          std::vector<int>::iterator it
            = std::find(new_asgn.begin(), new_asgn.end(), j);
          if(it == new_asgn.end()) { // not found
            new_asgn.push_back(j);
          } else
            continue;

          if(ind == max_dim - 1) {
            // A new assignement without unassigned costs is made here
            // allAsgn.push_back(new_asgn);

            // Construct assignments with unassigned costs
            std::queue<std::tuple<std::vector<int>, int>> unassignedQueue;
            unassignedQueue.push(std::make_tuple(new_asgn, 0));
            while(!unassignedQueue.empty()) {
              std::tuple<std::vector<int>, int> elem = unassignedQueue.front();
              unassignedQueue.pop();
              std::vector<int> elemVec = std::get<0>(elem);
              int elemInd = std::get<1>(elem);

              if(elemInd < (int)min_dim) {
                // Add unchanged assignment
                unassignedQueue.push(std::make_tuple(elemVec, elemInd + 1));

                // Add new assignment (with unassigned cost)
                std::vector<int> new_unasgn(elemVec);
                new_unasgn.push_back(new_unasgn[elemInd]);
                new_unasgn[elemInd] = max_dim;
                unassignedQueue.push(std::make_tuple(new_unasgn, elemInd + 1));
              } else {
                // A new assignment is made here
                // TODO there is probably be a better way to avoid duplicates
                std::sort(elemVec.begin() + min_dim, elemVec.end());
                std::string new_string = vectorToString(elemVec);
                auto it2 = unasgnAdded.find(new_string);
                if(it2 == unasgnAdded.end()) {
                  unasgnAdded.insert(vectorToString(elemVec));
                  allAsgn.push_back(elemVec);
                } // else
                  // std::cout << "collision " << vectorToString(elemVec) <<
                  // std::endl;
              }
            }

          } else
            queue.push(new_asgn);
        }
      }

      return allAsgn;
    }

    // std::vector<std::vector<int>>
    // AssignmentExhaustive::constructAssignments2(int min_dim, int max_dim){
    std::vector<std::vector<int>> constructAssignments2(unsigned int min_dim,
                                                        unsigned int max_dim) {
      std::vector<std::vector<int>> allAsgn;
      std::queue<std::tuple<std::vector<int>, std::vector<int>, bool>> queue;

      std::vector<int> asgn_temp, unasgn_temp;
      queue.push(std::make_tuple(asgn_temp, unasgn_temp, false));

      while(!queue.empty()) {
        auto queueElem = queue.front();
        queue.pop();
        std::vector<int> asgn = std::get<0>(queueElem);
        std::vector<int> unasgn = std::get<1>(queueElem);
        bool toUnasgn = std::get<2>(queueElem);
        unsigned int ind = asgn.size();
        for(unsigned int j = 0; j < max_dim + 1; ++j) {
          std::vector<int> new_asgn(asgn);
          std::vector<int> new_unasgn(unasgn);
          if(j < max_dim) {
            // TODO avoid high complexity of find in vector
            std::vector<int>::iterator it
              = std::find(new_asgn.begin(), new_asgn.end(), j);
            std::vector<int>::iterator it2
              = std::find(new_unasgn.begin(), new_unasgn.end(), j);
            if(it == new_asgn.end() and it2 == new_unasgn.end()) { // not found
              bool addIt = true;
              if(ind >= min_dim) {
                for(unsigned int cpt = min_dim; cpt < new_asgn.size(); ++cpt)
                  if(new_asgn[cpt] > (int)j)
                    addIt = false;
                for(unsigned int cpt = 0; cpt < new_unasgn.size(); ++cpt)
                  if(new_unasgn[cpt] > (int)j)
                    addIt = false;
              }
              if(addIt) {
                if(not toUnasgn)
                  new_asgn.push_back(j);
                else
                  new_unasgn.push_back(j);
              } else
                continue;
            } else
              continue;

            unsigned int new_ind = new_asgn.size();
            if(new_ind == max_dim) {
              // A new assignement is made here
              for(auto new_unasgn_elem : new_unasgn)
                new_asgn.push_back(new_unasgn_elem);
              // std::sort(new_asgn.begin()+min_dim, new_asgn.end());
              allAsgn.push_back(new_asgn);
            } else
              queue.push(std::make_tuple(new_asgn, new_unasgn, false));

          } else {
            if(not toUnasgn and ind < min_dim) {
              new_asgn.push_back(max_dim);
              queue.push(std::make_tuple(new_asgn, new_unasgn, true));
            }
          }
        }
      }

      return allAsgn;
    }

    // std::vector<std::vector<int>>
    // AssignmentExhaustive::constructAssignments3(int min_dim, int max_dim){
    std::vector<std::vector<int>> constructAssignments3(unsigned int min_dim,
                                                        unsigned int max_dim) {
      std::vector<std::vector<int>> allAsgn;
      std::queue<std::tuple<std::vector<int>, std::vector<int>, bool,
                            std::vector<bool>, int>>
        queue;

      std::vector<int> asgn_temp, unasgn_temp;
      std::vector<bool> done_temp(max_dim, false);
      queue.push(std::make_tuple(asgn_temp, unasgn_temp, false, done_temp, -1));

      while(!queue.empty()) {
        auto queueElem = queue.front();
        queue.pop();
        std::vector<int> asgn = std::get<0>(queueElem);
        std::vector<int> unasgn = std::get<1>(queueElem);
        bool toUnasgn = std::get<2>(queueElem);
        std::vector<bool> done = std::get<3>(queueElem);
        unsigned int maxDone = std::get<4>(queueElem);

        unsigned int ind = asgn.size();
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

            unsigned int new_ind = new_asgn.size();
            if(new_ind == max_dim) {
              // A new assignement is made here
              for(auto new_unasgn_elem : new_unasgn)
                new_asgn.push_back(new_unasgn_elem);
              // std::sort(new_asgn.begin()+min_dim, new_asgn.end());
              allAsgn.push_back(new_asgn);
            } else
              queue.push(std::make_tuple(
                new_asgn, new_unasgn, false, new_done, new_maxDone));

          } else {
            if(not toUnasgn and ind < min_dim) {
              new_asgn.push_back(max_dim);
              queue.push(std::make_tuple(
                new_asgn, new_unasgn, true, new_done, new_maxDone));
            }
          }
        }
      }

      return allAsgn;
    }

    void printAssignments(std::vector<std::vector<int>> &allAsgn) {
      printMsg(debug::Separator::L1, debug::Priority::VERBOSE);
      std::stringstream ss;
      ss << "{";
      for(auto vecTemp : allAsgn) {
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
  };

  template <typename dataType>
  dataType AssignmentExhaustive<dataType>::tryAssignment(
    std::vector<int> &asgn, std::vector<asgnMatchingTuple> &matchings) {
    unsigned int nRows = this->costMatrix.size() - 1;
    unsigned int nCols = this->costMatrix[0].size() - 1;
    // int max_dim = std::max(nRows, nCols);
    unsigned int min_dim = std::min(nRows, nCols);
    bool transpose = nRows > nCols;

    dataType cost = 0;
    for(unsigned int ind = 0; ind < asgn.size(); ++ind) {
      int indMatrix = std::min(ind, min_dim);
      int i = (!transpose) ? indMatrix : asgn[ind];
      int j = (!transpose) ? asgn[ind] : indMatrix;
      cost += this->costMatrix[i][j];
      matchings.push_back(std::make_tuple(i, j, this->costMatrix[i][j]));
    }
    return cost;
  }

  template <typename dataType>
  int AssignmentExhaustive<dataType>::run(
    std::vector<asgnMatchingTuple> &matchings) {
    int nRows = this->costMatrix.size() - 1;
    int nCols = this->costMatrix[0].size() - 1;
    int max_dim = std::max(nRows, nCols);
    int min_dim = std::min(nRows, nCols);

    // --- Construct all possible assignments
    std::vector<std::vector<int>> allAsgn;
    // We hard write the basic assignments to avoid the call of
    // constructAssignments These assignments are, of course, automatically
    // generated by constructAssignments
    if(min_dim == 1 and max_dim == 1)
      allAsgn = {{0}, {1, 0}};
    else if(min_dim == 1 and max_dim == 2)
      allAsgn = {{0, 1}, {2, 0, 1}, {1, 0}};
    else if(min_dim == 1 and max_dim == 3)
      allAsgn = {{0, 1, 2}, {3, 0, 1, 2}, {1, 0, 2}, {2, 0, 1}};
    else if(min_dim == 2 and max_dim == 2)
      allAsgn = {{0, 1}, {0, 2, 1}, {2, 1, 0}, {2, 2, 0, 1},
                 {1, 0}, {1, 2, 0}, {2, 0, 1}};
    else {
      std::stringstream ss;
      ss << min_dim << "_" << max_dim;
      std::string asgnName = ss.str();
      auto it = savedAsgn.find(asgnName);
      if(it == savedAsgn.end()) { // not found
        printMsg(asgnName, debug::Priority::VERBOSE);
        // allAsgn = constructAssignments(min_dim, max_dim);
        // allAsgn = constructAssignments2(min_dim, max_dim);
        allAsgn = constructAssignments3(min_dim, max_dim);
        savedAsgn[asgnName] = allAsgn;
        std::stringstream ss2;
        ss2 << allAsgn.size() << " done";
        printMsg(ss2.str(), debug::Priority::VERBOSE);
      } else
        allAsgn = savedAsgn[asgnName];
    }

    // --- Try these assignments and get the better
    dataType bestCost = std::numeric_limits<dataType>::max();
    std::vector<asgnMatchingTuple> bestMatching;
    for(std::vector<int> asgn : allAsgn) {
      std::vector<asgnMatchingTuple> tempMatching;
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

#endif
