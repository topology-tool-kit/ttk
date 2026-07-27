/// \ingroup base
/// \class ttk::AssignmentSolver
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
///
/// Assignment Problem Solver abstract class
///
/// For the unbalanced problem:
///   The cost matrix in input has a size of (n + 1) x (m + 1)
///   - n is the number of jobs, m the number of workers
///   - the nth row contains the cost of not assigning workers
///   - the mth column is the same but with jobs
///   - the last cell (costMatrix[n][m]) is not used

#pragma once

#include <Debug.h>
#include <PersistenceDiagramUtils.h>

namespace ttk {

  template <class dataType>
  class AssignmentSolver : virtual public Debug {

  public:
    AssignmentSolver() = default;

    ~AssignmentSolver() override = default;

    virtual int run(std::vector<MatchingType> &matchings) = 0;

    virtual inline void clear() {
      rowSize = 0;
      colSize = 0;
    }

    virtual inline void clearMatrix() {
      std::vector<std::vector<dataType>> C = getCostMatrix();
      for(int r = 0, rS0 = rowSize; r < rS0; ++r)
        for(int c = 0, cS0 = colSize; c < cS0; ++c)
          C[r][c] = 0.0;
    }

    virtual inline int setInput(std::vector<std::vector<dataType>> &C_) {
      rowSize = C_.size();
      colSize = C_[0].size();

      costMatrix = C_;

      setBalanced((this->rowSize == this->colSize));

      return 0;
    }

    virtual inline void setBalanced(bool balanced) {
      balancedAssignment = balanced;
    }

    virtual inline std::vector<std::vector<dataType>> getCostMatrix() {
      return costMatrix;
    }

    virtual inline std::vector<std::vector<dataType>> *getCostMatrixPointer() {
      return &costMatrix;
    }

    void printTableVector(std::vector<std::vector<dataType>> &table) {
      for(auto vecTemp : table) {
        const std::stringstream ss;
        for(auto valTemp : vecTemp) {
          ss << valTemp << " ";
        }
        printMsg(ss.str(), debug::Priority::VERBOSE);
      }
      printMsg(debug::Separator::L1, debug::Priority::VERBOSE);
    }

  protected:
    std::vector<std::vector<dataType>> costMatrix;

    int rowSize = 0;
    int colSize = 0;

    bool balancedAssignment;
  };
} // namespace ttk
