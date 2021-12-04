/// \ingroup base
/// \class ttk::AssignmentMunkres
/// \author Maxime Soler <soler.maxime@total.com>

#pragma once

#include <AssignmentSolver.h>

#include <cmath>
#include <iostream>
#include <limits>

namespace ttk {

  template <class dataType>
  class AssignmentMunkres : virtual public Debug,
                            public AssignmentSolver<dataType> {

  public:
    AssignmentMunkres() {
      this->setDebugMsgPrefix("AssignmentMunkres");
    }

    ~AssignmentMunkres() = default;

    int run(std::vector<asgnMatchingTuple> &matchings);

    inline void clear() {
      AssignmentSolver<dataType>::clear();
      pathCount = 0;
      createdZeros.clear();
    }

    inline int setInput(std::vector<std::vector<dataType>> &C_) {
      AssignmentSolver<dataType>::setInput(C_);

      createdZeros.clear();

      auto r = (unsigned long)this->rowSize;
      auto c = (unsigned long)this->colSize;

      rowCover.resize(r);
      colCover.resize(c);

      rowLimitsMinus.resize(r);
      rowLimitsPlus.resize(r);
      colLimitsMinus.resize(c);
      colLimitsPlus.resize(c);

      M.resize(r);
      for(int row = 0; row < this->rowSize; ++row)
        M[row].resize(c);

      int nbPaths = 1 + this->colSize + this->rowSize;
      path.resize((unsigned long)nbPaths);
      for(int p = 0; p < nbPaths; ++p)
        path[p].resize(2);

      resetMasks();

      return 0;
    }

    inline void showCostMatrix() {
      auto C = AssignmentSolver<dataType>::getCostMatrix();
      std::stringstream msg;
      for(int r = 0; r < this->rowSize; ++r) {
        msg << std::endl << "  ";
        for(int c = 0; c < this->colSize; ++c)
          msg << std::fixed << std::setprecision(3) << (*C)[r][c] << " ";
      }

      this->printMsg(msg.str(), debug::Priority::DETAIL);
    }

    inline void showMaskMatrix() {
      std::stringstream msg;
      msg << std::endl << "   | ";
      for(int c = 0; c < this->colSize; ++c) {
        msg << colCover[c] << "  ";
      }
      msg << std::endl << "   | ";
      for(int c = 0; c < this->colSize; ++c) {
        msg << "---";
      }
      msg << std::endl;
      for(int r = 0; r < this->rowSize; ++r) {
        msg << " " << rowCover[r] << " | ";
        for(int c = 0; c < this->colSize; ++c)
          msg << M[r][c] << "  ";
        msg << std::endl;
      }
      this->printMsg(msg.str(), debug::Priority::DETAIL);
    }

  private:
    std::vector<std::vector<int>> M;
    std::vector<bool> rowCover;
    std::vector<bool> colCover;

    std::vector<int> rowLimitsMinus;
    std::vector<int> rowLimitsPlus;
    std::vector<int> colLimitsMinus;
    std::vector<int> colLimitsPlus;

    std::vector<std::vector<int>> path;
    std::vector<std::pair<int, int>> createdZeros;

    int pathRow0;
    int pathCol0;
    int pathCount = 0;

    // internal methods
    int stepOne(int &step);

    int stepTwo(int &step);

    int stepThree(int &step);

    int stepFour(int &step);
    int findZero(int &row, int &col);
    int findStarInRow(int row);

    int stepFive(int &step);
    int findStarInCol(int col);
    int findPrimeInRow(int row);

    int stepSix(int &step);

    int stepSeven(int &step);

    int affect(std::vector<asgnMatchingTuple> &matchings,
               const std::vector<std::vector<dataType>> &C);

    int computeAffectationCost(const std::vector<std::vector<dataType>> &C);

    inline bool isZero(dataType t) {
      // return std::abs((double) t) < 1e-15;
      return t == 0.;
    }

    inline int resetMasks() {
      for(int r = 0; r < this->rowSize; ++r) {
        rowCover[r] = false;
        for(int c = 0; c < this->colSize; ++c)
          M[r][c] = 0;
      }
      for(int c = 0; c < this->colSize; ++c)
        colCover[c] = false;
      return 0;
    }

    inline int copyInputMatrix(std::vector<std::vector<dataType>> &saveInput) {
      auto C = AssignmentSolver<dataType>::getCostMatrixPointer();
      int rS = this->rowSize;
      int cS = this->colSize;

      for(int r = 0; r < rS; ++r)
        for(int c = 0; c < cS; ++c)
          saveInput[r][c] = (*C)[r][c];

      return 0;
    }
  };

} // namespace ttk

#include <AssignmentMunkresImpl.h>
