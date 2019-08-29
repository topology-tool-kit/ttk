/// \ingroup base
/// \class ttk::Munkres
/// \author Maxime Soler <soler.maxime@total.com>

#ifndef _MUNKRES_H
#define _MUNKRES_H

#include <Debug.h>
#include <PersistenceDiagram.h>
#include <cmath>
#include <iostream>
#include <limits>

namespace ttk {

  class Munkres : public Debug {

  public:
    Munkres(){};

    ~Munkres(){};

    template <typename dataType>
    int run(std::vector<matchingTuple> &matchings);

    template <typename dataType>
    inline void clear() {
      pathCount = 0;
      rowSize = 0;
      colSize = 0;
      createdZeros.clear();
    }

    template <typename dataType>
    inline void clearMatrix() {
      std::vector<std::vector<dataType>> C
        = *((std::vector<std::vector<dataType>> *)Cptr);
      for(int r = 0, rS0 = rowSize; r < rS0; ++r)
        for(int c = 0, cS0 = colSize; c < cS0; ++c)
          C[r][c] = 0.0;
    }

    inline int setInput(int rowSize_, int colSize_, void *C_) {
      rowSize = rowSize_;
      colSize = colSize_;

      createdZeros.clear();

      Cptr = C_;

      auto r = (unsigned long)rowSize_;
      auto c = (unsigned long)colSize_;

      rowCover.resize(r);
      colCover.resize(c);

      rowLimitsMinus.resize(r);
      rowLimitsPlus.resize(r);
      colLimitsMinus.resize(c);
      colLimitsPlus.resize(c);

      M.resize(r);
      for(int row = 0; row < rowSize_; ++row)
        M[row].resize(c);

      int nbPaths = 1 + colSize_ + rowSize_;
      path.resize((unsigned long)nbPaths);
      for(int p = 0; p < nbPaths; ++p)
        path[p].resize(2);

      resetMasks();

      return 0;
    }

    template <typename dataType>
    inline void showCostMatrix() {
      auto C = (std::vector<std::vector<dataType>> *)Cptr;
      std::stringstream msg;
      for(int r = 0; r < rowSize; ++r) {
        msg << std::endl << "  ";
        for(int c = 0; c < colSize; ++c)
          msg << std::fixed << std::setprecision(3) << (*C)[r][c] << " ";
      }

      msg << std::endl;
      dMsg(std::cout, msg.str(), infoMsg);
    }

    inline void showMaskMatrix() {
      std::stringstream msg;
      msg << std::endl << "   | ";
      for(int c = 0; c < colSize; ++c) {
        msg << colCover[c] << "  ";
      }
      msg << std::endl << "   | ";
      for(int c = 0; c < colSize; ++c) {
        msg << "---";
      }
      msg << std::endl;
      for(int r = 0; r < rowSize; ++r) {
        msg << " " << rowCover[r] << " | ";
        for(int c = 0; c < colSize; ++c)
          msg << M[r][c] << "  ";
        msg << std::endl;
      }
      msg << std::endl;
      dMsg(std::cout, msg.str(), infoMsg);
    }

  private:
    void *Cptr;

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

    int rowSize = 0;
    int colSize = 0;

    // internal methods
    template <typename dataType>
    int stepOne(int &step);

    template <typename dataType>
    int stepTwo(int &step);

    template <typename dataType>
    int stepThree(int &step);

    template <typename dataType>
    int stepFour(int &step);
    template <typename dataType>
    int findZero(int &row, int &col);
    template <typename dataType>
    int findStarInRow(int row);

    template <typename dataType>
    int stepFive(int &step);
    template <typename dataType>
    int findStarInCol(int col);
    template <typename dataType>
    int findPrimeInRow(int row);

    template <typename dataType>
    int stepSix(int &step);

    template <typename dataType>
    int stepSeven(int &step);

    template <typename dataType>
    int affect(std::vector<matchingTuple> &matchings,
               const std::vector<std::vector<dataType>> &C);

    template <typename dataType>
    int computeAffectationCost(const std::vector<std::vector<dataType>> &C);

    template <typename dataType>
    inline bool isZero(dataType t) {
      // return std::abs((double) t) < 1e-15;
      return t == 0.;
    }

    inline int resetMasks() {
      for(int r = 0; r < rowSize; ++r) {
        rowCover[r] = false;
        for(int c = 0; c < colSize; ++c)
          M[r][c] = 0;
      }
      for(int c = 0; c < colSize; ++c)
        colCover[c] = false;
      return 0;
    }

    template <typename dataType>
    inline int copyInputMatrix(std::vector<std::vector<dataType>> &saveInput) {
      auto C = (std::vector<std::vector<dataType>> *)Cptr;
      int rS = rowSize;
      int cS = colSize;

      for(int r = 0; r < rS; ++r)
        for(int c = 0; c < cS; ++c)
          saveInput[r][c] = (*C)[r][c];

      return 0;
    }
  };

// Include in namespace ttk
#include <MunkresImpl.h>

} // namespace ttk

#endif
