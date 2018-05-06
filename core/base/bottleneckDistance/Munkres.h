/// \ingroup base
/// \class ttk::Munkres 
/// \author Maxime Soler <soler.maxime@total.com>

#ifndef _MUNKRES_H
#define _MUNKRES_H

#ifndef matchingTuple
#define matchingTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::idVertex, dataType>
#endif

#include <cmath>
#include <limits>
#include <iostream>
#include <Debug.h>
#include <PersistenceDiagram.h>

namespace ttk {

  class Munkres : public Debug
  {

    public:

      Munkres() {};

      ~Munkres() {};

      template <typename dataType>
      int run(std::vector<matchingTuple> *matchings);

      template<typename dataType>
      inline void clear() {
        if (C_orig != nullptr) {
          for (int r = 0, rS0 = rowSize; r < rS0; ++r)
            delete[] ((dataType**)C_orig)[r];
          delete[] ((dataType**)C_orig);
          C_orig = nullptr;
        }

        if (M != nullptr) {
          for (int r = 0, rS0 = rowSize; r < rS0; ++r)
            delete[] M[r];
          delete[] M;
          M = nullptr;
        }

        if (path != nullptr) {
          for (int p = 0, nbp = pathCount; p < nbp; ++p)
            delete[] path[p];
          delete[] path;
          path = nullptr;
        }

        delete[] colCover;
        delete[] rowCover;

        pathCount = 0;
        rowSize = 0;
        colSize = 0;
      }

      template<typename dataType>
      inline void clearMatrix() {
        if (Cptr != nullptr) {
          for (int r = 0, rS0 = rowSize; r < rS0; ++r)
            delete[] ((dataType**)Cptr)[r];
          delete[] ((dataType**)Cptr);
          Cptr = nullptr;
        }
      }

      inline int setInput(int rowSize_, int colSize_, void* C_) {
        rowSize = rowSize_;
        colSize = colSize_;
        Cptr = C_;

        rowCover = new bool[rowSize_];
        colCover = new bool[colSize_];

        M = new int*[rowSize_];
        for (int r = 0; r < rowSize_; ++r)
          M[r] = new int[colSize_];

        int nbPaths = 1 + colSize_ + rowSize_;
        path = new int*[nbPaths];
        for (int p = 0; p < nbPaths; ++p)
          path[p] = new int[2];

        resetMasks();

        return 0;
      }

      template<typename dataType>
      inline void showCostMatrix() {
        auto** C = (dataType**) Cptr;
        std::stringstream msg;
        for (int r = 0; r < rowSize; ++r) {
          msg << std::endl << "  ";
          for (int c = 0; c < colSize; ++c)
            msg << std::scientific << C[r][c] << " ";
        }

        msg << std::endl;
        dMsg(std::cout, msg.str(), advancedInfoMsg);
      }

      inline void showMaskMatrix() {
        std::stringstream msg;
        msg << std::endl << "     ";
        for (int c = 0; c < colSize; ++c) {
          msg << colCover[c] << "  ";
        }
        for (int r = 0; r < rowSize; ++r) {
          msg << std::endl << "  " << rowCover[r] << "  ";
          for (int c = 0; c < colSize; ++c)
            msg << M[r][c] << "  ";
        }

        msg << std::endl;
        dMsg(std::cout, msg.str(), advancedInfoMsg);
      }

    private:

      void      *Cptr;
      void      *C_orig;

      int       **M;
      bool      *rowCover;
      bool      *colCover;

      int       **path;
      int       pathRow0;
      int       pathCol0;
      int       pathCount = 0;

      int       rowSize = 0;
      int       colSize = 0;

      // internal methods
      template <typename dataType>
      int stepOne(int& step);

      template <typename dataType>
      int stepTwo(int& step);

      template <typename dataType>
      int stepThree(int& step);

      template <typename dataType>
      int stepFour(int& step);
      template <typename dataType>
      int findZero(int& row, int& col);
      template <typename dataType>
      int findStarInRow(int row);

      template <typename dataType>
      int stepFive(int& step);
      template <typename dataType>
      int findStarInCol(int col, int& row);
      template <typename dataType>
      int findPrimeInRow(int row, int& col);

      template <typename dataType>
      int stepSix(int& step);

      template <typename dataType>
      int stepSeven(int& step);

      template <typename dataType>
      int affect(std::vector<matchingTuple> *matchings);

      template <typename dataType>
      int computeAffectationCost();

      inline int resetMasks()
      {
        for (int r = 0; r < rowSize; ++r) {
          rowCover[r] = false;
          for (int c = 0; c < colSize; ++c)
            M[r][c] = 0;
        }
        for (int c = 0; c < colSize; ++c)
          colCover[c] = false;
        return 0;
      }

      template<typename dataType>
      inline int copyInputMatrix() {

        auto** C = (dataType**) Cptr;
        int rS = rowSize;
        int cS = colSize;
        auto** C_orig_ = new dataType*[rS];
        for (int r = 0; r < rS; ++r)
        {
          C_orig_[r] = new dataType[cS];
          for (int c = 0; c < cS; ++c)
            C_orig_[r][c] = C[r][c];
        }
        C_orig = (void*) C_orig_;
        return 0;
      }

  };

}

#include <MunkresImpl.h>

#endif
