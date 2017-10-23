/// \ingroup base
/// \class ttk::Munkres 
/// \author Maxime Soler <soler.maxime@total.com>

#ifndef _MUNKRES_H
#define _MUNKRES_H

#include <cmath>
#include <limits>
#include <iostream>
#include <Debug.h>
#include <DataTypes.h>

namespace ttk {
  
  class Munkres : public Debug
  {
    
    public:
    
      Munkres();
    
      ~Munkres();

      template <typename T>
        int run(vector<tuple<idVertex, idVertex, T> > *matchings);

      template<typename T>
        inline void clear() {
          if (C_orig != NULL) {
            for (int r = 0, rS0 = rowSize; r < rS0; ++r)
              delete[] ((T**)C_orig)[r];
            delete[] ((T**)C_orig);
            C_orig = NULL;
          }

          // if (Cptr != NULL) {
          //   for (int r = 0, rS0 = rowSize; r < rS0; ++r)
          //     delete[] ((T**)Cptr)[r];
          //   delete[] ((T**)Cptr);
          //   Cptr = NULL;
          // }

          if (M != NULL) {
            for (int r = 0, rS0 = rowSize; r < rS0; ++r)
              delete[] M[r];
            delete[] M;
            M = NULL;
          }

          if (path != NULL) {
            for (int p = 0, nbp = pathCount; p < nbp; ++p)
              delete[] path[p];
            delete[] path;
            path = NULL;
          }

          if (colCover != NULL)
            delete[] colCover;
          if (rowCover != NULL)
            delete[] rowCover;

          pathCount = 0;
          rowSize = 0;
          colSize = 0;
        }

      inline int setInput(int colSize_, int rowSize_, void* C_) {
        colSize = colSize_;
        rowSize = rowSize_;
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

      template<typename T>
      inline void showCostMatrix() {
        T** C = (T**) Cptr;
        stringstream msg;
        for (int r = 0; r < rowSize; ++r) {
          msg << endl << "  ";
          for (int c = 0; c < colSize; ++c)
            msg << std::scientific << C[r][c] << " ";
        }

        msg << endl;
        dMsg(std::cout, msg.str(), advancedInfoMsg);
      }

      inline void showMaskMatrix() {
        stringstream msg;
        msg << endl << "     ";
        for (int c = 0; c < colSize; ++c) {
          msg << colCover[c] << "  ";
        }
        for (int r = 0; r < rowSize; ++r) {
          msg << endl << "  " << rowCover[r] << "  ";
          for (int c = 0; c < colSize; ++c)
            msg << M[r][c] << "  ";
        }

        msg << endl;
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
      int       pathCount;
    
      int       rowSize;
      int       colSize;
    
      // internal methods
      template <typename T>
        int stepOne(int& step);

      template <typename T>
       int stepTwo(int& step);

      template <typename T>
        int stepThree(int& step);

      template <typename T>
        int stepFour(int& step);
      template <typename T>
        int findZero(int& row, int& col);
      template <typename T>
        int findStarInRow(int row);

      template <typename T>
        int stepFive(int& step);
      template <typename T>
        int findStarInCol(int col, int& row);
      template <typename T>
        int findPrimeInRow(int row, int& col);

      template <typename T>
        int stepSix(int& step);

      template <typename T>
        int stepSeven(int& step);

      template <typename T>
        int affect(vector<tuple<idVertex, idVertex, T> > *matchings);

      template <typename T>
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

      template<typename T>
        inline int copyInputMatrix() {

          T** C = (T**) Cptr;
          int rS = rowSize;
          int cS = colSize;
          T** C_orig_ = new T*[rS];
          for (int r = 0; r < rS; ++r)
          {
            C_orig_[r] = new T[cS];
            for (int c = 0; c < cS; ++c)
              C_orig_[r][c] = C[r][c];
          }
          C_orig = (void*) C_orig_;
          return 0;
        }

  };

  template <typename T>
  int Munkres::run(vector<tuple<idVertex, idVertex, T> > *matchings)
  {
    int step = 1;
    int iter = 0;
    bool done = false;
    bool debug = false;
    Timer t;

    copyInputMatrix<T>();

    while (!done)
    {
      ++ iter;
      {
        stringstream msg;
        msg << "[BottleneckDistence] Munkres step " << step << ", Iteration "
            << iter << " in " << t.getElapsedTime() << endl;
        dMsg(std::cout, msg.str(), advancedInfoMsg);
      }

      if (iter % 20 == 0) {
        stringstream msg;
        double progress = std::round((double)iter)/5.0;
        msg << "[BottleneckDistance] Munkres progress " << progress
            << "%" << endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }

      if (iter > 500) {
        // debug = true;
        step = 7;
        // Abort. Still found something
        // albeit hardly optimized.
      }

      if (debug) {
        showCostMatrix<T>();
        showMaskMatrix();
      }

      switch(step)
      {
        case 1:
          stepOne<T>(step);
          break;
        case 2:
          stepTwo<T>(step);
          break;
        case 3:
          stepThree<T>(step);
          break;
        case 4:
          stepFour<T>(step);
          break;
        case 5:
          stepFive<T>(step);
          break;
        case 6:
          stepSix<T>(step);
          break;
        case 7:
          stepSeven<T>(step);
          done = true;
          break;
      }
    }

    computeAffectationCost<T>();
    affect<T>(matchings);
    clear<T>();

    return 0;

  }

  // Substract minimum value in every row.
  template <typename T>
  int Munkres::stepOne(int& step) {
    double minInRow;
    T** C = (T**) Cptr;

    for (int r = 0; r < rowSize; ++r) {
      minInRow = C[r][0];

      for (int c = 0; c < colSize; ++c)
        if (C[r][c] < minInRow) minInRow = C[r][c];

      for (int c = 0; c < colSize; ++c)
        C[r][c] -= minInRow;
    }

    step = 2;
    return 0;
  }

  // Find a zero in the matrix,
  // star it if it is the only one in its row and col.
  template <typename T>
  int Munkres::stepTwo(int& step) {
    T** C = (T**) Cptr;

    for (int r = 0; r < rowSize; ++r)
      for (int c = 0; c < colSize; ++c)
        if (C[r][c] == 0. && !rowCover[r] && !colCover[c]) {
          // TODO check == 0. criterium
          M[r][c] = 1;
          rowCover[r] = true;
          colCover[c] = true;
        }

    for (int r = 0; r < rowSize; ++r)
      rowCover[r] = false;
    for (int c = 0; c < colSize; ++c)
      colCover[c] = false;

    step = 3;
    return 0;
  }

  // Check column coverings.
  // If all columns are starred (1 star only per column is possible)
  // then the algorithm is terminated.
  template <typename T>
  int Munkres::stepThree(int& step) {
    // T** C = (T**) Cptr;

    for (int r = 0; r < rowSize; ++r)
      for (int c = 0; c < colSize; ++c)
        if (M[r][c] == 1)
          colCover[c] = true;

    int processedCols = 0;
    for (int c = 0; c < colSize; ++c)
      if (colCover[c]) ++processedCols;

    if (processedCols >= colSize || processedCols >= rowSize)
      step = 7; // end algorithm
    else
      step = 4; // follow prime scheme
    return 0;
  }

  // Find noncovered zero, prime it
  // . if current row has no starred zero -> step 5
  // . else, cover row and uncover the col with a star
  // Repeat until there are no uncovered zero left
  // Save smallest uncovered value then -> step 6
  template <typename T>
  int Munkres::stepFour(int& step) {
    // T** C = (T**) Cptr;

    int row = -1;
    int col = -1;
    bool done = false;

    while (!done) {

      findZero<T>(row, col);
      if (row == -1) {
        done = true;
        step = 6;
      }

      else {
        M[row][col] = 2;
        int colOfStarInRow = findStarInRow<T>(row);
        if (colOfStarInRow > -1) { // found a star
          rowCover[row] = true;
          colCover[colOfStarInRow] = false;
        }

        else {
          done = true;
          step = 5;
          pathRow0 = row;
          pathCol0 = col;
        }
      }
    }

    return 0;
  }

  template <typename T>
  int Munkres::findStarInRow(int row) {
    int col = -1;
    for (int c = 0; c < colSize; ++c)
      if (M[row][c] == 1) col = c;
    return col;
  }

  template <typename T>
  int Munkres::findZero(int& row, int& col) {
    T** C = (T**) Cptr;

    row = -1;
    col = -1;

    for (int r = 0; r < rowSize; ++r) {
      if (rowCover[r]) continue;

      for (int c = 0; c < colSize; ++c) {
        if (colCover[c]) continue;
        if (C[r][c] == (T) 0) {
          row = r;
          col = c;
          return 0;
        }
      }
    }

    stringstream msg;
    msg << "[Munkres] Zero not found" << std::endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);

    return 0;
  }

  // Make path of alternating primed and starred zeros
  // 1. uncovered primed found at step 4
  // 2. same column, starred (if any)
  // 3. same row, primed (always one)
  // 4. continue until a primed zero has no starred zero in its column
  // Unstar each starred zero in the series, star each primed zero
  // in the series,
  // erase all primes, uncover every line, return to step 3.
  template <typename T>
  int Munkres::stepFive(int& step) {
    {
      // T** C = (T**) Cptr;

      int r = -1;
      int c = -1;

      pathCount = 1;
      path[pathCount - 1][0] = pathRow0;
      path[pathCount - 1][1] = pathCol0;

      bool done = false;
      while (!done) {
        findStarInCol<T>(path[pathCount - 1][1], r);
        if (r == -1)
          done = true;

        else {
          ++pathCount;
          path[pathCount - 1][0] = r;
          path[pathCount - 1][1] = path[pathCount - 2][1];

          findPrimeInRow<T>(path[pathCount - 1][0], c);
          ++pathCount;
          path[pathCount - 1][0] = path[pathCount - 2][0];
          path[pathCount - 1][1] = c;
        }
      }
    }

    // process path
    for (int p = 0; p < pathCount; ++p) {
      if (M[path[p][0]][path[p][1]] == 1)
        M[path[p][0]][path[p][1]] = 0;
      else
        M[path[p][0]][path[p][1]] = 1;
    }

    // clear covers
    for (int r = 0; r < rowSize; ++r) rowCover[r] = false;
    for (int c = 0; c < colSize; ++c) colCover[c] = false;

    // erase primes
    for (int r = 0; r < rowSize; ++r)
      for (int c = 0; c < colSize; ++c)
        if (M[r][c] == 2) M[r][c] = 0;

    step = 3;
    return 0;
  }

  template <typename T>
  int Munkres::findStarInCol(int col, int& row) {
    // T** C = (T**) Cptr;

    row = -1;
    for (int r = 0; r < rowSize; ++r)
      if (M[r][col] == 1) row = r;

    return 0;
  }

  template <typename T>
  int Munkres::findPrimeInRow(int row, int& col) {
    // T** C = (T**) Cptr;

    for (int c = 0; c < colSize; ++c)
      if (M[row][c] == 2) col = c;

    return 0;
  }

  // Add value found in step 4 to every element of each covered row,
  // substract it from every element of each uncovered col.
  // Return to step 4 without altering any stars/primes/covers.
  template <typename T>
  int Munkres::stepSix(int& step) {
    T** C = (T**) Cptr;

    T min = std::numeric_limits<T>::max();

    // find smallest
    for (int r = 0; r < rowSize; ++r) {
      if (rowCover[r]) continue;

      for (int c = 0; c < colSize; ++c) {
        if (colCover[c]) continue;
        else if (C[r][c] < min) min = C[r][c];
      }
    }

    // add and substract
    for (int r = 0; r < rowSize; ++r) {
      for (int c = 0; c < colSize; ++c) {
        if (rowCover[r]) C[r][c] += min;
        if (!colCover[c]) C[r][c] -= min;
      }
    }

    step = 4;
    return 0;
  }

  template <typename T>
  int Munkres::stepSeven(int& step) {
    stringstream msg;
    msg << "[Munkres] Step 7 over." << endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);
    return 0;
  }

  template<typename T>
  int Munkres::affect(vector<tuple<idVertex, idVertex, T>> *matchings) {

    T** C = (T**) C_orig;
    int nbC = colSize;
    int nbR = rowSize;

    matchings->clear();

    for (int r = 0; r < nbR; ++r)
      for (int c = 0; c < nbC; ++c)
        if (M[r][c] == 1) {
          tuple<idVertex, idVertex, T> t = make_tuple(r, c, C[r][c]);
          matchings->push_back(t);
        }

    return 0;
  }

  template<typename T>
  int Munkres::computeAffectationCost() {

    T** C = (T**) C_orig;
    int nbC = colSize;
    int nbR = rowSize;

    T total = 0;
    for (int r = 0; r < nbR; ++r)
      for (int c = 0; c < nbC; ++c)
        if (M[r][c] == 1) {
          total += C[r][c];
        }

    stringstream msg;
    msg << "[Munkres] Total cost = " << total << endl;
    dMsg(std::cout, msg.str(), timeMsg);

    return 0;
  }

}

#endif
