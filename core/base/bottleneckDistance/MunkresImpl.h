#ifndef _MUNKRESIMPL_H
#define _MUNKRESIMPL_H

#ifndef matchingTuple
#define matchingTuple std::tuple<ftm::idVertex, ftm::idVertex, dataType>
#endif

template <typename dataType>
int ttk::Munkres::run(std::vector<matchingTuple> *matchings)
{
  int step = 1;
  int iter = 0;
  int maxIter = 5000;
  bool done = false;
  Timer t;

  copyInputMatrix<dataType>();

  while (!done)
  {
    ++ iter;
    {
      std::stringstream msg;
      msg << "[BottleneckDistence] Munkres step " << step << ", Iteration "
          << iter << " in " << t.getElapsedTime() << std::endl;
      dMsg(std::cout, msg.str(), advancedInfoMsg);
    }

    if (iter > 20 && (iter % (int) std::round((double) maxIter / 5.0) == 0)) {
      std::stringstream msg;
      double progress = std::round(100.0 * (double) iter / (double) maxIter);
      msg << "[BottleneckDistance] Munkres progress " << progress
          << "%" << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }

    if (iter > maxIter) {
      showCostMatrix<dataType>();
      showMaskMatrix();

      {
        std::stringstream msg;
        msg << "[BottleneckDistance] Munkres failed to converge after " << 
(maxIter)
            << " iterations. Aborting." << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
      }

      step = 7;
      // Abort. Still found something
      // albeit hardly optimized.
    }

    switch(step)
    {
      case 1:
        stepOne<dataType>(step);
        break;
      case 2:
        stepTwo<dataType>(step);
        break;
      case 3:
        stepThree<dataType>(step);
        break;
      case 4:
        stepFour<dataType>(step);
        break;
      case 5:
        stepFive<dataType>(step);
        break;
      case 6:
        stepSix<dataType>(step);
        break;
      case 7:
        stepSeven<dataType>(step);
        done = true;
        break;
      default:
        break;
    }
  }

  this->computeAffectationCost<dataType>();
  this->affect<dataType>(matchings);
  this->clear<dataType>();

  return 0;

}

// Substract minimum value in every row.
template <typename dataType>
int ttk::Munkres::stepOne(int& step) {
  double minInRow;
  auto** C = (dataType**) Cptr;

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
template <typename dataType>
int ttk::Munkres::stepTwo(int& step) {
  auto** C = (dataType**) Cptr;

  for (int r = 0; r < rowSize; ++r)
    for (int c = 0; c < colSize; ++c)
      // TODO Weakness: == 0. criterium
      if (C[r][c] == 0. && !rowCover[r] && !colCover[c]) {
        M[r][c] = 1;
        // [CORE] Important! here diagonal values shouldn't be discarded in 
        // rowCover.
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
template <typename dataType>
int ttk::Munkres::stepThree(int& step) {
  // [CORE] Important! Here, matchings to the diagonal are accounted for.
  for (int r = 0; r < rowSize; ++r)
    for (int c = 0; c < colSize; ++c)
      if (M[r][c] == 1)
        colCover[c] = true;

  int processedCols = 0;
  for (int c = 0; c < colSize-1; ++c)
    if (colCover[c]) ++processedCols;

  if (processedCols >= colSize-1 || processedCols >= rowSize-1)
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
template <typename dataType>
int ttk::Munkres::stepFour(int& step) {
  int row = -1;
  int col = -1;
  bool done = false;

  while (!done) {

    findZero<dataType>(row, col);
    if (row == -1) {
      done = true;
      step = 6;
    }

    else {
      M[row][col] = 2;
      int colOfStarInRow = findStarInRow<dataType>(row);
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

template <typename dataType>
int ttk::Munkres::findStarInRow(int row) {
  int col = -1;
  for (int c = 0; c < colSize; ++c)
    if (M[row][c] == 1) col = c;
  return col;
}

template <typename dataType>
int ttk::Munkres::findZero(int& row, int& col) {
  auto** C = (dataType**) Cptr;

  row = -1;
  col = -1;

  for (int r = 0; r < rowSize; ++r) {
    if (rowCover[r]) continue;

    for (int c = 0; c < colSize; ++c) {
      if (colCover[c]) continue;
      if (C[r][c] == (dataType) 0) {
        row = r;
        col = c;
        return 0;
      }
    }
  }

  std::stringstream msg;
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
template <typename dataType>
int ttk::Munkres::stepFive(int& step) {
  {
    int r = -1;
    int c = -1;

    pathCount = 1;
    path[pathCount - 1][0] = pathRow0;
    path[pathCount - 1][1] = pathCol0;

    bool done = false;
    while (!done) {
      findStarInCol<dataType>(path[pathCount - 1][1], r);
      if (r == -1)
        done = true;

      else {
        ++pathCount;
        path[pathCount - 1][0] = r;
        path[pathCount - 1][1] = path[pathCount - 2][1];

        findPrimeInRow<dataType>(path[pathCount - 1][0], c);
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

template <typename dataType>
int ttk::Munkres::findStarInCol(int col, int& row) {
  row = -1;
  for (int r = 0; r < rowSize; ++r)
    if (M[r][col] == 1) row = r;

  return 0;
}

template <typename dataType>
int ttk::Munkres::findPrimeInRow(int row, int& col) {
  for (int c = 0; c < colSize; ++c)
    if (M[row][c] == 2) col = c;

  return 0;
}

// Add value found in step 4 to every element of each covered row,
// substract it from every element of each uncovered col.
// Return to step 4 without altering any stars/primes/covers.
template <typename dataType>
int ttk::Munkres::stepSix(int& step) {
  auto** C = (dataType**) Cptr;

  dataType min = std::numeric_limits<dataType>::max();

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

template <typename dataType>
int ttk::Munkres::stepSeven(int& step) {
  std::stringstream msg;
  msg << "[Munkres] Step 7 over." << std::endl;
  dMsg(std::cout, msg.str(), advancedInfoMsg);
  return 0;
}

template<typename dataType>
int ttk::Munkres::affect(std::vector<matchingTuple> *matchings) {

  auto** C = (dataType**) C_orig;
  int nbC = colSize;
  int nbR = rowSize;

  matchings->clear();

  for (int r = 0; r < nbR; ++r)
    for (int c = 0; c < nbC; ++c)
      if (M[r][c] == 1) {
        matchingTuple t = std::make_tuple(r, c, C[r][c]);
        matchings->push_back(t);
      }

  return 0;
}

template<typename dataType>
int ttk::Munkres::computeAffectationCost() {

  auto** C = (dataType**) C_orig;
  int nbC = colSize;
  int nbR = rowSize;

  dataType total = 0;
  for (int r = 0; r < nbR; ++r)
    for (int c = 0; c < nbC; ++c)
      if (M[r][c] == 1) {
        total += C[r][c];
      }

  std::stringstream msg;
  msg << "[Munkres] Total cost = " << total << std::endl;
  dMsg(std::cout, msg.str(), timeMsg);

  return 0;
}

#endif
