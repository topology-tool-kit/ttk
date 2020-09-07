/// \ingroup base
/// \class ttk::RangeMinimumQuery
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief Class to answer range minimum queries in an array in constant time
/// after a linearithmic time preprocess.

#pragma once

#include <Debug.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace ttk {

  template <class DataType>
  class RangeMinimumQuery : virtual public Debug {
  public:
    RangeMinimumQuery();
    RangeMinimumQuery(std::vector<DataType> &input);

    inline void setVector(std::vector<DataType> &input) {
      input_ = input.data();
      input_end_ = input.data() + input.size();
    }

    int preprocess(const bool silent = false);
    int query(int i, int j) const;

  protected:
    // Input vector
    DataType *input_{};
    DataType *input_end_{};

    // Sparse Table
    std::vector<std::vector<int>> table_{};
  };

} // namespace ttk

/* Function definitions */

// Constructors
template <class DataType>
ttk::RangeMinimumQuery<DataType>::RangeMinimumQuery() {
  this->setDebugMsgPrefix("RangeMinimumQuery");
}
template <class DataType>
ttk::RangeMinimumQuery<DataType>::RangeMinimumQuery(
  std::vector<DataType> &input) {
  setVector(input);
}

// Preprocessing
template <class DataType>
int ttk::RangeMinimumQuery<DataType>::preprocess(const bool silent) {

  Timer t;

  // Compute the size of the matrix
  int sizeOfArray = static_cast<int>(input_end_ - input_);
  int numberOfBlocs = static_cast<unsigned int>(log2(sizeOfArray + 1)) + 1;

  // Init the matrix
  table_.resize(sizeOfArray);
  for(int i = 0; i < sizeOfArray; i++) {
    table_[i].resize(numberOfBlocs);
    fill(table_[i].begin(), table_[i].end(), -1);
  }

  // Initialize for blocs of size 1 (k==0)
  for(int i = 0; i < sizeOfArray; i++) {
    table_[i][0] = i;
  }
  // Compute other values recursively
  for(int k = 1; ((1 << k) - 1) < sizeOfArray; k++) {
    for(int i = 0; (i + (1 << k) - 1) < sizeOfArray; i++) {
      if(input_[table_[i][k - 1]]
         <= input_[table_[i + (1 << (k - 1))][k - 1]]) {
        table_[i][k] = table_[i][k - 1];
      } else {
        table_[i][k] = table_[i + (1 << (k - 1))][k - 1];
      }
    }
  }
  // Debug messages
  if(!silent) {
    this->printMsg("Preprocessed queries.", 1.0, t.getElapsedTime(), 1);
  }
  return 0;
}

template <class DataType>
int ttk::RangeMinimumQuery<DataType>::query(int i, int j) const {

#ifndef TTK_ENABLE_KAMIKAZE
  // i must be lower or equal to j, else swap values
  if(i > j) {
    std::swap(i, j);
  }
#endif

  // Compute size of blocs (2^k) to use
  int k = static_cast<int>(log2(j - i + 1));
  // Compute the range minimum
  if(input_[table_[i][k]] <= input_[table_[j - (1 << k) + 1][k]]) {
    return table_[i][k];
  } else {
    return table_[j - (1 << k) + 1][k];
  }
}
