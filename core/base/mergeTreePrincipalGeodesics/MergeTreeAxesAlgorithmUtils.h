#pragma once

#include <string>

namespace ttk {
  namespace axa {
    //----------------------------------------------------------------------------
    // Output Utils
    //----------------------------------------------------------------------------
    void zeroPadding(std::string &colName,
                     const size_t numberCols,
                     const size_t colIdx);

    std::string getTableCoefficientName(int noAxes, int axeNum);

    std::string getTableCoefficientNormName(int noAxes, int axeNum);

    std::string getTableVectorName(
      int noAxes, int axeNum, int vId, int vComp, bool isSecondInput = false);

    std::string getTableCorrelationName(int noAxes, int axeNum);

    std::string getTableCorrelationPersName(int noAxes, int axeNum);

    std::string getTableTreeName(int noTrees, int treeNum);
  } // namespace axa
} // namespace ttk
