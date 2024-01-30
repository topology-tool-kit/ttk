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

    std::string getTableCoefficientName(int noGeodesics, int geodesicNum);

    std::string getTableCoefficientNormName(int noGeodesics, int geodesicNum);

    std::string getTableVectorName(int noGeodesics,
                                   int geodesicNum,
                                   int vId,
                                   int vComp,
                                   bool isSecondInput = false);

    std::string getTableCorrelationName(int noGeodesics, int geodesicNum);

    std::string getTableCorrelationPersName(int noGeodesics, int geodesicNum);

    std::string getTableTreeName(int noTrees, int treeNum);
  } // namespace axa
} // namespace ttk
