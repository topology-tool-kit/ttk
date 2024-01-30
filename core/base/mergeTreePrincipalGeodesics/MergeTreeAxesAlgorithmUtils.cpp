#include <MergeTreeAxesAlgorithmUtils.h>

namespace ttk {
  namespace axa {
    //----------------------------------------------------------------------------
    // Output Utils
    //----------------------------------------------------------------------------
    void zeroPadding(std::string &colName,
                     const size_t numberCols,
                     const size_t colIdx) {
      std::string const max{std::to_string(numberCols - 1)};
      std::string const cur{std::to_string(colIdx)};
      std::string const zer(max.size() - cur.size(), '0');
      colName.append(zer).append(cur);
    }

    std::string getTableCoefficientName(int noGeodesics, int geodesicNum) {
      std::string name{"T"};
      zeroPadding(name, noGeodesics, geodesicNum);
      return name;
    }

    std::string getTableCoefficientNormName(int noGeodesics, int geodesicNum) {
      std::string name{"TNorm"};
      zeroPadding(name, noGeodesics, geodesicNum);
      return name;
    }

    std::string getTableVectorName(int noGeodesics,
                                   int geodesicNum,
                                   int vId,
                                   int vComp,
                                   bool isSecondInput) {
      std::string indexString{};
      zeroPadding(indexString, noGeodesics, geodesicNum);
      std::string const prefix{(isSecondInput ? "T2_" : "")};
      std::string name{prefix + "V" + indexString + "_" + std::to_string(vId)
                       + "_" + std::to_string(vComp)};
      return name;
    }

    std::string getTableCorrelationName(int noGeodesics, int geodesicNum) {
      std::string name{"Corr"};
      zeroPadding(name, noGeodesics, geodesicNum);
      return name;
    }

    std::string getTableCorrelationPersName(int noGeodesics, int geodesicNum) {
      std::string name{"CorrPers"};
      zeroPadding(name, noGeodesics, geodesicNum);
      return name;
    }

    std::string getTableTreeName(int noTrees, int treeNum) {
      std::string name{"Tree"};
      zeroPadding(name, noTrees, treeNum);
      return name;
    }
  } // namespace axa
} // namespace ttk
