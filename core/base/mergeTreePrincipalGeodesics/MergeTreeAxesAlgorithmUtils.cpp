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

    std::string getTableCoefficientName(int noAxes, int axeNum) {
      std::string name{"T"};
      zeroPadding(name, noAxes, axeNum);
      return name;
    }

    std::string getTableCoefficientNormName(int noAxes, int axeNum) {
      std::string name{"TNorm"};
      zeroPadding(name, noAxes, axeNum);
      return name;
    }

    std::string getTableVectorName(
      int noAxes, int axeNum, int vId, int vComp, bool isSecondInput) {
      std::string indexString{};
      zeroPadding(indexString, noAxes, axeNum);
      std::string const prefix{(isSecondInput ? "T2_" : "")};
      std::string name{prefix + "V" + indexString + "_" + std::to_string(vId)
                       + "_" + std::to_string(vComp)};
      return name;
    }

    std::string getTableCorrelationName(int noAxes, int axeNum) {
      std::string name{"Corr"};
      zeroPadding(name, noAxes, axeNum);
      return name;
    }

    std::string getTableCorrelationPersName(int noAxes, int axeNum) {
      std::string name{"CorrPers"};
      zeroPadding(name, noAxes, axeNum);
      return name;
    }

    std::string getTableTreeName(int noTrees, int treeNum) {
      std::string name{"Tree"};
      zeroPadding(name, noTrees, treeNum);
      return name;
    }
  } // namespace axa
} // namespace ttk
