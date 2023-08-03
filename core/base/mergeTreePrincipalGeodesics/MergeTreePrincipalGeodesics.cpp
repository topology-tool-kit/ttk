#include <MergeTreePrincipalGeodesics.h>

namespace ttk {
  double MergeTreePrincipalGeodesics::verifyOrthogonality(
    std::vector<std::vector<std::vector<double>>> &vS,
    std::vector<std::vector<std::vector<double>>> &v2s,
    bool doPrint) {
    unsigned int const geodesicNumber = vS.size() - 1;
    if(geodesicNumber == 0)
      return 0.0;
    if(doPrint)
      printMsg("Scalar Products:");
    double cost = 0.0;
    std::vector<std::vector<double>> sumVs;
    ttk::Geometry::multiAddVectorsFlatten(vS, v2s, sumVs);
    for(unsigned int i = 0; i < geodesicNumber; ++i) {
      double const scalarProd
        = ttk::Geometry::dotProduct(sumVs[i], sumVs[geodesicNumber]);
      cost += scalarProd * scalarProd;
      if(doPrint and geodesicNumber < 4) {
        std::stringstream ss;
        ss << " - " << i << " _ " << scalarProd;
        printMsg(ss.str());
      }
    }
    if(doPrint) {
      std::stringstream ss;
      ss << " - Ortho. Cost : " << cost;
      printMsg(ss.str());
    }
    return cost;
  }

  double MergeTreePrincipalGeodesics::verifyOrthogonality(
    std::vector<std::vector<std::vector<double>>> &vS,
    std::vector<std::vector<std::vector<double>>> &v2s,
    std::vector<std::vector<double>> &v,
    std::vector<std::vector<double>> &v2,
    bool doPrint) {
    std::vector<std::vector<std::vector<double>>> vSTemp = vS, v2sTemp = v2s;
    vSTemp.push_back(v);
    v2sTemp.push_back(v2);
    return verifyOrthogonality(vSTemp, v2sTemp, doPrint);
  }
} // namespace ttk
