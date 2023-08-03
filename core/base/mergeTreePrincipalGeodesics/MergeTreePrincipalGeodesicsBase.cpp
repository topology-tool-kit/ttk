#include <MergeTreePrincipalGeodesicsBase.h>

// other includes
#include <Geometry.h>
#include <vector>

namespace ttk {
  // ----------------------------------------------------------------------------
  // Vector Utils
  // ----------------------------------------------------------------------------
  void MergeTreePrincipalGeodesicsBase::callGramSchmidt(
    std::vector<std::vector<double>> &vS,
    std::vector<double> &v,
    std::vector<double> &newV) {
    std::vector<std::vector<double>> allVs = vS, uS;
    allVs.push_back(v);
    ttk::Geometry::gramSchmidt(allVs, uS);
    newV = uS[uS.size() - 1];
  }

  void
    MergeTreePrincipalGeodesicsBase::vectorToPointer(std::vector<double> &vec,
                                                     double *&pVec) {
    pVec = vec.data();
  }

  void MergeTreePrincipalGeodesicsBase::vectorsToPointers(
    std::vector<std::vector<double>> &vec, std::vector<double *> &pVec) {
    pVec.resize(vec.size());
    for(unsigned int i = 0; i < vec.size(); ++i)
      vectorToPointer(vec[i], pVec[i]);
  }

  void MergeTreePrincipalGeodesicsBase::vectorOfVectorsToPointers(
    std::vector<std::vector<std::vector<double>>> &vS,
    std::vector<std::vector<double *>> &pVS) {
    pVS.resize(vS.size());
    for(unsigned int i = 0; i < vS.size(); ++i)
      vectorsToPointers(vS[i], pVS[i]);
  }

  void MergeTreePrincipalGeodesicsBase::pointerToVector(
    double *pVec, size_t size, std::vector<double> &vec) {
    vec.resize(size);
    for(unsigned int i = 0; i < size; ++i)
      vec[i] = pVec[i];
  }

  void MergeTreePrincipalGeodesicsBase::pointersToVectors(
    std::vector<double *> &pVec,
    std::vector<size_t> sizes,
    std::vector<std::vector<double>> &vec) {
    vec.resize(pVec.size());
    for(unsigned int i = 0; i < pVec.size(); ++i)
      pointerToVector(pVec[i], sizes[i], vec[i]);
  }

  void MergeTreePrincipalGeodesicsBase::pointersToVectors(
    std::vector<double *> &pVec,
    size_t size,
    std::vector<std::vector<double>> &vec) {
    std::vector<size_t> const sizes(pVec.size(), size);
    pointersToVectors(pVec, sizes, vec);
  }
} // namespace ttk
