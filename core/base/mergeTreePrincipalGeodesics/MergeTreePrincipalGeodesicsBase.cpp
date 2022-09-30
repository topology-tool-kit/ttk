#include <MergeTreePrincipalGeodesicsBase.h>

// other includes
#include <Geometry.h>
#include <vector>

namespace ttk {
  // ----------------------------------------------------------------------------
  // Vector Utils
  // ----------------------------------------------------------------------------
  void
    MergeTreePrincipalGeodesicsBase::sumProjection(std::vector<double> &v1,
                                                   std::vector<double> &v2,
                                                   std::vector<double> &v1Out,
                                                   std::vector<double> &v2Out) {
    std::vector<double> sumV;
    ttk::Geometry::addVectors(v1, v2, sumV);
    ttk::Geometry::vectorProjection(v1, sumV, v1Out);
    ttk::Geometry::vectorProjection(v2, sumV, v2Out);
  }

  void MergeTreePrincipalGeodesicsBase::gramSchmidt(
    std::vector<std::vector<double>> &vS,
    std::vector<double> &v,
    std::vector<double> &newV) {
    std::vector<std::vector<double>> allVs = vS, uS;
    allVs.push_back(v);
    uS = allVs;
    for(unsigned int i = 1; i < allVs.size(); ++i) {
      std::vector<double> projecSum;
      ttk::Geometry::vectorProjection(allVs[i], uS[0], projecSum);
      for(unsigned int j = 1; j < i; ++j) {
        std::vector<double> projecTemp, projecSumTemp;
        ttk::Geometry::vectorProjection(allVs[i], uS[j], projecTemp);
        ttk::Geometry::addVectors(projecSum, projecTemp, projecSumTemp);
        projecSum = projecSumTemp;
      }
      ttk::Geometry::subtractVectors(allVs[i], projecSum, uS[i]);
    }
    newV = uS[uS.size() - 1];
  }

  void MergeTreePrincipalGeodesicsBase::multiSumVector(
    std::vector<std::vector<double>> &v1,
    std::vector<std::vector<double>> &v2,
    std::vector<std::vector<double>> &sumV) {
    sumV.resize(v1.size());
    for(unsigned int i = 0; i < v1.size(); ++i)
      ttk::Geometry::addVectors(v1[i], v2[i], sumV[i]);
  }

  void MergeTreePrincipalGeodesicsBase::multiSumVectorFlatten(
    std::vector<std::vector<std::vector<double>>> &v1,
    std::vector<std::vector<std::vector<double>>> &v2,
    std::vector<std::vector<double>> &sumV) {
    std::vector<std::vector<double>> v1_flatten, v2_flatten;
    multiFlatten(v1, v1_flatten);
    multiFlatten(v2, v2_flatten);
    multiSumVector(v1_flatten, v2_flatten, sumV);
  }

  void MergeTreePrincipalGeodesicsBase::flatten(
    std::vector<std::vector<double>> &v, std::vector<double> &newV) {
    newV.resize(v.size() * v[0].size());
    for(unsigned int i = 0; i < v.size(); ++i)
      for(unsigned int j = 0; j < v[0].size(); ++j)
        newV[i * v[0].size() + j] = v[i][j];
  }

  void MergeTreePrincipalGeodesicsBase::multiFlatten(
    std::vector<std::vector<std::vector<double>>> &v,
    std::vector<std::vector<double>> &newV) {
    newV.resize(v.size());
    for(unsigned int i = 0; i < v.size(); ++i)
      flatten(v[i], newV[i]);
  }

  void MergeTreePrincipalGeodesicsBase::unflatten(
    std::vector<double> &v, std::vector<std::vector<double>> &newV) {
    newV.resize(v.size() / 2);
    for(unsigned int i = 0; i < v.size(); i += 2)
      newV[i / 2] = std::vector<double>{v[i], v[i + 1]};
  }

  bool
    MergeTreePrincipalGeodesicsBase::isVectorUniform(std::vector<double> &v) {
    for(unsigned int i = 0; i < v.size() - 1; ++i)
      if(not(std::abs(v[i] - v[i + 1]) < 1e-6))
        return false;
    return true;
  }

  bool MergeTreePrincipalGeodesicsBase::isVectorNull(std::vector<double> &v) {
    for(unsigned int i = 0; i < v.size(); ++i)
      if(not(std::abs(v[i]) < 1e-6))
        return false;
    return true;
  }

  bool MergeTreePrincipalGeodesicsBase::isVectorNullFlatten(
    std::vector<std::vector<double>> &v) {
    std::vector<double> v_flat;
    flatten(v, v_flat);
    return isVectorNull(v_flat);
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
    std::vector<size_t> sizes(pVec.size(), size);
    pointersToVectors(pVec, sizes, vec);
  }

  // ----------------------------------------------------------------------------
  // Matrix Utils
  // ----------------------------------------------------------------------------
  void MergeTreePrincipalGeodesicsBase::matrixDot(
    std::vector<std::vector<double>> &m1,
    std::vector<std::vector<double>> &m2,
    std::vector<std::vector<double>> &newM) {
    newM.resize(m1.size(), std::vector<double>(m2[0].size(), 0.0));
    for(unsigned int i = 0; i < newM.size(); ++i)
      for(unsigned int j = 0; j < newM[i].size(); ++j)
        for(unsigned int k = 0; k < m1[i].size(); ++k)
          newM[i][j] += m1[i][k] * m2[k][j];
  }

  void MergeTreePrincipalGeodesicsBase::subMatrix(
    std::vector<std::vector<double>> &m1,
    std::vector<std::vector<double>> &m2,
    std::vector<std::vector<double>> &newM) {
    newM.resize(m1.size(), std::vector<double>(m1[0].size()));
    for(unsigned int i = 0; i < m1.size(); ++i)
      for(unsigned int j = 0; j < m1[0].size(); ++j)
        newM[i][j] = m1[i][j] - m2[i][j];
  }

  void MergeTreePrincipalGeodesicsBase::sumMatrix(
    std::vector<std::vector<double>> &m1,
    std::vector<std::vector<double>> &m2,
    std::vector<std::vector<double>> &newM) {
    newM.resize(m1.size(), std::vector<double>(m1[0].size()));
    for(unsigned int i = 0; i < m1.size(); ++i)
      for(unsigned int j = 0; j < m1[0].size(); ++j)
        newM[i][j] = m1[i][j] + m2[i][j];
  }

  void MergeTreePrincipalGeodesicsBase::multMatrix(
    std::vector<std::vector<double>> &m1,
    double mult,
    std::vector<std::vector<double>> &newM) {
    newM = m1;
    for(unsigned int i = 0; i < newM.size(); ++i)
      for(unsigned int j = 0; j < newM[i].size(); ++j)
        newM[i][j] *= mult;
  }

  void MergeTreePrincipalGeodesicsBase::transpose(
    std::vector<std::vector<double>> &m,
    std::vector<std::vector<double>> &newM) {
    std::vector<std::vector<double>> mTemp(
      m[0].size(), std::vector<double>(m.size()));
    for(unsigned int i = 0; i < m.size(); ++i)
      for(unsigned int j = 0; j < m[0].size(); ++j)
        mTemp[j][i] = m[i][j];
    newM = mTemp;
  }
} // namespace ttk
