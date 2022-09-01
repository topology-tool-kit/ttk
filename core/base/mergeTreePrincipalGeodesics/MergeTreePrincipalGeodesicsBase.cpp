#include <MergeTreePrincipalGeodesicsBase.h>

// other includes
#include <Geometry.h>
#include <vector>

namespace ttk {
  // ----------------------------------------------------------------------------
  // Vector Utils
  // ----------------------------------------------------------------------------
  double MergeTreePrincipalGeodesicsBase::distanceL2(std::vector<double> &v1,
                                                     std::vector<double> &v2) {
    return ttk::Geometry::distance<double>(v1.data(), v2.data(), v1.size());
  }

  double MergeTreePrincipalGeodesicsBase::distanceL2Flatten(
    std::vector<std::vector<double>> &v1,
    std::vector<std::vector<double>> &v2) {
    std::vector<double> v1_flatten, v2_flatten;
    flatten(v1, v1_flatten);
    flatten(v2, v2_flatten);
    return distanceL2(v1_flatten, v2_flatten);
  }

  double
    MergeTreePrincipalGeodesicsBase::scalarProduct(std::vector<double> &v1,
                                                   std::vector<double> &v2) {
    return ttk::Geometry::dotProduct(v1.data(), v2.data(), v1.size());
  }

  double MergeTreePrincipalGeodesicsBase::scalarProductFlatten(
    std::vector<std::vector<double>> &v1,
    std::vector<std::vector<double>> &v2) {
    std::vector<double> v1_flatten, v2_flatten;
    flatten(v1, v1_flatten);
    flatten(v2, v2_flatten);
    return scalarProduct(v1_flatten, v2_flatten);
  }

  double MergeTreePrincipalGeodesicsBase::norm(std::vector<double> &v) {
    return ttk::Geometry::magnitude(v.data(), v.size());
  }

  double MergeTreePrincipalGeodesicsBase::normFlatten(
    std::vector<std::vector<double>> &v) {
    std::vector<double> v_flatten;
    flatten(v, v_flatten);
    return norm(v_flatten);
  }

  void MergeTreePrincipalGeodesicsBase::vectorProjection(
    std::vector<double> &v1,
    std::vector<double> &v2,
    std::vector<double> &projec) {
    projec = std::vector<double>(v1.size(), 0.0);
    ttk::Geometry::vectorProjection(
      v2.data(), v1.data(), projec.data(), v1.size());
  }

  void
    MergeTreePrincipalGeodesicsBase::sumProjection(std::vector<double> &v1,
                                                   std::vector<double> &v2,
                                                   std::vector<double> &v1Out,
                                                   std::vector<double> &v2Out) {
    std::vector<double> sumV;
    sumVector(v1, v2, sumV);
    vectorProjection(sumV, v1, v1Out);
    vectorProjection(sumV, v2, v2Out);
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
      vectorProjection(uS[0], allVs[i], projecSum);
      for(unsigned int j = 1; j < i; ++j) {
        std::vector<double> projecTemp, projecSumTemp;
        vectorProjection(uS[j], allVs[i], projecTemp);
        sumVector(projecSum, projecTemp, projecSumTemp);
        projecSum = projecSumTemp;
      }
      subVector(allVs[i], projecSum, uS[i]);
    }
    newV = uS[uS.size() - 1];
  }

  void MergeTreePrincipalGeodesicsBase::sumVector(std::vector<double> &v1,
                                                  std::vector<double> &v2,
                                                  std::vector<double> &sumV) {
    if(sumV.size() != v1.size())
      sumV = std::vector<double>(v1.size());
    ttk::Geometry::addVectors(v1.data(), v2.data(), sumV.data(), v1.size());
  }

  void MergeTreePrincipalGeodesicsBase::subVector(std::vector<double> &v1,
                                                  std::vector<double> &v2,
                                                  std::vector<double> &subV) {
    subV = std::vector<double>(v1.size());
    ttk::Geometry::subtractVectors(
      v2.data(), v1.data(), subV.data(), v1.size());
  }

  void MergeTreePrincipalGeodesicsBase::multVectorByScalar(
    std::vector<double> &v, double scalar, std::vector<double> &multV) {
    multV = v;
    ttk::Geometry::scaleVector(v.data(), scalar, multV.data(), v.size());
  }

  void MergeTreePrincipalGeodesicsBase::multVectorByScalarFlatten(
    std::vector<std::vector<double>> &v,
    double scalar,
    std::vector<std::vector<double>> &multV) {
    std::vector<double> v_flat, multV_flat;
    flatten(v, v_flat);
    multVectorByScalar(v_flat, scalar, multV_flat);
    unflatten(multV_flat, multV);
  }

  double
    MergeTreePrincipalGeodesicsBase::sumVectorElements(std::vector<double> &v) {
    double sum = 0;
    for(auto e : v)
      sum += e;
    return sum;
  }

  double MergeTreePrincipalGeodesicsBase::sumVectorElementsFlatten(
    std::vector<std::vector<double>> &v) {
    std::vector<double> v_flatten;
    flatten(v, v_flatten);
    return sumVectorElements(v_flatten);
  }

  void MergeTreePrincipalGeodesicsBase::multiSumVector(
    std::vector<std::vector<double>> &v1,
    std::vector<std::vector<double>> &v2,
    std::vector<std::vector<double>> &sumV) {
    sumV = std::vector<std::vector<double>>(v1.size());
    for(unsigned int i = 0; i < v1.size(); ++i)
      sumVector(v1[i], v2[i], sumV[i]);
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
    newV = std::vector<double>(v.size() * v[0].size());
    for(unsigned int i = 0; i < v.size(); ++i)
      for(unsigned int j = 0; j < v[0].size(); ++j)
        newV[i * v[0].size() + j] = v[i][j];
  }

  void MergeTreePrincipalGeodesicsBase::multiFlatten(
    std::vector<std::vector<std::vector<double>>> &v,
    std::vector<std::vector<double>> &newV) {
    newV = std::vector<std::vector<double>>(v.size());
    for(unsigned int i = 0; i < v.size(); ++i)
      flatten(v[i], newV[i]);
  }

  void MergeTreePrincipalGeodesicsBase::unflatten(
    std::vector<double> &v, std::vector<std::vector<double>> &newV) {
    newV = std::vector<std::vector<double>>(v.size() / 2);
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
    pVec = std::vector<double *>(vec.size());
    for(unsigned int i = 0; i < vec.size(); ++i)
      vectorToPointer(vec[i], pVec[i]);
  }

  void MergeTreePrincipalGeodesicsBase::vectorOfVectorsToPointers(
    std::vector<std::vector<std::vector<double>>> &vS,
    std::vector<std::vector<double *>> &pVS) {
    pVS = std::vector<std::vector<double *>>(vS.size());
    for(unsigned int i = 0; i < vS.size(); ++i)
      vectorsToPointers(vS[i], pVS[i]);
  }

  void MergeTreePrincipalGeodesicsBase::pointerToVector(
    double *pVec, size_t size, std::vector<double> &vec) {
    vec = std::vector<double>(size);
    for(unsigned int i = 0; i < size; ++i)
      vec[i] = pVec[i];
  }

  void MergeTreePrincipalGeodesicsBase::pointersToVectors(
    std::vector<double *> &pVec,
    std::vector<size_t> sizes,
    std::vector<std::vector<double>> &vec) {
    vec = std::vector<std::vector<double>>(pVec.size());
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
    newM = std::vector<std::vector<double>>(
      m1.size(), std::vector<double>(m2[0].size(), 0.0));
    for(unsigned int i = 0; i < newM.size(); ++i)
      for(unsigned int j = 0; j < newM[i].size(); ++j)
        for(unsigned int k = 0; k < m1[i].size(); ++k)
          newM[i][j] += m1[i][k] * m2[k][j];
  }

  void MergeTreePrincipalGeodesicsBase::subMatrix(
    std::vector<std::vector<double>> &m1,
    std::vector<std::vector<double>> &m2,
    std::vector<std::vector<double>> &newM) {
    newM = std::vector<std::vector<double>>(
      m1.size(), std::vector<double>(m1[0].size()));
    for(unsigned int i = 0; i < m1.size(); ++i)
      for(unsigned int j = 0; j < m1[0].size(); ++j)
        newM[i][j] = m1[i][j] - m2[i][j];
  }

  void MergeTreePrincipalGeodesicsBase::sumMatrix(
    std::vector<std::vector<double>> &m1,
    std::vector<std::vector<double>> &m2,
    std::vector<std::vector<double>> &newM) {
    newM = std::vector<std::vector<double>>(
      m1.size(), std::vector<double>(m1[0].size()));
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

  // ----------------------------------------------------------------------------
  // Statistics Utils
  // ----------------------------------------------------------------------------
  double MergeTreePrincipalGeodesicsBase::mean(std::vector<double> &v) {
    double mean = 0.0;
    for(auto e : v)
      mean += e;
    return mean / v.size();
  }

  double MergeTreePrincipalGeodesicsBase::var(std::vector<double> &v) {
    return cov(v, v);
  }

  double MergeTreePrincipalGeodesicsBase::cov(std::vector<double> &v1,
                                              std::vector<double> &v2) {
    double cov = 0.0;
    double meanV1 = mean(v1);
    double meanV2 = mean(v2);
    for(unsigned int i = 0; i < v1.size(); ++i)
      cov += (v1[i] - meanV1) * (v2[i] - meanV2);
    return cov / (v1.size() - 1);
  }

  double MergeTreePrincipalGeodesicsBase::corr(std::vector<double> &v1,
                                               std::vector<double> &v2) {
    return cov(v1, v2) / (std::sqrt(var(v1)) * std::sqrt(var(v2)));
  }
} // namespace ttk
