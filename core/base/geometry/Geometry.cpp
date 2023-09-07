#include <Geometry.h>

#include <algorithm>

using namespace std;
using namespace ttk;

static const double PREC_DBL{Geometry::pow(10.0, -DBL_DIG)};
static const float PREC_FLT{powf(10.0F, -FLT_DIG)};
static const float PREC_FLT_1{powf(10.0F, -FLT_DIG + 1)};

template <typename T>
T Geometry::angle(const T *vA0, const T *vA1, const T *vB0, const T *vB1) {
  return M_PI
         - acos(dotProduct(vA0, vA1, vB0, vB1)
                / (magnitude(vA0, vA1) * magnitude(vB0, vB1)));
}

template <typename T>
bool Geometry::areVectorsColinear(const T *vA0,
                                  const T *vA1,
                                  const T *vB0,
                                  const T *vB1,
                                  std::array<T, 3> *coefficients,
                                  const T *tolerance) {

  int aNullComponents = 0, bNullComponents = 0;
  std::array<T, 3> a{}, b{};
  for(int i = 0; i < 3; i++) {
    a[i] = vA1[i] - vA0[i];
    if(fabs(a[i]) < PREC_FLT) {
      aNullComponents++;
    }
    b[i] = vB1[i] - vB0[i];
    if(fabs(b[i]) < PREC_FLT) {
      bNullComponents++;
    }
  }

  if((aNullComponents == 3) || (bNullComponents == 3)) {
    return true;
  }

  // check for axis aligned vectors
  if((aNullComponents > 1) || (bNullComponents > 1)) {
    if(aNullComponents == bNullComponents) {
      // only one non-null component for both vectors
      return true;
    }
  }

  bool useDenominatorA = false;
  T sumA = 0, sumB = 0;
  for(int i = 0; i < 3; i++) {
    sumA += fabs(a[i]);
    sumB += fabs(b[i]);
  }
  if(sumA > sumB) {
    useDenominatorA = true;
  }

  std::array<T, 3> k{};

  T maxDenominator = 0;
  int isNan = -1, maximizer = 0;
  for(int i = 0; i < 3; i++) {
    if(useDenominatorA) {
      if(fabs(a[i]) > PREC_FLT) {
        k[i] = b[i] / a[i];
      } else {
        isNan = i;
      }
    } else {
      if(fabs(b[i]) > PREC_FLT) {
        k[i] = a[i] / b[i];
      } else {
        isNan = i;
      }
    }

    if(!i) {
      maxDenominator = fabs(k[i]);
      maximizer = i;
    } else {
      if(fabs(k[i]) > maxDenominator) {
        maxDenominator = fabs(k[i]);
        maximizer = i;
      }
    }
  }

  T colinearityThreshold;

  colinearityThreshold = PREC_FLT;
  if(tolerance) {
    colinearityThreshold = *tolerance;
  }

  if(coefficients) {
    (*coefficients) = k;
  }

  if(isNan == -1) {

    if((fabs(1 - fabs(k[(maximizer + 1) % 3] / k[maximizer]))
        < colinearityThreshold)
       && (fabs(1 - fabs(k[(maximizer + 2) % 3] / k[maximizer]))
           < colinearityThreshold)) {
      return true;
    }
  } else {
    if(fabs(1 - fabs(k[(isNan + 1) % 3] / k[(isNan + 2) % 3]))
       < colinearityThreshold) {
      return true;
    }
  }

  k[0] = k[1] = k[2] = 0;

  return false;
}

template <typename T>
int Geometry::computeBarycentricCoordinates(const T *p0,
                                            const T *p1,
                                            const T *p,
                                            std::array<T, 2> &baryCentrics,
                                            const int &dimension) {

  if(dimension > 3) {
    return -1;
  }

  int bestI = 0;
  T maxDenominator = 0;

  for(int i = 0; i < dimension; i++) {

    T denominator = fabs(p0[i] - p1[i]);
    if(!i) {
      maxDenominator = denominator;
      bestI = i;
    } else {
      if(denominator > maxDenominator) {
        maxDenominator = denominator;
        bestI = i;
      }
    }
  }

  baryCentrics[0] = p0[bestI] - p1[bestI];
  baryCentrics[0] = (p[bestI] - p1[bestI]) / baryCentrics[0];

  baryCentrics[1] = 1 - baryCentrics[0];

  // check if the point lies in the edge
  std::array<T, 3> test{};
  for(int i = 0; i < dimension; i++) {
    test[i] = baryCentrics[0] * p0[i] + baryCentrics[1] * p1[i];
  }

  if((!((fabs(test[0] - p[0]) < PREC_FLT_1)
        && (fabs(test[1] - p[1]) < PREC_FLT_1)))) {
    for(int i = 0; i < 2; i++) {
      baryCentrics[i] = -baryCentrics[i];
    }
  }

  return 0;
}
template <typename T>
int Geometry::computeBarycentricCoordinates(const T *p0,
                                            const T *p1,
                                            const T *p2,
                                            const T *p,
                                            std::array<T, 3> &baryCentrics) {

  // find the pair of coordinates that maximize the sum of the denominators
  // (more stable computations)
  int bestI = 0, bestJ = 1;
  T maxDenominator = 0;

  for(int i = 0; i < 2; i++) {
    for(int j = i + 1; j < 3; j++) {

      baryCentrics[0]
        = (p1[j] - p2[j]) * (p0[i] - p2[i]) + (p2[i] - p1[i]) * (p0[j] - p2[j]);
      baryCentrics[1]
        = (p1[j] - p2[j]) * (p0[i] - p2[i]) + (p2[i] - p1[i]) * (p0[j] - p2[j]);

      T denominator = fabs(baryCentrics[0]);

      if(fabs(baryCentrics[1]) < denominator) {
        denominator = fabs(baryCentrics[1]);
      }

      if((i == 0) && (j == 1)) {
        maxDenominator = denominator;
      } else {
        if(denominator > maxDenominator) {
          maxDenominator = denominator;
          bestI = i;
          bestJ = j;
        }
      }
    }
  }

  baryCentrics[0] = (p1[bestJ] - p2[bestJ]) * (p0[bestI] - p2[bestI])
                    + (p2[bestI] - p1[bestI]) * (p0[bestJ] - p2[bestJ]);
  // (y1 - y2)*(x - x2) + (x2 - x1)*(y - y2)
  baryCentrics[0] = ((p1[bestJ] - p2[bestJ]) * (p[bestI] - p2[bestI])
                     + (p2[bestI] - p1[bestI]) * (p[bestJ] - p2[bestJ]))
                    / baryCentrics[0];

  // (y1 - y2)*(x0 - x2) + (x2 - x1)*(y0 - y2)
  baryCentrics[1] = (p1[bestJ] - p2[bestJ]) * (p0[bestI] - p2[bestI])
                    + (p2[bestI] - p1[bestI]) * (p0[bestJ] - p2[bestJ]);
  // (y2 - y0)*(x - x2) + (x0 - x2)*(y - y2)
  baryCentrics[1] = ((p2[bestJ] - p0[bestJ]) * (p[bestI] - p2[bestI])
                     + (p0[bestI] - p2[bestI]) * (p[bestJ] - p2[bestJ]))
                    / baryCentrics[1];

  baryCentrics[2] = 1 - baryCentrics[0] - baryCentrics[1];

  return 0;
}

template <typename T>
bool Geometry::computeSegmentIntersection(const T &xA,
                                          const T &yA,
                                          const T &xB,
                                          const T &yB,
                                          const T &xC,
                                          const T &yC,
                                          const T &xD,
                                          const T &yD,
                                          T &x,
                                          T &y) {

  T d = (xA - xB) * (yC - yD) - (yA - yB) * (xC - xD);

  if(fabs(d) < PREC_DBL) {
    return false;
  }

  x = ((xC - xD) * (xA * yB - yA * xB) - (xA - xB) * (xC * yD - yC * xD)) / d;

  y = ((yC - yD) * (xA * yB - yA * xB) - (yA - yB) * (xC * yD - yC * xD)) / d;

  if((x < std::min(xA, xB) - PREC_FLT) || (x > std::max(xA, xB) + PREC_FLT)) {
    return false;
  }

  if((x < std::min(xC, xD) - PREC_FLT) || (x > std::max(xC, xD) + PREC_FLT)) {
    return false;
  }

  return true;
}

template <typename T>
int Geometry::computeTriangleArea(const T *p0,
                                  const T *p1,
                                  const T *p2,
                                  T &area) {

  std::array<T, 3> cross{};

  crossProduct(p0, p1, p1, p2, cross);

  area = 0.5 * magnitude(cross.data());

  return 0;
}

template <typename T>
int Geometry::computeTriangleAreaFromSides(const T s0,
                                           const T s1,
                                           const T s2,
                                           T &area) {

  double const s = (s0 + s1 + s2) / 2.0;
  area = std::sqrt(s * (s - s0) * (s - s1) * (s - s2));

  return 0;
}

template <typename T>
int Geometry::computeTriangleAngles(const T *p0,
                                    const T *p1,
                                    const T *p2,
                                    std::array<T, 3> &angles) {

  angles[0] = angle(p0, p1, p1, p2);
  angles[1] = angle(p1, p2, p2, p0);
  angles[2] = angle(p2, p0, p0, p1);

  return 0;
}

template <typename T>
int Geometry::computeTriangleAngleFromSides(const T s0,
                                            const T s1,
                                            const T s2,
                                            T &angle) {

  angle = std::acos((s0 * s0 + s1 * s1 - s2 * s2) / (2.0 * s0 * s1));

  return 0;
}

template <typename T>
int Geometry::crossProduct(const T *vA0,
                           const T *vA1,
                           const T *vB0,
                           const T *vB1,
                           std::array<T, 3> &crossProduct) {

  std::array<T, 3> a{}, b{};

  for(int i = 0; i < 3; i++) {
    a[i] = vA1[i] - vA0[i];
    b[i] = vB1[i] - vB0[i];
  }

  for(int i = 0; i < 3; i++) {
    crossProduct[i]
      = a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3];
  }

  return 0;
}

template <typename T>
int Geometry::crossProduct(const T *vA, const T *vB, T *crossProduct) {
  crossProduct[0] = vA[1] * vB[2] - vA[2] * vB[1];
  crossProduct[1] = vA[2] * vB[0] - vA[0] * vB[2];
  crossProduct[2] = vA[0] * vB[1] - vA[1] * vB[0];
  return 0;
}

template <typename T>
T Geometry::distance(const T *p0, const T *p1, const int &dimension) {

  T distance = 0;

  for(int i = 0; i < dimension; i++) {
    distance += (p0[i] - p1[i]) * (p0[i] - p1[i]);
  }

  return sqrt(distance);
}

template <typename T>
T Geometry::distance(const std::vector<T> &p0, const std::vector<T> &p1) {
  return distance(p0.data(), p1.data(), p0.size());
}

template <typename T>
T Geometry::distanceFlatten(const std::vector<std::vector<T>> &p0,
                            const std::vector<std::vector<T>> &p1) {
  std::vector<T> p0_flatten, p1_flatten;
  flattenMultiDimensionalVector(p0, p0_flatten);
  flattenMultiDimensionalVector(p1, p1_flatten);
  return distance(p0_flatten, p1_flatten);
}

template <typename T>
T Geometry::dotProduct(const T *vA0, const T *vA1, const T *vB0, const T *vB1) {

  T dotProduct = 0;
  for(int i = 0; i < 3; i++) {
    dotProduct += (vA1[i] - vA0[i]) * (vB1[i] - vB0[i]);
  }

  return dotProduct;
}

template <typename T>
T Geometry::dotProduct(const T *vA, const T *vB, const int &dimension) {
  T dotProduct = 0;
  for(int i = 0; i < dimension; ++i)
    dotProduct += vA[i] * vB[i];
  return dotProduct;
}

template <typename T>
T Geometry::dotProduct(const std::vector<T> &vA, const std::vector<T> &vB) {
  return dotProduct(vA.data(), vB.data(), vA.size());
}

template <typename T>
T Geometry::dotProductFlatten(const std::vector<std::vector<T>> &vA,
                              const std::vector<std::vector<T>> &vB) {
  std::vector<T> vA_flatten, vB_flatten;
  flattenMultiDimensionalVector(vA, vA_flatten);
  flattenMultiDimensionalVector(vB, vB_flatten);
  return dotProduct(vA_flatten, vB_flatten);
}

template <typename T>
bool Geometry::isPointInTriangle(const T *p0,
                                 const T *p1,
                                 const T *p2,
                                 const T *p) {

  std::array<T, 3> barycentrics{};

  Geometry::computeBarycentricCoordinates(p0, p1, p2, p, barycentrics);

  for(int i = 0; i < static_cast<int>(barycentrics.size()); i++) {
    if(barycentrics[i] < -PREC_DBL) {
      return false;
    }
    if(barycentrics[i] > 1 + PREC_DBL) {
      return false;
    }
  }

  return true;
}

template <typename T>
bool Geometry::isPointOnSegment(
  const T &x, const T &y, const T &xA, const T &yA, const T &xB, const T &yB) {

  std::array<T, 2> pA{xA, yA}, pB{xB, yB}, p{x, y};
  return Geometry::isPointOnSegment(p.data(), pA.data(), pB.data(), 2);
}

template <typename T>
bool Geometry::isPointOnSegment(const T *p,
                                const T *pA,
                                const T *pB,
                                const int &dimension) {

  std::array<T, 2> baryCentrics{};

  Geometry::computeBarycentricCoordinates(pA, pB, p, baryCentrics, dimension);

  return (
    ((baryCentrics[0] > -PREC_DBL) && (baryCentrics[0] < 1 + PREC_DBL))
    && ((baryCentrics[1] > -PREC_DBL) && (baryCentrics[1] < 1 + PREC_DBL)));
}

template <typename T>
bool Geometry::isTriangleColinear(const T *p0,
                                  const T *p1,
                                  const T *p2,
                                  const T *tolerance) {

  bool maxDecision = false;
  T maxCoefficient = 0;
  std::array<T, 3> coefficients{};

  bool decision = areVectorsColinear(p0, p1, p1, p2, &coefficients, tolerance);
  maxDecision = decision;
  for(int i = 0; i < 3; i++) {
    if(!i) {
      maxCoefficient = fabs(coefficients[i]);
      maxDecision = decision;
    } else {
      if(fabs(coefficients[i]) > maxCoefficient) {
        maxCoefficient = fabs(coefficients[i]);
        maxDecision = decision;
      }
    }
  }

  decision = areVectorsColinear(p0, p2, p2, p1, &coefficients, tolerance);
  for(int i = 0; i < 3; i++) {
    if(fabs(coefficients[i]) > maxCoefficient) {
      maxCoefficient = fabs(coefficients[i]);
      maxDecision = decision;
    }
  }

  decision = areVectorsColinear(p1, p0, p0, p2, &coefficients, tolerance);
  for(int i = 0; i < 3; i++) {
    if(fabs(coefficients[i]) > maxCoefficient) {
      maxCoefficient = fabs(coefficients[i]);
      maxDecision = decision;
    }
  }

  return maxDecision;
}

template <typename T>
T Geometry::magnitude(const T *v, const int &dimension) {
  return sqrt(dotProduct(v, v, dimension));
}

template <typename T>
T Geometry::magnitude(const std::vector<T> &v) {
  return magnitude(v.data(), v.size());
}

template <typename T>
T Geometry::magnitudeFlatten(const std::vector<std::vector<T>> &v) {
  std::vector<T> v_flatten;
  flattenMultiDimensionalVector(v, v_flatten);
  return magnitude(v_flatten);
}

template <typename T>
T Geometry::magnitude(const T *o, const T *d) {

  T mag = 0;

  for(int i = 0; i < 3; i++) {
    mag += (o[i] - d[i]) * (o[i] - d[i]);
  }

  return sqrt(mag);
}

template <typename T>
int Geometry::subtractVectors(const T *a,
                              const T *b,
                              T *out,
                              const int &dimension) {
  for(int i = 0; i < dimension; ++i)
    out[i] = b[i] - a[i];
  return 0;
}

template <typename T>
int Geometry::subtractVectors(const std::vector<T> &a,
                              const std::vector<T> &b,
                              std::vector<T> &out) {
  out.resize(a.size());
  return subtractVectors(a.data(), b.data(), out.data(), a.size());
}

template <typename T>
int Geometry::addVectors(const T *a, const T *b, T *out, const int &dimension) {
  for(int i = 0; i < dimension; ++i)
    out[i] = b[i] + a[i];
  return 0;
}

template <typename T>
int Geometry::addVectors(const std::vector<T> &a,
                         const std::vector<T> &b,
                         std::vector<T> &out) {
  out.resize(a.size());
  return addVectors(a.data(), b.data(), out.data(), a.size());
}

template <typename T>
int Geometry::multiAddVectors(const std::vector<std::vector<T>> &a,
                              const std::vector<std::vector<T>> &b,
                              std::vector<std::vector<T>> &out) {
  out.resize(a.size());
  for(unsigned int i = 0; i < a.size(); ++i)
    addVectors(a[i], b[i], out[i]);
  return 0;
}

template <typename T>
int Geometry::multiAddVectorsFlatten(
  const std::vector<std::vector<std::vector<T>>> &a,
  const std::vector<std::vector<std::vector<T>>> &b,
  std::vector<std::vector<T>> &out) {
  std::vector<std::vector<T>> a_flatten, b_flatten;
  multiFlattenMultiDimensionalVector(a, a_flatten);
  multiFlattenMultiDimensionalVector(b, b_flatten);
  multiAddVectors(a_flatten, b_flatten, out);
  return 0;
}

template <typename T>
int Geometry::scaleVector(const T *a,
                          const T factor,
                          T *out,
                          const int &dimension) {
  for(int i = 0; i < dimension; ++i)
    out[i] = a[i] * factor;
  return 0;
}

template <typename T>
int Geometry::scaleVector(const std::vector<T> &a,
                          const T factor,
                          std::vector<T> &out) {
  out.resize(a.size());
  return scaleVector(a.data(), factor, out.data(), a.size());
}

template <typename T>
int Geometry::vectorProjection(const T *a,
                               const T *b,
                               T *out,
                               const int &dimension) {
  T dotProdBB = dotProduct(b, b, dimension);
  T dotProdAB;
  if(dotProdBB > PREC_DBL) {
    dotProdAB = dotProduct(a, b, dimension);
    dotProdAB /= dotProdBB;
  } else
    dotProdAB = 0; // Gram-Schmidt convention
  for(int i = 0; i < dimension; ++i)
    out[i] = b[i] * dotProdAB;
  return 0;
}

template <typename T>
int Geometry::vectorProjection(const std::vector<T> &a,
                               const std::vector<T> &b,
                               std::vector<T> &out) {
  out.resize(a.size(), 0.0);
  return vectorProjection(a.data(), b.data(), out.data(), a.size());
}

template <typename T>
void Geometry::addVectorsProjection(const std::vector<T> &a,
                                    const std::vector<T> &b,
                                    std::vector<T> &a_out,
                                    std::vector<T> &b_out) {
  std::vector<T> sumV;
  addVectors(a, b, sumV);
  vectorProjection(a, sumV, a_out);
  vectorProjection(b, sumV, b_out);
}

template <typename T>
void Geometry::gramSchmidt(const std::vector<std::vector<T>> &a,
                           std::vector<std::vector<T>> &out) {
  out.resize(a.size());
  out[0] = a[0];
  for(unsigned int i = 1; i < a.size(); ++i) {
    std::vector<T> projecSum;
    vectorProjection(a[i], out[0], projecSum);
    for(unsigned int j = 1; j < i; ++j) {
      std::vector<T> projecTemp, projecSumTemp;
      vectorProjection(a[i], out[j], projecTemp);
      addVectors(projecSum, projecTemp, projecSumTemp);
      projecSum = projecSumTemp;
    }
    subtractVectors(projecSum, a[i], out[i]);
  }
}

template <typename T>
bool Geometry::isVectorUniform(const std::vector<T> &a) {
  for(unsigned int i = 0; i < a.size() - 1; ++i)
    if(not(std::abs(a[i] - a[i + 1]) < PREC_DBL))
      return false;
  return true;
}

template <typename T>
bool Geometry::isVectorNull(const std::vector<T> &a) {
  for(unsigned int i = 0; i < a.size(); ++i)
    if(not(std::abs(a[i]) < PREC_DBL))
      return false;
  return true;
}

template <typename T>
bool Geometry::isVectorNullFlatten(const std::vector<std::vector<T>> &a) {
  std::vector<T> a_flatten;
  flattenMultiDimensionalVector(a, a_flatten);
  return isVectorNull(a_flatten);
}

template <typename T>
int Geometry::flattenMultiDimensionalVector(
  const std::vector<std::vector<T>> &a, std::vector<T> &out) {
  out.resize(a.size() * a[0].size());
  for(unsigned int i = 0; i < a.size(); ++i)
    for(unsigned int j = 0; j < a[0].size(); ++j)
      out[i * a[0].size() + j] = a[i][j];
  return 0;
}

template <typename T>
int Geometry::multiFlattenMultiDimensionalVector(
  const std::vector<std::vector<std::vector<T>>> &a,
  std::vector<std::vector<T>> &out) {
  out.resize(a.size());
  for(unsigned int i = 0; i < a.size(); ++i)
    flattenMultiDimensionalVector(a[i], out[i]);
  return 0;
}

template <typename T>
int Geometry::unflattenMultiDimensionalVector(const std::vector<T> &a,
                                              std::vector<std::vector<T>> &out,
                                              const int &no_columns) {
  if(a.size() % no_columns != 0)
    return -1;
  out.resize(a.size() / no_columns);
  for(unsigned int i = 0; i < out.size(); ++i) {
    out[i].resize(no_columns);
    for(unsigned int j = 0; j < out[i].size(); ++j)
      out[i][j] = a[i * no_columns + j];
  }
  return 0;
}

template <typename T>
void Geometry::matrixMultiplication(const std::vector<std::vector<T>> &a,
                                    const std::vector<std::vector<T>> &b,
                                    std::vector<std::vector<T>> &out) {
  out.resize(a.size(), std::vector<T>(b[0].size(), 0.0));
  for(unsigned int i = 0; i < out.size(); ++i)
    for(unsigned int j = 0; j < out[i].size(); ++j)
      for(unsigned int k = 0; k < a[i].size(); ++k)
        out[i][j] += a[i][k] * b[k][j];
}

template <typename T>
void Geometry::subtractMatrices(const std::vector<std::vector<T>> &a,
                                const std::vector<std::vector<T>> &b,
                                std::vector<std::vector<T>> &out) {
  out.resize(a.size(), std::vector<T>(a[0].size()));
  for(unsigned int i = 0; i < out.size(); ++i)
    for(unsigned int j = 0; j < out[0].size(); ++j)
      out[i][j] = b[i][j] - a[i][j];
}

template <typename T>
void Geometry::addMatrices(const std::vector<std::vector<T>> &a,
                           const std::vector<std::vector<T>> &b,
                           std::vector<std::vector<T>> &out) {
  out.resize(a.size(), std::vector<T>(a[0].size()));
  for(unsigned int i = 0; i < a.size(); ++i)
    for(unsigned int j = 0; j < a[0].size(); ++j)
      out[i][j] = a[i][j] + b[i][j];
}

template <typename T>
void Geometry::scaleMatrix(const std::vector<std::vector<T>> &a,
                           const T factor,
                           std::vector<std::vector<T>> &out) {
  out.resize(a.size(), std::vector<T>(a[0].size()));
  for(unsigned int i = 0; i < out.size(); ++i)
    for(unsigned int j = 0; j < out[i].size(); ++j)
      out[i][j] = a[i][j] * factor;
}

template <typename T>
void Geometry::transposeMatrix(const std::vector<std::vector<T>> &a,
                               std::vector<std::vector<T>> &out) {
  out.resize(a[0].size(), std::vector<T>(a.size()));
  for(unsigned int i = 0; i < a.size(); ++i)
    for(unsigned int j = 0; j < a[0].size(); ++j)
      out[j][i] = a[i][j];
}

#define GEOMETRY_SPECIALIZE(TYPE)                                              \
  template TYPE Geometry::angle<TYPE>(                                         \
    TYPE const *, TYPE const *, TYPE const *, TYPE const *);                   \
  template bool Geometry::areVectorsColinear<TYPE>(                            \
    TYPE const *, TYPE const *, TYPE const *, TYPE const *,                    \
    std::array<TYPE, 3> *, TYPE const *);                                      \
  template int Geometry::computeBarycentricCoordinates<TYPE>(                  \
    TYPE const *, TYPE const *, TYPE const *, std::array<TYPE, 2> &,           \
    int const &);                                                              \
  template int Geometry::computeBarycentricCoordinates<TYPE>(                  \
    TYPE const *, TYPE const *, TYPE const *, TYPE const *,                    \
    std::array<TYPE, 3> &);                                                    \
  template bool Geometry::computeSegmentIntersection<TYPE>(                    \
    TYPE const &, TYPE const &, TYPE const &, TYPE const &, TYPE const &,      \
    TYPE const &, TYPE const &, TYPE const &, TYPE &, TYPE &);                 \
  template int Geometry::computeTriangleAngles<TYPE>(                          \
    TYPE const *, TYPE const *, TYPE const *, std::array<TYPE, 3> &);          \
  template int Geometry::computeTriangleAngleFromSides<TYPE>(                  \
    TYPE const, TYPE const, TYPE const, TYPE &);                               \
  template int Geometry::computeTriangleArea<TYPE>(                            \
    TYPE const *, TYPE const *, TYPE const *, TYPE &);                         \
  template int Geometry::computeTriangleAreaFromSides<TYPE>(                   \
    TYPE const, TYPE const, TYPE const, TYPE &);                               \
  template int Geometry::crossProduct<TYPE>(TYPE const *, TYPE const *,        \
                                            TYPE const *, TYPE const *,        \
                                            std::array<TYPE, 3> &);            \
  template int Geometry::crossProduct<TYPE>(                                   \
    TYPE const *, TYPE const *, TYPE *);                                       \
  template TYPE Geometry::distance<TYPE>(                                      \
    TYPE const *, TYPE const *, int const &);                                  \
  template TYPE Geometry::distance<TYPE>(                                      \
    std::vector<TYPE> const &, std::vector<TYPE> const &);                     \
  template TYPE Geometry::distanceFlatten<TYPE>(                               \
    std::vector<std::vector<TYPE>> const &,                                    \
    std::vector<std::vector<TYPE>> const &);                                   \
  template TYPE Geometry::dotProduct<TYPE>(                                    \
    TYPE const *, TYPE const *, TYPE const *, TYPE const *);                   \
  template TYPE Geometry::dotProduct<TYPE>(                                    \
    TYPE const *, TYPE const *, int const &);                                  \
  template TYPE Geometry::dotProduct<TYPE>(                                    \
    std::vector<TYPE> const &, std::vector<TYPE> const &);                     \
  template TYPE Geometry::dotProductFlatten<TYPE>(                             \
    std::vector<std::vector<TYPE>> const &,                                    \
    std::vector<std::vector<TYPE>> const &);                                   \
  template bool Geometry::isPointInTriangle<TYPE>(                             \
    TYPE const *, TYPE const *, TYPE const *, TYPE const *);                   \
  template bool Geometry::isPointOnSegment<TYPE>(TYPE const &, TYPE const &,   \
                                                 TYPE const &, TYPE const &,   \
                                                 TYPE const &, TYPE const &);  \
  template bool Geometry::isPointOnSegment<TYPE>(                              \
    TYPE const *, TYPE const *, TYPE const *, int const &);                    \
  template bool Geometry::isTriangleColinear<TYPE>(                            \
    TYPE const *, TYPE const *, TYPE const *, TYPE const *);                   \
  template TYPE Geometry::magnitude<TYPE>(TYPE const *, int const &);          \
  template TYPE Geometry::magnitude<TYPE>(std::vector<TYPE> const &);          \
  template TYPE Geometry::magnitudeFlatten<TYPE>(                              \
    std::vector<std::vector<TYPE>> const &);                                   \
  template TYPE Geometry::magnitude<TYPE>(TYPE const *, TYPE const *);         \
  template int Geometry::subtractVectors<TYPE>(                                \
    TYPE const *, TYPE const *, TYPE *, int const &);                          \
  template int Geometry::subtractVectors<TYPE>(std::vector<TYPE> const &,      \
                                               std::vector<TYPE> const &,      \
                                               std::vector<TYPE> &);           \
  template int Geometry::addVectors<TYPE>(                                     \
    TYPE const *, TYPE const *, TYPE *, int const &);                          \
  template int Geometry::addVectors<TYPE>(std::vector<TYPE> const &,           \
                                          std::vector<TYPE> const &,           \
                                          std::vector<TYPE> &);                \
  template int Geometry::multiAddVectors<TYPE>(                                \
    std::vector<std::vector<TYPE>> const &,                                    \
    std::vector<std::vector<TYPE>> const &, std::vector<std::vector<TYPE>> &); \
  template int Geometry::multiAddVectorsFlatten<TYPE>(                         \
    std::vector<std::vector<std::vector<TYPE>>> const &,                       \
    std::vector<std::vector<std::vector<TYPE>>> const &,                       \
    std::vector<std::vector<TYPE>> &);                                         \
  template int Geometry::scaleVector<TYPE>(                                    \
    TYPE const *, TYPE const, TYPE *, int const &);                            \
  template int Geometry::scaleVector<TYPE>(                                    \
    std::vector<TYPE> const &, TYPE const, std::vector<TYPE> &);               \
  template int Geometry::vectorProjection<TYPE>(                               \
    TYPE const *, TYPE const *, TYPE *, int const &);                          \
  template int Geometry::vectorProjection<TYPE>(std::vector<TYPE> const &,     \
                                                std::vector<TYPE> const &,     \
                                                std::vector<TYPE> &);          \
  template void Geometry::addVectorsProjection<TYPE>(                          \
    std::vector<TYPE> const &, std::vector<TYPE> const &, std::vector<TYPE> &, \
    std::vector<TYPE> &);                                                      \
  template void Geometry::gramSchmidt<TYPE>(                                   \
    std::vector<std::vector<TYPE>> const &, std::vector<std::vector<TYPE>> &); \
  template bool Geometry::isVectorUniform<TYPE>(std::vector<TYPE> const &);    \
  template bool Geometry::isVectorNull<TYPE>(std::vector<TYPE> const &);       \
  template bool Geometry::isVectorNullFlatten<TYPE>(                           \
    std::vector<std::vector<TYPE>> const &);                                   \
  template int Geometry::flattenMultiDimensionalVector<TYPE>(                  \
    std::vector<std::vector<TYPE>> const &, std::vector<TYPE> &);              \
  template int Geometry::multiFlattenMultiDimensionalVector<TYPE>(             \
    std::vector<std::vector<std::vector<TYPE>>> const &,                       \
    std::vector<std::vector<TYPE>> &);                                         \
  template int Geometry::unflattenMultiDimensionalVector<TYPE>(                \
    std::vector<TYPE> const &, std::vector<std::vector<TYPE>> &, int const &); \
  template void Geometry::matrixMultiplication<TYPE>(                          \
    std::vector<std::vector<TYPE>> const &,                                    \
    std::vector<std::vector<TYPE>> const &, std::vector<std::vector<TYPE>> &); \
  template void Geometry::subtractMatrices<TYPE>(                              \
    std::vector<std::vector<TYPE>> const &,                                    \
    std::vector<std::vector<TYPE>> const &, std::vector<std::vector<TYPE>> &); \
  template void Geometry::addMatrices<TYPE>(                                   \
    std::vector<std::vector<TYPE>> const &,                                    \
    std::vector<std::vector<TYPE>> const &, std::vector<std::vector<TYPE>> &); \
  template void Geometry::scaleMatrix<TYPE>(                                   \
    std::vector<std::vector<TYPE>> const &, TYPE const,                        \
    std::vector<std::vector<TYPE>> &);                                         \
  template void Geometry::transposeMatrix<TYPE>(                               \
    std::vector<std::vector<TYPE>> const &, std::vector<std::vector<TYPE>> &);

// explicit specializations for float and double
GEOMETRY_SPECIALIZE(float);
GEOMETRY_SPECIALIZE(double);
