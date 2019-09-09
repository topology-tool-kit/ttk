#include <Geometry.h>

#include <algorithm>

// TODO:
// for faster computations, remove all pow10

using namespace std;
using namespace ttk;

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
                                  vector<T> *coefficients,
                                  const T *tolerance) {

  int aNullComponents = 0, bNullComponents = 0;
  vector<T> a(3), b(3);
  for(int i = 0; i < 3; i++) {
    a[i] = vA1[i] - vA0[i];
    if(fabs(a[i]) < pow10(-FLT_DIG)) {
      aNullComponents++;
    }
    b[i] = vB1[i] - vB0[i];
    if(fabs(b[i]) < pow10(-FLT_DIG)) {
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

  vector<T> k(3, 0);

  T maxDenominator = 0;
  int isNan = -1, maximizer = 0;
  for(int i = 0; i < 3; i++) {
    if(useDenominatorA) {
      if(fabs(a[i]) > pow10(-FLT_DIG)) {
        k[i] = b[i] / a[i];
      } else {
        isNan = i;
      }
    } else {
      if(fabs(b[i]) > pow10(-FLT_DIG)) {
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

  colinearityThreshold = pow10(-FLT_DIG);
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
                                            vector<T> &baryCentrics,
                                            const int &dimension) {

  baryCentrics.resize(2);

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
  vector<T> test(dimension);
  for(int i = 0; i < dimension; i++) {
    test[i] = baryCentrics[0] * p0[i] + baryCentrics[1] * p1[i];
  }

  if((!((fabs(test[0] - p[0]) < pow(10, -FLT_DIG + 1))
        && (fabs(test[1] - p[1]) < pow(10, -FLT_DIG + 1))))) {
    for(int i = 0; i < 2; i++) {
      baryCentrics[i] = -baryCentrics[i];
    }
  }

  return 0;
}
template <typename T>
int Geometry::computeBarycentricCoordinates(
  const T *p0, const T *p1, const T *p2, const T *p, vector<T> &baryCentrics) {

  baryCentrics.resize(3);

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

  if(fabs(d) < pow(10, -DBL_DIG)) {
    return false;
  }

  x = ((xC - xD) * (xA * yB - yA * xB) - (xA - xB) * (xC * yD - yC * xD)) / d;

  y = ((yC - yD) * (xA * yB - yA * xB) - (yA - yB) * (xC * yD - yC * xD)) / d;

  if((x < std::min(xA, xB) - pow10(-FLT_DIG))
     || (x > std::max(xA, xB) + pow10(-FLT_DIG))) {
    return false;
  }

  if((x < std::min(xC, xD) - pow10(-FLT_DIG))
     || (x > std::max(xC, xD) + pow10(-FLT_DIG))) {
    return false;
  }

  return true;
}

template <typename T>
int Geometry::computeTriangleArea(const T *p0,
                                  const T *p1,
                                  const T *p2,
                                  T &area) {

  vector<T> cross;

  crossProduct(p0, p1, p1, p2, cross);

  area = 0.5 * magnitude(cross.data());

  return 0;
}

template <typename T>
int Geometry::computeTriangleAngles(const T *p0,
                                    const T *p1,
                                    const T *p2,
                                    vector<T> &angles) {

  angles.resize(3);

  angles[0] = angle(p0, p1, p1, p2);
  angles[1] = angle(p1, p2, p2, p0);
  angles[2] = angle(p2, p0, p0, p1);

  return 0;
}

template <typename T>
int Geometry::crossProduct(const T *vA0,
                           const T *vA1,
                           const T *vB0,
                           const T *vB1,
                           vector<T> &crossProduct) {

  crossProduct.resize(3);

  vector<T> a(3), b(3);

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
T Geometry::dotProduct(const T *vA0, const T *vA1, const T *vB0, const T *vB1) {

  T dotProduct = 0;
  for(int i = 0; i < 3; i++) {
    dotProduct += (vA1[i] - vA0[i]) * (vB1[i] - vB0[i]);
  }

  return dotProduct;
}

template <typename T>
T Geometry::dotProduct(const T *vA, const T *vB) {
  return vA[0] * vB[0] + vA[1] * vB[1] + vA[2] * vB[2];
}

template <typename T>
int Geometry::getBoundingBox(const vector<vector<float>> &points,
                             vector<pair<T, T>> &bBox) {

  if(points.empty()) {
    return -1;
  }

  int dimension = points[0].size();

  bBox.resize(dimension);

  for(SimplexId i = 0; i < static_cast<SimplexId>(points.size()); i++) {

    if(i == 0) {
      for(int j = 0; j < dimension; j++) {
        bBox[j].first = points[i][j];
        bBox[j].second = points[i][j];
      }
    } else {
      for(int j = 0; j < dimension; j++) {
        if(points[i][j] < bBox[j].first) {
          bBox[j].first = points[i][j];
        }
        if(points[i][j] > bBox[j].second) {
          bBox[j].second = points[i][j];
        }
      }
    }
  }

  return 0;
}

template <typename T>
bool Geometry::isPointInTriangle(const T *p0,
                                 const T *p1,
                                 const T *p2,
                                 const T *p) {

  vector<T> barycentrics;

  Geometry::computeBarycentricCoordinates(p0, p1, p2, p, barycentrics);

  for(int i = 0; i < static_cast<int>(barycentrics.size()); i++) {
    if(barycentrics[i] < -pow10(-DBL_DIG)) {
      return false;
    }
    if(barycentrics[i] > 1 + pow10(-DBL_DIG)) {
      return false;
    }
  }

  return true;
}

template <typename T>
bool Geometry::isPointOnSegment(
  const T &x, const T &y, const T &xA, const T &yA, const T &xB, const T &yB) {

  vector<T> pA(2), pB(2), p(2);

  pA[0] = xA;
  pA[1] = yA;

  pB[0] = xB;
  pB[1] = yB;

  p[0] = x;
  p[1] = y;

  return Geometry::isPointOnSegment(p.data(), pA.data(), pB.data(), 2);
}

template <typename T>
bool Geometry::isPointOnSegment(const T *p,
                                const T *pA,
                                const T *pB,
                                const int &dimension) {

  vector<T> baryCentrics;

  Geometry::computeBarycentricCoordinates(pA, pB, p, baryCentrics, dimension);

  return (((baryCentrics[0] > -pow10(-DBL_DIG))
           && (baryCentrics[0] < 1 + pow10(-DBL_DIG)))
          && ((baryCentrics[1] > -pow10(-DBL_DIG))
              && (baryCentrics[1] < 1 + pow10(-DBL_DIG))));
}

template <typename T>
bool Geometry::isTriangleColinear(const T *p0,
                                  const T *p1,
                                  const T *p2,
                                  const T *tolerance) {

  bool maxDecision = false;
  T maxCoefficient = 0;
  vector<T> coefficients(3);

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
T Geometry::magnitude(const T *v) {

  T mag = 0;

  for(int i = 0; i < 3; i++) {
    mag += v[i] * v[i];
  }

  return sqrt(mag);
}

template <typename T>
T Geometry::magnitude(const T *o, const T *d) {

  T mag = 0;

  for(int i = 0; i < 3; i++) {
    mag += (o[i] - d[i]) * (o[i] - d[i]);
  }

  return sqrt(mag);
}

// explicit instantiations for double

template double Geometry::angle<double>(double const *,
                                        double const *,
                                        double const *,
                                        double const *);
template bool Geometry::areVectorsColinear<double>(double const *,
                                                   double const *,
                                                   double const *,
                                                   double const *,
                                                   std::vector<double> *,
                                                   double const *);
template int
  Geometry::computeBarycentricCoordinates<double>(double const *,
                                                  double const *,
                                                  double const *,
                                                  std::vector<double> &,
                                                  int const &);
// template int Geometry::computeBarycentricCoordinates<double>(
//     double const &, double const &, double const &, double const &,
//     double const &, double const &, std::vector<double> &);
template int
  Geometry::computeBarycentricCoordinates<double>(double const *,
                                                  double const *,
                                                  double const *,
                                                  double const *,
                                                  std::vector<double> &);
template bool Geometry::computeSegmentIntersection<double>(double const &,
                                                           double const &,
                                                           double const &,
                                                           double const &,
                                                           double const &,
                                                           double const &,
                                                           double const &,
                                                           double const &,
                                                           double &,
                                                           double &);
template int Geometry::computeTriangleAngles<double>(double const *,
                                                     double const *,
                                                     double const *,
                                                     std::vector<double> &);
template int Geometry::computeTriangleArea<double>(double const *,
                                                   double const *,
                                                   double const *,
                                                   double &);
template int Geometry::crossProduct<double>(double const *,
                                            double const *,
                                            double const *,
                                            double const *,
                                            std::vector<double> &);
template int
  Geometry::crossProduct<double>(double const *, double const *, double *);
template double
  Geometry::distance<double>(double const *, double const *, int const &);
template double Geometry::dotProduct<double>(double const *,
                                             double const *,
                                             double const *,
                                             double const *);
template double Geometry::dotProduct<double>(double const *, double const *);
template int
  Geometry::getBoundingBox<double>(std::vector<std::vector<float>> const &,
                                   std::vector<std::pair<double, double>> &);
template bool Geometry::isPointInTriangle<double>(double const *,
                                                  double const *,
                                                  double const *,
                                                  double const *);
template bool Geometry::isPointOnSegment<double>(double const &,
                                                 double const &,
                                                 double const &,
                                                 double const &,
                                                 double const &,
                                                 double const &);
template bool Geometry::isPointOnSegment<double>(double const *,
                                                 double const *,
                                                 double const *,
                                                 int const &);
template bool Geometry::isTriangleColinear<double>(double const *,
                                                   double const *,
                                                   double const *,
                                                   double const *);
template double Geometry::magnitude<double>(double const *);
template double Geometry::magnitude<double>(double const *, double const *);

// explicit instantiations for float

template float Geometry::angle<float>(float const *,
                                      float const *,
                                      float const *,
                                      float const *);
template bool Geometry::areVectorsColinear<float>(float const *,
                                                  float const *,
                                                  float const *,
                                                  float const *,
                                                  std::vector<float> *,
                                                  float const *);
template int
  Geometry::computeBarycentricCoordinates<float>(float const *,
                                                 float const *,
                                                 float const *,
                                                 std::vector<float> &,
                                                 int const &);
// template int Geometry::computeBarycentricCoordinates<float>(
//     float const &, float const &, float const &, float const &,
//     float const &, float const &, std::vector<float> &);
template int
  Geometry::computeBarycentricCoordinates<float>(float const *,
                                                 float const *,
                                                 float const *,
                                                 float const *,
                                                 std::vector<float> &);
template bool Geometry::computeSegmentIntersection<float>(float const &,
                                                          float const &,
                                                          float const &,
                                                          float const &,
                                                          float const &,
                                                          float const &,
                                                          float const &,
                                                          float const &,
                                                          float &,
                                                          float &);
template int Geometry::computeTriangleAngles<float>(float const *,
                                                    float const *,
                                                    float const *,
                                                    std::vector<float> &);
template int Geometry::computeTriangleArea<float>(float const *,
                                                  float const *,
                                                  float const *,
                                                  float &);
template int Geometry::crossProduct<float>(float const *,
                                           float const *,
                                           float const *,
                                           float const *,
                                           std::vector<float> &);
template int
  Geometry::crossProduct<float>(float const *, float const *, float *);
template float
  Geometry::distance<float>(float const *, float const *, int const &);
template float Geometry::dotProduct<float>(float const *,
                                           float const *,
                                           float const *,
                                           float const *);
template float Geometry::dotProduct<float>(float const *, float const *);
template int
  Geometry::getBoundingBox<float>(std::vector<std::vector<float>> const &,
                                  std::vector<std::pair<float, float>> &);
template bool Geometry::isPointInTriangle<float>(float const *,
                                                 float const *,
                                                 float const *,
                                                 float const *);
template bool Geometry::isPointOnSegment<float>(float const &,
                                                float const &,
                                                float const &,
                                                float const &,
                                                float const &,
                                                float const &);
template bool Geometry::isPointOnSegment<float>(float const *,
                                                float const *,
                                                float const *,
                                                int const &);
template bool Geometry::isTriangleColinear<float>(float const *,
                                                  float const *,
                                                  float const *,
                                                  float const *);
template float Geometry::magnitude<float>(float const *);
template float Geometry::magnitude<float>(float const *, float const *);
