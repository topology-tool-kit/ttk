/// \ingroup base
/// \class ttk::ContinuousScatterPlot
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date March 2016
///
/// \brief TTK processing package that computes the continuous scatterplot of
/// bivariate volumetric data.
///
/// \b Related \b publication \n
/// "Continuous Scatterplots" \n
/// Sven Bachthaler, Daniel Weiskopf \n
/// Proc. of IEEE VIS 2008.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2008.
///
/// \sa ttkContinuousScatterPlot.cpp %for a usage example.

#pragma once

#include <limits>

// base code includes
#include <Geometry.h>
#include <Triangulation.h>

namespace ttk {

  class ContinuousScatterPlot : virtual public Debug {

  public:
    ContinuousScatterPlot();
    ~ContinuousScatterPlot();

    template <typename dataType1,
              typename dataType2,
              class triangulationType = AbstractTriangulation>
    int execute(const dataType1 *,
                const dataType2 *,
                const triangulationType *) const;

    inline int setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    inline int setDummyValue(bool withDummyValue, double dummyValue) {
      if(withDummyValue) {
        withDummyValue_ = true;
        dummyValue_ = dummyValue;
      }
      return 0;
    }

    inline int setResolutions(const SimplexId &resolutionX,
                              const SimplexId &resolutionY) {
      resolutions_[0] = resolutionX;
      resolutions_[1] = resolutionY;
      return 0;
    }

    inline int setScalarMin(double *scalarMin) {
      scalarMin_ = scalarMin;
      return 0;
    }

    inline int setScalarMax(double *scalarMax) {
      scalarMax_ = scalarMax;
      return 0;
    }

    inline int setOutputDensity(std::vector<std::vector<double>> *density) {
      density_ = density;
      return 0;
    }

    inline int setOutputMask(std::vector<std::vector<char>> *mask) {
      validPointMask_ = mask;
      return 0;
    }

  protected:
    SimplexId vertexNumber_;
    bool withDummyValue_;
    double dummyValue_;
    SimplexId resolutions_[2];
    double *scalarMin_;
    double *scalarMax_;
    std::vector<std::vector<double>> *density_;
    std::vector<std::vector<char>> *validPointMask_;
  };
} // namespace ttk

template <typename dataType1, typename dataType2, class triangulationType>
int ttk::ContinuousScatterPlot::execute(
  const dataType1 *scalars1,
  const dataType2 *scalars2,
  const triangulationType *triangulation) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(!scalars1)
    return -1;
  if(!scalars2)
    return -2;
  if(!triangulation)
    return -3;
  if(!density_)
    return -4;

  if(triangulation->getNumberOfCells() <= 0) {
    this->printErr("no cells.");
    return -5;
  }

  if(triangulation->getCellVertexNumber(0) != 4) {
    this->printErr("no tetrahedra.");
    return -6;
  }
#endif

  Timer t;

  // helpers:
  const SimplexId numberOfCells = triangulation->getNumberOfCells();

  // rendering helpers:
  // constant ray direction (ortho)
  const double d[3]{0, 0, -1};
  const double delta[2]{
    scalarMax_[0] - scalarMin_[0], scalarMax_[1] - scalarMin_[1]};
  const double sampling[2]{
    delta[0] / resolutions_[0], delta[1] / resolutions_[1]};
  const double epsilon{0.000001};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId cell = 0; cell < numberOfCells; ++cell) {
    bool isDummy{};

    // get tetrahedron info
    SimplexId vertex[4];
    double data[4][3];
    float position[4][3];
    double localScalarMin[2]{};
    double localScalarMax[2]{};
    // for each triangle
    for(int k = 0; k < 4; ++k) {
      // get indices
      triangulation->getCellVertex(cell, k, vertex[k]);

      // get scalars
      data[k][0] = scalars1[vertex[k]];
      data[k][1] = scalars2[vertex[k]];
      data[k][2] = 0;

      if(withDummyValue_
         and (data[k][0] == dummyValue_ or data[k][1] == dummyValue_)) {
        isDummy = true;
        break;
      }

      // get local stats
      if(!k or localScalarMin[0] > data[k][0])
        localScalarMin[0] = data[k][0];
      if(!k or localScalarMin[1] > data[k][1])
        localScalarMin[1] = data[k][1];
      if(!k or localScalarMax[0] < data[k][0])
        localScalarMax[0] = data[k][0];
      if(!k or localScalarMax[1] < data[k][1])
        localScalarMax[1] = data[k][1];

      // get positions

      triangulation->getVertexPoint(
        vertex[k], position[k][0], position[k][1], position[k][2]);
    }
    if(isDummy)
      continue;

    // gradient:
    double g0[3];
    double g1[3];
    {
      double v12[3];
      double v13[3];
      double v14[3];
      double s12[3];
      double s13[3];
      double s14[3];
      for(int k = 0; k < 3; ++k) {
        v12[k] = position[1][k] - position[0][k];
        v13[k] = position[2][k] - position[0][k];
        v14[k] = position[3][k] - position[0][k];

        s12[k] = data[1][k] - data[0][k];
        s13[k] = data[2][k] - data[0][k];
        s14[k] = data[3][k] - data[0][k];
      }

      double a[3];
      double b[3];
      double c[3];
      Geometry::crossProduct(v13, v12, a);
      Geometry::crossProduct(v12, v14, b);
      Geometry::crossProduct(v14, v13, c);
      double det = Geometry::dotProduct(v14, a);
      if(det == 0.) {
        for(int k = 0; k < 3; ++k) {
          g0[k] = 0.0;
          g1[k] = 0.0;
        }
      } else {
        double invDet = 1.0 / det;
        for(int k = 0; k < 3; ++k) {
          g0[k] = (s14[0] * a[k] + s13[0] * b[k] + s12[0] * c[k]) * invDet;
          g1[k] = (s14[1] * a[k] + s13[1] * b[k] + s12[1] * c[k]) * invDet;
        }
      }
    }

    // volume:
    double volume;
    bool isLimit{};
    {
      double cp[3];
      Geometry::crossProduct(g0, g1, cp);
      volume = Geometry::magnitude(cp);
      if(volume == 0.)
        isLimit = true;
    }

    // classification:
    int index[4]{0, 1, 2, 3};
    bool isInTriangle{};
    if(Geometry::isPointInTriangle(data[0], data[1], data[2], data[3]))
      isInTriangle = true;
    else if(Geometry::isPointInTriangle(data[0], data[1], data[3], data[2])) {
      isInTriangle = true;
      index[0] = 0;
      index[1] = 1;
      index[2] = 3;
      index[3] = 2;
    } else if(Geometry::isPointInTriangle(data[0], data[2], data[3], data[1])) {
      isInTriangle = true;
      index[0] = 0;
      index[1] = 2;
      index[2] = 3;
      index[3] = 1;
    } else if(Geometry::isPointInTriangle(data[1], data[2], data[3], data[0])) {
      isInTriangle = true;
      index[0] = 1;
      index[1] = 2;
      index[2] = 3;
      index[3] = 0;
    }

    // projection:
    double density{};
    std::vector<std::vector<SimplexId>> triangles;
    double imaginaryPosition[3]{};
    std::vector<SimplexId> triangle(3);
    // class 0
    if(isInTriangle) {
      // mass density
      double massDensity{};
      {
        double A;

        Geometry::computeTriangleArea(
          data[index[0]], data[index[1]], data[index[2]], A);
        double invA = 1.0 / A;
        if(A == 0.) {
          invA = 0.0;
          isLimit = true;
        }

        double alpha, beta, gamma;

        Geometry::computeTriangleArea(
          data[index[1]], data[index[2]], data[index[3]], alpha);

        Geometry::computeTriangleArea(
          data[index[0]], data[index[2]], data[index[3]], beta);

        Geometry::computeTriangleArea(
          data[index[0]], data[index[1]], data[index[3]], gamma);

        alpha *= invA;
        beta *= invA;
        gamma *= invA;

        double p0[3];
        double p1[3];
        for(int k = 0; k < 3; ++k) {
          p0[k] = position[index[3]][k];

          p1[k] = alpha * position[index[0]][k] + beta * position[index[1]][k]
                  + gamma * position[index[2]][k];
        }
        massDensity = Geometry::distance(p0, p1);
      }

      if(isLimit)
        density = std::numeric_limits<decltype(density)>::max();
      else
        density = massDensity / volume;

      triangle[0] = vertex[index[3]];
      triangle[1] = vertex[index[0]];
      triangle[2] = vertex[index[1]];
      triangles.push_back(triangle);

      triangle[0] = vertex[index[3]];
      triangle[1] = vertex[index[0]];
      triangle[2] = vertex[index[2]];
      triangles.push_back(triangle);

      triangle[0] = vertex[index[3]];
      triangle[1] = vertex[index[1]];
      triangle[2] = vertex[index[2]];
      triangles.push_back(triangle);
    }
    // class 1
    else {
      double massDensity{};
      double p[3]{0, 0, 0};
      if(Geometry::computeSegmentIntersection(
           data[0][0], data[0][1], data[1][0], data[1][1], data[2][0],
           data[2][1], data[3][0], data[3][1], p[0], p[1])) {
        index[0] = 0;
        index[1] = 1;
        index[2] = 2;
        index[3] = 3;
      } else if(Geometry::computeSegmentIntersection(
                  data[0][0], data[0][1], data[2][0], data[2][1], data[1][0],
                  data[1][1], data[3][0], data[3][1], p[0], p[1])) {
        index[0] = 0;
        index[1] = 2;
        index[2] = 1;
        index[3] = 3;
      } else if(Geometry::computeSegmentIntersection(
                  data[0][0], data[0][1], data[3][0], data[3][1], data[1][0],
                  data[1][1], data[2][0], data[2][1], p[0], p[1])) {
        index[0] = 0;
        index[1] = 3;
        index[2] = 1;
        index[3] = 2;
      }

      double a = Geometry::distance(data[index[0]], p);
      double b = Geometry::distance(data[index[0]], data[index[1]]);
      double r0 = a / b;

      a = Geometry::distance(data[index[2]], p);
      b = Geometry::distance(data[index[2]], data[index[3]]);
      double r1 = a / b;

      double p0[3];
      double p1[3];
      for(int k = 0; k < 3; ++k) {

        p0[k] = position[index[0]][k]
                + r0 * (position[index[1]][k] - position[index[0]][k]);

        p1[k] = position[index[2]][k]
                + r1 * (position[index[3]][k] - position[index[2]][k]);
      }
      massDensity = Geometry::distance(p0, p1);

      if(isLimit)
        density = std::numeric_limits<decltype(density)>::max();
      else
        density = massDensity / volume;

      imaginaryPosition[0] = p[0];
      imaginaryPosition[1] = p[1];
      imaginaryPosition[2] = 0;

      // four triangles projection
      triangle[0] = -1; // new geometry
      triangle[1] = vertex[index[0]];
      triangle[2] = vertex[index[2]];
      triangles.push_back(triangle);

      triangle[1] = vertex[index[2]];
      triangle[2] = vertex[index[1]];
      triangles.push_back(triangle);

      triangle[1] = vertex[index[1]];
      triangle[2] = vertex[index[3]];
      triangles.push_back(triangle);

      triangle[1] = vertex[index[3]];
      triangle[2] = vertex[index[0]];
      triangles.push_back(triangle);
    }

    // rendering:
    // "Fast, Minimum Storage Ray/Triangle Intersection", Tomas Moller & Ben
    // Trumbore
    {
      const SimplexId minI
        = floor((localScalarMin[0] - scalarMin_[0]) / sampling[0]);
      const SimplexId minJ
        = floor((localScalarMin[1] - scalarMin_[1]) / sampling[1]);
      const SimplexId maxI
        = ceil((localScalarMax[0] - scalarMin_[0]) / sampling[0]);
      const SimplexId maxJ
        = ceil((localScalarMax[1] - scalarMin_[1]) / sampling[1]);

      for(SimplexId i = minI; i < maxI; ++i) {
        for(SimplexId j = minJ; j < maxJ; ++j) {
          // set ray origin
          const double o[3]{scalarMin_[0] + i * sampling[0],
                            scalarMin_[1] + j * sampling[1], 1};
          for(unsigned int k = 0; k < triangles.size(); ++k) {
            const auto &tr = triangles[k];

            // get triangle info
            double p0[3];
            if(isInTriangle) {
              p0[0] = scalars1[tr[0]];
              p0[1] = scalars2[tr[0]];
            } else {
              p0[0] = imaginaryPosition[0];
              p0[1] = imaginaryPosition[1];
            }
            p0[2] = 0;

            const double p1[3]{
              (double)scalars1[tr[1]], (double)scalars2[tr[1]], 0};
            const double p2[3]{
              (double)scalars1[tr[2]], (double)scalars2[tr[2]], 0};
            const double e1[3]{p1[0] - p0[0], p1[1] - p0[1], 0};
            const double e2[3]{p2[0] - p0[0], p2[1] - p0[1], 0};

            double q[3];
            Geometry::crossProduct(d, e2, q);
            const double a = Geometry::dotProduct(e1, q);
            if(a > -epsilon and a < epsilon)
              continue;

            const double f = 1.0 / a;
            const double s[3]{o[0] - p0[0], o[1] - p0[1], 1};
            const double u = f * Geometry::dotProduct(s, q);
            if(u < 0.0)
              continue;

            double r[3];
            Geometry::crossProduct(s, e1, r);
            const double v = f * Geometry::dotProduct(d, r);
            if(v < 0.0 or (u + v) > 1.0)
              continue;

              // triangle/ray intersection below
#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp atomic
#else
#pragma omp atomic update
#endif
#endif
            (*density_)[i][j] += (1.0 - u - v) * density;

#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp atomic
            (*validPointMask_)[i][j] += 1;
#else
#pragma omp atomic write
            (*validPointMask_)[i][j] = 1;
#endif
#else
            (*validPointMask_)[i][j] = 1;
#endif
            break;
          }
        }
      }
    }
  }

  {
    std::stringstream msg;
    msg << "Processed " << numberOfCells << " tetrahedra";
    this->printMsg(msg.str(), 1, t.getElapsedTime(), threadNumber_);
  }

  return 0;
}
