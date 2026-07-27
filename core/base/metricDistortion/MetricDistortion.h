/// \ingroup base
/// \class ttk::MetricDistortion
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2022.
///
/// This module defines the %MetricDistortion class that computes distance,
/// area and curvature information about a surface and an optional distance
/// matrix (giving the distance between the points of the surface in a metric
/// space).
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Geometry.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The MetricDistortion class provides methods to compute distance,
   * area and curvature information about a surface and an optional distance
   * matrix (giving the distance between the points of the surface in a metric
   * space).
   */
  class MetricDistortion : virtual public Debug {

  public:
    MetricDistortion();

    int preconditionTriangulation(AbstractTriangulation *triangulation) {
      // Pre-condition functions.
      if(triangulation) {
        triangulation->preconditionBoundaryVertices();
        triangulation->preconditionVertexNeighbors();
      }
      return 0;
    }

    template <class triangulationType, class tableDataType>
    void computeSurfaceCurvature(const triangulationType *triangulation,
                                 std::vector<tableDataType *> &distanceMatrix,
                                 std::vector<double> &surfaceCurvature,
                                 std::vector<double> &metricCurvature,
                                 std::vector<double> &diffCurvature) {
      unsigned int const dim = triangulation->getNumberOfVertices();
      surfaceCurvature = std::vector<double>(dim, std::nan(""));
      metricCurvature = std::vector<double>(dim, std::nan(""));
      diffCurvature = std::vector<double>(dim, std::nan(""));

      std::vector<std::vector<std::tuple<int, int>>> point2CellPoints(dim);
      unsigned int noTriangle = 0, noQuad = 0;
      for(int i = 0; i < triangulation->getNumberOfCells(); ++i) {
        auto cellNoPoints = triangulation->getCellVertexNumber(i);
        if(cellNoPoints < 3 or cellNoPoints > 4)
          continue;
        for(int j = 0; j < cellNoPoints; ++j) {
          auto first = (j + 1) % cellNoPoints;
          auto second = (j + (cellNoPoints - 1)) % cellNoPoints;

          SimplexId i0, i1;
          triangulation->getCellVertex(i, first, i0);
          triangulation->getCellVertex(i, second, i1);

          std::tuple<int, int> const tup = std::make_tuple(i0, i1);

          SimplexId ij;
          triangulation->getCellVertex(i, j, ij);
          point2CellPoints[ij].push_back(tup);
        }
        noTriangle += (cellNoPoints == 3);
        noQuad += (cellNoPoints == 4);
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < dim; ++i) {
        double sumAngleSurface = 0.0, sumAngleMetric = 0.0;

        for(unsigned int j = 0; j < point2CellPoints[i].size(); ++j) {
          auto tup = point2CellPoints[i][j];
          int const i0 = std::get<0>(tup);
          int const i1 = std::get<1>(tup);

          float p_i0[3], p_i1[3], p_i[3];
          triangulation->getVertexPoint(i0, p_i0[0], p_i0[1], p_i0[2]);
          triangulation->getVertexPoint(i1, p_i1[0], p_i1[1], p_i1[2]);
          triangulation->getVertexPoint(i, p_i[0], p_i[1], p_i[2]);

          double const dist_i_i0 = Geometry::distance(&p_i[0], &p_i0[0]);
          double const dist_i_i1 = Geometry::distance(&p_i[0], &p_i1[0]);
          double const dist_i0_i1 = Geometry::distance(&p_i0[0], &p_i1[0]);
          double angleSurface;
          Geometry::computeTriangleAngleFromSides(
            dist_i_i0, dist_i_i1, dist_i0_i1, angleSurface);
          sumAngleSurface += angleSurface;
          if(distanceMatrix.size() != 0) {
            double angleMetric;
            double const distMat_i_i0 = distanceMatrix[i][i0];
            double const distMat_i_i1 = distanceMatrix[i][i1];
            double const distMat_i0_i1 = distanceMatrix[i0][i1];
            Geometry::computeTriangleAngleFromSides(
              distMat_i_i0, distMat_i_i1, distMat_i0_i1, angleMetric);
            sumAngleMetric += angleMetric;
          }
        }

        unsigned int const cornerNoCell = (noTriangle > noQuad ? 2 : 1);
        double const coef
          = (point2CellPoints[i].size() <= cornerNoCell
               ? 0.5
               : (triangulation->isVertexOnBoundary(i) ? 1 : 2));
        surfaceCurvature[i] = coef * M_PI - sumAngleSurface;
        // surfaceCurvature[i] *= std::pow(coef, -1);

        if(distanceMatrix.size() != 0) {
          metricCurvature[i] = coef * M_PI - sumAngleMetric;
          // metricCurvature[i] *= std::pow(coef, -1);
          diffCurvature[i] = metricCurvature[i] - surfaceCurvature[i];
        }
      }
    }

    void initIndicators(double &min, double &max, double &avg) {
      min = std::numeric_limits<double>::max();
      max = std::numeric_limits<double>::lowest();
      avg = 0.0;
    }

    void updateIndicators(
      double &min, double &max, double &avg, double value, int noValue) {
      min = std::min(min, value);
      max = std::max(max, value);
      avg += value / noValue;
    }

    template <class triangulationType, class tableDataType>
    void computeSurfaceDistance(
      const triangulationType *triangulation,
      std::vector<tableDataType *> &distanceMatrix,
      std::vector<double> &surfaceDistance,
      std::vector<double> &metricDistance,
      std::vector<double> &ratioDistance,
      std::vector<std::array<double, 3>> &surfacePointDistance,
      std::vector<std::array<double, 3>> &metricPointDistance,
      std::vector<std::array<double, 3>> &ratioPointDistance) {
      unsigned int dim = triangulation->getNumberOfCells();
      surfaceDistance = std::vector<double>(dim, std::nan(""));
      metricDistance = std::vector<double>(dim, std::nan(""));
      ratioDistance = std::vector<double>(dim, std::nan(""));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < dim; ++i) {
        if(triangulation->getCellVertexNumber(i) != 2)
          continue;

        SimplexId i0, i1;
        triangulation->getCellVertex(i, 0, i0);
        triangulation->getCellVertex(i, 1, i1);

        float p0[3], p1[3];
        triangulation->getVertexPoint(i0, p0[0], p0[1], p0[2]);
        triangulation->getVertexPoint(i1, p1[0], p1[1], p1[2]);

        surfaceDistance[i] = Geometry::distance(&p0[0], &p1[0]);

        if(distanceMatrix.size() != 0) {
          metricDistance[i] = distanceMatrix[i0][i1];
          ratioDistance[i] = metricDistance[i] / surfaceDistance[i];
        }
      }

      dim = triangulation->getNumberOfVertices();
      surfacePointDistance = std::vector<std::array<double, 3>>(dim);
      metricPointDistance = std::vector<std::array<double, 3>>(dim);
      ratioPointDistance = std::vector<std::array<double, 3>>(dim);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < dim; ++i) {
        double minDistanceS, maxDistanceS, avgDistanceS;
        initIndicators(minDistanceS, maxDistanceS, avgDistanceS);
        double minDistanceM, maxDistanceM, avgDistanceM;
        initIndicators(minDistanceM, maxDistanceM, avgDistanceM);

        auto neighborNum = triangulation->getVertexNeighborNumber(i);
        for(int j = 0; j < neighborNum; ++j) {
          SimplexId neighbor;
          triangulation->getVertexNeighbor(i, j, neighbor);

          float p0[3], p1[3];
          triangulation->getVertexPoint(i, p0[0], p0[1], p0[2]);
          triangulation->getVertexPoint(neighbor, p1[0], p1[1], p1[2]);

          double const distance = Geometry::distance(&p0[0], &p1[0]);
          updateIndicators(
            minDistanceS, maxDistanceS, avgDistanceS, distance, neighborNum);

          if(distanceMatrix.size() != 0) {
            double const distanceM = distanceMatrix[i][neighbor];
            updateIndicators(
              minDistanceM, maxDistanceM, avgDistanceM, distanceM, neighborNum);
          }
        }
        surfacePointDistance[i][0] = minDistanceS;
        surfacePointDistance[i][1] = maxDistanceS;
        surfacePointDistance[i][2] = avgDistanceS;
        if(distanceMatrix.size() != 0) {
          metricPointDistance[i][0] = minDistanceM;
          metricPointDistance[i][1] = maxDistanceM;
          metricPointDistance[i][2] = avgDistanceM;
          for(unsigned int j = 0; j < 3; ++j)
            ratioPointDistance[i][j]
              = metricPointDistance[i][j] / surfacePointDistance[i][j];
        }
      }
    }

    template <class triangulationType, class tableDataType>
    void computeSurfaceArea(const triangulationType *triangulation,
                            std::vector<tableDataType *> &distanceMatrix,
                            std::vector<double> &surfaceArea,
                            std::vector<double> &metricArea,
                            std::vector<double> &ratioArea) {
      unsigned int const dim = triangulation->getNumberOfCells();
      surfaceArea = std::vector<double>(dim, std::nan(""));
      metricArea = std::vector<double>(dim, std::nan(""));
      ratioArea = std::vector<double>(dim, std::nan(""));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < dim; ++i) {
        auto cellNoPoints = triangulation->getCellVertexNumber(i);
        if(cellNoPoints < 3 or cellNoPoints > 4)
          continue;

        SimplexId i0, i1, i2;
        triangulation->getCellVertex(i, 0, i0);
        triangulation->getCellVertex(i, 1, i1);
        triangulation->getCellVertex(i, 2, i2);

        float p0[3], p1[3], p2[3];
        triangulation->getVertexPoint(i0, p0[0], p0[1], p0[2]);
        triangulation->getVertexPoint(i1, p1[0], p1[1], p1[2]);
        triangulation->getVertexPoint(i2, p2[0], p2[1], p2[2]);

        float area;
        Geometry::computeTriangleArea(&p0[0], &p1[0], &p2[0], area);
        surfaceArea[i] = area;

        if(distanceMatrix.size() != 0) {
          double areaMetric;
          double const distMat_i0_i1 = distanceMatrix[i0][i1];
          double const distMat_i1_i2 = distanceMatrix[i1][i2];
          double const distMat_i2_i0 = distanceMatrix[i2][i0];
          Geometry::computeTriangleAreaFromSides(
            distMat_i0_i1, distMat_i1_i2, distMat_i2_i0, areaMetric);
          metricArea[i] = areaMetric;
        }

        if(cellNoPoints == 4) {
          SimplexId i3;
          triangulation->getCellVertex(i, 3, i3);
          float p3[3];
          triangulation->getVertexPoint(i3, p3[0], p3[1], p3[2]);

          float areaSurface;
          Geometry::computeTriangleArea(&p1[0], &p2[0], &p3[0], areaSurface);
          surfaceArea[i] += areaSurface;

          if(distanceMatrix.size() != 0) {
            double areaMetric;
            double const distMat_i1_i2 = distanceMatrix[i1][i2];
            double const distMat_i2_i3 = distanceMatrix[i2][i3];
            double const distMat_i3_i1 = distanceMatrix[i3][i1];
            Geometry::computeTriangleAreaFromSides(
              distMat_i1_i2, distMat_i2_i3, distMat_i3_i1, areaMetric);
            metricArea[i] += areaMetric;
          }
        }

        if(distanceMatrix.size() != 0)
          ratioArea[i] = metricArea[i] / surfaceArea[i];
      }
    }

  }; // MetricDistortion class

} // namespace ttk
