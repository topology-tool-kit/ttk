/// \ingroup base
/// \class ttk::MetricDistortion
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2022.
///
/// This module defines the %MetricDistortion class that computes distance,
/// area and curvature information about a surface and an optionnal distance
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
   * area and curvature information about a surface and an optionnal distance
   * matrix (giving the distance between the points of the surface in a metric
   * space).
   */
  class MetricDistortion : virtual public Debug {

  public:
    MetricDistortion();

    void
      computeSurfaceCurvature(std::vector<std::vector<double>> &surfacePoints,
                              std::vector<std::vector<int>> &surfaceCells,
                              std::vector<std::vector<double>> &distanceMatrix,
                              std::vector<bool> &isPointBoundary,
                              std::vector<double> &surfaceCurvature,
                              std::vector<double> &metricCurvature,
                              std::vector<double> &diffCurvature) {
      unsigned int dim = surfacePoints.size();
      surfaceCurvature = std::vector<double>(dim, std::nan(""));
      metricCurvature = std::vector<double>(dim, std::nan(""));
      diffCurvature = std::vector<double>(dim, std::nan(""));

      std::vector<std::vector<std::tuple<int, int>>> point2CellPoints(dim);
      unsigned int noTriangle = 0, noQuad = 0;
      for(unsigned int i = 0; i < surfaceCells.size(); ++i) {
        auto cellNoPoints = surfaceCells[i].size();
        if(cellNoPoints < 3 or cellNoPoints > 4)
          continue;
        for(unsigned int j = 0; j < cellNoPoints; ++j) {
          auto first = (j + 1) % cellNoPoints;
          auto second = (j + (cellNoPoints - 1)) % cellNoPoints;
          std::tuple<int, int> tup
            = std::make_tuple(surfaceCells[i][first], surfaceCells[i][second]);
          point2CellPoints[surfaceCells[i][j]].push_back(tup);
        }
        noTriangle += (cellNoPoints == 3);
        noQuad += (cellNoPoints == 4);
      }

      for(unsigned int i = 0; i < dim; ++i) {
        double sumAngleSurface = 0.0, sumAngleMetric = 0.0;

        for(unsigned int j = 0; j < point2CellPoints[i].size(); ++j) {
          auto tup = point2CellPoints[i][j];
          int i0 = std::get<0>(tup);
          int i1 = std::get<1>(tup);
          auto p_i0 = surfacePoints[i0];
          auto p_i1 = surfacePoints[i1];
          auto p_i = surfacePoints[i];
          auto dist_i_i0 = Geometry::distance(&p_i[0], &p_i0[0]);
          auto dist_i_i1 = Geometry::distance(&p_i[0], &p_i1[0]);
          auto dist_i0_i1 = Geometry::distance(&p_i0[0], &p_i1[0]);
          sumAngleSurface += cosineLaw(dist_i_i0, dist_i_i1, dist_i0_i1);
          if(distanceMatrix.size() != 0)
            sumAngleMetric
              += cosineLaw(distanceMatrix[i][i0], distanceMatrix[i][i1],
                           distanceMatrix[i0][i1]);
        }

        unsigned int cornerNoCell = (noTriangle > noQuad ? 2 : 1);
        double coef = (point2CellPoints[i].size() <= cornerNoCell
                         ? 0.5
                         : (isPointBoundary[i] ? 1 : 2));
        surfaceCurvature[i] = coef * M_PI - sumAngleSurface;
        // surfaceCurvature[i] *= std::pow(coef, -1);

        if(distanceMatrix.size() != 0) {
          metricCurvature[i] = coef * M_PI - sumAngleMetric;
          // metricCurvature[i] *= std::pow(coef, -1);
          diffCurvature[i] = metricCurvature[i] - surfaceCurvature[i];
        }
      }
    }

    void
      computeSurfaceDistance(std::vector<std::vector<double>> &surfacePoints,
                             std::vector<std::vector<int>> &surfaceCells,
                             std::vector<std::vector<double>> &distanceMatrix,
                             std::vector<double> &surfaceDistance,
                             std::vector<double> &metricDistance,
                             std::vector<double> &ratioDistance) {
      unsigned int dim = surfaceCells.size();
      surfaceDistance = std::vector<double>(dim, std::nan(""));
      metricDistance = std::vector<double>(dim, std::nan(""));
      ratioDistance = std::vector<double>(dim, std::nan(""));
      for(unsigned int i = 0; i < dim; ++i) {
        if(surfaceCells[i].size() != 2)
          continue;

        int i0 = surfaceCells[i][0];
        int i1 = surfaceCells[i][1];

        surfaceDistance[i]
          = Geometry::distance(&surfacePoints[i0][0], &surfacePoints[i1][0]);

        if(distanceMatrix.size() != 0) {
          metricDistance[i] = distanceMatrix[i0][i1];
          ratioDistance[i] = metricDistance[i] / surfaceDistance[i];
        }
      }
    }

    void computeSurfaceArea(std::vector<std::vector<double>> &surfacePoints,
                            std::vector<std::vector<int>> &surfaceCells,
                            std::vector<std::vector<double>> &distanceMatrix,
                            std::vector<double> &surfaceArea,
                            std::vector<double> &metricArea,
                            std::vector<double> &ratioArea) {
      unsigned int dim = surfaceCells.size();
      surfaceArea = std::vector<double>(dim, std::nan(""));
      metricArea = std::vector<double>(dim, std::nan(""));
      ratioArea = std::vector<double>(dim, std::nan(""));
      for(unsigned int i = 0; i < dim; ++i) {
        if(surfaceCells[i].size() < 3 or surfaceCells[i].size() > 4)
          continue;

        int i0 = surfaceCells[i][0];
        int i1 = surfaceCells[i][1];
        int i2 = surfaceCells[i][2];

        Geometry::computeTriangleArea(&surfacePoints[i0][0],
                                      &surfacePoints[i1][0],
                                      &surfacePoints[i2][0], surfaceArea[i]);

        if(distanceMatrix.size() != 0) {
          Geometry::computeTriangleAreaFromSides(
            distanceMatrix[i0][i1], distanceMatrix[i1][i2],
            distanceMatrix[i2][i0], metricArea[i]);
        }

        if(surfaceCells[i].size() == 4) {
          int i3 = surfaceCells[i][3];
          double temp;
          Geometry::computeTriangleArea(&surfacePoints[i1][0],
                                        &surfacePoints[i2][0],
                                        &surfacePoints[i3][0], temp);
          surfaceArea[i] += temp;

          if(distanceMatrix.size() != 0) {
            double tempMetric;
            Geometry::computeTriangleAreaFromSides(
              distanceMatrix[i1][i2], distanceMatrix[i2][i3],
              distanceMatrix[i3][i1], tempMetric);
            metricArea[i] += tempMetric;
          }
        }

        if(distanceMatrix.size() != 0)
          ratioArea[i] = metricArea[i] / surfaceArea[i];
      }
    }

    //-------------------------------------------------------------------------
    // Utils
    //-------------------------------------------------------------------------
    double cosineLaw(double a, double b, double c) {
      // Get the angle opposite to edge c
      return std::acos((a * a + b * b - c * c) / (2.0 * a * b));
    }
  }; // MetricDistortion class

} // namespace ttk
