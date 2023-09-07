/// \ingroup base
/// \class ttk::ProjectionFromTable
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2022.
///
/// This module defines the %ProjectionFromTable class that
/// projects on a surface points in a vtkTable.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n

#pragma once

// ttk common includes
#include <Debug.h>
#include <Geometry.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The ProjectionFromTable class provides methods that
   * projects on a surface points in a vtkTable.
   */
  class ProjectionFromTable : virtual public Debug {

  protected:
  public:
    ProjectionFromTable();

    template <class triangulationType, class xDataType, class yDataType>
    void computeInputPoints(
      const triangulationType *triangulation,
      std::vector<std::tuple<int, double, double>> &surfaceValues,
      std::array<const long, 2> &surfaceDim,
      const xDataType *const tableXValues,
      const yDataType *const tableYValues,
      const size_t nTableValues,
      std::vector<std::vector<double>> &inputPoints) {
      unsigned int const noPoints = nTableValues;
      inputPoints = std::vector<std::vector<double>>(
        noPoints, std::vector<double>(3, 0.0));

      // Find quadrant of each point
      // 0 --- 2
      // |     |
      // 1 --- 3
      std::vector<std::array<int, 4>> quadPoints(noPoints);
      for(unsigned int i = 0; i < noPoints; ++i) {
        int i0 = 0, i1 = 0;
        for(unsigned int j = 0; j < surfaceDim[0]; ++j)
          if(std::get<1>(surfaceValues[j * surfaceDim[1]]) < tableXValues[i])
            i0 = j;
        for(unsigned int j = 0; j < surfaceDim[1]; ++j)
          if(std::get<2>(surfaceValues[j]) < tableYValues[i])
            i1 = j;
        quadPoints[i][0] = i0 * surfaceDim[1] + i1 + 1;
        quadPoints[i][1] = i0 * surfaceDim[1] + i1;
        quadPoints[i][2] = (i0 + 1) * surfaceDim[1] + i1 + 1;
        quadPoints[i][3] = (i0 + 1) * surfaceDim[1] + i1;
      }

      // Iterate through each point
      std::vector<std::array<double, 4>> coef(noPoints);
      for(unsigned int i = 0; i < noPoints; ++i) {
        std::vector<std::array<double, 2>> points(4);
        for(unsigned int j = 0; j < 4; ++j) {
          points[j][0] = std::get<1>(surfaceValues[quadPoints[i][j]]);
          points[j][1] = std::get<2>(surfaceValues[quadPoints[i][j]]);
        }

        // Find barycentric coordinates
        std::array<double, 3> tableValues{static_cast<double>(tableXValues[i]),
                                          static_cast<double>(tableYValues[i]),
                                          0};
        std::array<std::array<double, 3>, 3> trianglePoints;
        std::array<int, 3> indexes;
        std::array<double, 2> mid{(points[3][0] - points[1][0]) / 2.0,
                                  (points[0][1] - points[1][1]) / 2.0};
        if(tableXValues[i] > points[1][0] + mid[0]) {
          indexes[0] = 2;
          indexes[1] = 3;
          if(tableYValues[i] > points[1][1] + mid[1])
            indexes[2] = 0;
          else
            indexes[2] = 1;
        } else {
          indexes[0] = 0;
          indexes[1] = 1;
          if(tableYValues[i] > points[1][1] + mid[1])
            indexes[2] = 2;
          else
            indexes[2] = 3;
        }

        for(unsigned int j = 0; j < 3; ++j) {
          for(unsigned int k = 0; k < 2; ++k)
            trianglePoints[j][k] = points[indexes[j]][k];
          trianglePoints[j][2] = 0;
        }

        Geometry::computeTriangleArea(
          tableValues.data(), trianglePoints[0].data(),
          trianglePoints[1].data(), coef[i][indexes[2]]);
        Geometry::computeTriangleArea(
          tableValues.data(), trianglePoints[0].data(),
          trianglePoints[2].data(), coef[i][indexes[1]]);
        Geometry::computeTriangleArea(
          tableValues.data(), trianglePoints[1].data(),
          trianglePoints[2].data(), coef[i][indexes[0]]);

        double const sumArea
          = coef[i][indexes[2]] + coef[i][indexes[1]] + coef[i][indexes[0]];
        for(unsigned int j = 0; j < 3; ++j)
          coef[i][indexes[j]] /= sumArea;

        // Compute new point
        for(unsigned int k = 0; k < 4; ++k) {
          float surfacePoint[3];
          triangulation->getVertexPoint(
            std::get<0>(surfaceValues[quadPoints[i][k]]), surfacePoint[0],
            surfacePoint[1], surfacePoint[2]);
          for(unsigned int j = 0; j < 3; ++j)
            inputPoints[i][j] += coef[i][k] * surfacePoint[j];
        }
      }
    }
  }; // ProjectionFromTable class

} // namespace ttk
