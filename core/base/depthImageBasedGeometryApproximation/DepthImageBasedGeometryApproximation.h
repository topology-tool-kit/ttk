/// \ingroup base
/// \class ttk::DepthImageBasedGeometryApproximation
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.7.2018
///
/// \brief This module approximates the depicted geometry of a depth image.
///
/// This module approximates the depicted geometry of a depth image solely based
/// on the depth values and the corresponding camera parameter.
///
/// Related publication:
/// 'VOIDGA: A View-Approximation Oriented Image Database Generation Approach'
/// Jonas Lukasczyk, Eric Kinner, James Ahrens, Heike Leitte, and Christoph
/// Garth. IEEE 8th Symposium on Large Data Analysis and Visualization (LDAV),
/// 2018.

#pragma once

// base code includes
#include <Debug.h>

// std includes
#include <cmath>
#include <limits>

namespace ttk {

  class DepthImageBasedGeometryApproximation : virtual public Debug {

  public:
    DepthImageBasedGeometryApproximation() {
      this->setDebugMsgPrefix("DIBGA");
    }
    ~DepthImageBasedGeometryApproximation() = default;

    /**
     * This function computes for an input depth image and its corresponding
     * camera parameters the original 3D position of each pixel and then
     * connects neighbouring pixels via triangles. The function also returns a
     * scalar array (triangleDistortions) that records the distortion of each
     * triangle that can later be used to filter the resulting geometry.
     */
    template <class dataType, class idType>
    int execute(
      // output
      float *pointCoordinates,
      double *triangleDistortions,
      idType *connectivityList,
      idType *offsetArray,

      // input
      const dataType *depthValues,
      const double *camPos,
      const double *camDir,
      const double *camUp,
      const double *camNearFar,
      const double *camHeight,
      const double *resolution) const;
  };
} // namespace ttk

template <class dataType, class idType>
int ttk::DepthImageBasedGeometryApproximation::execute(
  // output
  float *pointCoordinates,
  double *triangleDistortions,
  idType *connectivityList,
  idType *offsetArray,

  // input
  const dataType *depthValues,
  const double *camPos,
  const double *camDir,
  const double *camUp,
  const double *camNearFar,
  const double *camHeight,
  const double *resolution) const {

  const double camDirMag = std::sqrt(
    camDir[0] * camDir[0] + camDir[1] * camDir[1] + camDir[2] * camDir[2]);
  const double camDirN[3]{
    camDir[0] / camDirMag, camDir[1] / camDirMag, camDir[2] / camDirMag};

  this->printMsg(ttk::debug::Separator::L2, ttk::debug::LineMode::NEW,
                 ttk::debug::Priority::DETAIL);
  this->printMsg({{"Resolution", std::to_string((int)resolution[0]) + "x"
                                   + std::to_string((int)resolution[1])},
                  {"CamPos", "[" + std::to_string(camPos[0]) + ","
                               + std::to_string(camPos[1]) + ","
                               + std::to_string(camPos[2]) + "]"},
                  {"CamDir", "[" + std::to_string(camDirN[0]) + ","
                               + std::to_string(camDirN[1]) + ","
                               + std::to_string(camDirN[2]) + "]"},
                  {"CamHeight", std::to_string(camHeight[0])},
                  {"CamNearFar", "[" + std::to_string(camNearFar[0]) + ","
                                   + std::to_string(camNearFar[1]) + "]"}},
                 ttk::debug::Priority::DETAIL);
  this->printMsg(ttk::debug::Separator::L1, ttk::debug::LineMode::NEW,
                 ttk::debug::Priority::DETAIL);

  Timer timer;

  const size_t resolutionST[2] = {(size_t)resolution[0], (size_t)resolution[1]};

  this->printMsg("Processing image (" + std::to_string(resolutionST[0]) + "x"
                   + std::to_string(resolutionST[1]) + ")",
                 0, -1, this->threadNumber_, ttk::debug::LineMode::REPLACE);
  // -------------------------------------------------------------------------
  // Compute Camera Vectors
  // -------------------------------------------------------------------------

  // Compute camera size
  const double camSize[2]
    = {resolution[0] / resolution[1] * camHeight[0], camHeight[0]};

  // Compute camRight = camDirN x CamUp
  double camRight[3] = {camDirN[1] * camUp[2] - camDirN[2] * camUp[1],
                        camDirN[2] * camUp[0] - camDirN[0] * camUp[2],
                        camDirN[0] * camUp[1] - camDirN[1] * camUp[0]};
  double temp = sqrt(camRight[0] * camRight[0] + camRight[1] * camRight[1]
                     + camRight[2] * camRight[2]);
  camRight[0] /= temp;
  camRight[1] /= temp;
  camRight[2] /= temp;

  // Compute true up vector
  double camUpTrue[3]
    = {camDirN[1] * (-camRight[2]) - camDirN[2] * (-camRight[1]),
       camDirN[2] * (-camRight[0]) - camDirN[0] * (-camRight[2]),
       camDirN[0] * (-camRight[1]) - camDirN[1] * (-camRight[0])};
  temp = sqrt(camUpTrue[0] * camUpTrue[0] + camUpTrue[1] * camUpTrue[1]
              + camUpTrue[2] * camUpTrue[2]);
  camUpTrue[0] /= temp;
  camUpTrue[1] /= temp;
  camUpTrue[2] /= temp;

  // -------------------------------------------------------------------------
  // Create Vertices
  // -------------------------------------------------------------------------
  {
    // Compute pixel size in world coordinates
    const double pixelWidthWorld = camSize[0] / resolution[0];
    const double pixelHeightWorld = camSize[1] / resolution[1];

    // Optimization: precompute half of the camera size to reduce the number of
    // operations in the for loop Include a half pixel offset (-0.5) to center
    // vertices at pixel centers
    const double camWidthWorldHalf = 0.5 * camSize[0] - 0.5 * pixelWidthWorld;
    const double camHeightWorldHalf = 0.5 * camSize[1] - 0.5 * pixelHeightWorld;

    // Compute depth delta
    const double delta = camNearFar[1] - camNearFar[0];

    // Optimization: reorient camera model to bottom left corner to reduce
    // operations in for loop
    const double camPosCorner[3] = {camPos[0] - camRight[0] * camWidthWorldHalf
                                      - camUpTrue[0] * camHeightWorldHalf,
                                    camPos[1] - camRight[1] * camWidthWorldHalf
                                      - camUpTrue[1] * camHeightWorldHalf,
                                    camPos[2] - camRight[2] * camWidthWorldHalf
                                      - camUpTrue[2] * camHeightWorldHalf};

// Compute vertex coordinates while parallizing over rows
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t y = 0; y < resolutionST[1]; y++) {
      const double v = ((double)y) * pixelHeightWorld;
      const double vTimesUp[3]
        = {v * camUpTrue[0], v * camUpTrue[1], v * camUpTrue[2]};

      const size_t yOffset = y * resolutionST[0];
      for(size_t x = 0; x < resolutionST[0]; x++) {
        const size_t pixelIndex = x + yOffset;

        // double d = (double)(depthValues[ pixelIndex ])*delta+camNearFar[0];
        const double depth = ((double)depthValues[pixelIndex]);
        // double d = depth > 0.98 ? 0 : depth * delta + camNearFar[0];
        const double d = depth * delta + camNearFar[0];
        const double u = ((double)x) * pixelWidthWorld;

        // compute vertex coordinate
        const size_t pointCoordinateOffset = pixelIndex * 3;
        pointCoordinates[pointCoordinateOffset]
          = camPosCorner[0] + u * camRight[0] + vTimesUp[0] + d * camDirN[0];
        pointCoordinates[pointCoordinateOffset + 1]
          = camPosCorner[1] + u * camRight[1] + vTimesUp[1] + d * camDirN[1];
        pointCoordinates[pointCoordinateOffset + 2]
          = camPosCorner[2] + u * camRight[2] + vTimesUp[2] + d * camDirN[2];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Create Triangles
  // -------------------------------------------------------------------------
  {
    auto absDiff = [](const dataType &a, const dataType &b) {
      return a > b ? a - b : b - a;
    };
    auto isNaN = [](const double &a) { return std::isnan(a) || a >= 1.0; };

    const size_t nTriangles = 2 * (resolution[0] - 1) * (resolution[1] - 1);

    for(size_t t = 0; t < nTriangles; t++)
      offsetArray[t] = t * 3;
    offsetArray[nTriangles] = 3 * nTriangles;

    /* Index Structure:
    0 - 1
    | / |
    2 - 3
    */
    const size_t xl = resolutionST[0] - 1;
    const size_t yl = resolutionST[1] - 1;
    const size_t trianglesPerRow = xl * 2;

    const double myNan = std::numeric_limits<double>::quiet_NaN();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
    for(size_t y = 0; y < yl; y++) {
      size_t yOffset = y * resolutionST[0];
      size_t triangleIndexOffset = y * trianglesPerRow * 3;
      size_t triangleDistortionOffset = y * trianglesPerRow;

      for(size_t x = 0; x < xl; x++) {
        size_t i0 = x + yOffset;
        size_t i1 = i0 + 1;
        size_t i2 = i0 + resolutionST[0];
        size_t i3 = i2 + 1;

        connectivityList[triangleIndexOffset++] = i0;
        connectivityList[triangleIndexOffset++] = i2;
        connectivityList[triangleIndexOffset++] = i1;

        connectivityList[triangleIndexOffset++] = i1;
        connectivityList[triangleIndexOffset++] = i2;
        connectivityList[triangleIndexOffset++] = i3;

        const double i0Depth = (double)depthValues[i0];
        const double i1Depth = (double)depthValues[i1];
        const double i2Depth = (double)depthValues[i2];
        const double i3Depth = (double)depthValues[i3];

        // wow
        triangleDistortions[triangleDistortionOffset++]
          = isNaN(i0Depth) || isNaN(i2Depth) || isNaN(i1Depth)
              ? myNan
              : std::max(
                absDiff(i0Depth, i1Depth),
                std::max(absDiff(i1Depth, i2Depth), absDiff(i0Depth, i2Depth)));

        triangleDistortions[triangleDistortionOffset++]
          = isNaN(i1Depth) || isNaN(i2Depth) || isNaN(i3Depth)
              ? myNan
              : std::max(
                absDiff(i1Depth, i3Depth),
                std::max(absDiff(i3Depth, i2Depth), absDiff(i2Depth, i1Depth)));
      }
    }
  }

  // Print performance
  this->printMsg("Processing image (" + std::to_string(resolutionST[0]) + "x"
                   + std::to_string(resolutionST[1]) + ")",
                 1, timer.getElapsedTime(), this->threadNumber_);

  return 1;
}
