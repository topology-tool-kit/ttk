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
#include <limits>

namespace ttk {

  class DepthImageBasedGeometryApproximation : virtual public Debug {

  public:
    DepthImageBasedGeometryApproximation() {
      this->setDebugMsgPrefix("DIBGA");
    };
    ~DepthImageBasedGeometryApproximation(){};

    /**
     * This function computes for an input depth image and its corresponding
     * camera parameters the original 3D position of each pixel and then
     * connects neighbouring pixels via triangles. The function also returns a
     * scalar array (triangleDistortions) that records the distortion of each
     * triangle that can later be used to filter the resulting geometry.
     */
    template <class dataType, class idType>
    int execute(
      // input
      const dataType *depthValues,
      const double *camPos,
      const double *camDir,
      const double *camUp,
      const double *camNearFar,
      const double *camHeight,
      const double *resolution,

      // output
      float *pointCoordinates,
      idType *connectivityList,
      double *triangleDistortions) const;
  };
} // namespace ttk

template <class dataType, class idType>
int ttk::DepthImageBasedGeometryApproximation::execute(
  // input
  const dataType *depthValues,
  const double *camPos,
  const double *camDir,
  const double *camUp,
  const double *camNearFar,
  const double *camHeight,
  const double *resolution,

  // output
  float *pointCoordinates,
  idType *connectivityList,
  double *triangleDistortions) const {

  this->printMsg(ttk::debug::Separator::L2, ttk::debug::LineMode::NEW,
                 ttk::debug::Priority::DETAIL);
  this->printMsg({{"Resolution", std::to_string((int)resolution[0]) + "x"
                                   + std::to_string((int)resolution[1])},
                  {"CamPos", "[" + std::to_string(camPos[0]) + ","
                               + std::to_string(camPos[1]) + ","
                               + std::to_string(camPos[2]) + "]"},
                  {"CamDir", "[" + std::to_string(camDir[0]) + ","
                               + std::to_string(camDir[1]) + ","
                               + std::to_string(camDir[2]) + "]"},
                  {"CamHeight", std::to_string(camHeight[0])},
                  {"CamNearFar", "[" + std::to_string(camNearFar[0]) + ","
                                   + std::to_string(camNearFar[1]) + "]"}},
                 ttk::debug::Priority::DETAIL);
  this->printMsg(ttk::debug::Separator::L1, ttk::debug::LineMode::NEW,
                 ttk::debug::Priority::DETAIL);

  Timer t;

  size_t resolutionST[2] = {(size_t)resolution[0], (size_t)resolution[1]};

  this->printMsg("Processing image (" + std::to_string(resolutionST[0]) + "x"
                   + std::to_string(resolutionST[1]) + ")",
                 0, -1, this->threadNumber_, ttk::debug::LineMode::REPLACE);
  // -------------------------------------------------------------------------
  // Compute Camera Vectors
  // -------------------------------------------------------------------------

  // Compute camera size
  double camSize[2]
    = {resolution[0] / resolution[1] * camHeight[0], camHeight[0]};

  // Compute camRight = camDir x CamUp
  double camRight[3] = {camDir[1] * camUp[2] - camDir[2] * camUp[1],
                        camDir[2] * camUp[0] - camDir[0] * camUp[2],
                        camDir[0] * camUp[1] - camDir[1] * camUp[0]};
  double temp = sqrt(camRight[0] * camRight[0] + camRight[1] * camRight[1]
                     + camRight[2] * camRight[2]);
  camRight[0] /= temp;
  camRight[1] /= temp;
  camRight[2] /= temp;

  // Compute true up std::vector
  double camUpTrue[3]
    = {camDir[1] * (-camRight[2]) - camDir[2] * (-camRight[1]),
       camDir[2] * (-camRight[0]) - camDir[0] * (-camRight[2]),
       camDir[0] * (-camRight[1]) - camDir[1] * (-camRight[0])};
  temp = sqrt(camUpTrue[0] * camUpTrue[0] + camUpTrue[1] * camUpTrue[1]
              + camUpTrue[2] * camUpTrue[2]);
  camUpTrue[0] /= temp;
  camUpTrue[1] /= temp;
  camUpTrue[2] /= temp;

  // Compute Index Map
  size_t n = resolutionST[0] * resolutionST[1];

  // -------------------------------------------------------------------------
  // Create Vertices
  // -------------------------------------------------------------------------
  {
    // Compute pixel size in world coordinates
    double pixelWidthWorld = camSize[0] / resolution[0];
    double pixelHeightWorld = camSize[1] / resolution[1];

    // Optimization: precompute half of the camera size to reduce the number of
    // operations in the for loop Include a half pixel offset (-0.5) to center
    // vertices at pixel centers
    double camWidthWorldHalf = 0.5 * camSize[0] - 0.5 * pixelWidthWorld;
    double camHeightWorldHalf = 0.5 * camSize[1] - 0.5 * pixelHeightWorld;

    // Compute depth delta
    double delta = camNearFar[1] - camNearFar[0];

    // Optimization: reorient camera model to bottom left corner to reduce
    // operations in for loop
    double camPosCorner[3] = {camPos[0] - camRight[0] * camWidthWorldHalf
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
      double v = ((double)y) * pixelHeightWorld;
      double vTimesUp[3]
        = {v * camUpTrue[0], v * camUpTrue[1], v * camUpTrue[2]};

      size_t yOffset = y * resolutionST[0];
      for(size_t x = 0; x < resolutionST[0]; x++) {
        size_t pixelIndex = x + yOffset;

        // double d = (double)(depthValues[ pixelIndex ])*delta+camNearFar[0];
        const double depth = ((double)depthValues[pixelIndex]);
        double d = depth > 0.98 ? 0 : depth * delta + camNearFar[0];
        double u = ((double)x) * pixelWidthWorld;

        // compute vertex coordinate
        size_t pointCoordinateOffset = pixelIndex * 3;
        pointCoordinates[pointCoordinateOffset]
          = camPosCorner[0] + u * camRight[0] + vTimesUp[0] + d * camDir[0];
        pointCoordinates[pointCoordinateOffset + 1]
          = camPosCorner[1] + u * camRight[1] + vTimesUp[1] + d * camDir[1];
        pointCoordinates[pointCoordinateOffset + 2]
          = camPosCorner[2] + u * camRight[2] + vTimesUp[2] + d * camDir[2];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Create Triangles
  // -------------------------------------------------------------------------
  {
    auto absDiff = [](dataType a, dataType b) { return a > b ? a - b : b - a; };

    /* Index Structure:
    0 - 1
    | / |
    2 - 3
    */
    size_t xl = resolutionST[0] - 1;
    size_t yl = resolutionST[1] - 1;
    size_t trianglesPerRow = xl * 2;

    double myNan = std::numeric_limits<double>::quiet_NaN();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
    for(size_t y = 0; y < yl; y++) {
      size_t yOffset = y * resolutionST[0];
      size_t triangleIndexOffset = y * trianglesPerRow * 4;
      size_t triangleDistortionOffset = y * trianglesPerRow;

      for(size_t x = 0; x < xl; x++) {
        size_t i0 = x + yOffset;
        size_t i1 = i0 + 1;
        size_t i2 = i0 + resolutionST[0];
        size_t i3 = i2 + 1;

        connectivityList[triangleIndexOffset++] = 3;
        connectivityList[triangleIndexOffset++] = i0;
        connectivityList[triangleIndexOffset++] = i2;
        connectivityList[triangleIndexOffset++] = i1;

        connectivityList[triangleIndexOffset++] = 3;
        connectivityList[triangleIndexOffset++] = i1;
        connectivityList[triangleIndexOffset++] = i2;
        connectivityList[triangleIndexOffset++] = i3;

        double i0Depth = depthValues[i0];
        double i1Depth = depthValues[i1];
        double i2Depth = depthValues[i2];
        double i3Depth = depthValues[i3];

        triangleDistortions[triangleDistortionOffset++]
          = i0Depth > 0.98 || i2Depth > 0.98 || i1Depth > 0.98
              ? myNan
              : std::max(
                absDiff(i0Depth, i1Depth),
                std::max(absDiff(i1Depth, i2Depth), absDiff(i0Depth, i2Depth)));

        triangleDistortions[triangleDistortionOffset++]
          = i1Depth > 0.98 || i2Depth > 0.98 || i3Depth > 0.98
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
                 1, t.getElapsedTime(), this->threadNumber_);

  return 0;
}
