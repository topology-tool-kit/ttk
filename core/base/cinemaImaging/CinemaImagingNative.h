/// \ingroup base
/// \class ttk::CinemaImagingNative
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \author Rosty Hnatyshyn <rostyslav.hnatyshyn@gmail.com>
/// \date 10.11.2020
///
/// \brief Native renderer that uses a bounding volume hierarchy for accelerated
/// raycasting.

#pragma once

#include "CinemaImaging.h"
#include <BoundingVolumeHierarchy.h>
#include <Ray.h>

namespace ttk {
  class CinemaImagingNative : public CinemaImaging {
  public:
    CinemaImagingNative() {
      this->setDebugMsgPrefix("CinemaImaging(Native)");
    }
    ~CinemaImagingNative() = default;

    template <typename IT>
    int renderImage(float *depthBuffer,
                    unsigned int *primitiveIds,
                    float *barycentricCoordinates,

                    const size_t &nVertices,
                    const float *vertexCoords,
                    const size_t &nTriangles,
                    const IT *connectivityList,

                    const BoundingVolumeHierarchy<IT> &bvh,

                    const double resolution[2],
                    const double camPos[3],
                    const double camDirRaw[3],
                    const double camUp[3],
                    const double &camHeight,
                    const bool &orthographicProjection,
                    const double &viewAngle) const;
  };

} // namespace ttk

template <typename IT>
int ttk::CinemaImagingNative::renderImage(
  float *depthBuffer,
  unsigned int *primitiveIds,
  float *barycentricCoordinates,
  const size_t &ttkNotUsed(nVertices),
  const float *vertexCoords,
  const size_t &ttkNotUsed(nTriangles),
  const IT *connectivityList,
  const BoundingVolumeHierarchy<IT> &bvh,
  const double resolution[2],
  const double camPos[3],
  const double camDirRaw[3],
  const double camUp[3],
  const double &camHeight,
  const bool &orthographicProjection,
  const double &viewAngle) const {
  ttk::Timer timer;
  int resX = resolution[0];
  int resY = resolution[1];

  this->printMsg("Rendering Image ("
                   + std::string(orthographicProjection ? "O" : "P") + "|"
                   + std::to_string(resX) + "x" + std::to_string(resY) + ")",
                 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

  // Compute camera size
  const double aspect = resolution[0] / resolution[1];
  const double camSize[2] = {aspect * camHeight, camHeight};

  const auto normalize = [](double out[3], const double in[3]) {
    double temp = sqrt(in[0] * in[0] + in[1] * in[1] + in[2] * in[2]);
    out[0] = in[0] / temp;
    out[1] = in[1] / temp;
    out[2] = in[2] / temp;
  };

  double camDir[3]{0, 0, 0};
  normalize(camDir, camDirRaw);

  // Compute camRight = camDir x CamUp
  double camRight[3]{camDir[1] * camUp[2] - camDir[2] * camUp[1],
                     camDir[2] * camUp[0] - camDir[0] * camUp[2],
                     camDir[0] * camUp[1] - camDir[1] * camUp[0]};
  normalize(camRight, camRight);

  // Compute true up std::vector
  double camUpTrue[3]{camDir[1] * (-camRight[2]) - camDir[2] * (-camRight[1]),
                      camDir[2] * (-camRight[0]) - camDir[0] * (-camRight[2]),
                      camDir[0] * (-camRight[1]) - camDir[1] * (-camRight[0])};
  normalize(camUpTrue, camUpTrue);

  // Compute pixel size in world coordinates
  double pixelWidthWorld = camSize[0] / resolution[0];
  double pixelHeightWorld = camSize[1] / resolution[1];

  // Optimization: precompute half of the camera size to reduce the number of
  // operations in the for loop. Include a half pixel offset (-0.5) to center
  // vertices at pixel centers
  double camWidthWorldHalf = 0.5 * camSize[0] - 0.5 * pixelWidthWorld;
  double camHeightWorldHalf = 0.5 * camSize[1] - 0.5 * pixelHeightWorld;

  // Optimization: reorient camera model to bottom left corner to reduce
  // operations in for loop
  double camPosCorner[3] = {camPos[0] - camRight[0] * camWidthWorldHalf
                              - camUpTrue[0] * camHeightWorldHalf,
                            camPos[1] - camRight[1] * camWidthWorldHalf
                              - camUpTrue[1] * camHeightWorldHalf,
                            camPos[2] - camRight[2] * camWidthWorldHalf
                              - camUpTrue[2] * camHeightWorldHalf};

  float nan = std::numeric_limits<float>::quiet_NaN();
  if(orthographicProjection) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
    for(int y = 0; y < resY; y++) {
      double v = ((double)y) * pixelHeightWorld;

      size_t pixelIndex = y * resX;
      size_t bcIndex = 2 * pixelIndex;

      for(int x = 0; x < resX; x++) {
        double u = ((double)x) * pixelWidthWorld;

        depthBuffer[pixelIndex] = nan;
        primitiveIds[pixelIndex] = CinemaImaging::INVALID_ID;
        barycentricCoordinates[bcIndex] = nan;
        barycentricCoordinates[bcIndex + 1] = nan;

        // set origin
        float org_x = camPosCorner[0] + u * camRight[0] + v * camUpTrue[0];
        float org_y = camPosCorner[1] + u * camRight[1] + v * camUpTrue[1];
        float org_z = camPosCorner[2] + u * camRight[2] + v * camUpTrue[2];

        float ray_origin[3] = {org_x, org_y, org_z};

        // set dir
        float dir_x = camDir[0];
        float dir_y = camDir[1];
        float dir_z = camDir[2];

        float ray_dir[3] = {dir_x, dir_y, dir_z};

        Ray ray(ray_dir, ray_origin);
        bool wasHit = false;
        int triIdx;
        float distance;
        wasHit = bvh.intersect(
          ray, connectivityList, vertexCoords, &triIdx, &distance);
        if(wasHit) {
          depthBuffer[pixelIndex] = distance;
          primitiveIds[pixelIndex] = triIdx;
          barycentricCoordinates[bcIndex] = ray.u;
          barycentricCoordinates[bcIndex + 1] = ray.v;
        }
        pixelIndex++;
        bcIndex += 2;
      }
    }
  } else {
    double factor = (viewAngle / 180.0 * 3.141592653589793) / resolution[0];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
    for(int y = 0; y < resY; y++) {
      double v = (y - resY * 0.5) * factor;
      size_t pixelIndex = y * resX;
      size_t bcIndex = 2 * pixelIndex;

      for(int x = 0; x < resX; x++) {
        double u = (x - resX * 0.5) * factor;

        depthBuffer[pixelIndex] = nan;
        primitiveIds[pixelIndex] = CinemaImaging::INVALID_ID;
        barycentricCoordinates[bcIndex] = nan;
        barycentricCoordinates[bcIndex + 1] = nan;

        // set origin
        float org_x = camPos[0];
        float org_y = camPos[1];
        float org_z = camPos[2];

        float ray_origin[3] = {org_x, org_y, org_z};
        // set dir
        float dir_x = camDir[0] + u * camRight[0] + v * camUpTrue[0];
        float dir_y = camDir[1] + u * camRight[1] + v * camUpTrue[1];
        float dir_z = camDir[2] + u * camRight[2] + v * camUpTrue[2];

        float ray_dir[3] = {dir_x, dir_y, dir_z};

        Ray ray(ray_dir, ray_origin);
        bool wasHit = false;
        int triIdx;
        float distance;
        wasHit = bvh.intersect(
          ray, connectivityList, vertexCoords, &triIdx, &distance);
        if(wasHit) {
          depthBuffer[pixelIndex] = distance;
          primitiveIds[pixelIndex] = triIdx;
          barycentricCoordinates[bcIndex] = ray.u;
          barycentricCoordinates[bcIndex + 1] = ray.v;
        }
        pixelIndex++;
        bcIndex += 2;
      }
    }
  }
  this->printMsg("Rendering Image ("
                   + std::string(orthographicProjection ? "O" : "P") + "|"
                   + std::to_string(resX) + "x" + std::to_string(resY) + ")",

                 1, timer.getElapsedTime(), this->threadNumber_);

  return 1;
}
