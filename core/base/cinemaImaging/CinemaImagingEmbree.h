/// \ingroup base
/// \class ttk::CinemaImagingEmbree
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.05.2020
///
/// \brief TTK %CinemaImagingEmbree processing package.
///
/// %CinemaImagingEmbree is a TTK processing package that

#pragma once

#include <CinemaImaging.h>

#if TTK_ENABLE_EMBREE
#include <embree3/rtcore.h>
#include <limits>
#include <string>
#endif

namespace ttk {

  class CinemaImagingEmbree : public CinemaImaging {
  public:
    CinemaImagingEmbree();
    ~CinemaImagingEmbree();

#if TTK_ENABLE_EMBREE

    int initializeDevice(RTCDevice &device) const;

    template <typename IT>
    int initializeScene(RTCScene &scene,

                        const RTCDevice &device,
                        const size_t &nVertices,
                        const float *vertexCoords,
                        const size_t &nTriangles,
                        const IT *connectivityList) const;

    int renderImage(float *depthBuffer,
                    unsigned int *primitiveIds,
                    float *barycentricCoordinates,

                    const RTCScene &scene,
                    const double resolution[2],
                    const double camCenter[3],
                    const double camDir[3],
                    const double camUp[3],
                    const double &camFactor, // either height or angle
                    const bool &orthographicProjection = true) const;

    int deallocateScene(RTCDevice &device, RTCScene &scene) const;

#endif
  };
} // namespace ttk

#if TTK_ENABLE_EMBREE

template <typename IT>
int ttk::CinemaImagingEmbree::initializeScene(
  RTCScene &scene,

  const RTCDevice &device,
  const size_t &nVertices,
  const float *vertexCoords,
  const size_t &nTriangles,
  const IT *connectivityList) const {
  ttk::Timer timer;
  this->printMsg("Initializing Scene (#v:" + std::to_string(nVertices)
                   + "|#t:" + std::to_string(nTriangles) + ")",
                 0, 0, ttk::debug::LineMode::REPLACE);

  scene = rtcNewScene(device);

  RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

  // vertices
  {
    rtcSetSharedGeometryBuffer(mesh, RTC_BUFFER_TYPE_VERTEX, 0,
                               RTC_FORMAT_FLOAT3, (const void *)vertexCoords, 0,
                               3 * sizeof(float), nVertices);
  }

  // triangles
  {
    // unfortunately embree does not support signed integer based indexing
    //   rtcSetSharedGeometryBuffer(
    //       mesh,
    //       RTC_BUFFER_TYPE_INDEX,
    //       0,
    //       RTC_FORMAT_LLONG3,
    //       static_cast<const void*>(connectivityList),
    //       0,
    //       3*sizeof(long long),
    //       nTriangles
    //   );

    unsigned int *indices = (unsigned int *)rtcSetNewGeometryBuffer(
      mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
      3 * sizeof(unsigned int), nTriangles);

    for(size_t t = 0, tn = nTriangles * 3; t < tn; t++)
      indices[t] = (unsigned int)connectivityList[t];

    rtcCommitGeometry(mesh);
    rtcAttachGeometry(scene, mesh);
    rtcReleaseGeometry(mesh);
  }

  rtcCommitScene(scene);

  this->printMsg("Initializing Scene (#v:" + std::to_string(nVertices)
                   + "|#t:" + std::to_string(nTriangles) + ")",
                 1, timer.getElapsedTime());
  return 1;
};

#endif
