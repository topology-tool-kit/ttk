#include <CinemaImagingEmbree.h>

ttk::CinemaImagingEmbree::CinemaImagingEmbree() {
  this->setDebugMsgPrefix("CinemaImaging(Embree)");
}
ttk::CinemaImagingEmbree::~CinemaImagingEmbree() {
}

#if TTK_ENABLE_EMBREE

int ttk::CinemaImagingEmbree::deallocateScene(RTCDevice &device,
                                              RTCScene &scene) const {
  ttk::Timer timer;
  this->printMsg("Deallocating Scene", 0, 0, ttk::debug::LineMode::REPLACE);

  rtcReleaseScene(scene);
  rtcReleaseDevice(device);

  this->printMsg("Deallocating Scene", 1, timer.getElapsedTime());

  return 1;
};

int ttk::CinemaImagingEmbree::initializeDevice(RTCDevice &device) const {
  ttk::Timer timer;
  this->printMsg("Initializing Device", 0, 0, ttk::debug::LineMode::REPLACE);

  device = rtcNewDevice("hugepages=1,threads=1");

  if(!device) {
    this->printErr("Unable to create device");
    this->printErr(std::to_string(rtcGetDeviceError(NULL)));
    return 0;
  }

  auto errorFunction
    = [](void *ttkNotUsed(userPtr), enum RTCError error, const char *str) {
        printf("error %d: %s\n", error, str);
      };

  rtcSetDeviceErrorFunction(device, errorFunction, NULL);

  this->printMsg("Initializing Device", 1, timer.getElapsedTime());

  return 1;
};

int ttk::CinemaImagingEmbree::renderImage(
  float *depthBuffer,
  unsigned int *primitiveIds,
  float *barycentricCoordinates,

  const RTCScene &scene,
  const double resolution[2],
  const double camPos[3],
  const double camDirRaw[3],
  const double camUp[3],
  const double &camFactor,
  const bool &orthographicProjection) const {

  ttk::Timer timer;
  const int resX = resolution[0];
  const int resY = resolution[1];

  this->printMsg("Rendering Image ("
                   + std::string(orthographicProjection ? "O" : "P") + "|"
                   + std::to_string(resX) + "x" + std::to_string(resY) + ")",
                 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

  struct RTCIntersectContext context;
  rtcInitIntersectContext(&context);

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

  struct RTCRayHit rayhit;
  rayhit.ray.dir_x = -1;
  rayhit.ray.dir_y = -1;
  rayhit.ray.dir_z = -1;

  size_t pixelIndex = 0;
  size_t bcIndex = 0;
  float nan = std::numeric_limits<float>::quiet_NaN();

  if(orthographicProjection) {

    // Compute camera size
    const double aspect = resolution[0] / resolution[1];
    const double camSize[2] = {aspect * camFactor, camFactor};

    // Compute pixel size in world coordinates
    const double pixelWidthWorld = camSize[0] / resolution[0];
    const double pixelHeightWorld = camSize[1] / resolution[1];

    // Optimization: precompute half of the camera size to reduce the number of
    // operations in the for loop. Include a half pixel offset (-0.5) to center
    // vertices at pixel centers
    const double camWidthWorldHalf = 0.5 * camSize[0] - 0.5 * pixelWidthWorld;
    const double camHeightWorldHalf = 0.5 * camSize[1] - 0.5 * pixelHeightWorld;

    // Optimization: reorient camera model to bottom left corner to reduce
    // operations in for loop
    const double camPosCorner[3] = {camPos[0] - camRight[0] * camWidthWorldHalf
                                      - camUpTrue[0] * camHeightWorldHalf,
                                    camPos[1] - camRight[1] * camWidthWorldHalf
                                      - camUpTrue[1] * camHeightWorldHalf,
                                    camPos[2] - camRight[2] * camWidthWorldHalf
                                      - camUpTrue[2] * camHeightWorldHalf};

    for(int y = 0; y < resY; y++) {
      double v = ((double)y) * pixelHeightWorld;

      for(int x = 0; x < resX; x++) {
        double u = ((double)x) * pixelWidthWorld;

        // set origin
        rayhit.ray.org_x = camPosCorner[0] + u * camRight[0] + v * camUpTrue[0];
        rayhit.ray.org_y = camPosCorner[1] + u * camRight[1] + v * camUpTrue[1];
        rayhit.ray.org_z = camPosCorner[2] + u * camRight[2] + v * camUpTrue[2];

        // set dir
        rayhit.ray.dir_x = camDir[0];
        rayhit.ray.dir_y = camDir[1];
        rayhit.ray.dir_z = camDir[2];

        // compute hit
        rayhit.ray.tnear = 0.01;
        rayhit.ray.tfar = INFINITY;
        rayhit.ray.mask = 0;
        rayhit.ray.flags = 0;
        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
        rtcIntersect1(scene, &context, &rayhit);

        // write depth
        bool hitPrimitive = rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID;
        if(hitPrimitive) {
          depthBuffer[pixelIndex] = std::max(0.0f, rayhit.ray.tfar);
          primitiveIds[pixelIndex] = rayhit.hit.primID;
          barycentricCoordinates[bcIndex++] = rayhit.hit.u;
          barycentricCoordinates[bcIndex++] = rayhit.hit.v;
        } else {
          depthBuffer[pixelIndex] = nan;
          primitiveIds[pixelIndex] = CinemaImaging::INVALID_ID;
          barycentricCoordinates[bcIndex++] = nan;
          barycentricCoordinates[bcIndex++] = nan;
        }
        pixelIndex++;
      }
    }
  } else {

    double factor = (camFactor / 180.0 * 3.141592653589793) / resolution[0];

    for(int y = 0; y < resY; y++) {
      double v = (y - resY * 0.5) * factor;

      for(int x = 0; x < resX; x++) {
        double u = (x - resX * 0.5) * factor;

        // set origin
        rayhit.ray.org_x = camPos[0];
        rayhit.ray.org_y = camPos[1];
        rayhit.ray.org_z = camPos[2];

        // // set dir
        rayhit.ray.dir_x = camDir[0] + u * camRight[0] + v * camUpTrue[0];
        rayhit.ray.dir_y = camDir[1] + u * camRight[1] + v * camUpTrue[1];
        rayhit.ray.dir_z = camDir[2] + u * camRight[2] + v * camUpTrue[2];

        // compute hit
        rayhit.ray.tnear = 0.01;
        rayhit.ray.tfar = INFINITY;
        rayhit.ray.mask = 0;
        rayhit.ray.flags = 0;
        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
        rtcIntersect1(scene, &context, &rayhit);

        // write depth
        bool hitPrimitive = rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID;
        if(hitPrimitive) {
          depthBuffer[pixelIndex] = std::max(0.0f, rayhit.ray.tfar);
          primitiveIds[pixelIndex] = rayhit.hit.primID;
          barycentricCoordinates[bcIndex++] = rayhit.hit.u;
          barycentricCoordinates[bcIndex++] = rayhit.hit.v;
        } else {
          depthBuffer[pixelIndex] = nan;
          primitiveIds[pixelIndex] = CinemaImaging::INVALID_ID;
          barycentricCoordinates[bcIndex++] = nan;
          barycentricCoordinates[bcIndex++] = nan;
        }
        pixelIndex++;
      }
    }
  }

  this->printMsg("Rendering Image ("
                   + std::string(orthographicProjection ? "O" : "P") + "|"
                   + std::to_string(resX) + "x" + std::to_string(resY) + ")",
                 1, timer.getElapsedTime(), this->threadNumber_);

  return 1;
};

#endif
