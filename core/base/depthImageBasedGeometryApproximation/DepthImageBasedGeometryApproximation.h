/// \ingroup base
/// \class ttk::DepthImageBasedGeometryApproximation
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.7.2018
///
/// \brief TTK %depthImageBasedGeometryApproximation processing package.
///
/// %DepthImageBasedGeometryApproximation is a TTK processing package that
/// approximates geomerty based on an input depth image and its corresponding
/// camera parameters.
///
/// Related publication:
/// 'VOIDGA: A View-Approximation Oriented Image Database Generation Approach'
/// Jonas Lukasczyk, Eric Kinner, James Ahrens, Heike Leitte, and Christoph
/// Garth. IEEE 8th Symposium on Large Data Analysis and Visualization (LDAV),
/// 2018.
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <Wrapper.h>

using namespace std;

namespace ttk {

  class DepthImageBasedGeometryApproximation : public Debug {

  public:
    DepthImageBasedGeometryApproximation(){};
    ~DepthImageBasedGeometryApproximation(){};

    // Execute the geometry approximation.
    template <class dataType>
    int execute(dataType *depthValues,
                double *camPos,
                double *camDir,
                double *camUp,
                double *camRes,
                double *camNearFar,
                double *camHeight,

                int subsampling,

                vector<size_t> &indicies,
                vector<tuple<double, double, double>> &vertices,
                vector<tuple<int, int, int>> &triangles,
                vector<double> &triangleDistortions) const;
  };
} // namespace ttk

template <class dataType>
int ttk::DepthImageBasedGeometryApproximation::execute(
  dataType *depthValues,
  double *camPos,
  double *camDir,
  double *camUp,
  double *camRes,
  double *camNearFar,
  double *camHeight,

  int subsampling,

  vector<size_t> &indicies,
  vector<tuple<double, double, double>> &vertices,
  vector<tuple<int, int, int>> &triangles,
  vector<double> &triangleDistortions) const {

  Timer t;
  size_t step = subsampling + 1;

  size_t camResST[2] = {(size_t)camRes[0], (size_t)camRes[1]};

  // -------------------------------------------------------------------------
  // Compute Camera Vectors
  // -------------------------------------------------------------------------

  // Compute camera size
  double camSize[2] = {camRes[0] / camRes[1] * camHeight[0], camHeight[0]};

  // Compute camRight = camDir x CamUp
  double camRight[3] = {camDir[1] * camUp[2] - camDir[2] * camUp[1],
                        camDir[2] * camUp[0] - camDir[0] * camUp[2],
                        camDir[0] * camUp[1] - camDir[1] * camUp[0]};
  double temp = sqrt(camRight[0] * camRight[0] + camRight[1] * camRight[1]
                     + camRight[2] * camRight[2]);
  camRight[0] /= temp;
  camRight[1] /= temp;
  camRight[2] /= temp;

  // Compute true up vector
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
  size_t n = camResST[0] * camResST[1];

  vector<int> pixelIndexVertexIndexMap;
  pixelIndexVertexIndexMap.resize(n);
  size_t numberNewVertices = 0;
  {
    for(size_t i = 0; i < n; i++)
      pixelIndexVertexIndexMap[i]
        = depthValues[i] > 0.99 ? -1 : numberNewVertices++;
  }

  // -------------------------------------------------------------------------
  // Create Vertices
  // -------------------------------------------------------------------------
  {
    // Compute pixel size in world coordinates
    double pixelWidthWorld = camSize[0] / camRes[0];
    double pixelHeightWorld = camSize[1] / camRes[1];

    // Optimization: precompute half of the camera size to reduce the number of
    // operations in the for loop Include a half pixel offset (-0.5) to center
    // vertices at pixel centers
    double camWidthWorldHalf = 0.5 * camSize[0] - 0.5 * pixelWidthWorld;
    double camHeightWorldHalf = 0.5 * camSize[1] - 0.5 * pixelHeightWorld;

    // Make room for new vertices
    vertices.resize(numberNewVertices);
    indicies.resize(numberNewVertices);

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

// Compute vertex positions and parallize over rows
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t y = 0; y < camResST[1]; y += step) {
      double v = ((double)y) * pixelHeightWorld;
      double vTimesUp[3]
        = {v * camUpTrue[0], v * camUpTrue[1], v * camUpTrue[2]};

      size_t yOffset = y * camResST[0];
      for(size_t x = 0; x < camResST[0]; x += step) {
        size_t pixelIndex = x + yOffset;
        int vertexIndex = pixelIndexVertexIndexMap[pixelIndex];
        if(vertexIndex < 0)
          continue;

        // double d = (double)(depthValues[ pixelIndex ])*delta+camNearFar[0];
        double d = ((double)depthValues[pixelIndex]) * delta + camNearFar[0];
        double u = ((double)x) * pixelWidthWorld;
        auto &vertex = vertices[vertexIndex];

        // Store pixel index of vertex
        indicies[vertexIndex] = pixelIndex;

        get<0>(vertex)
          = camPosCorner[0] + u * camRight[0] + vTimesUp[0] + d * camDir[0];
        get<1>(vertex)
          = camPosCorner[1] + u * camRight[1] + vTimesUp[1] + d * camDir[1];
        get<2>(vertex)
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
    size_t xl = camResST[0] - step;
    size_t yl = camResST[1] - step;
    size_t yD = step * camResST[0];

    for(size_t y = 0; y < yl; y += step) {
      for(size_t x = 0; x < xl; x += step) {
        size_t i0 = x + y * camResST[0];
        size_t i1 = i0 + step;
        size_t i2 = i0 + yD;
        size_t i3 = i2 + step;

        int i0Index = pixelIndexVertexIndexMap[i0];
        int i1Index = pixelIndexVertexIndexMap[i1];
        int i2Index = pixelIndexVertexIndexMap[i2];
        int i3Index = pixelIndexVertexIndexMap[i3];

        dataType i0Depth = depthValues[i0];
        dataType i1Depth = depthValues[i1];
        dataType i2Depth = depthValues[i2];
        dataType i3Depth = depthValues[i3];

        if(pixelIndexVertexIndexMap[i1] >= 0
           && pixelIndexVertexIndexMap[i2] >= 0) {
          // Check first triangle
          if(pixelIndexVertexIndexMap[i0] >= 0) {
            triangles.push_back(make_tuple(i0Index, i1Index, i2Index));

            dataType distortion
              = max(absDiff(i0Depth, i1Depth),
                    max(absDiff(i1Depth, i2Depth), absDiff(i0Depth, i2Depth)));

            triangleDistortions.push_back(distortion);
          }

          // Check second triangle
          if(pixelIndexVertexIndexMap[i3] >= 0) {
            triangles.push_back(make_tuple(i1Index, i3Index, i2Index));

            dataType distortion
              = max(absDiff(i3Depth, i1Depth),
                    max(absDiff(i1Depth, i2Depth), absDiff(i3Depth, i2Depth)));

            triangleDistortions.push_back(distortion);
          }
        }
      }
    }
  }

  // Print performance
  {
    stringstream msg;
    msg << "[ttkDepthImageBasedGeometryApproximation] Depth Image ("
        << camResST[0] << "x" << camResST[1] << ":" << step << ") processed in "
        << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  {
    stringstream msg;
    msg << "[ttkDepthImageBasedGeometryApproximation] Generated ("
        << vertices.size() << " vertices) and (" << triangles.size()
        << " triangles)." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}
