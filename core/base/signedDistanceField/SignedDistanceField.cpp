#include <SignedDistanceField.h>

#include <algorithm>
#include <cmath>

ttk::SignedDistanceField::SignedDistanceField() {
  this->setDebugMsgPrefix("SignedDistanceField");
}

void ttk::SignedDistanceField::findIntersections2D(
  const std::vector<std::array<float, 4>> &segments,
  const std::array<double, 3> &origin,
  std::vector<bool> &vertexIntersection,
  int *const isInterior) const {
  // the segments and the arrays below are expressed in the plane of the grid
  const auto dims = getPlaneDims();
  const long fullRes[3] = {xResolution_, yResolution_, zResolution_};
  const long fullStride[3] = {1, xResolution_, xResolution_ * yResolution_};
  const long res[2] = {fullRes[dims[0]], fullRes[dims[1]]};
  const long stride[2] = {fullStride[dims[0]], fullStride[dims[1]]};
  const double planeOrigin[2] = {origin[dims[0]], origin[dims[1]]};
  const double spacing[2] = {spacing_[dims[0]], spacing_[dims[1]]};
  // first coordinates of the crossings of each row with the curve, for the
  // parity test
  std::vector<std::vector<double>> rowCrossings(res[1]);

  for(const auto &s : segments) {
    // dim is the direction of the grid lines: rows (along the first axis of
    // the plane) then columns (along the second one)
    for(int dim = 0; dim < 2; ++dim) {
      const int other = 1 - dim;
      const double a0 = s[dim], a1 = s[2 + dim];
      const double b0 = s[other], b1 = s[2 + other];
      // segment parallel to these grid lines: its crossings are found by the
      // grid lines of the other direction
      if(b0 == b1)
        continue;

      // indices of the grid lines crossed by the segment
      long kMin
        = std::ceil((std::min(b0, b1) - planeOrigin[other]) / spacing[other]);
      long kMax
        = std::floor((std::max(b0, b1) - planeOrigin[other]) / spacing[other]);
      kMin = std::max(kMin, 0l);
      kMax = std::min(kMax, res[other] - 1);

      for(long k = kMin; k <= kMax; ++k) {
        const double b = planeOrigin[other] + k * spacing[other];
        const double a = a0 + (b - b0) / (b1 - b0) * (a1 - a0);
        // half-open rule: a curve vertex lying on a row is counted once
        if(dim == 0 and ((b0 <= b and b < b1) or (b1 <= b and b < b0)))
          rowCrossings[k].emplace_back(a);
        const double u = (a - planeOrigin[dim]) / spacing[dim];
        if(u < -1e-3 or u > res[dim] - 1 + 1e-3)
          continue;
        // mark both ends of the grid edge crossed by the segment
        long ind = std::floor(u);
        ind = std::max(0l, std::min(ind, res[dim] - 2));
        const long vertexId = ind * stride[dim] + k * stride[other];
        vertexIntersection[vertexId] = true;
        vertexIntersection[vertexId + stride[dim]] = true;
      }
    }
  }

  // a vertex is interior if the curve crosses its row an odd number of times
  // on its left
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(long k = 0; k < res[1]; ++k) {
    auto &crossings = rowCrossings[k];
    std::sort(crossings.begin(), crossings.end());
    size_t noCrossingsLeft = 0;
    for(long i = 0; i < res[0]; ++i) {
      const double a = planeOrigin[0] + i * spacing[0];
      while(noCrossingsLeft < crossings.size()
            and crossings[noCrossingsLeft] < a)
        ++noCrossingsLeft;
      isInterior[i * stride[0] + k * stride[1]] = noCrossingsLeft % 2;
    }
  }
}

float ttk::SignedDistanceField::distanceToSegments2D(
  const float x,
  const float y,
  const std::vector<std::array<float, 4>> &segments) const {
  float minDist2 = std::numeric_limits<float>::max();
  for(const auto &s : segments) {
    const float dx = s[2] - s[0], dy = s[3] - s[1];
    const float len2 = dx * dx + dy * dy;
    float t = (len2 > 0) ? ((x - s[0]) * dx + (y - s[1]) * dy) / len2 : 0;
    t = std::max(0.f, std::min(t, 1.f));
    const float px = s[0] + t * dx - x, py = s[1] + t * dy - y;
    minDist2 = std::min(minDist2, px * px + py * py);
  }
  return std::sqrt(minDist2);
}
