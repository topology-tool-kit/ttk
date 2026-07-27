/// \ingroup base
/// \class ttk::Ray
/// \authors Rosty Hnatyshyn <rostyslav.hnatyshyn@gmail.com>
/// \date 10.11.2020
///
/// \brief Data structure for a ray.

#pragma once

namespace ttk {
  class Ray {
  public:
    Ray(float *direction, float *origin) {
      m_direction = direction;
      m_origin = origin;
      this->distance = 0;
      this->u = 0;
      this->v = 0;
    }
    float *m_direction;
    float *m_origin;
    float distance;
    float u;
    float v;
  };
} // namespace ttk