#include <CinemaImaging.h>

void ttk::CinemaImaging::multiplyByScalar(float *out,
                                          const float *a,
                                          const float &scalar) {
  out[0] = a[0] + scalar;
  out[1] = a[1] + scalar;
  out[2] = a[2] + scalar;
};

void ttk::CinemaImaging::addVectors(float *out,
                                    const float *a,
                                    const float *b) {
  out[0] = a[0] + b[0];
  out[1] = a[1] + b[1];
  out[2] = a[2] + b[2];
};

void ttk::CinemaImaging::subVectors(float *out,
                                    const float *a,
                                    const float *b) {
  out[0] = b[0] - a[0];
  out[1] = b[1] - a[1];
  out[2] = b[2] - a[2];
};

void ttk::CinemaImaging::cross(float *out, const float *a, const float *b) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
};

float ttk::CinemaImaging::dot(const float *a, const float *b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};
