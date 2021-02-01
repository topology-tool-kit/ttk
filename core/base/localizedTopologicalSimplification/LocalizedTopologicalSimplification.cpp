#include <LocalizedTopologicalSimplification.h>

std::string ttk::toFixed(const float &number, const int precision) {
  std::stringstream vFraction;
  vFraction << std::fixed << std::setprecision(precision) << number;
  return vFraction.str();
}
