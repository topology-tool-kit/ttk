/// \ingroup base
/// \class ttk::OsCall
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2013.
///
/// \brief Os-specifics.

#ifndef _OS_H
#define _OS_H

#ifdef _WIN32
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#define drand48() (double(rand()) / RAND_MAX)
//  #define               isnan(x)      _isnan(x)
#define srand48(seed) srand(seed)
#endif // _WIN32

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector>

//#define SINGLE_PRECISION

#ifdef SINGLE_PRECISION
#define REAL_TYPE float
#define REAL_TYPE_STRING "float"
#define REAL_MAX FLT_MAX
#else
#define REAL_TYPE double
#define REAL_TYPE_STRING "double"
#define REAL_MAX DBL_MAX
#endif

#define DBL_SIGNIFICANT_DIGITS 14
#define FLT_SIGNIFICANT_DIGITS 7

#ifdef SINGLE_PRECISION
#define REAL_SIGNIFICANT_DIGITS FLT_SIGNIFICANT_DIGITS
#else
#define REAL_SIGNIFICANT_DIGITS DBL_SIGNIFICANT_DIGITS
#endif

#ifndef __APPLE__
#define M_PI 3.14159265358979323846
#endif

namespace ttk {

#ifdef SINGLE_PRECISION
  using real = float;
#else
  using real = double;
#endif

  class OsCall {
  public:
    static int getCurrentDirectory(std::string &directoryPath);

    static float getMemoryInstantUsage();

    static int getNumberOfCores();

    static std::vector<std::string>
      listFilesInDirectory(const std::string &directoryName,
                           const std::string &extension);

    static int mkDir(const std::string &directoryName);

    static int nearbyint(const double &x);

    static int rmDir(const std::string &directoryName);

    static int rmFile(const std::string &fileName);

    int static roundToNearestInt(const double &val);
  };

  inline int OsCall::roundToNearestInt(const double &val) {
    const double upperBound = ceil(val);
    const double lowerBound = floor(val);

    if(upperBound - val <= val - lowerBound) {
      return (int)upperBound;
    } else {
      return (int)lowerBound;
    }
  }

  class Memory {

  public:
    Memory() {
      initialMemory_ = OsCall::getMemoryInstantUsage();
    }

    inline float getInitialMemoryUsage() {
      return initialMemory_;
    }

    inline float getInstantUsage() {
      return OsCall::getMemoryInstantUsage();
    }

    inline float getElapsedUsage() {
      return OsCall::getMemoryInstantUsage() - initialMemory_;
    }

  protected:
    float initialMemory_;
  };

} // namespace ttk

#endif
