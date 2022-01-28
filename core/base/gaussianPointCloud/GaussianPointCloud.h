/// \ingroup base
/// \class ttk::GaussianPointCloud
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date February 2019.
///
/// \brief TTK %gaussianPointCloud processing package that generates a 1D, 2D or
/// 3D point cloud by randomly casting samples from a Gaussian distribution.
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <Debug.h>
#include <random>

namespace ttk {

  class GaussianPointCloud : virtual public Debug {
  public:
    GaussianPointCloud();

    template <class dataType>
    int castSample(const int &dimension,
                   dataType &x,
                   dataType &y,
                   dataType &z) const;

    // Execute the sampling
    template <class dataType>
    int generate(const int &dimension,
                 const int &numberOfSamples,
                 dataType *const outputData) const;
  };
} // namespace ttk

template <class dataType>
int ttk::GaussianPointCloud::castSample(const int &dimension,
                                        dataType &x,
                                        dataType &y,
                                        dataType &z) const {

  std::random_device rd{};
  std::mt19937 gen{rd()};

  dataType u, v, w;

  std::uniform_real_distribution<dataType> uniformDistribution(-1, 1);
  u = uniformDistribution(gen);
  v = (dimension >= 2 ? uniformDistribution(gen) : 0);
  w = (dimension == 3 ? uniformDistribution(gen) : 0);

  dataType norm = std::sqrt(u * u + v * v + w * w);

  std::normal_distribution<dataType> normalDistribution{0, 1};

  dataType gaussianScale = normalDistribution(gen);

  x = (u / norm) * gaussianScale;
  y = (v / norm) * gaussianScale;
  z = (w / norm) * gaussianScale;

  return 0;
}

template <class dataType>
int ttk::GaussianPointCloud::generate(const int &dimension,
                                      const int &numberOfSamples,
                                      dataType *const outputData) const {

  Timer t;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < numberOfSamples; i++) {
    castSample<dataType>(dimension, outputData[3 * i], outputData[3 * i + 1],
                         outputData[3 * i + 2]);
  }

  this->printMsg(std::vector<std::vector<std::string>>{
    {"#Samples", std::to_string(numberOfSamples)}});
  this->printMsg("Samples generated in " + std::to_string(dimension) + "D", 1.0,
                 t.getElapsedTime(), this->threadNumber_);

  return 0;
}
