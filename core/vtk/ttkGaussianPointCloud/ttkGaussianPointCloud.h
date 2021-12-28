/// \ingroup vtk
/// \class ttkGaussianPointCloud
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date February 2019.
///
/// \brief TTK VTK-filter that generates a 1D, 2D or 3D point cloud by randomly
/// casting samples from a Gaussian distribution.
///
/// VTK wrapping code for the ttk::GaussianPointCloud package.
///
/// \sa ttk::GaussianPointCloud

#pragma once

// VTK Module
#include <ttkGaussianPointCloudModule.h>

// TTK includes
#include <GaussianPointCloud.h>
#include <ttkAlgorithm.h>

class TTKGAUSSIANPOINTCLOUD_EXPORT ttkGaussianPointCloud
  : public ttkAlgorithm,
    protected ttk::GaussianPointCloud {

public:
  static ttkGaussianPointCloud *New();
  vtkTypeMacro(ttkGaussianPointCloud, ttkAlgorithm);

  vtkSetMacro(Dimension, int);
  vtkGetMacro(Dimension, int);

  vtkSetMacro(NumberOfSamples, int);
  vtkGetMacro(NumberOfSamples, int);

protected:
  ttkGaussianPointCloud();

  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int Dimension{2};
  int NumberOfSamples{1000};
};
