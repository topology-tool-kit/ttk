/// \ingroup vtk
/// \class ttkGaussianPointCloud
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date February 2019.
///
/// \brief TTK VTK-filter that generates a 1D, 2D or 3D point cloud by randomly
/// casting samples from a Gaussian distribution.
///
/// VTK wrapping code for the @GaussianPointCloud package.
///
/// \sa ttk::GaussianPointCloud

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridAlgorithm.h>

// TTK includes
#include <GaussianPointCloud.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkGaussianPointCloud
#else
class ttkGaussianPointCloud
#endif
  : public vtkUnstructuredGridAlgorithm,
    public ttk::Wrapper {

public:
  static ttkGaussianPointCloud *New();
  vtkTypeMacro(ttkGaussianPointCloud, vtkUnstructuredGridAlgorithm)

    vtkSetMacro(Dimension, int);
  vtkGetMacro(Dimension, int);

  vtkSetMacro(NumberOfSamples, int);
  vtkGetMacro(NumberOfSamples, int);

  // default ttk setters
  vtkSetMacro(debugLevel_, int);
  void SetThreads() {
    threadNumber_
      = !UseAllCores ? ThreadNumber : ttk::OsCall::getNumberOfCores();
    Modified();
  }
  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }
  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  int FillOutputPortInformation(int port, vtkInformation *info) override {

    if(!port) {
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    }

    return 1;
  }

protected:
  ttkGaussianPointCloud() {

    Dimension = 2;
    NumberOfSamples = 1000;

    UseAllCores = true;

    SetNumberOfInputPorts(0);
    SetNumberOfOutputPorts(1);
  }

  ~ttkGaussianPointCloud(){};

  bool UseAllCores;
  int ThreadNumber;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int Dimension;
  int NumberOfSamples;
  ttk::GaussianPointCloud gaussianPointCloud_;

  bool needsToAbort() override {
    return GetAbortExecute();
  };
  int updateProgress(const float &progress) override {
    UpdateProgress(progress);
    return 0;
  };
};
