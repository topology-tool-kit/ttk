/// \ingroup vtk
/// \class ttkDepthImageBasedGeometryApproximation
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.07.2018
///
/// \brief TTK VTK-filter that approximates the geomerty that is depicted by a
/// set of depth images.
///
/// VTK wrapping code for the @DepthImageBasedGeometryApproximation package.
///
/// This filter approximates the geometry that is depicted by a set of depth
/// images.
///
/// Related publication:
/// 'VOIDGA: A View-Approximation Oriented Image Database Generation Approach'
/// Jonas Lukasczyk, Eric Kinner, James Ahrens, Heike Leitte, and Christoph
/// Garth. IEEE 8th Symposium on Large Data Analysis and Visualization (LDAV),
/// 2018.
///
/// \param Input A vtkMultiBlockDataSet containing a set of depth images
/// represented by vtkImagedata objects. (vtkMultiBlockDataSet) \param
/// Subsampling The factor that controls the sampling rate (0: no Subsampling,
/// 1: skip every second sample, 2: skip every second and third sample...)
/// \param Output A set of unstructured grids where each grid corresponds to a
/// depth image (vtkMultiBlockDataSet)
///
/// \sa ttk::DepthImageBasedGeometryApproximation

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkMultiBlockDataSetAlgorithm.h>

// TTK includes
#include <DepthImageBasedGeometryApproximation.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkDepthImageBasedGeometryApproximation
#else
class ttkDepthImageBasedGeometryApproximation
#endif
  : public vtkMultiBlockDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkDepthImageBasedGeometryApproximation *New();
  vtkTypeMacro(ttkDepthImageBasedGeometryApproximation,
               vtkMultiBlockDataSetAlgorithm)

    vtkSetMacro(Subsampling, int);
  vtkGetMacro(Subsampling, int);

  vtkSetMacro(DepthScalarField, string);
  vtkGetMacro(DepthScalarField, string);

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

  int FillInputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
        break;
      default:
        return 0;
    }
    return 1;
  }

  int FillOutputPortInformation(int port, vtkInformation *info) override {
    switch(port) {
      case 0:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
        break;
      default:
        return 0;
    }
    return 1;
  }

protected:
  ttkDepthImageBasedGeometryApproximation() {
    Subsampling = 0;

    UseAllCores = false;
    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }
  ~ttkDepthImageBasedGeometryApproximation(){};

  bool UseAllCores;
  int ThreadNumber;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int Subsampling;
  string DepthScalarField;
  ttk::DepthImageBasedGeometryApproximation
    depthImageBasedGeometryApproximation_;

  bool needsToAbort() override {
    return GetAbortExecute();
  };
  int updateProgress(const float &progress) override {
    UpdateProgress(progress);
    return 0;
  };
};
