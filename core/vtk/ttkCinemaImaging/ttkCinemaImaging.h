/// \ingroup vtk
/// \class ttkCinemaImaging
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.9.2018
///
/// \brief TTK VTK-filter that generates images of a vtkDataSet.
///
/// This filter takes images of a vtkDataObject from positions specified on a
/// vtkPointSet. Each image will be a block of a vtkMultiBlockDataSet where
/// block order corresponds to point order. Each sample point can optionally
/// have vtkDoubleArrays to override the default rendering parameters, i.e, the
/// resolution, focus, clipping planes, and viewport height.
///
/// VTK wrapping code for the @CinemaImaging package.
///
/// \param Input vtkDataObject that will be depicted (vtkDataObject)
/// \param Input vtkPointSet that records the camera sampling locations
/// (vtkPointSet) \param Output vtkMultiBlockDataSet that represents a list of
/// images (vtkMultiBlockDataSet)

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkMultiBlockDataSetAlgorithm.h>

// TTK includes
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkCinemaImaging
#else
class ttkCinemaImaging
#endif
  : public vtkMultiBlockDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkCinemaImaging *New();
  vtkTypeMacro(ttkCinemaImaging, vtkMultiBlockDataSetAlgorithm)

    vtkSetVector2Macro(Resolution, int);
  vtkGetVector2Macro(Resolution, int);

  vtkSetVector2Macro(CamNearFar, double);
  vtkGetVector2Macro(CamNearFar, double);

  vtkSetVector3Macro(CamFocus, double);
  vtkGetVector3Macro(CamFocus, double);

  vtkSetMacro(CamHeight, double);
  vtkGetMacro(CamHeight, double);

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
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
        break;
      case 1:
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
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
  ttkCinemaImaging() {
    int res[2] = {256, 256};
    SetResolution(res);
    double nf[2] = {0.1, 2};
    SetCamNearFar(nf);
    double foc[3] = {0, 0, 0};
    SetCamFocus(foc);
    SetCamHeight(1);

    UseAllCores = false;

    SetNumberOfInputPorts(2);
    SetNumberOfOutputPorts(1);
  }
  ~ttkCinemaImaging(){};

  bool UseAllCores;
  int ThreadNumber;

  int Resolution[2];
  double CamNearFar[2];
  double CamFocus[3];
  double CamHeight;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool needsToAbort() override {
    return GetAbortExecute();
  };
  int updateProgress(const float &progress) override {
    UpdateProgress(progress);
    return 0;
  };
};
