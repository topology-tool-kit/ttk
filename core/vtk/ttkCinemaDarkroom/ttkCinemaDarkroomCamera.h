/// \ingroup vtk
/// \class ttkCinemaDarkroomCamera
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief This source generates a Cinema Darkroom Camera.
///
/// \param Output vtkUnstructuredGrid.
///
/// This source generates a single vertex with point data to represent a camera
/// that can be used for Cinema Darkroom rendering. The source can also be
/// synchronized with the camera of an active view port.
///
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkAlgorithm.h>
#include <ttkCinemaDarkroomModule.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomCamera : public ttkAlgorithm {

  double Position[3]{0, 0, 0};
  double Up[3]{0, 0, 1};
  double FocalPoint[3]{0, 0, 0};

public:
  static ttkCinemaDarkroomCamera *New();
  vtkTypeMacro(ttkCinemaDarkroomCamera, ttkAlgorithm);

  vtkSetVector3Macro(Position, double);
  vtkGetVector3Macro(Position, double);
  vtkSetVector3Macro(Up, double);
  vtkGetVector3Macro(Up, double);
  vtkSetVector3Macro(FocalPoint, double);
  vtkGetVector3Macro(FocalPoint, double);

  int SyncWithParaViewCamera();

protected:
  ttkCinemaDarkroomCamera();
  ~ttkCinemaDarkroomCamera() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};