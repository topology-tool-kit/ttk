/// \ingroup vtk
/// \class ttkCinemaDarkroomCamera
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief TODO
///
/// \param Output vtkPointSet.
///
///
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkAlgorithm.h>
#include <ttkCinemaDarkroomModule.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomCamera : public ttkAlgorithm {

  double Position[3]{0,0,1};
  double Up[3]{0,0,1};
  double Focus[3]{0,0,0};

public:
  static ttkCinemaDarkroomCamera *New();
  vtkTypeMacro(ttkCinemaDarkroomCamera, ttkAlgorithm);

  vtkSetVector3Macro(Position, double);
  vtkGetVector3Macro(Position, double);
  vtkSetVector3Macro(Up, double);
  vtkGetVector3Macro(Up, double);
  vtkSetVector3Macro(Focus, double);
  vtkGetVector3Macro(Focus, double);

protected:
  ttkCinemaDarkroomCamera();
  ~ttkCinemaDarkroomCamera() override;


  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

};