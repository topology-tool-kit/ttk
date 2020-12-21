/// \ingroup vtk
/// \class ttkCinemaDarkroomNoise
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief Noise
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// \sa ttkCinemaDarkroomShader

#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomNoise
  : public ttkCinemaDarkroomShader {
private:
public:
  static ttkCinemaDarkroomNoise *New();
  vtkTypeMacro(ttkCinemaDarkroomNoise, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomNoise();
  ~ttkCinemaDarkroomNoise() override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
