/// \ingroup vtk
/// \class ttkCinemaDarkroomFXAA
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief Fast Approximate Anti-Aliasing
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// \b Related \b Publication:
/// "Fast Approximate Anti-Aliasing (FXAA)".
/// Timothy Lottes.
/// NVIDIA Corporation. 2009
///
/// \sa ttkCinemaDarkroomShader

#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomFXAA
  : public ttkCinemaDarkroomShader {
public:
  static ttkCinemaDarkroomFXAA *New();
  vtkTypeMacro(ttkCinemaDarkroomFXAA, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomFXAA();
  ~ttkCinemaDarkroomFXAA() override;

  std::string GetVertexShaderCode() override;
  std::string GetFragmentShaderCode() override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
