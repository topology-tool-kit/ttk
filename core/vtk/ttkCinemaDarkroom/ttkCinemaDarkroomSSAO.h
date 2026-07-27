/// \ingroup vtk
/// \class ttkCinemaDarkroomSSAO
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief Screen Space Ambient Occlusion
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// \sa ttkCinemaDarkroomShader

#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomSSAO
  : public ttkCinemaDarkroomShader {
private:
  double Radius{6};
  double DiffArea{0.5};

public:
  vtkSetMacro(Radius, double);
  vtkGetMacro(Radius, double);
  vtkSetMacro(DiffArea, double);
  vtkGetMacro(DiffArea, double);

  static ttkCinemaDarkroomSSAO *New();
  vtkTypeMacro(ttkCinemaDarkroomSSAO, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomSSAO();
  ~ttkCinemaDarkroomSSAO() override;

  std::string GetFragmentShaderCode() override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
