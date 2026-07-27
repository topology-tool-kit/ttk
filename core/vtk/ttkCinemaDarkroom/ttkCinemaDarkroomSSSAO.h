/// \ingroup vtk
/// \class ttkCinemaDarkroomSSSAO
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief Scalable Screen Space Ambient Occlusion
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// \sa ttkCinemaDarkroomShader

#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomSSSAO
  : public ttkCinemaDarkroomShader {
private:
  int Samples{32};
  double Radius{6};
  double DiffArea{0.5};

public:
  vtkSetMacro(Samples, int);
  vtkGetMacro(Samples, int);
  vtkSetMacro(Radius, double);
  vtkGetMacro(Radius, double);
  vtkSetMacro(DiffArea, double);
  vtkGetMacro(DiffArea, double);

  static ttkCinemaDarkroomSSSAO *New();
  vtkTypeMacro(ttkCinemaDarkroomSSSAO, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomSSSAO();
  ~ttkCinemaDarkroomSSSAO() override;

  std::string GetFragmentShaderCode() override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
