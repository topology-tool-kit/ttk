/// \ingroup vtk
/// \class ttkCinemaDarkroomSSDoF
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief Screen Space Depth of Field
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// \b Related \b Publication:
/// "Efficiently Simulating the Bokeh of Polygonal Apertures in a Post-Process
/// Depth of Field Shader". L. McIntosh, B. E. Riecke and S. DiPaola. Computer
/// Graphics Forum. 2012.
///
/// \sa ttkCinemaDarkroomShader

#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomSSDoF
  : public ttkCinemaDarkroomShader {
private:
  double Radius{0.05};
  double MaxBlur{0.1};
  double Aperture{1.0};
  double FocalDepth{0.4};

public:
  vtkSetMacro(Radius, double);
  vtkGetMacro(Radius, double);
  vtkSetMacro(MaxBlur, double);
  vtkGetMacro(MaxBlur, double);
  vtkSetMacro(Aperture, double);
  vtkGetMacro(Aperture, double);
  vtkSetMacro(FocalDepth, double);
  vtkGetMacro(FocalDepth, double);

  static ttkCinemaDarkroomSSDoF *New();
  vtkTypeMacro(ttkCinemaDarkroomSSDoF, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomSSDoF();
  ~ttkCinemaDarkroomSSDoF() override;

  std::string GetFragmentShaderCode() override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
