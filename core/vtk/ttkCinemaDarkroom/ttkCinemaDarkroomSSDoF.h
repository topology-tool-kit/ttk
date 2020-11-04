#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomSSDoF : public ttkCinemaDarkroomShader {
private:

  double Radius{0.05};
  double MaxBlur{0.1};
  double Aperture{1.0};
  double Distance{0.4};

public:

  vtkSetMacro(Radius,double);
  vtkGetMacro(Radius,double);
  vtkSetMacro(MaxBlur,double);
  vtkGetMacro(MaxBlur,double);
  vtkSetMacro(Aperture,double);
  vtkGetMacro(Aperture,double);
  vtkSetMacro(Distance,double);
  vtkGetMacro(Distance,double);

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
