#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomIBS : public ttkCinemaDarkroomShader {
private:

  double Scale{1.0};
  double Strength{1.0};
  double Luminance{1.0};

public:

  vtkSetMacro(Scale,double);
  vtkGetMacro(Scale,double);
  vtkSetMacro(Strength,double);
  vtkGetMacro(Strength,double);
  vtkSetMacro(Luminance,double);
  vtkGetMacro(Luminance,double);

  static ttkCinemaDarkroomIBS *New();
  vtkTypeMacro(ttkCinemaDarkroomIBS, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomIBS();
  ~ttkCinemaDarkroomIBS() override;

  std::string GetFragmentShaderCode() override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
