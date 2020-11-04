#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomFXAA : public ttkCinemaDarkroomShader {
private:

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
