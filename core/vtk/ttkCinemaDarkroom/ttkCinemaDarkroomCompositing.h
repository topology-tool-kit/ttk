#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomCompositing : public ttkCinemaDarkroomShader {
private:

public:

  static ttkCinemaDarkroomCompositing *New();
  vtkTypeMacro(ttkCinemaDarkroomCompositing, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomCompositing();
  ~ttkCinemaDarkroomCompositing() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
