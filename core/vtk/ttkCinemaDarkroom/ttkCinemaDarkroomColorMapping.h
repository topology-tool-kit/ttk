#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class vtkPiecewiseFunction;

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomColorMapping : public ttkCinemaDarkroomShader {
private:

  static const std::vector<std::vector<double>> ColorMaps;

  int ColorMap{0};
  std::string ManualColorMap{""};
  double SolidColor[3]{0,0,0};
  double NANColor[3]{0,0,0};

public:

  vtkSetMacro(ColorMap, int);
  vtkGetMacro(ColorMap, int);
  vtkSetMacro(ManualColorMap, std::string);
  vtkGetMacro(ManualColorMap, std::string);
  vtkSetVector3Macro(NANColor, double);
  vtkGetVector3Macro(NANColor, double);
  vtkSetVector3Macro(SolidColor, double);
  vtkGetVector3Macro(SolidColor, double);

  static ttkCinemaDarkroomColorMapping *New();
  vtkTypeMacro(ttkCinemaDarkroomColorMapping, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomColorMapping();
  ~ttkCinemaDarkroomColorMapping() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};