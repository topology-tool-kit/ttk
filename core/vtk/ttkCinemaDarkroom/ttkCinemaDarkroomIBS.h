/// \ingroup vtk
/// \class ttkCinemaDarkroomIBS
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief Image Based Shading
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// \b Related \b Publication:
/// "Dynamic Nested Tracking Graphs".
/// Jonas Lukasczyk, Christoph Garth, Gunther H. Weber, Tim Biedert, Ross
/// Maciejewski, Heike Leitte. IEEE Transactions on Visualization and Computer
/// Graphics. 2019
///
/// \sa ttkCinemaDarkroomShader

#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomIBS
  : public ttkCinemaDarkroomShader {
private:
  double Strength{1.0};
  double Luminance{1.0};
  double Ambient{0.2};

public:
  vtkSetMacro(Strength, double);
  vtkGetMacro(Strength, double);
  vtkSetMacro(Luminance, double);
  vtkGetMacro(Luminance, double);
  vtkSetMacro(Ambient, double);
  vtkGetMacro(Ambient, double);

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
