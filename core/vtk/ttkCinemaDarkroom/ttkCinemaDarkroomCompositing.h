/// \ingroup vtk
/// \class ttkCinemaDarkroomCompositing
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief Composites multiple vtkImageData objects base on depth values.
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// Composites multiple vtkImageData objects passed to the first input port base
/// on depth point data arrays specified via
/// SetInputArrayToProcess(0,0,0,0,"DepthArrayName").
///
/// \b Related \b Publication:
/// "Cinema Database Specification - Dietrich Release v1.2".
/// D. Rogers, J. Woodring, J. Ahrens, J. Patchett, and J. Lukasczyk.
/// Technical Report LA-UR-17-25072, Los Alamos National Laboratory,
/// 2018.
///
/// \sa ttkCinemaDarkroomShader

#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomCompositing
  : public ttkCinemaDarkroomShader {
private:
public:
  static ttkCinemaDarkroomCompositing *New();
  vtkTypeMacro(ttkCinemaDarkroomCompositing, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomCompositing();
  ~ttkCinemaDarkroomCompositing() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
