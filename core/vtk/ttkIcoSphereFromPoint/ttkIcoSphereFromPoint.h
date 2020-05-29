/// \ingroup vtk
/// \class ttkIcoSphereFromPoint
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.09.2019
///
/// This filter creates an icosphere with a specified radius, center, and number
/// of subdivisions at each vertex of an input vtkPointset.
///
/// \sa ttk::IcoSphere
/// \sa ttk::ttkAlgorithm

#pragma once

// VTK Module
#include <ttkIcoSphereFromPointModule.h>

// VTK Includes
#include <ttkIcoSphere.h>

class TTKICOSPHEREFROMPOINT_EXPORT ttkIcoSphereFromPoint : public ttkIcoSphere {

private:
  bool CopyPointData{true};

public:
  vtkSetMacro(CopyPointData, bool);
  vtkGetMacro(CopyPointData, bool);

  static ttkIcoSphereFromPoint *New();
  vtkTypeMacro(ttkIcoSphereFromPoint, ttkIcoSphere);

protected:
  ttkIcoSphereFromPoint();
  ~ttkIcoSphereFromPoint();

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};