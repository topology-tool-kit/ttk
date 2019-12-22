/// \ingroup vtk
/// \class ttkIcoSphereFromPoint
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.09.2019
///
/// This filter creates an IcoSphere with a specified radius, center, and number
/// of subdivisions. Alternatively, by providing an optional input, the filter
/// will automatically determine the radius and center such that the resulting
/// IcoSphere encapsulates the input object.
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