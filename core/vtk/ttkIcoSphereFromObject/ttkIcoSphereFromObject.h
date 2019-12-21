/// \ingroup vtk
/// \class ttkIcoSphereFromObject
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
#include <ttkIcoSphereFromObjectModule.h>

// VTK Includes
#include <ttkIcoSphere.h>

class TTKICOSPHEREFROMOBJECT_EXPORT ttkIcoSphereFromObject
  : public ttkIcoSphere {

private:
  float Scale{1};

public:
  static ttkIcoSphereFromObject *New();
  vtkTypeMacro(ttkIcoSphereFromObject, ttkIcoSphere);

  vtkSetMacro(Scale, float);
  vtkGetMacro(Scale, float);

protected:
  ttkIcoSphereFromObject();
  ~ttkIcoSphereFromObject();

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};