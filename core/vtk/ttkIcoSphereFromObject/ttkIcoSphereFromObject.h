/// \ingroup vtk
/// \class ttkIcoSphereFromObject
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.09.2019
///
/// This filter creates an IcoSphere that encapsulates an input data object. The
/// icosphere can optionally be scaled.
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