/// \ingroup vtk
/// \class ttkIcosphereFromObject
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.09.2019
///
/// This filter creates an IcoSphere that encapsulates an input data object. The
/// icosphere can optionally be scaled.
///
/// \sa ttk::IcoSphere
/// \sa ttk::ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/geometryApproximation/">Geometry
///   Approximation example</a> \n

#pragma once

// VTK Module
#include <ttkIcosphereFromObjectModule.h>

// VTK Includes
#include <ttkIcosphere.h>

class TTKICOSPHEREFROMOBJECT_EXPORT ttkIcosphereFromObject
  : public ttkIcosphere {

private:
  float Scale{1};

public:
  static ttkIcosphereFromObject *New();
  vtkTypeMacro(ttkIcosphereFromObject, ttkIcosphere);

  vtkSetMacro(Scale, float);
  vtkGetMacro(Scale, float);

protected:
  ttkIcosphereFromObject();
  ~ttkIcosphereFromObject();

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};