/// \ingroup vtk
/// \class ttkIcoSphere
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
#include <ttkIcoSphereModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <IcoSphere.h>

class TTKICOSPHERE_EXPORT ttkIcoSphere : public ttkAlgorithm,
                                         protected ttk::IcoSphere {
private:
  int NumberOfSubdivisions{0};
  float Radius{1};

  // single ico sphere
  float Center[3]{0, 0, 0};

  // alternatively: multiple ico spheres
  int NumberOfIcoSpheres{1};
  float *Centers{nullptr};

public:
  static ttkIcoSphere *New();
  vtkTypeMacro(ttkIcoSphere, ttkAlgorithm);

  vtkSetMacro(NumberOfSubdivisions, int);
  vtkGetMacro(NumberOfSubdivisions, int);

  vtkSetVector3Macro(Center, float);
  vtkGetVector3Macro(Center, float);

  vtkSetMacro(Radius, float);
  vtkGetMacro(Radius, float);

  vtkSetMacro(NumberOfIcoSpheres, int);
  vtkGetMacro(NumberOfIcoSpheres, int);

  vtkSetMacro(Centers, float *);
  vtkGetMacro(Centers, float *);

protected:
  ttkIcoSphere();
  ~ttkIcoSphere();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};