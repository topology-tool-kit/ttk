/// \ingroup vtk
/// \class ttkIcoSphere
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.09.2019
///
/// This filter creates an IcoSphere with a specified radius, center, and number
/// of subdivisions.
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
  double Radius{1};
  bool ComputeNormals{false};

  // single ico sphere
  double Center[3]{0, 0, 0};

  // alternatvely create a sphere at each point
  vtkPoints* Centers{nullptr};

public:
  static ttkIcoSphere *New();
  vtkTypeMacro(ttkIcoSphere, ttkAlgorithm);

  vtkSetMacro(NumberOfSubdivisions, int);
  vtkGetMacro(NumberOfSubdivisions, int);

  vtkSetVector3Macro(Center, double);
  vtkGetVector3Macro(Center, double);

  vtkSetMacro(Radius, double);
  vtkGetMacro(Radius, double);

  vtkSetMacro(ComputeNormals, bool);
  vtkGetMacro(ComputeNormals, bool);

  vtkSetMacro(Centers, vtkPoints*);
  vtkGetMacro(Centers, vtkPoints*);

protected:
  ttkIcoSphere();
  ~ttkIcoSphere();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};