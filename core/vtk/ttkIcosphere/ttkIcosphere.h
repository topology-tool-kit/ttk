/// \ingroup vtk
/// \class ttkIcosphere
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.09.2019
///
/// This filter creates an Icosphere with a specified radius, center, and number
/// of subdivisions.
///
/// \sa ttk::Icosphere
/// \sa ttk::ttkAlgorithm

#pragma once

// VTK Module
#include <ttkIcosphereModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <Icosphere.h>

class vtkDataArray;

class TTKICOSPHERE_EXPORT ttkIcosphere : public ttkAlgorithm,
                                         protected ttk::Icosphere {
private:
  int NumberOfSubdivisions{0};
  double Radius{1};
  bool ComputeNormals{false};

  // single ico sphere
  double Center[3]{0, 0, 0};

  // alternatvely create a sphere at each point
  vtkDataArray *Centers{nullptr};

public:
  static ttkIcosphere *New();
  vtkTypeMacro(ttkIcosphere, ttkAlgorithm);

  vtkSetMacro(NumberOfSubdivisions, int);
  vtkGetMacro(NumberOfSubdivisions, int);

  vtkSetVector3Macro(Center, double);
  vtkGetVector3Macro(Center, double);

  vtkSetMacro(Radius, double);
  vtkGetMacro(Radius, double);

  vtkSetMacro(ComputeNormals, bool);
  vtkGetMacro(ComputeNormals, bool);

  vtkSetMacro(Centers, vtkDataArray *);
  vtkGetMacro(Centers, vtkDataArray *);

protected:
  ttkIcosphere();
  ~ttkIcosphere();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};