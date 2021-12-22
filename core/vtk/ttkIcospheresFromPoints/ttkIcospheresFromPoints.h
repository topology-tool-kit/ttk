/// \ingroup vtk
/// \class ttkIcospheresFromPoints
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.09.2019
///
/// This filter creates an icosphere with a specified radius, center, and number
/// of subdivisions at each vertex of an input vtkPointSet.
///
/// \sa ttk::IcoSphere
/// \sa ttk::ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/dragon/">Dragon
/// example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/harmonicSkeleton/">
///   Harmonic Skeleton example</a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/morseMolecule/">
/// Morse molecule example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n
///

#pragma once

// VTK Module
#include <ttkIcospheresFromPointsModule.h>

// VTK Includes
#include <ttkIcosphere.h>

class TTKICOSPHERESFROMPOINTS_EXPORT ttkIcospheresFromPoints
  : public ttkIcosphere {

private:
  bool CopyPointData{true};

public:
  vtkSetMacro(CopyPointData, bool);
  vtkGetMacro(CopyPointData, bool);

  static ttkIcospheresFromPoints *New();
  vtkTypeMacro(ttkIcospheresFromPoints, ttkIcosphere);

protected:
  ttkIcospheresFromPoints();
  ~ttkIcospheresFromPoints();

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
