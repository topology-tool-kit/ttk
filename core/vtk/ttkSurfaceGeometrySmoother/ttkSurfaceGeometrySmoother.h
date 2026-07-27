/// \ingroup vtk
/// \class ttkSurfaceGeometrySmoother
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2022.
///
/// \brief TTK VTK-filter for smoothing meshes on surfaces.
///
/// ttk::GeometrySmoother with a twist!
/// This class smoothes and projects a 1D or a 2D mesh onto a 2D
/// closed triangulated surface.
///
/// \param Input vtkDataSet
/// \param Output vtkDataSet
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The optional mask array can be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the mask array)
/// \note: To use this optional array, `ForceInputMaskScalarField` needs to be
/// enabled with the setter `setForceInputMaskScalarField()'.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \sa vtkGeometrySmoother
/// \sa ttk::SurfaceGeometrySmoother
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_casting/">Persistent
///   Generators Casting example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_fertility/">Persistent
///   Generators Fertility example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_skull/">Persistent
///   Generators Skull example</a> \n
///

#pragma once

// VTK Module
#include <ttkSurfaceGeometrySmootherModule.h>

// ttk code includes
#include <SurfaceGeometrySmoother.h>
#include <ttkAlgorithm.h>

class TTKSURFACEGEOMETRYSMOOTHER_EXPORT ttkSurfaceGeometrySmoother
  : public ttkAlgorithm,
    protected ttk::SurfaceGeometrySmoother {

public:
  static ttkSurfaceGeometrySmoother *New();

  vtkTypeMacro(ttkSurfaceGeometrySmoother, ttkAlgorithm);

  vtkSetMacro(NumberOfIterations, int);
  vtkGetMacro(NumberOfIterations, int);

  vtkSetMacro(UseMaskScalarField, bool);
  vtkGetMacro(UseMaskScalarField, bool);

  vtkSetMacro(ForceInputMaskScalarField, bool);
  vtkGetMacro(ForceInputMaskScalarField, bool);

  vtkSetMacro(ForceIdentifiersField, bool);
  vtkGetMacro(ForceIdentifiersField, bool);

protected:
  ttkSurfaceGeometrySmoother();
  ~ttkSurfaceGeometrySmoother() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int NumberOfIterations{1};
  bool UseMaskScalarField{true};
  bool ForceInputMaskScalarField{false};
  bool ForceIdentifiersField{false};
};
