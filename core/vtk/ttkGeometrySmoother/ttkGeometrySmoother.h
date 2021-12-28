/// \ingroup vtk
/// \class ttkGeometrySmoother
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter for geometry smoothing.
///
/// This filter is a dummy example for the development of TTK packages. It
/// smooths an input mesh by average the vertex locations on the link of each
/// vertex.
///
/// \param Input Input mesh (vtkDataSet)
/// \param Output Output mesh (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkScalarFieldSmoother
/// \sa ttk::ScalarFieldSmoother
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearning/">1-Manifold
///   Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearningCircles/">1-Manifold
///   Learning Circles example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/2manifoldLearning/">
///   2-Manifold Learning example</a> \n
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

#include <ScalarFieldSmoother.h>
#include <ttkAlgorithm.h>
#include <ttkGeometrySmootherModule.h>

class TTKGEOMETRYSMOOTHER_EXPORT ttkGeometrySmoother
  : public ttkAlgorithm,
    protected ttk::ScalarFieldSmoother {

private:
  int NumberOfIterations{1};
  int MaskIdentifier{0};
  bool ForceInputMaskScalarField{false};

public:
  vtkSetMacro(NumberOfIterations, int);
  vtkGetMacro(NumberOfIterations, int);

  vtkSetMacro(MaskIdentifier, int);
  vtkGetMacro(MaskIdentifier, int);

  vtkSetMacro(ForceInputMaskScalarField, bool);
  vtkGetMacro(ForceInputMaskScalarField, bool);

  vtkTypeMacro(ttkGeometrySmoother, ttkAlgorithm);
  static ttkGeometrySmoother *New();

protected:
  ttkGeometrySmoother();
  ~ttkGeometrySmoother();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
