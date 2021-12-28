/// \ingroup vtk
/// \class ttkContinuousScatterPlot
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date March 2016
///
/// \brief TTK VTK-filter that computes the continuous scatterplot of bivariate
/// volumetric data.
///
/// This filter produces a 2D vtkUnstructuredGrid with a scalar field (named
/// "Density") representing the continuous scatter plot (attached to the
/// 2D geometry as point data). A point mask is also attached to the 2D
/// geometry as point data.
///
/// The components of the input bivariate data must be specified as independent
/// scalar fields attached to the input geometry as point data.
///
/// \param Input Input bivariate volumetric data-set, either regular grids or
/// tetrahedral meshes (vtkDataSet)
/// \param Output Output 2D continuous scatter plot (vtkUnstructuredGrid)
///
/// This module respects the following convention regarding the order of the
/// input arrays to process (SetInputArrayToProcess()):
/// \param idx 0: first scalar array.
/// \param idx 1: second scalar array.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Continuous Scatterplots" \n
/// Sven Bachthaler, Daniel Weiskopf \n
/// Proc. of IEEE VIS 2008.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2008.
///
/// \sa ttk::ContinuousScatterPlot
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/builtInExample2/">
///   Builtin example 2</a> \n

#pragma once

// VTK Module
#include <ttkContinuousScatterPlotModule.h>

// ttk baseCode includes
#include <ContinuousScatterPlot.h>
#include <ttkAlgorithm.h>

class vtkDataArray;

class TTKCONTINUOUSSCATTERPLOT_EXPORT ttkContinuousScatterPlot
  : public ttkAlgorithm,
    protected ttk::ContinuousScatterPlot {

public:
  static ttkContinuousScatterPlot *New();

  vtkTypeMacro(ttkContinuousScatterPlot, ttkAlgorithm);

  vtkSetMacro(WithDummyValue, bool);
  vtkGetMacro(WithDummyValue, bool);

  vtkSetMacro(DummyValue, double);
  vtkGetMacro(DummyValue, double);

  vtkSetMacro(ProjectImageSupport, bool);
  vtkGetMacro(ProjectImageSupport, bool);

  void SetScatterplotResolution(int N, int M) {
    ScatterplotResolution[0] = N;
    ScatterplotResolution[1] = M;
    ScatterplotResolution[2] = 1;
    this->Modified();
  }

protected:
  ttkContinuousScatterPlot();
  ~ttkContinuousScatterPlot() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool WithDummyValue{false};
  double DummyValue{0};
  bool ProjectImageSupport{true};
  int ScatterplotResolution[3]{1920, 1080, 0};

  template <typename dataType1, class triangulationType>
  int dispatch(const dataType1 *scalars1,
               vtkDataArray *inputScalars2,
               const triangulationType *triangulation);
};
