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

#ifndef _TTK_CONTINUOUSSCATTERPLOT_H
#define _TTK_CONTINUOUSSCATTERPLOT_H

// VTK includes
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// ttk baseCode includes
#include <ContinuousScatterPlot.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkContinuousScatterPlot
#else
class ttkContinuousScatterPlot
#endif
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkContinuousScatterPlot *New();

  vtkTypeMacro(ttkContinuousScatterPlot, vtkDataSetAlgorithm);

  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }

  vtkSetMacro(ScalarField1, std::string);
  vtkGetMacro(ScalarField1, std::string);

  vtkSetMacro(ScalarField2, std::string);
  vtkGetMacro(ScalarField2, std::string);

  vtkSetMacro(UcomponentId, int);
  vtkGetMacro(UcomponentId, int);

  vtkSetMacro(VcomponentId, int);
  vtkGetMacro(VcomponentId, int);

  vtkSetMacro(WithVaryingConnectivity, int);
  vtkGetMacro(WithVaryingConnectivity, int);

  vtkSetMacro(WithDummyValue, int);
  vtkGetMacro(WithDummyValue, int);

  vtkSetMacro(DummyValue, double);
  vtkGetMacro(DummyValue, double);

  vtkSetMacro(ProjectImageSupport, bool);
  vtkGetMacro(ProjectImageSupport, bool);

  void SetScatterplotResolution(int N, int M) {
    ScatterplotResolution[0] = N;
    ScatterplotResolution[1] = M;
    ScatterplotResolution[2] = 1;
    Modified();
  }

  int getScalars(vtkDataSet *input);
  int getTriangulation(vtkDataSet *input);

protected:
  ttkContinuousScatterPlot();
  ~ttkContinuousScatterPlot();

  TTK_SETUP();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

private:
  bool WithVaryingConnectivity;
  bool WithDummyValue;
  double DummyValue;
  bool ProjectImageSupport;
  int ScatterplotResolution[3];
  int UcomponentId, VcomponentId;
  std::string ScalarField1;
  std::string ScalarField2;

  vtkDataArray *inputScalars1_;
  vtkDataArray *inputScalars2_;
  double scalarMin_[2];
  double scalarMax_[2];
  std::vector<std::vector<double>> density_;
  std::vector<std::vector<char>> validPointMask_;
  ttk::Triangulation *triangulation_;

  // output
  vtkSmartPointer<vtkUnstructuredGrid> vtu_;
  vtkSmartPointer<vtkPoints> pts_;
};

#endif // _TTK_CONTINUOUSSCATTERPLOT_H
