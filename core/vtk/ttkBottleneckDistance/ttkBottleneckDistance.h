/// \ingroup vtk
/// \class ttkBottleneckDistance
/// \author Maxime Soler <soler.maxime@total.com>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the bottleneckDistance processing package.
///
/// VTK wrapping code for the ttk::BottleneckDistance package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa ttk::BottleneckDistance
#pragma once

// TTK includes
#include <BottleneckDistance.h>
#include <ttkAlgorithm.h>
#include <ttkBottleneckDistanceModule.h>

class vtkUnstructuredGrid;

class TTKBOTTLENECKDISTANCE_EXPORT ttkBottleneckDistance
  : public ttkAlgorithm,
    protected ttk::BottleneckDistance {

public:
  static ttkBottleneckDistance *New();

  vtkTypeMacro(ttkBottleneckDistance, ttkAlgorithm);

  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);

  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);

  vtkSetMacro(PX, double);
  vtkGetMacro(PX, double);

  vtkSetMacro(PY, double);
  vtkGetMacro(PY, double);

  vtkSetMacro(PZ, double);
  vtkGetMacro(PZ, double);

  vtkSetMacro(PE, double);
  vtkGetMacro(PE, double);

  vtkSetMacro(PS, double);
  vtkGetMacro(PS, double);

  vtkSetMacro(UseOutputMatching, bool);
  vtkGetMacro(UseOutputMatching, bool);

  vtkSetMacro(UsePersistenceMetric, bool);
  vtkGetMacro(UsePersistenceMetric, bool);

  vtkSetMacro(WassersteinMetric, const std::string &);
  vtkGetMacro(WassersteinMetric, std::string);

  vtkSetMacro(DistanceAlgorithm, const std::string &);
  vtkGetMacro(DistanceAlgorithm, std::string);

  vtkSetMacro(PVAlgorithm, int);
  vtkGetMacro(PVAlgorithm, int);

  vtkSetMacro(UseGeometricSpacing, bool);
  vtkGetMacro(UseGeometricSpacing, bool);

  vtkSetMacro(Spacing, double);
  vtkGetMacro(Spacing, double);

  double Getresult() {
    return this->distance_;
  }

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation * /*request*/,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

protected:
  ttkBottleneckDistance();

private:
  bool UseOutputMatching{false};
  double Spacing{0.0};
  double Alpha{1.0};
  double Tolerance{1.0};
  double PX{0.0};
  double PY{0.0};
  double PZ{0.0};
  double PE{1.0};
  double PS{1.0};
  std::string DistanceAlgorithm{};
  std::string WassersteinMetric{"2"};
  bool UsePersistenceMetric{false};
  bool UseGeometricSpacing{false};
  int PVAlgorithm{-1};
};
