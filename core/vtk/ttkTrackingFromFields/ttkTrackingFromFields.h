#pragma once

#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <TrackingFromFields.h>
#include <ttkAlgorithm.h>
#include <ttkTrackingFromFieldsModule.h>

#include <algorithm>
#include <string>

class TTKTRACKINGFROMFIELDS_EXPORT ttkTrackingFromFields
  : public ttkAlgorithm,
    protected ttk::TrackingFromFields {

public:
  static ttkTrackingFromFields *New();

  vtkTypeMacro(ttkTrackingFromFields, ttkAlgorithm);

  vtkSetMacro(Sampling, int);
  vtkGetMacro(Sampling, int);

  vtkSetMacro(StartTimestep, int);
  vtkGetMacro(StartTimestep, int);

  vtkSetMacro(EndTimestep, int);
  vtkGetMacro(EndTimestep, int);

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

  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);

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
  vtkSetMacro(DoPostProc, bool);
  vtkGetMacro(DoPostProc, bool);

  vtkSetMacro(PostProcThresh, double);
  vtkGetMacro(PostProcThresh, double);

protected:
  ttkTrackingFromFields();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // Sampling config.
  int StartTimestep{0};
  int EndTimestep{-1};
  int Sampling{1};

  // Filtering config.
  double Tolerance{1};
  double PX{1};
  double PY{1};
  double PZ{0};
  double PE{0};
  double PS{0};

  // Bottleneck config.
  bool UseGeometricSpacing{false};
  bool DoPostProc{false};
  double PostProcThresh{0.0};
  double Spacing{1.0};
  double Alpha{1.0};
  std::string DistanceAlgorithm{"ttk"};
  int PVAlgorithm{-1};
  std::string WassersteinMetric{"2"};

  template <class dataType, class triangulationType>
  int trackWithPersistenceMatching(vtkDataSet *input,
                                   vtkUnstructuredGrid *output,
                                   unsigned long fieldNumber,
                                   const triangulationType *triangulation);
};
