/// \ingroup vtk
/// \class ttkUncertainDataEstimator
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \date August 2016.
///
/// \brief TTK VTK-filter that takes an input ensemble data set
/// (represented by a list of scalar fields) and which computes various
/// vertexwise statistics (PDF estimation, bounds, moments, etc.)
///
/// \param Input0 Input ensemble scalar field #0 (vtkDataSet)
/// \param Input1 Input ensemble scalar field #1 (vtkDataSet)\n
/// ...\n
/// \param InputN Input ensemble scalar field #%N (vtkDataSet)
/// \param Output0 Lower and upper bound fields (vtkDataSet)
/// \param Output1 Histogram estimations of the vertex probability density
/// functions (vtkDataSet)
/// \param Output2 Mean field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa vtkMandatoryCriticalPoints
/// \sa ttk::UncertainDataEstimator

#pragma once

// VTK Module
#include <ttkAlgorithm.h>
#include <ttkUncertainDataEstimatorModule.h>

// ttk code includes
#include <UncertainDataEstimator.h>

class TTKUNCERTAINDATAESTIMATOR_EXPORT ttkUncertainDataEstimator
  : public ttkAlgorithm,
    protected ttk::UncertainDataEstimator {

public:
  static ttkUncertainDataEstimator *New();

  vtkTypeMacro(ttkUncertainDataEstimator, ttkAlgorithm);

  vtkGetMacro(ComputeLowerBound, bool);
  vtkSetMacro(ComputeLowerBound, bool);

  vtkGetMacro(ComputeUpperBound, bool);
  vtkSetMacro(ComputeUpperBound, bool);

  vtkGetMacro(BinCount, int);
  vtkSetMacro(BinCount, int);

  void SetBoundToCompute(int value) {
    if(value == 0) {
      SetComputeLowerBound(true);
      SetComputeUpperBound(true);
    } else if(value == 1) {
      SetComputeLowerBound(true);
      SetComputeUpperBound(false);
    } else if(value == 2) {
      SetComputeLowerBound(false);
      SetComputeUpperBound(true);
    }
  }

protected:
  ttkUncertainDataEstimator();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
