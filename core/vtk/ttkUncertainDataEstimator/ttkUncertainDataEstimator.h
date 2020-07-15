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
/// \param InputN Input ensemble scalar field #N (vtkDataSet)
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

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKUNCERTAINDATAESTIMATOR_EXPORT ttkUncertainDataEstimator
  : public ttkAlgorithm,
    protected ttk::UncertainDataEstimator {

public:
  static ttkUncertainDataEstimator *New();

  vtkTypeMacro(ttkUncertainDataEstimator, vtkAlgorithm);

  // set-getters macros to define from each variable you want to access from
  // the outside (in particular from paraview) - to adapt.

  void ComputeLowerBound(bool state) {
    computeLowerBound_ = state;
  }

  void ComputeUpperBound(bool state) {
    computeUpperBound_ = state;
  }

  void BoundToCompute(int value) {
    switch(value) {
      case 0:
        ComputeLowerBound(true);
        ComputeUpperBound(true);
        break;
      case 1:
        ComputeLowerBound(true);
        ComputeUpperBound(false);
        break;
      case 2:
        ComputeLowerBound(false);
        ComputeUpperBound(true);
        break;
    }
    Modified();
  }

  void BinCount(int binCount) {
    binCount_ = binCount;
    Modified();
  }

protected:
  ttkUncertainDataEstimator();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool computeLowerBound_{true};
  bool computeUpperBound_{true};
  int binCount_{10};
};
