/// \ingroup vtk
/// \class ttkDiscreteGradient
/// \author Guillaume Favelier <guillaume.favelier@sorbonne-universite.fr>
/// \date April 2018.
///
/// \brief TTK VTK-filter that wraps the discreteGradient processing package.
///
/// VTK wrapping code for the ttk::dcg::DiscreteGradient package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The optional offset array can be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the offset array)
/// \note: To use this optional array, `ForceInputOffsetScalarField` needs to be
/// enabled with the setter `setForceInputOffsetScalarField()'.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa ttk::DiscreteGradient

#pragma once

// VTK Module
#include <ttkDiscreteGradientModule.h>

// ttk code includes
#include <DiscreteGradient.h>
#include <ttkAlgorithm.h>

class vtkPolyData;

class TTKDISCRETEGRADIENT_EXPORT ttkDiscreteGradient
  : public ttkAlgorithm,
    protected ttk::dcg::DiscreteGradient {

public:
  static ttkDiscreteGradient *New();
  vtkTypeMacro(ttkDiscreteGradient, ttkAlgorithm);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(ComputeGradientGlyphs, bool);
  vtkGetMacro(ComputeGradientGlyphs, bool);

protected:
  ttkDiscreteGradient();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  template <typename scalarType, typename triangulationType>
  int fillCriticalPoints(vtkPolyData *output,
                         vtkDataArray *const inputScalars,
                         const triangulationType &triangulation);

  template <typename triangulationType>
  int fillGradientGlyphs(vtkPolyData *const outputGradientGlyphs,
                         const triangulationType &triangulation);

  bool ForceInputOffsetScalarField{false};
  bool ComputeGradientGlyphs{true};
};
