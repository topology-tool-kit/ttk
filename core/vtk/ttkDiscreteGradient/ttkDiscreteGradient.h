/// \ingroup vtk
/// \class ttkDiscreteGradient
/// \author Guillaume Favelier <guillaume.favelier@sorbonne-universite.fr>
/// \date April 2018.
///
/// \brief TTK VTK-filter that wraps the discreteGradient processing package.
///
/// VTK wrapping code for the @DiscreteGradient package.
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
/// \sa ttk::DiscreteGradient

#pragma once

// VTK Module
#include <ttkDiscreteGradientModule.h>

// ttk code includes
#include <DiscreteGradient.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class vtkUnstructuredGrid;

class TTKDISCRETEGRADIENT_EXPORT ttkDiscreteGradient
  : public ttkAlgorithm,
    protected ttk::dcg::DiscreteGradient {

public:
  static ttkDiscreteGradient *New();
  vtkTypeMacro(ttkDiscreteGradient, ttkAlgorithm);

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(InputOffsetScalarFieldName, std::string);
  vtkGetMacro(InputOffsetScalarFieldName, std::string);

  vtkSetMacro(ComputeGradientGlyphs, bool);
  vtkGetMacro(ComputeGradientGlyphs, bool);

  vtkSetMacro(IterationThreshold, int);
  vtkGetMacro(IterationThreshold, int);

  vtkSetMacro(ScalarFieldId, int);
  vtkGetMacro(ScalarFieldId, int);

  vtkSetMacro(OffsetFieldId, int);
  vtkGetMacro(OffsetFieldId, int);

  int getScalars(vtkDataSet *input);
  int getOffsets(vtkDataSet *input);

protected:
  ttkDiscreteGradient();
  ~ttkDiscreteGradient() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  template <typename T>
  int dispatch(vtkUnstructuredGrid *outputCriticalPoints);

  std::string ScalarField{};
  std::string InputOffsetScalarFieldName{};
  bool ForceInputOffsetScalarField{false};
  bool ComputeGradientGlyphs{true};
  int IterationThreshold{-1};
  int ScalarFieldId{};
  int OffsetFieldId{-1};

  vtkDataArray *inputScalars_{};
  ttkSimplexIdTypeArray *offsets_{};
  vtkDataArray *inputOffsets_{};
  bool hasUpdatedMesh_{};
};
