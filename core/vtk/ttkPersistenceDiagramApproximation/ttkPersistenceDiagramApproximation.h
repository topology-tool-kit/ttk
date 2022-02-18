/// \ingroup vtk
/// \class ttkPersistenceDiagramApproximation
/// \author Jules Vidal <julien.tierny@lip6.fr>
/// \date 2022.
///
/// \brief TTK VTK-filter for the computation of an approximation of a
/// persistence diagram.
///
/// This filter computes an *approximation* of the persistence diagram of
/// the extremum-saddle pairs of an input scalar field.
/// The approximation comes with a user-controlled error on the relative
/// Bottleneck distance to the exact diagram.
///
/// \param Input Input scalar field, either 2D or 3D. Must be a regular grid.
//
/// \param Output The output of this filter is composed of:\n
/// 1. The approximation of the persistence diagram (VTU).
/// 2. The corresponding approximation of the scalar field (VTI).
/// 3. Uncertainty boumds on the correct localization of persistence pairs.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \b Related \b publication \n
/// "Fast Approximation of Persistence Diagrams with Guarantees" \n
/// Jules Vidal, Julien Tierny\n
/// IEEE Symposium on Large Data Visualization and Analysis (LDAV), 2021
///
/// \sa ttk::ApproximateTopology
/// \sa ttk::PersistenceDiagram
///
/// \b Online \b examples: \n
///

#pragma once

// VTK includes
#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkPersistenceDiagramApproximationModule.h>

// ttk code includes
#include <PersistenceDiagram.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class TTKPERSISTENCEDIAGRAMAPPROXIMATION_EXPORT
  ttkPersistenceDiagramApproximation : public ttkAlgorithm,
                                       protected ttk::PersistenceDiagram {

public:
  static ttkPersistenceDiagramApproximation *New();

  vtkTypeMacro(ttkPersistenceDiagramApproximation, ttkAlgorithm);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

  vtkSetMacro(ComputeSaddleConnectors, bool);
  vtkGetMacro(ComputeSaddleConnectors, bool);

  vtkSetMacro(ShowInsideDomain, bool);
  vtkGetMacro(ShowInsideDomain, bool);

  vtkGetMacro(StartingResolutionLevel, int);
  vtkSetMacro(StartingResolutionLevel, int);

  vtkGetMacro(StoppingResolutionLevel, int);
  vtkSetMacro(StoppingResolutionLevel, int);

  vtkGetMacro(IsResumable, bool);
  vtkSetMacro(IsResumable, bool);

  vtkGetMacro(TimeLimit, double);
  vtkSetMacro(TimeLimit, double);

  vtkGetMacro(Epsilon, double);
  vtkSetMacro(Epsilon, double);

protected:
  ttkPersistenceDiagramApproximation();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

private:
  template <typename scalarType, typename triangulationType>
  int dispatch(vtkUnstructuredGrid *outputCTPersistenceDiagram,
               vtkUnstructuredGrid *outputBounds,
               vtkDataArray *const inputScalarsArray,
               const scalarType *const inputScalars,
               scalarType *outputScalars,
               SimplexId *outputOffsets,
               int *outputMonotonyOffsets,
               const SimplexId *const inputOrder,
               const triangulationType *triangulation);

  template <typename scalarType, typename triangulationType>
  int setPersistenceDiagram(vtkUnstructuredGrid *outputCTPersistenceDiagram,
                            const std::vector<ttk::PersistencePair> &diagram,
                            vtkDataArray *inputScalarsArray,
                            const scalarType *const inputScalars,
                            const triangulationType *triangulation) const;

  template <typename scalarType, typename triangulationType>
  int drawBottleneckBounds(vtkUnstructuredGrid *outputBounds,
                           const std::vector<ttk::PersistencePair> &diagram,
                           vtkDataArray *inputScalarsArray,
                           const scalarType *const outputScalars,
                           const scalarType *const inputScalars,
                           const triangulationType *triangulation) const;

  bool ForceInputOffsetScalarField{false};
  bool ShowInsideDomain{false};
};
