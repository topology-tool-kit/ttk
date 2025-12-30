/// \ingroup vtk
/// \class ttkVectorWeightCurve
/// \author Tanner Finken
/// \author Josh A. Levine
/// \date 2024.
///
/// \brief TTK VTK-filter for the computation of vector field weight curves.
///
/// This filter takes a 2D vector field as input and computes the
/// number of pairs as a function of weight (i.e. the number of
/// pairs whose simplification weight is higher than a threshold).
///
/// These curves provide useful visual clues in order to fine-tune
/// simplification thresholds of a discrete vector field.
///
/// \param Input Input Vector Field
/// \param Output Table giving the number of all simplified pairs
/// as a function of weight(simplifying value) (vtkTable)
///
/// \b Related \b publication \n
/// "Localized Evaluation for Constructing Discrete Vector Fields" \n
/// Tanner Finken, Julien Tierny, Joshua A. Levine \n
/// IEEE VIS 2024.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/discreteVectorFieldTopology/">Discrete
///   Vector Field Topology example</a> \n
///

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkTable.h>

// VTK Module
#include <ttkVectorWeightCurveModule.h>

// ttk code includes
#include <VectorWeightCurve.h>
#include <ttkAlgorithm.h>

class TTKVECTORWEIGHTCURVE_EXPORT ttkVectorWeightCurve
  : public ttkAlgorithm,
    protected ttk::VectorWeightCurve {

public:
  static ttkVectorWeightCurve *New();

  vtkTypeMacro(ttkVectorWeightCurve, ttkAlgorithm);

  vtkTable *GetOutput();
  vtkTable *GetOutput(int);

  vtkSetMacro(DisplayOrbits, bool);
  vtkGetMacro(DisplayOrbits, bool);

  vtkSetMacro(DisplayExtrema, bool);
  vtkGetMacro(DisplayExtrema, bool);

  vtkSetMacro(ReverseFullOrbit, bool);
  vtkGetMacro(ReverseFullOrbit, bool);

protected:
  ttkVectorWeightCurve();
  ~ttkVectorWeightCurve() override = default;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;
};
