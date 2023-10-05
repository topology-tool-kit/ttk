/// \ingroup vtk
/// \class ttkPersistenceCurve
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date September 2016.
///
/// \brief TTK VTK-filter for the computation of persistence curves.
///
/// This filter takes a persistence diagram as input and computes the
/// number of pairs as a function of persistence (i.e. the number of
/// pairs whose persistence is higher than a threshold).
///
/// These curves provide useful visual clues in order to fine-tune persistence
/// simplification thresholds.
///
/// \param Input Input Persistence Diagram as computed by the
/// \ref ttkPersistenceDiagram module
/// \param Output0 Table giving the number of persistent minimum-saddle pairs
/// as a function of persistence (vtkTable)
/// \param Output1 Table giving the number of persistent saddle-saddle pairs
/// as a function of persistence (vtkTable)
/// \param Output2 Table giving the number of persistent saddle-maximum pairs
/// as a function of persistence (vtkTable)
/// \param Output3 Table giving the number of all persistent pairs as
/// a function of persistence (vtkTable)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttkMergeAndContourTreePP
/// \sa ttkPersistenceDiagram
/// \sa ttkTopologicalSimplification
/// \sa ttk::PersistenceCurve
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/BuiltInExample1/">BuiltInExample1
/// </a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/dragon/">Dragon
/// example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morsePersistence/">Morse
///   Persistence example</a> \n
///

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkTable.h>

// VTK Module
#include <ttkPersistenceCurveModule.h>

// ttk code includes
#include <PersistenceCurve.h>
#include <ttkAlgorithm.h>

class TTKPERSISTENCECURVE_EXPORT ttkPersistenceCurve
  : public ttkAlgorithm,
    protected ttk::PersistenceCurve {

public:
  static ttkPersistenceCurve *New();

  vtkTypeMacro(ttkPersistenceCurve, ttkAlgorithm);

  vtkTable *GetOutput();
  vtkTable *GetOutput(int);

protected:
  ttkPersistenceCurve();
  ~ttkPersistenceCurve() override = default;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;
};
