/// \ingroup vtk
/// \class ttkPersistentGenerators
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date February 2022.
///
/// \brief TTK VTK-filter for the computation of persistent generators.
///
/// \b Related \b publication \n
/// "Discrete Morse Sandwich: Fast Computation of Persistence Diagrams for
/// Scalar Data -- An Algorithm and A Benchmark" \n
/// Pierre Guillou, Jules Vidal, Julien Tierny \n
/// IEEE Transactions on Visualization and Computer Graphics, 2023.\n
/// arXiv:2206.13932, 2023.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_at/">Persistent
///   Generators AT example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_casting/">Persistent
///   Generators Casting example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_darkSky/">Persistent
///   Generators DarkSky example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_fertility/">Persistent
///   Generators Fertility example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_householdAnalysis/">Persistent
///   Generators Household Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_periodicPicture/">Persistent
///   Generators Periodic Picture example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_skull/">Persistent
///   Generators Skull example</a> \n
///

#pragma once

// VTK includes
#include <vtkDataArray.h>
#include <vtkPolyData.h>

// VTK Module
#include <ttkPersistentGeneratorsModule.h>

// ttk code includes
#include <PersistentGenerators.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class TTKPERSISTENTGENERATORS_EXPORT ttkPersistentGenerators
  : public ttkAlgorithm,
    protected ttk::PersistentGenerators {

public:
  static ttkPersistentGenerators *New();

  vtkTypeMacro(ttkPersistentGenerators, ttkAlgorithm);

  vtkSetMacro(PruneHandlesGenerators, bool);
  vtkGetMacro(PruneHandlesGenerators, bool);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

protected:
  ttkPersistentGenerators();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

private:
  template <typename triangulationType>
  int dispatch(vtkPolyData *output,
               vtkDataArray *const inputScalarsArray,
               const SimplexId *const inputOrder,
               const triangulationType &triangulation);

  bool ForceInputOffsetScalarField{false};
};
