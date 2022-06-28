/// \ingroup vtk
/// \class ttkPersistentGenerators
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date February 2022.
///
/// \brief TTK VTK-filter for the computation of persistent generators.

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
