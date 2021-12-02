/// \ingroup vtk
/// \class ttkScalarFieldSmoother
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief TTK VTK-filter for scalar field smoothing.
///
/// This class is a dummy example for the development of TTK filters. It
/// smooths an input scalar field by averaging the scalar values on the link
/// of each vertex.
///
/// \param Input vtkDataSet
/// \param Output vtkDataSet
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The optional mask array can be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the mask array)
/// \note: To use this optional array, `ForceInputMaskScalarField` needs to be
/// enabled with the setter `setForceInputMaskScalarField()'.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding standalone program for a usage example:
///   - standalone/ScalarFieldSmoother/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkGeometrySmoother
/// \sa ttk::ScalarFieldSmoother
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morsePersistence/">Morse
///   Persistence example</a> \n

#ifndef _TTK_SCALAR_FIELD_SMOOTHER_H
#define _TTK_SCALAR_FIELD_SMOOTHER_H

// VTK includes

// VTK Module
#include <ttkScalarFieldSmootherModule.h>

// ttk code includes
#include <ScalarFieldSmoother.h>
#include <ttkAlgorithm.h>

class TTKSCALARFIELDSMOOTHER_EXPORT ttkScalarFieldSmoother
  : public ttkAlgorithm,
    protected ttk::ScalarFieldSmoother {

public:
  static ttkScalarFieldSmoother *New();

  vtkTypeMacro(ttkScalarFieldSmoother, ttkAlgorithm);

  vtkSetMacro(NumberOfIterations, int);
  vtkGetMacro(NumberOfIterations, int);

  vtkSetMacro(ForceInputMaskScalarField, bool);
  vtkGetMacro(ForceInputMaskScalarField, bool);

protected:
  ttkScalarFieldSmoother();

  ~ttkScalarFieldSmoother() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int NumberOfIterations{1};
  bool ForceInputMaskScalarField{false};
};

#endif // _TTK_SCALAR_FIELD_SMOOTHER_H
