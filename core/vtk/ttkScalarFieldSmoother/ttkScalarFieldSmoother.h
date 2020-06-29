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
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// The name of the input scalar field to consider can be specified with the
/// standard VTK call SetInputArrayToProcess(), with the following parameters:
/// \param port 0
/// \param connection 0
/// \param fieldAssociation 0 (point data)
/// \param arrayName (const char* string representing the name of the VTK array)
///
/// This module respects the following convention regarding the order of the
/// input arrays to process (SetInputArrayToProcess()):
/// \param idx 0: input scalar field
///
/// See the corresponding standalone program for a usage example:
///   - standalone/ScalarFieldSmoother/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkGeometrySmoother
/// \sa ttk::ScalarFieldSmoother
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
  vtkDataArray *outputScalarField_{nullptr};
};

#endif // _TTK_SCALAR_FIELD_SMOOTHER_H
