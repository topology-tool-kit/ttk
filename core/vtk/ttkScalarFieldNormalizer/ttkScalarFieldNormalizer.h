/// \ingroup vtk
/// \class ttkScalarFieldNormalizer
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date May 2016.
///
/// \brief TTK VTK-filter that normalizes an input scalar field.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output normalized scalar field (vtkDataSet)
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
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkCellDataConverter
/// \sa vtkPointDataConverter
/// \sa vtkScalarFieldSmoother
#pragma once

// VTK includes -- to adapt

// VTK Module
#include <ttkScalarFieldNormalizerModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKSCALARFIELDNORMALIZER_EXPORT ttkScalarFieldNormalizer
  : public ttkAlgorithm {

public:
  static ttkScalarFieldNormalizer *New();

  vtkTypeMacro(ttkScalarFieldNormalizer, ttkAlgorithm);

protected:
  ttkScalarFieldNormalizer();

  ~ttkScalarFieldNormalizer() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int normalize(vtkDataArray *input, vtkDataArray *output) const;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
};
