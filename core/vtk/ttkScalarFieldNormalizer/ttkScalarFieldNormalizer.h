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
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa vtkCellDataConverter
/// \sa vtkPointDataConverter
/// \sa vtkScalarFieldSmoother
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/harmonicSkeleton/">
///   Harmonic Skeleton example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n
///

#pragma once

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
