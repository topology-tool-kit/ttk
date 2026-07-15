/// \ingroup vtk
/// \class ttkDiscreteVectorField
/// \author Tanner Finken <finkent@arizona.edu>
/// \author Joshua A. Levine <josh@cs.arizona.edu>
/// \date May 2024.
///
/// \brief TTK VTK-filter that wraps the discreteVectorField processing package.
///
/// VTK wrapping code for the ttk::dcvf::DiscreteVectorField package.
///
/// \b Related \b publication \n
/// "Localized Evaluation for Constructing Discrete Vector Fields" \n
/// Tanner Finken, Julien Tierny, Joshua A. Levine \n
/// IEEE Vis 2024.
///
/// \param Input Input vector field (vtkDataSet)
/// \param Output Output glyphs (vtkPolyData)
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \sa ttk::dcvf::DiscreteVectorField
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/discreteVectorFieldTopology/">Discrete
///   Vector Field Topology example</a> \n
///

#pragma once

// VTK Module
#include <ttkDiscreteVectorFieldModule.h>

// ttk code includes
#include <DiscreteVectorField.h>
#include <ttkAlgorithm.h>

class vtkPolyData;

class TTKDISCRETEVECTORFIELD_EXPORT ttkDiscreteVectorField
  : public ttkAlgorithm,
    protected ttk::dcvf::DiscreteVectorField {

public:
  static ttkDiscreteVectorField *New();
  vtkTypeMacro(ttkDiscreteVectorField, ttkAlgorithm);

  vtkSetMacro(ComputeVectorGlyphs, bool);
  vtkGetMacro(ComputeVectorGlyphs, bool);

protected:
  ttkDiscreteVectorField();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  template <typename scalarType, typename triangulationType>
  int fillCriticalPoints(vtkPolyData *output,
                         const triangulationType &triangulation);

  template <typename triangulationType>
  int fillVectorGlyphs(vtkPolyData *const outputVectorGlyphs,
                       const triangulationType &triangulation);

  bool ComputeVectorGlyphs{true};
};
