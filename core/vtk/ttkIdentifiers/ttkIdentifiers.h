/// \ingroup vtk
/// \class ttkIdentifiers
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date August 2015.
///
/// \brief TTK VTK-filter that computes the global identifiers for each vertex
/// and each cell as point data and cell data scalar fields.
///
/// This filter is useful to retrieve the global identifiers of vertices or
/// cells in subsequent filters throughout the VTK pipeline.
///
/// \param Input Input data-set (vtkDataSet)
/// \param Output Output data-set with identifier fields (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/imageProcessing/">Image
///   Processing example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tectonicPuzzle/">Tectonic
///   Puzzle example</a> \n
///

#pragma once

// VTK includes -- to adapt

// VTK Module
#include <ttkIdentifiersModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKIDENTIFIERS_EXPORT ttkIdentifiers : public ttkAlgorithm {

public:
  static ttkIdentifiers *New();

  vtkTypeMacro(ttkIdentifiers, ttkAlgorithm);

  vtkSetMacro(CellFieldName, const std::string &);
  vtkGetMacro(CellFieldName, std::string);

  vtkSetMacro(VertexFieldName, const std::string &);
  vtkGetMacro(VertexFieldName, std::string);

protected:
  ttkIdentifiers();

  ~ttkIdentifiers() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string CellFieldName{"CellIdentifiers"},
    VertexFieldName{ttk::VertexScalarFieldName};
};
