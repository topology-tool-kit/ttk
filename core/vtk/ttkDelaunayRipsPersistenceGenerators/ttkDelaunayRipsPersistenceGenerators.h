/// \ingroup vtk
/// \class ttkDelaunayRipsPersistenceGenerators
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date September 2025.
///
/// \brief TTK VTK-filter that wraps the ttk::DelaunayRipsPersistenceGenerators
/// module.
///
/// VTK wrapping code for the ttk::DelaunayRipsPersistenceGenerators package.
///
/// \param Input Input table (vtkTable)
/// \param Output PersistenceDiagram (vtkUnstructuredGrid)
/// \param Output 1-dimensional generators (vtkUnstructuredGrid)
/// \param Output 2-dimensional generators (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::DelaunayRipsPersistenceGenerators
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/delaunayRispPersistence/">DelaunayRips
///   Persistence example</a> \n

#pragma once

// VTK Module
#include <ttkDelaunayRipsPersistenceGeneratorsModule.h>

// VTK Includes
#include <ttkMacros.h>
#include <vtkUnstructuredGrid.h>

// TTK Includes
#include <DelaunayRipsPersistenceDiagram.h>
#include <ttkAlgorithm.h>

void GeneratorsToVTU(vtkUnstructuredGrid *vtu,
                     vtkPoints *inputPoints,
                     const std::vector<ttk::rpd::Generator2> &generators);

class TTKDELAUNAYRIPSPERSISTENCEGENERATORS_EXPORT
  ttkDelaunayRipsPersistenceGenerators
  : public ttkAlgorithm, // we inherit from the generic ttkAlgorithm class
    protected ttk::DelaunayRipsPersistenceDiagram { // and we inherit from the
                                                    // base
  // class
private:
  bool KeepAllDataArrays{true};
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  std::vector<std::string> ScalarFields{};

public:
  static ttkDelaunayRipsPersistenceGenerators *New();
  vtkTypeMacro(ttkDelaunayRipsPersistenceGenerators, ttkAlgorithm);

  void SetScalarFields(const std::string &s) {
    ScalarFields.push_back(s);
    Modified();
  }

  void ClearScalarFields() {
    ScalarFields.clear();
    Modified();
  }

  vtkSetMacro(KeepAllDataArrays, bool);
  vtkGetMacro(KeepAllDataArrays, bool);

  vtkSetMacro(SelectFieldsWithRegexp, bool);
  vtkGetMacro(SelectFieldsWithRegexp, bool);

  vtkSetMacro(RegexpString, const std::string &);
  vtkGetMacro(RegexpString, std::string);

protected:
  ttkDelaunayRipsPersistenceGenerators();
  ~ttkDelaunayRipsPersistenceGenerators() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
