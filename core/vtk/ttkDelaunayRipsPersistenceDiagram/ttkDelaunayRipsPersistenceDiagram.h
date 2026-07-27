/// \ingroup vtk
/// \class ttkDelaunayRipsPersistenceDiagram
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date September 2025.
///
/// \brief TTK VTK-filter that wraps the ttk::DelaunayRipsPersistenceDiagram
/// module.
///
/// VTK wrapping code for the ttk::DelaunayRipsPersistenceDiagram package.
///
/// \param Input Input table (vtkTable)
/// \param Output PersistenceDiagram (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::DelaunayRipsPersistenceDiagram
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/delaunayRispPersistence/">DelaunayRips
///   Persistence example</a> \n

#pragma once

// VTK Module
#include <ttkDelaunayRipsPersistenceDiagramModule.h>

// VTK Includes
#include <ttkMacros.h>
#include <vtkUnstructuredGrid.h>

// TTK Includes
#include <DelaunayRipsPersistenceDiagram.h>
#include <ttkAlgorithm.h>

class TTKDELAUNAYRIPSPERSISTENCEDIAGRAM_EXPORT ttkDelaunayRipsPersistenceDiagram
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
  static ttkDelaunayRipsPersistenceDiagram *New();
  vtkTypeMacro(ttkDelaunayRipsPersistenceDiagram, ttkAlgorithm);

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
  ttkDelaunayRipsPersistenceDiagram();
  ~ttkDelaunayRipsPersistenceDiagram() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
