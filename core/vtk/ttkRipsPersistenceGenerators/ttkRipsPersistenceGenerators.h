/// \ingroup vtk
/// \class ttkRipsPersistenceGenerators
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date June 2024.
///
/// \brief TTK VTK-filter that wraps the ttk::RipsPersistenceGenerators module.
///
/// VTK wrapping code for the ttk::RipsPersistenceGenerators package.
///
/// \param Input Input table (vtkTable)
/// \param Input Point set for geometric realization (vtkPointSet)
/// \param Output Generators (vtkUnstructuredGrid)
/// \param Output PersistenceDiagram (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::RipsPersistenceGenerators
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkRipsPersistenceGeneratorsModule.h>

// TTK Base Includes
#include <RipsPersistenceGenerators.h>
#include <ttkAlgorithm.h>

class TTKRIPSPERSISTENCEGENERATORS_EXPORT ttkRipsPersistenceGenerators
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  , protected ttk::RipsPersistenceGenerators // and we inherit from the base class
{
private:
  bool KeepAllDataArrays{true};
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  std::vector<std::string> ScalarFields{};

public:
  static ttkRipsPersistenceGenerators *New();
  vtkTypeMacro(ttkRipsPersistenceGenerators, ttkAlgorithm);

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

  vtkSetMacro(SimplexMaximumDiameter, double);
  vtkGetMacro(SimplexMaximumDiameter, double);

  vtkSetMacro(OutputCascade, bool);
  vtkGetMacro(OutputCascade, bool);

protected:
  ttkRipsPersistenceGenerators();
  ~ttkRipsPersistenceGenerators() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
