/// \ingroup vtk
/// \class ttkRipsPersistenceDiagram
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.
///
/// \brief TTK VTK-filter that wraps the ttk::RipsPersistenceDiagram module.
///
/// VTK wrapping code for the ttk::RipsPersistenceDiagram package.
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
/// \sa ttk::RipsPersistenceDiagram
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkRipsPersistenceDiagramModule.h>

// VTK Includes
#include <vtkUnstructuredGrid.h>
#include <ttkMacros.h>

// TTK Includes
#include <RipsPersistenceDiagram.h>
#include <ttkAlgorithm.h>

void DiagramToVTU(
  vtkUnstructuredGrid *vtu,
  const std::vector<ttk::rpd::Diagram> &diagram,
  double SimplexMaximumDiameter);

class TTKRIPSPERSISTENCEDIAGRAM_EXPORT ttkRipsPersistenceDiagram
  : public ttkAlgorithm, // we inherit from the generic ttkAlgorithm class
    protected ttk::RipsPersistenceDiagram { // and we inherit from the base
                                            // class
private:
  bool KeepAllDataArrays{true};
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  std::vector<std::string> ScalarFields{};

public:
  static ttkRipsPersistenceDiagram *New();
  vtkTypeMacro(ttkRipsPersistenceDiagram, ttkAlgorithm);

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

  vtkSetMacro(SimplexMaximumDimension, int);
  vtkGetMacro(SimplexMaximumDimension, int);

  vtkSetMacro(SimplexMaximumDiameter, double);
  vtkGetMacro(SimplexMaximumDiameter, double);

  vtkSetMacro(FieldOfCoefficients, int);
  vtkGetMacro(FieldOfCoefficients, int);

  vtkSetMacro(InputIsDistanceMatrix, int);
  vtkGetMacro(InputIsDistanceMatrix, int);

protected:
  ttkRipsPersistenceDiagram();
  ~ttkRipsPersistenceDiagram() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};